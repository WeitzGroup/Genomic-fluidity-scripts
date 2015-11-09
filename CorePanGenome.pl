#!/usr/bin/env perl

# Fluidity Estimator
# Copyright (C) 2011 Georgia Institute of Technology
# Fluidity Estimator comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions.
# Contact: Joshua Weitz <jsweitz@gatech.edu>
# Contact: Andrey Kislyuk <kislyuk@gmail.com>

package CorePanGenome;

my $settings = {
	min_ortho_coverage => 0.5,
	min_ortho_identity => 0.3,
#	ortho_coverage_step => 0.02,
#	ortho_identity_step => 0.02,
#	ortho_coverage_step => 0.1,
#	ortho_identity_step => 0.1,
	ortho_coverage_step => 1,
	ortho_identity_step => 1,
	tempdir => 'cptmp', keep => 1,
	xcdir => "xcdir", # must be mounted on all nodes; TODO: use PBS staging support
	use_pbs => 0,
	pbs_queue => 'topaz_main',
};
my $stats;

use strict;
use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use FindBin;
use lib "$FindBin::RealBin/lib";
use AKUtils qw(logmsg);
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
#use Statistics::Descriptive;

$0 = fileparse($0);

exit(main());

sub main() {
	die("Usage: $0 <genome1.(gff|fasta|gb)> <genome2> [...] [options]\n") if @ARGV < 2;

	my @cmd_options = ('identity=s', 'coverage=s', 'min_ortho_coverage=s', 'min_ortho_identity=s', 'use_pbs=s', 'pbs_queue=s', 'tempdir=s', 'xcdir=s', 'ortho_coverage_step=s', 'ortho_identity_step=s');
	GetOptions($settings, @cmd_options) or die;

	if ($$settings{identity}) {
		die("Arguments --identity and --coverage must be supplied together") unless $$settings{coverage};
		die("Identity setting $$settings{identity} out of bounds") if $$settings{identity} < 0.05 or $$settings{identity} > 1;
		die("Identity setting $$settings{coverage} out of bounds") if $$settings{coverage} < 0.05 or $$settings{coverage} > 1;
		$$settings{min_ortho_identity} = $$settings{identity};
		$$settings{min_ortho_coverage} = $$settings{coverage};
		$$settings{ortho_identity_step} = 1;
		$$settings{ortho_coverage_step} = 1;
	}

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	$$settings{tempdir} = File::Spec->rel2abs($$settings{tempdir});
	$$settings{xcdir} = File::Spec->rel2abs($$settings{xcdir});
	$$settings{outfile_prefix} = "$0";
	logmsg "Temporary directory is $$settings{tempdir}";

	mkdir($$settings{tempdir}) unless -d $$settings{tempdir};
	mkdir($$settings{xcdir}) unless -d $$settings{xcdir};

	my @input_files = @ARGV;
	die("No input files supplied") if @input_files < 1;
	foreach my $file (@input_files) { $file = File::Spec->rel2abs($file); }
	chdir($$settings{tempdir}) or die("Unable to change working directory to $$settings{tempdir}");

	my $checkpoint_passed;
	if (-f "$$settings{xcdir}/$0.checkpoint") {
		$checkpoint_passed = 1;
		open(CHK, '<', "$$settings{xcdir}/$0.checkpoint") or die;
		my @chk_files = <CHK>; chomp @chk_files;
		for (0..$#input_files) {
			$checkpoint_passed = 0 if $input_files[$_] ne $chk_files[$_];
#			warn "Line $_ failed checkpointing:\n$input_files[$_]\nvs.\n$chk_files[$_]\n" if $input_files[$_] ne $chk_files[$_];
		}
		close CHK;
	}

	if ($checkpoint_passed) {
		logmsg "Loading data from checkpoint in $$settings{xcdir}...";
	} else {
#		die("checkpoint failed; full calc temporarily disabled");
		unlink "$$settings{xcdir}/$0.checkpoint";
		my $invoke_string = "$FindBin::RealBin/CorePanGenome_blaster.pl @input_files";
		for (qw(min_ortho_coverage min_ortho_identity xcdir tempdir use_pbs pbs_queue outfile_prefix)) {
			$invoke_string .= " --$_=$$settings{$_}";
		}
		logmsg "Running $invoke_string...";
		system($invoke_string); die if $?;
		open(CHK, '>', "$$settings{xcdir}/$0.checkpoint") or die;
		print CHK join("\n", @input_files)."\n";
		close CHK;
	}
	runAnalysisJobs(\@input_files, $settings);

	return 0;
}

sub runAnalysisJobs($$) {
	my ($input_files, $settings) = @_;
	logmsg "Running analysis jobs...";

	$$settings{report_dir} = "$$settings{xcdir}/report";

	my @cmd_list;
	for (my $min_ortho_identity = $$settings{min_ortho_identity}; $min_ortho_identity <= 1; $min_ortho_identity += $$settings{ortho_identity_step}) {
		for (my $min_ortho_coverage = $$settings{min_ortho_coverage}; $min_ortho_coverage <= 1; $min_ortho_coverage += $$settings{ortho_coverage_step}) {
			my $invoke_string = "$FindBin::RealBin/CorePanGenome_analyzer.pl @$input_files";
			$invoke_string .= " --min_ortho_coverage=$min_ortho_coverage";
			$invoke_string .= " --min_ortho_identity=$min_ortho_identity";
			for (qw(xcdir outfile_prefix report_dir)) {
				$invoke_string .= " --$_=$$settings{$_}";
			}

#			unlink "$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out.done";
			push(@cmd_list, $invoke_string);
		}
	}

	print "C: @cmd_list\n";

	if ($$settings{use_pbs}) {
		AKUtils::runPBSjobs(\@cmd_list, undef, {%$settings, wait_on_pbs_jobs=>1});
	} else {
		for (my $i=0; $i<@cmd_list; $i += 16) {
			my @cmd_batch = @cmd_list[$i..$i+15];
			foreach my $cmd (@cmd_batch) {
				next unless $cmd;
				logmsg "Running $cmd...";
				system("$cmd"); die if $?;
				#system("$cmd &"); die if $?;
			}
			#my $pacer_cmd = $cmd_list[$i];
			#system($pacer_cmd); die if $?;
		}
	}
	logmsg "Output is in $$settings{report_dir}";
}
