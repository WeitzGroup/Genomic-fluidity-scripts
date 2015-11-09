#!/usr/bin/env perl

# Fluidity Estimator
# Copyright (C) 2011 Georgia Institute of Technology
# Fluidity Estimator comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions.
# Contact: Joshua Weitz <jsweitz@gatech.edu>
# Contact: Andrey Kislyuk <kislyuk@gmail.com>

package CorePanGenome;

my $settings = {};
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

use CorePanGenome;

$0 = fileparse($0);

exit(main());

sub main() {
	die("Usage: $0 <genome1.(gff|fasta|gb)> <genome2> [...] [options]\n") if @ARGV < 2;

	my @cmd_options = ('min_ortho_coverage=s', 'min_ortho_identity=s', 'use_pbs=s', 'pbs_queue=s', 'tempdir=s', 'xcdir=s', 'outfile_prefix=s');
	GetOptions($settings, @cmd_options) or die;

	for (qw(min_ortho_coverage min_ortho_identity xcdir outfile_prefix)) {
		die("$0: argument $_ is required") unless defined $$settings{$_};
	}

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	$$settings{tempdir} = File::Spec->rel2abs($$settings{tempdir});
	$$settings{xcdir} = File::Spec->rel2abs($$settings{xcdir});
	logmsg "Temporary directory is $$settings{tempdir}";

	my @input_files = @ARGV;
	foreach my $file (@input_files) { $file = File::Spec->rel2abs($file); }
	chdir($$settings{tempdir}) or die("Unable to change working directory to $$settings{tempdir}");

	my $genes = CorePanGenome::loadGenesFromInput(\@input_files);

	runAlignments($genes, $settings);
	logmsg "Output is in $$settings{xcdir}";
	
	return 0;
}

sub runAlignments($$) {
	my ($genes, $settings) = @_;
	logmsg "Preparing protein data...";

	my @orgs = keys %$genes;
	foreach my $i (0..$#orgs) { # Print protein sequences to files for the worker to run BLAST on
		my $org1 = $orgs[$i];
		my $blast_db_h = Bio::SeqIO->new(-file => ">$$settings{xcdir}/$org1.aa.fasta", -format => 'Fasta');
		foreach my $gene (values %{$$genes{$org1}}) {
			$blast_db_h->write_seq(Bio::Seq->new(-seq => $gene->seq->translate()->seq(),
				-display_id => ($gene->get_tag_values('locus_tag'))[0]));
		}
	}

	my @cmd_list;
	foreach my $i (0..$#orgs) { # Form commands for the worker
		my $org1 = $orgs[$i]; # org1 is designated as database, org2 as query		
		foreach my $j ($i..$#orgs) { # Every pair, plus against self
			my $org2 = $orgs[$j]; # org1 is designated as database, org2 as query
			my $cmd = "$FindBin::RealBin/CorePanGenome_worker.pl";
			$cmd .= " --run_id='$$settings{outfile_prefix}.$org2.$org1'";
			$cmd .= " --outfile='$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out'";
			$cmd .= " --db_mfa='$$settings{xcdir}/$org1.aa.fasta' --query_mfa='$$settings{xcdir}/$org2.aa.fasta'";
			for ('min_ortho_coverage', 'min_ortho_identity') {
				$cmd .= " --$_=$$settings{$_}";
			}
			
			if (-f "$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out.done"
				and (stat("$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out"))[7] > 0) {
#				warn "would skip job $$settings{outfile_prefix}.$org2.$org1"; next;
			}

			unlink "$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out.done";
			push(@cmd_list, $cmd);
		}
	}

	if ($$settings{use_pbs}) {
		AKUtils::runPBSjobs(\@cmd_list, undef, {%$settings, wait_on_pbs_jobs=>1});
	} else {
		foreach my $cmd (@cmd_list) {
			logmsg "Running $cmd...";
			system($cmd); die if $?;
		}
	}
}
