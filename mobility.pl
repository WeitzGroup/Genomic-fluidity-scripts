#!/usr/bin/env perl

# Fluidity Estimator
# Copyright (C) 2011 Georgia Institute of Technology
# Fluidity Estimator comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions.
# Contact: Joshua Weitz <jsweitz@gatech.edu>
# Contact: Andrey Kislyuk <kislyuk@gmail.com>

package Mobility;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my %settings = (
	logfile => '-',
	datadir => $FindBin::RealBin,
	blast_db => "/storage/db/nt/nt",
);
my %stats;

use strict;
use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use Class::Struct;
use IPC::Open2;
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/lib";
use AKUtils;

sub main();

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

my $usage = "$0 input.mfa [options]
Options can be abbreviated. Options are:
";

exit(main());

sub main() {
	die("Usage: $usage") if @ARGV < 1;

	my ($input_file_full) = @ARGV;
	$input_file_full = File::Spec->rel2abs($input_file_full);
	my ($input_seq_file, $input_seq_dir) = fileparse($input_file_full);

	my @cmd_options = ('outfile=s', 'keep', 'frag_length=s');
	GetOptions(\%settings, @cmd_options) or die("Error while parsing command line switches. Check the input options. Usage:".$usage);

	my @uinfo; # = getpwuid($>);
	warn "$0: Version ".$VERSION." started ".localtime()." by ".$uinfo[0]."\n";

	$settings{outfile} ||= "$input_seq_dir/$input_seq_file.mob";
	$settings{logfile} ||= "$settings{outfile}.log";
	open(FH, '>>', $settings{outfile}) or die "Unable to open file $settings{outfile} for writing: $!"; close FH;
	if ($settings{logfile} eq '-') {
		$FSFind::LOG = *STDERR;
	} else {
		open($FSFind::LOG, '>', $settings{logfile}) or die "Unable to open file $settings{logfile} for writing: $!";
		warn "Logging to $settings{logfile}. Use \"-logfile=-\" to direct logging to stderr.\n";
	}
	$settings{keep} = 1;
	$settings{tempdir} = tempdir($settings{tempdir} or File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($settings{keep}));
	logmsg "Temporary directory is $settings{tempdir}";

	my $seqs = AKUtils::readMfa($input_file_full, {replace_spaces => 1});
	my %frags;

	exit if length((values %$seqs)[0]) < 1e6;
	$settings{frag_num} =  int(length((values %$seqs)[0])/4000);
	die("set frag_length") unless $settings{frag_length};
#	$settings{frag_length} = 500;
	
die if scalar(keys %$seqs) > 1;
	foreach my $seqname (keys %$seqs) {
		
#		die("Sequence $seqname too short") if length($$seqs{$seqname}) < $settings{frag_length} * $settings{frag_num};
		foreach my $i (1..$settings{frag_num}) {
			my $start_pos = int(rand(length($$seqs{$seqname}) - $settings{frag_length}));
			my $frag = substr($$seqs{$seqname}, $start_pos, $settings{frag_length});
			#print "$seqname $i $start_pos ".length($frag)."\n";
			$frags{$seqname.":".$start_pos} = $frag;
		}
	}

	open(BLAT_IN, '>', "$settings{tempdir}/blat-in.mfa");
	foreach my $frag_name (keys %frags) {
		my $seq = $frags{$frag_name};
		$seq =~ s/(.{80})/$1\n/g;
		$seq .= "\n" unless $seq =~ /\n$/;
		print BLAT_IN ">$frag_name\n$seq";
	}
	close BLAT_IN;

	$settings{min_blat_score} = 60;
	system("blat  -maxIntron=20  -minScore=$settings{min_blat_score} $settings{tempdir}/blat-in.mfa $settings{tempdir}/blat-in.mfa $settings{tempdir}/blat-out.psl");
	my $hits = AKUtils::readPSL2("$settings{tempdir}/blat-out.psl");
#	my $hits = AKUtils::readPSL2("tmp/blat-out.psl");

my $map;
my ($tot_hits, $tot_seqs, $tot_hit_bp);

	$settings{min_contextlen} = 90;
#	$settings{min_hitlen} = 60;
	foreach my $seq (keys %$hits) {
		foreach my $hit (@{$$hits{$seq}}) {
			next if $$hit{qName} eq $$hit{tName};



			my ($q_lflank, $q_rflank) = ($$hit{qStart}, $$hit{qSize}-$$hit{qEnd});
			my ($t_lflank, $t_rflank) = ($$hit{tStart}, $$hit{tSize}-$$hit{tEnd});
			my $l_dissim = 1 if (($$hit{strand} eq '+' and $q_lflank > $settings{min_contextlen} and $t_lflank > $settings{min_contextlen})
				or ($$hit{strand} eq '-' and $q_lflank > $settings{min_contextlen} and $t_rflank > $settings{min_contextlen}));
			my $r_dissim = 1 if (($$hit{strand} eq '+' and $q_rflank > $settings{min_contextlen} and $t_rflank > $settings{min_contextlen})
				or ($$hit{strand} eq '-' and $q_rflank > $settings{min_contextlen} and $t_lflank > $settings{min_contextlen}));

#			print "L $q_lflank:\t" if $l_dissim;
#			print "R $q_rflank:\t" if $r_dissim;
			#print join(":", values %$hit)."\n" if $l_dissim or $r_dissim;
#			print "$$hit{qName}:$$hit{qStart}..$$hit{qEnd} (".($$hit{qEnd}-$$hit{qStart}).") :: $$hit{tName}:$$hit{tStart}..$$hit{tEnd} (".($$hit{tEnd}-$$hit{tStart}).") (".sprintf("%.1f", 100*$$hit{matches}/($$hit{matches}+$$hit{misMatches}))."% id)\n"
#				if $l_dissim or $r_dissim;

			push(@{$$map{$$hit{qName}}}, $$hit{tName});

			$tot_hit_bp += $$hit{qEnd}-$$hit{qStart};
#			push(@{$$map{$$hit{tName}}}, $$hit{qName});
#			print(($$hit{matches}/($$hit{matches}+$$hit{misMatches}))."\n");
#			if ($q_lflank > $settings{min_contextlen} and $t_lflank > $settings{min_contextlen}) { print "w00t L\n"; }
#			if ($q_rflank > $settings{min_contextlen} and $t_rflank > $settings{min_contextlen}) { print "w00t R\n"; }
		}
	}

	foreach my $seq (sort keys %$map) {
#		print "$seq:\n";
		$tot_seqs++;
		foreach my $hit_seq (@{$$map{$seq}}) {
#			print "\t$hit_seq\n";
			$tot_hits++;
		}
	}
#	print "hits/seq: ".sprintf("%.2f", $tot_hits/$tot_seqs)."\n";


	open(BLAT_IN, '>', "$settings{tempdir}/mobile.blat-in.mfa");
	foreach my $seq_name (sort keys %$map) {
		foreach my $hit_seq_name (@{$$map{$seq_name}}) {
			my $seq = $frags{$hit_seq_name}; die unless $seq;
			$seq =~ s/(.{80})/$1\n/g;
			$seq .= "\n" unless $seq =~ /\n$/;
			print BLAT_IN ">$hit_seq_name\n$seq";
		}
	}
	close BLAT_IN;

	my %db_hits;
	foreach my $db_name (qw(plasmids transposons viruses)) {
		system("blat -maxIntron=20 -minScore=$settings{min_blat_score} $settings{tempdir}/mobile.blat-in.mfa $settings{datadir}/db/$db_name.fasta $settings{tempdir}/blat-out.$db_name.psl");
		$db_hits{$db_name} = AKUtils::readPSL2("$settings{tempdir}/blat-out.$db_name.psl");
	}

	print "query\tdb\thits\ttot length\tself hits\tself hit length\tsimulated frags\tfrag length\tqlength\n";
	my %hit_counts;
	my %match_counts;
	foreach my $db_name (keys %db_hits) {
		foreach my $seq_name (keys %{$db_hits{$db_name}}) {
			foreach my $hit (@{$db_hits{$db_name}->{$seq_name}}) {
				$hit_counts{$db_name}++;
				$match_counts{$db_name} += $$hit{matches};
#matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts

my $qalr = $$hit{matches}/($$hit{qEnd}-$$hit{qStart});
my $talr = $$hit{matches}/($$hit{tEnd}-$$hit{tStart});
			}
		}
		print((keys %$seqs)[0]."\t$db_name\t$hit_counts{$db_name}\t$match_counts{$db_name}\t$tot_hits\t$tot_hit_bp\t$tot_seqs\t$settings{frag_length}\t".length((values %$seqs)[0])."\n");
	}

	return 0;
}
