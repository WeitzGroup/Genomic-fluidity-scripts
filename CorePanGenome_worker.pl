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
#	tempdir => 'cptmp', keep => 1,
	min_ortho_coverage => 0.5,
	min_ortho_identity => 0.3,
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

$0 = fileparse($0);

exit(main());

sub main() {
	warn "$0: ".`uname -a`."\nCMD: $0 @ARGV\n";
	my @cmd_options = ('outfile=s', 'run_id=s', 'query_mfa=s', 'db_mfa=s', 'min_ortho_coverage=s', 'min_ortho_identity=s');
	GetOptions($settings, @cmd_options) or die;
	die("Arguments outfile and run_id must be supplied") unless $$settings{outfile} and $$settings{run_id};
	die("Arguments query_mfa and db_mfa must be supplied") unless $$settings{query_mfa} and $$settings{db_mfa};
	die("File $$settings{query_mfa} does not exist") unless -f $$settings{query_mfa};
	die("File $$settings{db_mfa} does not exist") unless -f $$settings{db_mfa};

	my $query_seqs = AKUtils::readMfa($$settings{query_mfa});
	my $db_seqs = AKUtils::readMfa($$settings{db_mfa});

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

	my $db_loc = AKUtils::formatBLASTdb($$settings{db_mfa}, {%$settings, formatdb_protein_db => 1});
	$ENV{BLASTDB} = (fileparse($db_loc))[1];
	my $bf = Bio::Tools::Run::StandAloneBlast->new(database => $db_loc,
		outfile => "$$settings{tempdir}/$0.$$.blast_out", # $settings{outfile},
		program => 'blastp',
		a => AKUtils::getNumCPUs());

	logmsg "Running BLAST on $$settings{query_mfa} vs. $$settings{db_mfa} [$db_loc]...";
	my $report = $bf->blastall($$settings{query_mfa});
#	my $report = new Bio::SearchIO(-file => $report_file); #, -format => 'blasttable');
	my %reported_hits;

	open(OUT, '>', $$settings{outfile}) or die;
	while (my $result = $report->next_result) {
		while (my $hit = $result->next_hit) {
			HSP: while (my $hsp = $hit->next_hsp) {
				my ($id1, $id2) = ($hsp->seq('hit')->id, $hsp->seq('query')->id);
				next if $id1 eq $id2; # discard self hits. TODO: if same ids occur in diff orgs, we have a problem

				# WARNING: the line below skips pairs already seen, assuming reciprocity
				next if $reported_hits{$id1}->{$id2};

#foreach my $gene (values %{$$genes{$org1}}) {
#$blast_db_h->write_seq(Bio::Seq->new(-seq => $gene->seq->translate()->seq(), -display_id => ($gene->get_tag_values('locus_tag'))[0]));
# WARNING: WRONG
#				my $q_coverage = $hsp->length('query') / length($hsp->query_string);
#				my $t_coverage = $hsp->length('hit') / length($hsp->hit_string);
# END WARNING
				my $l1 = length($$db_seqs{$id1}); die("Internal error: sequence $id1 not found") unless $l1;
				my $l2 = length($$query_seqs{$id2}); die("Internal error: sequence $id2 not found") unless $l2;
				my $q_coverage = $hsp->length('query') / $l2;
				my $t_coverage = $hsp->length('hit') / $l1;

#				warn("$id1, $id2: $q_coverage = ".$hsp->length('query')." / $l2; $t_coverage = ".$hsp->length('hit')." / $l1\nD:$$db_seqs{$id1}\nQ:$$query_seqs{$id2}\nQS:".$hsp->query_string."\nDS:".$hsp->hit_string."\n"); #[L1Q $l1q L1D $l1d L2Q $l2q L2D $l2d] ".$hsp->length('query')." / ".length($hsp->query_string)." = $q_coverage; ".$hsp->length('hit')." / ".length($hsp->hit_string)." = $t_coverage\n".$hsp->query_string."\n".$hsp->hit_string."\n");
				die("Internal error") if $q_coverage > 1 or $t_coverage > 1;

				next if $hsp->percent_identity < $$settings{min_ortho_identity} * 100;
				

				if ($q_coverage > $$settings{min_ortho_coverage}
					and $t_coverage > $$settings{min_ortho_coverage}
					and $hsp->percent_identity >= $$settings{min_ortho_identity} * 100) {
=head1
					# Retain only the best hit for each gene1
					foreach my $alt_hit_id2 (keys %{$ortholog_hits{$org1}->{$org2}->{$id1}}) {
						my $alt_hsp = $ortholog_hits{$org1}->{$org2}->{$id1}->{$alt_hit_id2}->{hsp};
						if ($alt_hsp->evalue < $hsp->evalue) { # better hit already found
							next HSP;
						} else { # discard previous hit
							delete $ortholog_hits{$org1}->{$org2}->{$id1}->{$alt_hit_id2};
						}
					}
=cut
#							warn if $ortholog_hits{$org1}->{$org2}->{$id1}->{$id2}; # TODO: skip if evalue > cur
							# In the HSP, gene1 is "hit", gene2 is "query"
##							$ortholog_hits{$org1}->{$org2}->{$id1}->{$id2}
##								= {gene1 => $$genes{$org1}->{$id1}, gene2 => $$genes{$org2}->{$id2}, hsp => $hsp};
#							print OUT "$org1\t$id1\t$org2\t$id2\n";
					# TODO: may need to report hit rank, other details
					# WARNING: the line below marks pairs already seen, assuming reciprocity
                    $reported_hits{$id1}->{$id2} = 1;

                    # BUGCHECK: this is egregiously wrong
					# $reported_hits{$id2}->{$id1} = 1;
                    my @l = ($id1, $id2, $hsp->evalue, sprintf("%.4f", $q_coverage), sprintf("%.4f", $t_coverage), sprintf("%.4f", $hsp->percent_identity));
                    print OUT join("\t", @l)."\n";
				}
			}
		}
	}
	close OUT;
#	my $report = $bf->blastall($$settings{query_mfa});
	system("touch $$settings{outfile}.done");
	logmsg "Report is in $$settings{outfile}";
	return 0;
}
