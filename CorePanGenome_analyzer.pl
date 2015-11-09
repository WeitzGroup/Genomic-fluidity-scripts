#!/usr/bin/env perl

# Fluidity Estimator
# Copyright (C) 2011 Georgia Institute of Technology
# Fluidity Estimator comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions.
# Contact: Joshua Weitz <jsweitz@gatech.edu>
# Contact: Andrey Kislyuk <kislyuk@gmail.com>

package CorePanGenome;

my $settings = {max_sample_size => 200};
my $stats = {};

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

use CorePanGenome;

$0 = fileparse($0);

exit(main());

sub main() {
	warn "$0: ".`uname -a`."\nCMD: $0 @ARGV\n";
	my @cmd_options = ('xcdir=s', 'min_ortho_coverage=s', 'min_ortho_identity=s', 'outfile_prefix=s', 'report_dir=s');
	GetOptions($settings, @cmd_options) or die;

	for (qw(min_ortho_coverage min_ortho_identity xcdir outfile_prefix)) {
		die("$0: argument $_ is required") unless defined $$settings{$_};
	}

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

	my @input_files = @ARGV;
	die("No input files supplied") if @input_files < 1;
	foreach my $file (@input_files) { $file = File::Spec->rel2abs($file); }
	chdir($$settings{tempdir}) or die("Unable to change working directory to $$settings{tempdir}");

	my $genes = CorePanGenome::loadGenesFromInput(\@input_files);

	logmsg "Analyzing set at i=$$settings{min_ortho_identity}, c=$$settings{min_ortho_coverage}";

	$$settings{report_dir} ||= "$$settings{xcdir}/report";
	$$settings{report_dir} .= "/$$settings{min_ortho_identity}i_$$settings{min_ortho_coverage}c";
	$$settings{manifest_file} = "$$settings{report_dir}/set.list";
	mkpath($$settings{report_dir}) or warn("Unable to create $$settings{report_dir}");
	$$settings{ortholog_pairwise_shared_file} = "$$settings{report_dir}/pair_shared.csv";
	$$settings{ortholog_pairwise_totals_file} = "$$settings{report_dir}/pair_total.csv";
	$$settings{rarefaction_stats_file} = "$$settings{report_dir}/rarefaction.csv";
	$$settings{fluidity_stats_file} = "$$settings{report_dir}/fluidity.csv";
	$$settings{pairwise_fluidity_stats_file} = "$$settings{report_dir}/pw_fluidity.csv";
	$$settings{cog_stats_file} = "$$settings{report_dir}/chogs.csv";
	$$settings{report_file} = "$$settings{report_dir}/report.txt";

	my $hits = loadOrthologHits($genes, $settings);

	my @flat_hits;
	foreach my $org1 (keys %$hits) {
		foreach my $org2 (keys %{$$hits{$org1}}) {
			foreach my $id1 (keys %{$$hits{$org1}->{$org2}}) {
				foreach my $id2 (keys %{$$hits{$org1}->{$org2}->{$id1}}) {
					my $hit = $$hits{$org1}->{$org2}->{$id1}->{$id2};
					push(@flat_hits, $hit);
				}
			}
		}
	}
	my @sorted_flat_hits = sort {$$a{evalue} <=> $$b{evalue}} @flat_hits;

	my ($cogs, $genes2cogs, $cog_matrix) = groupOrthologs($genes, $hits, \@sorted_flat_hits, undef, $settings);
	delete $$settings{report_file};
	my ($pairwise_shared_cogs, $pairwise_total_cogs) = computePairStats($genes, $hits, \@sorted_flat_hits, $settings);
=head1
	my $pairwise_shared_cogs = {}; my @orgs = sort keys %$genes;
	foreach my $i (0..$#orgs) {
		foreach my $j ($i..$#orgs) {
			my ($org1, $org2) = ($orgs[$i], $orgs[$j]);
			foreach my $c (0..$#$cogs) {
				if ($$cog_matrix[$c]->{$org1} > 0 and $$cog_matrix[$c]->{$org2} > 0) {
					$$pairwise_shared_cogs{$org1}->{$org2}++;
					$$pairwise_shared_cogs{$org2}->{$org1}++ if $org1 ne $org2; # avoid overcounting own cogs
				}
			}
		}
	}
=cut

	printOrthologStats($genes, $pairwise_shared_cogs, $pairwise_total_cogs, $settings);
	printCOGStats($genes, $cogs, $genes2cogs, $cog_matrix, $settings);

	my $coll_stats;
	$coll_stats = compileRarefaction($genes, $cogs, $genes2cogs, $cog_matrix, $settings);

	printRarefactionStats($coll_stats, $cogs, $genes, $pairwise_shared_cogs, $pairwise_total_cogs, $settings);

	open(FH, '>', $$settings{manifest_file}) or die;
	my @orgs = sort keys %$genes;
	print FH join("\n", @orgs)."\n";
	close FH;

=head1
# bugcheck
	my %org_keys;
	foreach my $i (0..$#$cogs) {
		foreach my $member (@{$$cogs[$i]}) {
			$org_keys{$member->seq_id}++;
		}
	}

	foreach my $org (keys %org_keys) {
		if ($org_keys{$org} > scalar(keys(%{$$genes{$org}}))) {
			die "Number of cogs for org $org ($org_keys{$org}) exceeds number of genes in that org (".scalar(keys(%{$$genes{$org}})).")";
		}
	}
# end bugcheck
# todo: sanity check: every gene should be a member of a cog
=cut

	return 0;
}

sub loadOrthologHits($$) {
	my ($genes, $settings) = @_;
	logmsg "Loading alignment outputs...";
	my %ortholog_hits;
	my @orgs = keys %$genes;
	foreach my $i (0..$#orgs) {
		my $org1 = $orgs[$i]; # org1 is designated as database, org2 as query
		foreach my $j ($i..$#orgs) { # Every pair, plus against self
			my $org2 = $orgs[$j]; # org1 is designated as database, org2 as query

			my ($report_file, $orgs_reversed);
			if (-f "$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out") {
				$report_file = "$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out";
			} elsif (-f "$$settings{xcdir}/$$settings{outfile_prefix}.$org1.$org2.out") {
				$report_file = "$$settings{xcdir}/$$settings{outfile_prefix}.$org1.$org2.out";
				$orgs_reversed = 1;
			} else {
				die("Input file not found: tried \"$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out\""
					. "and \"$$settings{xcdir}/$$settings{outfile_prefix}.$org2.$org1.out\"");
			}
			open(IN, '<', $report_file) or die("Unable to open file $report_file for reading: $!");
			my $l;
			while (<IN>) {
				$l++;
				chomp;
				my ($id1, $id2, $evalue, $q_coverage, $t_coverage, $percent_id) = split /\t/;
				
				next unless $q_coverage > $$settings{min_ortho_coverage}
				    and $t_coverage > $$settings{min_ortho_coverage}
				    and $percent_id >= $$settings{min_ortho_identity} * 100;

				if ($orgs_reversed) {
					die("Unable to find gene record for $org1:$id2") unless defined $$genes{$org1}->{$id2};
					die("Unable to find gene record for $org2:$id1") unless defined $$genes{$org2}->{$id1};
					$ortholog_hits{$org1}->{$org2}->{$id1}->{$id2}
						= {gene1 => $$genes{$org1}->{$id2}, gene2 => $$genes{$org2}->{$id1},
							id1 => $id2, id2 => $id1, evalue => $evalue}; #, hsp => $hsp};
				} else {
					die("Unable to find gene record for $org1:$id1") unless defined $$genes{$org1}->{$id1};
					die("Unable to find gene record for $org2:$id2") unless defined $$genes{$org2}->{$id2};
					$ortholog_hits{$org1}->{$org2}->{$id1}->{$id2}
						= {gene1 => $$genes{$org1}->{$id1}, gene2 => $$genes{$org2}->{$id2},
							id1 => $id1, id2 => $id2, evalue => $evalue}; #, hsp => $hsp};
				}
			}
			warn("No ortholog hits reported in $report_file") if $l < 1;
			close IN;
		}
	}
	return \%ortholog_hits;
}

# A gene can join a cog if it has hits above the orthology threshold to all members of the cog.
sub canJoinCog($$$;$) {
	my ($gene, $cog, $hits, $settings) = @_;
	my ($gene_id, $gene_org) = (($gene->get_tag_values('locus_tag'))[0], $gene->seq_id);

	my $ortho_count;
	foreach my $member (@$cog) {
		my ($member_id, $member_org) = (($member->get_tag_values('locus_tag'))[0], $member->seq_id);
		return 0 if $member_id eq $gene_id and $member_org eq $gene_org; # this gene is already in cog - check the other gene in the hit instead
#		print "$gene_id ($gene_org) vs. $member_id ($member_org): ".($$hits{$gene_org}->{$member_org}->{$gene_id}->{$member_id}
#			or $$hits{$member_org}->{$gene_org}->{$member_id}->{$gene_id} ? "hit found" : "no hit found")."\n";
		$ortho_count++ if $$hits{$gene_org}->{$member_org}->{$gene_id}->{$member_id}
			or $$hits{$member_org}->{$gene_org}->{$member_id}->{$gene_id};
	}
#	print "gene $gene_id org $gene_org: $ortho_count/".@$cog."\n";
	return 0 if $ortho_count < @$cog;
	return $ortho_count;
}

# Input: names of orgs and sets of hits representing orthologous relationships
# Output:
#  array of all clusters of orthologous genes (cogs): [[geneid1, geneid2, ...], [geneid3, geneid4, ...], ...]
#  hash of gene ids -> groups: { speciesid->{geneid->group, ...}, ...}
# org_mask: if defined, construct orgs from only the orgs whose ids are keys in the org_mask hash
sub groupOrthologs($$$$$) {
	my ($genes, $hits, $sorted_flat_hits, $org_mask, $settings) = @_;
	my @cogs;
	my %geneid2cogs;

	foreach my $hit (@$sorted_flat_hits) {
		die unless defined $hit;
		my ($org1, $id1) = ($hit->{gene1}->seq_id, ($hit->{gene1}->get_tag_values('locus_tag'))[0]);
		my ($org2, $id2) = ($hit->{gene2}->seq_id, ($hit->{gene2}->get_tag_values('locus_tag'))[0]);
		next if defined $org_mask and not (defined $$org_mask{$org1} and defined $$org_mask{$org2});
		print "$org1:$id1 :: $org2:$id2 $$hit{evalue}\n" if $$settings{verbose};
		if (defined $geneid2cogs{$id1} and defined $geneid2cogs{$id2}) {
			print " attempting to join 2 cogs\n" if $$settings{verbose};
			# find the bigger cog and try to join the smaller cog's contents to it one by one.
			# if successful, enroll the smaller cog's members in the bigger cog and destroy the smaller cog.
			die("Internal error") if @{$geneid2cogs{$id1}} > 1 or @{$geneid2cogs{$id2}} > 1; # multicog membership not currently supported
			next if $geneid2cogs{$id1}->[0] eq $geneid2cogs{$id2}->[0];
			my ($big_cog, $small_cog);
			if (scalar(@{$geneid2cogs{$id1}->[0]}) > scalar(@{$geneid2cogs{$id2}->[0]})) {
				($big_cog, $small_cog) = ($geneid2cogs{$id1}->[0], $geneid2cogs{$id2}->[0]);
			} else {
				($big_cog, $small_cog) = ($geneid2cogs{$id2}->[0], $geneid2cogs{$id1}->[0]);
			}
			my $cogs_joinable = 1;
			foreach my $member (@$small_cog) {
				unless (canJoinCog($member, $big_cog, $hits, $settings)) { $cogs_joinable = 0; last; }
			}
			if ($cogs_joinable) { # destroy smaller cog, enroll its members in the bigger cog
				print "  cogs of size ".@$big_cog.", ".@$small_cog." joined\n" if $$settings{verbose};
				push(@$big_cog, @$small_cog);
				foreach my $member (@$small_cog) {
					my $member_id = ($member->get_tag_values('locus_tag'))[0];
					$geneid2cogs{$member_id}->[0] = $big_cog;
				}
			}
		} elsif (defined $geneid2cogs{$id1}) {
			print " attempting to enroll $id2 in cog of $id1\n" if $$settings{verbose};
			# try to join id2 to id1's cog. if can't, make a singleton cog. destroy all other hits encountered
			die("Internal error") if @{$geneid2cogs{$id1}} > 1; # multicog membership not currently supported
			my $cog = $geneid2cogs{$id1}->[0];
			if (canJoinCog($hit->{gene2}, $cog, $hits, $settings)) {
				push(@$cog, $hit->{gene2});
				$geneid2cogs{$id2}->[0] = $cog;
			} else { # new singleton cog
				$geneid2cogs{$id2}->[0] = [$hit->{gene2}];
			}
		} elsif (defined $geneid2cogs{$id2}) {
			print " attempting to enroll $id1 in cog of $id2\n" if $$settings{verbose};
			# try to join id1 to id2's cog. if can't, make a singleton cog. destroy all other hits encountered
			die("Internal error") if @{$geneid2cogs{$id2}} > 1; # multicog membership not currently supported
			my $cog = $geneid2cogs{$id2}->[0];
			if (canJoinCog($hit->{gene1}, $cog, $hits, $settings)) {
				push(@$cog, $hit->{gene1});
				$geneid2cogs{$id1}->[0] = $cog;
			} else { # new singleton cog
				$geneid2cogs{$id1}->[0] = [$hit->{gene1}];
			}
		} else {
			print " form new cog\n" if $$settings{verbose};
			# new cog
			my $cog = [$hit->{gene1}, $hit->{gene2}];
			$geneid2cogs{$id1}->[0] = $cog;
			$geneid2cogs{$id2}->[0] = $cog;
		}
	}

	# Count true singletons (genes not seen in any alignment)
	foreach my $org (keys %$genes) {
		my ($s, $m);
		foreach my $gene_id (keys %{$$genes{$org}}) {
			next if defined $geneid2cogs{$gene_id};
			my $cog = [$$genes{$org}->{$gene_id}];
			$geneid2cogs{$gene_id}->[0] = $cog;
#			if (not defined $geneid2cogs{$gene}) { $m++; } else { $s++; }
		}
		
	}

	my %c;
	foreach my $id (keys %geneid2cogs) {
		$c{$geneid2cogs{$id}->[0]} = $geneid2cogs{$id}->[0];
	}
	@cogs = values %c;
	# TODO: form @cogs from geneid2cogs here

	my @cog_matrix; # $cog_matrix[$cog_i]->{org_id} = N (number of participant genes)
	my %member_matrix;
	my %org_keys;
	foreach my $i (0..$#cogs) {
		foreach my $member (@{$cogs[$i]}) {
			$member_matrix{$member}++;
			my $id = ($member->get_tag_values('locus_tag'))[0];
			die("gene $id is present in $member_matrix{$member} cogs") if $member_matrix{$member} > 1;
			$cog_matrix[$i]->{$member->seq_id}++;
			$org_keys{$member->seq_id}++;
		}
	}

	if (defined $$settings{report_file}) {
		open(REPORT_FH, '>', $$settings{report_file}) or die;
		print REPORT_FH "Genes per org:\n"; print REPORT_FH "\t$_\t".keys(%{$$genes{$_}})."\n" for keys %$genes;
		my $c; $c += scalar(keys(%{$$genes{$_}})) for keys %$genes;
		print REPORT_FH "$c genes total\n";
		print REPORT_FH "CHoGs per org:\n"; print REPORT_FH "\t$_\t$org_keys{$_}\n" for keys %org_keys;
		print REPORT_FH @cogs." cogs total\n";
		my %cshist; for my $c (@cogs) { $cshist{scalar(@$c)}++; }
		print REPORT_FH "Cluster sizes\n";
		print REPORT_FH "\t$_\t$cshist{$_}\t".int(100*$cshist{$_}/scalar(@cogs))."%\n" for sort {$a<=>$b} keys %cshist;
=head1
	foreach my $org (keys %$genes) {
		my ($s, $m);
		foreach my $gene (keys %{$$genes{$org}}) {
			if ($geneid2cogs{$gene}) { $m++; } else { $s++; }
		}
		print REPORT_FH "Org $org, $s singleton genes, $m cog participant genes\n";
	}
=cut
		close REPORT_FH;
	}

#  fixme: sanity check: number of cogs (membership) per org can never exceed number of genes

	return (\@cogs, \%geneid2cogs, \@cog_matrix);
}

# Re-compute COG tables for just 2 genomes at a time, to avoid break-up interference from other genomes in the set.
sub computePairStats($$$$) {
	my ($genes, $hits, $sorted_flat_hits, $settings) = @_;
	my @orgs = sort keys %$genes;

	my (%pairwise_shared_cogs, %pairwise_total_cogs);
	foreach my $i (0..$#orgs) {
		foreach my $j ($i..$#orgs) {
			my ($org1, $org2) = ($orgs[$i], $orgs[$j]);
			my ($cogs, $genes2cogs, $cog_matrix) = groupOrthologs($genes, $hits, $sorted_flat_hits, {$org1 => 1, $org2 => 1}, $settings);
			foreach my $c (0..$#$cogs) {
				if ($$cog_matrix[$c]->{$org1} > 0 and $$cog_matrix[$c]->{$org2} > 0) {
					$pairwise_shared_cogs{$org1}->{$org2}++;
					$pairwise_shared_cogs{$org2}->{$org1}++ if $org1 ne $org2; # avoid overcounting own cogs
				}
				if ($$cog_matrix[$c]->{$org1} > 0 or $$cog_matrix[$c]->{$org2} > 0) {
					$pairwise_total_cogs{$org1}->{$org2}++;
					$pairwise_total_cogs{$org2}->{$org1}++ if $org1 ne $org2; # avoid overcounting own cogs
				}
			}
		}
	}
	return (\%pairwise_shared_cogs, \%pairwise_total_cogs);
}

# In a collector's run, add a new species' genes to the collection.
sub addSpeciesToSet($$$$) {
#	my ($genes, $genes2cogs, $permutation, $timestep, $org_id, $collector) = @_;
	my ($permutation, $timestep, $cog_matrix, $collector) = @_;

	my (@cogs_seen_in_orgs, @cogs_members_in_orgs);
	foreach my $i (0..$#$cog_matrix) {
		foreach my $collected_org (@$permutation[0..$timestep]) {
			$cogs_seen_in_orgs[$i]++ if $$cog_matrix[$i]->{$collected_org};
			$cogs_members_in_orgs[$i] += $$cog_matrix[$i]->{$collected_org};
		}
		next unless $cogs_seen_in_orgs[$i];
		$$collector{$timestep}->{$cogs_seen_in_orgs[$i]}->{presence}++;
#		$$collector{$timestep}->{$cogs_seen_in_orgs[$i]}->{total} += $cogs_members_in_orgs[$i]; #CHECK
	}

=head1
	foreach my $gene_id (keys %{$$genes{$org_id}}) {
		# "This gene has been seen in X previously examined orgs"
		my ($seen_in_org_count, $total_in_org_count) = (0, 0);
		if ($$genes2cogs{$gene_id}) {
			my (%seen_in_orgs, %total_in_orgs);
			foreach my $cog (@{$$genes2cogs{$gene_id}}) {
				foreach my $member (@$cog) {
					$seen_in_orgs{$member->seq_id} = 1;
					$total_in_orgs{$member->seq_id}++;
				}
			}
			foreach my $prev_org (@$permutation[0..$timestep-1]) {
				$seen_in_org_count++ if $seen_in_orgs{$prev_org};
				$total_in_org_count += $total_in_orgs{$prev_org};
			}
		} else { # this gene is a singleton
			$seen_in_org_count = 0;
		}
		# now update with counts from current org
		# FIXME: These counts reflect only genes seen in the current org, but the histogram must include genes seen in other orgs
		$$collector{$timestep}->{$seen_in_org_count}->{presence}++;
#		$$collector{$timestep}->{$seen_in_org_count}->{total} += $total_in_org_count;
	}
=cut
#	for each gene in org, check all cogs that gene is in.
#	count # orgs which genes from those cogs originate from.
#	[ updates for previously accounted cogs -> ? ]
#	return ($set, $num_new_genes, $num_core_genes);
	return $collector;
}

# Given BLAST hits, output numbers of pairwise hits for every pair of orgs
sub printOrthologStats($$$$) {
	my ($genes, $pairwise_shared_cogs, $pairwise_total_cogs, $settings) = @_;
	my @orgs = sort keys %$genes;
	open(FH, '>', $$settings{ortholog_pairwise_shared_file}) or die;
	foreach my $org1 (@orgs) {
		my @l;
		foreach my $org2 (@orgs) {
			my $value = $$pairwise_shared_cogs{$org1}->{$org2}; $value ||= 0;
			push(@l, $value);
		}
		print FH join(",", @l)."\n";
	}
	close FH;
	open(FH, '>', $$settings{ortholog_pairwise_totals_file}) or die;
	foreach my $org1 (@orgs) {
		my @l;
		foreach my $org2 (@orgs) {
			my $value = $$pairwise_total_cogs{$org1}->{$org2}; $value ||= 0;
			push(@l, $value);
		}
		print FH join(",", @l)."\n";
	}
	close FH;
}

# Print cog_matrix to a file.
sub printCOGStats($$$$$) {
	my ($genes, $cogs, $genes2cogs, $cog_matrix, $settings) = @_;

	my @orgs = sort keys %$genes;

	open(FH, '>', $$settings{cog_stats_file}) or die;
	print FH join(",", ("cog", @orgs))."\n";
	foreach my $i (0..$#$cogs) {
		my @line = ($i);
		foreach my $org (@orgs) {
			# $cog_matrix[$cog_i]->{org_id} = N (number of participant genes)
			my $value = $$cog_matrix[$i]->{$org}; $value ||= 0;
			push(@line, $value);
		}
		print FH join(",", @line)."\n";
	}
	close FH;
}

sub getPermutations($$) {
	my ($genes, $settings) = @_;
	my @species_permutations;
	my $max_exhaustive_set_size = 1;
	my $k = 1;
	while ($k < $$settings{max_sample_size}) { $max_exhaustive_set_size++; $k *= $max_exhaustive_set_size; }

	if (scalar(keys %$genes) < $max_exhaustive_set_size) { # exhaustive search
		@species_permutations = AKUtils::permutations(keys %$genes);
	} else {
		my @set = keys %$genes;
		my $set_size = scalar(@set);
		for (1..$$settings{max_sample_size}) {
			my $p = AKUtils::sampleWithoutReplacement(\@set, scalar(@set));
			push(@species_permutations, $p);
		}
	}
	return \@species_permutations;
}

# Given BLAST hits, run collector's curves for all possible species combinations
# Find distributions (min, max, avg, ...) of all cluster size deltas (0=>new, N=>core, ...) for each increment position
# Fit a rarefaction curve with estimator of choice
sub compileRarefaction($$$$$) {
	my ($genes, $cogs, $genes2cogs, $cog_matrix, $settings) = @_;

	# collect_series = {permutation => {timestep => {cog size => {(total|presence) => count}}}}
	# Depending on definition, core genes count is average of presence counts for cog size = # species over all permutations,
	my %collect_series;
	# collect_counts = {timestep => {cog size => [count, count, ...]}} (presence count in all permutations)
	my %collect_counts;

	logmsg "Compiling rarefaction data...";
	my $species_permutations = getPermutations($genes, $settings);

	foreach my $i (0..$#$species_permutations) {
		my $permutation = $$species_permutations[$i];
		my %collector;

		warn "[$i/".@$species_permutations."] P: [@$permutation]\n";
		foreach my $timestep (0..$#$permutation) {
			my $org_id = $$permutation[$timestep];

#			addSpeciesToSet($genes, $genes2cogs, $permutation, $timestep, $org_id, \%collector);
			addSpeciesToSet($permutation, $timestep, $cog_matrix, \%collector);

##			print "\t$timestep/[@$permutation]\n";
			foreach my $bin (sort keys %{$collector{$timestep}}) {
##				print "\t\t$bin\t$collector{$timestep}->{$bin}->{presence}\n";
				push(@{$collect_counts{$timestep}->{$bin}}, $collector{$timestep}->{$bin}->{presence});
			}
#			foreach my $gene_id (keys %{$$genes{$org_id}}) {
#				foreach my $cog (@{$$genes2cogs{$gene_id}}) {
#					print "P: [@$permutation]\t$org_id\t$gene_id\t[@$cog]\n";
#				}
#			}
		}
#		$collect_series{$permutation} = \%collector;
	}

	# Compute means, stdevs, min, med, max, 5%, 95% bounds
	my %coll_stats;
	foreach my $timestep (sort {$a <=> $b} keys %collect_counts) {
		foreach my $bin (sort {$a <=> $b} keys %{$collect_counts{$timestep}}) {
			die("Internal error: missing collector data") if @{$collect_counts{$timestep}->{$bin}} < 1;
#			my $stat = Statistics::Descriptive::Full->new();
#			$stat->add_data(@{$collect_counts{$timestep}->{$bin}});
#			$coll_stats{$timestep}->{$bin} = $stat;
			my @values = @{$collect_counts{$timestep}->{$bin}};
			$coll_stats{$timestep}->{$bin}->{min} = min(@values);
			$coll_stats{$timestep}->{$bin}->{max} = max(@values);
			$coll_stats{$timestep}->{$bin}->{avg} = sprintf("%.8f", sum(@values)/scalar(@values));
			$coll_stats{$timestep}->{$bin}->{stdev} = sprintf("%.8f", AKUtils::stdev(@values));
			$coll_stats{$timestep}->{$bin}->{median} = AKUtils::median(@values);
#			$coll_stats{$timestep}->{$bin}->{p5} = 
#			$coll_stats{$timestep}->{$bin}->{p95} = 
		}
	}
	return \%coll_stats;
}

sub printRarefactionStats($$$$$$) {
	my ($coll_stats, $cogs, $genes, $pairwise_shared_cogs, $pairwise_total_cogs, $settings) = @_;

	if (defined $coll_stats) {
		open(REPORT, '>', $$settings{rarefaction_stats_file}) or die;
		foreach my $timestep (sort {$a <=> $b} keys %$coll_stats) {
			foreach my $bin (sort {$a <=> $b} keys %{$$coll_stats{$timestep}}) {
	##			print REPORT "TS $timestep\tBin $bin: @{$collect_counts{$timestep}->{$bin}}\n\tm=", $coll_stats{$timestep}->{$bin}->{avg},
	##				"\tsd=", $coll_stats{$timestep}->{$bin}->{stdev},
	#				"\t[5,95]=", $stat->percentile(5), ",", $stat->percentile(95),
	##				"\n";
				my @l = ($timestep, $bin);
				for (qw(avg stdev min max median)) {
					push(@l, $$coll_stats{$timestep}->{$bin}->{$_});
				}
				print REPORT join(",", @l)."\n";
			}
		}
		close REPORT;
		
	#	my $species_permutations = getPermutations($genes, $settings);

		open(REPORT, '>', $$settings{fluidity_stats_file}) or die;
		# S_Chao1 = S_obs + (n1 (n1-1)) / (2(n2 + 1))
		# S_obs = observed #genes ("species") -- WRONG
		# S_obs = observed #OTUs (CHoGs)
		# n1 = number of OTUs with 1 sequence [cogs of size 1 or represented in only 1 sp.], n2 = same w/2 sequences
		my $s_obs = scalar(@$cogs);
	#	$s_obs += scalar(keys %{$$genes{$_}}) for keys %$genes;
		my $n1 = $$coll_stats{scalar(keys %$genes)-1}->{1}->{avg};
		my $n2 = $$coll_stats{scalar(keys %$genes)-1}->{2}->{avg};
		my $s_chao1 = $s_obs + (($n1 * ($n1 - 1)) / (2*($n2 + 1)));
		print REPORT "S_CHAO1 = $s_chao1 [$s_obs $n1 $n2]\n";

		# Fluidity:
		# nu_hat = (1/K) \sum_{m=2}^{M}{(m(m-1)/(M(M-1))) G_m}
		# K = average number of genes/genome = avg_genes
		# M = total number of genomes = tot_genomes
		# G_m = number of CHoGs present in m of M genomes
		my $tot_genomes = scalar(keys %$genes); # TODO: count cogs, not genes
		my ($tot_genes, $avg_genes);
		foreach my $org (keys %$genes) {
			$tot_genes += scalar(keys %{$$genes{$org}});
		}
		$avg_genes = $tot_genes / $tot_genomes;

		# TODO: this sampling method is incorrect
		#foreach my $tot_genomes (2..scalar(keys %$genes)) {
		foreach my $tot_genomes (scalar(keys %$genes)) {
			my $fluidity;

			foreach my $m (2..$tot_genomes) {
				# die if $$coll_stats{$tot_genomes-1}->{$m}->{stdev} != 0;
				my $G_m = $$coll_stats{$tot_genomes-1}->{$m}->{avg}; # / scalar(@$cogs);
				my $weight = ($m * ($m - 1)) / ($tot_genomes * ($tot_genomes - 1));
				$fluidity += ($weight * $G_m);
				# print REPORT "tot_genomes=$tot_genomes; tot_genes=$tot_genes; avg_genes=$avg_genes; m=$m; G_m=$G_m; weight=$weight; running_f=$fluidity\n";
			}
			$fluidity = 1 - ($fluidity / $avg_genes);
			print REPORT "$tot_genomes,$fluidity\n";
		}
		close REPORT;
	}

	# Pairwise fluidity:
	# nu = 1 - avg(sum_{i,j}{ 2 S_{ij} / (M_i + M_j) })

	open(REPORT, '>', $$settings{pairwise_fluidity_stats_file}) or die;
	my $species_permutations = getPermutations($genes, $settings);

	my %pairwise_fluidities; # fluidities{2} = [x, y, z]
	foreach my $permutation (@$species_permutations) {
		foreach my $subgroup_size (2..@$permutation) {
			my $pairwise_fluidity;
			my $k;
			foreach my $i (0..$subgroup_size-1) {
				foreach my $j ($i+1..$subgroup_size-1) {
					$k++;
					my ($org1, $org2) = ($$permutation[$i], $$permutation[$j]);
					my $shared_cogs = $$pairwise_shared_cogs{$org1}->{$org2};
					my $org1_cogs = $$pairwise_total_cogs{$org1}->{$org1};
					my $org2_cogs = $$pairwise_total_cogs{$org2}->{$org2};
					$pairwise_fluidity += ((2 * $shared_cogs) / ($org1_cogs + $org2_cogs));
#					my $c_pwf = 1 - ($pairwise_fluidity / $k);
#					print REPORT "shared cogs=$shared_cogs; org1 cogs=$org1_cogs; org2 cogs=$org2_cogs; unc_pwf=$pairwise_fluidity, pwf=$c_pwf\n";
				}
			}
			$pairwise_fluidity = 1 - ($pairwise_fluidity / $k);
			push(@{$pairwise_fluidities{$subgroup_size}}, sprintf("%.8f", $pairwise_fluidity));
		}
	}
	foreach my $subgroup_size (sort {$a<=>$b} keys %pairwise_fluidities) {
		my @values = @{$pairwise_fluidities{$subgroup_size}};
		my %stats;
		$stats{avg} = sprintf("%.8f", sum(@values)/scalar(@values));
		$stats{stdev} = sprintf("%.8f", AKUtils::stdev(@values));
		$stats{min} = min(@values);
		$stats{max} = max(@values);
		$stats{median} = AKUtils::median(@values);
		print REPORT join(',', $subgroup_size, $stats{avg}, $stats{stdev}, $stats{min}, $stats{max}, $stats{median})."\n";
	}
	close REPORT;
}
