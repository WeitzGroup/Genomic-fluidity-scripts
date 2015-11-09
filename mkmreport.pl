#!/usr/bin/env perl

# Fluidity Estimator
# Copyright (C) 2011 Georgia Institute of Technology
# Fluidity Estimator comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions.
# Contact: Joshua Weitz <jsweitz@gatech.edu>
# Contact: Andrey Kislyuk <kislyuk@gmail.com>


use strict;
use List::Util qw(min max sum shuffle);

my $workdir = $ARGV[0]; $workdir ||= ".";

my $r_prefix="report";

open(R, '>', "$workdir/$r_prefix.m") or die;

my (%cogs, %shared, %pw_shared_cogs, %pw_total_cogs, %raref, %gpo, %cpo, %cs, %pw_fluidity);

foreach my $dir (glob "$workdir/*i_*c") {
	open(IN, '<', "$dir/report.txt") or die;
	my $rtxt; { local $/; $rtxt = <IN>; }
	close IN;
	$rtxt =~ /Genes per org:\n(.+)\n(\d+) genes total\nCOGs per org:\n(.+)\n(\d+) cogs total\nCluster sizes\n(.+)/s;
	my ($gpo, $gt, $cpo, $ct, $cs) = ($1, $2, $3, $4, $5);
	for (split /\n/, $gpo) { /\s(.+)\s(\d+)/; $gpo{$1}=$2; }
	for (split /\n/, $cpo) { /\s(.+)\s(\d+)/; $cpo{$1}=$2; }
	for (split /\n/, $cs) { /\s([^\s]+)\s(\d+)/; $cs{$1}=$2; }

	open(IN, '<', "$dir/chogs.csv") or die;
	my $ch=<IN>; chomp $ch; my @ch=split(/,/, $ch); shift @ch;
	my @cogs;
	while (<IN>) { chomp; my @l=split(/,/); shift @l; push(@cogs, \@l); }
	close IN;

=head1
	my @pairwise_total_cogs;
	foreach my $cog (@cogs) {
		foreach my $i (0..scalar(keys %gpo)-1) {
			foreach my $j (0..scalar(keys %gpo)-1) {
				$pairwise_total_cogs[$i]->[$j]++ if $$cog[$i] > 0 or $$cog[$j] > 0;
			}
		}
	}
=cut

	my (@pairwise_shared_cogs, @pairwise_total_cogs);
	open(IN, '<', "$dir/pair_shared.csv") or die("Could not open $dir/pair_shared.csv for reading: $!");
	while (<IN>) { chomp; my @l=split(/,/); push(@pairwise_shared_cogs, \@l); }
	close IN;
	open(IN, '<', "$dir/pair_total.csv") or die("Could not open $dir/pair_total.csv for reading: $!");
	while (<IN>) { chomp; my @l=split(/,/); push(@pairwise_total_cogs, \@l); }
	close IN;

	my %trial_raref;
	open(IN, '<', "$dir/rarefaction.csv") or warn("Could not open $dir/rarefaction.csv for reading: $!");
	while (<IN>) {
		chomp;
		my ($timestep, $bin, $mean, $stdev, $min, $max, $median) = split(/,/);
		$trial_raref{$timestep}->{$bin} = {mean => $mean, stdev => $stdev, min => $min, max => $max, median => $median};
	}
	close IN;

	open(IN, '<', "$dir/pw_fluidity.csv") or die;
	my %trial_pwf;
	while (<IN>) {
		chomp;
		my ($bin, $mean, $stdev, $min, $max, $median) = split(/,/);
		$trial_pwf{$bin} = {mean => $mean, stdev => $stdev, min => $min, max => $max, median => $median};
	}
	close IN;

	$dir =~ /(\d.*)i_(.+)c/; my ($i, $c) = ($1, $2);
	$cogs{$i}->{$c} = \@cogs;
	$pw_shared_cogs{$i}->{$c} = \@pairwise_shared_cogs;
	$raref{$i}->{$c} = \%trial_raref;
	$pw_total_cogs{$i}->{$c} = \@pairwise_total_cogs;
	$pw_fluidity{$i}->{$c} = \%trial_pwf;
}


print R "cg_set.name='$workdir';\n";
print R "cg_set.genome_ids={";
for (sort keys %gpo) { print R "'$_'," }
print R "};\n";
print R "cg_set.gene_counts=[";
for (sort keys %gpo) { print R "$gpo{$_} " }
print R "];\n";

print R "cg_set.data = {};\n";
foreach my $i (sort keys %cogs) {
	foreach my $c (sort keys %{$cogs{$i}}) {
		print R "cg_set.data{length(cg_set.data)+1} = struct('min_aln_identity', $i, 'min_aln_coverage', $c, ...\n";
		print R "\t'cog_counts', [";
		for (sort keys %cpo) { print R "$cpo{$_}," }
		print R "\t], ...\n";
=head1
		print R "\t'rarefaction', [ ...\n";
		foreach my $ts (sort {$a <=> $b} keys %{$raref{$i}->{$c}}) {
			my @bin_lines;
			foreach my $bin (sort {$a <=> $b} keys %{$raref{$i}->{$c}->{$ts}}) {
				my @field_lines = ("'timestep'", $ts, "'bin'", $bin);
				foreach my $field (qw(mean stdev min max median)) {
					push(@field_lines, "'$field'", "$raref{$i}->{$c}->{$ts}->{$bin}->{$field}");
				}
				push(@bin_lines, "struct(".join(", ", @field_lines).")");
			}
			print R "\t\t[".join(",\n\t\t\t", @bin_lines)."]; ...\n";
		}
		print R "\t], ...\n";
=cut

		if (defined $raref{$i}->{$c}) {
		my $max_ts = max(keys %{$raref{$i}->{$c}});
		print R "\t'rarefaction', struct( ...\n\t\t";
		my @raref_lines;
		foreach my $field (qw(mean stdev min max median)) {
			push(@raref_lines, "'$field', [");
			foreach my $ts (sort {$a <=> $b} keys %{$raref{$i}->{$c}}) {
				my @l;
#				foreach my $bin (sort {$a <=> $b} keys %{$raref{$i}->{$c}->{$ts}}) {
				foreach my $bin (1..$max_ts+1) {
					my $value;
					if (defined $raref{$i}->{$c}->{$ts}->{$bin}) {
						$value = $raref{$i}->{$c}->{$ts}->{$bin}->{$field};
					} else {
						$value = -1;
					}
					push(@l, $value);
				}
				push(@raref_lines, join(" ", @l).";");
			}
			push(@raref_lines, "]");
		}
		print R join(", ...\n\t\t", @raref_lines)."), ...\n";
		}

		print R "\t'pw_fluidity', struct( ...\n\t\t";
		my @pwf_lines;
		foreach my $field (qw(mean stdev min max median)) {
			push(@pwf_lines, "'$field', [");
			my @l;
			foreach my $bin (1..max(keys %{$pw_fluidity{$i}->{$c}})) {
				push(@l, $pw_fluidity{$i}->{$c}->{$bin}->{$field});
			}
			push(@pwf_lines, join(" ", @l).";");
			push(@pwf_lines, "]");
		}
		print R join(", ...\n\t\t", @pwf_lines)."), ...\n";

		print R "\t'pairwise_shared_cogs', [ ...\n";
		foreach my $row (@{$pw_shared_cogs{$i}->{$c}}) {
			print R "\t\t@$row; ...\n";
		}
		print R "\t], ...\n";
		print R "\t'pairwise_total_cogs', [ ...\n";
		foreach my $row (@{$pw_total_cogs{$i}->{$c}}) {
			print R "\t\t@$row; ...\n";
		}

		print R "\t], ...\n";

		print R "\t'cogs', [ ...\n";

		foreach my $row (@{$cogs{$i}->{$c}}) {
			print R "\t\t@$row; ...\n";
		}

		print R "]);\n";
	}
}

close R;

# system("cd $workdir; octave --eval \"$r_prefix; save('$r_prefix.mat', 'cg_set');\" > /dev/null");

system("cd $workdir; matlab -r \"$r_prefix; save('$r_prefix.mat', 'cg_set'); exit;\" > /dev/null");
die if $?;
