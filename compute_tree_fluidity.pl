#!/usr/bin/env perl

# Fluidity Estimator
# Copyright (C) 2011 Georgia Institute of Technology
# Fluidity Estimator comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions.
# Contact: Joshua Weitz <jsweitz@gatech.edu>
# Contact: Andrey Kislyuk <kislyuk@gmail.com>

my $settings = {
	min_ortho_coverage => 0.5,
	min_ortho_identity => 0.5,
	
};
my $stats;

use strict;
use File::Basename;
use File::Path;
use Bio::SeqIO;
use Bio::TreeIO;
use FindBin;
use lib "$FindBin::RealBin/lib";

$0 = fileparse($0);

exit(main());

sub main() {
	die("Usage: $0 tree.newick data_dir xcdir\n") if @ARGV < 3;

	my ($tree_file, $data_dir, $xcdir) = @ARGV;
	$$settings{xcdir} = $xcdir;
	$$settings{xcdir} = File::Spec->rel2abs($$settings{xcdir});
	$$settings{data_dir} = File::Spec->rel2abs($$settings{data_dir});
	$$settings{outfile_prefix} = "CorePanGenome.pl";
	$$settings{report_dir} = $xcdir."/report";
	die("Tree file $tree_file does not exist") unless -f $tree_file;
	my $tree = Bio::TreeIO->new(-format => 'newick', -file => $tree_file)->next_tree;
	die("Unable to parse a tree from file $tree_file") unless $tree;

	my %node_manifests;
	foreach my $node ($tree->get_nodes(-order => 'breadth')) {
		next if $node->is_Leaf();
		foreach my $child ($node->get_all_Descendents) {
			next unless $child->is_Leaf();
			my $child_filename = "$data_dir/".$child->id.".gbk.fasta.reannot.gbk";
			my $child_filename = "$data_dir/".$child->id.".gbk";
			die("file $child_filename not found") unless -f $child_filename;
			$node_manifests{$node}->{$child} = $child_filename;
			# warn "Node $node manifest: $child_filename\n";
		}
		next if scalar(keys %{$node_manifests{$node}}) < 2;
		
		my @files = values %{$node_manifests{$node}};
		my $invoke_string = "$FindBin::RealBin/CorePanGenome_analyzer.pl @files";
		for (qw(min_ortho_coverage min_ortho_identity xcdir outfile_prefix)) { # report_dir)) {
			$invoke_string .= " --$_=$$settings{$_}";
		}
		
#		print "$invoke_string\n";
		system("$invoke_string");
		warn "MANIFEST: @files\n";
		system("tail -n 1 $$settings{report_dir}/$$settings{min_ortho_identity}i_$$settings{min_ortho_coverage}c/pw_fluidity.csv");
		
		rmtree($$settings{report_dir});
	}
}
