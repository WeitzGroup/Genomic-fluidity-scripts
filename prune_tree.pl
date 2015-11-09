#!/usr/bin/env perl

# Fluidity Estimator
# Copyright (C) 2011 Georgia Institute of Technology
# Fluidity Estimator comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions.
# Contact: Joshua Weitz <jsweitz@gatech.edu>
# Contact: Andrey Kislyuk <kislyuk@gmail.com>

my $settings = {
};
my $stats;

use strict;
use File::Basename;
use Bio::SeqIO;
use Bio::TreeIO;

$0 = fileparse($0);

exit(main());

sub main() {
	die("Usage: $0 tree.newick node1 node2 node3 ...\n") if @ARGV < 2;

	my ($tree_file, @nodes) = @ARGV;
	die("Tree file $tree_file does not exist") unless -f $tree_file;
	my $tree = Bio::TreeIO->new(-format => 'newick', -file => $tree_file)->next_tree;
	die("Unable to parse a tree from file $tree_file") unless $tree;

	my %nodes_to_keep;
	$nodes_to_keep{$_} = 1 for @nodes;

	pruneTree($tree, \%nodes_to_keep);
	
	my $out_tree_io = new Bio::TreeIO(-file => ">$tree_file.pruned.tree", -format => 'newick');
	$out_tree_io->write_tree($tree);

	my $out_tree_io = new Bio::TreeIO(-file => ">$tree_file.pruned.svg", -format => 'svggraph');
	$out_tree_io->write_tree($tree);
#	use Bio::Tree::Draw::Cladogram;
#	my $c = Bio::Tree::Draw::Cladogram->new(-bootstrap => 1, -tree => $tree, -compact => 0);
#	$c->print(-file => "$tree_file.cladogram.eps");
}

sub pruneTree($$) {
	my ($tree, $nodes_to_keep) = @_;
	NODE: foreach my $node ($tree->get_nodes()) {
		my $has_kept_children;
		foreach my $child ($node->get_all_Descendents) {
			next NODE if $$nodes_to_keep{$child->id};
		}
		$tree->splice($node) unless $$nodes_to_keep{$node->id} or $tree->get_root_node == $node;
		# TODO: why doesn't this work?
		# $tree->splice(-remove_id => [$node->id], -preserve_lengths => 1) unless $$nodes_to_keep{$node->id} or $tree->get_root_node == $node;
	}
	return $tree;
}
