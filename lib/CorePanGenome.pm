#!/usr/bin/env perl

# Author: Andrey Kislyuk (kislyuk@gmail.com)

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
use lib "$FindBin::RealBin/../lib";
use AKUtils qw(logmsg);
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
#use Statistics::Descriptive;

1;

sub loadGenesFromInput($) {
	my ($input_files) = @_;
	my %genes;
	foreach my $file (@$input_files) {
		logmsg "Loading data from $file...";
		my ($org_filename) = fileparse($file); $org_filename =~ s/\.(gb|gbk)$//;
		if ($genes{$org_filename}) {
			warn "WARNING: duplicate organism id $org_filename, skipping file $file";
			next;
		}
		my $gbh = Bio::SeqIO->new(-format => 'genbank', -file => $file);
		while (my $gb_ent = $gbh->next_seq()) {
			logmsg "Loading GenBank entry ".$gb_ent->id." from $org_filename";
			# warn unless $gb_ent->species;
			for my $cds ($gb_ent->get_SeqFeatures) {
				next if $cds->primary_tag ne 'CDS';
				my $locus_name = ($cds->get_tag_values('locus_tag'))[0];
				# WARNING: this may have nasty repercussions
				#$locus_name =~ s/\s//g;
				die("Error: space in sequence name \"$locus_name\"") if $locus_name =~ /\s/;
				if ($genes{$org_filename}->{$locus_name}) {
					warn "WARNING: duplicate locus id $locus_name in $org_filename, skipping";
					next;
				}
				$cds->seq_id($org_filename);
#				print $cds->seq_id." id \n";
				$genes{$org_filename}->{$locus_name} = $cds;
			}
			logmsg keys(%{$genes{$org_filename}})." genes loaded";
		}
	}
	return \%genes;
}
