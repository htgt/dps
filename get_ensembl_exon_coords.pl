#!/usr/bin/perl

use Text::CSV;
use Bio::EnsEMBL::Registry;
use feature qw/ say /;
use strict;
use warnings;
use Getopt::Long;

my $input_file_param;
my $output_file_param;
my $species_param;
my $pre_region = 150;
my $post_region = 150;

GetOptions(
    'input=s'         	=> \$input_file_param,
    'output=s'			=> \$output_file_param,
    'species=s'       	=> \$species_param,
    'pre_region=i'		=> \$pre_region,
    'post_region=i'		=> \$post_region,
) or die usage_message();

sub usage_message {
    print <<'END_USAGE';
usage: get_ensembl_coords.pl
	--input=<input.csv>
	--output=<output.csv>
	--species=[ Human | Mouse ]
	--pre_region=150
	--post_region=150
END_USAGE
return;
}

my $doomed;

unless ($input_file_param) {
	print "Specify input_file\n";
	$doomed++;
}
unless ($output_file_param) {
	print "Specify output_file\n";
	$doomed++;
}
unless ($species_param) {
	print "Specify species\n";
	$doomed++;
}
unless ($pre_region) {
	print "Specify pre_region\n";
	$doomed++;
}
unless ($post_region) {
	print "Specify post_region\n";
	$doomed++;
}
if ($doomed) {
	usage_message();
	die;
}

my @gene_array;
my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
                or die "Cannot use CSV: ".Text::CSV->error_diag ();
 
open my $fh, "<:encoding(utf8)", $input_file_param or die $input_file_param . ": $!";
while ( my $gene_row = $csv->getline( $fh ) ) {
    next if $gene_row->[0] eq 'Gene ID'; # deals with the header line
    push @gene_array, $gene_row;
}
$csv->eof or $csv->error_diag();
close $fh;

say 'Read ' . scalar @gene_array . ' rows from ' . $input_file_param;

say 'Connecting to the ensembl registry';

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
        -host => $ENV{LIMS2_ENSEMBL_HOST} || 'ensembldb.internal.sanger.ac.uk',
        -user => $ENV{LIMS2_ENSEMBL_USER} || 'anonymous'
    );

my $exon_adaptor = $registry->get_adaptor('Mouse', 'Core', 'Exon');

my @out_array;
push @out_array,
	[	'Gene Id',
		'Exon Id',
		'chr_start',
		'chr_end',
		'pre_chr_name',
		'pre_chr_start',
		'pre_chr_end',
		'post_chr_name',
		'post_chr_start',
		'post_chr_end',
		];

foreach my $gene_row ( @gene_array ) {
    my $gene_id = $gene_row->[0];
    my $exon_id = $gene_row->[1];
$DB::single=1;
    my $exon = $exon_adaptor->fetch_by_stable_id( $exon_id );
    if ( ! $exon ) {
    	say 'No exon stable id found: ' . $exon_id . ' (' . $gene_id . ')';
    	next; 
    }
    my ($pre_chr_name, $pre_chr_start, $pre_chr_end);
    my ($post_chr_name, $post_chr_start, $post_chr_end);
    my ($exon_chr_start, $exon_chr_end);

	say 'Exon: ' . $exon_id;

 	$pre_chr_name = $exon->slice->seq_region_name();
 	$exon_chr_start = $exon->start;
    $pre_chr_start = $exon->start - 151;
    $pre_chr_end = $exon->start - 1;

    $exon_chr_end = $exon->end;
    $post_chr_name = $exon->slice->seq_region_name;
    $post_chr_start = $exon->end + 1;
    $post_chr_end = $exon->end + 150;

    
    push @out_array, [
    		$gene_id,
    		$exon_id,
    		$exon_chr_start,
    		$exon_chr_end,
    		$pre_chr_name,
    		$pre_chr_start,
    		$pre_chr_end,
    		$post_chr_name,
    		$post_chr_start,
    		$post_chr_end,
    		];
}

$csv->eol("\n");
open $fh, ">:encoding(utf8)", $output_file_param or die $output_file_param . ": $!";
$csv->print ($fh, $_) for @out_array;
close $fh or die $output_file_param . ": $!";

=head1 NAME

get_ensembl_coords.pl - find the ensembl coordinates for defined regions before and after the input exon

=head1 SYNOPSIS

get_ensembl_coords.pl
	--input=<input.csv>
	--output=<output.csv>
	--species=[ Human | Mouse ]
	--pre_region=150
	--post_region=150

Given an input csv file of form 'Gene ID', 'Exon ID' produce an output csv file:
'pre_chr_name','pre_chr_start','pre_chr_end','post_chr_name','post_chr_start','post_chr_end'


=head1 DESCRIPTION

As synopsis. All command line options are required.

=cut