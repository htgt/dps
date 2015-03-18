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
my $pre_region_param = 150; # unless overridden below
my $post_region_param = 150;

GetOptions(
    'input=s'         	=> \$input_file_param,
    'output=s'			=> \$output_file_param,
    'species=s'       	=> \$species_param,
    'pre_region=i'		=> \$pre_region_param,
    'post_region=i'		=> \$post_region_param,
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
unless ($pre_region_param) {
	print "Specify pre_region\n";
	$doomed++;
}
unless ($post_region_param) {
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
        -host => $ENV{LIMS2_ENSEMBL_HOST} || 'ensembldb.ensembl.org',
        -user => $ENV{LIMS2_ENSEMBL_USER} || 'anonymous'
    );

my $exon_adaptor = $registry->get_adaptor('Mouse', 'Core', 'Exon');

my $slice_adaptor = $registry->get_adaptor('Mouse', 'Core', 'Slice');

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

    my $ov = {};

	say 'Gene: ' . $gene_id;
    say 'Exon: ' . $exon_id;

    $ov->{gene_id} = $gene_id;
    $ov->{exon_id} = $exon_id;
 	$ov->{pre_chr_name} = $exon->slice->seq_region_name();
 	$ov->{exon_chr_start} = $exon->start;
    $ov->{pre_chr_start} = $exon->start - ($pre_region_param + 1); # starts at 1 not zero
    $ov->{pre_chr_end} = $exon->start - 1; # parameterize this?

    $ov->{exon_chr_end} = $exon->end;
    $ov->{post_chr_name} = $exon->slice->seq_region_name;
    $ov->{post_chr_start} = $exon->end + 1; # parameterize this?
    $ov->{post_chr_end} = $exon->end + $post_region_param;

    # Now get a slice that covers the whole region we have defined and check if there are any other exons
    # that overlap that region.

    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $ov->{pre_chr_name}, $ov->{pre_chr_start}, $ov->{post_chr_end} );

    my $transcripts = $slice->get_all_Transcripts();

    foreach my $trans ( @{$transcripts} ) {
        say '-- transcript: ' . $trans->display_id;
        my $exons = $trans->get_all_Exons();

        foreach my $sample_exon ( @{$exons}){
            say '- sample_exon id: ' . $sample_exon->display_id;
            if ( $sample_exon->display_id eq $ov->{exon_id} ) {
                say '++ This exon is the critical exon in the same transcript';
                next;
            }
            my $trans_exon = $sample_exon->transform('chromosome');
            if (my $new_left = is_left_overlap( $trans_exon, $ov))
             {
                $DB::single=1;
                say '+ Found pre_exon overlap with ' . $trans_exon->display_id;
                say '+ old pre_chr_start was: ' . $ov->{pre_chr_start};
                $ov->{pre_chr_start} = $new_left + 1;
                say '+ new pre_chr_start is:  ' . $ov->{pre_chr_start};
            } elsif ( my $new_right = is_right_overlap( $trans_exon, $ov) ) {
                $DB::single=1;
                say ' Found post_exon overlap ' . $trans_exon->display_id;
                say '+ old post_chr_end was: ' . $ov->{post_chr_end};
                $ov->{post_chr_end} = $new_right - 1;
                say '+ new post_chr_end is:  ' . $ov->{post_chr_end};
            }
        }
        say '---';
    }
    
    push @out_array, [
    		$ov->{gene_id},
    		$ov->{exon_id},
    		$ov->{exon_chr_start},
    		$ov->{exon_chr_end},
    		$ov->{pre_chr_name},
    		$ov->{pre_chr_start},
    		$ov->{pre_chr_end},
    		$ov->{post_chr_name},
    		$ov->{post_chr_start},
    		$ov->{post_chr_end},
    		];
}

$csv->eol("\n");
open $fh, ">:encoding(utf8)", $output_file_param or die $output_file_param . ": $!";
$csv->print ($fh, $_) for @out_array;
close $fh or die $output_file_param . ": $!";



# Needs a diagram!
#
#                                a                            f
#                    ||||TE|||||||                            |||||TE|||||||
#                             c          b            d           e
#                             <----------|||||CE|||||||----------->
#
# for left overlap following must both be true:
# (b - a) > 0
# (a - c) > 0
#
# for right overlap following must both be true:
# (f - d) > 0
# (e - f) > 0
#
#


sub is_left_overlap {
    my $test_exon = shift;
    my $ce_data = shift;

    my $a = $test_exon->end;
    my $b = $ce_data->{exon_chr_start};
    my $c = $ce_data->{pre_chr_start};

    my $ret_val;

    $ret_val = $a if ( (($b - $a ) > 0) and (($a - $c) > 0 ) );

    return $ret_val;
}

sub is_right_overlap {
    my $test_exon = shift;
    my $ce_data = shift;

    my $f = $test_exon->start;
    my $d = $ce_data->{exon_chr_end};
    my $e = $ce_data->{post_chr_end};

    my $ret_val;

    $ret_val = $f if ( (($f - $d) > 0 ) and (($e - $f) > 0) );

    return $ret_val;
}


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

=head1 AUTHOR

D J Parry-Smith
Copyright Welcome Trust Genome Institute 2015


=head1 DESCRIPTION

As synopsis. All command line options are required.

=cut