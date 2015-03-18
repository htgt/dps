#!/usr/bin/env perl
use strict;
use warnings;

use WGE::Model::DB;
use LIMS2::Model;
use Data::Dumper;
use Path::Class;
use LIMS2::Util::EnsEMBL;
use Text::CSV;
use Log::Log4perl ':easy';
use feature qw( say );
use YAML::Any;
use List::Util qw( sum min max );

Log::Log4perl->easy_init( { level => $DEBUG } );

die "Usage: rank_crisprs_human.pl <exons.csv> <species>" unless @ARGV == 2;

#should be pair or crispr

#i couldn't be bothered to think of a clever way to do this
#basically how to distribute the crisprs over the exons
# this array of array refs - each arrayref represents the dist of crisprs amongst the
# number of 
my @dists = (
    [], # zero critical exons - so get no crisprs
    [2], # 1 critical exon - I want 2 crisprs from one critical exon
    [1, 1], # 2 critical exons - I want 1 crispr from one and 1 crispr from the other
);

my $total_snps = 0;

use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $sa = $reg->get_adaptor("human", "core", "slice");
my $vfa = $reg->get_adaptor("human", "variation", "variationfeature");

my $w = WGE::Model::DB->new;

my $species = $ARGV[1];
my $SPECIES_ID = $w->resultset('Species')->find( { id => $species } )->numerical_id;

die "Couldn't find species $species" unless $SPECIES_ID;

my $csv = Text::CSV->new( { eol => "\n" });
my $fh = file( $ARGV[0] )->openr or die( "Can't open $ARGV[0] " . $! );

# input csv file must have column headers
$csv->column_names( @{ $csv->getline( $fh ) } );

my @original_fields = $csv->column_names;
my $error_fh = file( 'human_genes_needing_crisprs.csv' )->openw;
$csv->print( $error_fh, [ @original_fields ] );

my @crispr_fields = qw(
    strand
    id
    exon_id
    num_exons
    rank
    transcript_midpoint
    chr_name
    chr_start
    chr_end
    seq
    pam_right
    off_target_summary
);

say join ",", "gene_id", @crispr_fields;

my $total = 0;
my $current_gene;
my @current_crisprs;
my @current_rows;
my %problem_genes;
my $num_crisprs = 50;
my $num_too_far = 0;
my $min_crisprs = 2;
while ( my $data = $csv->getline_hr( $fh ) ) {
    if ( ! $current_gene ) {
        $current_gene = $data->{gene_id};
        say STDERR "Processing $current_gene";
    }
    elsif ( $current_gene ne $data->{gene_id} ) {

        my @missing = grep { ! $_->{off_target_summary} } @current_crisprs;
        if ( @missing ) {
            say $_->{id} for @missing;
            say STDERR "Have " . scalar( @current_crisprs ) . " crisprs:";
        }

        process_crisprs();

        $current_gene = $data->{gene_id};
        @current_crisprs = ();
        @current_rows = ();

        say STDERR "Processing $current_gene";
    }

    #keep track of all of this genes crisprs
    push @current_rows, $data;

    next if @current_crisprs >= $num_crisprs;

    next unless $data->{chr} and $data->{coding_start} and $data->{coding_end};

    #only allow crisprs up to the midpoint
    # if ( $data->{strand} == 1 ) {
    #     #dont bother looking if the exon is past halfway
    #     next if $data->{coding_start} > $data->{transcript_midpoint};
    # }
    # elsif ( $data->{strand} == -1 ) {
    #     next if $data->{coding_end} < $data->{transcript_midpoint};
    # }
    # else {
    #     die "Invalid strand";
    # }

    my $crisprs = crisprs_for_region( {
        chr_name  => $data->{chr},
        chr_start => $data->{coding_start},
        chr_end   => $data->{coding_end},
        strand    => $data->{strand},
    } );

    #data to add into the crispr hash for easy printing
    my $exon_data = {
        num_exons           => $data->{num_exons},
        rank                => $data->{rank},
        exon_id             => $data->{exon_id},
        transcript_midpoint => $data->{transcript_midpoint},
        strand              => $data->{strand},
    };

    my @valid_crisprs =
        grep { crispr_in_target_region( $_, $data->{coding_start}, $data->{coding_end}, $exon_data ) }
            map { $_->as_hash }
                @{ $crisprs };

    my $remaining = $num_crisprs - @current_crisprs;

    #add up to 15 crisprs for this gene
    if ( @valid_crisprs >= $remaining ) {
        push @current_crisprs, @valid_crisprs[0 .. ($remaining-1)];
    }
    else {
        push @current_crisprs, @valid_crisprs;
    }
}

say $_->{id} for grep { ! $_->{off_target_summary} } @current_crisprs;

#process last line as while loop will have finished before its processed
process_crisprs();

sub process_crisprs {
    $total++;
    say STDERR "Found " . scalar( @current_crisprs ) . " crisprs";
    #say STDERR $_->{id} for @current_crisprs;
    my @crisprs = grep { crispr_unique( $_ ) } @current_crisprs;

    if ( @crisprs < $min_crisprs ) {
        #say STDERR "Found " . scalar( @crisprs ) . " unique crisprs";
        for my $row ( @current_rows ) {
            #$problem_genes{ $data->{gene_id} } = scalar( @ranked );
            $csv->print( $error_fh, [ @{ $row }{ @original_fields } ] );
        }

        return;
    }
    elsif ( @crisprs == 2 ) {
        for my $crispr ( @crisprs ) {
            say join ",", $current_gene, @{ $crispr }{ @crispr_fields };
        }

        return;
    }

    #we have enough crisprs so try and choose the best ones

    #should group them when we add them tbh
    my %grouped;
    for my $crispr ( @crisprs ) {
        my $rank = $crispr->{rank};

        #at position x add this crispr
        push @{ $grouped{$rank} }, $crispr;
    }

    my @crisprs_per_exon;

    my $num_exons = min 2, scalar( keys %grouped );
    die "0 exons??" unless $num_exons;

    my $dist = $dists[$num_exons];

    #remainder keeps track of if we didnt hit the wanted target
    my ( $remainder, $i ) = ( 0, 0 );
    for my $rank ( sort keys %grouped ) {
        my $exon_crisprs = $grouped{$rank};
        my $crisprs_wanted = $dist->[$i++] + $remainder;
        say STDERR "Finding $crisprs_wanted for exon with rank $rank";
        my $chosen = crisprs_for_exon( $exon_crisprs, $crisprs_wanted );

        #set remainder for the next loop to pick up
        $remainder = $crisprs_wanted - @{ $chosen };

        #now we should just print out these bad boys
        for my $crispr ( @{ $chosen } ) {
            unless ( @{ $crispr }{ @crispr_fields } ) {
                die $current_gene . "\n" . Dumper( $crispr );
            }
            say join ",", $current_gene, @{ $crispr }{ @crispr_fields };
        }

        #we skip any remaining ones after the fifth exon (5 due to 0-4)
        last if $i >= 2;
    }
}

say STDERR "$total genes total";
say STDERR "Skipped a total of $total_snps crisprs due to SNPs";

=head crisprs_for_region

Find all the single crisprs in and around the target region.

=cut
sub crisprs_for_region {
    my ( $params ) = @_;

    DEBUG "Getting crisprs for " . region_str( $params );

    my $order_by;

    if ( $params->{strand} == 1 ) {
        $order_by = { -asc => 'chr_start' };
    }
    elsif ( $params->{strand} == -1 ) {
        $order_by = { -desc => 'chr_start' };
    }
    else {
        die "Invalid strand";
    }

    # we use 90 because the spaced between the crisprs in a pair can be 50 bases.
    # 50 + the size of 2 crisprs is around 90
    # that should bring back all the possible crisprs we want ( and some we do not want
    # which we must filter out )
    my @crisprs = $w->resultset('Crispr')->search(
        {
            'species_id'  => $SPECIES_ID,
            'chr_name'    => $params->{chr_name},
            # need all the crisprs starting with values >= start_coord
            # and whose start values are <= end_coord
			# flanks here are 25
            'chr_start'   => {
                -between => [ $params->{chr_start} - 25, $params->{chr_end} + 25 ],
            },
        },
        { order_by => $order_by }
    );

    #say $_->id for @crisprs;

    return \@crisprs;
}

=head crispr_in_target_region

Crispr is in target region if at least one base of the pam site is within
the target region.

=cut
sub crispr_in_target_region {
    my ( $crispr, $start, $end, $data_to_add ) = @_;

    for my $key ( keys %{ $data_to_add } ) {
        $crispr->{$key} = $data_to_add->{$key};
    }

    #only allow crisprs up to the midpoint
    # if ( $data_to_add->{strand} == 1 ) {
    #     #dont bother looking if the exon is past halfway
    #     return if $start > $data_to_add->{transcript_midpoint};
    #     $end = min $end, $data_to_add->{transcript_midpoint};
    # }
    # elsif ( $data_to_add->{strand} == -1 ) {
    #     return if $end < $data_to_add->{transcript_midpoint};
    #     $start = max $start, $data_to_add->{transcript_midpoint};
    # }
    # else {
    #     die "Invalid strand";
    # }

    # if ( $crispr->{pam_right} ) {
    #     if ( $crispr->{chr_end} > $start && $crispr->{chr_end} <= $end ) {
    #         return 1;
    #     }
    # }
    # else {
    #     if ( $crispr->{chr_start} >= $start && $crispr->{chr_start} < $end ) {
    #         return 1;
    #     }
    # }

    #must be completely contained
    if ( $crispr->{chr_start} >= $start && $crispr->{chr_end} <= $end ) {
        return 1;
    }

    return;
}

sub rank_crisprs {
    my ( $crisprs ) = @_;

    my @ranked;
    for my $crispr ( @{ $crisprs } ) {
        die $crispr->{id} unless $crispr->{off_target_summary};
        my $score = score_off_target_summary( $crispr->{off_target_summary} );

        push @ranked, { score => $score, crispr => $crispr };
    }

    return sort { $a->{score} <=> $b->{score} } @ranked;
}

sub score_off_target_summary {
    my ( $text ) = @_;

    die "off_target_summary is null" unless $text;

    #load yaml string
    my $summary = Load( $text );

    #the higher the number of mismatches the less weight it should have
    my @weights = ( 100, 50, 10, 1 );

    #get score by adding 0 + 1 + 2 + 3 fields
    return sum map { $summary->{$_} * $weights[$_] } 0 .. 3;
}

sub crispr_unique {
    my ( $crispr ) = @_;

    die "off_target_summary is null for " . Dumper( $crispr ) unless $crispr->{off_target_summary};

    return Load( $crispr->{off_target_summary} )->{0} <= 1;

    #return if its not unique
    return if Load( $crispr->{off_target_summary} )->{0} > 1;

    my $slice = $sa->fetch_by_region('chromosome', $crispr->{chr_name}, $crispr->{chr_start}, $crispr->{chr_end});

    #return false if it has a snp
    if ( slice_has_snps( $slice ) ) {
        ++$total_snps;
        DEBUG "Found SNP";
        return;
    }
    else {
        return 1;
    }
}

sub region_str {
    my $crispr = shift;

    unless ( $crispr->{chr_name} && $crispr->{chr_start} && $crispr->{chr_end} ) {
        die Dumper( $crispr );
    }

    return $crispr->{chr_name} . ':' . $crispr->{chr_start} . '-' . $crispr->{chr_end};
}

#attempt to extract x crisprs 10bp apart from a group of crisprs
sub crisprs_for_exon {
    my ( $crisprs, $crisprs_wanted ) = @_;

    die "crisprs_wanted is undefined" unless $crisprs_wanted;

    my $max_distance = 10;

    return $crisprs if @{ $crisprs } <= $crisprs_wanted;
    return [ $crisprs->[0] ] if $crisprs_wanted == 1; #just return the first crispr

    my @too_close;

    #first crispr we pick is the very first one in the exon
    my @chosen = shift @{ $crisprs };
    for my $crispr ( @{ $crisprs } ) {
        return \@chosen if @chosen == $crisprs_wanted;
        #say STDERR $crispr->{chr_start} . " " . $chosen[-1]->{chr_end};
        my $distance;
        if ( $crispr->{strand} == 1 ) {
            $distance = ( $crispr->{chr_start} - $chosen[-1]->{chr_end} ) - 1;
        }
        elsif ( $crispr->{strand} == -1 ) {
            #if its -ve the chosen crispr will always be further along
            $distance = ( $chosen[-1]->{chr_start} - $crispr->{chr_end} ) - 1;
        }
        else {
            die "Invalid strand";
        }

        if ( $distance >= $max_distance ) {
            #say STDERR "Crispr chosen with distance $distance";
            push @chosen, $crispr;
            #if we got as many crisprs as we need then we're done
        }
        else {
            #say STDERR "Crispr with distance $distance skipped";
            push @too_close, { distance => $distance, crispr => $crispr };
        }
    }

    say STDERR "Couldn't get $crisprs_wanted crisprs spaced nicely";

    #sort by distance descending
    #this isn't that great cause the distance to the NEXT crispr could be 0
    #but its good enough i think
    @too_close = sort { $b->{distance} <=> $a->{distance} } @too_close;
    my $remaining = ( $crisprs_wanted - @chosen ) - 1;
    push @chosen, map { $_->{crispr} } @too_close[0 .. $remaining];

    return \@chosen;
}

sub slice_has_snps {
    my ( $slice ) = @_;

    my @variation_features = @{ $vfa->fetch_all_by_Slice( $slice ) };

    #see if there is a 1000 genomes entry in all the variation features and their alleles
    for my $vf ( @variation_features ) {
        for my $allele ( @{ $vf->variation->get_all_Alleles } ) {
            if ( $allele->population && $allele->population->name =~ /1000GENOMES.*YRI/ ) {
                return 1;
            }
        }
    }

    return;
}

=head1 NAME

rank_crisprs_human.pl - find groups of crisprs pairs in loxp regions

=head1 SYNOPSIS

rank_crisprs_human.pl <exons.csv> <species>

get the "best" crisprs for a given list of exons

=head1 DESCRIPTION

Iterate through a csv of exons with co-ordinate data, pull up all crisprs for that region.
Loop through the CRISPRs, throw away those with low off target scores and SNPs,
attempt to take the best 2 crisprs that are left

The exons csv must have the following fields:
ensembl_exon_id
gene_id
chr
coding_start
coding_end
num_exons
rank
exon_id
transcript_midpoint
strand

=cut