#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );
use Data::Dumper;
use Math::Round qw( round );

use LIMS2::Util::EnsEMBL;

die "Usage: $0 <species> <genes.csv>" unless @ARGV == 2;

my $species = shift;
my $e = LIMS2::Util::EnsEMBL->new( species => $species );

<>;
my @genes;
while ( my $line = <> ) {
    chomp $line;
    my ( $ms, $ens, $exons ) = split /,/, $line;

    next if $ms =~ /^#/; #allow commented lines
    say STDERR "Processing ($ens) $ms";
    my $gene = $e->gene_adaptor->fetch_by_stable_id($ens);
    say STDERR "$ms ens_gene is " . $gene->external_name;
    push @genes, $e->gene_adaptor->fetch_by_stable_id($ens);
}

#my @genes = ( $e->gene_adaptor->fetch_by_stable_id("ENSSSCG00000017529") );

#say $gene->canonical_transcript->slice->coord_system->name;

# ensembl human 'my $g = $e->gene_adaptor->fetch_by_stable_id("ENSG00000139618"); my $t = $g->canonical_transcript; my $start = $g->canonical_transcript->cdna_coding_start; my $end = $g->canonical_transcript->cdna_coding_end; say $start; say $end; my (@locs) = $g->canonical_transcript->cdna2genomic($start+100, $end); say join ", ", map { $_->start."-".$_->end } @locs; say "";say $g->canonical_transcript->coding_region_start; say $t->five_prime_utr->length'

say join ",", qw( marker_symbol gene_id num_exons strand transcript_id transcript_midpoint exon_id rank chr coding_start coding_end seq num_crisprs );

my $count;

#my $genes = $e->gene_adaptor->fetch_all_by_biotype( 'protein_coding' );
#while ( my $gene = shift @{ $genes } ) {
for my $gene ( @genes ) {
    say STDERR "Processed $count genes" if ++$count % 1000 == 0;

    if ( $gene->seq_region_name !~ /^(?:\d+|X|Y|MT)$/ ) {
        say STDERR "Skipping gene on " . $gene->seq_region_name;
        next;
    }

    my $transcript = $gene->canonical_transcript;

    my ( $cdna_start, $cdna_end ) = ( $transcript->cdna_coding_start, $transcript->cdna_coding_end );

    die "Canonical transcript of " . $gene->stable_id . " isn't protein coding ???"
        unless $cdna_start and $cdna_end;

    if ( $cdna_start+100 >= $cdna_end ) {
        say STDERR "Gene " . $gene->stable_id . " is too short, skipping.";
        say join ",", $gene->external_name,
                      $gene->stable_id,
                      scalar( @{ $transcript->get_all_Exons } ),
                      $gene->seq_region_strand,
                      $transcript->stable_id,
                      ("") x 7;
        next;
    }

    my $transcript_data = build_transcripts_hash( $gene->get_all_Transcripts );

    #die Dumper($transcript_data);

    unless ( %{$transcript_data} ) {
        say STDERR "Skipping " . $gene->stable_id . " as it has no coding transcripts";
        next;
    }

    #get the mid point in cdna co-ordinates, this is the same as
    #$cdna_start + ( ( $t->cdna_coding_end - $t->cdna_coding_start ) / 2 )
    my $cdna_mid = $cdna_start + round( length($transcript->translateable_seq)/2 );

    #get the chromosome coordinate for our cdna middle
    my ( $middle ) = $transcript->cdna2genomic($cdna_mid, $cdna_mid);


    #TODO
    # add the real canonical transcript check

    my @locs = $transcript->cdna2genomic($cdna_start+100, $cdna_end);

    my ( $rank, $location, $num_found ) = (0) x 3;
    my @exons = @{ $transcript->get_all_Exons };
    for my $exon ( @exons ) {
        #say $exon->stable_id;
        ++$rank;

        my $start = $exon->coding_region_start( $transcript );
        my $end   = $exon->coding_region_end( $transcript );
        next unless $start and $end; #skip non coding

        #get the location we're considering
        my $region = $locs[$location];

        #say STDERR $gene->stable_id . " has " . scalar(@exons) . " exons";
        #say STDERR "Pulled region " . ($location);
        #say STDERR Dumper( $region );

        #if this region is not contained by the exon then skip it as it
        #is outside that 100bp buffer we set.
        #we will try again on the next one
        next unless $region->start >= $start and $region->end <= $end;

        ++$location;

        next unless exon_is_constitutive( $transcript_data, $exon->stable_id );

        #next unless $exon->is_constitutive; #skip non constitutive

        $num_found++;

        if ( $region->start > $region->end ) {
            say STDERR $gene->stable_id;
            say STDERR $exon->stable_id;
            die Dumper($region);
        }

        #we are on chromosome coord system so this is fine
        #say STDERR join ",", $gene->seq_region_name, $start, $end;
        my $slice = $e->slice_adaptor->fetch_by_region(
            "chromosome",
            $gene->seq_region_name,
            $region->start,
            $region->end
        );

        #get naive crispr count by counting the number of GG/CCs
        my $seq = $slice->seq;
        my $num_crisprs = () = $seq =~ /(GG|CC)/g;

        #say $exon->stable_id;
        #say $start, ", ", $end;
        #say $exon->coord_system_name;

        say join ",", $gene->external_name,
                      $gene->stable_id,
                      scalar( @exons ),
                      $gene->seq_region_strand,
                      $transcript->stable_id,
                      $middle->start, #display the halfway transcript point
                      $exon->stable_id,
                      $rank,
                      $gene->seq_region_name,
                      $region->start,
                      $region->end,
                      $seq,
                      $num_crisprs;

        #die if $gene->seq_region_strand == -1;
    }

    unless ( $num_found ) {
        say STDERR "BLANK ROW";
        say join ",", $gene->external_name,
                      $gene->stable_id,
                      scalar( @exons ),
                      $gene->seq_region_strand,
                      $transcript->stable_id,
                      ("") x 7;
    }

    #die "Couldn't find any exons for " . $gene->stable_id unless $num_found;
}

sub exon_is_constitutive {
    my ( $transcript_data, $exon_id ) = @_;

    while ( my ( $transcript, $exons ) = each %{ $transcript_data } ) {
        return 0 unless exists $exons->{$exon_id};
    }

    return 1;
}

sub build_transcripts_hash {
    my ( $transcripts ) = @_;

    my %all_transcripts;

    for my $transcript ( @{ $transcripts } ) {
        #skip all non coding transcripts
        next unless $transcript->translation;
        next if $transcript->biotype eq 'nonsense_mediated_decay';

        $all_transcripts{$transcript->stable_id} = { map { $_->stable_id => 1 }
                                                    @{ $transcript->get_all_Exons } };
    }

    return \%all_transcripts;
}

sub valid_coding_transcript {
    my ( $transcript ) = @_;
    my $id = $transcript->stable_id;

    TRACE( "$id biotype: " . $transcript->biotype );
    if ( !$transcript->translation ) {
        TRACE("Transcript $id is non protein coding");
        return 0;
    }

    if ( $transcript->biotype eq 'nonsense_mediated_decay') {
        TRACE("Transcript $id is NMD");
        return 0;
    }

    # CDS incomplete check, both 5' and 3'
    if ( _get_transcript_attribute( $transcript, 'cds_end_NF' ) ) {
        TRACE("Transcript $id has incomplete CDS end");
        return 0;
    }

    if ( _get_transcript_attribute( $transcript, 'cds_start_NF' ) ) {
        TRACE("Transcript $id has incomplete CDS start");
        return 0;
    }

    TRACE("Transcript $id is VALID");
    return 1;
}

sub _get_transcript_attribute {
    my ( $transcript, $code ) = @_;

    my ( $attr ) = @{ $transcript->get_all_Attributes($code) };
    if ( $attr ) {
        return $attr->value();
    }
    return 0;
}

1;
