#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Set::IntervalTree;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Usage:

perl FEATURE_SATURATION_CURVE.pl -g file.gff3 -s samfile.sam -f gene -n 1000000

Description:

Purpose:       Determine the saturation curve for any feature type.
 
Input:         -s -- A sorted or unsorted SAM alignment file.  
               -g -- A sorted or unsorted valid GFF3 file.
               -f -- A feature type.  The feature type MUST match what is 
               in column 3.  GENE will not match/find gene.
               -n -- The lowest sampling frequency.  100 = 1/100 reads sampled; 
               1000 = 1/1000 reads sampled;

Output:        Output writes two columns to standard out.  The first column 
               is the probability of a read being sampled. The second is the proportion 
               of 'type/feature' covered. Plotting: column 1 = x and column 2 = y.  


";


my ($help);
my $gff3;
my $sam;
my $feature_type;
my $lower;
my $opt_success = GetOptions('help'            => \$help,
			     'gff3=s'          => \$gff3,
			     'sam=s'           => \$sam,
			     'feature_type=s'  => \$feature_type,
			     'lower=s'         => \$lower
    );

die $usage if $help || ! $opt_success;
die $usage unless $sam && $gff3 && $feature_type && $lower;


my ($search_tree, $features) = read_gff($gff3, $feature_type);
my $coverage_by_id           = read_sam($sam, $search_tree, $features);
generate_dist($coverage_by_id, $features);

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
sub read_gff {
  
    my %CHR_TREE;
    my %FEATURES;

    my ($file, $feature_type) = @_;
    open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";
    GFF3_LINE: while(<$IN>){
	chomp;
	next GFF3_LINE if $_ =~/^\#/;
	my @line = split /\t/, $_;
	next GFF3_LINE unless $line[2] =~ /$feature_type/;
	my @attributes = split /;/, $line[8];
	$CHR_TREE{$line[0]} = Set::IntervalTree->new if !defined $CHR_TREE{$line[0]};
	$CHR_TREE{$line[0]}->insert($attributes[0], $line[3], $line[4]);
	$FEATURES{$attributes[0]} = 0;
    }
    close $IN;
    return (\%CHR_TREE, \%FEATURES);
}
#-----------------------------------------------------------------------------
sub read_sam{  
    my ($file, $tree) = @_;
    open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";
    
    my %TREE = %{$tree};
    my %COV_HASH;
    
  SAM_LINE: while(<$IN>){
      chomp;
      next SAM_LINE if $_ =~ /^@/;
      my @l = split /\t/, $_;
      my @ids = @{$TREE{$l[2]}->fetch($l[3], $l[3] + length($l[9]))};
      
      for(my $i = $lower; $i >= .01; $i = $i /10){
	  my $rand = int(rand($i));
	  if ($rand == 1){
	      foreach my $id (@ids){
		  $COV_HASH{$i}{$id} = 1;
	      }
	  }
      }
  }
    close $IN;
    return \%COV_HASH;
}
#-----------------------------------------------------------------------------
sub generate_dist {
    my ($counts, $features) = @_;

    my $n_features = scalar (keys %{$features});
    while(my($key, $class) = each %{$counts}){
	my $count = 0;
	while(my ($classv, $value) = each %{$class}){
	    $count++ if $counts->{$key}{$classv};
	}
	my $fraction_features_covered = 0;
	my $fraction_reads_sampled    = 0;
	$fraction_features_covered    = $count / $n_features if $count ne 0;
	$fraction_reads_sampled       = 1 / $key if $key ne 0;
	print "$fraction_reads_sampled\t$fraction_features_covered\n";
    }
}
