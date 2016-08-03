#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

multiFasta_size_split.pl my.fasta -b 0 -e 100 -o seqs_0_100.fa

Description:

Right closes size extraction.

Options for single size range:
 -f multi fasta file
 -b start of size fraction (inclusive)
 -e end of size fraction (exclusive)
 -o output name, optional; otherwise [seqs_b_e.fa]

Options for automatic size range:
 -m max size
 -s step
";


my ($help);
my $file;
my $start;
my $end;
my $output;
my $max;
my $step;
my $opt_success = GetOptions('help'     => \$help,
                             "file=s"   => \$file,
                             "begin=s"  => \$start,
                             "end=s"    => \$end,
                             "output=s" => \$output,
                             "max=s"    => \$max,
                             "step=s"   => \$step,
			      );

die $usage if $help || ! $opt_success;

die $usage unless defined $file;
die $usage unless (defined $max && defined $step)
    || (defined $start && defined $end);




my $auto = defined $max && defined $step;

if($auto){
    for(my $i = 0; $i < ($max - $step); $i += $step){
        $start = $i;
        $end   = $start + $step;
        print STDERR "$start\t$end\n";
        $output = "seq\_$start\_$end.fasta";
        getSizeFraction();
    }
}
else{
    getSizeFraction();
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub getSizeFraction{

    if(! defined $output){
        $output = "seq\_$start\_$end";
    }
    my $IN  = Bio::SeqIO->new(-file  => $file);
    my $OUT = Bio::SeqIO->new(-file => ">$output", -format=>'fasta');
        while (my $seq = $IN->next_seq) {
            my $l = $seq->length;
                print STDERR "$start - $end : seq len : $l\n";
            if($l >= $start && $l < $end){

                $OUT->write_seq($seq);
            }
        }
}
