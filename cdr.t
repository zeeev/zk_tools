#!/usr/bin/perl

use lib '/home/zkronenb/tools/zk_tools/';
use CDR;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

cdr.t my.cdr

Description:

For development of a CDR object oriented code.

";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success;

my $file = shift;
die $usage unless $file;

my $cdr = CDR->new('file' => $file);
my $whatever = $cdr->Get_Seqids();
my $pragma = $cdr->Get_Pragma();

my $s = "chr1"; 
$cdr->Query_Range($s);    





#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}

