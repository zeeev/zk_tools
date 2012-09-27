#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Math::BigInt;
use lib "/home/zkronenb/tools/zk_tools";
use fet;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

PARSE_SNVer.pl SNVer.pooled.filtered.vcf -t 1 -b 2 -d 10 > output.txt

Description:

file:            First argument to pass in.  This should be the direct output from SNVer. 
-t/--target:     The column(s) in the file corrisponding to the libraries to use as the target.
-b/--background: The column(s) in the file corrisponding to the libraries to use as the background.
-d/--depth:      The depth cutoff to use.  Both the target and background depths must be greater than -d.

";

my ($help);
my $t;
my $b;
my $depth;
my $opt_success = GetOptions('help'        => \$help,
			     'target=s'    => \$t,
			     'depth=s'     => \$depth,
			     'backround=s' => \$b,
);

die $usage if $help || ! $opt_success;

my $file = shift;
die $usage unless $file && $t && $b && $depth;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

my %PRAGMA;

$b = Parse_group($b);
$t = Parse_group($t);


print "CHRM\tPOS\tMA_COUNT_POP1\tMA_COUNT_POP2\tCOUNT_POP1\tCOUNT_POP2\tQ_FREQ_POP1\tQ_FREQ_POP2\tHET_POP1\tHET_POP2\tFST\n";

LINE: while (my $l = <$IN>) {
    chomp $l;
    if($l =~ /^#/){
    }else{
	my $line_dat = Parse_Genotype_Lines($l);
	my $t_dat    = Group_Column($line_dat, $t);
	my $b_dat    = Group_Column($line_dat, $b);

	next LINE if @{$b_dat}[1] < $depth;
	next LINE if @{$t_dat}[1] < $depth;
	
	
	my $res = Fst($b_dat, $t_dat);
#	my $g_stat = G_test($res);

	print "$line_dat->{CHRM}\t$line_dat->{POS}\t$res->{b}{a}\t$res->{t}{a}\t$res->{b}{count}\t$res->{t}{count}\t$res->{b}{bq}\t$res->{t}{tq}\t$res->{b}{het}\t$res->{t}{het}\t$res->{fst}\n";
    }
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
sub Parse_group{
    my $unsplit = shift;
    my @split;
    if($unsplit =~ /,/){
	@split  = split /,/, $unsplit; 
    }else{
	$split[0] = $unsplit;
    }
    return \@split;
}
#-----------------------------------------------------------------------------

sub Parse_Pragma {
    my $l = shift;
    $l /\#//g;
    map {$PRAGMA{$_} = $_} split /,|=/;
}

#-----------------------------------------------------------------------------

sub Parse_Genotype_Lines {
    
    my $l = shift;    
    my %line_dat;
    my @VCF_DAT = split "\t", $l;
    $line_dat{CHRM}   = shift @VCF_DAT;
    $line_dat{POS}    = shift @VCF_DAT;
    $line_dat{ID}     = shift @VCF_DAT;
    $line_dat{REF}    = shift @VCF_DAT;
    $line_dat{ALT}    = shift @VCF_DAT;
    $line_dat{QUAL}   = shift @VCF_DAT;
    $line_dat{FILTER} = shift @VCF_DAT;
    $line_dat{INFO}   = shift @VCF_DAT;
    $line_dat{FORMAT} = shift @VCF_DAT;
    Parse_Genotypes(\@VCF_DAT, \%line_dat);
    return \%line_dat;
}
#-----------------------------------------------------------------------------
sub Parse_Genotypes{
    
    my ($sample, $data_struct) = @_;
    
    my $id = 1;
    foreach my $pool (@{$sample}){
	my @ac_dp = split ":", $pool;
	$data_struct->{genotypes}{$id} = \@ac_dp;
	$id++;
    }
}

#-----------------------------------------------------------------------------

sub Fst{

# John H. Gillespie pg 132                                                                                        
    my ($background, $target) = @_;

    my @background = @{$background};
    my @target     = @{$target};

    my %STRUCT;

    my $b_q = $background[0] / $background[1];
    my $t_q = $target[0]     / $target[1];

    my $b_p = $background[1] - $background[0];
    my $t_p = $target[1] - $target[0]; 

    my $b_2pq = 2 * (1 - $b_q) * $b_q;
    my $t_2pq = 2 * (1 - $t_q) * $t_q;

    my $n_tot       = $background[1] + $target[1];
    my $tot_q_count = $background[0] + $target[0]; 
    my $tot_q       = $tot_q_count   / $n_tot;

    my $Gs = (($background[1] / $n_tot) * ($b_q**2 + (1 - $b_q)**2)) + (($target[1] / $n_tot) * ($t_q**2 + (1 - $t_q)**2));
    my $Gt = ((1 - $tot_q)**2) + ($tot_q**2);

    my $Fst = 'NA';
    $Fst = 0 if $b_q - $t_q == 0;
    if ($Gt == 0 || ($Gs - $Gt) == 0){
    }else{
	$Fst = ($Gs - $Gt) / (1 - $Gt);
    }
    

    my $fishers = fet::fet($b_q, $b_p, $t_q, $t_p, 1);

    $STRUCT{fst}      = $Fst;
    $STRUCT{fishers}  = $fishers;
    $STRUCT{b}{het}   = $b_2pq;
    $STRUCT{t}{het}   = $t_2pq;
    $STRUCT{t}{tq}    = $t_q;
    $STRUCT{b}{bq}    = $b_q;
    $STRUCT{b}{count} = $background[1];
    $STRUCT{t}{count} = $target[1];
    $STRUCT{b}{a}     = $background[0];
    $STRUCT{t}{a}     = $target[0];

    return \%STRUCT;

}

#-----------------------------------------------------------------------------
sub Group_Column{

    my($line_dat, $group_numbers) = @_;
        
    my @counts = (0,0);

    foreach my $g (@{$group_numbers}){
	$counts[1] += @{$line_dat->{genotypes}{$g}}[1] if @{$line_dat->{genotypes}{$g}}[1] ne 'NA';
	$counts[0] += @{$line_dat->{genotypes}{$g}}[0] if @{$line_dat->{genotypes}{$g}}[0] ne 'NA';
    }
    return \@counts;
}
#-----------------------------------------------------------------------------

sub G_test{

    my $data = shift;

#       $STRUCT{fst}      = $Fst;                                                                 
#       $STRUCT{b}{het}   = $b_2pq;                                                               
#       $STRUCT{t}{het}   = $t_2pq;                                                               
#       $STRUCT{t}{tq}    = $t_q;                                                                 
#       $STRUCT{b}{bq}    = $b_q;                                                                 
#       $STRUCT{b}{count} = $background[1];                                                       
#       $STRUCT{t}{count} = $target[1];                                                           

    my $q_obs_b = $data->{b}{bq};
    my $p_obs_b = 1 - $data->{b}{bq};
    my $q_exp_b = $data->{t}{tq};
    my $p_exp_b = 1 - $data->{t}{tq};

    my $q_obs_t = $data->{t}{tq};
    my $p_obs_t = 1 - $data->{t}{tq};
    my $q_exp_t = $q_obs_b;
    my $p_exp_t = $p_obs_b;

    my $sum = 0;

    $sum += log10(($q_obs_b / $q_exp_b)) if $q_obs_b > 0 && $q_exp_b > 0;
    $sum += log10(($p_obs_b / $p_exp_b)) if $p_obs_b > 0 && $p_exp_b > 0;
    $sum += log10(($q_obs_t / $q_exp_t)) if $q_obs_t > 0 && $q_exp_t > 0;
    $sum += log10(($p_obs_t / $p_exp_t)) if $p_obs_t > 0 && $p_exp_t > 0;
    
    my $Gstat = 'NA';
    $Gstat = 2 * $sum if $sum > 0;
    
    
    return $Gstat;
}

#-----------------------------------------------------------------------------

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
