#!/usr/bin/perl
use Inline C;
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

print "CHRM\tPOS\tMA_COUNT_POP1\tMA_COUNT_POP2\tCOUNT_POP1\tCOUNT_POP2\tQ_FREQ_POP1\tQ_FREQ_POP2\tHET_POP1\tHET_POP2\tFST\tFISHER\'S_EXACT\tPERMUTATION\n";

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

	print "$line_dat->{CHRM}\t$line_dat->{POS}\t$res->{b}{a}\t$res->{t}{a}\t$res->{b}{count}\t$res->{t}{count}\t$res->{b}{bq}\t$res->{t}{tq}\t$res->{b}{het}\t$res->{t}{het}\t$res->{fst}\t$res->{fishers}\t$res->{per}\n";
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

    my $b_qcount = $background[0];
    my $t_qcount = $target[0];

    my $b_pcount = $background[1] - $background[0];
    my $t_pcount = $target[1] - $target[0]; 

    my $b_2pq = 2 * (1 - $b_q) * $b_q;
    my $t_2pq = 2 * (1 - $t_q) * $t_q;

    my $n_tot       = $background[1] + $target[1];
    my $tot_q_count = $background[0] + $target[0]; 
    my $tot_q       = $tot_q_count   / $n_tot;

#watch out for pool sizes that are teh same;

    my $Gs = (0.5 * ($b_q**2 + (1 - $b_q)**2)) + (0.5 * ($t_q**2 + (1 - $t_q)**2));
    #my $Gs = (($background[1] / $n_tot) * ($b_q**2 + (1 - $b_q)**2)) + (($target[1] / $n_tot) * ($t_q**2 + (1 - $t_q)**2));
    my $Gt = ((1 - $tot_q)**2) + ($tot_q**2);

    my $Fst = 'NA';
    $Fst = 0 if $b_q - $t_q == 0;
    if ($Gt == 0 || ($Gs - $Gt) == 0){
    }else{
	$Fst = ($Gs - $Gt) / (1 - $Gt);
    }
    
#    print "$t_qcount\t$t_pcount\t$b_qcount\t$b_pcount\n";

    #two sided
    my $fishers = fet::fet($t_qcount,$t_pcount,$b_qcount,$b_pcount, 1);




#    double add(int t_mac, int t_tot, int b_mac, int b_tot) {

    my $per = permute_counts($target[0], $target[1], $background[0], $background[1]); 
	


    $STRUCT{per}      = $per;
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

#-----------------------------------------------------------------------------
__END__
__C__

    double permute_counts(int t_mac, int t_tot, int b_mac, int b_tot) {
	

/*	printf("target_mac: %d\t", t_mac);      */
/*	printf("target_tot: %d\t", t_tot); 	*/
/*	printf("background_mac: %d\t", b_mac);  */
	



	int t_makeup[t_tot];
	int b_makeup[b_tot];
	
	int total_count = t_tot + b_tot;

/*	printf("total: %d\n", total_count); */

	int total[total_count];
	
	int i;
	for(i = 0; i < t_mac; i++){
	    t_makeup[i] = 1;
	}
	int j;
	for(j = t_mac; j < t_tot; j++){
	    t_makeup[j] = 0;
	}

/*	printf("target_makeup:\n");       */
/*	int v;				  */
/*	for(v = 0; v < t_tot; v++){	  */
/*            printf("%d", t_makeup[v]);  */
/*        }				  */
/*	printf("\n");                     */



	int k;
	for(k = 0; k < b_mac; k++){
	    b_makeup[k] = 1;
	}
	int l;
	for(l = b_mac; l < b_tot; l++){
	    b_makeup[l] = 0;
	}


/*	printf("background_makeup:\n");   */ 
/*	int g;				  */
/*	for(g = 0; g < b_tot; g++){	  */
/*            printf("%d", b_makeup[g]);  */
/*        }				  */
/*        printf("\n");                   */

	int n;
	for(n = 0; n < t_tot; n++){
            total[n] = t_makeup[n]; 
	}
        int o;
	for(o = 0; o < b_tot; o++){
	    int bspot    = t_tot + o;
	    total[bspot] = b_makeup[o]; 
	}

/*	printf("total_makeup:\n");         */
/*	int x;				   */
/*	for(x = 0; x < total_count; x++){  */
/*	    printf("%d", total[x]);	   */
/*	}				   */
/*	printf("\n");                      */

	
	long double n_total = 0;
        long double n_per   = 0;	
        while(1 == 1){ 	
	    ++n_per;
	    double n_hits = 0;     
	    if (n_total >= 10){
		break;
	    }
	    if (n_per >= 10000000){
		break;
	    }
	    shuffle(total, total_count);
	    int z;
	    for(z = 0; z < t_tot; z++){
           /*   printf("%d", total[z]); */
		if(total[z] == 1){
		    ++n_hits;
		}
	    }
	   /* printf("\n"); */
	    if(n_hits >= t_mac){
		++n_total;
	    }
	}
 	
	
	long double p_val = n_total / n_per;
/*	printf("pvalue: %f\n", p_val);  */
	return p_val;
}
    
static int rand_int(int n) {
    int limit = RAND_MAX - RAND_MAX % n;
    int rnd;

    do {
	rnd = rand();
    } while (rnd >= limit);
    return rnd % n;
}

/* I changed to shuffle as needed */
/* n - 1 -> n */
void shuffle(int *array, int n) {
    int i, j, tmp;

    for (i = n - 1; i > 0; i--) {
	j = rand_int(i + 1);
	tmp = array[j];
	array[j] = array[i];
	array[i] = tmp;
    }
}
