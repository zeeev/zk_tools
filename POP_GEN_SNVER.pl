#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Tabix;
use Statistics::Descriptive;
use Math::BigFloat;
use POSIX;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "

Synopsis:

BIN_VALUE_BY_REGION_W_BUILD -b build.file.txt -d data_file.gz -w window -t 1,2 -b 3,4 -segregating 40 -depth 10

Description:

I bin non overlapping windows.  If a feature is shorter than the window size 
I will just use the whole feature.  use the value flag to tell me what data you want worked on (zero indexed).
Seg is the lower limit for the number of sites withing the windowsize.  less than 40 means I will through out that window.
Depth is the lower average bound over the window.
";

my ($help);
my $build;
my $t;
my $b;
my $window;
my $data;
my $call; 
my $seg;
my $depth;
my $opt_success = GetOptions('help'      => \$help,
			     'call=s'    => \$call,
			     'data=s'    => \$data,
			     'build=s'   => \$build,
			     'window=s'  => \$window,
			     'target=s'    => \$t,
                             'background=s' => \$b,
			     'segregating=s' => \$seg,
			     'depth=s'       => \$depth,
    );

my %BUILD_RECORD;
my %N_CHOOSE;

die $usage if $help || ! $opt_success;
die $usage unless $data && $build && $t && $b && $seg;

my %PRAGMA;

Parse_Build();
$t = Parse_group($t);
$b = Parse_group($b);

print "SEQID\tSTART\tSTOP\tGENOME_START\tGENOME_STOP\tSEG_SITES_IN_WINDOW\tW_t\tW_b\tF_t\tF_b\tP_t\tP_b\tH_t\tH_b\tW-test\tF-test\tP-test\tH-test\tD_FST\tMD\n";

my $regions = BIN_REGIONS($build, $window);
$regions    = PARSE_REGIONS($data, $regions, $t, $b);
my $results = SELECTION_REGIONS($regions);

#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub BIN_REGIONS{
    
    my %REGION_STRUCT;

    my ($build_file, $window) = @_;
    
    open (my $IN, '<', $build_file) or die "Can't open $build_file for reading\n$!\n";
    
    LINE: while (my $line = <$IN>) {
	chomp $line;
	my @line_contents = split /\s+/, $line;
	next LINE if ! defined $line_contents[0];

	if($line_contents[2] > $window + 1){
	    my $start  = 1;
	    my $end    = $start + $window;
	    my $slide  = int($window / 10);
	    while ($end  <= $line_contents[2]){
		$end   += $slide;
		$start += $slide;
		$REGION_STRUCT{$line_contents[0]}{"$start:$end"} = 'NA';
	    }
	}
	else{
       
	$REGION_STRUCT{$line_contents[0]}{"1:$line_contents[2]"} = 'NA';
	}
    }
    return \%REGION_STRUCT;
}
#-----------------------------------------------------------------------------
sub PARSE_REGIONS {
    
    my ($data_file, $regions, $t, $b) = @_;    
    my $tabix_obj = Tabix->new(-data=>$data_file); 
    my @data_file_seqs = $tabix_obj->getnames;     	
    SEQID: foreach my $seqid (keys %{$regions}){
	next SEQID if ! grep {/$seqid$/} @data_file_seqs;
	print STDERR "working on $seqid\n";
	POSTION: foreach my $pos (keys %{$regions->{$seqid}}){
	    $regions->{$seqid}{$pos} = QUERY_RANGE($tabix_obj, $seqid, $pos, $t, $b);
	}
    }
    return $regions;
}

#-----------------------------------------------------------------------------

sub QUERY_RANGE { 

    my ($data, $seqid, $pos, $t, $b) = @_ ;
    my %COUNTS;
    my %NON_REF;
    my %DEPTH;
    my %PI;
    $PI{TARGET} = 0;
    $PI{BETWEEN} = 0;
    my @positions = split /:/, $pos;
    my $iter      = $data->query($seqid, $positions[0], $positions[1]);
    my %gvf_dat;
    
    LINE: while(my $l = $data->read($iter)){
	chomp $l;
	my @l = split /\t/, $l;
	my $line_dat = Parse_Genotype_Lines($l);
        my $t_dat    = Group_Column($line_dat, $t);
        my $b_dat    = Group_Column($line_dat, $b);

	next LINE if $b_dat->[1] < 5;
	next LINE if $t_dat->[1] < 5;

	my $ref_t = $t_dat->[1] - $t_dat->[0];
	my $total_non_ref =  $b_dat->[0] + $t_dat->[0];
	my $total_ref     = ($b_dat->[1] + $t_dat->[1]) - $total_non_ref;
	my $total         = ($b_dat->[1] + $t_dat->[1]);

	$PI{TARGET}  += (choose($t_dat->[0],2)    + choose($ref_t, 2))     / choose($t_dat->[1],2);
	$PI{BETWEEN} += (choose($total_non_ref,2) + choose($total_ref, 2)) / choose($total, 2);
	
	$NON_REF{BACK}{@{$b_dat}[0]}++;
	$NON_REF{TARG}{@{$t_dat}[0]}++;
        $DEPTH{BACK}{@{$b_dat}[1]}++;
	$COUNTS{BACK}{total}++;
	$DEPTH{TARG}{@{$t_dat}[1]}++; 
	$COUNTS{TARG}{total}++;
    }
    
    $PI{TARGET} =  1 - ($PI{TARGET}  / ($positions[1] - $positions[0]));
    $PI{BETWEEN} = 1 - ($PI{BETWEEN} / ($positions[1] - $positions[0]));
    
    my @result_vect = (\%NON_REF, \%DEPTH, \%COUNTS, \%PI);
    
    return \@result_vect;
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
sub SELECTION_REGIONS{

my $regions = shift;

  SEQID: foreach my $seqid (keys %{$regions}){
    POS: foreach my $pos (keys %{$regions->{$seqid}}){
	my $data = $regions->{$seqid}{$pos};
	next POS if $data eq 'NA';
	next POS if ! defined $data;
	
	my $theta_target = THETA_ESTIMATES($data, 'TARG');
	my $theta_backgr = THETA_ESTIMATES($data, 'BACK');
	my $total_sites  = $data->[2]{TARG}{total};
	
	$total_sites = 0 if ! defined $total_sites;
	
	my $md = 'NA';
	$md = MEAN_DIFF($data, $total_sites) if $total_sites = 50;

	my @pos = split /:/, $pos;
	my $S_W = 'NA';
	my $S_F = 'NA';
	my $S_P = 'NA';
	my $S_H = 'NA';
	
# [$Watterson_theta, $Fu_theta, $Pi_theta, $H__theta];

	$S_W = log($theta_target->[0] / $theta_backgr->[0]) if defined $theta_target->[0] && defined $theta_backgr->[0];
	$S_F = log($theta_target->[1] / $theta_backgr->[1]) if defined $theta_target->[1] && defined $theta_backgr->[1];
	$S_P = log($theta_target->[2] / $theta_backgr->[2]) if defined $theta_target->[2] && defined $theta_backgr->[2];
	$S_H = log($theta_target->[3] / $theta_backgr->[3]) if defined $theta_target->[3] && defined $theta_backgr->[3];
	
	my $directional_FST = 1 - ($data->[3]{TARGET} / $data->[3]{BETWEEN});
	my $g_start = $BUILD_RECORD{$seqid} + $pos[0];
	my $g_end   = $BUILD_RECORD{$seqid} + $pos[1];

#	print "$seqid\t$pos[0]\t$pos[1]\t$g_start\t$g_end\t$total_sites\t$S_F\t$S_W\t$S_P\t$S_H\t$directional_FST\t$md\n";
	print "$seqid\t$pos[0]\t$pos[1]\t$g_start\t$g_end\t$total_sites\t";
	print "$theta_target->[0]\t$theta_backgr->[0]\t";
	print "$theta_target->[1]\t$theta_backgr->[1]\t";
	print "$theta_target->[2]\t$theta_backgr->[2]\t";
	print "$theta_target->[3]\t$theta_backgr->[3]\t";
	print "$S_W\t$S_F\t$S_P\t$S_H\t$directional_FST\t$md\n";
    }
  }
}

#-----------------------------------------------------------------------------
sub THETA_ESTIMATES{
    my ($data, $group) = @_;

    return if ! defined $data->[2]{$group}{total};
    my $total_sites = $data->[2]{$group}{total};
    
    return if $total_sites < $seg;
    
    my $average_depth = 0;
    my $max_depth     = 0;
    
    foreach my $j (keys %{$data->[1]{$group}}){
	my $j_i = $data->[1]{$group}{$j};
	$max_depth = $j if $j > $max_depth;
	$average_depth += ($j_i * $j) / $total_sites;
    }
    
    return if $average_depth < $depth;
    
    my $a_n = 0;

    for (my $i = 1; $i < $average_depth; $i++){
	$a_n += 1 / $i;
    }
    my $num          = 0;
    my $fu_sum_theta = 0;
    my $H_sum_theta  = 0;
    my $PI_sum_theta = 0;

    NON_REF: foreach my $i (keys %{$data->[0]{$group}}){
	my $i_i     = $data->[0]{$group}{$i};
	my $Xi      = $i_i / $total_sites;
	my $Theta_i = $i * $Xi;
	
	$H_sum_theta  += $Theta_i * $i;
	$PI_sum_theta += $Theta_i * ($max_depth - $i);
	$num          += $Theta_i * (1 / $i) if $i > 0;
	$fu_sum_theta += $Theta_i;	
    }
    
    return if $num == 0;
    return if $fu_sum_theta == 0;
    
    my $two_over_chr = 2/($average_depth * ($average_depth  - 1));
    
    my $Pi_theta        = $PI_sum_theta * $two_over_chr;
    my $Fu_theta        = (1 / $average_depth) * $fu_sum_theta;
    my $H__theta        = $H_sum_theta * $two_over_chr;
    my $Watterson_theta = $num / $a_n;
    
    return [$Watterson_theta, $Fu_theta, $Pi_theta, $H__theta];
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
sub Parse_Build{
    my $running_position = 0;
    open (my $IN, '<', $build) or die "Can't open $build for reading\n$!\n";
    while (my $l = <$IN>) {
        chomp $l;
        my @l = split /\t/, $l;
	$BUILD_RECORD{$l[0]} = $running_position;
        $running_position += $l[2];
    }
}
#-----------------------------------------------------------------------------                                                          
sub factorial{
    my $n = shift;
    my $s = 1;
    $s *= $n-- while $n > 0;
    return $s
}
#-----------------------------------------------------------------------------                     
sub choose{
    my ($n, $k) = @_;
    return $N_CHOOSE{$n} if defined $N_CHOOSE{$n};
    my $results = Math::BigInt->new();
    my $x    = Math::BigInt->new($n);
    my $y    = Math::BigInt->new($k);
    $results =  $x->Math::BigInt::bnok($y);
    $N_CHOOSE{$n} =  $results->{value}[0];
    return $N_CHOOSE{$n};
}
#-----------------------------------------------------------------------------                                                                         
sub MEAN_DIFF{

    my ($data, $total) = @_;
    my @g = ('TARG', 'BACK');
    
    my $md_sum = 0;

    foreach my $non_ref (keys %{$data->[0]{$g[0]}}){
	my $non_ref_count = $data->[0]{$g[1]}{$non_ref};
	foreach my $non_ref_2 (keys %{$data->[0]{$g[1]}}){
	    $md_sum += abs($non_ref_count - $non_ref_2) * $non_ref_2;
	}
    }
    
    return (1 / ($total * ($total -1))) * $md_sum;
}
