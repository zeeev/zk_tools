#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Tabix;
use Statistics::Descriptive;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

BIN_VALUE_BY_REGION_W_BUILD -b build.file.txt -d data_file.gz -w window -v value_column -c stats descriptive call


Description:

I bin non overlapping windows.  If a feature is shorter than the window size 
I will just use the whole feature.  use the value flag to tell me what data you want worked on (zero indexed).

";


my ($help);
my $build;
my $t;
my $b;
my $window;
my $data;
my $call; 
my $opt_success = GetOptions('help'      => \$help,
			     'call=s'    => \$call,
			     'data=s'    => \$data,
			     'build=s'   => \$build,
			     'window=s'  => \$window,
			     'target=s'    => \$t,
                             'backround=s' => \$b,

    );



die $usage if $help || ! $opt_success;
die $usage unless $data && $build && $t && $b;

my %PRAGMA;

$t = Parse_group($t);
$b = Parse_group($b);

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

	my $tail = 1;
	my $cur  = 1;
	if($line_contents[2] > $window){
	    while ($tail + $window <= $line_contents[2]){
		my $cur = $tail + $window;
		$REGION_STRUCT{$line_contents[0]}{"$tail:$cur"} = 'NA';
		$tail = $cur;
	    }
	}
	$tail = 0 if $tail == $window + 1;
       
	$REGION_STRUCT{$line_contents[0]}{"$tail:$line_contents[2]"} = 'NA';

    }
    return \%REGION_STRUCT;
}


#-----------------------------------------------------------------------------
sub PARSE_REGIONS {
    
    my ($data_file, $regions, $genomic_column, $value_column) = @_;    
    my $tabix_obj = Tabix->new(-data=>$data_file); 
    my @data_file_seqs = $tabix_obj->getnames;     	
    SEQID: foreach my $seqid (keys %{$regions}){
	print STDERR "working on $seqid\n";
	POSTION: foreach my $pos (keys %{$regions->{$seqid}}){
	    next SEQID if ! grep {/$seqid$/} @data_file_seqs;
	    $regions->{$seqid}{$pos} = QUERY_RANGE($tabix_obj, $seqid, $pos);
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
    my @positions = split /:/, $pos;
    my $iter    = $data->query($seqid, $positions[0], $positions[1]);
    my %gvf_dat;
    
    LINE: while(my $l = $data->read($iter)){
	chomp $l;
	my @l = split /\t/, $l;
	my $line_dat = Parse_Genotype_Lines($l);
        my $t_dat    = Group_Column($line_dat, $t);
        my $b_dat    = Group_Column($line_dat, $b);

        $NON_REF{BACK}{@{$b_dat}[0]}++;
        $DEPTH{BACK}{@{$b_dat}[1]}++;
	$COUNTS{BACK}{total}++;
	$NON_REF{TARG}{@{$t_dat}[1]}++; 
	$DEPTH{TARG}{@{$t_dat}[1]}++; 
	$COUNTS{TARG}{total}++;
    }
    
    my @result_vect = (\%NON_REF, \%DEPTH, \%COUNTS);
    
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
	my $w_theta_target = WATTERSON_THETA($data, 'TARG');
	my $w_theta_target = WATTERSON_THETA($data, 'BACK');
	my $S_W = log($w_theta_target / $w_theta_target);
    }
  }
}

#-----------------------------------------------------------------------------

sub WATTERSON_THETA{
    my ($data, $group) =  @_;

#    $NON_REF{BACK}{@{$b_dat}[0]}++;
#    $DEPTH{BACK}{@{$b_dat}[1]}++;
#    $COUNTS{BACK}{total}++;
#    $NON_REF{TARG}{@{$t_dat}[1]}++;
#    $DEPTH{TARG}{@{$t_dat}[1]}++;
#    $COUNTS{TARG}{total}++;
# (\%NON_REF, \%DEPTH, \%COUNTS);
    
    
    my $total_sites = @{$data}[2]->{$group}{total};
    
    #calculating the theta sum;

    my $theta_i_sum = 0;
    
    foreach my $i (keys %{@{$data}[0]->{$group}}){
	my $n_i = @{$data}[0]->{$group}{$i};
	my $theta_i = $i * ($n_i / $total_sites); 
	$theta_i_sum += (1/$i) * $theta_i;
    }
    
    #calculating the avearge depth

    my $average_depth = 0;
    
    foreach my $j (keys %{@{$data}[1]->{$group}}){
	my $j_i = @{$data}[0]->{$group}{$i}
        $average_depth += ($j_i * $j) / $total_sites;
    }
    
    my $Watterson_theta = (1/$average_depth) * $theta_i_sum;
    return $Watterson_theta;
}


