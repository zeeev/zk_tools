#!/usr/bin/perl

package gamma_est;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(camel);
our @EXPORT_OK = qw($weight);
our @VERSION = 1.00;

use Math::GammaFunction qw/:all/;

#my $file = $ARGV[0];
#open(FIN, $file);
#my @a;
#while(my $line = <FIN>){
#	$line =~ s/\s//g;
#	if((length($line) > 0) && ($line>0)){
#		push(@a,$line);
#	}
#}

my $data_len=0;
my $x_sum=0;
my $ln_xsum=0;
my $h = .00001;

sub log_likelihood{
	my $k = shift;
	($k-1)*$ln_xsum - $data_len*$k - $data_len*$k*log($x_sum/($data_len*$k)) - $data_len*log(gamma($k));		
}

sub first_derivative{
	my $k = shift;
	return (log_likelihood($k+$h)-log_likelihood($k-$h))/(2*$h);
}

sub second_derivative{
	my $k = shift;
	return (log_likelihood($k-$h)-2*log_likelihood($k)+log_likelihood($k+$h))/$h**2.0;
}

sub get_cutoff{
	my $k = shift;
	my $t = shift;
	my $per_cut = shift;
	my $x = 0;
	my $step = .001;
	my $cum_area = 0;
	while($cum_area < $per_cut){
		my $tmp_x = $x+.5*$step;
		$cum_area+=$step*(($tmp_x**($k-1))*(exp(-$tmp_x/$t)/(($t**$k)*gamma($k))));
		$x+=$step;
	}
	return $x-$step;
}
	

#note: this module requires a data array: a, and a decimal percent cutoff: percent_cutoff
sub newton_method{
	my $a = shift;
	my $percent_cutoff = shift;
	foreach my $data (@{$a}){
        #foreach my $data (@a){
		if($data > 0){
			$x_sum += $data;
        		$ln_xsum += log($data);
			$data_len++;
		}
	}
	if($data_len == 0){print "All data was zero\n"; die;}
	my $k = 1.3;
	my $i = 0;
	my $diff = 5;
	my $tolerance = .00001;
	while(($diff > $tolerance) && ($i<10)){
		my $tmp_k = $k - first_derivative($k)/second_derivative($k);
		$diff = abs($tmp_k - $k);
		$k = $tmp_k;
		$i+=1;
		#print $k."\n";
	}

	my $t = $x_sum/($k*$data_len);
	my $cutoff = get_cutoff($k,$t,$percent_cutoff);
	return $cutoff	

}



