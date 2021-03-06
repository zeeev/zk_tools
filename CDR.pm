package CDR;
use Tabix;
use strict;
use warnings;
use Data::Dumper;
use PDL;
use PDL::Stats::Basic;
use Math::CDF;
use Set::IntSpan::Fast::XS;
use List::MoreUtils qw(uniq);
use Math::BigFloat;
use Bit::Vector;

#-----------------------------------------------------------------------------
#--------------------------------- Constructor -------------------------------
#-----------------------------------------------------------------------------

=head1 new

     Title   : new
     Usage   : CDR->new();
     Function: Creates a CDR object;
     Returns : A CDR object
     Args    :
=cut

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    $self->_initialize_args(@_);
    $self->_Parse_Pragma;
    if($self->{'file'} !~ /\.gz$/){
	system ("bgzip $self->{'file'}");
	system ("tabix -f -s 1 -b 2 -e 3 $self->{'file'}\.gz");
	$self->{'file'} = "$self->{'file'}\.gz";
    }

    die "no defined pragma" unless defined $self->{'pragma'};
    die "no defined filekeys" unless defined $self->{'file_keys'};
    $self->_Load_Seqids();
    return $self;
}

#-----------------------------------------------------------------------------

sub _initialize_args {
    my $self = shift;

    my $args = $self->prepare_args(@_);
    # Set valid class attributes here
    my %valid_attributes = map {$_ => 1} qw(file file_pragma);
    for my $attribute (keys %{$args}) {
        if (! exists $valid_attributes{$attribute}) {
            my $message = ("$attribute is an invalid attribute and will be" .
                           "skipped");
            $self->warn('invalid_attribute', $message);
        }
        else {
            $self->$attribute($args->{$attribute})
	}
    }
}

#-----------------------------------------------------------------------------
#-------------------------------- Attributes ---------------------------------
#-----------------------------------------------------------------------------

=head1  ATTRIBUTES

All attributes can be supplied as parameters to the constructor as a
list (or referenece) of key value pairs.

=head2 file

  Title   : file
  Usage   : $file = $self->file($cdr_filename);
  Function: Get/set the value of the CDR file (to be) indexed.
  Returns : The CDR file name.
  Args    : The CDR file name.

=cut

#-----------------------------------------------------------------------------

sub file {
    my ($self, $filename) = @_;

    if ($filename) {
        if (! -r $filename) {
            my $message = "$filename is not readable";
            $self->throw('file_not_readable', $message);
        }
        $self->{file} = $filename;
    }
    return $self->{file};
}

#----------------------------------------------------------------------------- 

sub file_pragma {
    my ($self, $filename) = @_;

    if ($filename) {
        if (! -r $filename) {
            my $message = "$filename is not readable";
            $self->throw('file_not_readable', $message);
        }
	$self->{file_pragma} = $filename;
    }
    return $self->{file_pragma};
}

#-----------------------------------------------------------------------------   

=head2 sparse

  Title   : sparse
  Usage   : $sparse = $self->sparse(1);
  Function: Get/set the value of sparse.  A true value for sparse will
            cause all features to be returned as a 9 element array ref
            corresponding to the columns in GFF3.  A false value for sparse
            (the default) will cause features to be returned as hash refs.
  Returns : 1 or 0
  Args    : 1 or 0

=cut

sub sparse {
    my ($self, $value) = @_;

    if (defined $value) {
        $self->{sparse} = $value ? 1 : 0;
    }
    return $self->{sparse};
}

#-----------------------------------------------------------------------------

=head2 prepare_args

 Title   : prepare_args
 Usage   : $args = $self->prepare_args(@_);
 Function: Take a list of key value pairs that may be an array, hash or ref
           to either and return them as a hash or hash reference depending on
           calling context.
 Returns : Hash or hash reference
 Args    : An array, hash or reference to either

=cut

sub prepare_args {

    my ($self, @args) = @_;

    my %args_hash;

    if (scalar @args == 1 && ref $args[0] eq 'ARRAY') {
	%args_hash = @{$args[0]};
    }
    elsif (scalar @args == 1 && ref $args[0] eq 'HASH') {
	%args_hash = %{$args[0]};
    }
    elsif (scalar @args % 2 == 0) {
	%args_hash = @args;
    }
    else {
	my $class = ref($self);
	my $err_code = 'invalid_arguments_to_prepare_args';
                my $err_msg  = ("Bad arguments passed to $class. A list "   .
                                "of key value pairs or a reference to "     .
                                "such a list was expected, But we got:\n"   .
                                join ' ', @args);
                $self->throw(message => $err_msg,
                             code    => $err_code);
    }

    return wantarray ? %args_hash : \%args_hash;
}

#-----------------------------------------------------------------------------

sub throw {
    shift->handle_message('FATAL', @_);
}

#-----------------------------------------------------------------------------    

=head2 handle_message                                                             
                                                                                  
 Title   : handle_message                                                         
 Usage   : $base->handle_message($level, $code, $message);                        
 Function: Handle a message and report to the user appropriately.                 
 Returns : None                                                                   
 Args    : level:   FATAL, WARN, INFO                                             
           code:    $info_code # single_word_code_for_info                        
           message: $info_msg  # Free text description of info                    
                                                                                  
=cut                                                                              

sub handle_message {
    my ($self, $level, $code, $message) = @_;

    my $caller = ref($self);
    $level ||= 'UNKNOWN';
    $message ||= $caller;
    if (! $code) {
	$code = 'unspecified_code';
          $message = join '', ("VAAST::Base is handeling a message " .
                                "for $caller without an error code.  "  .
			       "Complain to the author of $caller!");
    }
    chomp $message;
    $message .= "\n";
 
   if ($level eq 'FATAL') {
	$message = join ' : ', ('FATAL', $code, $message);
	croak $message;
    }
    elsif ($level eq 'WARN') {
	$message = join ' : ', ('WARN', $code, $message);
	print STDERR $message;
    }
    elsif ($level eq 'INFO') {
	$message = join ' : ', ('INFO', $code, $message);
	print STDERR $message;
    }
    else {
          $message = join '', ("VAAST::Base is handeling a message " .
                               "for $caller without an error level.  "  .
                               "Complain to the author of $caller!\n");
          chomp $message;
          $message = join ' : ', ('UNKNOWN', $code, $message);
          croak $message;
    }
}
#-----------------------------------------------------------------------------   

sub _Load_Seqids{
    my $self   = shift;
    my $t      = Tabix->new(-data => $self->{'file'});
    my @seqids = $t->getnames(); 
    $self->{'seqids'} = \@seqids;
}

#-----------------------------------------------------------------------------   

sub _Parse_Pragma{  
    my $self = shift;    
	my @pragma;
	my %file_key;
	my $file = $self->{file_pragma};
	open(FH, '<', $file) || die "cannot open $file\n";
	while(my $l = <FH>){
	    chomp $l;
	    next unless $l =~ /##/;
   	    push @pragma, $l;
	    if($l =~ /FILE-INDEX/){
		my @l = split /\t/, $l;
		$file_key{$l[2]} = $l[3];
	    }
	}
	close FH;
	$self->{'pragma'} = \@pragma;
	$self->{'file_keys'} = \%file_key;
}

#-----------------------------------------------------------------------------   

sub Get_Pragma { 
    my $self = shift;      
    my @pragma = @{$self->{'pragma'}};
    return \@pragma;
}

#-----------------------------------------------------------------------------   

sub Get_Seqids {
    my $self = shift;
    return $self->{'seqids'}; 
}

#-----------------------------------------------------------------------------   

sub _Parse_Line{
    my $self = shift;
    $self->_Parse_Raw_Line;
    $self->_Parse_Reference_Genotypes;
    $self->_Parse_Genotype;
}

#-----------------------------------------------------------------------------   
sub start{
    my $self = shift;
    $self->_Parse_Raw_Line if ! $self->{'line'}{'refined'}{'start'};
    return $self->{'line'}{'refined'}{'start'};
}
#-----------------------------------------------------------------------------   
sub end{
    my $self = shift;
    $self->_Parse_Raw_Line if ! $self->{'line'}{'refined'}{'start'};
    return $self->{'line'}{'refined'}{'start'};
}
#-----------------------------------------------------------------------------     
sub _type{
    my $self = shift;
    $self->_Parse_Raw_Line if ! $self->{'line'}{'refined'}{'type'};
    return $self->{'line'}{'refined'}{'type'};
}
#-----------------------------------------------------------------------------     

sub _Parse_Raw_Line{

    my $self = shift;
    my %lc;
    my @l = split /\t/, $self->{'line'}{'raw'};
    $lc{'seqid'}  = shift @l;
    $lc{'start'}  = shift @l;
    $lc{'end'}    = shift @l;
    $lc{'type'}   = shift @l;
    $lc{'effect'} = shift @l;
    $lc{'ref'}    = _Parse_Ref(shift @l);
    $self->{'line'}{'refined'} = \%lc;
    $self->{'line'}{'raw_genotypes'} = \@l;

}	       

#-----------------------------------------------------------------------------   

sub _Parse_Ref{
    my $ref = shift;
    my @ref;
    if ($ref =~ /\|/){
	push @ref, split /\|/, $ref;
    }
    esle{
	push @ref, $ref;
    }
    return \@ref;
}

#-----------------------------------------------------------------------------   

sub _Parse_Genotype{

    my $self = shift;
    my %genotype;

    for my $g (@{$self->{'line'}{'raw_genotypes'}}){
	my @g = split /\|/, $g; 

	my @indvs = @{get_tag_from_string($g[0])};

	INDV: foreach my $i (@indvs){
#	    die "overlapping genotypes : $self->{'line'}{'raw'}\n" if exists $genotype{$i}{'genotype'};
	    $self->{line}{genotypes}{$i}{'genotype'} = $g[1];
	    $self->{line}{genotypes}{$i}{'effect'}   = $g[2];
	}
    }
}

#-----------------------------------------------------------------------------   
sub _Parse_Reference_Genotypes{
    
    my $self  = shift;
    if (defined $self->{'line'}{'refined'}{'ref'}[1]){
	
	while(my($key, $file) = each %{$self->{'file_keys'}}){
	    $self->{'line'}{'genotypes'}{$key}{'genotype'}  = "$self->{'line'}{'refined'}{'ref'}[0]:$self->{'line'}{'refined'}{'ref'}[0]";
	    $self->{'line'}{'genotypes'}{$key}{'effect'}  = "$self->{'line'}{'refined'}{'ref'}[1]:$self->{'line'}{'refined'}{'ref'}[1]";	
	}
    }
    else{
	while(my($key, $file) = each %{$self->{'file_keys'}}){
            $self->{'line'}{'genotypes'}{$key}{'genotype'}  = "$self->{'line'}{'refined'}{'ref'}[0]:$self->{'line'}{'refined'}{'ref'}[0]";
            $self->{'line'}{'genotypes'}{$key}{'effect'}  = undef;
	}
    }
}

    
#-----------------------------------------------------------------------------   
sub _Parse_Alleles{

    my $genotypes = shift;
    my %genotype_dat;
    
    while(my ($indv, $g) = each %{$genotypes}){
	map {$genotype_dat{$_}++} split /:/, $g->{genotype};
    }
    
    my @alleles = keys %genotype_dat;
    return \@alleles;
}
#-----------------------------------------------------------------------------   

sub _Count_Genotypes{

    my $info_hash = shift;
  
    my %allele_counts;
    my %genotype_counts;
    
    $allele_counts{'^'}     = 0;
    $genotype_counts{'^:^'} = 0;

    while( my($info, $value) = each %{$info_hash}){
	while(my($key, $value_2) = each %{$value}){
	    if ($key =~ /genotype/){
		my @gen = split /:/, $value_2;
		map {$allele_counts{$_}++} @gen; 
		$genotype_counts{$value_2}++;
	    }
	}
    }
    
    my %total;

    $total{'allele_counts'}{'nocall'} = 0;
    $total{'allele_counts'}{'called'} = 0;
    $total{'genotype_counts'}{'called'} = 0;
    $total{'genotype_counts'}{'nocall'} = 0;
    
    while( my ($key, $value) = each %genotype_counts){
	$total{'genotype_counts'}{'nocall'} += $value if $key =~ /\^/;
	$total{'genotype_counts'}{'called'} += $value if $key !~ /\^/;
    }

    while( my ($key, $value) = each %allele_counts){
	$total{'allele_counts'}{'nocall'} += $value if $key =~ /\^/;
	$total{'allele_counts'}{'called'} += $value if $key !~ /\^/;
    }

    return (\%total, \%allele_counts, \%genotype_counts);
}

#-----------------------------------------------------------------------------   

sub HWE_Departure{
    my ($self, $chr) = @_;
    my $tabix = Tabix->new(-data => $self->{'file'});
    my $it = $tabix->query($chr);
    return if ! defined $it;

  LINE: while(my $l = $tabix->read($it)){
      $self->{line}{raw} = $l;
      $self->_Parse_Line();
      my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
      my @counts  = _Count_Genotypes($self->{'line'}{'genotypes'});
      my $ref = @{$self->{line}{refined}{ref}}[0];
      my ($alt) = grep {!/\^|$ref/} @alleles;
      return if scalar (grep {!/\^|$ref/} @alleles) != 1;
      
      my $total_alleles   = $counts[0]->{'allele_counts'}{'called'};
      my $total_genotypes = $counts[0]->{'genotype_counts'}{'called'};
      
      my ($p, $q) = (0,0);
	      
      $p = $counts[1]->{$ref} / $total_alleles if defined $counts[1]->{$ref} && $counts[1]->{$ref} > 0;
      $q = $counts[1]->{$alt} / $total_alleles if defined $counts[1]->{$alt} && $counts[1]->{$alt} > 0;
                  
      my $exp_pp = int ($p**2  * $total_genotypes);
      my $exp_qq = int ($q**2  * $total_genotypes);
      my $exp_pq = int (2 * $p * $q * $total_genotypes); 
      
      my ($AA, $AB, $BB, $AN, $BN) = (0,0,0,0,0);
      
      $AB = $counts[2]->{"$alleles[0]:$alleles[1]"} if defined $counts[2]->{"$alleles[0]:$alleles[1]"};
      $AA = $counts[2]->{"$ref:$ref"} if defined $counts[2]->{"$ref:$ref"};
      $BB = $counts[2]->{"$alt:$alt"} if defined $counts[2]->{"$alt:$alt"};
      $AN = $counts[2]->{"$ref:\^"} if defined $counts[2]->{"$ref:\^"};
      $BN = $counts[2]->{"$alt:\^"} if defined $counts[2]->{"$alt:\^"};
      
      my $p_pp =  _Chi_Lookup($exp_pp, $AA);
      my $p_pq =  _Chi_Lookup($exp_pq, $AB);
      my $p_qq =  _Chi_Lookup($exp_qq, $BB);
      
 
      my $chi = $p_pp + $p_pq + $p_qq;
      my $p_value   =  1 - Math::CDF::pchisq($chi, 1);  
    
      #print "$self->{line}{refined}{start}\t$p_pp\t$AA\t$AB\t$BB\t$exp_pp\t$exp_pq\t$exp_qq\t$p\t$q\t$chi\t$p_value\n";
      print "$self->{line}{refined}{start}\t$AA\t$AB\t$BB\t$AN\t$BN\n";
  } 
}
#-----------------------------------------------------------------------------   

sub _Chi_Lookup{

    my ($exp, $obs) = @_;
    my $chi;
   
    if ($exp != 0 ){
	$chi = (($obs - $exp)**2) / $exp;
    }
    else{
	$chi = 0;
    }
    return $chi;
}
    
#-----------------------------------------------------------------------------   
sub Weir_Cockerham{

#the equation was taken from weir and cockerham 1984 equation (1) pg 1359
#right now it only allows for two subpopulations since I am only doing fst
#genome scans for two populations.  All subroutines used for this calculation
#begin with wc_.  This has been mimiced after Pegas [R]'s equations.  The pigeon
#headcrest snp has validated this stie. 
    
#library(pegas)
#my.dat<-as.data.frame(cbind(c(rep("T/T", 8), rep("T/C",2), rep("C/C", 23)),c(rep(1, 8), rep(2,2), rep(2, 23))))
#colnames(my.dat)<-c("genotypes", "population")
#my.loci<-as.loci(my.dat)
    
    my ($self, $groups, $scaff) = @_;
    my $t = Tabix->new(-data => $self->{'file'});
    my $it = $t->query($scaff);
    next SCAFF if ! defined $it;
    
    my @DAT;

  LINE: while(my $l = $t->read($it)){
        $self->{line}{raw} = $l;
        $self->_Parse_Line();
	$self->_Group($groups);
	my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
	my @gens = grep {!/\^/} @alleles;
	next LINE if (scalar grep {!/\^/} @alleles) != 2;
	my @tot          = _Count_Genotypes($self->{'line'}{'genotypes'});
	next LINE if $tot[0]->{allele_counts}{nocall} > $tot[0]->{allele_counts}{called};
	my @g_a          = _Count_Genotypes($self->{'line'}{'group_A'});
	print STDERR Dumper(@g_a) if ! defined $g_a[0]->{allele_counts}{called};
	next LINE if $g_a[0]->{allele_counts}{called} == 0; 	
	my @g_b          = _Count_Genotypes($self->{'line'}{'group_B'});
	next LINE if $g_b[0]->{allele_counts}{called} == 0;
	my ($n_bar, $r)  = wc_n_bar([\@g_a, \@g_b]);

	my $N = $n_bar*$r;

	next LINE if $N == 0;

	my $n_c          = wc_n_c($N, $r, [\@g_a, \@g_b]);
	my $p_bar        = wc_p_bar($N, \@gens, [\@g_a, \@g_b]);
	my $s_2          = wc_s_2($n_bar, $r, $p_bar, \@gens,[\@g_a, \@g_b]);
	my $h_bar        = wc_h_bar($N, \@gens, [\@g_a, \@g_b]);
	my $a            = wc_a($n_bar, $r, $n_c, @{$p_bar}[0], @{$s_2}[0], $h_bar);
	my $b            = wc_b($n_bar, @{$p_bar}[0], $r, @{$s_2}[0], $h_bar);
	my $c            = wc_c($h_bar);
	next LINE if $a == 0 || $c == 0;
	my $fit =  1 - ($c/($a + $b + $c));
	my $fst =  $a/($a + $b + $c);
	my $fis =  1 - ($c/($b + $c));
	push @DAT,  "$self->{line}{refined}{seqid}\t$self->{line}{refined}{start}\t$fit\t$fst\t$fis";
    }
    return \@DAT;
  }

#-----------------------------------------------------------------------------   
#This may be a wonky way to do it.  But it allows for multiple groups.

sub wc_n_bar{
    my $grouped_counts = shift;
    my $r = scalar @{$grouped_counts};
    my $n_bar = 0;
    foreach my $g (@{$grouped_counts}){
	my $n = @{$g}[0]->{'genotype_counts'}{'called'};
	$n_bar += $n / $r;
    }
    return $n_bar, $r;
}
#-----------------------------------------------------------------------------   
sub wc_n_c{
    my ($N, $r, $grouped_counts) = @_;
    my $sum_exp = 0;
    foreach my $g (@{$grouped_counts}){
	$sum_exp += (@{$g}[0]->{'genotype_counts'}{'called'}**2 / $N);
    }
    my $nc =  ($N - $sum_exp)/($r - 1);
    return $nc;
}
#-----------------------------------------------------------------------------   
sub wc_p_bar{
    
    #you are summing across alleles NOT groups FML.
    
    my ($N, $alleles, $grouped_counts) = @_;
    my @p_bars = ();
    foreach my $a (@{$alleles}){
	my $p_bar = 0;
	foreach my $g (@{$grouped_counts}){
	    if(defined @{$g}[1]->{$a}){
		my $ptild = @{$g}[1]->{$a} / @{$g}[0]->{allele_counts}{called};
		$p_bar   += @{$g}[0]->{genotype_counts}{called} * $ptild / $N;  
	    }
	}				      
	push @p_bars, $p_bar;
    }
    return \@p_bars;
}
#-----------------------------------------------------------------------------   
sub wc_s_2{
    my ($n_bar, $r, $p_bar, $alleles, $grouped_counts) = @_;
    my @s_2s = ();
    my $dem = ($r - 1) * $n_bar;
    
    my $count=0;
    foreach my $a (@{$alleles}){
	my $s_2 = 0;
	foreach my $g (@{$grouped_counts}){
	    my $n = defined @{$g}[0]->{'allele_counts'}{'called'} ? @{$g}[0]->{'allele_counts'}{'called'} : 0;
	    my $p_tild = defined @{$g}[1]->{$a} ?  @{$g}[1]->{$a} / $n : 0;
	    $s_2 += $n != 0 ? ((@{$g}[0]->{'genotype_counts'}{'called'} *($p_tild - @{$p_bar}[$count])**2) / $dem) : 0;
	}
	$count++;
	push @s_2s, $s_2;
    }
    return \@s_2s;
}
#-----------------------------------------------------------------------------   
sub wc_h_bar{
    my ($N, $alleles, $grouped_counts) = @_;
    my $h_bar = 0;
    my $het = join ":", sort {$a cmp $b} @{$alleles};
    
    foreach my $g (@{$grouped_counts}){
	$h_bar += defined @{$g}[2]->{$het} ?  @{$g}[2]->{$het}  / $N : 0;
    }
    return $h_bar;
}
#-----------------------------------------------------------------------------   
sub wc_a{
    my ($n_bar, $r, $n_c, $p_bar, $s_2, $h_bar) = @_;
    my $a = 0; 
    my $A = $p_bar * (1 - $p_bar) - ($r - 1) * ($s_2/$r);
    $a = $n_bar * ($s_2 - ($A - $h_bar/4)/($n_bar - 1))/$n_c;
    return $a;
}
#-----------------------------------------------------------------------------   
sub wc_b{
    my ($n_bar, $p_bar, $r, $s_sq, $h_bar) = @_;
    my $b = 0;
    my $inner = $p_bar * (1 - $p_bar) - (($r - 1)/$r) * $s_sq - (2 * $n_bar - 1)/(4 * $n_bar) * $h_bar;
    $b = ($n_bar / ($n_bar - 1)) * $inner;
    return $b;
}
#-----------------------------------------------------------------------------   
sub wc_c{
    my $h_bar= shift;
    my $c = 0;
    $c = 1/2*$h_bar;
    }
#-----------------------------------------------------------------------------   
sub _Group{

    #Groups the loci into two groups.  You can use ! to exclude an indvidual from
    #both groups.  If !-10 for example 10 will be removed from all the genotypes
    #including $self->{line}{genotypes}.  Total exclusion. 
    

    my %groupA;
    my %group_disgard;

    my %groupA_DAT;
    my %groupB_DAT;

    my ($self, $groups) = @_;
    
    foreach my $i (@{$groups}){
	if($i !~ /!/){
	    $groupA{$i} = 1; 
	}
	else{
	    my @is = split /-/, $i;
	    $group_disgard{$is[1]} = 1;
	}
    }
    
    GROUP: while(my($indv, $indv_hash) = each %{$self->{'line'}{'genotypes'}}){
	if (exists $groupA{$indv}){
	    $groupA_DAT{$indv} = $indv_hash;
	}
	elsif(exists $group_disgard{$indv}){
	    next GROUP;
	}
	else{
	    $groupB_DAT{$indv} = $indv_hash;
	}
    }

    foreach my $key (keys %group_disgard){
	delete $self->{line}{genotypes}{$key};
    }
    
    $self->{line}{group_A} = \%groupA_DAT;
    $self->{line}{group_B} = \%groupB_DAT; 
}

#-----------------------------------------------------------------------------   

sub _Scrub_No_Call{

    my $self    = shift;
    my $alleles = _Parse_Alleles($self->{'line'}{'genotypes'});
    return 0 if scalar @{$alleles} == 1 && grep {/\^/} @{$alleles};
    return 1;
}
	
#-----------------------------------------------------------------------------   

sub _Print_Line{
    my $self = shift;
    print "$self->{'line'}{'raw'}\n";

}

#-----------------------------------------------------------------------------   

sub _Print_Pragma{
    my $self = shift;   
    my @pragma = @{$self->{'pragma'}};
    foreach my $p (@pragma){
	print "$p\n";
    }
}


#-----------------------------------------------------------------------------   
sub _sort_by_increasing_vals{
    
    my $hashref = shift;
    my %hash = %{$hashref};
    my @s_keys = sort {$hash{$a} <=> $hash{$b}} keys %hash; 
    
    my $lval;
    RID: while($lval = shift @s_keys){
	last RID if $lval ne '^' 
    }
    return $lval;

}

#-----------------------------------------------------------------------------   
sub _Print_MAP{

    my $self = shift;
    return if $self->{line}{refined}{type} ne 'SNV';

    print "$self->{line}{refined}{seqid}\t";
    print "$self->{line}{refined}{seqid}:";
    print "$self->{line}{refined}{start}:";
    print "@{$self->{line}{refined}{ref}}[0]\t";
    print "0\t";
    print "$self->{line}{refined}{start}\n";

}

#-----------------------------------------------------------------------------   

sub _Print_Beagle{

    my $self = shift;


    return if $self->{line}{refined}{type} ne 'SNV';
    my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
    my $as       = join ",", @alleles;

    print "M\t";
    print "$self->{line}{refined}{seqid}:";
    print "$self->{line}{refined}{start}:";
    print "@{$self->{line}{refined}{ref}}[0]:";
    print "$as\t";

    my ($flag) = sort {$b <=> $a} keys %{$self->{line}{genotypes}};
    
    foreach my $key (sort {$a <=> $b} keys %{$self->{line}{genotypes}}){
	my @gens = split /:/, $self->{line}{genotypes}{$key}{genotype};
	print "$gens[0]\t";
	print "$gens[1]";
	$key == $flag ? print "\n" : print "\t";
    }
	
}

#-----------------------------------------------------------------------------   
sub Print_Beagle_Header{

    my $self = shift;
    print "I\t";
    print "id\t";
    
    my ($flag) = sort {$b <=> $a} keys %{$self->{file_keys}};
    foreach my $key (sort {$a <=> $b} keys %{$self->{file_keys}}){
	print "$key\t";
	print "$key";
	$key == $flag ? print "\n" : print "\t";
    }
}

#-----------------------------------------------------------------------------   
sub LD{
    my ($self, $indvs, $args) = @_;
    my $count = 0; 
    my %markers;
    my @pdl;
    
    return if ! defined $args;
    
    my $tabix = Tabix->new(-data => $self->{'file'});
    my $it = $tabix->query($args);
    
  
    return if ! defined $it;
    
    my $n_indv = scalar @{$indvs};
    my $n_alleles = 2 * $n_indv;
    
  LINE: while(my $l = $tabix->read($it)){
      $self->{line}{raw} = $l;
      $self->_Parse_Line();
      
      next LINE unless $self->{line}{refined}{type} eq 'SNV';
      
      my $ref = $self->{line}{refined}{ref}[0];
      $self->_Group($indvs);
      my @counts  = _Count_Genotypes($self->{line}{group_A});
      my @counts2 = _Count_Genotypes($self->{line}{group_B});

      my @alleles = sort {$a cmp $b} keys %{$counts[1]};
      
      next LINE if scalar (grep {!/\^|$ref/} @alleles) != 1;
      next LINE if $counts[0]->{genotype_counts}{nocall}  > 0;
      next LINE if $counts2[0]->{genotype_counts}{nocall}  > 0;    
      next LINE if ($counts[1]->{$alleles[0]} / $n_alleles) > 0.70 || ($counts[1]->{$alleles[0]} / $n_alleles) < 0.30;
      $markers{$count} = $self->{line}{refined}{start};
      
      my @vals;
      
# AA = 1
# AB = 2 
# BB = 3
      
      my @non_ref = grep {!/$ref|\^/} @alleles;
      
      my $hets = "$alleles[0]:$alleles[1]";
      my $homr = "$ref:$ref";
      my $homn = "$non_ref[0]:$non_ref[0]";

    GENOTYES: foreach my $key (@{$indvs}){
	push @vals, 1 if $self->{line}{genotypes}{$key}{genotype} eq $homr;
	push @vals, 2 if $self->{line}{genotypes}{$key}{genotype} eq $hets;
	push @vals, 3 if $self->{line}{genotypes}{$key}{genotype} eq $homn;
    }
      
      $pdl[$count++] = pdl @vals;    
  }
    
# A PDL data stucture. 
    
    return if scalar @pdl <= 1; 
    my $piddle = pdl(@pdl);
    
    my %r_data;
    my $n_row =  $piddle->slice('0,:')->nelem;
    
  OUTER: for (my $i = 0; $i < $n_row; $i++){
      my $loc_a = $piddle->PDL::slice(":,$i");
    INNER: for (my $j = 0; $j < $n_row; $j++){
	next INNER if defined $r_data{$i}{$j};
	next OUTER if abs($markers{$i} - $markers{$j}) >= 1_000_000;
	my $loc_b = $piddle->PDL::slice(":,$j");
	my $cor =  $loc_a->corr($loc_b)->at(0);
	$r_data{$i}{$j} = $cor**2;
    }
  }
    
    my @print_dat;
    
    while(my($loc1, $loc2_dat) = each %r_data){
	my $pos1 = $markers{$loc1};
	while(my ($loc2, $cov) = each %{$loc2_dat}){
	    my $pos2 = $markers{$loc2};
	    my $p_diff = abs($pos1 - $pos2);
	    push @print_dat, "$args\t$pos1\t$pos2\t$p_diff\t$cov"; 
	}
    }
    return @print_dat;
}
#-----------------------------------------------------------------------------   
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

#-----------------------------------------------------------------------------   
sub Distance_Matrix{
    my ($self, $indvs, $scaffs) = @_;
    my $count = 0;
    my %markers;
    my @pdl;

    print join "\t", @{$indvs};
    print "\n";

  SCAFF: foreach my $scaff (@{$scaffs}){

      my $t = Tabix->new(-data => $self->{'file'});
      my $it = $t->query($scaff);
      
      next SCAFF if ! defined $it;
      
      my $n_indv = scalar @{$indvs};
      my $n_alleles = 2 * $n_indv;
      
    LINE: while(my $l = $t->read($it)){
	$self->{line}{raw} = $l;
	$self->_Parse_Line();
	my $ref = $self->{line}{refined}{ref}[0];
	my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
	next LINE unless $self->{line}{refined}{type} eq 'SNV';
	next LINE if scalar (grep {!/\^|$ref/} @alleles) != 1;
	$self->_Group($indvs);
	my @counts  = _Count_Genotypes($self->{line}{group_A});
	next LINE if $counts[0]->{genotype_counts}{nocall} > 0; 
	$markers{$count} = $self->{line}{refined}{start};
	my @vals;
	
      GENOTYES: foreach my $key (@{$indvs}){
	  my @genotype = split /:/, $self->{line}{genotypes}{$key}{genotype};
	  if($genotype[0] ne $ref || $genotype[1] ne $ref){
	      push @vals, '1';
	  }
	  else{
	      push @vals, '0';
	  }
      }
	#  $pdl[$count++] = pdl @vals;
	print join "\t", @vals;
	print "\n";
#ending all loops to load data      
    }
  }
# A PDL data stucture.           

#    #return if there aren't sufficent sites.
#    return if scalar @pdl <= 1000;
#    my $piddle = pdl(@pdl);
#    
#    my @dist;
#    my @distance_matrix;
#    
#    my $n_col =  $piddle->slice(':,0')->nelem;
#    
#  OUTER: for (my $i = 0; $i < $n_col; $i++){
#      my $col_a = $piddle->PDL::slice("$i,:");
#    INNER: for (my $j = 0; $j < $n_col; $j++){
#        next INNER if defined $dist[$i][$j];
#        my $col_b = $piddle->PDL::slice("$j,:");  
#	my $sum = PDL::sum($col_a - $col_b);
#        $dist[$i][$j] = $sum;
#    }
#  }
#    return \@dist
}



#-----------------------------------------------------------------------------   
sub Print_PLINK{

    my ($self, $indvs, $fam, $pre, $features) = @_;
    my $count = 0;
    my %PED;
    my %MAP;
    my @pdl;

    my $count2 = 0;

    foreach my $fkey (@{$indvs}){
	$PED{$fkey} = [$fam, $self->{file_keys}{$fkey}, '0', '0', 'unknown', '0'];
    }


    my $n_indv = scalar @{$indvs};
    my $n_alleles = 2 * $n_indv;

    my $t = Tabix->new(-data => $self->{'file'});
    foreach my $f (@{$features}){
	print STDERR "INFO working on: $f\n";
	my $it = $t->query($f);
	
	#create data structure PLINK format
	
	my $lpos = 0;	
	
      LINE: while(my $l = $t->read($it)){
	  $self->{line}{raw} = $l;
	  $self->_Parse_Line();
	  next LINE if $self->{line}{refined}{seqid} =~ /^C/;  #contigs 
	  next LINE unless $self->{line}{refined}{type} eq 'SNV';
	  if($self->{line}{refined}{start} - $lpos > 100_000){
	      my $ref = $self->{line}{refined}{ref}[0];
	      $self->_Group($indvs);
	      my @counts  = _Count_Genotypes($self->{line}{group_A});
	      next LINE if $counts[0]->{genotype_counts}{nocall} > 0;
	      next LINE if scalar  (grep {!/\^|$ref/} keys %{$counts[1]}) != 1;
	      $lpos = $self->{line}{refined}{start};
	      $MAP{$count2} = [$self->{line}{refined}{seqid}, "$self->{line}{refined}{seqid}:$self->{line}{refined}{start}:$self->{line}{refined}{start}", '0', $self->{line}{refined}{start}];
	      $count2++;
	      foreach my $indv (@{$indvs}){
		  push @{$PED{$indv}}, (split /:/, $self->{'line'}{'genotypes'}{$indv}{genotype});
	      }
	      $lpos = $self->{line}{refined}{start};
	  }
      }
    }
    
#PRINT PED
    open(FH, '>', "$pre.ped") || die "cannot open $pre.per for writing";
    
    foreach my $key (sort {$a <=> $b} keys %PED){
	foreach my $dat (@{$PED{$key}}){
	    print FH "$dat\t";
	}
	print FH "\n";
    }
    
    close(FH);
    
#PRINT MAP
    
    open(FH2, '>', "$pre.map") || die "cannot open $pre.map for writing";
    
    foreach my $m (sort {$a <=> $b} keys %MAP){
	foreach my $dat (@{$MAP{$m}}){
	    print FH2 "$dat\t";
	}
	print FH2 "\n";
    }
    close(FH2);
}


#-----------------------------------------------------------------------------   

#This sub is designed to generate counts of nocall by groups.  Lets say you 
#want to know if your groupsize is 10 how many loci contain at least one NC.
#prints:

#seqid group_size Y_nc_count N_nc_count 

sub cout_nocall_by_group_size{
    my ($self, $indvs, $args) = @_;
    my $tabix = Tabix->new(-data => $self->{'file'});
    my $it = $tabix->query($args);
    return if ! defined $it;
  LINE: while(my $l = $tabix->read($it)){
      $self->{line}{raw} = $l;
      $self->_Parse_Line();
      for(my $i = 1; $i < scalar(@{$indvs}); $i++){
	  my @tmp_indvs = @{$indvs};
	  fisher_yates_shuffle(\@tmp_indvs);
	  my @temp_group;
	  for(my $j = 1; $j <= $i; $j++){
	      push @temp_group, shift @tmp_indvs;
	  }
	  
	  $self->_Group(\@temp_group);
	  my @counts  = _Count_Genotypes($self->{'line'}{'group_A'});
	  my $res = $counts[2]->{'^:^'} > 0 ? 1 : 0;  
	  print "$i\t$res\n";
      }
  }      
}

#-----------------------------------------------------------------------------   

sub Gillespie_Weighted_FST{
    
    my ($self, $groups, $scaffs) = @_;
    my $t = Tabix->new(-data => $self->{'file'});
    my $it = $t->query($scaffs);

    my @dat;
    
  LINE: while(my $l = $t->read($it)){
	$self->{line}{raw} = $l;
	$self->_Parse_Line();
	$self->_Group($groups);
	my @alleles = grep {!/\^/} @{_Parse_Alleles($self->{'line'}{'genotypes'})};
	next LINE if (scalar @alleles) != 2;
	
	my @tot          = _Count_Genotypes($self->{'line'}{'genotypes'});
	my @g_a          = _Count_Genotypes($self->{'line'}{'group_A'});
	my @g_b          = _Count_Genotypes($self->{'line'}{'group_B'});
		
	my $group_sums = sum_heterozygosity(\@tot, \@g_a, \@g_b, $alleles[0], 1);
	my $total_sum  = sum_heterozygosity(\@tot, \@tot, $alleles[0], 0);
	my $fst = ($total_sum - $group_sums) / $total_sum;
	push @dat,"$self->{line}{refined}{seqid}\t$self->{line}{refined}{start}\t$fst";
    }       
    return @dat;
}


#-----------------------------------------------------------------------------   
sub sum_heterozygosity{

    my @groups = @_;
    my $tot = shift @groups;
    my $flag = pop @groups;
    my $allele = pop @groups;
    my $weight = 1;
   
    my $sum = 0; 
    
    foreach my $g (@groups){
	$weight = defined @{$g}[0]->{allele_counts}{called} ? @{$g}[0]->{allele_counts}{called} / @{$tot}[0]->{allele_counts}{called} : 0 if $flag == 1; 
	my $p = defined @{$g}[1]->{$allele} ?  @{$g}[1]->{$allele} / @{$g}[0]->{allele_counts}{called} : 0;
	my $q = (1 - $p);
	$sum += $p > 0 && $q > 0 ? $weight*2*$p*$q : 0;
    }
    return $sum;
} 

#-----------------------------------------------------------------------------   

#This this script takes the entire cdr and loads up all indviduals as a single
#string that contains binary data for rapid intersections to measure deminishing
#returns.  This is some very heavy lifting.

sub RETURN_BIT_HASH{
    
    my ($self, $indvs, $scaffs) = @_;
        
    $self->{bit}{SNV}{$_} = ()       for @{$indvs};
    $self->{bit}{insertion}{$_} = () for @{$indvs};
    $self->{bit}{deletion}{$_} = ()  for @{$indvs};


    my $t = Tabix->new(-data => $self->{'file'});
    
  SCAFF: foreach my $f (@{$scaffs}){
      print STDERR "INFO working on: $f\n";
      my $it = $t->query($f);
      
    LINE: while(my $l = $t->read($it)){
	$self->{line}{raw} = $l;
	$self->_Parse_Line();
	$self->_Group($indvs);
	my @alleles =  grep {!/@{$self->{line}{refined}{ref}}[0]|\^/} @{_Parse_Alleles($self->{line}{group_A})};
	next LINE if ! defined $alleles[0];
	my $type = $self->{line}{refined}{type};

	my %lookup;
	$lookup{$_} = 0 for @alleles;
	
      OUTER: while( my($key, $value) = each %{$self->{line}{group_A}}){
	  my %cp_lookup = %lookup;
	  $cp_lookup{$_} = 1 for split /:/, $self->{line}{group_A}{$key}{genotype};
	  
	ALLELE: foreach my $a (@alleles){
	    push @{$self->{bit}{$type}{$key}}, $cp_lookup{$a};
	}
      }
    }
  }
    print STDERR "INFO building vectors\n";
    my $vect_data = $self->vectorize_data();
    return $vect_data;     
}
#-----------------------------------------------------------------------------   
#This call Bit::Vector, vectorizes the data and returns the data structure

sub vectorize_data{

    my $self = shift;

    my %DATA_STRUCT;

    while(my($type, $indv_hash) = each %{$self->{bit}}){
	foreach my $indv (keys %{$self->{bit}{$type}}){
	    my $bits = scalar @{$self->{bit}{$type}{$indv}};
	    my $vector = Bit::Vector->new($bits);
	    $vector->from_Bin(join "", @{$self->{bit}{$type}{$indv}});
	    $DATA_STRUCT{$type}{$indv} = $vector;
	    delete $self->{bit}{$type}{$indv};
	}
    }
    return \%DATA_STRUCT;
}
#-----------------------------------------------------------------------------   
sub LINKER{
    my ($self, $indvs, $args) = @_;
    my $count = 0;
    my %markers;
    my @pdl; 
    my %BFH;
    
  SEQ: foreach my $seq (@{$args}){
      
      my $tabix = Tabix->new(-data => $self->{'file'});
      
      my @iters;
      
      push @iters, $tabix->query($seq, 0, 100000);
      push @iters, $tabix->query($seq, 2400000, 2600000);

      my $start_end = 0;  

    START_END: foreach my $it (@iters){ 
	my $n_indv = scalar @{$indvs};
	my $n_alleles = 2 * $n_indv;
	my @pos;
      LINE: while(my $l = $tabix->read($it)){
	  
	  $self->{line}{raw} = $l;
	  $self->_Parse_Line();
	  
	  next LINE unless $self->{line}{refined}{type} eq 'SNV';
	  
	  my $ref = $self->{line}{refined}{ref}[0];
	  $self->_Group($indvs);
	  my @counts  = _Count_Genotypes($self->{line}{group_A});
	  my @counts2 = _Count_Genotypes($self->{line}{group_B});
	  
	  my @alleles = sort {$a cmp $b} keys %{$counts[1]};
	  
	  next LINE if scalar (grep {!/\^|$ref/} @alleles) != 1;
	  next LINE if $counts[0]->{genotype_counts}{nocall}  > 0;
	  next LINE if $counts2[0]->{genotype_counts}{nocall}  > 0;
	  next LINE if ($counts[1]->{$alleles[0]} / $n_alleles) > 0.95 || ($counts[1]->{$alleles[0]} / $n_alleles) < 0.05;
	  $markers{$count} = $self->{line}{refined}{start};
	  
	  my @vals;
	  my @non_ref = grep {!/$ref|\^/} @alleles;

	  my $hets = "$alleles[0]:$alleles[1]";
	  my $homr = "$ref:$ref";
	  my $homn = "$non_ref[0]:$non_ref[0]";

	GENOTYES: foreach my $key (@{$indvs}){
	    push @vals, 1 if $self->{line}{genotypes}{$key}{genotype} eq $homr;
	    push @vals, 2 if $self->{line}{genotypes}{$key}{genotype} eq $hets;
	    push @vals, 3 if $self->{line}{genotypes}{$key}{genotype} eq $homn;
	}

	  next LINE if ! defined $vals[0];
	  push @pos,  \@vals;
      }
	$start_end++;	
	$BFH{$seq}{$start_end} = pdl @pos;
    }    
  }

    my @keys = (1,2);

    my %no_repeat;
    
  FIRST_SEQ: foreach my $seq1 (sort {$a cmp $b} keys %BFH){
      F_E1: foreach my $n1 (@keys){
	  my $pdla = $BFH{$seq1}{$n1};
	SECOND_SEQ: foreach my $seq2 (sort {$a cmp $b} keys %BFH){
	    next SECOND_SEQ if $seq1 eq $seq2;
	    F_E2: foreach my $n2 (@keys){
		
		my $pdlb = $BFH{$seq2}{$n2};
		my ($s1, $offA) = split /_/, $seq1;
		my ($s2, $offB) = split /_/, $seq2;

		next SECOND_SEQ if $no_repeat{$seq1}{$seq2}{$n1}{$n2} == 1;
		
		$no_repeat{$seq1}{$seq2}{$n1}{$n2} = 1;
		$no_repeat{$seq2}{$seq1}{$n1}{$n2} = 1;
		
		my $n_rowa =  $pdla->slice('0,:')->nelem;
		my $n_rowb =  $pdlb->slice('0,:')->nelem;
		
		my $total = 0;
		my $sum = 0;

		my %r_data;
				
	      OUTER: for (my $i = 0; $i < $n_rowa; $i++){
		  my $loc_a = $pdla->PDL::slice(":,$i");
		  next OUTER if $n_rowa == 1;
		  next OUTER if $n_rowb == 1;
		INNER: for (my $j = 0; $j < $n_rowb; $j++){
		    next INNER if defined $r_data{$i}{$j};
		    my $loc_b = $pdlb->PDL::slice(":,$j");
		    my $cor =  $loc_a->corr($loc_b)->at(0);
		    $r_data{$i}{$j} = 1;
		    $sum++ if $cor**2 > 0.20; 
		    $total++;
		}
	      }
		print "$offA\t$offB\t$sum\t$total\t$n1\t$n2\n";
	    }
	}
      }
  }   
}

#-----------------------------------------------------------------------------   
# To generate basic statistics about files.
# Primarly counting
sub COUNT{
    use Data::Tabular::Dumper;
    
    my ($self, $args) = @_; 
    my %COUNT_DAT;
    my $tabix = Tabix->new(-data => $self->{'file'});
    
  SEQ: foreach my $seq (@{$args}){
      my $it = $tabix->query($seq);
    LINE: while(my $l = $tabix->read($it)){
	$self->{line}{raw} = $l;
	$self->_Parse_Line();
	my $ref = @{$self->{line}{refined}{ref}}[0];
	my $type = $self->{line}{refined}{type};
      INDV: foreach my $k (keys %{$self->{line}{genotypes}}){
	  my @gen = split /:/, $self->{line}{genotypes}{$k}{genotype};
	  if($gen[0] ne $gen[1]){
	      $COUNT_DAT{$k}{het_site}{$type}++;
	      $COUNT_DAT{total}{het_site}{$type}++;
	  }
	  foreach my $g (@gen){
	      if ($g eq '^'){
		  $COUNT_DAT{$k}{nocall_site}{$type}++;
		  $COUNT_DAT{total}{nocall_site}{$type}++;
		  next INDV;
	      }
	      if($g ne $ref){
		  $COUNT_DAT{$k}{non_ref_site}{$type}++;
		  $COUNT_DAT{total}{non_ref_site}{$type}++;
		  next INDV;
	      }
	  }
      }
    }
  }

    while(my($k1, $sh1) = each %COUNT_DAT){
	while(my($k2, $sh2) = each %{$COUNT_DAT{$k1}}){
	    foreach my $k3 (keys %{$COUNT_DAT{$k1}{$k2}}){
		my $val = $COUNT_DAT{$k1}{$k2}{$k3};
		print "$k1\t$k2\t$k3\t$val\n";
	    }
	}
    }

}    


#-----------------------------------------------------------------------------   
# per loci tajima's d.  The supportive subroutines are labled with a _TD
sub Tajimas_D{

    my ($self, $groups, $scaff) = @_;

    my $t = Tabix->new(-data => $self->{'file'});
    my $it = $t->query($scaff);

    next SCAFF if ! defined $it;

    my @DAT;

  LINE: while(my $l = $t->read($it)){
      $self->{line}{raw} = $l;
      $self->_Parse_Line();
      $self->_Group($groups);
      my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
      my @gens = grep {!/\^/} @alleles;
      next LINE if (scalar grep {!/\^/} @alleles) != 2;
  }
}
#-----------------------------------------------------------------------------   
sub A_1_2_TD{
    my $n_hap = shift;

    my $a1 = 0;
    my $a2 = 0;
    for(my $i = 1; $i <= ($n_hap -1); $i++){
	$a1 += (1 / $i);
	$a2 += (1 / $i**2);
    }
    return ($a1, $a2);
}
#-----------------------------------------------------------------------------   
sub B_1_2_TD{
    my $n_hap = shift;

    my $b1 = ($n_hap + 1) / (3*($n_hap - 1));
    my $b2 = (2*($n_hap**2 + $n_hap + 3))/(9*$n_hap*($n_hap -1));
    return $b1, $b2;
}
#-----------------------------------------------------------------------------   
sub C_1_2_TD{
    my ($a1, $a2, $b1, $b2, $n_hap) = @_;
   
    my $c1 = $b1 - (1 / $a1);
    my $c2 = $b2 - (($n_hap + 2)/($a1*$n_hap)) + ($a2/$a1**2);
    return $c1, $c2;
}
#-----------------------------------------------------------------------------   
sub E_1_2_TD{
    my ($a1, $a2, $c1, $c2) = @_;
    
    my $e1 = $c1/$a1;
    my $e2 = $c2/($a1**2+$a2);
    return $e1, $e2;
}
#-----------------------------------------------------------------------------   

sub COAL_T{

my ($self, $groups, $scaff) = @_;
my $t = Tabix->new(-data => $self->{'file'});
my $it = $t->query($scaff);

next SCAFF if ! defined $it;

my @DAT;

 LINE: while(my $l = $t->read($it)){
     $self->{line}{raw} = $l;
     $self->_Parse_Line();
     $self->_Group($groups);
     my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
     my @gens = grep {!/\^/} @alleles;
     next LINE if (scalar grep {!/\^/} @alleles) != 2;
     my @tot          = _Count_Genotypes($self->{'line'}{'genotypes'});
     next LINE if $tot[0]->{allele_counts}{nocall} > $tot[0]->{allele_counts}{called};
     my @g_a          = _Count_Genotypes($self->{'line'}{'group_A'});

 }
}

#-----------------------------------------------------------------------------   
# Taken from Hao in VAAST

sub get_tag_from_string {
    my $item = shift;
    
    my @tags;
    return \@tags unless defined $item;
    
    my @contiguous = split /,/, $item;
    foreach my $one (@contiguous) {
	my ($a, $b) = split /-/, $one;
	if (defined $b) {
	    for (my $i = $a; $i<=$b; $i++) {
		push @tags, $i;
	    }
	}
	else {  push @tags, $a; }
    }

    return \@tags;
}

#----------------------------------------------------------------------------- 
sub crap2{

    my ($self, $groups, $scaff, $window_len, $step) = @_;
    my $t = Tabix->new(-data => $self->{'file'});
    my $it = $t->query($scaff);
    return if ! defined $it;
    
    my %NON_REF_ALLELES;
    my $non_ref = \%NON_REF_ALLELES;
    my $current_window = $window_len;
    
  LINE: while(my $l = $t->read($it)){
      $self->{line}{raw} = $l;
      $self->_Parse_Line();
      $self->_Group($groups);
      next LINE unless $self->_type eq 'SNV';
      my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
      my @gens = grep {!/\^/} @alleles;
      next LINE if (scalar grep {!/\^/} @alleles) != 2;
      my $ref = @{$self->{line}{refined}{ref}}[0];
      my ($alt) = grep {!/\^|$ref/} @alleles;
      
      # genotype counts 

#0  HASH(0x5f7fbf0)   'allele_counts' => HASH(0x6358000)
#      'called' => 16
#      'nocall' => 66
#   'genotype_counts' => HASH(0x6357fe0)
#      'called' => 8
#      'nocall' => 33
#1  HASH(0x5f7f8f0)
#   'A' => 5
#   'G' => 11
#   '^' => 66
#2  HASH(0x5f7f920)
#   'A:A' => 2
#   'A:G' => 1
#   'G:G' => 5
#   '^:^' => 33
      
      my @all        = _Count_Genotypes($self->{'line'}{'genotypes'});
      my @target     = _Count_Genotypes($self->{'line'}{'group_A'});
      my @background = _Count_Genotypes($self->{'line'}{'group_B'});
      
      # position

      my $loci = $self->start;
      
      #number of non reference alleles at that given position

      $non_ref->{$loci}{alt}{all}        = $all[1]->{$alt};
      $non_ref->{$loci}{alt}{target}     = $target[1]->{$alt};
      $non_ref->{$loci}{alt}{background} = $background[1]->{$alt};
      
      # number of called chromosomes

      $non_ref->{$loci}{called}{all}        = $all[0]->{allele_counts}{called};
      $non_ref->{$loci}{called}{target}     = $target[0]->{allele_counts}{called};
      $non_ref->{$loci}{called}{background} = $background[0]->{allele_counts}{called};

     ($current_window, $non_ref) = _process_window($non_ref, $current_window, $window_len, $loci, $step) if $loci >= $current_window;
  }
}
#----------------------------------------------------------------------------- 
sub _process_window{
    
    my ($allelic_count, $current_window, $window_len, $current_pos, $offset, $chr_targ, $chr_back) = @_;
    my $step = $current_window + $offset;

    my %spectra;

    while(my ($loci_pos, $type_hash) = each %{$allelic_count}){
	while(my ($set, $count) = each %{$type_hash->{alt}}){
	    $count = 0 if !defined $count;
	    $spectra{$set}{$count}++;
	    delete $allelic_count->{$loci_pos}{$set} if $loci_pos < $step;
	}	
    }

    my %total_seg;
  
    while(my ($loci_pos, $type_hash) = each %{$allelic_count}){
	$total_seg{target}{loci}++;
	$total_seg{background}{loci}++;
	$total_seg{all}{loci}++;
	while(my ($set, $count) = each %{$type_hash->{called}}){   
	    if (defined $allelic_count->{$loci_pos}{called}{$count}){
		$total_seg{$set}{called_alleles} += $allelic_count->{$loci_pos}{called}{$count};
		$total_seg{$set}{max} = $allelic_count->{$loci_pos}{called}{$count} if $allelic_count->{$loci_pos}{called}{$count} > $total_seg{max};
		delete $allelic_count->{$loci_pos}{$set} if $loci_pos < $step;
	    }
	}
    }
    
    my $T = _CDF($spectra{target}, 100);
    my $B = _CDF($spectra{background}, 100);
    my $A = _CDF($spectra{all}, 100);
    
    Theta($spectra{target},$total_seg{target}) if $total_seg{target}{loci} > 10;

    my $D = _KS_TEST($T, $B, $current_window);
    
    $current_window += $window_len; 
    return($current_window, $allelic_count);
}
#----------------------------------------------------------------------------- 
sub _KS_TEST{

    my ($target_cdf, $background_cdf, $current_window) = @_;
    my $max_diff = 0;
    foreach my $non_ref (keys %{$target_cdf}){
	my $diff = abs($target_cdf->{$non_ref} - $background_cdf->{$non_ref});
	$max_diff = $diff if $diff > $max_diff;
    }
    print "$max_diff\t$current_window\n";
    return($max_diff);
}
#----------------------------------------------------------------------------- 
sub _CDF{
    
    my ($spectra, $max) = @_;
    my $sum = 0;
  NON_REF: for(my $i = 0; $i <= $max; $i++){
      next NON_REF if ! exists $spectra->{$i};
      $sum += $spectra->{$i};	
  }
    my %cdf;
    my $last = 0;
  NON_REF: for(my $i = 0; $i <= $max; $i++){
      $cdf{$i} = $last;
      next NON_REF if! defined $spectra->{$i};
      my $current = ($spectra->{$i} + $last) / $sum;
      $cdf{$i} += $current;
      $last = $current;
  }
    
    return(\%cdf);
}
#----------------------------------------------------------------------------- 
sub Theta {
    my ($data, $total_seg) = @_;
    my $a_n = 0;

    my $avg_chrs = $total_seg->{called_alleles} / $total_seg->{loci};

    for (my $i = 1; $i < 2 * $avg_chrs; $i++){
        $a_n += 1 / $i;
    }
 
    my $num          = 0;
    my $fu_sum_theta = 0;
    my $H_sum_theta  = 0;
    my $PI_sum_theta = 0;

  NON_REF: foreach my $i (keys %{$data}){
      my $i_i     = $i;
      my $Xi      = $i_i / $total_seg->{loci};
      my $Theta_i = $i * $Xi;

      $H_sum_theta  += $Theta_i * $i;
      $PI_sum_theta += $Theta_i * ($total_seg->{max} - $i);
      $num          += $Theta_i * (1 / $i) if $i > 0;
      $fu_sum_theta += $Theta_i;


  }  
      my $two_over_chr = 2/($avg_chrs * ($avg_chrs  - 1));    
      my $Pi_theta        = $PI_sum_theta * $two_over_chr;
      my $Fu_theta        = (1 / $avg_chrs) * $fu_sum_theta;
      my $H__theta        = $H_sum_theta * $two_over_chr;
      my $Watterson_theta = $num / $a_n;
   
    return [$Watterson_theta, $Fu_theta, $Pi_theta, $H__theta];

}
#----------------------------------------------------------------------------- 

sub MAF_WINDOWS{
    my ($self, $scaff_info) = @_;
    my @scaffold_info = @{$scaff_info};
    my %MAF_COUNTS;
    my $n_windows;
    my @RESULTS;

  CHRS: while (@scaffold_info){
      
      my $seqid = shift @scaffold_info;
      my $total_len = shift @scaffold_info;      
      my $t = Tabix->new(-data => $self->{'file'});

      my $last_pos = 0;

      print STDERR "working on feature: $seqid\n";
      
    WINDOW: for(my $i = 10_000; $i < $total_len; $i += 10_000){

	my $n_seg    = 0;
	
	my $it = $t->query($seqid, $last_pos, $i);
	return if ! defined $it;
	
	my $window_maf_count = 0;
	
      LINE: while(my $l = $t->read($it)){
	  $self->{line}{raw} = $l;
	  $self->_Parse_Line();
#	  $self->_Group($groups);
	  next LINE unless $self->_type eq 'SNV';
	  my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
	  my @gens = grep {!/\^/} @alleles;
	  next LINE if (scalar grep {!/\^/} @alleles) != 2;
	  my $ref        = @{$self->{line}{refined}{ref}}[0];
	  my ($alt)      = grep {!/\^|$ref/} @alleles;
	  my @all        = _Count_Genotypes($self->{'line'}{'genotypes'});
	  $n_seg++;
	  my $maf = $all[1]{$alt} / $all[0]{'allele_counts'}{'called'};
	  $window_maf_count++ if $maf > .5;
      }	
	$last_pos = $i;    
	$MAF_COUNTS{$window_maf_count}++;
	$n_windows++;	    
	push @RESULTS, $window_maf_count;	
    }
  }
    while(my ($count, $number) = each %MAF_COUNTS){
	print "COUNT_RESULTS:$count\t$number\n";
    }
    foreach my $results (@RESULTS){
	print "ALL:$results\n";
    }
}

1;
#----------------------------------------------------------------------------- 
sub Print_genotype{
    my ($self, $groups, $scaff, $start, $end) = @_;
    my $t = Tabix->new(-data => $self->{'file'});
    my $it = $t->query($scaff, $start, $end);
    return if ! defined $it;
     
  LINE: while(my $l = $t->read($it)){
      $self->{line}{raw} = $l;
      $self->_Parse_Line();
      $self->_Group($groups);
      next LINE unless $self->_type eq 'SNV';
      my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
      my @gens = grep {!/\^/} @alleles;
    next LINE if (scalar grep {!/\^/} @alleles) != 2;
      my $ref = @{$self->{line}{refined}{ref}}[0];
      my ($alt) = grep {!/\^|$ref/} @alleles;
      foreach my $tindv (sort {$a <=> $b} keys %{$self->{'line'}{'group_A'}}){
	  print "$tindv\t$self->{'line'}{'group_A'}{$tindv}{genotype}\n";
      }
      foreach my $bindv (sort {$a <=> $b} keys %{$self->{'line'}{'group_B'}}){
          print"$bindv\t$self->{'line'}{'group_B'}{$bindv}{genotype}\n";
      }
  }
}
