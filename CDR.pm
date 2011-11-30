package CDR;
use Tabix;
use strict;
use warnings;
use Math::CDF;
use Set::IntSpan::Fast;
use List::MoreUtils qw(uniq);

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
    system ("bgzip $self->{'file'}");
    system ("tabix -f -s 1 -b 2 -e 3 $self->{'file'}\.gz");
    $self->{'file'} = "$self->{'file'}\.gz";
    $self->_Load_Seqids();
    return $self;
}

#-----------------------------------------------------------------------------

sub _initialize_args {
    my $self = shift;

    my $args = $self->prepare_args(@_);
    # Set valid class attributes here
    my %valid_attributes = map {$_ => 1} qw(file fh sparse);
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
	my $file = $self->{'file'};
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

sub Query_Range {
    my ($self, $args) = @_;
    my $t = Tabix->new(-data => $self->{'file'});
    my $i = $t->query($args);
    LINE: while(my $l = $t->read($i)){
	next LINE if ! defined $l;
	$self->{'line'}{'raw'} = $l;
	$self->_Parse_Line;
	next LINE if !defined $self->{'line'}{'genotypes'};
	my $flag = $self->_Scrub_No_Call;
	$flag = 0 if $self->{'line'}{'refined'}{'type'} ne 'SNV';
	next LINE if $flag == 0;
	my @group = (40, 0, 3, 31, 35, 19, 24);
	$self->FST(\@group);
	#$self->HWE_Departure;	
    }
}


#-----------------------------------------------------------------------------   

sub _Parse_Line{
    my $self = shift;
    $self->_Parse_Raw_Line;
    $self->_Parse_Genotype;
    $self->_Parse_Reference_Genotypes;
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
	my $set = Set::IntSpan::Fast->new();
        $set->add_from_string($g[0]);
	foreach my $i ($set->as_array){
	    print STDERR "overlapping genotypes : $self->{'line'}{'raw'}\n" if exists $genotype{$i}{'genotype'};
	    return if exists $genotype{$i}{'genotype'};
	    $genotype{$i}{'genotype'} = join ":", sort {$a cmp $b} split /:/, $g[1];
	    $genotype{$i}{'effect'} = $g[2];
	}
    }
    $self->{'line'}{'genotypes'} = \%genotype;
}

#-----------------------------------------------------------------------------   
sub _Parse_Reference_Genotypes{
   
    my $self  = shift;
    while(my($key, $file) = each %{$self->{'file_keys'}}){
	if(! defined $self->{'line'}{'genotypes'}{$key}){
	    $self->{'line'}{'genotypes'}{$key}{'genotype'}  = "$self->{'line'}{'refined'}{'ref'}[0]:$self->{'line'}{'refined'}{'ref'}[0]";
	    if (defined $self->{'line'}{'refined'}{'ref'}[1]){
		$self->{'line'}{'genotypes'}{$key}{'effect'}  = "$self->{'line'}{'refined'}{'ref'}[1]:$self->{'line'}{'refined'}{'ref'}[1]";	
	    }
	    else{
		$self->{'line'}{'genotypes'}{$key}{'effect'}  = undef;
	    }
	}
    }
}
    
#-----------------------------------------------------------------------------   

sub _Parse_Alleles{

    my $genotypes = shift;
    my @genotypes;


    while(my ($indv, $g) = each %{$genotypes}){
	while(my ($info, $d) = each %{$g}){
	    push @genotypes, split /:/, $d if $info eq 'genotype';
	}
    }

    my @uniq = sort {$a cmp $b} uniq @genotypes;
    return \@uniq;
}

#-----------------------------------------------------------------------------   

sub _Count_Genotypes{

    my $info_hash = shift;

    my %allele_counts;
    my %genotype_counts;
 
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

    my $self = shift;

    my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
    my @counts  = _Count_Genotypes($self->{'line'}{'genotypes'});
    

    return if scalar (grep {!/\^/} @alleles) > 2;
    return if scalar (grep {!/\^/} @alleles) == 1;

  #  return if ! defined $counts[2]->{"$alleles[0]:$alleles[1]"};     
  #  return if ! defined $counts[2]->{"$alleles[0]:$alleles[0]"};
  #  return if ! defined $counts[2]->{"$alleles[1]:$alleles[1]"};
  #  
  #  my $total_alleles   = $counts[0]->{'allele_counts'}{'called'};
  #  my $total_genotypes = $counts[0]->{'genotype_counts'}{'called'};
  #  my $p               = 0;
  #  my $q               = 0;
  #  
  #  $p = $counts[1]->{$alleles[0]} / $total_alleles if defined $counts[1]->{$alleles[0]};
  #  $q = $counts[1]->{$alleles[1]} / $total_alleles if defined $counts[1]->{$alleles[1]};
  #  
  #  my $exp_pp = int $p**2  * $total_genotypes;
  #  my $exp_qq = int $q**2  * $total_genotypes;
  #  my $exp_pq = int 2 * $p * $q * $total_genotypes; 

    my $AA = 0;
    my $AB = 0; 
    my $BB = 0;


    $AB = $counts[2]->{"$alleles[0]:$alleles[1]"} if defined $counts[2]->{"$alleles[0]:$alleles[1]"};
    $AA = $counts[2]->{"$alleles[0]:$alleles[0]"} if defined $counts[2]->{"$alleles[0]:$alleles[0]"};
    $BB = $counts[2]->{"$alleles[1]:$alleles[1]"} if defined $counts[2]->{"$alleles[1]:$alleles[1]"};

     

    print "$self->{line}{refined}{start}\t$AA\t$AB\t$BB\n";


#    my $p_pp =  _Chi_Lookup($exp_pp, $AA);
#    my $p_pq =  _Chi_Lookup($exp_pq, $AB);
#    my $p_qq =  _Chi_Lookup($exp_qq, $BB);
#
#    my $chi = $p_pp + $p_pq + $p_qq;
#    my $p_value   = Math::CDF::pchisq($chi, 1);
#    print "$p_value\n";
} 

#-----------------------------------------------------------------------------   

sub _Chi_Lookup{

    my ($exp, $obs) = @_;
    my $chi;
   
    if ($exp != 0 ){
	$chi = ($exp - $obs)**2 / $exp;
    }
    else{
	$chi = 0;
    }
    return $chi;
}
    
#-----------------------------------------------------------------------------   

sub FST{

    my ($self, $groups) = @_;
   

    my @alleles = @{_Parse_Alleles($self->{'line'}{'genotypes'})};
    return if scalar (grep {!/\^/} @alleles) == 1;
    return if scalar (grep {!/\^/} @alleles) > 2;

    $self->_Group($groups);
    my ($a_t_counts, $a_a_counts)      = _Count_Genotypes($self->{'line'}{'group_A'});
    my ($b_t_counts, $b_a_counts)      = _Count_Genotypes($self->{'line'}{'group_B'});
    my ($total_counts, $total_alleles) = _Count_Genotypes($self->{'line'}{'genotypes'});

    my ($A,$B) = keys %{$total_alleles};
    my $a_c    = $total_alleles->{$A};
    my $b_c    = $total_alleles->{$B};
    
    my $minor_allele = $a_c < $b_c ? $A : $B;
    my $n            = $total_counts->{'allele_counts'}{'called'};
    my $x            = $total_alleles->{$minor_allele} / $n;
    my $n_j          = $total_counts->{'allele_counts'}{'called'};
    $n_j             = $total_counts->{'allele_counts'}{'called'} + $total_counts->{'allele_counts'}{'nocall'} 
    if defined $total_counts->{'allele_counts'}{'nocall'};
    
    my  $x_ja = defined $a_a_counts->{$A} || defined $a_a_counts->{$B} ? 
        $a_a_counts->{$minor_allele} / ($a_a_counts->{$B} + $a_a_counts->{$A}) : 0;
    my  $x_jb = defined $b_a_counts->{$A} || defined $b_a_counts->{$B} ? 
	$b_a_counts->{$minor_allele} / ($b_a_counts->{$A} + $b_a_counts->{$B}) : 0; 
    q
    my $n_ja = $a_t_counts->{'allele_counts'}{'called'}; 
    my $n_jb = $b_t_counts->{'allele_counts'}{'called'};

    my $dem = 2 * $n / ($n - 1) * $x * (1 - $x);
    
    my $main_num_a = 2 * $n_ja/($n_ja -1) * $x_ja * (1 - $x_ja);
    my $main_num_b = 2 * $n_jb/($n_jb -1) * $x_jb * (1 - $x_jb);
    
    my $n_ja_2 = $n_ja * ($n_ja -1 )/2;
    my $n_jb_2 = $n_jb * ($n_jb -1 )/2;
    my $n_jt_2 = $n_ja_2 + $n_jb_2;

    my $a_final = $n_ja_2 * $main_num_a;
    my $b_final = $n_jb_2 * $main_num_b;

    my $FST =  1 - ((($a_final + $b_final) / $n_jt_2 ) / $dem);
    print "$self->{'line'}{'refined'}{'seqid'}\t$self->{'line'}{'refined'}{'start'}\t$FST\n";
}

#-----------------------------------------------------------------------------   

sub _Group{

    my %groupA;
    my %groupB;

    my %groupA_DAT;
    my %groupB_DAT;

    my ($self, $groups) = @_;
    
    foreach my $i (@{$groups}){
	$groupA{$i} = 1;
    }
    
    while(my($indv, $indv_hash) = each %{$self->{'line'}{'genotypes'}}){
	if (exists $groupA{$indv}){
	    $groupA_DAT{$indv} = $indv_hash;
	}
	else{
	    $groupB_DAT{$indv} = $indv_hash;
	} 
    }
    $self->{'line'}{'group_A'} = \%groupA_DAT;
    $self->{'line'}{'group_B'} = \%groupB_DAT; 
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



1;


