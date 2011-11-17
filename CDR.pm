package CDR;
use Tabix;
use strict;
use warnings;
use Set::IntSpan::Fast;


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
    my $seqids = $t->getnames(); 
    $self->{'seqids'} = $seqids;
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
    while(my $l = $t->read($i)){
	$self->{'line'}{'raw'} = $l;
	$self->_Parse_Line;
	$self->Genotype_Counting;
	$self->Expected_Hetrozygosity;
	#$self->HWE_Departure;
    }
}

#-----------------------------------------------------------------------------   
sub _Parse_Line{

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

    $self->_Parse_Genotype;
    $self->_Load_Reference_Genotypes;
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
	    die "overlapping genotypes : $self->{'line'}{'raw'}" if exists $genotype{$i}{'genotype'};
	    $genotype{$i}{'genotype'} = join ":", sort {$a cmp $b} split /:/, $g[1];
	    $genotype{$i}{'effect'} = $g[2];
	}
    }
    $self->{'line'}{'genotypes'} = \%genotype;
}

#-----------------------------------------------------------------------------   
sub _Load_Reference_Genotypes{
   
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

sub HapMap_FST{

    my $self = shift;
    my @alleles = sort {$a cmp $b} keys %{$self->{'allele_counts'}};
    return if scalar grep {!/^/} @alleles >  2;
    return if scalar grep {!/^/} @alleles == 1;
    my $total_alleles = $self->{'allele_counts'}{$alleles[0]} + $self->{'allele_counts'}{$alleles[1]};
    
}

#-----------------------------------------------------------------------------   

sub Genotype_Counting{

    my $self = shift;
    my %allele_counts;
    my %genotype_counts;


    while(my($indv, $info_hash) = each %{$self->{'line'}{'genotypes'}}){
	while( my($info, $value) = each %{$info_hash}){
	    if ($info =~ /genotype/){
		my @gen = split /:/, $value;
		map {$allele_counts{$_}++} @gen;
		$genotype_counts{$value}++;
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
    
    $self->{'genotype_counts'} = \%genotype_counts;
    $self->{'allele_counts'}   = \%allele_counts;
    $self->{'total_counts'}    = \%total;
    my @alleles = sort {$a cmp $b} keys %{$self->{'allele_counts'}};
    $self->{'alleles'} = \@alleles;
}
#-----------------------------------------------------------------------------   
sub Expected_Hetrozygosity{
    my $self = shift;
    
    my @alleles = @{$self->{'alleles'}};

    return if !defined $self->{'total_counts'}{'allele_counts'}{'called'};
    return if scalar (grep {!/\^/} @alleles) > 2;
    if(scalar (grep {!/\^/} @alleles) < 2){
	
	my $expected_hets = 0;
	my $observed_hets = 0;

	if(defined $self->{'genotype_counts'}{"$alleles[0]:$alleles[1]"}){
	    $observed_hets = $self->{'genotype_counts'}{"$alleles[0]:$alleles[1]"};
	}
	print "$expected_hets\t$observed_hets\n";
	$self->{'expected_hets'}  = $expected_hets;
        $self->{'observed_hets'}  = $observed_hets;
    }
    else{
	my $total_alleles = $self->{'total_counts'}{'allele_counts'}{'called'};
	my $p = $self->{'allele_counts'}{$alleles[0]} / $total_alleles;
	my $q = $self->{'allele_counts'}{$alleles[1]} / $total_alleles;
	my $expected_hets = 2 * $p * $q * $total_alleles / 2;
	$expected_hets = int $expected_hets;
	my $observed_hets = 0;
	if (defined $self->{'genotype_counts'}{"$alleles[0]:$alleles[1]"}){
	    $observed_hets = $self->{'genotype_counts'}{"$alleles[0]:$alleles[1]"};
	}
	$self->{'expected_hets'}  = $expected_hets;
	$self->{'observed_hets'}  = $observed_hets;
	print "$expected_hets\t$observed_hets\n";
    }
}
#-----------------------------------------------------------------------------   

sub HWE_Departure{

    my $self = shift;
    my @alleles = @{$self->{'alleles'}};

    return if !defined $self->{'total_counts'}{'allele_counts'}{'called'};
    return if scalar (grep {!/\^/} @alleles) > 2;

    my $total_alleles = $self->{'total_counts'}{'allele_counts'}{'called'};
    my $total_indvs   = $self->{'total_counts'}{'genotype_counts'}{'called'};
 
    my $p = 0;
    my $q = 0;

    $p = $self->{allele_counts}{$alleles[0]} / $total_alleles if defined $self->{allele_counts}{$alleles[0]};
    $q = $self->{allele_counts}{$alleles[1]} / $total_alleles if defined $self->{allele_counts}{$alleles[1]};

    my $exp_pp = $p**2 * $total_indvs ;
    my $exp_qq = $q**2 * $total_indvs ;
    my $exp_pq = 2 * $q * $p;  
    
    my  $obs_pp = $self->{'genotype_counts'}{"$alleles[0]:$alleles[0]"} ;
    my  $obs_qq = $self->{'genotype_counts'}{"$alleles[1]:$alleles[1]"} ; 
    my  $obs_pq = $self->{'genotype_counts'}{"$alleles[0]:$alleles[1]"} ;

    $obs_pp = 0 if ! defined $obs_pp;
    $obs_qq = 0 if ! defined $obs_qq;
    $obs_pq = 0 if ! defined $obs_pq;

    my $chi = ($obs_pp - $exp_pp)**2 + ($obs_pq - $exp_pq)**2  + ($obs_qq - $exp_qq)**2;
    if ($chi > 3.84){
	print "1\n";
    }
    else{
	print "0\n";
    }
} 


1;

#-----------------------------------------------------------------------------   

sub Rid_No_Call_Only{

    my $self = shift;

}

#-----------------------------------------------------------------------------   

sub Group{





}
