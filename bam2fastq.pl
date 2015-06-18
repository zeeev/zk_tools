#!/usr/bin/perl

use forks;
use forks::shared;

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil :sys_wait_h);
use Compress::Zlib;
use File::Copy;
use File::Spec;
use File::Which;
use File::Basename;
use File::Temp qw(tempfile);
use IPC::Shareable qw(:lock);
use Storable qw(store retrieve);

our $LFS;
BEGIN {
    binmode(STDIN);
    binmode(STDOUT);
    select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
    select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

    $SIG{INT} = sub {exit(130)};

    mkdir("$ENV{HOME}/.Inline") if(! -d "$ENV{HOME}/.Inline");

    $LFS = File::Which::which('lfs');
}

#use Inline (Config => DIRECTORY => File::Spec->tmpdir());
use Inline qw(C);

my ($exe) = $0 =~ /([^\/]+)$/;
my $usage = "
Usage:

     $exe <bam_file>

     Converts a BAM file to FASTQ format. Files sorted on query name convert
     faster but are not required. Paired end reads with missing pairs will be
     ignored but produce a warning. You can optionally save these missing pairs
     using the -fq3 option. Single end reads mixed in with paired end reads
     produces a fatal error unless you provide a location to store single end
     reads using the -fq3 option. Setting -fq3 to /dev/null will supress any
     warnings and errors related to both single end reads and missing pairs.

Options:
     id        <STRING>   Auto-generate file names using this ID. Will be
                          generated from readgroup info if -id is supplied with
                          an empty string or -id is a path to a directory (this
                          option assumes paired end reads)

     fq         <PATH>    Output FASTQ file (STDOUT otherwise)

     fq2        <PATH>    Optionally split second pair into separate output file

     fq3        <PATH>    Optional file for unmatched paired and single end reads

     gzip|z               Compress fastq output files (gzip readable)

     validate|v           Validates the ISIZE and CRC32 values in bam headers.
                          This is like doing an an MD5 checksum for each file.
                          Integrity of compressed output will also be checked.

     q64                  Fix Phred+64 BAM qaulity scores (should be Phred+33)

     solexa               Fix Solexa odds based qaulity scores (should be Phred+33)

     restore    <TAG>     Restore original quality string using values from the
                          given BAM tag

     cpus|c     <INT>     CPUs to use for conversion

     stripes    <INT>     Lustre stipes to use. Only valid for Lustre file
                          systems. Defaults to 4.

     IO         <DIR>     Hack to reduce IO on CHPC by writing intermediate files to
                          this location (preferablely an in-memory location)

     range|b  <INT,INT>   Only convert part of file (start, end). Zero based
                          (Not yet implemented)

     help|?               Prints this usage statement

";

my $cpus = 1;
my $stripes = 4;
my $range;
my $outfile1; #paired end
my $outfile2; #2nd pair
my $outfile3; #unmatched pairs
my $gzip;
my $q64;
my $solexa;
my $restore;
my $validate;
my $io_dir;
my $id;

#get options from command line
my @argv = @ARGV; #backup
GetOptions("cpus|c=i" => \$cpus,
	   "fq=s" => \$outfile1,
	   "fq2=s" => \$outfile2,
	   "fq3=s" => \$outfile3,
	   "gzip|z" => \$gzip,
	   "q64" => \$q64,
	   "solexa" => \$solexa,
	   "restore=s" => \$restore,
	   "validate|v" => \$validate,
	   "IO=s" => \$io_dir,
	   "id:s" => \$id,
	   "stripes=i" => \$stripes,
           "range|b=s" => \$range,
           "help|?" => sub{print $usage; exit(0)});

$q64 = 2 if($solexa); #solexa goes in the $q64 value

#make sure the id option didn't eat the file name
if($id && !@ARGV && $id =~ /\.bam$/){
    push(@ARGV, $id);
    $id = '';
}

my $file = $ARGV[0];
if(!$file){
    print $usage;
    exit(0);
}
die "ERROR: The file $file does not exist\n" if(! -f $file);
die "ERROR: Cannot specify fq2 without fq\n" if($outfile2 && !$outfile1);

#generate needed file names
if(defined($id)){
    if(length($id) > 0 && ! -d $id){
	#make file names
	my $name = "$id";
	$outfile1 ||= "./$name\_1.fq";
	$outfile2 ||= "./$name\_2.fq";
	$outfile3 ||= "./$name.fq";
    }
    else{
	#get read group from header
	my $header = get_header($file, $validate);
	my @rgs = grep {/^\@RG\t/} split(/\n/, $header->{text});

	die "ERROR: No read groups specified in BAM\n" if(!@rgs);
	warn "WARNING: Multiple read groups specified in BAM. I will just use the first.\n" if(@rgs > 1);

	my %rg = map {/^([^\:]+)\:(.*)$/} grep {!/^\@RG/} split(/\t/, $rgs[0]);
	my $hex = unpack('H4', pack('S<', crc32($rgs[0]))); #helps force uniqueness

	#fill in missing values
	if(!defined($rg{PU})){
	    $rg{PU} = 'L1';
	}
	if(!defined($rg{SM})){
	    $rg{SM} = $rg{ID};
	}

	#values for file name
	my $id   = $rg{SM};
	my $lane = $rg{PU};

	#fix lane values
	if($rg{PL}){
	    $lane =~ s/^${rg{PL}}[\-_](\d+)/L$1/;
	}

	#fix for WashU data
	if(defined($rg{SM}) && defined($rg{LB}) && $rg{LB} =~ /^\"$rg{SM}\-.*lib\d+\"$/){
	    my ($washu_prefix) = $rg{SM} =~ /^(H_[A-Z]{2}\-)[^\-]+\-[^\-]+$/;
	    if($washu_prefix){
		$washu_prefix = quotemeta($washu_prefix);
		$id =~ s/^$washu_prefix//;
		$lane =~ s/^.*\.(\d+)$/L$1/;
	    }
	}

	#make file names
	my $name = "$id\_$hex\_$lane";
	$outfile1 ||= "./$name\_1.fq";
	$outfile2 ||= "./$name\_2.fq";
	$outfile3 ||= "./$name.fq";
    }

    #let user know where the output will be
    print STDERR "##NOTE: Output will be sent to $outfile1, $outfile2, and $outfile3\n";
}

#io hack
my %io_files;
if($io_dir){
    die "ERROR: The directory $io_dir does not exist\n" if(!-d $io_dir);

    if($outfile1){
	$io_files{OUTFILE1} = $outfile1;
	(undef, $outfile1) = tempfile('BAM2FASTQ_XXXXX', DIR => $io_dir, SUFFIX => '.fq');
    }
    if($outfile2){
	$io_files{OUTFILE2} = $outfile2;
	(undef, $outfile2) = tempfile('BAM2FASTQ_XXXXX', DIR => $io_dir, SUFFIX => '.fq');
    }
    if($outfile3){
	$io_files{OUTFILE3} = $outfile3;
	 (undef, $outfile3) = tempfile('BAM2FASTQ_XXXXX', DIR => $io_dir, SUFFIX => '.fq');
    }
}

#make sections of right size
my $size = (stat($file))[7];
my $count = ceil($size/33554432); #32Mb sections
my $window = ceil($size/$count);
my @sections = (1..$count);
share(@sections); #share is non destructive for forks::shared

#prepare output file(s)
my $osize = ($outfile2) ? 2*$size : 4*$size;
$osize = int($osize/3) if($gzip);
foreach my $outfile ($outfile1, $outfile2, $outfile3, @io_files{qw(OUTFILE1 OUTFILE2 OUTFILE3)}){
    next if(!$outfile || -c $outfile);

    $outfile .= '.gz' if($gzip && $outfile !~ /\.gz$/); #add extension
    my (undef, $dir) = fileparse($outfile);

    #optimize for lustre file system stripes
    if($LFS && (my $max = `$LFS osts $dir 2> /dev/null | wc -l`) > 0){
	chomp($max);
	unlink($outfile);
	my $count = ($stripes < $max) ? $stripes : $max-1; #max includes header
	system("$LFS setstripe -c $count -S 32m $outfile");
    }

    get_handle('>', $outfile); #initialize
    close_handles(); #clear handle
    truncate($outfile, $osize); #estimate that outfile needs to be 4x infile
}

#make shared values I need
our $IPCHANDLE = tie(my $ipc_share0, 'IPC::Shareable', undef, { destroy => 1 });
tie(my $ipc_share1, 'IPC::Shareable', undef, { destroy => 1 });
tie(my $ipc_share2, 'IPC::Shareable', undef, { destroy => 1 });
tie(my $ipc_share3, 'IPC::Shareable', undef, { destroy => 1 });

$ipc_share0 = 0; #output buffer offset
$ipc_share1 = 0; #output buffer offset
$ipc_share2 = 0; #output buffer offset
$ipc_share3 = 0; #output buffer offset

#associate files and shared variables
our %IPC4FILE = ( PIPE => \$ipc_share0);
$IPC4FILE{$outfile1} = \$ipc_share1 if($outfile1);
$IPC4FILE{$outfile2} = \$ipc_share2 if($outfile2);
$IPC4FILE{$outfile3} = \$ipc_share3 if($outfile3);

#run all threads on data
my @threads;
my %param = (SIZE => $size, COUNT => $count, WINDOW => $window, SECTIONS => \@sections,
	     VALIDATE => $validate, Q64 => $q64, RESTORE => $restore, INFILE => $file,
	     OUTFILE1 => $outfile1, OUTFILE2 => $outfile2, OUTFILE3 => $outfile3);
for(my $i = 1; $i < $cpus; $i++){
    my $thr = threads->new({'context' => 'scalar'}, \&thread_run, \%param);
    push(@threads, $thr);
}
my $buffer1 = thread_run(\%param); #let main process join in

#gather and try and to colapse remaining buffers
my $fq1_ref;
my $fq2_ref;
my $gz1_ref;
my $gz2_ref;
{my $fastq1 = ''; $fq1_ref = \$fastq1;} #set reference
{my $fastq2 = ''; $fq2_ref = ($outfile2) ? \$fastq2 : $fq1_ref;} #set reference
{my $gz1 = ''; $gz1_ref = \$gz1;} #set reference
{my $gz2 = ''; $gz2_ref = \$gz2;} #set reference
foreach my $thr (@threads){
    my $tfile = $thr->join; #wait on threads
    my $buffer2 = retrieve($tfile);
    for (my $i = 0; $i < @$buffer2; $i++){
	$buffer1->[$i] ||= {}; #just in case

	my $hash1 = $buffer1->[$i];
	my $hash2 = $buffer2->[$i];
	foreach my $read_name (keys %$hash2){
	    if(!exists($hash1->{$read_name})){ #copy over to local hash
		$hash1->{$read_name} = $hash2->{$read_name};
		next;
	    }

	    if($hash1->{$read_name}[1] == 1){ #first in pair
		$$fq1_ref .= delete($hash1->{$read_name})->[0];
		$$fq2_ref .= $hash2->{$read_name}->[0];
	    }
	    else{ #second in pair
		$$fq1_ref .= $hash2->{$read_name}->[0];
		$$fq2_ref .= delete($hash1->{$read_name})->[0];
	    }
	}
	undef $buffer2->[$i]; #free memory
	
	#compress
	if($gzip){
	    #max ISIZE is 65536
	    deflate(\ (substr($$fq1_ref, 0, 65536, '')), $gz1_ref) while(length($$fq1_ref) >= 65536);
	    deflate(\ (substr($$fq2_ref, 0, 65536, '')), $gz2_ref) while(length($$fq2_ref) >= 65536);
	}
	else{
	    $gz1_ref = $fq1_ref;
	    $gz2_ref = $fq2_ref;
	}
	
	print_async([$gz1_ref, $gz2_ref], [$outfile1, $outfile2], 0) if(length($$gz1_ref) >= 33554432 ||
									length($$gz2_ref) >= 33554432);
    }
    unlink($tfile);
}

#gather reads with no mates
my $fq3_ref;
my $fq4_ref;
my $gz3_ref;
my $gz4_ref;
{my $fastq3 = ''; $fq3_ref = \$fastq3;} #set reference
{my $fastq4 = ''; $fq4_ref = \$fastq4;} #set reference
{my $gz3 = ''; $gz3_ref = \$gz3;} #set reference
{my $gz4 = ''; $gz4_ref = \$gz4;} #set reference
for (my $i = 0; $i < @$buffer1; $i++){
    my $hash1 = $buffer1->[$i];
    foreach my $read_name (keys %$hash1){
	if($hash1->{$read_name}[1] == 0){ #unpaired (should never happen at this stage)
	    $$fq4_ref .= delete($hash1->{$read_name})->[0];
	}
	else{ #missing mate
	    if($outfile3){
		$$fq3_ref .= delete($hash1->{$read_name})->[0];
	    }
	    else{
		print STDERR "WARNING: $read_name did not have a mate\n";
		delete($hash1->{$read_name});
	    }
	}
    }
    
    #compress
    if($gzip){
	#max ISIZE is 65536
	deflate(\ (substr($$fq3_ref, 0, 65536, '')), $gz3_ref) while(length($$fq3_ref) >= 65536);
	deflate(\ (substr($$fq4_ref, 0, 65536, '')), $gz4_ref) while(length($$fq4_ref) >= 65536);
    }
    else{
	$gz3_ref = $fq3_ref;
	$gz4_ref = $fq4_ref;
    }

    print_async([$gz3_ref, $gz4_ref], [$outfile3, $outfile3], 0) if(length($$gz3_ref) >= 33554432 ||
								    length($$gz4_ref) >= 33554432);
}

#compress what's left
if($gzip){
    #max ISIZE is 65536
    deflate(\ (substr($$fq1_ref, 0, 65536, '')), $gz1_ref) while(length($$fq1_ref));
    deflate(\ (substr($$fq2_ref, 0, 65536, '')), $gz2_ref) while(length($$fq2_ref));
    deflate(\ (substr($$fq3_ref, 0, 65536, '')), $gz3_ref) while(length($$fq3_ref));
    deflate(\ (substr($$fq4_ref, 0, 65536, '')), $gz4_ref) while(length($$fq4_ref));
}
else{
    $gz1_ref = $fq1_ref;
    $gz2_ref = $fq2_ref;
    $gz3_ref = $fq3_ref;
    $gz4_ref = $fq4_ref;
}

#write remaining data
print_async([$gz1_ref, $gz2_ref, $gz3_ref, $gz4_ref],
	    [$outfile1, $outfile2, $outfile3, $outfile3], 0);

#add end of file block
if($gzip){
    if($ipc_share0){ #output on STDOUT needs an EOF block
	$ipc_share0 += length(eof_block());
	syswrite(STDOUT, eof_block());
    }
    if($outfile1){
	my $loc1 = $ipc_share1;
	$ipc_share1 += length(eof_block());
	my $OUT = get_handle('+<', $outfile1);
	binmode($OUT);
	sysseek($OUT, $loc1, 0);
	syswrite($OUT, eof_block());
    }
    if($outfile2){
	my $loc2 = $ipc_share2;
	$ipc_share2 += length(eof_block());
	my $OUT = get_handle('+<', $outfile2);
	binmode($OUT);
	sysseek($OUT, $loc2, 0);
	syswrite($OUT, eof_block());
    }
    if($outfile3){
	my $loc3 = $ipc_share3;
	$ipc_share3 += length(eof_block());
	my $OUT = get_handle('+<', $outfile3);
	binmode($OUT);
	sysseek($OUT, $loc3, 0);
	syswrite($OUT, eof_block());
    }
}
close_handles(); #clear handles

truncate($outfile1, $ipc_share1) if($outfile1 && -f $outfile1); #make file proper size
truncate($outfile2, $ipc_share2) if($outfile2 && -f $outfile2); #make file proper size
truncate($outfile3, $ipc_share3) if($outfile3 && -f $outfile3); #make file proper size

#io hack
if($io_dir){
    move($outfile1, $io_files{OUTFILE1}) if($outfile1);
    move($outfile2, $io_files{OUTFILE2}) if($outfile2);
    move($outfile3, $io_files{OUTFILE3}) if($outfile3);
}

#zero length third file so just delete it
if($outfile3){
    unlink($outfile3) if($ipc_share3 == 0 || ($gzip && $ipc_share3 == 28));
}

#validate that bgzip files are not corrupt
if($validate && $gzip){
    foreach my $file ($outfile1, $outfile2, $outfile3){
	next if(!$file || ! -f $file);

	#make sections of right size
	my $size = (stat($file))[7];
	my $count = ceil($size/33554432); #32Mb sections
	my $window = ceil($size/$count);
	my @sections = (1..$count);
	share(@sections); #share is non destructive for forks::shared

	#validate end of file
	my $IN = get_handle('<', $file);
	binmode($IN);
	seek($IN, -28, 2); #last 28 bytes
	my $last;
	my $stat = read($IN, $last, 28, 0);
	die "ERROR: Output file is missing the BGZF EOF marker: $file\n" if($stat != 28 || $last ne eof_block());
	close_handles();

	#run all threads on data
	my @threads;
	my %param = (SIZE => $size, COUNT => $count, WINDOW => $window, SECTIONS => \@sections, INFILE => $file);
	for(my $i = 1; $i < $cpus; $i++){
	    my $thr = threads->new({'context' => 'scalar'}, \&thread_validate, \%param);
	    push(@threads, $thr);
	}
	thread_validate(\%param); #let main process join in
	
	#wait on threads
	foreach my $thr (@threads){
	    $thr->join; #wait on threads
	}
    }
}

#finished
exit(0);

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------

#get header from bam
sub get_header {
    my $file = shift;
    my $validate = shift;

    #declare variables
    my $header;
    my $header_off;
    my $header_block;
    
    #open input file
    my $IN = get_handle('<', $file);
    binmode($IN);
    
    #get header
    inflate(readblock($IN), \$header_block, $validate);
    
    my $magic = substr($header_block, 0, 4);
    $header_off += 4;
    die "ERROR: Magic string mismatch. This does not appear to be a bam file\n"
	if($magic ne "BAM\1");
    
    #grow header block to needed size
    my $l_text = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    while((length($header_block)-$header_off < $l_text+4)){
	inflate(readblock($IN), \$header_block, $validate);
    }
    my $text = substr($header_block, $header_off, $l_text);
    $header_off += $l_text;
    $text =~ s/\0$//g; #remove null padding
    $header->{text} = $text;
    
    my $n_ref = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    $header->{n_ref} = $n_ref;

    #get each reference
    for(my $i = 0; $i < $n_ref; $i++){
	while(length($header_block)-$header_off < 4){ #grow if needed
	    inflate(readblock($IN), \$header_block, $validate);
	}
	my $l_name = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	
	while(length($header_block)-$header_off < $l_name+4){ #grow if needed
	    inflate(readblock($IN), \$header_block, $validate);
	}
	my $name = substr($header_block, $header_off, $l_name);
	$header_off += $l_name;
	$name =~ s/\0$//g; #remove null padding
	$header->{ref}[$i]{name} = $name;
	
	my $l_ref = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	$header->{ref}[$i]{l_ref} = $l_ref;
    }
    substr($header_block,0,$header_off,''); #chop off header from block
    $header_off = tell($IN); #make offset be position of first block after header

    $header->{offset} = $header_off; #offset of first block imediately following the header
    $header->{block}  = $header_block; #piece of alignment accidentally stuck inside header block

    return $header;
}

#convenience method to see if I am in a thread or not
{my $tid; BEGIN{$tid = threads->tid if exists $INC{'threads.pm'}};
sub is_thread {
    my $is_thread = 1 if(defined($tid) && $tid != threads->tid);
    return $is_thread;
}}

#main fastq conversion done here by each process
sub thread_run {
    my $param = shift;

    #fix thread signalling
    if(is_thread){
	binmode(STDIN);
	binmode(STDOUT);
	select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
	select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

	$SIG{'__DIE__'} = sub {print STDERR  "$_[0]\n"; exit(255)};
	$SIG{INT} = sub {exit(130)};
    }

    #load parameters
    my $size = $param->{SIZE};
    my $count = $param->{COUNT};
    my $window = $param->{WINDOW};
    my $file = $param->{INFILE};
    my $q64      = $param->{Q64};
    my $restore  = $param->{RESTORE};
    my $validate = $param->{VALIDATE};
    my $outfile1 = $param->{OUTFILE1};
    my $outfile2 = $param->{OUTFILE2};
    my $outfile3 = $param->{OUTFILE3};
    my $sections = $param->{SECTIONS};

    #declare variables
    my @forks;
    my @buffer; #paired read buffer
    my $header;
    my $header_off;
    my $header_block;
    my $paired;
    my ($fq1_ref,$fq2_ref,$fq3_ref,$fq4_ref);
    my ($gz1_ref,$gz2_ref,$gz3_ref,$gz4_ref);
    {my $fastq1 = ''; $fq1_ref = \$fastq1;} #set reference
    {my $fastq2 = ''; $fq2_ref = ($outfile2) ? \$fastq2 : $fq1_ref} #set reference
    {my $fastq3 = ''; $fq3_ref = \$fastq3;} #set reference
    {my $fastq4 = ''; $fq4_ref = \$fastq4;} #set reference
    {my $gz1 = ''; $gz1_ref = \$gz1;} #set reference
    {my $gz2 = ''; $gz2_ref = \$gz2;} #set reference
    {my $gz3 = ''; $gz3_ref = \$gz3;} #set reference
    {my $gz4 = ''; $gz4_ref = \$gz4;} #set reference
    
    #open input file
    my $IN = get_handle('<', $file);
    binmode($IN);
    
    #get header
    inflate(readblock($IN), \$header_block, $validate);
    
    my $magic = substr($header_block, 0, 4);
    $header_off += 4;
    die "ERROR: Magic string mismatch. This does not appear to be a bam file\n"
	if($magic ne "BAM\1");
    
    #grow header block to needed size
    my $l_text = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    while((length($header_block)-$header_off < $l_text+4)){
	inflate(readblock($IN), \$header_block, $validate);
    }
    my $text = substr($header_block, $header_off, $l_text);
    $header_off += $l_text;
    $text =~ s/\0$//g; #remove null padding
    $header->{text} = $text;
    
    my $n_ref = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    $header->{n_ref} = $n_ref;
    
    #get each reference
    for(my $i = 0; $i < $n_ref; $i++){
	while(length($header_block)-$header_off < 4){ #grow if needed
	    inflate(readblock($IN), \$header_block, $validate);
	}
	my $l_name = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	
	while(length($header_block)-$header_off < $l_name+4){ #grow if needed
	    inflate(readblock($IN), \$header_block, $validate);
	}
	my $name = substr($header_block, $header_off, $l_name);
	$header_off += $l_name;
	$name =~ s/\0$//g; #remove null padding
	$header->{ref}[$i]{name} = $name;
	
	my $l_ref = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	$header->{ref}[$i]{l_ref} = $l_ref;
    }
    substr($header_block,0,$header_off,''); #chop off header from block
    $header_off = tell($IN); #make offset be position of first block after header

    #now read alignment sections
    while(my $section = shift @$sections){
	my $start = ($section-1) * $window;
	my $end   = $section * $window;
	$start = $header_off if($start < $header_off);
	$end = $size if($end > $size); #don't go after end

	#adjust to start of next block
	seek_next_block($IN, $start, 0);
	
	#get first block in section
	my $block = ($section == 1) ? $header_block : ''; #piece of alignment stuck at end of header
	inflate(readblock($IN), \$block, $validate);

	#find first alignment
	my $offset = tell_next_alignment(\$block, $header);
	while(! defined($offset) && tell($IN) < $start+65536 && tell($IN) < $size){ #get past short and empty blocks
	    inflate(readblock($IN), \$block, $validate);
	    $offset = tell_next_alignment(\$block, $header);
	}
	die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($block));
	substr($block, 0, $offset) = ''; #chop off leading data

	#get all blocks in section
	if($gzip){
	    while(tell($IN) < $end){
                inflate(readblock($IN), \$block, $validate);
                bam2fastq(\$block, $fq1_ref, $fq2_ref, $fq4_ref, \@buffer, $q64, $restore);

		#compress
		#max ISIZE is 65536
                deflate(\ substr($$fq1_ref, 0, 65536, ''), $gz1_ref) while(length($$fq1_ref) >= 65536);
                deflate(\ substr($$fq2_ref, 0, 65536, ''), $gz2_ref) while(length($$fq2_ref) >= 65536);
                deflate(\ substr($$fq4_ref, 0, 65536, ''), $gz4_ref) while(length($$fq4_ref) >= 65536);
	    }
	}
	else{
	    while(tell($IN) < $end){
                inflate(readblock($IN), \$block, $validate);
                bam2fastq(\$block, $fq1_ref, $fq2_ref, $fq4_ref, \@buffer, $q64, $restore);
	    }
	}

	#get last block and adjust split terminating line
	my $tail = '';
	inflate(readblock($IN), \$tail, $validate);
	$offset = tell_next_alignment(\$tail, $header);
	while(! defined($offset) && tell($IN) < $end+65536 && tell($IN) < $size){ #get past short and empty blocks
	    inflate(readblock($IN), \$tail, $validate);
	    $offset = tell_next_alignment(\$tail, $header);
	}
	die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($tail));
	$block .= substr($tail, 0, $offset) if($offset); #chop off trailing data

	#convert to fastq
	bam2fastq(\$block, $fq1_ref, $fq2_ref, $fq4_ref, \@buffer, $q64, $restore);
	
	#compress
	if($gzip){
	    #max ISIZE is 65536
	    deflate(\ (substr($$fq1_ref, 0, 65536, '')), $gz1_ref) while(length($$fq1_ref));
	    deflate(\ (substr($$fq2_ref, 0, 65536, '')), $gz2_ref) while(length($$fq2_ref));
	    deflate(\ (substr($$fq4_ref, 0, 65536, '')), $gz4_ref) while(length($$fq4_ref));
	}
	else{
	    $gz1_ref = $fq1_ref;
	    $gz2_ref = $fq2_ref;
	    $gz4_ref = $fq4_ref;
	}

	#minimal validation
	if(!$outfile3){
	    $paired |= 1 if(length($$gz4_ref));
	    $paired |= 2 if(length($$gz1_ref) || length($$gz2_ref));
	    if($paired & 1 && $outfile2){
		die "ERROR: The input BAM contains single end reads. You indicated\n".
		    "that you expect paired end reads by supplying option -fq2, but\n".
		    "you failed to specify a file for single end reads using -fq3.\n";
	    }
	    if($paired & 1 && $paired & 2){
		die "ERROR: The input BAM contains both paired and single end reads,\n".
		    "but you failed to specify a file for single end reads using -fq3\n";
	    }
	}

	print_async([$gz1_ref, $gz2_ref, $gz4_ref],
		    [$outfile1, $outfile2, $outfile3], 1) if(length($$gz1_ref) >= 33554432 ||
							     length($$gz2_ref) >= 33554432 ||
							     length($$gz4_ref) >= 33554432);
    }

    #compress remainder
    if($gzip){
	#max ISIZE is 65536)
	deflate(\ (substr($$fq1_ref, 0, 65536, '')), $gz1_ref) while(length($$fq1_ref));
	deflate(\ (substr($$fq2_ref, 0, 65536, '')), $gz2_ref) while(length($$fq2_ref));
	deflate(\ (substr($$fq4_ref, 0, 65536, '')), $gz4_ref) while(length($$fq4_ref));
    }
    else{
	$gz1_ref = $fq1_ref;
	$gz2_ref = $fq2_ref;
	$gz4_ref = $fq4_ref;
    }

    #print leftovers
    print_async([$gz1_ref, $gz2_ref, $gz4_ref],
		[$outfile1, $outfile2, $outfile3], 0);

    #clear filehandles
    close_handles();

    #thread stores buffer to avoid high memory usage
    if(is_thread){
	my (undef, $tfile) = tempfile("bam2fastq\_$$\_XXXXX", CLEANUP => 0, DIR => File::Spec->tmpdir);
	store(\@buffer, $tfile);
	return $tfile;
    }

    return \@buffer;
}

our %HANDLES;
sub get_handle {
    my $mode = shift;
    my $file = shift || 'PIPE';

    if(!$file){
	return \*STDOUT if($mode eq '>' || $mode eq '>>' || $mode eq '+<');
	return \*STDIN if($mode eq '<');
    }
    if(!$HANDLES{$mode}{$file}){
	open(my $FH, $mode, $file);
	$HANDLES{$mode}{$file} = $FH
    }

    return $HANDLES{$mode}{$file};
}
sub close_handles {
    close($_) foreach(map {values %$_} values %HANDLES);
    undef %HANDLES;
}

sub print_async {
    my $refs     = shift;
    my $outfiles = shift;
    my $async    = shift;

    #lock and change file size if necessary
    my $ipc_handle = $IPCHANDLE;
    if($async){
	return unless($ipc_handle->shlock(LOCK_EX|LOCK_NB));
    }
    else{
	sleep 0.1 while(!$ipc_handle->shlock(LOCK_EX));
    }

    my @locs;
    for(my $i = 0; $i < @$refs; $i++){
	my $ref = $refs->[$i];
	my $outfile = $outfiles->[$i];
	my $ipc_share = ($outfile) ? $IPC4FILE{$outfile} : $IPC4FILE{PIPE};

	#output first pair or interleaved pairs
	$locs[$i] = $$ipc_share;
	$$ipc_share += length($$ref);
	if(!$outfile){
	    syswrite(STDOUT, $$ref);
	    $$ref = ''; #reset
	}
	else{
	    truncate($outfile, $$ipc_share) if(-f $outfile && $$ipc_share > $osize); #grow
	}
    }
    $ipc_handle->shunlock;

    #write data once lock is no longer needed
    for(my $i = 0; $i < @$refs; $i++){
	my $ref = $refs->[$i];
	my $outfile = $outfiles->[$i];
	my $ipc_share = ($outfile) ? $IPC4FILE{$outfile} : $IPC4FILE{PIPE};

	#output first pair or interleaved pairs
	next if(! length($$ref));

	my $loc = $locs[$i];
	next if(!$outfile);
	
	my $OUT = get_handle('+<', $outfile);
	binmode($OUT);
	sysseek($OUT, $loc, 0);
	syswrite($OUT, $$ref);
	$$ref = ''; #reset
    }
}
    
{#predeclare values for efficienty
my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen,
    $bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,@seq,$seq,@qual,$qual);
my ($tag, $val_type, $value, $val_len, $sub_type, $sub_len);
my ($cl_ref,$cl_query,$match);
my ($substr_off, $data_len, $block_offset);
my ($data, $fq1_ref, $fq2_ref, $buffer, $id, $q64);

my @seq_index;
BEGIN {@seq_index = qw(= A C M G R S V T W Y H K D B N);} #force it to exist

sub bam2fastq {
    $data = shift; #ref
    $fq1_ref = shift; #ref
    $fq2_ref = shift; #ref
    $fq4_ref = shift; #ref
    $buffer  = shift; #ref
    $q64     = shift || 0;

    #collect and process bam alignments
    $data_len = length($$data);
    $block_offset = 0;
    while($block_offset < $data_len - 4){ #continue until near end of block
	$substr_off = $block_offset;
	
	($block_size) = unpack('l<', substr($$data, $substr_off, 4));
	$substr_off += 4;
	last if($substr_off+$block_size > $data_len); #current alignment is spit across next block
	$block_offset = $substr_off+$block_size; #end of current alignmnet
	
	($refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen)
	    = unpack('l<l<L<L<l<l<l<l<', substr($$data, $substr_off, 32));
	$substr_off += 32;
	
	#bin_mq_nl processing into sub values
	#$bin = ($bin_mq_nl>>16);
	#$mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1 
	$l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1
	
	#flag_nc processing into sub values
	$flag = ($flag_nc>>16); #flag is probably not usful here
	$n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
	next if($flag & 2816); #skip secondary, supplemental, and vendor failed alignments
	
	#get read name
	$read_name = substr($$data, $substr_off, $l_read_name);
	$substr_off += $l_read_name;
	chop($read_name); #removes trailing null
	
	#process cigar into sub values
	if($n_cigar_op){
	    @cigar = unpack("L<"x$n_cigar_op, substr($$data, $substr_off, $n_cigar_op*4));
	    $substr_off += $n_cigar_op*4;
	    @cigar = map {[($_ & 15), ($_>>4)]} @cigar; #mask is (1<<16)-1
	    $cl_ref = 0;
	    $cl_query = 0;
	    foreach (@cigar){
		$match = ((0x3C1A7 >> ($_->[0] << 1)) & 3);
		$cl_query += $_->[1] if($match & 1);
		$cl_ref += $_->[1] if($match & 2);
	    }
	    if($pos >= 0 && $cl_query != $l_seq){ #hard masking.  I've lost sequence
		warn "WARNING: $read_name has lost sequence because of hard masking of reads\n";
	    }
	}
	
	#process seq string (slow)
	my $count = int(($l_seq+1)/2);
	$seq = substr($$data, $substr_off, $count);
	$substr_off += $count;
	convert_seq(\$seq, $l_seq, ($flag & 16));
	
	#process quality string (slow)
	my $qual = substr($$data, $substr_off, $l_seq);
	$substr_off += $l_seq;
	convert_qual(\$qual, $l_seq, ($flag & 16), $q64);

	#process tags to get original quality values
	if($restore){
	    while($substr_off < $block_offset){
		$tag =  substr($$data, $substr_off, 2);
		$substr_off += 2;
		$val_type = substr($$data, $substr_off, 1);
		$substr_off += 1;

		$val_len = 0;
		if   ($val_type eq 'A'){ $val_len = 1 }
		elsif($val_type eq 'c'){ $val_len = 1 }
		elsif($val_type eq 'C'){ $val_len = 1 }
		elsif($val_type eq 's'){ $val_len = 2 }
		elsif($val_type eq 'S'){ $val_len = 2 }
		elsif($val_type eq 'i'){ $val_len = 4 }
		elsif($val_type eq 'I'){ $val_len = 4 }
		elsif($val_type eq 'f'){ $val_len = 4 }
		elsif($val_type eq 'd'){ $val_len = 8 }
		elsif($val_type eq 'Z' || $val_type eq 'H'){
		    $val_len = (index($$data,"\0",$substr_off) - $substr_off)+1; #plus 1 for null
		}
		elsif($val_type eq 'B'){
		    $sub_type = substr($$data, $substr_off, 1);
		    $substr_off += 1;
		    $sub_len = unpack('l<', substr($$data, $substr_off, 4));
		    $substr_off += 4;
		    
		    if   ($sub_type eq 'c'){ $val_len = 1*$sub_len }
		    elsif($sub_type eq 'C'){ $val_len = 1*$sub_len }
		    elsif($sub_type eq 's'){ $val_len = 2*$sub_len }
		    elsif($sub_type eq 'S'){ $val_len = 2*$sub_len }
		    elsif($sub_type eq 'i'){ $val_len = 4*$sub_len }
		    elsif($sub_type eq 'I'){ $val_len = 4*$sub_len }
		    elsif($sub_type eq 'f'){ $val_len = 4*$sub_len }
		}

		if($tag eq $restore){
		    die "ERROR: Wrong datatype for restored qaulity value\n" if($val_type ne 'Z');
		    die "ERROR: Restored quality value does not match sequence length\n" if($val_len-1 != $l_seq);
		    $value = substr($$data, $substr_off, $val_len-1); #-1 to ignore null
		    $substr_off += $val_len;
		    last; #short circuit loop
		}
		else{ #skip past value since it is not the right one
		    $substr_off += $val_len;
		}
	    }

	    #restore value
	    if(!$value){
		warn "WARNING: Original quality value not found for $read_name\n";
	    }
	    else{
		fix_qual(\$value, $l_seq, ($flag & 16), $q64) if(($flag & 16) || $q64);
		$qual = $value;
	    }
	}

	#group pairs together
	if(!($flag & 1)){ #read not paired
	    $$fq4_ref .= "\@$read_name\n$seq\n+\n$qual\n";
	}
	else{ #pairs
	    $id = ($flag & 64) ? $refID : $next_refID;
	    if(exists($buffer->[$id+1]{$read_name})){
                if($flag & 64){ #first in pair
                    $$fq1_ref .= "\@$read_name/1\n$seq\n+\n$qual\n";
		    $$fq2_ref .= delete($buffer->[$id+1]{$read_name})->[0];
                }
                else{
                    $$fq1_ref .= delete($buffer->[$id+1]{$read_name})->[0];
		    $$fq2_ref .= "\@$read_name/2\n$seq\n+\n$qual\n";
                }
	    }
	    else{
		if($flag & 64){ #first in pair
		    $buffer->[$id+1]{$read_name} = ["\@$read_name/1\n$seq\n+\n$qual\n", 1];
		}
		else{
		    $buffer->[$id+1]{$read_name} = ["\@$read_name/2\n$seq\n+\n$qual\n", 2];
		}
	    }
        }
    }
    
    $$data = substr($$data, $block_offset); #remove what I processed off of the block

    return;
}}

sub seek_next_block {
    my $FH = shift; #filehandle
    my $start = shift; #where to start looking from
    my $wence = shift || 0;

    #go to position
    seek($FH, $start, $wence);
    my $pos = tell($FH);
    my $size = (stat($FH))[7]; #file size

    #not enough space left for a block
    return seek($FH, $size, 0) if($size-$start < 28);

    #read in some data to check
    my $data;
    read($FH, $data, 65536, 0); #max block

    #find block keys in data
    my $offset = 0;
    my $key1 = pack('H8', '1f8b0804'); #4 byte static header key
    my $key2 = pack('H12', '060042430200'); #key1, 6 non-static bytes, then key2 (6 bytes)
    while(1){
	#match first static key
	$offset = index($data, $key1, $offset);
	return seek($FH, $size, 0) if($offset == -1);
	
	#match second static key
	my $offset2 = index($data, $key2, $offset);
	return seek($FH, $size, 0) if($offset2 == -1);

	#second match must be 10 bytes from first offset (4 byte key1 and 6 byte gap)
	if($offset2-$offset == 10){
	    last;
	}
	else{
	    $offset += 4; #jump over key1 match and try again
	}
    }
    $pos += $offset; #adjust file position with offset

    return seek($FH, $pos, 0);
}

#assumes file position is set to start of block
sub readblock {
    my $FH = shift;

    #read BGZF header
    my $data;
    my $stat = read($FH, $data, 18, 0);
    return undef if($stat == 0);
    die "ERROR: Failure to read BAM stream\n" if($stat != 18);

    #validate header
    my ($id1,$id2,$cm,$flg,$mtime,$xfl,$os,
	$xlen, $si1,$si2,$slen,$bsize) = unpack('CCCCVCCvCCvv', $data);
    die "ERROR: Does not appear to be a BGZF file\n"
	if($id1 != 31 || $id2 != 139 || $cm != 8 || $flg != 4 ||
	   $xlen != 6 || $si1 != 66 || $si2 != 67 || $slen != 2);
    
    #read compression block and footer
    my $c_data_size = $bsize-$xlen-19; #compression block size
    $stat = read($FH, $data, $c_data_size + 8, length($data)); #the +8 is for the footer 
    die "ERROR: Could not read compression block\n"
	if($stat != $c_data_size + 8);

    return \$data;
}

sub inflate {
    my $buffer = shift; #reference
    my $ref    = shift; #reference
    my $validate = shift;

    return if(!$buffer);

    my ($i_obj, $stat) = inflateInit(-WindowBits => -15, -Bufsize => 65536);
    die "ERROR: Failed to create zlib inflation object with status: $stat\n" if($stat);

    substr($$buffer, 0, 18, '');
    $$ref .= scalar($i_obj->inflate($buffer));

    #validate the length and crc32
    if($validate){
	die "ERROR: Trailing garbage in compression block\n" if(length($$buffer) != 8);

        my ($crc32, $isize) = unpack('L<L<', $$buffer);
	
        if($isize != $i_obj->total_out()){ #size does not match
            die "ERROR: The expected ISIZE of the uncompressed block does not match\n";
        }
        if($crc32 != crc32(substr($$ref, -$isize, $isize))){ #crc32 does not match
            die "ERROR: The expected CRC32 of the uncompressed block does not match\n";
        }
    }

    return;
}

sub deflate {
    my $buffer = shift; #reference
    my $ref = shift; #reference

    my ($d_obj, $stat) = deflateInit(-WindowBits => -15, -Bufsize => 65536, -Level => 4);
    die "ERROR: Failed to create zlib deflation object with status: $stat\n" if($stat);

    my $offset = length($$ref);
    $$ref .= pack('H36', '1f8b08040000000000ff0600424302000000'); #header (BSIZE=0)
    $$ref .= scalar($d_obj->deflate($buffer)).scalar($d_obj->flush); #compressed data
    $$ref .= pack('V', crc32($buffer)); #CRC32
    $$ref .= pack('V', $d_obj->total_in()); #ISIZE
    substr($$ref, $offset+16, 2, pack('v', $d_obj->total_out+6+19)); #set final BSIZE

    return;
}

#tells the virtual offset of the next alignment record
sub tell_next_alignment {
    my $data = shift; #ref
    my $header = shift;

    return if(!$$data);

    #predeclare values for efficienty
    my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,
	$tlen,$bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,@seq,@qual);
    my ($substr_off, $l_ref, $match, $cl_ref, $cl_query);
    my $n_ref = $header->{n_ref} if($header);

    #find correct virtual offset
    for(my $i = 0; $i < length($$data)-44; $i++){
	$substr_off = $i; #shift 1 byte to the right each time through loop
	
	($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,
	 $next_refID,$next_pos,$tlen) = unpack('l<l<l<L<L<l<l<l<l<', substr($$data, $substr_off, 36));
	$substr_off += 36;
	
	next unless(44 <= $block_size && $block_size < 65536); #reasonable assumption?
	next unless( 1 <= $l_seq && $l_seq < 10000); #reasonable assumption?
	
	#bin_mq_nl processing into sub values
	$bin = ($bin_mq_nl>>16);
	$mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1 
	$l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1
	next unless($bin <= 37450); #reasonable assumption?
	
	#flag_nc processing into sub values
	$flag = ($flag_nc>>16); #flag is probably not usful here
	$n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
	
	#get read name
	$read_name = substr($$data, $substr_off, $l_read_name);
	$substr_off += $l_read_name;
	next unless(substr($read_name, -1) eq "\0"); #last character must be null
	
	#process cigar into sub values
	@cigar = unpack("L<"x$n_cigar_op, substr($$data, $substr_off, $n_cigar_op*4));
	$substr_off += $n_cigar_op*4;
	@cigar = map {[($_ & 15), ($_>>4)]} @cigar; #mask is (1<<16)-1
	$cl_ref = 0;
	$cl_query = 0;
	foreach (@cigar){
	    $match = ((0x3C1A7 >> ($_->[0] << 1)) & 3);
	    $cl_query += $_->[1] if($match & 1);
	    $cl_ref += $_->[1] if($match & 2);
	}
	next unless(!@cigar || $cl_query == $l_seq); #reasonable assumption?
	next unless(reg2bin($pos, $pos+$cl_ref) == $bin || ($pos == -1 && $bin == 0));

	#removed because TLEN is niether consistent nor calculable without both reads
	#next unless($refID != $next_refID || $pos < 0 ||
	#	    $next_pos < 0 || $next_pos - ($pos+$cl_ref) == $tlen);
	
	#additionally validate contig and positions info (necessary?)
	if($header && $header->{text}){
	    next unless(-1<= $refID && $refID < $n_ref);
	    $l_ref = ($refID != -1) ? $header->{ref}[$refID]{l_ref} : 0;
	    next unless(-1 <= $pos && $pos < $l_ref);
	    next unless(-1 <= $next_refID && $next_refID < $n_ref);
	    $l_ref = ($next_refID != -1) ? $header->{ref}[$next_refID]{l_ref} : 0;
	    next unless(-1 <= $next_pos && $next_pos < $l_ref);
	}
	
	#process seq string
	
	#process quality string
	
	#ignore aux values for now
	
	return $i; #seems to be a valid bam line
    }

    return;
}

sub eof_block {
    return pack('H56', '1f8b08040000000000ff0600424302001b0003000000000000000000');
}

sub thread_validate {
    my $param = shift;

    #fix thread signalling
    if(is_thread){
        binmode(STDIN);
        binmode(STDOUT);
        select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
        select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

        $SIG{'__DIE__'} = sub {print STDERR  "$_[0]\n"; exit(255)};
        $SIG{INT} = sub {exit(130)};
    }

    #load parameters
    my $size = $param->{SIZE};
    my $count = $param->{COUNT};
    my $window = $param->{WINDOW};
    my $file = $param->{INFILE};
    my $sections = $param->{SECTIONS};
    my $validate = 1;

    #open input file
    my $IN = get_handle('<', $file);
    binmode($IN);

    #now read alignment sections
    while(my $section = shift @$sections){
        my $start = ($section-1) * $window;
        my $end   = $section * $window;
        $end = $size if($end > $size); #don't go after end

        #adjust to start of next block
        seek_next_block($IN, $start, 0) if($start > 0);

        #validate all blocks in section
        while(tell($IN) < $end){
            my $block;
            inflate(readblock($IN), \$block, $validate);
        }

        #validate one more (avoids weirdness with bad blocks at boundaries)
        if(tell($IN) < $size){
            my $block;
            inflate(readblock($IN), \$block, $validate);
        }
    }

    #clear any open file handles
    close_handles();

    return;
}

#the C code to implement MPI calls from perl
 __END__
__C__

#include <stdint.h>
#include <math.h>

#define bam1_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

int convert_seq(SV* ref, int len, int reverse) {
    SV* scalar = SvRV(ref);
    uint8_t* seq = (uint8_t*)SvPV_nolen(scalar);

    char* buf = (char*)malloc((len+1)*sizeof(char));
    buf[len] = 0;

    int i;
    if (reverse != 0) {
	for (i = 0; i < len; ++i){
	    buf[len - i - 1] = bam_nt16_rev_table[seq_comp_table[bam1_seqi(seq, i)]];
	}
    }
    else{
	for (i = 0; i < len; ++i) {
	    buf[i] = bam_nt16_rev_table[bam1_seqi(seq, i)];
	}
    }

    sv_setpvn(scalar, buf, len);
    free(buf);
}

int fix_qual(SV* ref, int len, int reverse, int q64) {    
    //note these qualities are already expected to be phred+33
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);

    //char* buf = (char*)malloc((len+1)*sizeof(char));
    //buf[len] = 0;

    //recalculate quality shift of 64
    int i;
    if(q64 == 2){ //recalculate solexa encoding
	int qsol;
	for (i = 0; i < len; ++i){
	    if(qual[i] < 59) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
	    qsol = qual[i]-64;
	    qual[i] = 33 + qsol + 10 * log10(1 + pow(10, qsol/-10)); //change scaling equation
	}
    }
    else if(q64 != 0){ //other illumina encoding
	for (i = 0; i < len; ++i){
	    if(qual[i] < 64) croak("ERROR: Qaulity values are too low to be in Phred+64 format\n");
	    qual[i] -= 31;
	}
    }

    if (reverse != 0) {
	uint8_t b;
	for (i = 0; i < len/2; ++i){
	    if(qual[i] > 75) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	    if(qual[len-1-i] > 75) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");

	    //swap
	    b = qual[i];
	    qual[i] = qual[len-1-i];
	    qual[len-1-i] = b;
	}
	sv_setpvn(scalar, (char*)qual, len);
    }
    else{
	for (i = 0; i < len; ++i) {
	    if(qual[i] > 75) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	}
    }
}

int convert_qual(SV* ref, int len, int reverse, int q64) {    
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);

    char* buf = (char*)malloc((len+1)*sizeof(char));
    buf[len] = 0;

    //recalculate quality shift of 64
    int i;
    if(q64 == 2){ //recalculate solexa encoding
	int qsol;
	for (i = 0; i < len; ++i){
	    if(qual[i] < 26) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
	    qsol = qual[i]-31;
	    qual[i] = qsol + 10 * log10(1 + pow(10, qsol/-10)); //change scaling equation
	}
    }
    else if(q64 != 0){ //other illumina encoding
	for (i = 0; i < len; ++i){
	    if(qual[i] < 31) croak("ERROR: Qaulity values are too low to be in Phred+64 format\n");
	    qual[i] -= 31;
	}
    }

    if (reverse != 0) {
	for (i = 0; i < len; ++i){
	    if(qual[i] > 42) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	    buf[len - i - 1] = 33 + qual[i];
	}
    }
    else{
	for (i = 0; i < len; ++i) {
	    if(qual[i] > 42) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	    buf[i] = 33 + qual[i];
	}
    }

    sv_setpvn(scalar, buf, len);
    free(buf);
}

short reg2bin(int beg, int end){
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
    return 0;
}
