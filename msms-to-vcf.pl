#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Math::Random;

use Data::Dumper;


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

perl msms-to-vcf.pl  --msms basic.msms.txt --genotypeErrorRate 0.1 --nocallRate 0.1 --lowQualityFraction 0.1

Description:

conversts msms to a VCF 4.1 file.

";

#global parameters
my ($help);
my $msms;
my %msms_data ;
my $nocallRate = 0;
my $genotypeErrorRate  = 0;
my $lowQualityFraction = 0;

$msms_data{haplotypes} = ();
$msms_data{positions}  = ();

my $opt_success = GetOptions('help'                  => \$help,
			     'nocallRate=s'          => \$nocallRate,
			     'genotypeErrorRate=s'   => \$genotypeErrorRate,
			     'lowQualityFraction=s'  => \$lowQualityFraction,
			     'msms=s'                => \$msms
    );

die $usage if $help || ! $opt_success || ! $msms;

loadMsmsData();

if(scalar (($msms_data{haplotypes}) % 2) != 0){
    print "Fatal: must simulate even number of haplotypes for diploids!\n";
    die ;
}
printVCFHeader();
printVCFGenotypes();
#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
sub loadMsmsData{

    open (my $IN, '<', $msms) or die "Can't open $msms for reading\n$!\n";

    my $lcounter = 0;

    while (my $line = <$IN>) {

	chomp $line;

	if($lcounter == 0){
	    $msms_data{command} = $line;
	}
	if($lcounter == 1){
	    $msms_data{unique} = $line;
	}
	if($lcounter == 4){
	    ($_, $msms_data{segsites}) = split /:/, $line;
	}
	if($lcounter == 5){
	    my @pos = split /:|\s+/, $line;
	    shift @pos;
	    shift @pos;
	    #megabase
	    my @physical;
	    foreach my $p (@pos){
		push @{$msms_data{positions}}, ($p * 100000);
	    }
	}
	if($lcounter > 4 && $line =~ /\/\//){
	    print "Fatal:  More than one MSMS simulation found in the file\n";
	    print "Info :  Found // more than once in the file.           \n";
	}
	if($lcounter > 6 && $line =~ /\s+/){
	    print "Info : Read in all of MSMS simulation \n";
	    last;
	}

	if($lcounter > 5){
	    if($line =~ /^[01]/){
		push @{$msms_data{haplotypes}}, $line;
	    }
	}
	$lcounter++;
    }
    close $IN;
}
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub printVCFHeader {

    my $chromLine ;
    $chromLine .= "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    for(my $i = 0; $i < scalar @{$msms_data{haplotypes}}/2; $i++){
	my $id = sprintf ("%08X", rand(0xffffffff));
	$chromLine .= "\t$id";
    }
    $chromLine .= "\n";
    
    my $header;

    $header .= "##fileformat=VCFv4.1\n";
    $header .= "##source=msms\n";
    $header .= "##msmsCommand=$msms_data{command}\n";
    $header .= "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
    $header .= "##INFO=<ID=HF,Number=1,Type=Float,Description=\"Heterozygous genotype Frequency\">\n";
    $header .= "##FORMAT=<ID=GT,Number=G,Type=String,Description=\"Genotype\">\n";
    $header .= "##FORMAT=<ID=TR,Number=1,Type=String,Description=\"True phased genotype\">\n";
    $header .= "##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Simulated genotype likelihoods AA,AB,BB\">\n";
    $header .= "##FORMAT=<ID=QT,Number=1,Type=String,Description=\"Quality distribution simulated high\|low \">\n";
    $header .= "##FORMAT=<ID=GE,Number=1,Type=Integer,Description=\"Marks genotyping errors 1=true 0=false\">\n";
    $header .= $chromLine;
    print $header;
    
}
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub genotypeLikelihood{

    my $genotype       = shift;
    my $lowQualityFlag = shift;

    my $stick = 1;

    my @logGenotypeLikelihoods;
    my @GL = ('NA', 'NA', 'NA');

    my $alpha = 2;
    my $beta  = 1;

    if($lowQualityFlag != 1){
	$alpha = 100;
	$beta  = 1;
    }

    my $firstBreak =  Math::Random::random_beta(1, $alpha, $beta);
    push @logGenotypeLikelihoods,  log($firstBreak);
    $stick -= $firstBreak;
    my $seconBreak = Math::Random::random_beta(1, 2, 1);
    push @logGenotypeLikelihoods,  log($seconBreak*$stick);
    $stick -= ($seconBreak*$stick);
    push @logGenotypeLikelihoods, log($stick);
    

    my @sortedllGenotyps = sort {$b <=> $a} @logGenotypeLikelihoods;

    if($genotype == 0){
	$GL[0] = $sortedllGenotyps[0];
	$GL[1] = $sortedllGenotyps[1];
	$GL[2] = $sortedllGenotyps[2];
    }
    if($genotype == 1){
	$GL[1] = $sortedllGenotyps[0];
	if(rand() > 0.5){
	    $GL[0] = $sortedllGenotyps[0];
	    $GL[2] = $sortedllGenotyps[1];
	}
	else{
	    $GL[0] = $sortedllGenotyps[1];
	    $GL[2] = $sortedllGenotyps[0];
	}
    }
    if($genotype == 2){
	$GL[2] = $sortedllGenotyps[0];
	$GL[1] = $sortedllGenotyps[1];
	$GL[0] = $sortedllGenotyps[2];
    }

    my $glString ;
    
    $glString .=  sprintf("%.4f", $GL[0]);
    $glString .= ",";
    $glString .=  sprintf("%.4f", $GL[1]);
    $glString .= ",";
    $glString .=  sprintf("%.4f", $GL[2]);

    return $glString;
}

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub printVCFGenotypes{

    my %seen_pos;

    my $snpIndex = -1;

    my @GTS = ("0/0", "0/1", "1/1");

    POS: foreach my $pos (@{$msms_data{positions}}){

	if(exists $seen_pos{$pos}){
	    next POS;
	}

	$seen_pos{$pos} = 1;

	$snpIndex += 1;
	
	my $vcfLine;
	
	$vcfLine .= "chrSim\t";
	$vcfLine .= "$pos\t";
	$vcfLine .= ".\t";
	$vcfLine .= "A\t";
	$vcfLine .= "T\t";
	$vcfLine .= ".\t";
	$vcfLine .= "PASS\t";
	
	my $nalt  = 0;
	my $nref  = 0;
	my $nhet  = 0;
	my $nhomr = 0;
	my $nhoma = 0;
	my $ngeno = 0;

	my $genotypeField;
	


	for(my $hapIndex = 0; $hapIndex < scalar @{$msms_data{haplotypes}} -1; $hapIndex+=2){

	    my $qflag = 0;
	    if(rand() < $lowQualityFraction){
		$qflag = 1;
	    }

	    my $g1 = substr ($msms_data{haplotypes}[$hapIndex   ], $snpIndex, 1);
	    my $g2 = substr ($msms_data{haplotypes}[$hapIndex +1], $snpIndex, 1);
	    
	    my $genotype;
	    my $truth   ; 
	    my $gls     ;
	    my $qt  = "high";
	    my $ge  = 0;

	    $qt = "low"  if $qflag == 1;
	    
	    $ngeno++    ;

	    if($g1 == 1 && $g2 == 1){
		$truth    = "1|1";
		$genotype = "1/1";
		$gls = genotypeLikelihood(2, $qflag);
		$nalt  += 2;
		$nhoma += 1;
	    }
	    if($g1 == 0 && $g2 == 0){
		$truth     = "0|0";
		$genotype  = "0/0";
		$gls = genotypeLikelihood(0, $qflag);
		$nref  += 2;
		$nhomr += 1;
	    }
	    if($g1 == 1 && $g2 == 0){
		$truth     = "1|0";
		$genotype  = "0/1";
		$gls = genotypeLikelihood(1, $qflag);
		$nref  += 1;
		$nalt  += 1;
		$nhet  += 1;
	    }
	    if($g1 == 0 && $g2 == 1){
		$truth    = "0|1";
		$genotype = "0/1";
		$gls = genotypeLikelihood(1, $qflag);
		$nref  += 1;
		$nalt  += 1;
		$nhet  += 1;
	    }

	    if(rand() < $nocallRate){
		$genotype = "./.";
		$gls = ".,.,.";
		$qt  = ".";
	    }
	    if(rand() < $genotypeErrorRate){
		my $newG = int(rand(2));
		$genotype = $GTS[$newG];
		$gls = genotypeLikelihood(1, 1);
		$ge = 1;
	    }

	    $genotypeField .= "\t$genotype:$truth:$gls:$qt:$ge";
	}
	
	my $af = $nalt / ($nalt + $nref);
	my $hf = $nhet / $ngeno;
	$af = sprintf("%.5f", $af);
	$hf = sprintf("%.5f", $hf);
	
	my $infoField = "AF=$af;HF=$hf;";
	my $formatField = "\tGT:TR:GL:QT:GE";

      
	$vcfLine .= $infoField;
	$vcfLine .= $formatField;
	$vcfLine .= $genotypeField;
	
	print $vcfLine, "\n";

    }
}

