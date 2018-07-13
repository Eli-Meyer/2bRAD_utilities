#!/usr/bin/perl
# written by E Meyer, eli.meyer@oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Determines SNP genotypes from nucleotide frequencies. Input file contains nucleotide frequencies 
from multiple samples. Output file lists the genotypes called from those frequencies. 
Usage: $scriptname -i input -o output <OPTIONS>
Required arguments:
	-i input	Input file, tab delimited text file of nucleotide frequencies (output from NtFrequences.pl)
			column = tag, column 2 = locus, column 3 = reference genotype
			subsequent columns = nucleotide frequences in each sample, as A/C/G/T (e.g. 0/0/10/12)
	-o output	A name for the output file (tab delimited text)
Options:
	-c coverage	Minimum coverage required to determine genotypes. Lower coverage loci wil be discarded.
			Default: 10
	-e ends		y: exclude terminal positions in alignments where errors may arise during ligation. (default)
			n: do not exclude terminal positions. 
	-m method	"nf": nucleotide frequencies (classic method; the default). 
			  This method determines genotypes directly from nucleotide frequencies, using thresholds
			  defined by the user. If minor allele frequency (MAF) <= x at a locus, the genotype is called 
			  homozygous for the major allele at that locus. If MAF >= n, genotype is called heterozygous.
			  Genotypes are not called at intermediate MAF (if x > MAF > n) where errors are likely.  
			"pgf": NF informed by population genotype frequencies. (an update on the classic method)
			  This method first identifies valid alleles at each locus based on their frequency in the 
			  population (the two most common alleles with frequencies >=y times in >= q individuals),
			  then applies relaxed nucleotide frequency thresholds for those alleles (using y instead
			  of n for valid alleles). 
			"bgc" = Bayesian Genotype Caller
			  This method calls the BGC software, which implements a maximum-likelihood (ML) method for 
			  calling genotypes that incorporates prior population-level information on genotype 
			  frequencies and error rates from a genotype-frequency estimator. For more details see 
			  (Maruki & Lynch, [doi: 10.1534/g3.117.039008], and cite that paper if using this option.	
	Options for method "nf" or "pgf":
	-x max_MAF	Maximum frequency of the minor allele you're willing to ignore and call the position 
			homozygous for the major allele (0-1). Default: 0.01
	-n min_MAF	Minimum frequency of the minor allele you're willing to accept as evidence of 
			heterozygosity, and call the locus heterozygous (0-1). Default: 0.25
	-r min_reads	Because low frequencies translate into 1 or fewer reads at low coverage, the script
			also imposes a minimum read number for detection of heterozygotes. (default: 2) 
	Options for method "pgf":
	-y frequency	Minimum frequency a second allele must be detected to be considered valid.
			(default: 0.05)
	-q samples	Each allele must present in at least q samples to be considered valid.
			(default: 2)	
	Options for method "bgc":
	-p p-value	Critical p-value for the chi-square polymorphism test (BGC)
			(default: 0.05)
	-v maxcov	Coverage at which the pipeline switches from BGC (for low coverage data) to
			HGC (for high coverage data). Default=80 (i.e. HGC above 80). 
Examples:
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab		# basic usage
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab -c 20 		# increase coverage threshold
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab -m bgc		# use Bayesian Genotype Caller
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab -m pgf -y 0.05	# use population method
USAGE
if ($#ARGV < 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SeqIO;
$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use File::Which;

# get variables from input
getopts('i:o:e:c:m:x:n:y:q:p:r:v:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
$ftabfile = $opt_i;	# input filename
$outfile = $opt_o;	# output filename
if ($opt_m)	{$method = $opt_m;} else	{$method = "nf";}
if ($method ne "bgc" && $method ne "nf" && $method ne "pgf") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($method eq "bgc")
	{
	$dep1 = "GFE";
	unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}
	$dep1 = "BGC";
	unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}
	$dep1 = "HGC";
	unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}
	}
if ($opt_c) {$mincov = $opt_c;} else {$mincov = 10;}		# default min coverage to be genotyped
if ($opt_x) {$maxhom = $opt_x;} else {$maxhom = 0.01;}		# default maximum MAF ignored and called homozygous for major allele
if ($opt_n) {$minhet = $opt_n;} else {$minhet = 0.25;}		# default minimum MAF considered heterozygous
if ($opt_r) {$minreads = $opt_r;} else {$minreads = 2;}		# default minimum MA read depth for heterozygous calls
if ($opt_y) {$minall = $opt_y;} else {$minall = 0.05;}		# default minimum coverage per allele to define valid alleles
if ($opt_q) {$minsamp = $opt_q;} else {$minsamp = 2;}		# default minimum number of samples to define valid alleles
if ($opt_p) {$bgcp = $opt_p;} else {$bgcp = 0.05;}		# default critical p value for polymorphism test (BGC)
if ($opt_e) {$eopt = $opt_e;} else {$eopt = "y";}		# exclude terminal positions from alignments
if ($opt_v) {$maxcov = $opt_v;} else {$maxcov = 80;}		# max coverage for BGC

system("date");

# Nucleotide Frequency Method - the classic approach
if ($method eq "nf")
{
# make a big hash of all base calls
print "Reading $ftabfile ...\n";
my %bh;
@boa = qw{A C G T};
open(FTAB, $ftabfile);
$rowno = 0;
while(<FTAB>)
	{chomp;
	$rowno++;
	@cols = split("\t", $_);
	if ($rowno == 1) 
		{
		@header = @cols;
		next;
		}
	for ($i=3; $i<@cols; $i++)
		{
		@bases = split("/", $cols[$i]);
		$jno = 0;
		foreach $j (@bases)
			{
			$bh{$cols[0]}{$cols[1]}{$header[$i]}{$boa[$jno]} = $j;
			$jno++;
			}
		}
	}
print "Finished. \n";
print $rowno-1, " loci with nucleotide frequency data. \n";
system("date");
close(FTAB);		

@samlabs = @header[3..@cols-1];

# loop through all loci and samples, applying genotyping rules using threshold defined above
# in this method, genotyping is based directly on nucleotide frequency thresholds
print "Determining genotypes from nucleotide frequencies ...\n";
open(OUT,">$outfile");
print OUT "scaffold", "\t", "position", "\t", join ("\t", @samlabs), "\n";
@sites = sort(keys(%bh));
$countgt = $gtlocno = $sampno = 0;
foreach $s (@sites)
	{
	%siteh = %{$bh{$s}};
	@loca = sort{$a<=>$b}(keys(%siteh));
	foreach $l (@loca)
		{
		$sgt = 0;
		if ($l == 1 | $l == 36)
			{
			if($eopt eq "y")
				{
				next;
				}
			}
		%loch = %{$siteh{$l}};
		@gta = ();
		foreach $s (@samlabs)
			{
			%samh = %{$loch{$s}};
			%tmph = ();
			foreach $b (@boa) {$tmph{$b} = $samh{$b};}
			@saa = sort{$tmph{$b}<=>$tmph{$a}}(keys(%tmph));
			$majmin = $tmph{$saa[0]}+$tmph{$saa[1]};
			$minreject = $majmin * $maxhom;
			$gti = ""; $rati = "";
#			print $saa[0], ":", $tmph{$saa[0]}, ",", $saa[1], ":", $tmph{$saa[1]}, "[", $majmin, "]";
# first conditions where we arent interested in or cannot compare nucleotide frequencies
# exclude loci with >2 alleles observed, and loci with low coverage
			if (($tmph{$saa[2]}>=$minreject)||($majmin<$mincov))
				{$gti = 0;}	 
# then conditions where we can compare nucleotide frequencies
			else
				{
				$rati = $tmph{$saa[1]}/($majmin);
#				print "(", $rati, ")";
				if ($rati<$maxhom) {$gti = $saa[0];}
				elsif ($rati>=$minhet) {@rsa = sort(@saa[0..1]); $gti = "@rsa";}
				elsif (($rati >= $maxhom && $rati < $minhet)||($tmph{$saa[1]} < $minreads)) {$gti = 0;}
				}
#			print "\t";
			push @gta, $gti;
			if ($gti =~ /[ACGT]/) {$sgt++;}
			}
#		print join("\t", @gta), "\n";
		$countgt += $sgt;
		if ($sgt > 0) 
			{
			$gtlocno++;
			print OUT $s, "\t", $l, "\t", join("\t", @gta), "\n";
			}
		}
	}
print "Finished. \n";
$sampno = @samlabs;
print "Called genotypes at $gtlocno loci in $sampno samples ($countgt genotypes altogether)\n";
system("date");
}

# Nucleotide Frequencies Informed by Population Genotype Frequencies 
# in this method, valid alleles are defined at each locus by counting the number of times each 
# allele occurs at high enough coverage in enough samples. Each sample is then genotyped only for those
# alleles based on observed numbers and user defined nucleotide frequency thresholds. The increased
# confidence in the valid alleles provided by population analysis makes it possible to relax MAF
# threshold for heterozygous calls, determining genotypes at loci with intermediate MAFs that could 
# otherwise not be determined directly from nucleotide frequencies. 
if ($method eq "pgf")
{
# make a big hash of all base calls
print "Reading $ftabfile ...\n";
my %bh;
@boa = qw{A C G T};
open(FTAB, $ftabfile);
$rowno = 0;
while(<FTAB>)
	{chomp;
	$rowno++;
	@cols = split("\t", $_);
	if ($rowno == 1) 
		{
		@header = @cols;
		next;
		}
	for ($i=3; $i<@cols; $i++)
		{
		@bases = split("/", $cols[$i]);
		$jno = 0;
		foreach $j (@bases)
			{
			$bh{$cols[0]}{$cols[1]}{$header[$i]}{$boa[$jno]} = $j;
			$jno++;
			}
		}
	}
print "Finished. \n";
print $rowno-1, " loci with nucleotide frequency data. \n";
system("date");
close(FTAB);		

@samlabs = @header[3..@cols-1];

# loop through all loci and samples, applying genotyping rules using thresholds defined above.
print "Determining genotypes from population genotype frequencies and sample nucleotide frequencies ...\n";
open(OUT,">$outfile");
print OUT "scaffold", "\t", "position", "\t", join ("\t", @samlabs), "\n";
@sites = sort(keys(%bh));
foreach $s (@sites)
	{
	%siteh = %{$bh{$s}};
	@loca = sort{$a<=>$b}(keys(%siteh));
# filter out terminal bases
	foreach $l (@loca)
		{
		$sgt = 0;
		if ($l == 1 | $l == 36)
			{
			if($eopt eq "y")
				{
				next;
				}
			}
		%loch = %{$siteh{$l}};
		@gta = (); %vah = ();
# first identify valid alleles
		foreach $s (@samlabs)
			{
			%samh = %{$loch{$s}};
			%tmph = (); $covi = 0; 
			foreach $b (@boa) {$tmph{$b} = $samh{$b}; $covi += $samh{$b};}
			@saa = sort{$tmph{$b}<=>$tmph{$a}}(keys(%tmph));
			if ($covi > 0)
				{
				if ((($tmph{$saa[0]}/$covi) >= $minall)&&($tmph{$saa[0]} >= $minreads)) {$vah{$saa[0]}++;}
				if ((($tmph{$saa[1]}/$covi) >= $minall)&&($tmph{$saa[1]} >= $minreads)) {$vah{$saa[1]}++;}
				}
			}
		@tmp = sort{$vah{$b}<=>$vah{$a}}(keys(%vah));
		%sva = (); $nva = 0;
		foreach $t (@tmp)
			{
			if ($vah{$t} >= $minsamp) {$sva{$t}++; $nva++;}
			}
#		print $s, "\t", $l, "\t", $nva, " valid alleles\n";
# filter out positions with more than two alleles
# then genotype samples based on frequency of those alleles in each sample
# taking into account the population genotype frequencies determined above
		foreach $s (@samlabs)
			{
			%samh = %{$loch{$s}};
			%tmph = ();
			foreach $b (@boa) {$tmph{$b} = $samh{$b};}
			@saa = sort{$tmph{$b}<=>$tmph{$a}}(keys(%tmph));
			$majmin = $tmph{$saa[0]}+$tmph{$saa[1]};
			$minreject = $majmin * $maxhom;
			$gti = ""; $rati = "";
#			print $saa[0], ":", $tmph{$saa[0]}, ",", $saa[1], ":", $tmph{$saa[1]}, "[", $majmin, "]";
# first conditions where we arent interested in or cannot compare nucleotide frequencies
# exclude loci with >2 alleles observed, and loci with low coverage
			if (($tmph{$saa[2]}>=$minreject)||($majmin<$mincov))
				{
				$gti = 0;
#				print "r";
				}		
# then conditions where we can compare nucleotide frequencies
			else
				{
				$rati = $tmph{$saa[1]}/($majmin);
#				print "(", $rati, ")";
				if ($rati<$maxhom) {$gti = $saa[0];}
				elsif ($rati>=$minhet) {@rsa = sort(@saa[0..1]); $gti = "@rsa";}
				elsif (($rati >= $maxhom && $rati < $minhet)||($tmph{$saa[1]} < $minreads)) 
					{
					if(exists($sva{$saa[1]}))	# if minor allele in this sample is one 
						{			# of top 2 in pgf, handle this way
						if (($rati >= $minall) && ($tmph{$saa[1]} >= $minreads)) 
							{@rsa = sort(@saa[0..1]); $gti = "@rsa";}
						else
							{$gti = 0;}
						}
					else				# if minor allele is not one of the top
						{			# 2 in pgf handle it this way
						$gti = 0;
						}
					}
				}
#			print "\t";
			push @gta, $gti;
			if ($gti ne 0) {$sgt++;}
			}
#		print join("\t", @gta), "\t", $sgt, "\n";
		$countgt += $sgt;
		if ($sgt > 0) 
			{
			$gtlocno++;
			print OUT $s, "\t", $l, "\t", join("\t", @gta), "\n";
			}
		}
	}
print "Finished. \n";
$sampno = @samlabs;
print "Called genotypes at $gtlocno loci in $sampno samples ($countgt genotypes altogether)\n";
system("date");
}

# Bayesian Genotyper Caller Method
if ($method eq "bgc")
{
# estimate genotype frequencies in population
$gfe_out = $ftabfile.".gfe.out";
print "Running GFE ...\n";
system "GFE -cv $bgcp -mode c -as 1 -min_cov $mincov -max_cov $maxcov -in $ftabfile -out $gfe_out > gfe_log.txt";
print "Finished.\n";
system("date");

# call genotypes at low coverage sites using ML approach informed by prior population estimates
$bgc_out = $ftabfile.".bgc.out";
print "Running BGC ...\n";
system "BGC -as 1 -min_cov $mincov -max_cov $maxcov -in $gfe_out -out $bgc_out > bgc_log.txt";
print "Finished.\n";
system("date");

# call genotypes at high coverage sites using ML approach informed by prior population estimates
$hgc_out = $ftabfile.".hgc.out";
print "Running HGC ...\n";
HGC:
system "HGC -in $ftabfile -out $hgc_out -min_cov $maxcov > hgc_log.txt";
$numrows=`wc -l $hgc_out`;
$numrows =~ s/\s.+//;
chomp($numrows);
if ($numrows eq 1) {goto HGC;}		# because HGC sometimes exits without throwing an error. This catches it.
print "Finished.\n";
system("date");

# integrate output from BGC and HGC
$rownum=0;
open(IN, $hgc_out);
while(<IN>)
	{
	chomp;
	$rownum++;
	@cols = split("\t", $_);
	if ($rownum==1) {@header=@cols; next;}
	$len=@cols; $end=$len-8;
	for ($a=11; $a<$end; $a++)
		{
		unless ($cols[$a] eq "NA")
			{
			$fh{$cols[0]}{$cols[1]}{$header[$a]} = $cols[$a];
			}
		}
	}
close(IN);

open(IN, $bgc_out);
open(OUT, ">integrated_bgc.tab");
$rownum=0;
while(<IN>)
	{
	chomp;
	$rownum++;
	@cols = split("\t", $_);
	if ($rownum==1) {print OUT $_, "\n"; @header=@cols; next;}
	$end=@cols-1;
	for ($a=5; $a<$end; $a++)
		{
		if(exists($fh{$cols[0]}{$cols[1]}{$header[$a]}))
			{
			$cols[$a] = $fh{$cols[0]}{$cols[1]}{$header[$a]};
			}
		}
	print OUT join("\t", @cols), "\n";
	}
close(IN);

# convert output to same format as the other methods
print "Converting output file...\n";
open(OUT, ">$outfile");
foreach $a (qw{0 1 2 3 4})
	{$exc_col{$a}++;}
foreach $a (qw{2 3 4})
	{$exa_col{$a}++;}
open(IN, "integrated_bgc.tab");
$rownum = 0;
while(<IN>)
	{
	$rownum++;
	chomp;
	@cols = split("\t", $_);
	$cno = 0;
	if ($rownum==1)
		{
		foreach $c (@cols)
			{
			unless (exists($exa_col{$cno})) 
				{
				print OUT $c, "\t";
				}
			$cno++;
			}
		print OUT "\n";
		next;
		}

	$nonmd = 0; $cno = 0;
	foreach $c (@cols)
		{
		unless (exists($exc_col{$cno})) 
			{
			if ($c ne "NA")	{$nonmd++;}
			}
		$cno++;
		}
	if ($nonmd < 1) {next;}

	if ($eopt eq "y")
		{
		if ($cols[1] == 1 || $cols[1] == 36)
			{
			next;
			}
		}

	$cno = 0; $gtino = 0;
	foreach $c (@cols)
		{
		if ($cno < 2)
			{
			print OUT $c, "\t";
			}
		unless (exists($exc_col{$cno})) 
			{
			@alls = split("", $c);
			if ($c eq "NA")
				{
				print OUT "0\t";
				}
			elsif ($alls[0] eq $alls[1])
				{
				$countgt++; $gtino++;
				print OUT $alls[0], "\t";
				}
			else
				{
				$countgt++; $gtino++;
				print OUT $alls[0], " " , $alls[1], "\t";
				}
			}
		$cno++;
		}
	if ($gtino>0) {$gtloc++;}
	print OUT "\n";
	}	
print "Finished.\n";
$sampno = @cols - 5;
print "Genotyped $gtloc loci in $sampno samples ($countgt genotypes altogether)\n";
system("date");
}


