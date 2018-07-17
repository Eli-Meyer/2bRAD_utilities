#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Excludes loci containing too few classes of genotypes or numbers of alleles (keeps polymorphic loci).
Usage: $scriptname -i input <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
Options:
	-n genotypes	minimum number of unique genotypes (for -m g) or alleles (for -m a) 
			required to consider a locus polymorphic (default=2)
	-m mode		g=genotypes (default). 
			  Input file contains genotypes from individuals.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          subsequent columns contain genotypes for each sample
			a=allele counts.
			  Input file contains allele counts from pooled samples.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          column 3 = major allele, column 4 = minor allele
                          subsequent pairs of columns contain allele counts 
			  (major then minor) for each sample
	-v min_cov	(for -m a) minimum coverage required to consider an allele present
			(default=2)
	-s min_samp	minimum number of samples in which an allele must be present to be counted.
			(default=1)
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:n:p:s:v:m:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_n) {$minall = $opt_n;} else {$minall = 2;}	# min genotypes or alleles to consider polymorphic
if ($opt_p) {$pvar = $opt_p;} else {$pvar = "n";}	# whether to print loci passing filter
if ($pvar eq "y" && !$opt_o) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_o) {$ovar = $opt_o;} else {$ovar = "null";}	# name for output file
if ($opt_m) {$mode = $opt_m;} else {$mode = "g";}	# mode
if ($mode ne "g" && $mode ne "a") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_v) {$mincov = $opt_v;} else {$mincov = 2;}	# min cov, for -m a
if ($opt_s) {$minfreq = $opt_s;} else {$minfreq = 1;}	# min samps, for -m a
my $opt = $pvar;

# mode g
if ($mode eq "g")
{
open(IN, $opt_i);
open(OUT, ">$ovar");
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	%ahi = ();
	for ($a=2; $a<$ncols; $a++)
		{
		if ($cols[$a] eq "0") {next;}
		$ahi{$cols[$a]}++;
		
		}
	@alla = sort(keys(%ahi));
	$passno = 0;
	foreach $z (@alla)
		{
		if ($ahi{$z} >= $minfreq) {$passno++;}
		}
	$nalla = @alla;
	if ($passno < $minall) {$nopoly++; next;}
	if ($opt eq "y") {print OUT $_, "\n";}
	}
close(IN); close(OUT);
print STDERR $rowcount-1, " loci in input file\n";
print STDERR $nopoly, " loci with fewer than ", $minall, " genotypes\n";
print STDERR $rowcount-1-$nopoly, " polymorphic loci (at least ", $minall, " genotypes in $minfreq or more samples)\n";
}

# mode a
elsif ($mode eq "a")
{
open(IN, $opt_i);
open(OUT, ">$ovar");
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	$gt_no = $alt_no = $alt_freq = 0;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if ($cols[$a] eq "NA") {next;}
		$gt_no++;
		if ($cols[$a+1] >= $mincov) 
			{
			$alt_no++;
			}
		}
	$alt_freq = $alt_no / $gt_no;
	if ($alt_freq eq 0) {$mono++; next;}
	if ($alt_freq < $minfreq) {$toolow++; next;}
	if ($opt eq "y") {print OUT $_, "\n";}
	}
close(IN); close(OUT);
print STDERR $rowcount-1, " loci in input file\n";
print STDERR $toolow, " loci with minor alleles present in fewer than $minfreq of sample\n";
print STDERR $mono, " monomorphic loci\n";
print STDERR $rowcount-1-$toolow-$mono, " polymorphic loci\n";
}

if (-e "null")
        {
        system("rm null");
        }
