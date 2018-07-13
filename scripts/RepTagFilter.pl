#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Excludes tags containing too many SNPs, suggesting repetive regions of the genome
Usage: $scriptname -i input -n max_snps <OPTIONS>
Required arguments:
	-i input	name of the SNP input file, a matrix of genotypes or allele counts
			(see -m for format)
			note: the script assumes the input only includes polymorphic loci
	-n max_snps	all SNPs from tags containing more than this number of SNPs will be excluded
Options:
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
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:n:p:m:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_n || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_p) {$pvar = $opt_p;} else {$pvar = "n";}	# whether to print loci passing filter
if ($pvar eq "y" && !$opt_o) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_o) {$ovar = $opt_o;} else {$ovar = "null";}	# name for output file
if ($opt_m) {$mode = $opt_m;} else {$mode = "g";}	# name for output file
if ($mode ne "g" && $mode ne "a") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
# note that choice of mode does nothing in this script
# its only included for consistency with other filtering scripts

my $maxsnps = $opt_n;
my $opt = $pvar;

open(IN, $opt_i);
open(OUT, ">$ovar");

while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	$ch{$cols[0]}++;
	}
close(IN);

$rowcount = 0;
open(IN, $opt_i);
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	if ($ch{$cols[0]}>$maxsnps) {$toomany++; next;}
	if ($opt eq "y") {print OUT $_, "\n";}
	}
close(IN); close(OUT);
if (-e "null")
	{
	system("rm null");
	}
print STDERR $rowcount-1, " loci\n";
print STDERR $toomany, " loci with more than ", $maxsnps, " SNPs per tag\n";
print STDERR $rowcount-1-$toomany, " loci with no more than more than ", $maxsnps, " SNPS per tag\n";
