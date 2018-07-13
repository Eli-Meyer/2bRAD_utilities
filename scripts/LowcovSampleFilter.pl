#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Excludes samples with too much missing data (genotypes called at too few loci)
Usage: $scriptname -i input -n min_data <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
	-n min_data	samples in which fewer than this number of loci were genotyped will be excluded
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

my $mindat = $opt_n;
my $opt = $pvar;

# mode g
if ($mode eq "g")
{

open(OUT, ">$ovar");
open(IN, $opt_i); while(<IN>) {$irowcount++;} close(IN);
$maxmiss = $irowcount - 1 - $mindat;

open(IN, $opt_i);
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	if($rowcount==1){next;}
	$ncols = @cols;
	for ($a=2; $a<$ncols; $a++)
		{
		if ($cols[$a]eq "0") {$ch{$a}++;}
		}
	}
close(IN);

$toofew = 0;
open(IN, $opt_i);
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	if($rowcount==1)
		{
		if($opt eq "y") 
			{
			@hdr = @cols;
			print OUT $hdr[0], "\t", $hdr[1], "\t";
			for ($a=2; $a<$ncols; $a++)
				{
				if ($ch{$a}<=$maxmiss) {print OUT $cols[$a], "\t";}
				}
			print OUT "\n";
			} 
		next;
		}
	$ncols = @cols;
	if ($opt eq "y") {print OUT $cols[0], "\t", $cols[1], "\t";}
	for ($a=2; $a<$ncols; $a++)
		{
		if($opt eq "y") {if ($ch{$a}<=$maxmiss) {print OUT $cols[$a], "\t";}}
		}
	if ($opt eq "y") {print OUT "\n";}
	}
close(IN); close(OUT); 
if (-e "null")
	{
	system("rm null");
	}

foreach $k (sort(keys(%ch)))
	{
	if ($ch{$k}>$maxmiss) {$toofew++;}
	}

print STDERR $ncols-2, " samples\n";
print STDERR $toofew, " samples with more than ", $maxmiss, " missing data\n";
print STDERR $ncols-2-$toofew, " samples with no more than more than ", $maxmiss, " missing data (i.e. at least $mindat loci genotyped)\n";

}

# mode a
elsif ($mode eq "a")
{

open(OUT, ">$ovar");
open(IN, $opt_i); while(<IN>) {$irowcount++;} close(IN);
$maxmiss = $irowcount - 1 - $mindat;

open(IN, $opt_i); 
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	if($rowcount==1){next;}
	$ncols = @cols;
	$nsamps = ($ncols-4)/2;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if (($cols[$a] eq "0" && $cols[$a+1] eq "0")||($cols[$a] eq "NA" && $cols[$a+1] eq "NA")) 
			{
			$ch{$a}++;
			}
		}
	}
close(IN);

$toofew = 0;
open(IN, $opt_i);
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	if($rowcount==1)
		{
		if($opt eq "y") 
			{
			@hdr = @cols;
			print OUT $hdr[0], "\t", $hdr[1], "\t", $hdr[2], "\t", $hdr[3], "\t";
			for ($a=4; $a<$ncols; $a+=2)
				{
				if ($ch{$a}<=$maxmiss) {print OUT $cols[$a], "\t", $cols[$a+1], "\t";}
				}
			print OUT "\n";
			} 
		next;
		}
	$ncols = @cols;
	if ($opt eq "y") {print OUT $cols[0], "\t", $cols[1], "\t", $cols[2], "\t", $cols[3], "\t";}
	for ($a=4; $a<$ncols; $a+=2)
		{
		if($opt eq "y") {if ($ch{$a}<=$maxmiss) {print OUT $cols[$a], "\t", $cols[$a+1], "\t";}}
		}
	if ($opt eq "y") {print OUT "\n";}
	}
close(IN);

foreach $k (sort(keys(%ch)))
	{
	if ($ch{$k}>$maxmiss) {$toofew++;}
	}

print STDERR $nsamps, " samples\n";
print STDERR $toofew, " samples with more than ", $maxmiss, " missing data\n";
print STDERR $nsamps-$toofew, " samples with no more than more than ", $maxmiss, " missing data (i.e. at least $mindat loci called)\n";

}
if (-e "null")
	{
	system("rm null");
	}
