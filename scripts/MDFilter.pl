#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Excludes loci containing too many missing data (genotyped in too few samples)
Usage: $scriptname -i input -n min_data <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
	-n min_data	loci that were genotyped in fewer samples than this will be excluded
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
open(IN, $opt_i);
open(OUT, ">$ovar");
$toofew = 0;
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	if ($rowcount==2)
		{
		$maxmiss = $ncols - 2 - $mindat;
		}
	$missi = 0;
	for ($a=2; $a<$ncols; $a++)
		{
		if ($cols[$a]eq "0") {$missi++;}
		}
	if ($missi>$maxmiss) {$toofew++; next;}
	if ($opt eq "y") {print OUT $_, "\n";}
	}
close(IN); close(OUT);

print STDERR $rowcount-1, " loci\n";
print STDERR $toofew, " loci with more than ", $maxmiss, " missing data (i.e. genotyped in fewer than $mindat samples)\n";
print STDERR $rowcount-1-$toofew, " loci with no more than more than ", $maxmiss, " missing data\n";
}

# mode a
elsif ($mode eq "a")
{
open(IN, $opt_i);
open(OUT, ">$ovar");
$toofew = 0;
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	if ($rowcount==2)
		{
		$maxmiss = ($ncols - 4)/2 - $mindat;
		}
	$missi = 0;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if (($cols[$a] eq 0 && $cols[$a+1] eq 0)||($cols[$a] eq "NA" && $cols[$a+1] eq "NA")) {$missi++;}
		}
	if ($missi>$maxmiss) {$toofew++; next;}
	if ($opt eq "y") {print OUT $_, "\n";}
	}
close(IN); close(OUT);
print STDERR $rowcount-1, " loci\n";
print STDERR $toofew, " loci with more than ", $maxmiss, " missing data (i.e. genotyped in fewer than $mindat samples)\n";
print STDERR $rowcount-1-$toofew, " loci with no more than more than ", $maxmiss, " missing data\n";
}
if (-e "null")
        {
        system("rm null");
        }

