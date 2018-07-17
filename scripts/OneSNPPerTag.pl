#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Selects a single SNP from each tag in a matrix or genotypes or allele counts. 
Chooses the locus with the least missing data. 
Usage: $scriptname -i input <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
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
getopts('i:o:p:m:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_p) {$pvar = $opt_p;} else {$pvar = "n";}	# whether to print loci passing filter
if ($pvar eq "y" && !$opt_o) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_o) {$ovar = $opt_o;} else {$ovar = "null";}	# name for output file
if ($opt_m) {$mode = $opt_m;} else {$mode = "g";}	# name for output file
if ($mode ne "g" && $mode ne "a") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $opt = $pvar;

if ($mode eq "g")
{
open(IN, $opt_i);
open(OUT, ">$ovar");
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	$mdcount = 0;
	foreach $a (@cols[2..$ncols]) {if ($a eq "0") {$mdcount++;}}
	$ch{$cols[0]}{$cols[1]} = $mdcount;
	}
close(IN);
foreach $tag (sort(keys(%ch)))
	{
	%tagh = %{$ch{$tag}};
	@sorted = sort{$tagh{$a}<=>$tagh{$b}}(keys(%tagh));
	$useh{$tag} = $sorted[0];
	}
open(IN, $opt_i);
$rowcount = 0;
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	if($useh{$cols[0]} ne $cols[1]) {next;}
	if ($opt eq "y") {print OUT $_, "\n";}
	$kept++;
	}
close(IN); close(OUT);
print STDERR $rowcount-1, " loci\n";
print STDERR $kept, " SNPs selected.\n";
}

elsif ($mode eq "a")
{
open(IN, $opt_i);
open(OUT, ">$ovar");
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	$mdcount = 0;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if (($cols[$a] eq "0" && $cols[$a+1] eq "0")||($cols[$a] eq "NA" && $cols[$a+1] eq "NA"))
			{$mdcount++;}
		}
	$ch{$cols[0]}{$cols[1]} = $mdcount;
	}
close(IN);

foreach $tag (sort(keys(%ch)))
	{
	%tagh = %{$ch{$tag}};
	@sorted = sort{$tagh{$a}<=>$tagh{$b}}(keys(%tagh));
	$useh{$tag} = $sorted[0];
	}

$rowcount = 0;
open(IN, $opt_i);
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	if($useh{$cols[0]} ne $cols[1]) {next;}
	if ($opt eq "y") {print OUT $_, "\n";}
	$kept++;
	}
close(IN);
print STDERR $rowcount-1, " loci\n";
print STDERR $kept, " SNPs selected.\n";
}

if (-e "null")
	{
	system("rm null");
	}
