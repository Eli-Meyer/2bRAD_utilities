#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Excludes loci containing too many alleles.
Usage: $scriptname -i input -n max_alleles <options>
Required arguments:
	-i input	name of the input file, a matrix of genotypes.
			Input format: rows=loci and columns=samples.
                        row 1 = column label, column 1 = tag, column 2 = position
                        subsequent columns contain genotypes for each sample
	-n max_alleles	maximum number of alleles allowed. Loci with more than this 
			number of allels will be excluded. 
Options:
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:n:p:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_n || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_p) {$pvar = $opt_p;} else {$pvar = "n";}	# whether to print loci passing filter
if ($pvar eq "y" && !$opt_o) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_o) {$ovar = $opt_o;} else {$ovar = "null";}	# name for output file

my $maxall = $opt_n;
my $opt = $pvar;

open(IN, $opt_i);
open(OUT, ">$ovar");
$toomany = 0;
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print OUT $_, "\n";} next;}
	@cols = split("\t", $_);
	%allh = ();
	foreach $c (@cols[2..@cols])
		{
		if ($c =~ /0/) {next;}
		@alls = split(" ", $c);
		foreach $a (@alls) {$allh{$a}++;}
		}
	@aa = sort(keys(%allh));
	if (@aa > $maxall) 
		{
		$toomany++; 
		next;
		}
	if ($opt eq "y") {print OUT $_, "\n";}
	}
close(IN); close(OUT);
print STDERR $rowcount-1, " loci\n";
print STDERR $toomany, " loci with more than ", $maxall, " alleles\n";
print STDERR $rowcount-1-$toomany, " loci with no more than ", $maxall, " alleles\n";

if (-e "null")
	{
	system("rm null");
	}
