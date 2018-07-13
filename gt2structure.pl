#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a genotype matrix (loci x samples) into the appropriate input format for STRUCTURE.
Usage: $scriptname -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. STRUCTURE input format.
USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

$zh{"A"} = 1;
$zh{"C"} = 2;
$zh{"G"} = 3;
$zh{"T"} = 4;
$zh{"0"} = 0;

open(IN, $opt_i);
open(OUT, ">$opt_o");
while(<IN>)
	{
	chomp;
	$rowcount++;
	if ($rowcount==1)
		{
		@noma = split("\t", $_);
		$nom = @noma;
		@noma = @noma[2..$nom];
		next;
		}
	@cols = split("\t", $_);
	$nom = @cols;
	for ($a=2; $a<$nom; $a++)
		{
		$gh{$noma[$a-2]}{$rowcount-1} = $cols[$a];
		}
	$oh{$noma[$a-2]} = $rowcount;
	}

foreach $s (@noma)
#foreach $s (sort(keys(%gh)))
	{
	if ($s eq "") {last;}
	%sh = %{$gh{$s}};
	print OUT $s,"-1\t";
	for ($b=1; $b<$rowcount; $b++)
		{
		$gi = $sh{$b};
		if ($gi =~ /\s/) {$gi =~ s/ .+//;}
		print OUT $zh{$gi}, "\t";
		}
	print OUT "\n";
	print OUT $s,"-2\t";
	for ($b=1; $b<$rowcount; $b++)
		{
		$gi = $sh{$b};
		if ($gi =~ /\s/) {$gi =~ s/.+ //;}
		print OUT $zh{$gi}, "\t";
		}
	print OUT "\n";
	}

