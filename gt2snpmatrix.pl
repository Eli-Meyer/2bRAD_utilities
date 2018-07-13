#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a genotype matrix (loci x samples) to a snp matrix, as described in the manual
for the R package diveRsity. This snp matrix can be converted to genepop format
using diveRsity's snp2gen function.
Usage: $scriptname -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. SNP matrix format, input for snp2gen.
USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

open(IN, $opt_i);
open(KEY, ">snpkey.$opt_o");
open(OUT, ">$opt_o");
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	$nom = @cols;

	if ($rowcount==1)
		{
		print OUT "SNP_ID", "\t";
		print OUT join("\t", @cols[2..$nom]), "\n";
		next;
		}
	$tagcount++;
	print OUT "SNP", $tagcount, "\t";
	print KEY $tagcount, "\t", $cols[0], "\t", $cols[1], "\n";
	for ($a=2; $a<$nom; $a++)
		{
		if ($cols[$a] =~ / /)
			{
			$cols[$a] =~ s/ //;
			print OUT $cols[$a], "\t";
			}
		elsif ($cols[$a] eq 0)
			{
			print OUT "--", "\t";
			}
		else
			{
			print OUT $cols[$a], $cols[$a], "\t";	
			}
		}
	print OUT "\n";
	}

