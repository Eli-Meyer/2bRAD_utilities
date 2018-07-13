#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a genotype matrix (loci x samples) to a PHYLIP-formatted alignment.
Usage: $scriptname -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. PHYLIP alignment format.
USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

$ch{"A C"} = "M"; $ch{"C A"} = "M";
$ch{"A G"} = "R"; $ch{"G A"} = "R";
$ch{"A T"} = "W"; $ch{"T A"} = "W";
$ch{"C G"} = "S"; $ch{"C G"} = "S";
$ch{"C T"} = "Y"; $ch{"C T"} = "Y";
$ch{"G T"} = "K"; $ch{"G T"} = "K";

$window = 10;

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
		if ($cols[$a] =~ / /)
			{
			$cols[$a] = $ch{$cols[$a]};
			}
		if ($cols[$a] eq 0)
			{
			$cols[$a] = "-";
			}
		$gh{$noma[$a-2]}{$rowcount-1} = $cols[$a];
		}
	}


print OUT " " x 4, $nom, " " x 3, $rowcount-1, "\n";
$bsize = 50;
$noblocks = int(($rowcount-1)/$bsize+0.5);
foreach $a (1..$noblocks)
	{

foreach $s (sort(keys(%gh)))
	{
	%sh = %{$gh{$s}};
	if ($a==1) 
		{print OUT $s;
		$lens = length($s);
		$diff = $window-$lens;
		print OUT " " x $diff;
		}
	else	{print OUT " " x $window;}
	print OUT " ";
	for ($c=0; $c<5; $c++)
		{
	for ($b=0; $b<10; $b++)
		{
		$cni = $b + ($c * 10) + ($bsize * ($a-1));
		print OUT $sh{$cni+1};
		}
		print OUT " ";
		}
	print OUT "\n";
	}
	print OUT "\n";
	}
