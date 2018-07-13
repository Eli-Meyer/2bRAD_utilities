#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a genotype matrix (loci x samples) to FSTAT format.
Usage: $scriptname -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. FSTAT format. corresponding locus_key 
			and sample_key files are also produced.
USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

$zh{"A"} = "01";
$zh{"C"} = "02";
$zh{"G"} = "03";
$zh{"T"} = "04";
$zh{"0"} = "0";

open(IN, $opt_i);
open(OK, ">locus_key.$opt_o");
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
	$snpno = $rowcount-1;
	print OK "SNP-", $snpno, "\t", $cols[0], "\t", $cols[1], "\n";
	$snpname = "SNP-".$snpno;
	$kh{$rowcount-1} = $snpname;
	for ($a=2; $a<$nom; $a++)
		{
		$gh{$noma[$a-2]}{$rowcount-1} = $cols[$a];
		}
	}


open(OF, ">sample_key.$opt_o");
open(OUT, ">$opt_o");
my $rowno = 0;
print OUT $nom-2, " ", $rowcount-1, " ", 4, " ", 2, "\n";
foreach $k (sort(keys(%kh)))	{print OUT $kh{$k}, "\n";}
foreach $s (sort(keys(%gh)))
	{
	$rowno++;
	%sh = %{$gh{$s}};
	print OF $rowno, "\t", $s, "\n";
	print OUT $rowno," ";
	for ($b=1; $b<$rowcount; $b++)
		{
		$gi = $sh{$b};
		if ($gi !~ /\s/)
			{
			print OUT $zh{$gi}, $zh{$gi}, " ";
			}
		else
			{
			@alla = split(" ", $gi);
			print OUT $zh{$alla[0]}, $zh{$alla[1]}, " ";
			}
		}
	print OUT "\n";
	}
close(OF);
