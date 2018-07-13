#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a SNP genotype matrix (loci x samples) produced from 2bRAD genotyping
into the format required for the BLUPF90 family of programs for mixed models
and quantitative genetic analysis. See BLUPF90 manual for details of that format.
Usage: $scriptname -i input -o output
Required arguments:
	-i input	Name of the input file, from CallGenotypes.pl.
			(rows=loci, columns=samples, columns 1 & 2 show tag name and position in tag)
	-o output	A name for the output file, in the format expected by BLUPF90 programs, e.g.
				sample0   02221022511020101020
				sample100 12221222221222200010
USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

open(IN, $opt_i);
open(OUT, ">$opt_o");
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	$nom = @cols;

	if ($rowcount==1)
		{
		@samples = @cols;
		next;
		}
	%ch = ();
	for ($a=2; $a<$nom; $a++)
		{
		@obsall = split(" ", $cols[$a]);
		foreach $o (@obsall) {unless ($o eq 0) {$ch{$o}++;}}
		}
	@saa = sort{$ch{$b}<=>$ch{$a}}(keys(%ch));
	if (@saa > 2) 
		{
		print STDERR "Warning -- > 2 alleles detected at $cols[0] position $cols[1]\n";
		next;
		}
	$maja = $saa[0];
	$mina = $saa[1];
	@recoded = ();
	for ($a=2; $a<$nom; $a++)
		{
		if ($cols[$a] eq 0)	{push @recoded, "5"; $gth{$samples[$a]}{$cols[0]}{$cols[1]} = 5;}
		elsif ($cols[$a] =~ $maja && $cols[$a] =~ $mina)	{push @recoded, "1"; $gth{$samples[$a]}{$cols[0]}{$cols[1]} = 1;}
		elsif ($cols[$a] =~ $maja)	{push @recoded, "2"; $gth{$samples[$a]}{$cols[0]}{$cols[1]} = 2;}
		elsif ($cols[$a] =~ $mina)	{push @recoded, "0"; $gth{$samples[$a]}{$cols[0]}{$cols[1]} = 0;}
		else	{print STDERR "Warning -- foreign allele detected at $cols[0] position $cols[1]\n";}
		}
	}

# check length of sample labels
$maxlen = 0;
foreach $s (sort(keys(%gth)))
	{
	$leni = length($s);
	if ($leni>$maxlen) {$maxlen=$leni;}
	}
$bufflen = $maxlen+1;

# loop through samples and combine genotypes from all loci into a single vector
# output in the format expected by BLUPF90 programs
foreach $s (sort(keys(%gth)))
	{
	%sgth = %{$gth{$s}};
	print OUT $s, " "x($bufflen-length($s));
	foreach $t (sort(keys(%sgth)))
		{
		%tsgth = %{$sgth{$t}};
		foreach $l (sort(keys(%tsgth)))
			{
			print OUT $tsgth{$l};
			}
		}
	print OUT "\n";
	}
