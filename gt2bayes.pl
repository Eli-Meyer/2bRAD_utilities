#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a 2bRAD genotype matrix into the input format required by BayeScan.
Usage: $scriptname input pop.file output
Where:
	input:		tab-delimited genotype matrix (columns=samples, rows=loci) 
			produced by RADGenotyper.pl or NFGenotyper.pl
	pop.file:	a tab-delimited text file showing which population each
			sample was drawn from. Formatted as: SampleName \"\\t\" PopName \"\\n\"
			Note -- make sure that sample names in this file are exactly identical
			to those shown in the first row of the genotype matrix.
	output: 	a name for the BayeScan formatted output file.
USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
# use this block if checking for executable dependencies
# copy the block and edit to check for additional Perl modules required by the script
#$mod1="File::Which";
#unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
#use File::Which;

# use this block and edit to check for executables required by the script
#$dep1 = "blastx";
#unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}

# -- define variables from input
$infile = $ARGV[0];
$popfile = $ARGV[1];
$outfile = $ARGV[2];

# -- read and store sample-population association provided by user in popfile
open(IN, $popfile);
$sampno = $popno = 0;
while (<IN>)
	{
	chomp;
	$sampno++;
	($samp, $pop) = split("\t", $_);
	$poph{$pop}++;
	$samph{$samp} = $pop;
	}
@pops = sort(keys(%poph));
$popno = @pops;
print "Input file contains ", $sampno, " samples from ", $popno, " populations (", join(", ", @pops), ").\n";

# -- count the number of loci in input file
$locno = `wc -l $infile` - 1;

# -- initialize output file with locus and population information
open(OUT, ">$outfile");
print OUT "[loci]=", $locno, "\n\n";
print OUT "[populations]=", $popno, "\n";

# -- identify the colums containing samples from each population
open(GT, $infile);
$rowno = 0;
while(<GT>)
	{
	$rowno++;
	if ($rowno>1) {last;}
	chomp;
	@cols = split("\t", $_);
	foreach $pop (@pops)
		{
		$colno = -1;
		foreach $c (@cols)
			{
			$colno++;
			if (exists($samph{$c}))
				{
				if ($samph{$c} eq $pop)
					{
					$keyh{$pop}{$colno}++;
					}
				}
			}
		%tmppoph = %{$keyh{$pop}};
		@usevec = sort{$a<=>$b}(keys(%tmppoph));
#		print $pop, "\t", "@usevec", "\n";
		}
	}

# -- for each population, read through SNP matrix and output BayeScan formatted data for that population
$popcount = 0;
foreach $pop (@pops)
	{
	$popcount++;
	%tmppoph = %{$keyh{$pop}};
	@usevec = sort{$a<=>$b}(keys(%tmppoph));
	print OUT "\n", "[pop]=", $popcount, "\n";
	open(GT, $infile);
	$rowno = 0; $locnum = 0;
	while(<GT>)
		{
		$rowno++; if ($rowno eq 1) {next;}
		$locnum++;
		print OUT $locnum, "\t";
		chomp;
		@cols = split("\t", $_);
		$ncols = @cols;
		@subcols = @cols[@usevec];
# -- count missing data for each population and locus		
		$mdl = $gtl = 0;
		foreach $s (@subcols)
			{
			if ($s eq 0) {$mdl++;}
			else	{$gtl++;}	
			}
		$numgenes = $gtl*2;
		print OUT $numgenes, "\t";
# -- count alleles observed at each locus
		%allh = ();
		foreach $c (@cols[2..$ncols])
			{
			if ($c eq 0) {next;}
			@alls = split(" ", $c);
			foreach $a (@alls) {$allh{$a}++;}
			}
		@obsalls = sort(keys(%allh));
		$numall = @obsalls;
		print OUT $numall, "\t";
#		print $cols[0], "\t", $cols[1], "\t", "$numall", "\n";

# -- count alleles observed for each population and locus
		%popallh = ();
		foreach $s (@subcols)
			{
			if ($s eq 0) {next;}
			if ($s !~ / /)
				{
				$popallh{$s}++; $popallh{$s}++; 
				}
			else
				{
				@popalls = split(" ", $s);
				foreach $p (@popalls)
					{
					$popallh{$p}++;
					}
				}
			}
		foreach $a (@obsalls)
			{
			if (exists($popallh{$a}))	{print OUT $popallh{$a}, " ";}
			else	{print OUT "0 ";}
			}
		print OUT "\n";
		}
	}
