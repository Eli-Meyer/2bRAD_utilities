#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a 2bRAD genotype matrix into the csv input format required by R/qtl.
Usage: $scriptname -i input -t traits -o output
Where:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-t traits	tab-delimited file of data on  traits, as
				"sample1	trait1	...	traitN"
			(note that sample names must match column headers in snps file)
	-m map		tab-delimited file of map positions as
				"marker	  LG	position"
	-o output	a name for the csv formatted output file.
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

# -- define variables from input
$infile = $opt_i;
$traitfile = $opt_t;
$mapfile = $opt_m;
$outfile = $opt_o;

# -- read and store trait data 
open(IN, $traitfile);
$sampno = $popno = 0;
while (<IN>)
	{
	chomp;
	$sampno++;
	@cols = split("\t", $_);
	if ($sampno == 1) {@headers = @cols; next;}
	for ($c=1; $c<@cols; $c++)
		{
		$th{$cols[0]}{$headers[$c]} = $cols[$c];
		}
	}
@samps = sort(keys(%th));
$nosamps = @samps;
@traits = @headers[1..$c];
$notraits = $c - 1 ;
print STDERR "Input file $traitfile contains $nosamps samples, with data on $notraits traits.\n";

# -- build a hash of genotypes for each sample and locus
open(GT, $infile);
$rowno = 0;
while(<GT>)
	{
	chomp;
	$rowno++;
	@cols = split("\t", $_);
	if ($rowno eq 1) {@header = @cols; next;}
	$ncols = @cols;
	%allh = ();
	for ($n=2; $n<$ncols; $n++)
		{
		@alls = split(" ", $cols[$n]);
		foreach $a (@alls) {unless ($a eq 0) {$allh{$a}++;}}
		}
	@sall = sort(keys(%allh));
	$sname = $cols[0].".".$cols[1];
	$snph{$sname}++;
	for ($n=2; $n<$ncols; $n++)
		{
		if ($cols[$n] eq 0)
			{
			$gh{$header[$n]}{$sname} = "-";
			}
		elsif ($cols[$n] eq $sall[0])
			{
			$gh{$header[$n]}{$sname} = "A";
			}
		elsif ($cols[$n] eq $sall[1])
			{
			$gh{$header[$n]}{$sname} = "B";
			}
		elsif ($cols[$n] =~ / /)
			{
			$gh{$header[$n]}{$sname} = "H";
			}
		}	
	}
close(GT);

# -- read and store map
open(IN, $mapfile);
$sampno = 0;
while (<IN>)
	{
	chomp;
	$sampno++;
	@cols = split("\t", $_);
	if (!exists($snph{$cols[0]})) {next;}
	if ($sampno == 1) {@headers = @cols; next;}
	$mh{$cols[0]}{"LG"} = $cols[1];	
	$mh{$cols[0]}{"pos"} = $cols[2];	
	}
@marks = sort(keys(%mh));
$nomarks = @marks;
print STDERR "Input file $mapfile contains $nomarks markers.\n";

# -- write output in CSV format for R/qtl
open (OUT, ">$outfile");
print OUT join("\t",@traits), join("\t", @marks), "\n";

print OUT "\t" x $notraits;
foreach $m (@marks)
	{
	print OUT $mh{$m}{"LG"}, "\t";
	}
print OUT "\n";

print OUT "\t" x $notraits;
foreach $m (@marks)
	{
	print OUT $mh{$m}{"pos"}, "\t";
	}
print OUT "\n";

foreach $s (@samps)
	{
	if ($s eq "Locus") {next;}
	foreach $t (@traits)
		{
		if(exists($th{$s}{$t}))
			{print OUT $th{$s}{$t}, "\t";}
		}	
	foreach $m (@marks)
		{
		print OUT $gh{$s}{$m}, "\t";		
		}
	print OUT "\n";
	}
