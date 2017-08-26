#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Converts a 2bRAD genotype matrix into the csv input format required by R/qtl.
Usage: $scriptname snps traits output
Where:
	snps:		tab-delimited genotype matrix (columns=samples, rows=loci) 
			produced by RADGenotyper.pl or NFGenotyper.pl
	traits:		tab-delimited file of data on  traits
			arranged as sample/trait1/..traitN
			note that sample names must match column headers in snps file
	map:		tab-delimited file of map positions as
			marker/LG/position
	output: 	a name for the csv formatted output file.
USAGE
if ($#ARGV != 3 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- define variables from input
$infile = $ARGV[0];
$traitfile = $ARGV[1];
$mapfile = $ARGV[2];
$outfile = $ARGV[3];

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
