#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Counts nucleotide frequencies at each locus in a 2bRAD sequence data set.
Usage: $scriptname -i input -r reference -o <OPTIONS>
Required arguments:
	-i input		input alignments, SAM format
	-r reference		reference used to generate the input alignments, FASTA format
	-o output		a name for the output file (tab delimited text)
Options:
	-c coverage		loci with lower coverage are discarded (default: 3)
USAGE
if ($#ARGV < 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:r:o:c:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_r || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_c) {$thd = $opt_c;} else {$thd = 3;}	
$infile = $opt_i;
$reffile = $opt_r;
$outfile = $opt_o;
# load reference into memory
system("date");
print "beginning analysis of nucleotide frequencies in $infile...\n";
print "Reading reference...\n";
open(REF,$reffile);
while(<REF>)
	{
	chomp;
	if ($_ =~ /\>/)
		{
		$_ =~ s/\>//;
		$refi = $_;
		}
	elsif ($_ =~ /[ACGT]/)
		{
		@bits = split("", $_);
		$bcount = 0;
		foreach $b (@bits)
			{
			$bcount++;
			$rrh{$refi}{$bcount} = $b;
			}
		}
	}
close(REF);
print "Finished loading reference. into memory.\n";

# define positions to be ignored (recognition sites, which are invariant)
foreach $p (13, 14, 15, 22, 23, 24) {$exch{$p}++;}		# for AlfI

# read in sam file output line by line
print "Reading alignments and counting nucleotides at each position...\n";
open(IN, $infile);
my %nfh;
while(<IN>)
	{
	if ($_ =~ /^@/) {next;}
	chomp;
	$rowi = $_;
	$rowi =~ s/^>//;
	@cols = split("\t", $rowi);
	$ncols = @cols;

	if ($cols[2] eq "*") {next;}
	$rawmap++;

# -- build a key for each read matching each position to a reference position
	$cig = $cols[5];
	@ciga = split(/(?<=\D)/, $cig);
	@codea = (); $refloc = $cols[3];
	foreach $c (@ciga)
		{
		if ($c =~ /M/)
			{
			$c =~ s/M//;
			for ($d=1; $d<=$c; $d++) {push @codea, $refloc; $refloc++;}
			}
		elsif ($c =~ /D/)
			{
			$c =~ s/D//;
			for ($d=1; $d<=$c; $d++) {$refloc++;}
			}
		elsif ($c =~ /I/)
			{
			$c =~ s/I//;
			for ($d=1; $d<=$c; $d++) {push @codea, 0;}
			}
		}

# -- count each base toward nucleotide frequencies at the corresponding reference position
	@reada = split("",$cols[9]);

	if (@reada ne @codea) {print "Something is wrong.\n @reada \n @codea\n"; exit;}
	$readpos = -1;
	foreach $b (@reada)
		{
		$readpos++;
		if ($codea[$readpos] eq 0) {next;}
		if ($codea[$readpos] eq "") {next;}
		$nfh{$cols[2]}{$codea[$readpos]}{$b}++;
		}
	}
print "Finished reading and counting alignments.\n";

# -- write output, applying coverage threshold per locus
print "Writing output...\n";
open(OUT, ">$outfile");
print OUT "Tag\tLocus\tRef\tA\tC\tG\tT\n";
@nbins = qw{A C G T};
$rownum = 0;
foreach $ref (sort(keys(%nfh)))
	{
	%refh = %{$nfh{$ref}};
	foreach $loc (sort{$a<=>$b}(keys(%refh)))
		{
		%loch = %{$refh{$loc}};
		$rownum++;
		$covi = 0;
		foreach $n (@nbins)
			{
			if (exists($loch{$n})) {$covi += $loch{$n};}
			}
		if(exists($exch{$loc})) {next;}
		if ($covi < $thd) {next;}
		print OUT $ref, "\t", $loc, "\t", $rrh{$ref}{$loc}, "\t";
		foreach $n (@nbins)
			{
			if (exists($loch{$n})) {print OUT $loch{$n}, "\t";}
			else	{print OUT "0\t";}
			}
		print OUT "\n";
		}
	}
print "Finished with analysis of nucleotide frequencies.\n";
system("date");
print "\n";
