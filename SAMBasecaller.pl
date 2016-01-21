#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Counts nucleotide frequencies at each locus in a 2bRAD sequence data set.
Usage: $scriptname input.sam reference.fasta cov_threshold output.tab
Where:
	input.sam:		input alignments, SAM format
	reference.fasta:	reference used to generate the input alignments, FASTA format
	cov_threshold:		loci with lower coverage are discarded
	output.tab:		a name for the output file (tab delimited text)
USAGE
if ($#ARGV != 3 || $ARGV[0] eq "-h") {print "\n", $usage, "\n"; exit;}

# define variables from input arguments
$infile = $ARGV[0];
$ref = $ARGV[1];
$thd = $ARGV[2];
$outfile = $ARGV[3];

# define positions to be exluded from basecalling (recognition sites, which are invariant)
foreach $p (13, 14, 15, 22, 23, 24) {$exch{$p}++;}		# for AlfI

# read in sam file output line by line
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

# -- check which strand the read matches
	$strand = 0;
	if ($cols[1] eq 16)	{$strand = -1;}
	else	{$strand = 1;}

# -- count each base toward nucleotide frequencies at the corresponding reference position
	@reada = split("",$cols[9]);

#	print $cols[0], "\t", $cols[2], "\t", "$cig", "\t", $strand, "\n";
#	print join("\t", @codea), "\n", join("\t", @reada), "\n";

	if (@reada ne @codea) {print "Something is wrong.\n @reada \n @codea\n"; exit;}
#	print "@reada", "\n";
	$readpos = -1;
	foreach $b (@reada)
		{
		$readpos++;
#		print $readpos, "\t", $codea[$readpos], "\t", $b, "\n";
		if ($codea[$readpos] eq 0) {next;}
		if ($codea[$readpos] eq "") {next;}
		$nfh{$cols[2]}{$codea[$readpos]}{$b}++;
#		if ($strand eq 1) {$nfh{$cols[2]}{$codea[$readpos]}{$b}++;}
#		elsif ($strand eq -1) {$nfh{$cols[2]}{$codea[$readpos]}{$rch{$b}}++;}
		}
	}

# -- write output, applying coverage threshold per locus
open(OUT, ">$outfile");
print OUT "\tTag\tLocus\tA\tC\tG\tT\tN\n";
@nbins = qw{A C G T N};
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
		print OUT $rownum, "\t", $ref, "\t", $loc, "\t";
		foreach $n (@nbins)
			{
			if (exists($loch{$n})) {print OUT $loch{$n}, "\t";}
			else	{print OUT "0\t";}
			}
		print OUT "\n";
		}
	}
