#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Converts a SNP matrix (produced from CallGenotypes.pl) into the format required by
the software DADI, described at: https://bitbucket.org/gutenkunstlab/dadi/wiki/DataFormats
Usage: $scriptname -i input -k key -r reference -o output
Where:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-k key		a tab-delimited text file associating each sample in the input
			with a population label. Alleles will be counted and reported
			by the population labels assigned in this file. Formated as:
			Sample_name	Population_name
	-r reference	Name of the reference file from which these SNPs were called (FASTA format)
	-o output	a name for the output file. Tab delimited text in DADI format.

USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:k:r:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_k || !$opt_r || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# build variables from user input
$infile = $opt_i;
$keyfile = $opt_k;
$reffile = $opt_r;
$outfile = $opt_o;

# read population info from key file
open(KEY,$keyfile);
while(<KEY>)
	{
	chomp;
	@cols = split("\t", $_);
	$kh{$cols[0]} = $cols[1];
	$poph{$cols[1]}++;
	}
close(KEY);
@poplist = sort(keys(%poph));

# read sequence info from reference file
open(REF,$reffile);
while(<REF>)
	{
	chomp;
	if ($_ =~ /^>/)
		{
		if ($sstr ne "")
			{
			$sh{$sid} = $sstr;
			$sid = $sstr = "";
			}
		$sid = $_;
		$sid =~ s/>//;
		$sid =~ s/\s.*//;
		$sstr = "";
		next;
		}
	else
		{
		$sstr = $sstr.$_;
		}
	}
$sh{$sid} = $sstr;
close(REF);
@refs = sort(keys(%sh));
$nrefs = @refs;
#print "$nrefs reference sequences in $reffile.\n";
#print "e.g.\n";
#print $refs[0], "\t", $sh{$refs[0]}, "\n";
#exit;

# main loop
open(GT,$infile);
open(OUT,">$outfile");
print OUT "NAME\tOUT\tAllele1\t";
print OUT join ("\t", @poplist), "\t";
print OUT "Allele2\t";
print OUT join ("\t", @poplist), "\t";
print OUT "Sequence\tPosition\n";

$rownum = 0;
while(<GT>)
	{
	chomp;
	@cols = split("\t", $_);
	$ncols = @cols;
	$rownum++;
# set up sample key
	if ($rownum eq 1)
		{
		@headers = @cols;
		next;
		}
# scan all genotypes to find major and minor allele
	%allh = (); $maja = $mina = "";
	for ($n=2; $n<$ncols; $n++)
		{
		$gti = $cols[$n];
		if ($gti eq 0) {next;}
		if ($gti !~ / /) {$allh{$gti}+=2;}
		else
			{
			@alls = split(" ", $gti);
			$allh{$alls[0]}++;
			$allh{$alls[1]}++;
			}
		}
	@salls = sort{$allh{$b}<=>$allh{$a}}(keys(%allh));
	$maja = @salls[0];
# kick out monomorphic and multiallelic loci
	if (!exists($salls[1])) {$mono++; next;}
	if (exists($salls[2])) {$multi++; next;}
	$mina = $salls[1];
#	print $maja, "\t", $mina, "\n";

# read genotypes
	%afh = ();
	for ($n=2; $n<$ncols; $n++)
		{
		$pop = $kh{$headers[$n]};
		$gti = $cols[$n];
		if ($gti eq 0 ) {next;}
# count alleles by populations
		elsif ($gti !~ / /) {$afh{$pop}{$gti}+=2;}
		else
			{
			@alls = split(" ", $gti);
			$afh{$pop}{$alls[0]}++;
			$afh{$pop}{$alls[1]}++;
			}
		}
#	print $cols[0], "\t", $cols[1], "\t";
#	foreach $p (@poplist)
#		{
#		print $p, ": ";
#		print $maja, " ";
#		if (exists($afh{$p}{$maja})) {print $afh{$p}{$maja};}
#		else	{print "0";}
#		print "; ";
#		print $mina, " ";
#		if (exists($afh{$p}{$mina})) {print $afh{$p}{$mina};}
#		else	{print "0";}
#		print "; ";
#		}
#	print "\n";

# write output
	$refi = $cols[0];
	$posi = $cols[1];
# exclude postions 1 and 36, since flanking sequences cant be extracted for SNPs at those positions
	if ($posi eq 1 || $posi eq 36) {next;}
	$trimer = substr ($sh{$refi}, $posi-2, 3);
	
	print OUT $trimer, "\t", $trimer, "\t";
	print OUT $maja, "\t";
	foreach $p (@poplist)
		{
		if (exists($afh{$p}{$maja})) {print OUT $afh{$p}{$maja};}
		else	{print OUT "0";}
		print OUT "\t";
		}
	print OUT $mina, "\t";
	foreach $p (@poplist)
		{
		if (exists($afh{$p}{$mina})) {print OUT $afh{$p}{$mina};}
		else	{print OUT "0";}
		print OUT "\t";
		}
	@parts = split("_", $refi);
	$snploc = $posi + $parts[1] - 1;
	print OUT $parts[0], "\t";
	print OUT $snploc, "\n";
	}
