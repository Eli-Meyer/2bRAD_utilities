#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Counts observations of alleles at each locus in a collection of base counts from 2bRAD 
(the output from SAMBaseCounts.pl). 

This script identifies the major and minor allele at each locus, and combines all samples
into a single file containing the number of times each of these alleles was observed
in each sample (two columns per sample, for major and minor allele respectively).

Output format: columns 1=tag, 2=position, 3=major allele, 4=minor allele,
5=major allele counts for sample A, 6=minor allele counts for sample A, etc.

Missing data are shown as "NA" for both alleles, and the minor allele is reported as 
"NA" for monomorphic loci.

Usage: $scriptname <options> file_1 file_2 ... file_n > output_file
Required arguments:
	files 1-n:	nucleotide frequencies (output from SAMBaseCounts.pl) for each sample
	output_file:	a name for the output; tab-delimited text
Options:
	-a max_alleles	maximum number of alleles allowed at each locus. Loci with more than this
			number of alleles will be excluded. (default=2)
        -v min_cov      minimum coverage required to consider an allele present. (default=2)
        -s min_samp     minimum number of samples in which an allele must be present
                        (default=1)

USAGE
if ($#ARGV < 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('a:v::s:h');	# in this example a is required, b is optional, h is help
if ($opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_a) {$avar = $opt_a;} else {$avar = 2;}	# max number of alleles allowed
if ($opt_v) {$vvar = $opt_v;} else {$vvar = 2;}	# max number of alleles allowed
if ($opt_s) {$svar = $opt_s;} else {$svar = 1;}	# max number of alleles allowed

my %bigh;

print "Tag\tLocus\tMajor\tMinor\t";
foreach $argi (0..$#ARGV)
	{
	$tmpname = $ARGV[$argi];
	$tmpname =~ s/\/[a-zA-Z0-9_]+\.tab//;
	$tmpname =~ s/.+\///g;
	print $tmpname.".maj", "\t", $tmpname.".min", "\t";
	print STDERR "reading ", $ARGV[$argi], "\n";
	open(TAB, $ARGV[$argi]);
	$rownum = 0;
	while(<TAB>)
		{
		chomp;
		$rownum++; if ($rownum==1) {next;}
		@cols = split("\t", $_);
		$bigh{$cols[0]}{$cols[1]}{$argi}{"A"} = $cols[3];		
		$bigh{$cols[0]}{$cols[1]}{$argi}{"C"} = $cols[4];		
		$bigh{$cols[0]}{$cols[1]}{$argi}{"G"} = $cols[5];		
		$bigh{$cols[0]}{$cols[1]}{$argi}{"T"} = $cols[6];		
		}
	}
print "\n";

foreach $r (sort(keys(%bigh)))
	{
	%ch = %{$bigh{$r}};
	foreach $l (sort{$a<=>$b}(keys(%ch)))
		{
		my %th;
		foreach $argi (0..$#ARGV)
			{
			if ($bigh{$r}{$l}{$argi}{"A"} >= $vvar) {$th{"A"} += $bigh{$r}{$l}{$argi}{"A"};}
			if ($bigh{$r}{$l}{$argi}{"C"} >= $vvar) {$th{"C"} += $bigh{$r}{$l}{$argi}{"C"};}
			if ($bigh{$r}{$l}{$argi}{"G"} >= $vvar) {$th{"G"} += $bigh{$r}{$l}{$argi}{"G"};}
			if ($bigh{$r}{$l}{$argi}{"T"} >= $vvar) {$th{"T"} += $bigh{$r}{$l}{$argi}{"T"};}
			}
		@stha = sort{$th{$b}<=>$th{$a}}(keys(%th));
		$passno = 0;
		foreach $s (@stha)
			{
			if ($th{$stha[$s]} >= $svar) {$passno++;}
			}
		if ($passno >= $avar) {$rejmulti++; next;}
		print $r, "\t", $l, "\t", $stha[0], "\t";
		if ($th{$stha[1]}<1) {print "NA", "\t";}
		else	{print $stha[1], "\t";}
		foreach $argi (0..$#ARGV)
			{
			if (exists($bigh{$r}{$l}{$argi}{$stha[0]}) || exists($bigh{$r}{$l}{$argi}{$stha[1]})) 
				{
				print $bigh{$r}{$l}{$argi}{$stha[0]}, "\t", $bigh{$r}{$l}{$argi}{$stha[1]}, "\t";
				}
			else {print "NA\tNA\t";}
			}
		print "\n";
		}
	}
exit;
