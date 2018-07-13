#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Counts the number of times each allele was observed, for each locus, in a collection of 
2bRAD data describing nucleotide frequencies for each locus and sample (the output from
SAMBaseCounts.pl). 

Output format: columns 1=tag, 2=position, 3=reference allele,
5=allele counts for sample 1 (A/C/G/T), 6=for sample 2, etc..
Missing data are shown as "NA" for all alleles.

Usage: $scriptname file_1 file_2 ... file_n > output_file
Where:
	files 1-n:	nucleotide frequencies (output from SAMBasecaller.pl) for each sample
	output_file:	a name for the output; tab-delimited text
USAGE
if ($#ARGV < 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my %bigh;
my $maxthird = 2;

print "Tag\tLocus\tRef\t";
$countrow = 0;
foreach $argi (0..$#ARGV)
	{
	$countrow++;
	if ($countround==1) {next;}
	$tmpname = $ARGV[$argi];
	$tmpname =~ s/\/*[a-zA-Z0-9_]+\.tab//;
	$tmpname =~ s/.+\///g;
	print $tmpname, "\t";
	print STDERR "reading ", $ARGV[$argi], "\n";
	open(TAB, $ARGV[$argi]);
	$rownum = 0;
	while(<TAB>)
		{
		chomp;
		$rownum++; if ($rownum==1) {next;}
		@cols = split("\t", $_);
		$bigh{$cols[0]}{$cols[1]}{"ref"} = $cols[2];		
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
		print $r, "\t", $l, "\t", $bigh{$r}{$l}{"ref"}, "\t";
		foreach $argi (0..$#ARGV)
			{
			if (exists($bigh{$r}{$l}{$argi}))
				{
				foreach $base (qw{A C G T})
					{
					print $bigh{$r}{$l}{$argi}{$base};
					unless ($base eq "T") {print "/";}
					}
				}
			else {print "0/0/0/0";}
			print "\t";
			}
		print "\n";
		}
	}
exit;
