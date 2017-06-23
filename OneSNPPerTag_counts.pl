#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Selects a single SNP from each tag in a 2bRAD genotype matrix to minimize missing data.
Usage: OneSNPPerTag_counts.pl input print_option
Where:
                input:          allele counts matrix where rows=loci and columns=samples.
                                row 1 = column label, column 1 = tag, column 2 = position
                                column 3 = reference allele, column 4 = alternate allele
                                subsequent pairs of columns contain allele counts (major then minor) for each sample
                print_option:   y = print allele counts and summary, n = only summary
USAGE
if ($#ARGV != 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

open(IN, "$ARGV[0]");
my $opt = $ARGV[1];

while(<IN>)
	{
	chomp;
	$rowcount++;
#	print $_, "\n";
	if($rowcount==1){next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	$mdcount = 0;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if (($cols[$a] eq "0" && $cols[$a+1] eq "0")||($cols[$a] eq "NA" && $cols[$a+1] eq "NA"))
			{$mdcount++;}
		}
#	print $cols[0], "\t", $cols[1], "\t", $mdcount, "\n";
	$ch{$cols[0]}{$cols[1]} = $mdcount;
	}
close(IN);
#exit;

foreach $tag (sort(keys(%ch)))
	{
	%tagh = %{$ch{$tag}};
	@sorted = sort{$tagh{$a}<=>$tagh{$b}}(keys(%tagh));
	$useh{$tag} = $sorted[0];
#	foreach $s (@sorted)
#		{
#		print $tag, "\t", $s, "\t", $tagh{$s}, "\n";
#		}
	}

$rowcount = 0;
open(IN, "$ARGV[0]");
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print $_, "\n";} next;}
	@cols = split("\t", $_);
	if($useh{$cols[0]} ne $cols[1]) {next;}
	if ($opt eq "y") {print $_, "\n";}
	$kept++;
	}
close(IN);
print STDERR $rowcount-1, " loci\n";
print STDERR $kept, " SNPs selected.\n";
