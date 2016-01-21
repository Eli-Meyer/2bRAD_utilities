#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
$scriptname=$0; $scriptname =~ s/.+\///g;

unless ($#ARGV==1) 
	{
	print "\nSelects a single SNP from each tag in a 2bRAD genotype matrix to minimize missing data.\n";
	print "Usage: $scriptname input print_option\n";
	print "Where:\n";
	print "\t\tinput: \t\tgenotypes matrix where rows=loci and columns=samples.\n";
	print "\t\t\t\trow 1 = column label, column 1 = tag, column 2 = position\n";
	print " \t\tprint_option: \ty = print genotypes and summary, n = only summary\n\n"; 
	exit;
	}

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
	foreach $a (@cols[2..$ncols]) {if ($a eq "0") {$mdcount++;}}
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
