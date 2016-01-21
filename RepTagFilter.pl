#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
$scriptname=$0; $scriptname =~ s/.+\///g;

unless ($#ARGV==2) 
	{
	print "\nExcludes tags containing too many SNPs, suggesting repetive regions of the genome\n";
	print "Usage: $scriptname input max_snps print_option\n";
	print "Where:\n";
	print "\t\tinput: \t\tgenotypes matrix where rows=loci and columns=samples.\n";
	print "\t\t\t\trow 1 = column label, column 1 = tag, column 2 = position\n";
	print "\t\tmax_snps: \ttags containing more than this number of SNPs will be excluded\n";
	print " \t\tprint_option: \ty = print genotypes and summary, n = only summary\n\n"; 
	exit;
	}

open(IN, "$ARGV[0]");
my $maxsnps = $ARGV[1];
my $opt = $ARGV[2];

while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){next;}
	@cols = split("\t", $_);
	$ncols = @cols;
#	print$cols[0], "\n";
	$ch{$cols[0]}++;
	}
close(IN);

$rowcount = 0;
open(IN, "$ARGV[0]");
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	if ($ch{$cols[0]}>$maxsnps) {$toomany++; next;}
	if ($opt eq "y") {print $_, "\n";}
	}
close(IN);
print STDERR $rowcount-1, " loci\n";
print STDERR $toomany, " loci with more than ", $maxsnps, " SNPs per tag\n";
print STDERR $rowcount-1-$toomany, " loci with no more than more than ", $maxsnps, " SNPS per tag\n";
