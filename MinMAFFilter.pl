#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Excludes loci that are monomorphic in all samples, or loci where the alternate
allele was never present in any sample above a user specified minimum allele frequency.
Usage: $scriptname input min_maf print_option
Where:
                input:          matrix of allele frequencies from CombineMAFs.pl, where
                                row 1 = column label, column 1 = tag, column 2 = position, 
				column 3 = major allele, column 4 = minor allele, 
				columns 5 - n: frequency of the alternate allele in each sample
                min_maf:    	loci where minor allele was never present at a greater allele frequency than this will be excluded
                print_option:   y = print data and summary, n = only summary

USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

open(IN, "$ARGV[0]");
my $minmaf = $ARGV[1];
my $opt = $ARGV[2];

$rowcount = $mono = $lowfreq = 0;
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print $_, "\n";} next;}
	@cols = split("\t", $_);
	%allh = ();
	if ($cols[3] eq "NA") {$mono++; next;}
	foreach $c (@cols[4..@cols])
		{
		if ($c =~ /"NA"/) {next;}
		$allh{$c}++;
		}
	@aa = sort{$b<=>$a}(keys(%allh));
	if ($aa[0] < $minmaf) 
		{
		$lowfreq++;
		next;
		}
	if ($opt eq "y") {print $_, "\n";}
	}
close(IN);
print STDERR $rowcount-1, " loci\n";
print STDERR $mono, " loci were monomorphic across all samples\n";
print STDERR $lowfreq, " loci had a second allele but its frequency never exceeded $minmaf\n";
print STDERR $rowcount-1-$mono-$lowfreq, " loci with a second allele present at acceptable frequencies (>= $minmaf) in one or more samples\n";
