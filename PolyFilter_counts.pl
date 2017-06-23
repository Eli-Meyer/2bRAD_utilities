#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Excludes monomorphic loci and loci where an alternate allele was observed
in too few samples to be useful. 
Note that this script counts alleles, not genotypes.
Usage: $scriptname input min_alt min_freq print_option
Where:
                input:          genotypes matrix where rows=loci and columns=samples.
                                row 1 = column label, column 1 = tag, column 2 = position, 
				column 3 = reference allele, column 4 = alternate allele
                min_alt:	alternate allele must be detected at least this many times in a sample to count as present
		min_freq:	alternate allele must be present in at least this number of samples to pass the filter
                print_option:   y = print allele counts and summary, n = only summary
USAGE
if ($#ARGV != 3 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

open(IN, "$ARGV[0]");
my $minall = $ARGV[1];
my $minfreq = $ARGV[2];
my $opt = $ARGV[3];

while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	$gt_no = $alt_no = $alt_freq = 0;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if ($cols[$a] eq "NA") {next;}
		$gt_no++;
		if ($cols[$a+1] >= $minall) 
			{
			$alt_no++;
			}
		}
	$alt_freq = $alt_no / $gt_no;
	if ($alt_freq eq 0) {$mono++; next;}
	if ($alt_freq < $minfreq) {$toolow++; next;}
	if ($opt eq "y") {print $_, "\n";}
	}
close(IN);
print STDERR $rowcount-1, " loci in input file\n";
print STDERR $toolow, " loci with minor alleles present in fewer than $minfreq of sample\n";
print STDERR $mono, " monomorphic loci\n";
print STDERR $rowcount-1-$toolow-$mono, " polymorphic loci\n";
