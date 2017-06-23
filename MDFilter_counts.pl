#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Excludes loci containing too many missing data (genotyped in too few samples)
Usage: $scriptname input min_data print_option
Where:
                input:          allele counts matrix where rows=loci and columns=samples.
                                row 1 = column label, column 1 = tag, column 2 = position
                                column 3 = reference allele, column 4 = alternate allele
                                subsequent pairs of columns contain allele counts (major then minor) for each sample
                min_data:       loci that were genotyped in fewer samples than this will be excluded
                print_option:   y = print genotypes and summary, n = only summary
USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

open(IN, "$ARGV[0]");
my $mindat = $ARGV[1];
my $opt = $ARGV[2];

$toofew = 0;
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	if ($rowcount==2)
		{
		$maxmiss = ($ncols - 4)/2 - $mindat;
		}
	$missi = 0;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if (($cols[$a] eq 0 && $cols[$a+1] eq 0)||($cols[$a] eq "NA" && $cols[$a+1] eq "NA")) {$missi++;}
		}
	if ($missi>$maxmiss) {$toofew++; next;}
	if ($opt eq "y") {print $_, "\n";}
	}
close(IN);
print STDERR $rowcount-1, " loci\n";
print STDERR $toofew, " loci with more than ", $maxmiss, " missing data (i.e. genotyped in fewer than $mindat samples)\n";
print STDERR $rowcount-1-$toofew, " loci with no more than more than ", $maxmiss, " missing data\n";
