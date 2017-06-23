#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Excludes samples with too much missing data (genotypes called at too few loci)
Usage: $scriptname input min_data print_option
Where:
                input:          allele counts matrix where rows=loci and columns=samples.
                                row 1 = column label, column 1 = tag, column 2 = position
                                column 3 = reference allele, column 4 = alternate allele
				subsequent pairs of columns contain allele counts (reference then alternate) for each sample
                min_data:       samples in which fewer loci than specified here were called will be excluded
                print_option:   y = print allele counts and summary, n = only summary
USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $mindat = $ARGV[1];
my $opt = $ARGV[2];

open(IN, "$ARGV[0]"); while(<IN>) {$irowcount++;} close(IN);
$maxmiss = $irowcount - 1 - $mindat;
#print $maxmiss, "\n";

open(IN, "$ARGV[0]");
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	if($rowcount==1){next;}
	$ncols = @cols;
	$nsamps = ($ncols-4)/2;
	for ($a=4; $a<$ncols; $a+=2)
		{
		if (($cols[$a] eq "0" && $cols[$a+1] eq "0")||($cols[$a] eq "NA" && $cols[$a+1] eq "NA")) 
			{
			$ch{$a}++;
#			print $a, "\t", $_, "\n";
			}
		}
	}
close(IN);

#foreach $k (sort(keys%ch))	{print $k, "\t", $ch{$k}, "\n";}
#exit;

$toofew = 0;
open(IN, "$ARGV[0]");
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	if($rowcount==1)
		{
		if($opt eq "y") 
			{
			@hdr = @cols;
			print $hdr[0], "\t", $hdr[1], "\t", $hdr[2], "\t", $hdr[3], "\t";
			for ($a=4; $a<$ncols; $a+=2)
				{
				if ($ch{$a}<=$maxmiss) {print $cols[$a], "\t", $cols[$a+1], "\t";}
				}
			print "\n";
			} 
		next;
		}
	$ncols = @cols;
	if ($opt eq "y") {print $cols[0], "\t", $cols[1], "\t", $cols[2], "\t", $cols[3], "\t";}
	for ($a=4; $a<$ncols; $a+=2)
		{
		if($opt eq "y") {if ($ch{$a}<=$maxmiss) {print $cols[$a], "\t", $cols[$a+1], "\t";}}
		}
	if ($opt eq "y") {print "\n";}
	}
close(IN);

foreach $k (sort(keys(%ch)))
	{
	if ($ch{$k}>$maxmiss) {$toofew++;}
	}

print STDERR $nsamps, " samples\n";
print STDERR $toofew, " samples with more than ", $maxmiss, " missing data\n";
print STDERR $nsamps-$toofew, " samples with no more than more than ", $maxmiss, " missing data (i.e. at least $mindat loci called)\n";

