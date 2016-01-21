#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
$scriptname=$0; $scriptname =~ s/.+\///g;

unless ($#ARGV==2) 
	{
	print "\nExcludes samples with too much missing data (genotypes called at too few loci)\n";
	print "Usage: $scriptname input min_data print_option\n";
	print "Where:\n";
	print "\t\tinput: \t\tgenotypes matrix where rows=loci and columns=samples.\n";
	print "\t\t\t\trow 1 = column label, column 1 = tag, column 2 = position\n";
	print "\t\tmin_data: \tsamples in which fewer loci than specified here were genotyped will be excluded\n";
	print " \t\tprint_option: \ty = print genotypes and summary, n = only summary\n\n"; 
	exit;
	}

my $mindat = $ARGV[1];
my $opt = $ARGV[2];

open(IN, "$ARGV[0]"); while(<IN>) {$irowcount++;} close(IN);
$maxmiss = $irowcount - 1 - $mindat;

open(IN, "$ARGV[0]");
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	if($rowcount==1){next;}
	$ncols = @cols;
	for ($a=2; $a<$ncols; $a++)
		{
		if ($cols[$a]eq "0") {$ch{$a}++;}
		}
	}
close(IN);

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
			print $hdr[0], "\t", $hdr[1], "\t";
			for ($a=2; $a<$ncols; $a++)
				{
				if ($ch{$a}<=$maxmiss) {print $cols[$a], "\t";}
				}
			print "\n";
			} 
		next;
		}
	$ncols = @cols;
	if ($opt eq "y") {print $cols[0], "\t", $cols[1], "\t";}
	for ($a=2; $a<$ncols; $a++)
		{
		if($opt eq "y") {if ($ch{$a}<=$maxmiss) {print $cols[$a], "\t";}}
		}
	if ($opt eq "y") {print "\n";}
	}
close(IN);

foreach $k (sort(keys(%ch)))
	{
	if ($ch{$k}>$maxmiss) {$toofew++;}
	}

print STDERR $ncols-2, " samples\n";
print STDERR $toofew, " samples with more than ", $maxmiss, " missing data\n";
print STDERR $ncols-2-$toofew, " samples with no more than more than ", $maxmiss, " missing data (i.e. at least $mindat loci genotyped)\n";
