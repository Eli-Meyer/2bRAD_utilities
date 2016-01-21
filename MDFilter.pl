#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
$scriptname=$0; $scriptname =~ s/.+\///g;

unless ($#ARGV==2) 
	{
	print "\nExcludes loci containing too many missing data (genotyped in too few samples)\n";
	print "Usage: $scriptname input min_data print_option\n";
	print "Where:\n";
	print "\t\tinput: \t\tgenotypes matrix where rows=loci and columns=samples.\n";
	print "\t\t\t\trow 1 = column label, column 1 = tag, column 2 = position\n";
	print "\t\tmin_data: \tloci that were genotyped in fewer samples than this will be excluded\n";
	print " \t\tprint_option: \ty = print genotypes and summary, n = only summary\n\n"; 
	exit;
	}
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
		$maxmiss = $ncols - 2 - $mindat;
		}
	$missi = 0;
	for ($a=2; $a<$ncols; $a++)
		{
		if ($cols[$a]eq "0") {$missi++;}
		}
	if ($missi>$maxmiss) {$toofew++; next;}
	if ($opt eq "y") {print $_, "\n";}
	}
close(IN);
print STDERR $rowcount-1, " loci\n";
print STDERR $toofew, " loci with more than ", $maxmiss, " missing data (i.e. genotyped in fewer than $mindat samples)\n";
print STDERR $rowcount-1-$toofew, " loci with no more than more than ", $maxmiss, " missing data\n";
