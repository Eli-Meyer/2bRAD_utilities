#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
$scriptname=$0; $scriptname =~ s/.+\///g;

unless ($#ARGV==2) 
	{
	print "\nExcludes loci containing too few genotypes (keeps polymorphic loci)\n";
	print "Note that this counts genotypes, not alleles (e.g. AA, AB, BB = 3 genotypes)\n";
	print "Usage: $scriptname input min_gt print_option\n";
	print "Where:\n";
	print "\t\tinput: \t\tgenotypes matrix where rows=loci and columns=samples.\n";
	print "\t\t\t\trow 1 = column label, column 1 = tag, column 2 = position\n";
        print "\t\tmin_gt: \tminimum number of genotypes required to consider a locus polymorphic\n";
	print " \t\tprint_option: \ty = print genotypes and summary, n = only summary\n\n"; 
	exit;
	}

open(IN, "$ARGV[0]");
my $minall = $ARGV[1];
my $opt = $ARGV[2];

while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print $_, "\n";} next;}
	@cols = split("\t", $_);
	$ncols = @cols;
	%ahi = ();
	for ($a=2; $a<$ncols; $a++)
		{
		if ($cols[$a] eq "0") {next;}
		$ahi{$cols[$a]}++;
		
		}
	@alla = sort(keys(%ahi));
	$nalla = @alla;
	if ($nalla < $minall) {$nopoly++; next;}
	if ($opt eq "y") {print $_, "\n";}
	}
close(IN);
print STDERR $rowcount-1, " loci in input file\n";
print STDERR $nopoly, " loci with fewer than ", $minall, " genotypes\n";
print STDERR $rowcount-1-$nopoly, " polymorphic loci (at least ", $minall, " genotypes)\n";
