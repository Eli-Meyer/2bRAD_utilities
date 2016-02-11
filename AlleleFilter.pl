#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Excludes loci containing too many alleles
Usage: $scriptname input max_alleles print_option
Where:
                input:          genotypes matrix where rows=loci and columns=samples.
                                row 1 = column label, column 1 = tag, column 2 = position
                max_alleles:    loci with more than this number of alleles are excluded
                print_option:   y = print genotypes and summary, n = only summary

USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

open(IN, "$ARGV[0]");
my $maxall = $ARGV[1];
my $opt = $ARGV[2];

$toomany = 0;
while(<IN>)
	{
	chomp;
	$rowcount++;
	if($rowcount==1){if($opt eq "y") {print $_, "\n";} next;}
	@cols = split("\t", $_);
	%allh = ();
	foreach $c (@cols[2..@cols])
		{
		if ($c =~ /0/) {next;}
		@alls = split(" ", $c);
		foreach $a (@alls) {$allh{$a}++;}
		}
	@aa = sort(keys(%allh));
	if (@aa > $maxall) 
		{
		$toomany++; 
#		print "@aa", "\n", $_, "\n"; 
		next;
		}
	if ($opt eq "y") {print $_, "\n";}
	}
close(IN);
print STDERR $rowcount-1, " loci\n";
print STDERR $toomany, " loci with more than ", $maxall, " alleles\n";
print STDERR $rowcount-1-$toomany, " loci with no more than ", $maxall, " alleles\n";
