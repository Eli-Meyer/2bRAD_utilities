#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h") 
	{
	print "\nConverts a genotype matrix (loci x samples) into the appropriate input format for STRUCTURE.\n";
	print "Usage: $scriptname input > output\n"; 
	print "Where:\tinput:\ttab-delimited genotype matrix, with rows=loci and columns=samples.\n";
	print "\t\tFirst two columns indicate tag and position respectively.\n";
	print "\t\tThis format is the output from RADGenotyper.pl.\n";
	print "\toutput:\ta name for the output file. STRUCTURE input format.\n\n";
	exit;
	}

$zh{"A"} = 1;
$zh{"C"} = 2;
$zh{"G"} = 3;
$zh{"T"} = 4;
$zh{"0"} = 0;

open(IN, $ARGV[0]);
while(<IN>)
	{
	chomp;
	$rowcount++;
	if ($rowcount==1)
		{
		@noma = split("\t", $_);
		$nom = @noma;
		@noma = @noma[2..$nom];
		next;
		}
	@cols = split("\t", $_);
	$nom = @cols;
	for ($a=2; $a<$nom; $a++)
		{
		$gh{$noma[$a-2]}{$rowcount-1} = $cols[$a];
		}
	}


foreach $s (sort(keys(%gh)))
	{
	%sh = %{$gh{$s}};
	print $s,"-1\t";
	for ($b=1; $b<$rowcount; $b++)
		{
		$gi = $sh{$b};
		if ($gi =~ /\s/) {$gi =~ s/ .+//;}
		print $zh{$gi}, "\t";
		}
	print "\n";
	print $s,"-2\t";
	for ($b=1; $b<$rowcount; $b++)
		{
		$gi = $sh{$b};
		if ($gi =~ /\s/) {$gi =~ s/.+ //;}
		print $zh{$gi}, "\t";
		}
	print "\n";
	}
