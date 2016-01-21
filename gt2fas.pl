#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h") 
	{
	print "\nConverts a genotype matrix (loci x samples) to a FASTA-formatted alignment.\n";
	print "Usage: $scriptname input > output\n"; 
	print "Where:\tinput:\ttab-delimited genotype matrix, with rows=loci and columns=samples.\n";
	print "\t\tFirst two columns indicate tag and position respectively.\n";
	print "\t\tThis format is the output from RADGenotyper.pl.\n";
	print "\toutput:\ta name for the output file. FASTA alignment format.\n\n";
	exit;
	}

$ch{"A C"} = "M"; $ch{"C A"} = "M";
$ch{"A G"} = "R"; $ch{"G A"} = "R";
$ch{"A T"} = "W"; $ch{"T A"} = "W";
$ch{"C G"} = "S"; $ch{"C G"} = "S";
$ch{"C T"} = "Y"; $ch{"C T"} = "Y";
$ch{"G T"} = "K"; $ch{"G T"} = "K";

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
		if ($cols[$a] =~ / /)
			{
			$cols[$a] = $ch{$cols[$a]};
			}
		if ($cols[$a] eq 0)
			{
			$cols[$a] = "-";
			}
		$gh{$noma[$a-2]}{$rowcount-1} = $cols[$a];
		}
	}


foreach $s (sort(keys(%gh)))
	{
	%sh = %{$gh{$s}};
	print ">",$s,"\n";
	for ($b=1; $b<$rowcount; $b++)
		{
		print $sh{$b};
		}
	print "\n";
	}
