#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h") 
	{
	print "\nConverts a genotype matrix (loci x samples) to a PHYLIP-formatted alignment.\n";
	print "Usage: $scriptname input > output\n"; 
	print "Where:\tinput:\ttab-delimited genotype matrix, with rows=loci and columns=samples.\n";
	print "\t\tFirst two columns indicate tag and position respectively.\n";
	print "\t\tThis format is the output from RADGenotyper.pl.\n";
	print "\toutput:\ta name for the output file. PHYLIP alignment format.\n\n";
	exit;
	}

$ch{"A C"} = "M"; $ch{"C A"} = "M";
$ch{"A G"} = "R"; $ch{"G A"} = "R";
$ch{"A T"} = "W"; $ch{"T A"} = "W";
$ch{"C G"} = "S"; $ch{"C G"} = "S";
$ch{"C T"} = "Y"; $ch{"C T"} = "Y";
$ch{"G T"} = "K"; $ch{"G T"} = "K";

$window = 10;

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


print " " x 4, $nom, " " x 3, $rowcount-1, "\n";
$bsize = 50;
$noblocks = int(($rowcount-1)/$bsize+0.5);
foreach $a (1..$noblocks)
	{

foreach $s (sort(keys(%gh)))
	{
	%sh = %{$gh{$s}};
	if ($a==1) 
		{print $s;
		$lens = length($s);
		$diff = $window-$lens;
		print " " x $diff;
		}
	else	{print " " x $window;}
	print " ";
	for ($c=0; $c<5; $c++)
		{
	for ($b=0; $b<10; $b++)
		{
		$cni = $b + ($c * 10) + ($bsize * ($a-1));
#		print $cni, "\n";
		print $sh{$cni+1};
		}
		print " ";
		}
	print "\n";
	}
	print "\n";
	}
