#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h") 
	{
	print "\nConverts a genotype matrix (loci x samples) to FSTAT format.\n";
	print "Usage: $scriptname input > output\n"; 
	print "Where:\tinput:\ttab-delimited genotype matrix, with rows=loci and columns=samples.\n";
	print "\t\tFirst two columns indicate tag and position respectively.\n";
	print "\t\tThis format is the output from RADGenotyper.pl.\n";
	print "\toutput:\ta name for the output file. FSTAT format.\n\n";
	exit;
	}

$zh{"A"} = "01";
$zh{"C"} = "02";
$zh{"G"} = "03";
$zh{"T"} = "04";
$zh{"0"} = "0";

open(IN, $ARGV[0]);
open(OK, ">$ARGV[0].locuskey.tab");
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
	$snpno = $rowcount-1;
	print OK "SNP-", $snpno, "\t", $cols[0], "\t", $cols[1], "\n";
	$snpname = "SNP-".$snpno;
	$kh{$rowcount-1} = $snpname;
	for ($a=2; $a<$nom; $a++)
		{
		$gh{$noma[$a-2]}{$rowcount-1} = $cols[$a];
		}
	}


open(OF, ">$ARGV[0].samplekey.tab");
my $rowno = 0;
print $nom-2, " ", $rowcount-1, " ", 4, " ", 2, "\n";
foreach $k (sort(keys(%kh)))	{print $kh{$k}, "\n";}
foreach $s (sort(keys(%gh)))
	{
	$rowno++;
	%sh = %{$gh{$s}};
	print OF $rowno, "\t", $s, "\n";
	print $rowno," ";
	for ($b=1; $b<$rowcount; $b++)
		{
		$gi = $sh{$b};
		if ($gi !~ /\s/)
			{
			print $zh{$gi}, $zh{$gi}, " ";
			}
		else
			{
			@alla = split(" ", $gi);
			print $zh{$alla[0]}, $zh{$alla[1]}, " ";
			}
		}
	print "\n";
	}
close(OF);
