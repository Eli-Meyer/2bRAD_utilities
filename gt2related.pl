#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h") 
	{
        print "\nConverts a genotype matrix (loci x samples) into the appropriate input format for the R package related.\n";
	print "Usage: $scriptname input > output\n"; 
	print "Where:\tinput:\ttab-delimited genotype matrix, with rows=loci and columns=samples.\n";
	print "\t\tFirst two columns indicate tag and position respectively.\n";
	print "\t\tThis format is the output from RADGenotyper.pl.\n";
	print "\toutput:\ta name for the output file. RELATED input format.\n\n";
	exit;
	}

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
	$locname = $cols[0]."-".$cols[1];
	$loch{$rowcount-1} = $locname;
	for ($a=2; $a<$nom; $a++)
		{
		$gh{$noma[$a-2]}{$rowcount-1} = $cols[$a];
		}
	}

$lkph{"A"}=1;
$lkph{"C"}=2;
$lkph{"G"}=3;
$lkph{"T"}=4;

#print "SampleID\t";
#foreach $s (sort{$a<=>$b}(keys(%loch))) {print $loch{$s}, "-1\t", $loch{$s}, "-2\t";}
#print "\n";

foreach $s (sort(keys(%gh)))
	{
	%sh = %{$gh{$s}};
	print $s,"\t";
	for ($b=1; $b<$rowcount; $b++)
		{
		$gi = $sh{$b};
		if ($gi eq 0)				# missing data
			{
			print "NA\tNA\t";
			next;
			}
		elsif ($gi !~ /\s/)			#homozygous
			{
			print $lkph{$gi}, "\t", $lkph{$gi}, "\t";
			next;
			}				
		elsif ($gi =~ /\s/)			#heterozygous
			{
			@alla = split(" ", $gi);
			print $lkph{$alla[0]}, "\t", $lkph{$alla[1]}, "\t";
			next;
			}
		}
	print "\n";
	}
