#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h") 
	{
	print "\nConverts a genotype matrix (loci x samples) to a snp matrix, as described in the manual\n";
	print "for the R package diveRsity. This snp matrix can be converted to genepop format\n";
	print "using diveRsity's snp2gen function.\n";
	print "Usage: $scriptname input > output\n"; 
	print "Where:\tinput:\ttab-delimited genotype matrix, with rows=loci and columns=samples.\n";
	print "\t\tFirst two columns indicate tag and position respectively.\n";
	print "\t\tThis format is the output from RADGenotyper.pl.\n";
	print "\toutput:\ta name for the output file. SNP matrix format.\n\n";
	exit;
	}

open(IN, $ARGV[0]);
open(KEY, ">$ARGV[0].snpkey");
while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	$nom = @cols;

	if ($rowcount==1)
		{
		print "SNP_ID", "\t";
		print join("\t", @cols[2..$nom]), "\n";
		next;
		}
	$tagcount++;
	print "SNP", $tagcount, "\t";
	print KEY $tagcount, "\t", $cols[0], "\t", $cols[1], "\n";
	for ($a=2; $a<$nom; $a++)
		{
		if ($cols[$a] =~ / /)
			{
			$cols[$a] =~ s/ //;
			print $cols[$a], "\t";
			}
		elsif ($cols[$a] eq 0)
			{
			print "--", "\t";
			}
		else
			{
			print $cols[$a], $cols[$a], "\t";	
			}
		}
	print "\n";
	}
