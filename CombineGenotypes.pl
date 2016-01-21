#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
use Bio::SeqIO;

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV<1) 
	{
	print "\nCombines output files produced by RAD_Genotyper.pl\n";
	print "Usage:\t$scriptname multiple_input_files > output_file\n"; 
        print "Arguments:\n";
        print "\t multiple_input_files\t output files from RAD_genotyper.pl (wildcards can be used)\n";
        print "\t output_file\t\t a name for the output file\n";
        print "\n"; exit;
	}

my %bigh;

print "Chr\tLocus\t";
foreach $argi (0..$#ARGV)
	{
	print $ARGV[$argi], "\t";
	open(TAB, $ARGV[$argi]);
	while(<TAB>)
		{
		chomp;
		@cols = split("\t", $_);
		$bigh{$cols[1]}{$cols[2]}{$argi} = $cols[8];		
		}
	}
print "\n";

foreach $r (sort(keys(%bigh)))
	{
	%ch = %{$bigh{$r}};
	foreach $l (sort{$a<=>$b}(keys(%ch)))
		{
		print $r, "\t", $l, "\t";
		foreach $argi (0..$#ARGV)
			{
			if(exists($bigh{$r}{$l}{$argi}))
				{print $bigh{$r}{$l}{$argi}, "\t";
				}
			else	{print 0, "\t";}
			}
		print "\n";
		}
	}
