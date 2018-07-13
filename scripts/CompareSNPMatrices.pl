#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Compare two matrices of SNP genotypes (e.g. produced from RAD data) from the same set of 
samples to evaluate overlap in genotyped loci, and the level of agreement and disagreement in genotypes.
This is useful for comparing different genotyping algorithms. 
Input files are formatted as the output from NFGenotyper or BGCGenotyper: tab-delimited text, 
rows are loci, columns 1-2 are tag and locus and subsequent columns are samples, homozygotes 
shown as e.g. "A", heterozygotes as e.g. "A C", and missing data as "0". 

Usage: $scriptname -f file1 -s file2 <OPTIONS>
Required arguments:
	-f file1	name of the first SNP matrix (tab delimited text)
	-s file2	name of the second SNP matrix (tab delimited text)
Options:
	-o option	0: (default) don't print any detailed info on disagreements
			1: show detailed info on loci called different homozygous genotypes in each file
			2: show detailed info on loci called different heterozygous genotypes in each file
			3: show detailed info on loci called homozygous in file1 and heterozygous in file2 
			4: show detailed info on loci called homozygous in file2 and heterozygous in file1 
			5: show detailed info on loci called in file 1 but not in file 2
	-b counts	(required if -o > 0) the file of allele counts from which genotypes were called.
USAGE
if ($#ARGV < 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('f:s:o:b:h');	# in this example a is required, b is optional, h is help
if (!$opt_f || !$opt_s || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_o && !$opt_b ) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_o) {$eopt = $opt_o;} else {$eopt = 0;}	
if ($opt_b) {$bfile = $opt_b;}

$file1 = $opt_f;
$file2 = $opt_s;

open(IN, $file1);
$rowno = 0;
while(<IN>)
	{
	chomp;
	$rowno++;
	@cols = split("\t", $_);
	if($rowno eq 1) 
		{
		@labels = @cols;
		next;
		}
	for ($a=2; $a<@cols-1; $a++)
		{
		if ($cols[$a] eq 0) {next;}
		$gt1h{$cols[0]}{$cols[1]}{$labels[$a]} = $cols[$a];
		}	
	}
close(IN);


open(IN, $file2);
$rowno = 0;
while(<IN>)
	{
	chomp;
	$rowno++;
	@cols = split("\t", $_);
	if($rowno eq 1) 
		{
		@labels = @cols;
		next;
		}
	for ($a=2; $a<@cols-1; $a++)
		{
		if ($cols[$a] eq 0) {next;}
		$gt2h{$cols[0]}{$cols[1]}{$labels[$a]} = $cols[$a];
		}	
	}
close(IN);

print "\n";
$gt1 = $also2 = $samein2 = $het1het2 = $het1hom2 = $hom1het2 = $hom1hom2 = 0;
foreach $tag (sort(keys(%gt1h)))
	{
	%tagh = %{$gt1h{$tag}};
	foreach $loc (sort{$a<=>$b}(keys(%tagh)))
		{
		%loch = %{$tagh{$loc}};
		foreach $samp (sort(keys(%loch)))
			{
			if ($loch{$samp} eq 0) {next;}
			$gt1++;
			if (exists($gt2h{$tag}{$loc}{$samp}))
				{
				$also2++;
				if ($gt1h{$tag}{$loc}{$samp} eq $gt2h{$tag}{$loc}{$samp})
					{
					$samein2++;
					}
				else
					{
					if ($gt1h{$tag}{$loc}{$samp} =~ / /)
						{
						if ($gt2h{$tag}{$loc}{$samp} =~ / /)
							{
							$het1het2++;
							if ($eopt eq 2)
								{
								print $tag, "\t", $loc, "\n";
								system("grep -P '$tag\t$loc\t' $file1");
								system("grep -P '$tag\t$loc\t' $file2");
								system("grep -P '$tag\t$loc\t' $bfile");
								}
							}
						else
							{
							$het1hom2++;
							if ($eopt eq 4)
								{
								print $tag, "\t", $loc, "\n";
								system("grep -P '$tag\t$loc\t' $file1");
								system("grep -P '$tag\t$loc\t' $file2");
								system("grep -P '$tag\t$loc\t' $bfile");
								}
							}
						}
					else
						{
						if ($gt2h{$tag}{$loc}{$samp} =~ / /)
							{
							$hom1het2++;
							if ($eopt eq 3)
								{
								print $tag, "\t", $loc, "\n";
								system("grep -P '$tag\t$loc\t' $file1");
								system("grep -P '$tag\t$loc\t' $file2");
								system("grep -P '$tag\t$loc\t' $bfile");
								}
							}
						else
							{
							$hom1hom2++;
							if ($eopt eq 1)
								{
								print $tag, "\t", $loc, "\n";
								system("grep -P '$tag\t$loc\t' $file1");
								system("grep -P '$tag\t$loc\t' $file2");
								system("grep -P '$tag\t$loc\t' $bfile");
								}
							}
						}
					}
				}
			else
				{
				if ($eopt eq 5)
					{
					print $tag, "\t", $loc, "\n";
					system("grep -P '$tag\t$loc\t' $file1");
					system("grep -P '$tag\t$loc\t' $file2");
					system("grep -P '$tag\t$loc\t' $bfile");
					}
				}
			}
		}
	}
print "\n";
print "$gt1 genotypes in $file1.\n";
$pboth = $also2 / $gt1 * 100;
print "$also2 of these also genotyped in $file2. ($pboth %)\n";
$diff = $gt1 - $also2;
print "$diff genotypes in $file1 but not in $file2\n";
$psame = $samein2 / $also2 * 100;
print "$samein2 of these were the same genotype in both files. ($psame %)\n";
$difftot = $also2 - $samein2;
print "$difftot of these were different genotypes in these files.\n";
print "$hom1hom2 were different homozygous genotypes in each file.\n";
print "$het1het2 were different heterozygous genotypes in each file.\n";
print "$hom1het2 were homozygous in $file1 and heterozygous in $file2.\n";
print "$het1hom2 were heterozygous in $file1 and homozygous in $file2.\n";
print "\n";

