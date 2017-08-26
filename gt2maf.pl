#!/usr/bin/perl

$infile = $ARGV[0];

$rowcount = 0;
open(IN,$infile);
while(<IN>)
	{
	chomp;
	@cols = split("\t", $_);
	$ncols = @cols;
	$rowcount++;
	if($rowcount==1){next;}
	%ah = ();
	for ($a=2; $a<$ncols; $a++)
		{
		if ($cols[$a] =~ /0/) {next;}
		if ($cols[$a]!~ / /) {$ah{$cols[$a]}++; $ah{$cols[$a]}++;}
		else 
			{
			@alla = split(" ", $cols[$a]);
			foreach $a (@alla) {$ah{$a}++;}
			}
		}
	@alleles = sort{$ah{$b}<=>$ah{$a}}(keys(%ah));
	foreach $a (@alleles)
		{
#		print $cols[0], "\t", $cols[1], "\t", $a, "\t", $ah{$a}, "\n";
		}
	$nall = @alleles;
	print $cols[0], "\t", $cols[1], "\t";
	if ($nall<2) {print 0, "\n";}
	else
		{
		$mafi = $ah{$alleles[1]}/($ah{$alleles[1]}+$ah{$alleles[0]});
		print $mafi, "\n";
		}
	}
