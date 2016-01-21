#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 3 || $ARGV[0] eq "-h") 
	{
	print "\nDetermines genotypes from nucleotide frequencies based on user defined frequency thresholds.\n";
	print "Usage: $scriptname input max_MAF_homo min_MAF_hetero min_cov\n"; 
	print "Where:\tinput:\t\tInput file, tab delimited text file of nucleotide frequencies (output from AlfIBasecaller.pl)\n";
	print "\tmax_MAF_homo:\tLoci at which the minor allele frequency (MAF) is below this threshold will be called\n";
	print "\t\t\thomozygous for the major allele (i.e. the minor allele is considered seuqencing error and ignored).\n";
	print "\tmin_MAF_hetero:\tLoci at which MAF exceeds this frequency will be called heterozygous for the two observed alleles.\n";
	print "\t\t\t(Note: loci with MAF between these two thresholds will be called \"unknown\" and not reported).\n";
	print "\tmin_cov:\tMinimum coverage required to determine genotypes. Lower coverage loci wil be discarded.\n\n";
	exit;
	}
use Bio::SeqIO;

# define variables from user input
my $ftabfile = $ARGV[0];
my $lethd = $ARGV[1];
my $uethd = $ARGV[2];
my $mincov = $ARGV[3];

# make a big hash of all base calls
my %bh;
@boa = qw{A C G T N};
open(FTAB, $ftabfile);
while(<FTAB>)
	{chomp;
	if($_ =~ /\tRef\t/) {next;}
	@cols = split("\t", $_);
	for ($i=3; $i<=7; $i++)	{$bh{$cols[1]}{$cols[2]}{$boa[$i-3]} = $cols[$i];}
	}
close(FTAB);		

#print "@boa\n";
#foreach $k (sort(keys(%bh)))
#	{%kh = %{$bh{$k}};
#	foreach $l (sort(keys(%kh)))
#		{print $k, "\t", $l, "\t", $kh{$l}{"A"}, "\t";
#		print $kh{$l}{"C"}, "\t", $kh{$l}{"G"}, "\t", $kh{$l}{"T"}, "\n";}}

# loop through all sites and identify limits of any terminal gaps
@sites = sort(keys(%bh));
$gtl = 0;
foreach $s (@sites)
	{
	%siteh = %{$bh{$s}};
	@loca = sort{$a<=>$b}(keys(%siteh));
	$lcount = 0; $status = 0; $bgap = 0;
	foreach $l (@loca)
		{$lcount++;
		if ($siteh{$l}{"N"}>0 && $status==0) {$bgap++; next;}
		elsif ($siteh{$l}{"N"}==0) {$status++; $egap = 0; next;}
		else {$egap++;}
		}
#	print $lcount, "\t", $bgap, "\t", $egap, "\n";
	if ($bgap>0) {$th{$s}{"bgap"} = $loca[$bgap-1];} else {$th{$s}{"bgap"} = $loca[0];}
	if ($egap>0) {$th{$s}{"egap"} = $loca[$lcount-$egap-1];} else {$th{$s}{"egap"} = $loca[-1];}
#	print $lcount, "\t", $th{$s}{"bgap"}, "\t", $th{$s}{"egap"}, "\n";
	}

# loop through all sites again and apply genotyping rules	
@cba = qw{A C G T};
foreach $s (@sites)
	{
	%siteh = %{$bh{$s}};
	@loca = sort{$a<=>$b}(keys(%siteh));

	foreach $l (@loca)
		{
		if ($l<=$th{$s}{"bgap"}) {next;}
		if ($l>=$th{$s}{"egap"}) {next;}

		%loch = %{$siteh{$l}};
		%tmph = ();
		foreach $c (@cba) {$tmph{$c} = $loch{$c};}
		@saa = sort{$tmph{$b}<=>$tmph{$a}}(keys(%tmph));
#		foreach $s (@cba) {print $s, " ", $tmph{$s}, "\t";} print "\n";

		if (($tmph{$saa[0]}+$tmph{$saa[1]})<$mincov) {next;}	# exclude low coverage
		if ($tmph{$saa[2]}>0) {next;}				# exclude loci with >2 alleles
		$rati = $tmph{$saa[1]}/($tmph{$saa[0]}+$tmph{$saa[1]});
		if ($rati<$lethd)
			{
			$fh{$s}{$l}{"type"} = "homo";	
			$fh{$s}{$l}{"call"} = $saa[0];	
			}
		elsif ($rati >= $lethd && $rati < $uethd)
			{
			$fh{$s}{$l}{"type"} = "uncertain";	
			@rsa = sort(@saa[0..1]);
			$fh{$s}{$l}{"call"} = "@rsa";	
			}
		elsif ($rati>=$uethd)
			{
			$fh{$s}{$l}{"type"} = "hetero";	
			@rsa = sort(@saa[0..1]);
			$fh{$s}{$l}{"call"} = "@rsa";	
			}
		if ($fh{$s}{$l}{"type"} eq "uncertain") {next;}

		$gtl++;	
		print $gtl, "\t", $s, "\t", $l, "\t";
		foreach $s (@cba) {print $tmph{$s}, "\t";} 
		print $fh{$s}{$l}{"type"}, "\t", $fh{$s}{$l}{"call"}, "\n", 
	
		}
	}
