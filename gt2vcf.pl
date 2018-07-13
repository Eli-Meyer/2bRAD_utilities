#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Converts a genotype matrix (loci x samples) to a VCF file.
Usage: $scriptname -i input -r reference -o output <options>
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. FASTA alignment format.
	-r reference	Complete path to the reference file used to generate these genotypes (FASTA). 
Options:
	-f filters	a text file described filters applied to the genotypes. this information
			will be included in the VCF file. e.g.	
			  "MD	removed loci genotyped in <20 samples"

USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:r:f:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || !$opt_r || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# open input file
open(IN, $opt_i);
open(OUT, ">$opt_o");


# add metadata lines to output file
print "##fileformat=VCFv4.2\n";
@ta = localtime(time());
@months = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12");
@days = ("01", "02", "03", "04", "05", "06", "07", "08", "09", 10..30);
print "##fileDate=", 1900+$ta[5], $months[$ta[4]], $days[$ta[3]-1], "\n";
print "##source=",$scriptname, "\n";
if (defined($rfile)) {print "##reference=", $rfile, "\n";}

# INFO lines
print "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";

# FILTER lines
if (defined($ffile))
	{
	open (FF, $ffile);
	while (<FF>)
		{
		chomp;
		@bits = split("\t", $_);
		print "##FILTER=<ID=", $bits[0], ",Description=\"", $bits[1], "\">\n";
		}
	}

# FORMAT lines
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

while(<IN>)
	{
	chomp;
	$rowcount++;
	@cols = split("\t", $_);
	$nom = @cols;

	if ($rowcount==1)
		{
		print join("\t", @cols[2..$nom]), "\n";
		next;
		}

	$tagcount++;
	print $cols[0], "\t", $cols[1], "\t", ".", "\t";

# identify and print reference allele
	$refstr = `grep $cols[0] -A 1 -m 2 $rfile | grep -v ">" | grep -v "-"`;
	@refa = split(//,$refstr);
	$refcall = @refa[$cols[1]-1];
#	print $refcall, "\t";
	
# identify and print a list of all alternate alleles
	%aah = (); $gtni = 0;
	for ($z=2; $z<$nom; $z++)
		{
		unless ($cols[$z] =~ /0/) {$gtni++;}
		@ai = split(" ",$cols[$z]);
		foreach $i (@ai) {unless ($i eq 0 || $i eq $refa[$cols[1]-1]) {$aah{$i}++;}}
		}
	@aaa = sort(keys(%aah));
	$naa = @aaa;
	if ($naa eq 0) {push @aaa, ".";}
	if ($refcall eq "N")
		{
		$refcall = $aaa[0];
		print $refcall, "\t";
		shift @aaa;
		print join(",", @aaa), "\t";
		}
	else
		{
		print $refcall, "\t";
		print join(",", @aaa), "\t";
		}
# make a hash of alleles
	%ach = (); $allno = 0;
	$ach{$refcall} = 0;
	foreach $a (@aaa)
		{
		$allno++;
		$ach{$a} = $allno;
		}

# print empty genotype quality score
	print ".", "\t";

# print filter status (all passed; assumption of this script)
	print "PASS\t";

# print number of samples with genotypes called at this locus
	print "NS=", $gtni, "\t";

# print FORMAT 
	print "GT", "\t";

# print genotype data for each sample

	for ($a=2; $a<$nom; $a++)
		{

		if ($cols[$a] =~ / /)
			{
			@obsall = split(" ", $cols[$a]);
			print $ach{$obsall[0]}, "/", $ach{$obsall[1]}, "\t";
			}
		elsif ($cols[$a] eq 0)
			{
			print "./.", "\t";
			}
		else
			{
			print $ach{$cols[$a]}, "/", $ach{$cols[$a]}, "\t";	
			}
		}
	print "\n";
	}

