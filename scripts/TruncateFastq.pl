#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Truncates a set of short reads in FASTQ format to keep the region specified
Usage:   $scriptname -i input -s start -e end -o output
Required arguments:
	-i input	file of short reads to be filtered, fastq format
	-s start	beginning of the region to keep, nucleotide position
	-e end		end of the region to keep, nucleotide position
	-o output	a name for the output file (fastq format)
USAGE
if ($#ARGV < 3 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:s:e:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i ||!$opt_s ||!$opt_e ||!$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $seqfile = $opt_i;		# raw reads, fastq format
my $startloc = $opt_s-1;	# beginning of region to keep, base 1
my $endloc = $opt_e-1;		# end of region to keep, base 1
my $outfile = $opt_o;		# name for output file, fastq format

# loop through fastq file and truncate sequences and quality scores
open (IN, $seqfile);
open (OUT, ">$outfile");
my $switch = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 4);}
	if ($_ =~ /^$ss/) 
		{
		print OUT $_, "\n";
		next;
		}
	if ($_ =~ /^\+$/) 
		{
		print OUT "$_", "\n";
		next;
		}
	else
		{
		$ssi = $_;
		$slen = length($ssi);
		$subi = substr ($ssi, $startloc, $endloc-$startloc+1);
		print OUT $subi, "\n";
		}
	}
close(IN);
