#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Counts and extracts type IIb restriction fragments from a set of DNA sequences.
Output:  a fasta file of those sites, named by position.
Usage:   $scriptname -i input -o output
Required arguments:
         -i input	a fasta file containing the sequences to be searched
	 -e enzyme	choice of enzyme (AlfI, BsaXI, BcgI)
         -o output	name for the output file, a fasta file of those sites
USAGE
if ($#ARGV < 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:e:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || !$opt_e || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $seqfile = $opt_i;
open(SEQ, $seqfile);
if ($opt_e =~ /AlfI/i)	
	{
	$patt = ".{12}GCA.{6}TGC.{12}";
	$rcpatt = ".{12}GCA.{6}TGC.{12}";
	$size = 36;
	$type="pal";
	}
elsif ($opt_e =~ /BcgI/i)	
	{
	$patt = ".{10}CGA.{6}TGC.{12}";
	$rcpatt = ".{12}GCA.{6}TCG.{10}";
	$size = 36;
	$type="non"
	}
elsif ($opt_e =~ /BsaXI/i)	
	{
	$patt = ".{12}AC.{5}CTCC.{10}";
	$rcpatt = ".{10}GGAG.{5}GT.{12}";
	$size = 33;
	$type="non";
	}
else {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

print "Forward sequence of recognition site: ", $patt, "\n";
if ($type eq "non") {print "Reverse compliment sequence of recognition site: ", $rcpatt, "\n";}
my $outfile = $opt_o;
open(OUT, ">$outfile");

my $found = 0;
my $nseqs = 0;
$ss = "";
while(<SEQ>)
	{
	chomp;
	if ($_ =~ />/)
		{
		if ($ss ne "")
			{
			while($ss =~ /(?=$patt)/ig) 
				{
				$loci = pos($ss)+$size;
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci;
				if ($type eq "non") {print OUT "_F",;}
				print OUT "\n";
				print OUT substr($ss, $loci-$size, $size), "\n";
				$found++;
				}
			if (type eq "non")
			{
			while($ss =~ /(?=$rcpatt)/ig) 
				{
				$loci = pos($ss)+$size;
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_R",;
				print OUT "\n";
				print OUT substr($ss, $loci-$size, $size), "\n";
				$found++;
				}
			}
			$ss = "";
			}
		$nseqs++;
		$idi = $_; $idi =~ s/>//; $idi =~ s/\s+.*//;
		}
	else
		{
		$ss = $ss.$_;
		}
	}
if ($ss ne "")
	{
	while($ss =~ /(?=$patt)/ig) 
		{
		$loci = pos($ss)+$size;
		print OUT ">", $idi, "_", $loci-$size+1, "_", $loci;
		if ($type eq "non") {print OUT "_F",;}
		print OUT "\n";
		print OUT substr($ss, $loci-$size, $size), "\n";
		$found++;
		}
	if (type eq "non")
	{
	while($ss =~ /(?=$rcpatt)/ig) 
		{
		$loci = pos($ss)+$size;
		print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_R",;
		print OUT "\n";
		print OUT substr($ss, $loci-$size, $size), "\n";
		$found++;
		}
	}
	$ss = "";
	}

print $nseqs, " sequences searched\n";
print $found, " recognition sites found altogether\n";


