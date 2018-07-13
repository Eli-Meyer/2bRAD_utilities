#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Draws the specified number of sequences randomly from a FASTQ sequence file.
Usage: $scriptname -i input -n num_seq -o output
Required arguments:
	-i input	name of the input file from which sequences will be randomly drawn.
	-n num_seq	number of sequences to draw
	-o output	a name for the output file (FASTQ)
USAGE
if ($#ARGV < 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:n:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_n || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# name variables from user input
$infile = $opt_i;
$todraw = $opt_n;
$outfile = $opt_o;

# count number of sequences in input
$nolines = `wc -l $infile`;
$nolines =~ s/\s.+//g;
chomp($nolines);
$noseqs = $nolines / 4;

if ($todraw >= $noseqs)
	{
	print "\nThere are only $noseqs sequences in $infile and you've asked to randomly draw $todraw.\n\n";
	exit;
	}

# generate list of sequence numbers to extract
%rh = ();
for ($a=1; $a<=$todraw; $a++)
	{
	$randi = 0;
	while ($randi eq 0 || exists($rh{$randi}))
		{
		$randi = int(rand($noseqs+1));
		}
	$rh{$randi}++;
#	print $randi, "\n";
	}

# loop through fastq file and print out randomly chosen sequence numbers
open (IN, $infile);
open (OUT, ">$outfile");
$state = 0; $seqcount = 0; $count = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 4);}
	if ($_ =~ /^$ss/) 
		{
		$seqcount++;
		if (exists($rh{$seqcount}))
			{$state = 1;}
		else 	{$state = 0;}
		}
	if ($state eq 1) {print OUT $_, "\n";}
	}
close(IN);

