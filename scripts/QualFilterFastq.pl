#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Removes reads containing too many low quality basecalls from a set of short sequences
Output:  high-quality reads in FASTQ format
Usage:   $scriptname -i input -m min_score -x max_LQ -o output
Required arguments:
	-i input	raw input reads in FASTQ format
	-m min_score	quality scores below this are considered low quality (LQ)
	-x max_LQ	reads with more than this many LQ bases are excluded
	-o output	name for ourput file of HQ reads in FASTQ format
USAGE
if ($#ARGV < 3 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SeqIO;

# get variables from input
getopts('i:m:x:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i ||!$opt_m ||!$opt_x ||!$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $fastqfile = $opt_i;
my $lowq = $opt_m;
my $minlq = $opt_x;
my $outfqfile = $opt_o;

#my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq-illumina");
my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq");

my %sh; my $scount = 0;
while ($seq = $inseqs->next_seq) 
	{
	$scount++;
	$qo = new Bio::Seq::Quality(-accession_number=>$seq->display_id, -qual=>$seq->qual,
				-verbose=>-1);
	$qot = $qo->qual_text;
	@qoa = split(" ", $qot);
	$qid = $qo->accession_number;
	$lqcount = 0;
	foreach $q (@qoa)
		{	
		if ($q - $const < $lowq) {$lqcount++;}
		}
	if ($lqcount > $minlq) {$toolow++; next;}
	$gh{$qid}++;
	}

print "Output from ", $scriptname, "\n";
print "Checked ", $scount, " reads.\n";
print $toolow, " failed.\n";
print $scount - $toolow, " passed.\n";
print $toolow/$scount, " rejection rate.\n";

print "Writing sequences to output...\n";
#my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq-illumina");
#my $outseqs = new Bio::SeqIO(-file=>">$outfqfile", -format=>"fastq-illumina");
my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq");
my $outseqs = new Bio::SeqIO(-file=>">$outfqfile", -format=>"fastq");
$ocount = 0;
while ($seq = $inseqs->next_seq) 
	{
	$qo = new Bio::Seq::Quality(-accession_number=>$seq->display_id, -qual=>$seq->qual,
				-verbose=>-1);
	$qid = $qo->accession_number;
	if (exists($gh{$qid})) {$outseqs->write_seq($seq); $ocount++;};
	}

print "Done.\n";
print $ocount, " sequences written to output.\n";
system("date");
