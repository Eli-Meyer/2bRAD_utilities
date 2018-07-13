#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Evaluates the uniqueness of type IIb restriction fragments in a FASTA file.
e.g. a collection of 36-bp AlfI fragments extracted from a genome sequence 
using AlfI_Extract.pl
Usage: $scriptname input.fasta 
Where:
	input.fasta:	a collection of 36-bp sequences (FASTA)
USAGE
if ($#ARGV != 0 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
# use this block if checking for executable dependencies
# copy the block and edit to check for additional Perl modules required by the script
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SeqIO;

my $infile = $ARGV[0];
my $seqs = new Bio::SeqIO(-file=>$infile, -format=>"fasta");
print "\n$infile\n", "-"x25, "\n";
my %srh;
my $maxrep = 0;
my $siteno = 0;
my $ncount = 0;
while($seq = $seqs->next_seq)
	{
	$siteno++;
	if ($seq->seq =~ /N/i) {$ncount++; next;}
	$srh{$seq->seq}++;
	$rco = $seq->revcom();
	$srh{$rco->seq}++;
	}
print $siteno, " sites altogether\n";
print $ncount, " excluded for ambiguous bases\n";

@usa = keys(%srh);
$nus = @usa / 2;
print $nus, " different sequences\n";

my %dh;
foreach $k (keys(%srh))
	{
	$dh{$srh{$k}}++;
	}
@nba = keys(%dh);
$nbn = @nba;

print "Csize\tNoC\t%Unique\t%Cov\n";
foreach $bin (sort{$a<=>$b}(keys(%dh)))
	{
	print $bin, "\t", $dh{$bin}/2, "\t";
	print int($dh{$bin}/2/$nus*100+0.5), "\t";
	print int($dh{$bin}/2*$bin/$nus*100+0.5), "\n";
	}
print "-"x25, "\n\n";

