#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Summarizes length statistics for a set of DNA sequences.
Usage: $scriptname -i input -o output
Required arguments:
	-i input	name of the input file (FASTA)
	-o output	a name for the output file (TXT)
USAGE
if ($#ARGV < 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SeqIO;

# get variables from input
getopts('i:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $infile = $opt_i;
my $outfile = $opt_o;
my $seqs = new Bio::SeqIO(-file=>$infile, -format=>"fasta");

open(OUT, ">$outfile");
print OUT "\n";
print OUT $infile, "\n", "-"x25, "\n";
$countseq = 0;
$totbp = 0;
$ambigbp = 0;
$maxlen = 0;
$minlen = 1000000000000;
$countn = 0;
my %slh;
while($seq = $seqs->next_seq)
	{
	$lenbp = 0;
	$countseq++;
	$ss = $seq->seq;
	while($ss =~ /\w/g) {$lenbp++;}
	$slh{$seq->display_id} = $lenbp;
	foreach(split("",$ss)) 
		{
		if($_ !~ /[ACGT]/i){$ambigbp++;}
		if($_ =~ /N/i){$countn++;}
		}
	$totbp += $lenbp;		
	if($lenbp>$maxlen){$maxlen = $lenbp;}
	if($lenbp<$minlen){$minlen = $lenbp;}
	}
$meanbp = int($totbp/$countseq+0.5);
@ssia = sort{$slh{$a}<=>$slh{$b}}(keys(%slh));
my $tmpsum = 0;
foreach ($i=@ssia; $i>=0; $i--)
	{
	$tmpsum += $slh{$ssia[$i]};
	if ($tmpsum+$slh{$ssia[$i-1]}>=$totbp*0.5) {$n50=$slh{$ssia[$i]}; last;}
	}

print OUT $countseq, " sequences.\n";
print OUT $meanbp, " average length.\n";
print OUT $maxlen, " maximum length.\n";
print OUT $minlen, " minimum length.\n";
print OUT "N50 = ", $n50, "\n";

if($totbp<1000)
	{
	print OUT $totbp, " bp altogether.\n";
	print OUT $ambigbp, " ambiguous bp. ";
	print OUT "(", int($ambigbp/$totbp*1000+0.5)/10, "%)\n";
	print OUT $countn, " of those are Ns. ";
	print OUT "(", int($countn/$totbp*1000+0.5)/10, "%)\n";
	}
elsif ($totbp<1000000)	
	{
	print OUT int($totbp/100+0.5)/10, " kb altogether ($totbp bp).\n";
	print OUT int($ambigbp/100+0.5)/10, " ambiguous kb. ";
	print OUT "(", $ambigbp, " bp, ";
	print OUT int($ambigbp/$totbp*1000+0.5)/10, "%)\n";
	print OUT int($countn/100+0.5)/10, " kb of Ns. ";
	print OUT "(", $countn, " bp, ";
	print OUT int($countn/$totbp*1000+0.5)/10, "%)\n";
	}
elsif ($totbp>=1000000)
	{
	print OUT int($totbp/100000+0.5)/10, " Mb altogether ($totbp bp).\n";
	print OUT int($ambigbp/100000+0.5)/10, " ambiguous Mb. ";
	print OUT "(", $ambigbp, " bp, ";
	print OUT int($ambigbp/$totbp*1000+0.5)/10, "%)\n";
	print OUT int($countn/100000+0.5)/10, " Mb of Ns. ";
	print OUT "(", $countn, " bp, ";
	print OUT int($countn/$totbp*1000+0.5)/10, "%)\n";
	}
print OUT "-"x25, "\n\n";


