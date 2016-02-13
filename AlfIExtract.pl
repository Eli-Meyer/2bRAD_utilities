#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Counts and extracts AlfI restriction fragments from a set of DNA sequences.
Output:  a fasta file of those sites, named by position.
Usage:   $scriptname sequences output
Where:
         sequences      a fasta file containing the sequences to be searched
         output         a fasta file of those sites
USAGE
if ($#ARGV != 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $seqfile = $ARGV[0];
open(SEQ, $seqfile);
$patt = ".{12}GCA.{6}TGC.{12}";
$size = 36;

print "Forward sequence of recognition site: ", $patt, "\n";
print "Reverse compliment sequence of recognition site: ", $rcpatt, "\n";
my $outfile = $ARGV[1];
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
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_F\n";
				print OUT substr($ss, $loci-$size, $size), "\n";
				$found++;
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
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_F\n";
				print OUT substr($ss, $loci-$size, $size), "\n";
				$found++;
		}
	$ss = "";
	}

print $nseqs, " sequences searched\n";
print $found, " recognition sites found altogether\n";


