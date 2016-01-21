#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
use Bio::SeqIO;

$scriptname=$0; $scriptname =~ s/.+\///g;

# -- program description and required arguments
unless ($#ARGV == 1)
        {print "\nCounts the number of occurences of the AlfI recognition site\n";
	print "or sequence motif in a set of DNA sequences.\n";
        print "Output:\t a fasta file of those sites, named by position.\n";
        print "Usage:\t $scriptname sequences output\n";
        print "Arguments:\n";
        print "\t sequences\t a fasta file containing the sequences to be searched \n";
        print "\t output\t\t a fasta file containing the extracted sites \n";
        print "\n"; exit;
        }
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
			while($ss =~ /(?=$patt)/g) 
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
	while($ss =~ /(?=$patt)/g) 
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

