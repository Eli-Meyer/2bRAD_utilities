#!/usr/bin/env perl

# -- program name
print "-"x60, "\n";
print "BsaXI_Extract.pl v 3.0\n";
print "Created 02 Mar 2010 E Meyer\n";
print "Last modified 05 Jan 2011\n";
print "-"x60, "\n";

# -- program description and required arguments
unless ($#ARGV == 1)
        {print "Counts the number of occurences of the BsaXI recognition site\n";
	print "or sequence motif in a set of DNA sequences.\n";
        print "Output:\t a fasta file of those sites, namedd by position.\n";
        print "Usage:\t script sequences output\n";
        print "Arguments:\n";
        print "\t sequences\t a fasta file containing the sequences to be searched \n";
        print "\t output\t a fasta file of those sites \n";
        print "\n"; exit;
        }
my $seqfile = $ARGV[0];
open(SEQ, $seqfile);

$patt = ".{10}CGA.{6}TGC.{12}";
$rcpatt = ".{12}GCT.{5}ACG.{10}";
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
				$subi = substr($ss, $loci-$size, $size);
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_F\n";
				print OUT $subi, "\n";
				$found++;
				}
			while($ss =~ /(?=$rcpatt)/g) 
				{
				$loci = pos($ss)+$size;
				$subi = substr($ss, $loci-$size, $size);
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_R\n";
				print OUT $subi, "\n";
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
		$subi = substr($ss, $loci-$size, $size);
		print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_F\n";
		print OUT $subi, "\n";
		$found++;
		}
	while($ss =~ /(?=$rcpatt)/g) 
		{
		$loci = pos($ss)+$size;
		$subi = substr($ss, $loci-$size, $size);
		print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "_R\n";
		print OUT $subi, "\n";
		$found++;
		}
	$ss = "";
	}

print $nseqs, " sequences searched\n";
print $found, " recognition sites found altogether\n";

