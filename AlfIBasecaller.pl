#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
use Bio::SeqIO;

$scriptname=$0; $scriptname =~ s/.+\///g;

# -- program description and required arguments
unless ($#ARGV == 2)
        {print "\nCounts the number of alleles detected at each position in the specified set\n";
	print "of mapped short reads (output from SHRiMP prettyprint).\n";
        print "Output:\t A table listing the number of calls at each position, with the positions\n";
	print "\tcorresponding to the recognition site (AlfI) excluded.\n";
        print "Usage:\t usage: $scriptname reference shrimp.out min_coverage\n";
        print "Arguments:\n";
        print "\t reference \t the reference sequences, a fasta file \n";
        print "\t shrimp.out \t the output from SHRiMP (gmapper -> probcalc -> prettyprint)\n";
        print "\t min_coverage \t positions with lower coverage are ignored\n";
        print "\n"; exit;
        }

my $reffile = $ARGV[0];
my $refs = new Bio::SeqIO(-file=>$reffile, -format=>'fasta');
my $infile = $ARGV[1];
open(IN, $infile);
my $critcov = $ARGV[2];

# -- build a hash of reference sequence orientation
# -- note that this is not necessary for palindromic sites but is for non-palindromic sites
my %orih;
my $nref = 0;
my $fref = 0;
my $rref = 0;
SEQ: while($seq = $refs->next_seq)
	{$nref++;
	$str = $seq->seq;
	$str =~ tr/[a,c,g,t]/[A,C,G,T]/;
	@loci = split("", $str);
	for($i=-1; $i<2; $i++)	
		{
		$main5 = join("", @loci[12+$i..14+$i]); $main3 = join("", @loci[21+$i..23+$i]);
#		print $seq->display_id,  "\t", $str, "\t",$i, "\t", $main5, "\t", $main3, "\t", $alt5, "\t", $alt3, "\n";
		if (($main5 eq "GCA") && ($main3 eq "TGC")) {$orih{$seq->display_id} = "F"; $fref++; next SEQ;}
		}
	}

# -- loop through prettyprint output and build a hash of hashes of arrays
# where the first level of hash keys are reference sequences, the next level
# are positions within each reference, and the values are arrays of bases
# called for each position.
my %covh; my %nh; my $nosite = 0;
while (<IN>)
	{chomp($_);
	if ($_ =~ /^>/) 
		{
		@bits = split("\t", $_);
		$rid = $bits[0]; $rid =~ s/^>//;
		$fid = $bits[1]; $fid =~ s/\s.+//;
		$std = $bits[2]; 
		$nh{$fid}++;
		if (!exists($orih{$fid})) {$nosite++; next;}
		}
	if ($_ =~ /^G:/) 
		{@rloc = split(" ", $_); 
		$ri = $rloc[1]; 
		$rf = $rloc[3]; 
		$rfstr = $rloc[2];
		@rfarr = split("", $rfstr);
		$nfc = @rfarr;
		}
	if ($_ =~ /^R:/) 
		{@rloc = split(" ", $_); 
		$rdstr = $rloc[2]; 
		$rdstr =~ tr/acgt/N/;			
		@rdarr = split("", $rdstr);
#		print $fid, "\t", $rfstr, "\n";
#		print $rid, "\t", $rdstr, "\n";

		$bothpre = $rfpre = $rdpre = $bothpost = $rfpost = $rdpost = $sofar = 0;
		for ($i=1; $i<=$nfc; $i++) 
			{
			if ($sofar == 0)	{
			if ($rfarr[$i-1] eq "-" && $rdarr[$i-1] eq "-") {$bothpre++;}
			if ($rfarr[$i-1] eq "-" && $rdarr[$i-1] ne "-") {$rfpre++;}
			if ($rfarr[$i-1] ne "-" && $rdarr[$i-1] eq "-") {$rdpre++;}
			if ($rfarr[$i-1] ne "-" && $rdarr[$i-1] ne "-")	{$sofar++;}
						}
			else			{
			if ($rfarr[$i-1] eq "-" && $rdarr[$i-1] eq "-") {$bothpost++;}
			if ($rfarr[$i-1] eq "-" && $rdarr[$i-1] ne "-") {$rfpost++;}
			if ($rfarr[$i-1] ne "-" && $rdarr[$i-1] eq "-") {$rdpost++;}
			if ($rfarr[$i-1] ne "-" && $rdarr[$i-1] ne "-")
			 	{
				$bothpost = $rfpost = $rdpost = 0;
				}
						}
			}
#		print join("\t", ($bothpre, $rfpre, $rdpre, $bothpost, $rfpost, $rdpost)), "\n";
		$preoff = $bothpre + $rfpre + $rdpre;
		$postoff = $bothpost + $rfpost + $rdpost;
#		print $preoff, "\t", $postoff, "\n";
#		print $preoff, "\t", $postoff, "\n";

		if ($std eq "+")	{
		@locarr = (); 
		$poscount = 0; $igaps = 0;
#		print "\n";
		for ($i=$ri-$preoff; $i<=$rf+$postoff+$igaps; $i++) 
			{$poscount++;
			if ($rfarr[$poscount-1] eq "-" && $rdarr[$poscount-1] ne "-") 
				{$igaps++; 
#				print "- "; 
				push @locarr, "-"; 
				next;}
#			print $i-$igaps, " ";
			push @locarr, $i-$igaps;
			}
#		print "\n";
					}
		if ($std eq "-")	{
		@locarr = (); 
		$poscount = 0; $igaps = 0;
#		print "\n";
		for ($i=$ri+$preoff; $i>=$rf-$postoff-$igaps; $i--) 
			{$poscount++;
			if ($rfarr[$poscount-1] eq "-" && $rdarr[$poscount-1] ne "-") 
				{$igaps++; 
#				print "- "; 
				push @locarr, "-"; 
				next;}
#			print $i-$igaps, " ";
			push @locarr, $i-$igaps;
			}
#		print "\n";
					}

#		print "@locarr", "\n";

		%gth = (); $poscount = 0;
		foreach $l (@locarr) 
			{$poscount++; 
#			print $l, "\t", $rdarr[$poscount-1], "\n";
			if ($l eq "-") {next;}
			$gth{$l} = $rdarr[$poscount-1];
			}
		if ($std eq "+")
			{
			if ($orih{$fid} eq "F")	{
			foreach $l (sort{$a<=>$b}(keys(%gth)))
				{
				if ($l eq "-") {next;}
				if ($l==13||$l==14||$l==15||$l==22||$l==23||$l==24) {next;}
				if ($l<1||$l>36) {next;}
				$tcb = $gth{$l};	
				unless ($tcb eq "-") {push @{$covh{$fid}{$l}}, $tcb;}
				if ($tcb eq "-") {push @{$covh{$fid}{$l}}, "N";}
				}
						}

			if ($orih{$fid} eq "R")	{
			foreach $l (sort{$a<=>$b}(keys(%gth)))
				{
				if ($l eq "-") {next;}
				if ($l==13||$l==14||$l==15||$l==22||$l==23||$l==24) {next;}
				if ($l<1||$l>36) {next;}
				$tcb = $gth{$l};	
				unless ($tcb eq "-") {push @{$covh{$fid}{$l}}, $tcb;}
				if ($tcb eq "-") {push @{$covh{$fid}{$l}}, "N";}
				}
						}
			}
		if ($std eq "-")
			{
			if ($orih{$fid} eq "F")	{
			foreach $l (sort{$a<=>$b}(keys(%gth)))
				{
				if ($l eq "-") {next;}
				if ($l==13||$l==14||$l==15||$l==22||$l==23||$l==24) {next;}
				if ($l<1||$l>36) {next;}
				$tcb = $gth{$l};	
				$tcb =~ tr/ACGT/TGCA/;
				unless ($tcb eq "-") {push @{$covh{$fid}{$l}}, $tcb;}
				if ($tcb eq "-") {push @{$covh{$fid}{$l}}, "N";}
				}
						}
			if ($orih{$fid} eq "R")	{
			foreach $l (sort{$a<=>$b}(keys(%gth)))
				{
				if ($l eq "-") {next;}
				if ($l==13||$l==14||$l==15||$l==22||$l==23||$l==24) {next;}
				if ($l<1||$l>36) {next;}
				$tcb = $gth{$l};	
				$tcb =~ tr/ACGT/TGCA/;
				unless ($tcb eq "-") {push @{$covh{$fid}{$l}}, $tcb;}
				if ($tcb eq "-") {push @{$covh{$fid}{$l}}, "N";}
				}
						}
			}
		}
	}

#for $refi (keys %nh) {print $refi, "\t", $nh{$refi}, "\n";}

# -- loop through the hash of hashes and print out table
$lc = 0;
print "\tRef\tPosition\tA\tC\tG\tT\tN\n";
for $ref (keys %covh)
	{
	if ($nh{$ref} < $critcov) {next;}
	for $pos (sort{$a<=>$b}(keys %{$covh{$ref}}))
		{
		$lc++;
		$nA = 0; $nC = 0; $nG = 0; $nT = 0; $nN = 0; 
		print $lc, "\t", $ref, "\t", $pos, "\t";
		@ostr = @{$covh{$ref}{$pos}};
		for (@ostr) {if($_ eq "A") {$nA++;}}
		for (@ostr) {if($_ eq "C") {$nC++;}}
		for (@ostr) {if($_ eq "G") {$nG++;}}
		for (@ostr) {if($_ eq "T") {$nT++;}}
		for (@ostr) {if($_ eq "N") {$nN++;}}
		print $nA, "\t", $nC, "\t", $nG, "\t", $nT, "\t", $nN, "\n";
		}
	}	

