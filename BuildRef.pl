#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Builds a reference for de novo analysis of 2bRAD sequences from samples
lacking a sequenced genome. Attempts to reconstruct the set of loci
(recognition sites) represented in a collection of raw 2bRAD sequences.
Usage: $scriptname input.fasta input.qual output.fasta 
Where:
	input.fasta:	10-20 million truncated, high-quality reads from a
			representative subset of the samples
	input.qual:	quality scores for the input reads
	output.fasta:	a name for the output file, which will serve as a 
			reference for mapping and genotyping
USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
# use this block if checking for executable dependencies
# copy the block and edit to check for additional Perl modules required by the script
$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use File::Which;
$mod2="Bio::TreeIO";
unless(eval("require $mod2")) {print "$mod2 not found. Exiting\n"; exit;}
use Bio::TreeIO;

# use this block and edit to check for executables required by the script
$dep1 = "raxmlHPC";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}
$dep2 = "cd-hit-est";
unless (defined(which($dep2))) {print $dep2, " not found. Exiting.\n"; exit;}

# define variables
my $seqfile = $ARGV[0];	# name of the input sequence file (fasta format)
my $qualfile = $ARGV[1];# name of the input quality score file (qual format)
my $outfile = $ARGV[2];	# a name for the final output file (cluster derived reference to be used for mapping)
my $mindepth = 0.028;	# min branch length and depth to select clades. 0.028 corresponds to 1 bp distance
my $initc = 0.944;	# 0.944 = 2 bp distance for initial clustering of alleles and paralogs
my $mincov = 2;		# min coverage needed in cluster to count a valid allele
my $site = "\\w{12}GCA\\w{6}TGC\\w{12}";	# recognition site in Perl regexp
my $qthd = 30;		# minimum allowed score; anything lower counts as low quality (LQ)
my $maxlq = 2;		# maximum allowed number of LQ positions

#print $site, "\n";

# -- screen for Ns in input reads
open(SF, $seqfile);
while(<SF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; next;}
	else	{
		if ($_ =~ /N/i)
			{
			$nhash{$name}++;
			}
		}
	}

# -- stringently filter and truncate raw colorspace reads and scores
my $nqual = $lqn = 0;
open(QF, $qualfile);
open(OQ, ">VHQ.qual");
while(<QF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; $inqual++; next;}
	else	{@qa = split(" ", $_);
		$lqn = 0; $poscount = 0;
		foreach $q (@qa) 
			{
			$poscount++; 
			if ($poscount == 13 || $poscount == 14 || $poscount == 15 || $poscount == 22 || $poscount == 23 || $poscount == 24)
				{
				next;
				} 
			if ($q < $qthd) {$lqn++;} 
			}
		if (exists($nhash{$name})) {$nqual++; next;}
		if ($lqn>$maxlq) {$nlq++; next;}
		elsif ($lqn <= $maxlq)
			{print OQ ">", $name, "\n";
			print OQ join(" ", @qa), "\n";
			$gh{$name}++;
			$outqual++;
#			print $name, "\n";
			}
		}
	}
close(QF);

open(SF, $seqfile);
open(OS, ">VHQ.fasta");
while(<SF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; $inseq++; next;}
	else	{
#		print $name, "\n";
		if (exists($gh{$name})) 
			{print OS ">", $name, "\n";
			print OS $_, "\n";
			$outseq++;
			}
		}
	}
close(SF);
close(OS);
print $inseq, " sequences in.\n";
print $inqual, " scores in.\n";
print $nqual, " sequences with Ns excluded.\n";
print $nlq, " low quality sequences excluded.\n";
print $outseq, " sequences out.\n";
print $outqual, " scores out.\n";
print "Done filtering raw reads.\n";

# -- Filter for perfect match to restriction site
system("grep -P \"$site\" VHQ.fasta -B 1 > matches.fasta");
system("perl -pi -e \"s/^--\n//g\" matches.fasta");
print "Done filtering for site matches. Matches found:\n";
system("grep '>' matches.fasta -c");

# cluster reads at 1.0 and extract valid alleles based on depth
print "Beginning process...\n";
system("date");
system("cd-hit-est -i matches.fasta -M 0 -d 0 -c 1 -o va >va.log");
print "Finished with initial clustering. Parsing...\n";
open(IN, "va.clstr");
while(<IN>)
	{
	chomp;
	if ($_ =~ /^>/)
		{
		$cname = $_; 
		$cname =~ s/>//;
		$cname =~ s/ /_/;
		$cnum++;
		}
	else
		{
		$ch{$cname}++;
		$allreads++;
		}
	}
close(IN);

foreach $c (sort(keys(%ch)))
	{
	if ($ch{$c}<$mincov)	{next;}
#	print $c, "\t", $ch{$c}, "\n";	
	$pass++;
	$passh{$c}++;
	}

open(IN, "va.clstr");
$status = 0;
while(<IN>)
	{
	chomp;
	if ($_ =~ /^>/)
		{
		$cname = $_;
		$cname =~ s/ /_/;
		$cname =~ s/>//;
		if (exists($passh{$cname}))	{$status = 1;}
		}
	elsif ($status > 0 && $_ =~ /\*/)
		{
		$repname = $_;
		$repname =~ s/.*>//;
		$repname =~ s/\..+//;
		$reph{$repname}++;
#		print $repname, "\n";
		$status = 0;
		}
	}
close(IN);

open(IN,"va");
open(OUT, ">filtered_va.fasta");
$status = 0;
while(<IN>)
	{
	chomp;
	if ($_ =~ />/)
		{
		$rname = $_;
		$rname =~ s/>//;
		$rname =~ s/ .+//;
		if(exists($reph{$rname})) {print OUT $_, "\n"; $status++;}
		}
	elsif ($status > 0)
		{
		print OUT $_, "\n";
		$status = 0;		
		}
	}
close(IN);
close(OUT);
print "Finished parsing initial clusters. \n";
print "Identified $pass putative alleles (sequences observed at least $mincov times).\n";
system("date");
print "Beginning second clustering...\n";

# cluster valid alleles to identify groups of related alleles that may represent paralogs or alleles
system("cd-hit-est -i filtered_va.fasta -M 0 -d 0 -c $initc -o gra >gra.log");
print "Finished clustering. Parsing cluster output...\n";

# parse clusters to record cluster membership of each sequence
open(IN, "gra.clstr");
$status = 0; $gra_c = 0;
while(<IN>)
	{
	chomp;
	if ($_ =~ /^>/)
		{
		$cname = $_;
		$cname =~ s/ /_/;
		$cname =~ s/>//;
		$gra_c++;
		}
	else
		{
		$memname = $_;
		$memname =~ s/.*>//;
		$memname =~ s/\..+//;
		if ($_ =~ /\*/)
			{
			$grareps{$cname} = $memname;
			}
		$grah{$cname}{$memname}++;
#		print $memname, "\t", $cname, "\n";
		}
	}
print "Finished parsing.\n";
print $gra_c, " clusters (groups of related alleles) identified.\n";
system("date");

print "Testing for clusters containing multiple loci...\n";
# build trees for each cluster to distinguish between paralogs and alleles
foreach $c (sort(keys(%grah)))
	{
#	print $c, "\n";
	%cih = %{$grah{$c}};
	@cmems = sort(keys(%cih));
	$ncmems = @cmems;
#	print $ncmems, "\t", "@cmems", "\n";
# if fewer than 4 alleles observed, not possible to use sequence relationships among
# alleles to infer the number of loci in the cluster. Assign representative sequence
# for those clusters
	if ($ncmems < 4)
		{
#		print $c, " not tested, too few sequences in the cluster.\n";
		$outh{$grareps{$c}}++;
		if ($ncmems < 2) {$mono_loc++;}
		else	{$poly_one++;}
		next;
		}
# extract sequences of all cluster members
	$c4tree++;
	$found = 0;
	open(SEQ, "filtered_va.fasta");
	open(TMP, ">tmp.fasta");
	while(<SEQ>)
		{
		chomp;
		if ($_ =~ />/)
			{
			$status = 0;
			$_ =~ s/>//;
			$_ =~ s/\s.*//;
			if (exists($cih{$_}))
				{
#				print $_, "\n";
				$outname = $_;
				$outname =~ s/\:/-/g;
				print TMP ">", $outname, "\n";
				$status++;
				$found++;
				next;
				}
			}
		elsif ($status > 0)
			{
			print TMP $_, "\n";
			$status = 0;
			}
		if ($found >= $ncmems) {last;}
		}
	close(SEQ); close(TMP);
#	system("cat tmp.fasta");

# tree building code here

# prepare for tree building
	if(glob("*.tree"))	{system("rm *.tree");}
	%cladeh = %deph = %blenh = ();

# build tree describing relationships among members of the cluster
	system("raxmlHPC -s tmp.fasta -p 123 -m GTRCAT -n out.tree > ml.log");
#	system("cat RAxML_bestTree.out.tree");
#	system("cat RAxML_parsimonyTree.out.tree");

# parse tree to identify clades that differ by at least the critical 
# minimum difference ($mindepth variable)
	$treeobj = new Bio::TreeIO(-file=>"RAxML_bestTree.out.tree", -format=>"newick");
	$tree = $treeobj->next_tree;
	$rootnode = $tree->get_root_node;

# store id, depth, and branch length for each node
	foreach $node ($rootnode->get_all_Descendents)
		{
		if (defined($node->id))
			{
#			print $node->id, "\t";
			@lineage = $tree->get_lineage_nodes($node);
			foreach $l (@lineage)
				{
#				print $l->internal_id, "\t";
				}
#			print "\n";
#			print $lineage[-1]->internal_id, "\t";
#			print $lineage[-1]->height, "\t";
#			print $lineage[-1]->branch_length, "\t";
#			print $lineage[-1]->depth, "\n";
			$cladeh{$lineage[-1]->internal_id}{$node->id}++;
			$deph{$lineage[-1]->internal_id} = $lineage[-1]->depth;
			$blenh{$lineage[-1]->internal_id} = $lineage[-1]->branch_length;
			}
		}

# identify number of loci in each cluster and representative sequence for each locus
#	print "Clade\tRepresentative\n";
	%tmpfateh = ();
	foreach $tclade (sort(keys(%cladeh)))
		{
		$cladetest++; $cci = 0;
		%th = %{$cladeh{$tclade}};
		@tha = sort(keys(%th));
		$ntha = @tha;
# store a representative of clade at root 
		if ($deph{$tclade}==0)
			{
			$outh{$tha[0]}++;
			$tmpfateh{$tclade}++;
#			print $tclade, "\t", $tha[0], "\t", $ntha, "\t", "@tha", "\n";
			}	
# store a representative of each clade judged to be sufficiently different to
# represent a seperate locus
		elsif ($deph{$tclade}>=$mindepth && $blenh{$tclade}>=$mindepth)
			{
#			print $tclade, "\t", $tha[0], "\t", $ntha, "\t", "@tha", "\n";
			$outh{$tha[0]}++;
			$tmpfateh{$tclade}++;
			}
		}
# end tree building code
	@tmpkeys = keys(%tmpfateh);
	$nkeys = @tmpkeys;
	if ($nkeys > 1) {$splitclust++; $newclust += $nkeys;}
	else	{$remainsingle++;}
	}

print "Finished tree analysis.\n";
system("date");
print $mono_loc, " singleton clusters\n";
print $poly_one, " clusters made of two few sequences to build useful trees\n";
print $c4tree, " cluster with sufficient diversity to build trees\n";
print $remainsingle, " clusters remained as individual loci after tree analysis\n";
print $splitclust, " clusters were split into two or more loci\n";
print "Those clusters were split into $newclust new clusters based on tree analysis.\n";

@outlocs = sort(keys(%outh));
$nout = @outlocs;
print $nout, " loci identified altogether.\n";

# write out valid sites to a de novo reference file
open(IN, "filtered_va.fasta");
open(OUT, ">$outfile");
$status = 0; $found = 0;
while(<IN>)
	{
	chomp;
	if ($_ =~ /^>/)
		{
		$sid = $_;
		$sid =~ s/>//;
		$sid =~ s/ .*//g;
		$altsid = $sid;
		$altsid =~ s/\:/-/g;
		if (exists($outh{$sid}) || exists($outh{$altsid}))
			{
			$status = 1;
			$found++;
			print OUT ">denovoLocus", $found, "\n";
			next;
			}
		else
			{
			next;
			}
		}
	else	
		{
		if ($status == 1)
			{
			print OUT $_, "\n";
			}
		$status = 0;
		}
	}
close(IN);
close(OUT);
print "Wrote $found sequences to output file $outfile.\n";
print "This is now ready to use as a reference for mapping and genotyping.\n";
system("date");

