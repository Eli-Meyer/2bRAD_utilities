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

# define variables. edit these if you understand the consequences.
my $seqfile = $ARGV[0];	# name of the input sequence file (fasta format)
my $qualfile = $ARGV[1];# name of the input quality score file (qual format)
my $outfile = $ARGV[2];	# a name for the final output file (cluster derived reference to be used for mapping)
my $mindepth = 0.028;	# min branch length and depth to select clades. 0.028 corresponds to 1 bp distance
my $initc = 0.944;	# 0.944 = 2 bp distance for initial clustering of alleles and paralogs
my $mincov = 2;		# min coverage needed in cluster to count a valid allele
my $site = "\\w{12}GCA\\w{6}TGC\\w{12}";	# recognition site in Perl regexp
my $qthd = 30;		# minimum allowed score; anything lower counts as low quality (LQ)
my $maxlq = 0;		# maximum allowed number of LQ positions
my $ow_opt = 0;		# decide whether to overwrite existing files. 1 = force overwrite; 0 = use existing files
my $minmems = 4;	# minimum number of sequences in a clade to consider tree building. Don't set below 4.
my $maxmems = 100;	# maximum number of sequences in a clade to consider tree building.

# fix sequence name characters because raxmlHP is picky
print "Checking input sequence names for \":\" characters...\n";
$checkseq = `head -n 1 $seqfile | grep \"\:\" -m 1`;
$checkqual = `head -n 1 $qualfile | grep \"\:\" -m 1`;
print "Finished checking input files.\n";
if (length($checkseq)>0) 
	{
	print "Found \":\" in $seqfile. Correcting...\n";
	system("perl -pi -e \"s/\:/-/g\" $seqfile");
	print "Finished.\n";
	}
if (length($checkqual)>0) 
	{
	print "Found \":\" in $qualfile. Correcting...\n";
	system("perl -pi -e \"s/\:/-/g\" $qualfile");
	print "Finished.\n";
	}

# -- counting input
$inseq = `grep \">\" -c  $seqfile`;
$inqual = `grep \">\" -c  $qualfile`;
chomp($inseq); chomp($inqual);

# -- screen for Ns in input reads
print "Checking for Ns in input...\n";
open(SF, $seqfile);
while(<SF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; next;}
	else	{
		if ($_ =~ /N/i)
			{
			$nhash{$name}++;
			$withn++;
			}
		}
	}
print "Finished.\n";

# -- designate positions to exclude from quality filtering
@excvec = qw{13 14 15 16 22 23 24};
foreach $e (@excvec) {$exch{$e}++;}

# -- stringently filter and truncate raw reads and scores
if ($ow_opt eq 0)
	{
	if (-e "VHQ.qual" && -e "VHQ.fasta")
		{
		print "VHQ.fasta and VHQ.qual found. Using existing files...\n";
		$outqual = `grep -c \">\" VHQ.qual`;
		$outseq = `grep -c \">\" VHQ.fasta`;
		chomp($outseq); chomp($outqual);
		$nlq = $inseq - $outseq - $withn;
		goto VHQ;
		}
	}
print "Conducting stringent quality filtering...\n";
my $nqual = $lqn = 0;
open(QF, $qualfile);
open(OQ, ">VHQ.qual");
while(<QF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; next;}
	else	{@qa = split(" ", $_);
		$lqn = 0; $poscount = 0;
		foreach $q (@qa) 
			{
			$poscount++; 
			if (exists($exch{$poscount})) {next;}
			if ($q < $qthd) {$lqn++;} 
			}
		if (exists($nhash{$name})) {next;}
		if ($lqn>$maxlq) {$nlq++; next;}
		elsif ($lqn <= $maxlq)
			{print OQ ">", $name, "\n";
			print OQ join(" ", @qa), "\n";
			$gh{$name}++;
			$outqual++;
			}
		}
	}
close(QF);

open(SF, $seqfile);
open(OS, ">VHQ.fasta");
while(<SF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; next;}
	else	{
		if (exists($gh{$name})) 
			{print OS ">", $name, "\n";
			print OS $_, "\n";
			$outseq++;
			}
		}
	}
close(SF);
close(OS);
print "Finished quality filtering.\n";
VHQ:					# end of quality filtering section. 	

print $inseq, " sequences in.\n";
print $inqual, " scores in.\n";
print $withn, " sequences with Ns excluded.\n";
print $nlq, " low quality sequences excluded.\n";
print $outseq, " sequences out.\n";
print $outqual, " scores out.\n";

# -- Filter for perfect match to restriction site
print "Filtering for perfect matches to restriction site...\n";
if ($ow_opt eq 0)
	{
	if (-e "matches.fasta")
		{
		print "matches.fasta found. Using existing file...\n";
		goto MATCHES;
		}
	}
system("grep -P \"$site\" VHQ.fasta -B 1 > matches.fasta");
system("perl -pi -e \"s/^--\n//g\" matches.fasta");
MATCHES:				# end of restriction site filtering section

print "Done filtering for site matches. Matches found:\n";
system("grep '>' matches.fasta -c");

# -- this step identifies sequences observed repeatedly in the VHQ dataset
# -- we call these "valid alleles"
# -- cluster reads at 1.0 and extract valid alleles based on depth
print "Beginning initial clustering step...\n";
if ($ow_opt eq 0)
	{
	if (-e "va.clstr")
		{
		print "va.clstr found. Using existing file...\n";
		goto CLSTR1;
		}
	}
system("date");
system("cd-hit-est -i matches.fasta -M 0 -d 0 -c 1 -o va >va.log");
print "Finished with initial clustering.\n";
CLSTR1:					# end of initial clustering section

print "Parsing initial clusters...\n";
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

# -- apply coverage filter to identify sequences identified more than $mincov times in the VHQ data
foreach $c (sort(keys(%ch)))
	{
	if ($ch{$c}<$mincov)	{next;}
	$pass++;
	$passh{$c}++;
	}

# -- identify representative sequences from valid clusters
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
		$status = 0;
		}
	}
close(IN);

# -- extract representative sequences from valid clusters
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

# cluster valid alleles to identify groups of related alleles that may represent paralogs or alleles
print "Beginning second clustering step...\n";
if ($ow_opt eq 0)
	{
	if (-e "gra.clstr")
		{
		print "gra.clstr found. Using existing file...\n";
		goto CLSTR2;
		}
	}
system("cd-hit-est -i filtered_va.fasta -M 0 -d 0 -c $initc -o gra >gra.log");
print "Finished second clustering step.\n";
CLSTR2:					# end of second clustering step
print "Parsing cluster output...\n";

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
		}
	}
print "Finished parsing.\n";
print $gra_c, " clusters (groups of related alleles) identified.\n";
system("date");

# -- build trees for each cluster in an effort to distinguish between paralogs and alleles
print "Testing for clusters containing multiple loci...\n";
foreach $c (sort(keys(%grah)))
	{			# count members in each cluster
	%cih = %{$grah{$c}};
	@cmems = sort(keys(%cih));
	$ncmems = @cmems;
# If fewer than 4 alleles were observed, it's not possible to use sequence relationships
# among alleles to infer the number of loci in the cluster. 
	if ($ncmems < $minmems)
		{
#		print $c, " not tested, too few sequences in the cluster.\n";
# pick representative seq from this cluster and add it to reference as a single locus
		$outh{$grareps{$c}}++;
		if ($ncmems < 2) {$mono_loc++;}
		else	{$poly_one++;}
		next;
		}
# If more than the maximum number of sequences were observed, this indicates a highly 
# repetitive cluster that we'd like to exclude because such loci will be error prone. 
	if ($ncmems > $maxmems)
		{
#		print $c, " not tested, too many sequences in the cluster.\n";
# pick representative seq from this cluster and add it to reference as a single locus
		$outh{$grareps{$c}}++;
		$poly_max++;
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
	RAXML:
	system("raxmlHPC -s tmp.fasta -p 123 -m GTRCAT -n out.tree > ml.log");
	system("sleep 2");

# parse tree to identify clades that differ by at least the critical 
# minimum difference ($mindepth variable)
	if (-e "RAxML_bestTree.out.tree") 
		{
		$treeobj = new Bio::TreeIO(-file=>"RAxML_bestTree.out.tree", -format=>"newick");
		$tree = $treeobj->next_tree;
		$rootnode = $tree->get_root_node;
		}
	else	{
		print "Tree not found. Re-building...\n";
		goto RAXML;
		}

# store id, depth, and branch length for each node
	foreach $node ($rootnode->get_all_Descendents)
		{
		if (defined($node->id))
			{
			@lineage = $tree->get_lineage_nodes($node);
			$cladeh{$lineage[-1]->internal_id}{$node->id}++;
			$deph{$lineage[-1]->internal_id} = $lineage[-1]->depth;
			$blenh{$lineage[-1]->internal_id} = $lineage[-1]->branch_length;
			}
		}

# identify number of loci in each cluster and representative sequence for each locus
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

