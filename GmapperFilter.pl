#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 4 || $ARGV[0] eq "-h") 
	{
	print "\nFilters alignments produced by SHRiMP's gmapper and probcalc to exclude\n";
	print "weak or short, partial matches based on user specified criteria.\n";
        print "Ambiguous mapping is evaluated at the level of sequences in your reference.\n"; 
	print "If a read matches two sequences equally well, that read is discarded.\n\n";
	print "Usage: $scriptname input_alignments min_matching min_align output_alignments output_counts\n"; 
	print "Where: input_alignments:\tinput file, the output from gmapper -> probcalc.\n";
	print "\tmin_matching:\t\tMinimum number of matching bases required in each alignment.\n";
	print "\tmin_align:\t\tMinimum length of each alignment required.\n";
	print "\toutput_alignments:\tA name for the output file containing all alignments passing the above criteria.\n";
	print "\t(This file is further processed for 2bRAD genotyping)\n";
	print "\toutput_counts:\t\tA name for the output file documenting coverage for each reference sequence.\n";
	print "\t(This file is used for RNASeq gene expression analysis)\n\n";
	exit;
	}

my $infile = $ARGV[0];
my $mthd = $ARGV[1];
my $athd = $ARGV[2];
my $outfile = $ARGV[3];
my $countfile = $ARGV[4];

# read in gmapper output and build a hash, counting raw mappings
open(IN, $infile);
my %maph;
while(<IN>)
	{
	if ($_ =~ /^#/) {next;}
	chomp;
	$rowi = $_;
	$rowi =~ s/^>//;
	@cols = split("\t", $rowi);
	$rawmap++;

	$al = $cols[6] - $cols[5] + 1;
	if ($al < $mthd) {$tooshort++; next;}
        $numstr = $cols[9];
        @chars = split("", $numstr);
        my $mmno = 0;
        foreach $c (@chars)
                {
		if ($insw == 0)	{
		if ($c =~ /\d/) {next;}
		if ($c =~ /[xACGT-]/i) {$mmno++; next;}
		if ($c =~ /\(/) {$insw = 1; next;}
				}
		if ($insw == 1)	{
		if ($c =~ /[xACGT-]/i) {$mmno++; next;}
		if ($c =~ /\)/) {$insw = 0; next;}
				}	
                }
	$matchi = $al - $mmno;
	if ($matchi < $athd) {$tooshort++; next;}

	$maph{$cols[0]}{$cols[1]}{"count"}++;
	$maph{$cols[0]}{$cols[1]}{"string"} = $_;
	$maph{$cols[0]}{$cols[1]}{"score"} = $cols[8];
	$maph{$cols[0]}{$cols[1]}{"match"} = $matchi;
	$maph{$cols[0]}{$cols[1]}{"align"} = $al;
	}
# count number of reads with one or more matches
@mra = keys(%maph);
$nmr = @mra;
$ambig = 0;

# select final set of unique mappings
open(OUT, ">$outfile");
my %refch;
foreach $r (@mra)
	{
	%rh = %{$maph{$r}};
	@ma = keys(%rh);
	$nma = @ma;
	if ($nma>1)
		{
		@moa = sort{$rh{$b}{"score"}<=>$rh{$a}{"score"}}(keys(%rh));
		if ($rh{$moa[0]}{"score"}==$rh{$moa[1]}{"score"}) 
			{
			$ambig++; 
			next;
			}
		}
	if ($nma == 1) {@moa = @ma;}
	$unimap++;
	$refch{$moa[0]}++;
	print OUT $rh{$moa[0]}{"string"}, "\n";
	}

open(CTS, ">$countfile");
@sref = sort{$refch{$b}<=>$refch{$a}}(keys(%refch));
foreach $s (@sref)
	{
	print CTS $s, "\t", $refch{$s}, "\n";
	}

print $rawmap, " raw mappings altogether.\n";
print $nmr, " reads had one or more matches\n";
print $ambig, " excluded for ambiguous matches\n";
print $tooshort, " excluded for short partial matches\n";
print $unimap, " unique mappings remained\n";

