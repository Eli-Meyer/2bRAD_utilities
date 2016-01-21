#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without restrictions or guarantees

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 4 || $ARGV[0] eq "-h") 
	{
	print "\nFilters the output from mapping short reads against a reference,\n";
	print "excluding short, weak, and ambiguous matches.\n";
	print "NOTE: make sure you've chosen settings for your mapper that output multiple\n";
	print "alignments for reads that match multiple reference sequences, because this is\n";
	print "required to test for ambiguous matches.\n";
	print "Usage:\t$scriptname input.sam min_matching min_aligned output.sam counts.tab\n";
	print "\tinput.sam: \tOutput from any short read mapper, SAM format.\n";
	print "\tmin_matching: \tMinimum number of matching bases.\n";
	print "\tmin_aligned: \tMinimum length of aligned region (bp).\n";
	print "\toutput.sam: \tName for output file containing alignments passing this filter, SAM format.\n";
	print "\tcounts.tab: \tName for the count output (reads per gene), tab delimited text.\n\n";
	exit;
	}
my $infile = $ARGV[0];
my $mthd = $ARGV[1];
my $athd = $ARGV[2];
my $outfile = $ARGV[3];
my $countfile = $ARGV[4];
my $ambig = 0;
my $tooshort = 0;

# read in sam file output and build a hash, counting raw mappings
open(IN, $infile);
my %maph;
while(<IN>)
	{
	if ($_ =~ /^@/) {next;}
	chomp;
	$rowi = $_;
	$rowi =~ s/^>//;
	@cols = split("\t", $rowi);
	$ncols = @cols;

	if ($cols[2] eq "*") {next;}
	$rawmap++;

# -- extract alignment length and apply alignment threshold
        $numstr = $cols[5];
        @chars = split("M", $numstr);
	$aligni = $mismatchi = 0;
#	print $numstr, "\t", "@chars", "\t";
	foreach $c (@chars)
		{
		$c =~ s/.+\D//g;
		if ($c > 0) {$aligni+=$c;}
#		print $c, " ";
		}
#	print "\n";
	if ($aligni < $athd) {$tooshort++; next;}


# -- extract mismatches and apply mismatch threshold
	foreach $flag (@cols[11..$ncols])
		{
		if ($flag =~ /NM\:i/)
			{
			$flag =~ s/NM\:i\://;
			$mismatchi = $flag;
			}
		}
	$matchi = $aligni - $mismatchi;
	if ($matchi < $mthd) {$tooweak++; next;}

# -- add extracted data to hash
	$maph{$cols[0]}{$cols[2]}{"count"}++;
	$maph{$cols[0]}{$cols[2]}{"string"} = $_;
	$maph{$cols[0]}{$cols[2]}{"align"} = $aligni;
	$maph{$cols[0]}{$cols[2]}{"match"} = $matchi;
#	print $_, "\n";
	}

# count number of reads with one or more matches passing thresholds
@mra = keys(%maph);
$nmr = @mra;

# select final set of unique mappings
open(OUT, ">$outfile");
open(IN, $infile);
while(<IN>)
	{
	chomp;
	if ($_ =~ /^@/) {print OUT $_, "\n";}
	else {last;}
	}

my %refch;
foreach $r (@mra)
	{
	%rh = %{$maph{$r}};
	@ma = keys(%rh);
	$nma = @ma;
	if ($nma>1)
		{
		@moa = sort{$rh{$b}{"match"}<=>$rh{$a}{"match"}}(keys(%rh));
		if ($rh{$moa[0]}{"match"}==$rh{$moa[1]}{"match"}) 
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
print $tooshort, " excluded for short alignments\n";
print $tooweak, " excluded for weak matches\n";
print $ambig, " excluded for ambiguous matches\n";
print $unimap, " unique mappings remained\n";

