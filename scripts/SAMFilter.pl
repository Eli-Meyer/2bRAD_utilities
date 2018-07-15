#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Filters the alignments produced by mapping short reads against a reference,
excluding ambiguous, short, and weak matches.
NOTE: make sure that when a read matches multiple reference sequences (ambigous)
your mapper reports at least two alignments in the output. This is NOT the default 
behavior for some mappers, but is required to exclude ambiguous matches before genotyping.

Usage:  $scriptname -i input -m matches -o output <options>
Required arguments:
	-i input	Output from any short read mapper, in SAM format.
	-m matches	Minimum number of matching bases required to consider an alignment valid. 
	-o output	A name for the filtered output (SAM format). 
Options:
	-c option	1: Report the number of reads matching each reference sequence
			in a separate output files "counts.tab". 0: Don't produce this file (default).
	-l length	Minimum length of aligned region (match, mismatch, + gaps) required to consider 
			an alignment valid. Only relevant if your mapper uses local alignment. For global
			alignments, this is set equal to -m. 
USAGE
if ($#ARGV < 3 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:m:o:c:l:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_m || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_c) {$cprint = 1;} else {$cprint = 0;}
if ($opt_l) {$athd = $opt_l;} else {$a_thd = $opt_m;}
my $infile = $opt_i;
my $mthd = $opt_m;
my $outfile = $opt_o;
$ambig = $tooshort = 0;

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
	foreach $c (@chars)
		{
		$c =~ s/.+\D//g;
		if ($c > 0) {$aligni+=$c;}
		}
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

if ($cprint eq 1)
{
open(CTS, ">counts.tab");
@sref = sort{$refch{$b}<=>$refch{$a}}(keys(%refch));
foreach $s (@sref)
	{
	print CTS $s, "\t", $refch{$s}, "\n";
	}
}

print $rawmap, " raw mappings altogether.\n";
print $nmr, " reads had one or more matches\n";
print $tooshort, " excluded for short alignments\n";
print $tooweak, " excluded for weak matches\n";
print $ambig, " excluded for ambiguous matches\n";
print $unimap, " unique mappings remained\n";

