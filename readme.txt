A collection of scripts for analysis of 2bRAD sequence data. 

-------------------------
AlleleFilter.pl
-------------------------

------------------------------------------------------------
AlleleFilter.pl
Excludes loci containing too many alleles.
Usage: AlleleFilter.pl -i input -n max_alleles <options>
Required arguments:
	-i input	name of the input file, a matrix of genotypes.
			Input format: rows=loci and columns=samples.
                        row 1 = column label, column 1 = tag, column 2 = position
                        subsequent columns contain genotypes for each sample
	-n max_alleles	maximum number of alleles allowed. Loci with more than this 
			number of allels will be excluded. 
Options:
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

------------------------------------------------------------

-------------------------
BuildRef.pl
-------------------------

------------------------------------------------------------
BuildRef.pl
Builds a reference for de novo analysis of 2bRAD sequences from samples
without a sequenced genome. The script filters, clusters, and compares 
similar sequences to infer the set of loci present in the species of
interest, using a subset of reads from the samples themselves.
Usage: BuildRef.pl -i input -o output <OPTIONS>
Required arguments:
	-i input	The set of processed (truncated, HQ) reads (FASTQ) to be used for reference
			development. Ideally this should include 10-20 million reads spanning the range
			of known diversity (e.g. from all populations in your study). Prepare this 
			ahead of time by concatenating together a subset of reads from your samples. 
	-o output	a name for the output file, to be used as a reference in mapping and genotyping
Options:
	-v overwrite	0=do not overwrite existing files; use them for analysis. (default) 
			1=do not use existing files; overwrite them with new files. 
	-n mincov	Minimum depth to qualify as a valid allele. (default=2)
	-q threshold	Quality scores below this threshold are low quality (default=30)
	-x number	Maximum number of low quality bases allowed for reference construction (default=0)
	-m mismatches	Maximum number of mismatches allowed in clustering of related alleles (default=2)
	-d distance	Minimum number of bases required to resolve sub-clusters (default=1)
	-a haplotypes	For very large clusters containing more than this number of unique sequences, 
			do not attempt to resolve sub clusters. These indicate repetitive sequences 
			which are not useful for genotyping anyway, and resolving these large clusters
			is computationally intensive. (default=32) 
------------------------------------------------------------

-------------------------
CallGenotypes.pl
-------------------------

------------------------------------------------------------
CallGenotypes.pl

Determines SNP genotypes from nucleotide frequencies. Input file contains nucleotide frequencies 
from multiple samples. Output file lists the genotypes called from those frequencies. 
Usage: CallGenotypes.pl -i input -o output <OPTIONS>
Required arguments:
	-i input	Input file, tab delimited text file of nucleotide frequencies (output from NtFrequences.pl)
			column = tag, column 2 = locus, column 3 = reference genotype
			subsequent columns = nucleotide frequences in each sample, as A/C/G/T (e.g. 0/0/10/12)
	-o output	A name for the output file (tab delimited text)
Options:
	-c coverage	Minimum coverage required to determine genotypes. Lower coverage loci wil be discarded.
			Default: 10
	-e ends		y: exclude terminal positions in alignments where errors may arise during ligation. (default)
			n: do not exclude terminal positions. 
	-m method	"nf": nucleotide frequencies (classic method; the default). 
			  This method determines genotypes directly from nucleotide frequencies, using thresholds
			  defined by the user. If minor allele frequency (MAF) <= x at a locus, the genotype is called 
			  homozygous for the major allele at that locus. If MAF >= n, genotype is called heterozygous.
			  Genotypes are not called at intermediate MAF (if x > MAF > n) where errors are likely.  
			"pgf": NF informed by population genotype frequencies. (an update on the classic method)
			  This method first identifies valid alleles at each locus based on their frequency in the 
			  population (the two most common alleles with frequencies >=y times in >= q individuals),
			  then applies relaxed nucleotide frequency thresholds for those alleles (using y instead
			  of n for valid alleles). 
			"bgc" = Bayesian Genotype Caller
			  This method calls the BGC software, which implements a maximum-likelihood (ML) method for 
			  calling genotypes that incorporates prior population-level information on genotype 
			  frequencies and error rates from a genotype-frequency estimator. For more details see 
			  (Maruki & Lynch, [doi: 10.1534/g3.117.039008], and cite that paper if using this option.	
	Options for method "nf" or "pgf":
	-x max_MAF	Maximum frequency of the minor allele you're willing to ignore and call the position 
			homozygous for the major allele (0-1). Default: 0.01
	-n min_MAF	Minimum frequency of the minor allele you're willing to accept as evidence of 
			heterozygosity, and call the locus heterozygous (0-1). Default: 0.25
	-r min_reads	Because low frequencies translate into 1 or fewer reads at low coverage, the script
			also imposes a minimum read number for detection of heterozygotes. (default: 2) 
	Options for method "pgf":
	-y frequency	Minimum frequency a second allele must be detected to be considered valid.
			(default: 0.05)
	-q samples	Each allele must present in at least q samples to be considered valid.
			(default: 2)	
	Options for method "bgc":
	-p p-value	Critical p-value for the chi-square polymorphism test (BGC)
			(default: 0.05)
	-v maxcov	Coverage at which the pipeline switches from BGC (for low coverage data) to
			HGC (for high coverage data). Default=80 (i.e. HGC above 80). 
Examples:
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab		# basic usage
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab -c 20 		# increase coverage threshold
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab -m bgc		# use Bayesian Genotype Caller
  CallGenotypes.pl -i allele_counts.tab -o genotypes.tab -m pgf -y 0.05	# use population method
------------------------------------------------------------

-------------------------
CombineAlleleCounts.pl
-------------------------

------------------------------------------------------------
CombineAlleleCounts.pl
Counts observations of alleles at each locus in a collection of base counts from 2bRAD 
(the output from SAMBaseCounts.pl). 

This script identifies the major and minor allele at each locus, and combines all samples
into a single file containing the number of times each of these alleles was observed
in each sample (two columns per sample, for major and minor allele respectively).

Output format: columns 1=tag, 2=position, 3=major allele, 4=minor allele,
5=major allele counts for sample A, 6=minor allele counts for sample A, etc.

Missing data are shown as "NA" for both alleles, and the minor allele is reported as 
"NA" for monomorphic loci.

Usage: CombineAlleleCounts.pl <options> file_1 file_2 ... file_n > output_file
Required arguments:
	files 1-n:	nucleotide frequencies (output from SAMBaseCounts.pl) for each sample
	output_file:	a name for the output; tab-delimited text
Options:
	-a max_alleles	maximum number of alleles allowed at each locus. Loci with more than this
			number of alleles will be excluded. (default=2)
        -v min_cov      minimum coverage required to consider an allele present. (default=2)
        -s min_samp     minimum number of samples in which an allele must be present
                        (default=1)

------------------------------------------------------------

-------------------------
CombineBaseCounts.pl
-------------------------

------------------------------------------------------------
CombineBaseCounts.pl
Counts the number of times each allele was observed, for each locus, in a collection of 
2bRAD data describing nucleotide frequencies for each locus and sample (the output from
SAMBaseCounts.pl). 

Output format: columns 1=tag, 2=position, 3=reference allele,
5=allele counts for sample 1 (A/C/G/T), 6=for sample 2, etc..
Missing data are shown as "NA" for all alleles.

Usage: CombineBaseCounts.pl file_1 file_2 ... file_n > output_file
Where:
	files 1-n:	nucleotide frequencies (output from SAMBasecaller.pl) for each sample
	output_file:	a name for the output; tab-delimited text
------------------------------------------------------------

-------------------------
CompareSNPMatrices.pl
-------------------------

------------------------------------------------------------
CompareSNPMatrices.pl

Compare two matrices of SNP genotypes (e.g. produced from RAD data) from the same set of 
samples to evaluate overlap in genotyped loci, and the level of agreement and disagreement in genotypes.
This is useful for comparing different genotyping algorithms. 
Input files are formatted as the output from NFGenotyper or BGCGenotyper: tab-delimited text, 
rows are loci, columns 1-2 are tag and locus and subsequent columns are samples, homozygotes 
shown as e.g. "A", heterozygotes as e.g. "A C", and missing data as "0". 

Usage: CompareSNPMatrices.pl -f file1 -s file2 <OPTIONS>
Required arguments:
	-f file1	name of the first SNP matrix (tab delimited text)
	-s file2	name of the second SNP matrix (tab delimited text)
Options:
	-o option	0: (default) don't print any detailed info on disagreements
			1: show detailed info on loci called different homozygous genotypes in each file
			2: show detailed info on loci called different heterozygous genotypes in each file
			3: show detailed info on loci called homozygous in file1 and heterozygous in file2 
			4: show detailed info on loci called homozygous in file2 and heterozygous in file1 
			5: show detailed info on loci called in file 1 but not in file 2
	-b counts	(required if -o > 0) the file of allele counts from which genotypes were called.
------------------------------------------------------------

-------------------------
EvalFrags.pl
-------------------------

------------------------------------------------------------
EvalFrags.pl
Evaluates the uniqueness of type IIb restriction fragments in a FASTA file.
e.g. a collection of 36-bp AlfI fragments extracted from a genome sequence 
using AlfI_Extract.pl
Usage: EvalFrags.pl input.fasta 
Where:
	input.fasta:	a collection of 36-bp sequences (FASTA)
------------------------------------------------------------

-------------------------
ExtractSites.pl
-------------------------

------------------------------------------------------------
ExtractSites.pl
Counts and extracts type IIb restriction fragments from a set of DNA sequences.
Output:  a fasta file of those sites, named by position.
Usage:   ExtractSites.pl -i input -o output
Required arguments:
         -i input	a fasta file containing the sequences to be searched
	 -e enzyme	choice of enzyme (AlfI, BsaXI, BcgI)
         -o output	name for the output file, a fasta file of those sites
------------------------------------------------------------

-------------------------
FastaStats.pl
-------------------------

------------------------------------------------------------
FastaStats.pl
Summarizes length statistics for a set of DNA sequences.
Usage: FastaStats.pl -i input -o output
Required arguments:
	-i input	name of the input file (FASTA)
	-o output	a name for the output file (TXT)
------------------------------------------------------------

-------------------------
gt2bayes.pl
-------------------------

------------------------------------------------------------
gt2bayes.pl
Converts a 2bRAD genotype matrix into the input format required by BayeScan.
Usage: gt2bayes.pl -i input -p pop.file -o output
Required arguments:
        -i input        tab-delimited genotype matrix, with rows=loci and columns=samples.
                        First two columns indicate tag and position respectively.
                        This format is the output from CallGenotypes.pl.
	-p pop.file	a tab-delimited text file showing which population each
			sample was drawn from. Formatted as: SampleName "\t" PopName "\n"
			Note -- make sure that sample names in this file are exactly identical
			to those shown in the first row of the genotype matrix.
	-o output 	a name for the BayeScan formatted output file.
------------------------------------------------------------

-------------------------
gt2colony.pl
-------------------------

------------------------------------------------------------
gt2colony.pl
Converts a genotype matrix (loci x samples) into the appropriate input format for COLONY.
Usage: gt2colony.pl -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. COLONY input format.
------------------------------------------------------------

-------------------------
gt2dadi.pl
-------------------------

------------------------------------------------------------
gt2dadi.pl

Converts a SNP matrix (produced from CallGenotypes.pl) into the format required by
the software DADI, described at: https://bitbucket.org/gutenkunstlab/dadi/wiki/DataFormats
Usage: gt2dadi.pl -i input -k key -r reference -o output
Where:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-k key		a tab-delimited text file associating each sample in the input
			with a population label. Alleles will be counted and reported
			by the population labels assigned in this file. Formated as:
			Sample_name	Population_name
	-r reference	Name of the reference file from which these SNPs were called (FASTA format)
	-o output	a name for the output file. Tab delimited text in DADI format.

------------------------------------------------------------

-------------------------
gt2fasta.pl
-------------------------

------------------------------------------------------------
gt2fasta.pl
Converts a genotype matrix (loci x samples) to a FASTA-formatted alignment.
Usage: gt2fasta.pl -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. FASTA alignment format.
------------------------------------------------------------

-------------------------
gt2fstat.pl
-------------------------

------------------------------------------------------------
gt2fstat.pl
Converts a genotype matrix (loci x samples) to FSTAT format.
Usage: gt2fstat.pl -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. FSTAT format. corresponding locus_key 
			and sample_key files are also produced.
------------------------------------------------------------

-------------------------
gt2phy.pl
-------------------------

------------------------------------------------------------
gt2phy.pl
Converts a genotype matrix (loci x samples) to a PHYLIP-formatted alignment.
Usage: gt2phy.pl -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. PHYLIP alignment format.
------------------------------------------------------------

-------------------------
gt2related.pl
-------------------------

------------------------------------------------------------
gt2related.pl
Converts a genotype matrix (loci x samples) into the appropriate input format for the R package related
Usage: gt2related.pl -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. RELATED format.
------------------------------------------------------------

-------------------------
gt2remlf90.pl
-------------------------

------------------------------------------------------------
gt2remlf90.pl
Converts a SNP genotype matrix (loci x samples) produced from 2bRAD genotyping
into the format required for the BLUPF90 family of programs for mixed models
and quantitative genetic analysis. See BLUPF90 manual for details of that format.
Usage: gt2remlf90.pl -i input -o output
Required arguments:
	-i input	Name of the input file, from CallGenotypes.pl.
			(rows=loci, columns=samples, columns 1 & 2 show tag name and position in tag)
	-o output	A name for the output file, in the format expected by BLUPF90 programs, e.g.
				sample0   02221022511020101020
				sample100 12221222221222200010
------------------------------------------------------------

-------------------------
gt2Rqtl.pl
-------------------------

------------------------------------------------------------
gt2Rqtl.pl
Converts a 2bRAD genotype matrix into the csv input format required by R/qtl.
Usage: gt2Rqtl.pl -i input -t traits -o output
Where:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-t traits	tab-delimited file of data on  traits, as
				"sample1	trait1	...	traitN"
			(note that sample names must match column headers in snps file)
	-m map		tab-delimited file of map positions as
				"marker	  LG	position"
	-o output	a name for the csv formatted output file.
------------------------------------------------------------

-------------------------
gt2snpmatrix.pl
-------------------------

------------------------------------------------------------
gt2snpmatrix.pl
Converts a genotype matrix (loci x samples) to a snp matrix, as described in the manual
for the R package diveRsity. This snp matrix can be converted to genepop format
using diveRsity's snp2gen function.
Usage: gt2snpmatrix.pl -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. SNP matrix format, input for snp2gen.
------------------------------------------------------------

-------------------------
gt2structure.pl
-------------------------

------------------------------------------------------------
gt2structure.pl
Converts a genotype matrix (loci x samples) into the appropriate input format for STRUCTURE.
Usage: gt2structure.pl -i input -o output
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. STRUCTURE input format.
------------------------------------------------------------

-------------------------
gt2vcf.pl
-------------------------

------------------------------------------------------------
gt2vcf.pl

Converts a genotype matrix (loci x samples) to a VCF file.
Usage: gt2vcf.pl -i input -r reference -o output <options>
Required arguments:
	-i input	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from CallGenotypes.pl.
	-o output	a name for the output file. FASTA alignment format.
	-r reference	Complete path to the reference file used to generate these genotypes (FASTA). 
Options:
	-f filters	a text file described filters applied to the genotypes. this information
			will be included in the VCF file. e.g.	
			  "MD	removed loci genotyped in <20 samples"

------------------------------------------------------------

-------------------------
LowcovSampleFilter.pl
-------------------------

------------------------------------------------------------
LowcovSampleFilter.pl

Excludes samples with too much missing data (genotypes called at too few loci)
Usage: LowcovSampleFilter.pl -i input -n min_data <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
	-n min_data	samples in which fewer than this number of loci were genotyped will be excluded
Options:
	-m mode		g=genotypes (default). 
			  Input file contains genotypes from individuals.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          subsequent columns contain genotypes for each sample
			a=allele counts.
			  Input file contains allele counts from pooled samples.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          column 3 = major allele, column 4 = minor allele
                          subsequent pairs of columns contain allele counts 
			  (major then minor) for each sample
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

------------------------------------------------------------

-------------------------
MDFilter.pl
-------------------------

------------------------------------------------------------
MDFilter.pl

Excludes loci containing too many missing data (genotyped in too few samples)
Usage: MDFilter.pl -i input -n min_data <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
	-n min_data	loci that were genotyped in fewer samples than this will be excluded
Options:
	-m mode		g=genotypes (default). 
			  Input file contains genotypes from individuals.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          subsequent columns contain genotypes for each sample
			a=allele counts.
			  Input file contains allele counts from pooled samples.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          column 3 = major allele, column 4 = minor allele
                          subsequent pairs of columns contain allele counts 
			  (major then minor) for each sample
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

------------------------------------------------------------

-------------------------
OneSNPPerTag.pl
-------------------------

------------------------------------------------------------
OneSNPPerTag.pl

Selects a single SNP from each tag in a matrix or genotypes or allele counts. 
Chooses the locus with the least missing data. 
Usage: OneSNPPerTag.pl -i input <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
Options:
	-m mode		g=genotypes (default). 
			  Input file contains genotypes from individuals.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          subsequent columns contain genotypes for each sample
			a=allele counts.
			  Input file contains allele counts from pooled samples.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          column 3 = major allele, column 4 = minor allele
                          subsequent pairs of columns contain allele counts 
			  (major then minor) for each sample
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

------------------------------------------------------------

-------------------------
PolyFilter.pl
-------------------------

------------------------------------------------------------
PolyFilter.pl
Excludes loci containing too few classes of genotypes or numbers of alleles (keeps polymorphic loci).
Usage: PolyFilter.pl -i input <OPTIONS>
Required arguments:
	-i input	name of the input file, a matrix of genotypes or allele counts
			(see -m for format)
Options:
	-g genotypes	minimum number of unique genotypes (for -m g) or alleles (for -m a) 
			required to consider a locus polymorphic (default=2)
	-m mode		g=genotypes (default). 
			  Input file contains genotypes from individuals.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          subsequent columns contain genotypes for each sample
			a=allele counts.
			  Input file contains allele counts from pooled samples.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          column 3 = major allele, column 4 = minor allele
                          subsequent pairs of columns contain allele counts 
			  (major then minor) for each sample
	-v min_cov	(for -m a) minimum coverage required to consider an allele present
			(default=2)
	-s min_samp	minimum number of samples in which an allele must be present
			(default=1)
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

------------------------------------------------------------

-------------------------
QualFilterFastq.pl
-------------------------

------------------------------------------------------------
QualFilterFastq.pl

Removes reads containing too many low quality basecalls from a set of short sequences
Output:  high-quality reads in FASTQ format
Usage:   QualFilterFastq.pl -i input -m min_score -x max_LQ -o output
Required arguments:
	-i input	raw input reads in FASTQ format
	-m min_score	quality scores below this are considered low quality (LQ)
	-x max_LQ	reads with more than this many LQ bases are excluded
	-o output	name for ourput file of HQ reads in FASTQ format
------------------------------------------------------------

-------------------------
RandomFastq.pl
-------------------------

------------------------------------------------------------
RandomFastq.pl
Draws the specified number of sequences randomly from a FASTQ sequence file.
Usage: RandomFastq.pl -i input -n num_seq -o output
Required arguments:
	-i input	name of the input file from which sequences will be randomly drawn.
	-n num_seq	number of sequences to draw
	-o output	a name for the output file (FASTQ)
------------------------------------------------------------

-------------------------
RepTagFilter.pl
-------------------------

------------------------------------------------------------
RepTagFilter.pl

Excludes tags containing too many SNPs, suggesting repetive regions of the genome
Usage: RepTagFilter.pl -i input -n max_snps <OPTIONS>
Required arguments:
	-i input	name of the SNP input file, a matrix of genotypes or allele counts
			(see -m for format)
			note: the script assumes the input only includes polymorphic loci
	-n max_snps	all SNPs from tags containing more than this number of SNPs will be excluded
Options:
	-m mode		g=genotypes (default). 
			  Input file contains genotypes from individuals.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          subsequent columns contain genotypes for each sample
			a=allele counts.
			  Input file contains allele counts from pooled samples.
			  Input format: rows=loci and columns=samples.
                          row 1 = column label, column 1 = tag, column 2 = position
                          column 3 = major allele, column 4 = minor allele
                          subsequent pairs of columns contain allele counts 
			  (major then minor) for each sample
	-p option	y=print filtered loci and summary; n=only print summary
			(default=n)
	-o output	a name for the output file (loci passing this filter) (required if -p y)

------------------------------------------------------------

-------------------------
SAMBaseCounts.pl
-------------------------

------------------------------------------------------------
SAMBaseCounts.pl
Counts nucleotide frequencies at each locus in a 2bRAD sequence data set.
Usage: SAMBaseCounts.pl -i input -r reference -o <OPTIONS>
Required arguments:
	-i input		input alignments, SAM format
	-r reference		reference used to generate the input alignments, FASTA format
	-o output		a name for the output file (tab delimited text)
Options:
	-c coverage		loci with lower coverage are discarded (default: 3)
------------------------------------------------------------

-------------------------
SAMFilter.pl
-------------------------

------------------------------------------------------------
SAMFilter.pl
Filters the alignments produced by mapping short reads against a reference,
excluding ambiguous, short, and weak matches.
NOTE: make sure that when a read matches multiple reference sequences (ambigous)
your mapper reports at least two alignments in the output. This is NOT the default 
behavior for some mappers, but is required to exclude ambiguous matches before genotyping.

Usage:  SAMFilter.pl -i input -m matches -o output <options>
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
------------------------------------------------------------

-------------------------
TruncateFastq.pl
-------------------------

------------------------------------------------------------
TruncateFastq.pl
Truncates a set of short reads in FASTQ format to keep the region specified
Usage:   TruncateFastq.pl -i input -s start -e end -o output
Required arguments:
	-i input	file of short reads to be filtered, fastq format
	-s start	beginning of the region to keep, nucleotide position
	-e end		end of the region to keep, nucleotide position
	-o output	a name for the output file (fastq format)
------------------------------------------------------------

