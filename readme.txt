-------------------------
AlfIExtract.pl
-------------------------

------------------------------------------------------------
AlfIExtract.pl
Counts and extracts AlfI restriction fragments from a set of DNA sequences.
Output:  a fasta file of those sites, named by position.
Usage:   AlfIExtract.pl sequences output
Where:
         sequences      a fasta file containing the sequences to be searched
         output         a fasta file of those sites
------------------------------------------------------------

-------------------------
AlleleFilter.pl
-------------------------

------------------------------------------------------------
AlleleFilter.pl
Excludes loci containing too many alleles
Usage: AlleleFilter.pl input max_alleles print_option
Where:
                input:          genotypes matrix where rows=loci and columns=samples.
                                row 1 = column label, column 1 = tag, column 2 = position
                max_alleles:    loci with more than this number of alleles are excluded
                print_option:   y = print genotypes and summary, n = only summary

------------------------------------------------------------

-------------------------
BcgIExtract.pl
-------------------------

------------------------------------------------------------
BcgIExtract.pl
Counts and extracts BcgI restriction fragments from a set of DNA sequences.
Output:  a fasta file of those sites, named by position.
Usage:   BcgIExtract.pl sequences output
Where:
         sequences      a fasta file containing the sequences to be searched
         output         a fasta file of those sites
------------------------------------------------------------

-------------------------
BsaXIExtract.pl
-------------------------

------------------------------------------------------------
BsaXIExtract.pl
Counts and extracts BsaXI restriction fragments from a set of DNA sequences.
Output:  a fasta file of those sites, named by position.
Usage:   BsaXIExtract.pl sequences output
Where:
         sequences      a fasta file containing the sequences to be searched
         output  	a fasta file of those sites
------------------------------------------------------------

-------------------------
BuildRef.pl
-------------------------

------------------------------------------------------------
BuildRef.pl
Builds a reference for de novo analysis of 2bRAD sequences from samples
lacking a sequenced genome. Attempts to reconstruct the set of loci
(recognition sites) represented in a collection of raw 2bRAD sequences.
Usage: BuildRef.pl input.fasta input.qual output.fasta 
Where:
	input.fasta:	10-20 million truncated, high-quality reads from a
			representative subset of the samples
	input.qual:	quality scores for the input reads
	output.fasta:	a name for the output file, which will serve as a 
			reference for mapping and genotyping
------------------------------------------------------------

-------------------------
CombineGenotypes.pl
-------------------------

Combines output files produced by RAD_Genotyper.pl
Usage:	CombineGenotypes.pl multiple_input_files > output_file
Arguments:
	 multiple_input_files	 output files from RAD_genotyper.pl (wildcards can be used)
	 output_file		 a name for the output file

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
gt2bayes.pl
-------------------------

------------------------------------------------------------
gt2bayes.pl
Converts a 2bRAD genotype matrix into the input format required by BayeScan.
Usage: gt2bayes.pl input pop.file output
Where:
	input:		tab-delimited genotype matrix (columns=samples, rows=loci) 
			produced by RADGenotyper.pl or NFGenotyper.pl
	pop.file:	a tab-delimited text file showing which population each
			sample was drawn from. Formatted as: SampleName "\t" PopName "\n"
			Note -- make sure that sample names in this file are exactly identical
			to those shown in the first row of the genotype matrix.
	output: 	a name for the BayeScan formatted output file.
------------------------------------------------------------

-------------------------
gt2colony.pl
-------------------------

Converts a genotype matrix (loci x samples) into the appropriate input format for COLONY.
Usage: gt2colony.pl input > output
Where:	input:	tab-delimited genotype matrix, with rows=loci and columns=samples.
		First two columns indicate tag and position respectively.
		This format is the output from RADGenotyper.pl.
	output:	a name for the output file. COLONY input format.

-------------------------
gt2dadi.pl
-------------------------

------------------------------------------------------------
gt2dadi.pl
Converts a SNP matrix (produced from NFGenotyper.pl) into the format required by
the software DADI, described at: https://bitbucket.org/gutenkunstlab/dadi/wiki/DataFormats
Usage: gt2dadi.pl input key reference output
Where:
	input:		name of the input file (tab-delimited text)
	key:		a tab-delimited text file associating each sample in the input
			with a population label. Alleles will be counted and reported
			by the population labels assigned in this file. Formated as:
			Sample_name	Population_name
	reference:	Name of the reference file from which these SNPs were called (FASTA format)
	out:		a name for the output file. Tab delimited text in DADI format.
------------------------------------------------------------

-------------------------
gt2fas.pl
-------------------------

Converts a genotype matrix (loci x samples) to a FASTA-formatted alignment.
Usage: gt2fas.pl input > output
Where:	input:	tab-delimited genotype matrix, with rows=loci and columns=samples.
		First two columns indicate tag and position respectively.
		This format is the output from RADGenotyper.pl.
	output:	a name for the output file. FASTA alignment format.

-------------------------
gt2fstat.pl
-------------------------

Converts a genotype matrix (loci x samples) to FSTAT format.
Usage: gt2fstat.pl input > output
Where:	input:	tab-delimited genotype matrix, with rows=loci and columns=samples.
		First two columns indicate tag and position respectively.
		This format is the output from RADGenotyper.pl.
	output:	a name for the output file. FSTAT format.

-------------------------
gt2phy.pl
-------------------------

Converts a genotype matrix (loci x samples) to a PHYLIP-formatted alignment.
Usage: gt2phy.pl input > output
Where:	input:	tab-delimited genotype matrix, with rows=loci and columns=samples.
		First two columns indicate tag and position respectively.
		This format is the output from RADGenotyper.pl.
	output:	a name for the output file. PHYLIP alignment format.

-------------------------
gt2remlf90.pl
-------------------------

------------------------------------------------------------
gt2remlf90.pl
Converts a SNP genotype matrix (loci x samples) produced from 2bRAD genotyping
into the format required for the BLUPF90 family of programs for mixed models
and quantitative genetic analysis. See BLUPF90 manual for details of that format.
Usage: gt2remlf90.pl input > output
Where	input:	SNP matrix produced from CombineGenotypes.pl 
		(rows=loci, columns=samples, columns 1 & 2 show tag name and position in tag)
	output: format expected by BLUPF90 programs
		e.g.
		sample0   02221022511020101020
		sample100 12221222221222200010
------------------------------------------------------------

-------------------------
gt2snpmatrix.pl
-------------------------

Converts a genotype matrix (loci x samples) to a snp matrix, as described in the manual
for the R package diveRsity. This snp matrix can be converted to genepop format
using diveRsity's snp2gen function.
Usage: gt2snpmatrix.pl input > output
Where:	input:	tab-delimited genotype matrix, with rows=loci and columns=samples.
		First two columns indicate tag and position respectively.
		This format is the output from RADGenotyper.pl.
	output:	a name for the output file. SNP matrix format.

-------------------------
gt2structure.pl
-------------------------

Converts a genotype matrix (loci x samples) into the appropriate input format for STRUCTURE.
Usage: gt2structure.pl input > output
Where:	input:	tab-delimited genotype matrix, with rows=loci and columns=samples.
		First two columns indicate tag and position respectively.
		This format is the output from RADGenotyper.pl.
	output:	a name for the output file. STRUCTURE input format.

-------------------------
gt2vcf.pl
-------------------------

------------------------------------------------------------
gt2vcf.pl
Converts a genotype matrix (loci x samples) to a VCF file.
Usage: gt2vcf.pl input reference filters > output
Where:  input:  	tab-delimited genotype matrix, with rows=loci and columns=samples.
                	First two columns indicate tag and position respectively.
                	This format is the output from RADGenotyper.pl.
	reference:	Complete path to the reference file used to generate these genotypes (FASTA). 
	filters:	(OPTIONAL) a text file described filters applied to the genotypes, as:
				(code	description). e.g.	"MD	removed loci genotyped in <20 samples"
        output: a name for the output file. VCF format.
------------------------------------------------------------

-------------------------
LowcovSampleFilter.pl
-------------------------

Excludes samples with too much missing data (genotypes called at too few loci)
Usage: LowcovSampleFilter.pl input min_data print_option
Where:
		input: 		genotypes matrix where rows=loci and columns=samples.
				row 1 = column label, column 1 = tag, column 2 = position
		min_data: 	samples in which fewer loci than specified here were genotyped will be excluded
 		print_option: 	y = print genotypes and summary, n = only summary

-------------------------
MDFilter.pl
-------------------------

Excludes loci containing too many missing data (genotyped in too few samples)
Usage: MDFilter.pl input min_data print_option
Where:
		input: 		genotypes matrix where rows=loci and columns=samples.
				row 1 = column label, column 1 = tag, column 2 = position
		min_data: 	loci that were genotyped in fewer samples than this will be excluded
 		print_option: 	y = print genotypes and summary, n = only summary

-------------------------
NFGenotyper.pl
-------------------------

Determines genotypes from nucleotide frequencies based on user defined frequency thresholds.
Usage: NFGenotyper.pl input max_MAF_homo min_MAF_hetero min_cov
Where:	input:		Input file, tab delimited text file of nucleotide frequencies (output from AlfIBasecaller.pl)
	max_MAF_homo:	Loci at which the minor allele frequency (MAF) is below this threshold will be called
			homozygous for the major allele (i.e. the minor allele is considered seuqencing error and ignored).
	min_MAF_hetero:	Loci at which MAF exceeds this frequency will be called heterozygous for the two observed alleles.
			(Note: loci with MAF between these two thresholds will be called "unknown" and not reported).
	min_cov:	Minimum coverage required to determine genotypes. Lower coverage loci wil be discarded.

-------------------------
OneSNPPerTag.pl
-------------------------

Selects a single SNP from each tag in a 2bRAD genotype matrix to minimize missing data.
Usage: OneSNPPerTag.pl input print_option
Where:
		input: 		genotypes matrix where rows=loci and columns=samples.
				row 1 = column label, column 1 = tag, column 2 = position
 		print_option: 	y = print genotypes and summary, n = only summary

-------------------------
PolyFilter.pl
-------------------------

Excludes loci containing too few genotypes (keeps polymorphic loci)
Note that this counts genotypes, not alleles (e.g. AA, AB, BB = 3 genotypes)
Usage: PolyFilter.pl input min_gt print_option
Where:
		input: 		genotypes matrix where rows=loci and columns=samples.
				row 1 = column label, column 1 = tag, column 2 = position
		min_gt: 	minimum number of genotypes required to consider a locus polymorphic
 		print_option: 	y = print genotypes and summary, n = only summary

-------------------------
RepTagFilter.pl
-------------------------

Excludes tags containing too many SNPs, suggesting repetive regions of the genome
Usage: RepTagFilter.pl input max_snps print_option
Where:
		input: 		genotypes matrix where rows=loci and columns=samples.
				row 1 = column label, column 1 = tag, column 2 = position
		max_snps: 	tags containing more than this number of SNPs will be excluded
 		print_option: 	y = print genotypes and summary, n = only summary

-------------------------
SAMBasecaller.pl
-------------------------

Counts nucleotide frequencies at each locus in a 2bRAD sequence data set.
Usage: SAMBasecaller.pl input.sam reference.fasta cov_threshold output.tab
Where:
	input.sam:		input alignments, SAM format
	reference.fasta:	reference used to generate the input alignments, FASTA format
	cov_threshold:		loci with lower coverage are discarded
	output.tab:		a name for the output file (tab delimited text)

-------------------------
SamFilter.pl
-------------------------

Filters the output from mapping short reads against a reference,
excluding short, weak, and ambiguous matches.
NOTE: make sure you've chosen settings for your mapper that output multiple
alignments for reads that match multiple reference sequences, because this is
required to test for ambiguous matches.
Usage:	SamFilter.pl input.sam min_matching min_aligned output.sam counts.tab
	input.sam: 	Output from any short read mapper, SAM format.
	min_matching: 	Minimum number of matching bases.
	min_aligned: 	Minimum length of aligned region (bp).
	output.sam: 	Name for output file containing alignments passing this filter, SAM format.
	counts.tab: 	Name for the count output (reads per gene), tab delimited text.

