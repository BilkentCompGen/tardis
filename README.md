tardis
======

Toolkit for Automated and Rapid DIscovery of Structural variants

Requirements
============

 * zlib   (http://www.zlib.net)
 * mrfast (https://github.com/BilkentCompGen/mrfast)
 * htslib (included as submodule; http://htslib.org/)

Fetching tardis
===============

	git clone https://github.com/BilkentCompGen/tardis.git --recursive

Compilation
===========

Type:

	make libs
	make
	cp tardis /path/to/your/favorite/binaries


Auxiliary files
===============

GRCh37 annotations available under aux/

 * build37.dups.bed: Segmental duplication coordinates.
 * build37.gaps.bed: Assembly gap coordinates.
 * build37.reps.bed: RepeatMasker annotations (as described below). This file is provided as compressed. Unzip it before use.

Also download the reference genome from the UCSC Genome Browser. For GRCh37, this file is at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Merge all FASTA files into a single file. Make sure that the same reference was used to align the reads beforehand (BAM file).

** Reference genome and its annotations should use the EXACT same names for the chromosomes. The example provided below use the 1000 Genomes Project naming convention. **

Building the repeats file
=========================

Download the RepeatMasker out files from the UCSC Genome Browser. For GRCh37 (hg19), this file is at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromOut.tar.gz

Create a temp directory to make things easier to cleanup:

	mkdir rmasker
	tar zxf chromOut.tar.gz -C rmasker

Concetanate all out files and get rid of unneeded information. Make sure that the chromosome names in the reference genome and the repeats file matches.

	cat `find rmasker -name *.out` \
		| grep -v "matching\|begin"\
		| awk '{OFS="\t"; if ($1 ~ /./) print $5,$6-1,$7,$9,$10,$11}'\
		| sed s/chr// | sed s/Un_// | sed s/_random// | sed s/gl/GL/ | grep -v hap \
		| sed s/\.\._// | sed s/._// \
		| awk '{OFS="\t"; if ($1=="M") $1="MT"; if ($1 ~ /GL/) $1=$1".1"; print $0}'\
		| sort -k 1,1 -k 2,2n > build37.reps.bed

Remove unnecessary files:

	rm chromOut.tar.gz
	rm -fr rmasker


Running tardis
==============

	tardis -i myinput.bam --ref human_g1k_v37.fasta --gaps build37.gap.bed \
		--reps build37.reps.bed --dups build37_dups.bed --vh \
		--out myoutput

Additional parameters, helpful when debugging:

	--skip-fastq --skip-sort --skip-remap

All parameters
==============

	--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.
	--input [BAM files]        : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.
	--out   [output prefix]    : Prefix for the output file names.Output will be 3 files:
					<out>.vcf: Contains the predicted SV's
					<out>.name: Contains read names of all read pairs. Each line a read name.
					<out>.clusters: Contains all initially created SV clusters. Each line is a cluster.
	--ref   [reference genome] : Reference genome in FASTA format. Use the same reference you used for aligning in the BAM file. Make sure it contains all chromosomal contigs as those present in the gaps bed file (like chr11_gl000202_random or chr7_gl000195_random). The reference genome MUST have a mrFAST index!
	--gaps  [gaps file]        : Assembly gap coordinates in BED3 format. Choose the same version as that of the reference genome (hg19 or hg38).
	--dups  [dups file]        : Segmental duplication coordinates in BED3 format. Choose the same version as that of the reference genome (hg19 or hg38).
	--reps  [reps file]        : RepeatMasker annotation coordinates in BED6 format.
	--mei   ["Alu:L1Hs:SVA"]   : List of mobile element names.
	--skip-fastq               : Skip FASTQ dump for discordants. Use this only if you are regenerating the calls.
	--skip-sort                : Skip FASTQ sort for discordants. Use this only if you are regenerating the calls.
	--skip-remap               : Skip FASTQ remapping for discordants. Use this only if you are regenerating the calls.
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.
	--sensitive		   : Use multiple remappings of discordant reads to apply set cover on them.
	--output-hs		   : Regardless of whether we are using 10x homogeneity score (HS), output the selected clusters HS scores to the vcf file
	--threads		   : An interger. The number of threads that multithreaded mrFAST will use
	--rp			   : The minimum read-pair support required for an SV cluster to be considered

Converting output VCF file to BED
==============

	awk '! /\#/' out.vcf |\
	awk '{print $1"\t"($2-1)"\t"(substr($8,match($8,/SVLEN=[0-9]+/)+length("SVLEN="),RLENGTH-length("SVLEN="))+$2-1)}' > out.bed

