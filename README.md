tardis
======

Toolkit for Automated and Rapid DIscovery of Structural variants

Requirements
============

 * zlib   (http://www.zlib.net)
 * mrfast (https://github.com/BilkentCompGen/mrfast)
 * htslib (included as submodule; http://htslib.org/)
 * sonic  (included as submodule; https://github.com/calkan/sonic)

Fetching TARDIS
===============

	git clone https://github.com/BilkentCompGen/tardis.git --recursive

Compilation
===========

Type:

	make libs
	make
	cp tardis /path/to/your/favorite/binaries


SONIC file (annotations container)
==================================

SONIC files are available under aux/

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
 	* Also download the reference genome at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz. 
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * ucsc_hg38.sonic: SONIC file for the human reference genome build 38.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.

Building the SONIC file
=======================

Please refer to the SONIC development repository: https://github.com/calkan/sonic/

However, you can still build the SONIC file using TARDIS:

	tardis --ref human_g1k_v37.fasta --make-sonic human_g1k_v37.sonic \
		--dups human_g1k_v37.segmental_duplications.bed \
		--gaps human_g1k_v37.assembly_gaps.bed \
		--reps human_g1k_v37.repeatmasker.out 

	

Running TARDIS - QUICK mode
===========================

	tardis -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--out myoutput

Additional parameters in SENSITIVE mode, helpful when debugging:

	--skip-fastq --skip-sort --skip-remap

All parameters
==============

	--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.
	--input [BAM files]        : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.
	--out   [output prefix]    : Prefix for the output file names.
	--ref   [reference genome] : Reference genome in FASTA format.
	--sonic [sonic file]       : SONIC file that contains assembly annotations.
	--mei   ["Alu:L1:SVA"]     : List of mobile element names.
	--no-soft-clip             : Skip soft clip remapping.

	Additional parameters for sensitive mode:

	--sensitive                : Sensitive mode that uses all map locations. Requires mrFAST remapping.
	--skip-fastq               : Skip FASTQ dump for discordants. Use this only if you are regenerating the calls. Sensitive mode only.
	--skip-sort                : Skip FASTQ sort for discordants. Use this only if you are regenerating the calls. Sensitive mode only. 
	--skip-remap               : Skip FASTQ remapping for discordants. Use this only if you are regenerating the calls. Sensitive mode only
	--threads                  : Number of threads for mrFAST to remap discordant reads.

	Additional parameters to build SONIC file within TARDIS:

	--make-sonic [sonic file]  : SONIC file that will contain the assembly annotations.
	--sonic-info [\"string\"]  : SONIC information string to be used as the reference genome name.
	--gaps  [gaps file]        : Assembly gap coordinates in BED3 format.
	--dups  [dups file]        : Segmental duplication coordinates in BED3 format.
	--reps  [reps file]        : RepeatMasker annotation coordinates in RepeatMasker format. See manual for details.

	Information:
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.


Converting output VCF file to BED
==============

	awk '! /\#/' out.vcf |\
	awk '{print $1"\t"($2-1)"\t"(substr($8,match($8,/SVLEN=[0-9]+/)+length("SVLEN="),RLENGTH-length("SVLEN="))+$2-1)}' > out.bed

Alternatively, use VCFlib: https://github.com/vcflib/vcflib


