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

	git clone https://github.com/calkan/tardis.git --recursive

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

	Additional parameters for 10X Genomics Linked Reads (under development):

	--10x                      : Enable 10X Genomics Linked Reads mode.
	--output-hs                : Output the selected clusters homogeneity scores to the VCF file.

	Information:
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.


Converting output VCF file to BED
==============

	awk '! /\#/' out.vcf |\
	awk '{print $1"\t"($2-1)"\t"(substr($8,match($8,/SVLEN=[0-9]+/)+length("SVLEN="),RLENGTH-length("SVLEN="))+$2-1)}' > out.bed

Alternatively, use VCFlib: https://github.com/vcflib/vcflib


Running TARDIS for 10x Data (A Partial Manual)
==============
For running TARDIS on 10x data, I start by running this:

	./tardis \
	--input /mnt/storage1/projects/giab/10x/NA12878_hg19/NA12878_hg19_phased_possorted_bam.bam \
	--input /mnt/storage1/projects/giab/10x/NA12891_hg19/NA12891_hg19_phased_possorted_bam.bam \
	--input /mnt/storage1/projects/giab/10x/NA12892_hg19/NA12892_hg19_phased_possorted_bam.bam \
	--ref /home/calkan/projects/iman-tardis/genome.fa \
	--gaps aux/hg19.gap.bed \
	--dups aux/hg19.dups.bed \
	--reps aux/hg19.reps.bed \
	--sensitive \
	--rd 2
	--out sensitive/yay \
	--10x \
	--output-hs \
	--threads 32

This will run TARDIS in 10x mode, utilizing the long range information of the reads.

Before this command finishes execution, and once mrFAST finishes creation of all DIVET files (there should be as many DIVET files as there are libraries in your BAM files), I run this command:

	./tardis \
	--input /mnt/storage1/projects/giab/10x/NA12878_hg19/NA12878_hg19_phased_possorted_bam.bam \
	--input /mnt/storage1/projects/giab/10x/NA12891_hg19/NA12891_hg19_phased_possorted_bam.bam \
	--input /mnt/storage1/projects/giab/10x/NA12892_hg19/NA12892_hg19_phased_possorted_bam.bam \
	--ref /home/calkan/projects/iman-tardis/genome.fa \
	--gaps aux/hg19.gap.bed \
	--dups aux/hg19.dups.bed \
	--reps aux/hg19.reps.bed \
	--sensitive \
	--rd 2
	--out sensitive/yay \
	--skip-fastq \
	--skip-remap \
	--output-hs
	
This will run TARDIS without considering 10x barcode information. I use to compare against the result of using 10x barcode information. If 10x mode is to be used,  it is essential to have the DIVET files created by TARDIS that is running in 10x mode. This is because only in 10x mode does TARDIS write the 10x barcode information to the DIVET files.

If for some reason running with 10x flag was interrupted, but you nonetheless still have the DIVET files of TARDIS that was running in 10x mode, you can run this command that would save you having to recreate DIVET files:

	./tardis \
	--input /mnt/storage1/projects/giab/10x/NA12878_hg19/NA12878_hg19_phased_possorted_bam.bam \
	--input /mnt/storage1/projects/giab/10x/NA12891_hg19/NA12891_hg19_phased_possorted_bam.bam \
	--input /mnt/storage1/projects/giab/10x/NA12892_hg19/NA12892_hg19_phased_possorted_bam.bam \
	--ref /home/calkan/projects/iman-tardis/genome.fa \
	--gaps aux/hg19.gap.bed \
	--dups aux/hg19.dups.bed \
	--reps aux/hg19.reps.bed \
	--sensitive \
	--rd 2
	--out sensitive/yay \
	--skip-fastq \
	--skip-remap \
	--output-hs \
	--10x

Developer's Documentation On 10x Pipeline
==============
So TARDIS takes in bam files as input. In its sensitive mode (which are using here) it will read through the bam file, and find all reads that are discordants. Deciding which read pairs are discordant is done by:

	common.h:int is_concordant( bam1_core_t bam_alignment_core, int min, int max);

Then TARDIS will write the left and right pairs in a FASTQ file. The FASTQ file will have the following name:

	<sample_name>.<library_name>_remap_1.fastq.gz ---> the first of the pair
	<sample_name>.<library_name>_remap_2.fastq.gz ---> the second of the pair

The two files will be sorted and any mate-less pair will be discarded.
During FASTQ files creation, if 10x flag (--10x) is passed in the command line arguments, the 10x barcode (BX field) will be read by TARDIS. To save it in the FASTQ file, TARDIS first encodes the barcode into an unsigned long using:

	common.h:unsigned long encode_ten_x_barcode(char* source);

Since all barcodes are exactly 16 nucleotides long, they fit in the first 4 bytes of the unsigned long. If a read does not have a BX field, it will be given:

	(unsigned long) -1 = 18446744073709551615

value. In either case, this unsigned long value will be appended to the read name with width of 20 digits, zero padded. The read name will look like this:

	<qname><zeros padding><encoded barcode>

This is done at:

	if (params->ten_x == 1){
	    ten_x_barcode = encode_ten_x_barcode(bam_aux_get(bam_alignment, "BX"));
	    sprintf(qname + strlen(qname), "%020lu\0", ten_x_barcode);
	}
			
We use 20 digit padding since the largest unsigned long is 20 digit long.
Once all FASTQ files are created and sorted, TARDIS runs mrFAST for each pair FASTQ files.The number of threads that mrFAST will use is specified by the --threads argument passed to TARDIS command line. 

Note: If you are debugging, you might wanna save yourself time by NOT deleting FASTQ files. Use --skip-fastq command line flag if the FASTQ files are already there. If this flag is on, TARDIS will get infer the FASTQ file names from the BAM files, but skip their actual creation.
mrFAST will create a bunch of files, most pertinent among them to TARDIS are the DIVET files. Each pair of FASTQ files (i.e. each library) will have one DIVET file. Each DIVET line is discordant alignment location that suggests an SV event.
Once mrFAST aligning is complete, TARDIS will start loading those alignments from the DIVET files to apply set cover over them. When it loads a line, it passes it for parsing to:

	vh/vh_divethandler.h:struct DivetRow *vh_loadDivetRowFromString (struct ReadName *hash[], char *line, struct LibraryInfo *libInfo, int id);

Part of the parsing done here, is that the 10x barcode is extracted from the read name:

	sscanf(readName + strlen(readName) - 20, "%020lu%0", &ten_x_barcode);

Then, the library's ID, libInfo->libId , is then stored in the last byte of the unsigned long representing the 10x barcode for this read:

	if (ten_x_barcode != (unsigned long)-1){
	    ten_x_barcode = ten_x_barcode | ((unsigned long)libInfo->libId << (sizeof(unsigned long)-1)*8);
	}

This is done to ensure that no two libraries share a barcode (sharing barcode tag does not mean anything useful in terms of 10x barcode generation).
