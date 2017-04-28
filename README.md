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

	git clone https://github.com/calkan/tardis.git --recursive

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
	--10x			   : Read 10x barcodes (BX field) and use the homogeneity function as a modifier to the cluster weight function
	--output-hs		   : Regardless of whether we are using 10x homogeneity score (HS), output the selected clusters HS scores to the vcf file
	--threads		   : An interger. The number of threads that multithreaded mrFAST will use
	--rp			   : The minimum read-pair support required for an SV cluster to be considered

Converting output VCF file to BED
==============

	awk '! /\#/' out.vcf |\
	awk '{print $1"\t"($2-1)"\t"(substr($8,match($8,/SVLEN=[0-9]+/)+length("SVLEN="),RLENGTH-length("SVLEN="))+$2-1)}' > out.bed

Running Tardis for 10x Data (A Partial Manual)
==============
For running Tardis on 10x data, I start by running this:

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

This will run Tardis in 10x mode, utilizing the long range information of the reads.

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
	
This will run Tardis without considering 10x barcode information. I use to compare against the result of using 10x barcode information. If 10x mode is to be used,  it is essential to have the DIVET files created by Tardis that is running in 10x mode. This is because only in 10x mode does Tardis write the 10x barcode information to the DIVET files.

If for some reason running with 10x flag was interrupted, but you nonetheless still have the DIVET files of Tardis that was running in 10x mode, you can run this command that would save you having to recreate DIVET files:

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
So Tardis takes in bam files as input. In its sensitive mode (which are using here) it will read through the bam file, and find all reads that are discordants. Deciding which read pairs are discordant is done by:

	common.h:int is_concordant( bam1_core_t bam_alignment_core, int min, int max);

Then Tardis will write the left and right pairs in a FASTQ file. The FASTQ file will have the following name:

	<sample_name>.<library_name>_remap_1.fastq.gz ---> the first of the pair
	<sample_name>.<library_name>_remap_2.fastq.gz ---> the second of the pair

The two files will be sorted and any mate-less pair will be discarded.
During FASTQ files creation, if 10x flag (--10x) is passed in the command line arguments, the 10x barcode (BX field) will be read by Tardis. To save it in the FASTQ file, Tardis first encodes the barcode into an unsigned long using:

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
Once all FASTQ files are created and sorted, Tardis runs mrFAST for each pair FASTQ files.The number of threads that mrFAST will use is specified by the --threads argument passed to Tardis command line. 

Note: If you are debugging, you might wanna save yourself time by NOT deleting FASTQ files. Use --skip-fastq command line flag if the FASTQ files are already there. If this flag is on, Tardis will get infer the FASTQ file names from the BAM files, but skip their actual creation.
mrFAST will create a bunch of files, most pertinent among them to Tardis are the DIVET files. Each pair of FASTQ files (i.e. each library) will have one DIVET file. Each DIVET line is discordant alignment location that suggests an SV event.
Once mrFAST aligning is complete, Tardis will start loading those alignments from the DIVET files to apply set cover over them. When it loads a line, it passes it for parsing to:

	vh/vh_divethandler.h:struct DivetRow *vh_loadDivetRowFromString (struct ReadName *hash[], char *line, struct LibraryInfo *libInfo, int id);

Part of the parsing done here, is that the 10x barcode is extracted from the read name:

	sscanf(readName + strlen(readName) - 20, "%020lu%0", &ten_x_barcode);

Then, the library's ID, libInfo->libId , is then stored in the last byte of the unsigned long representing the 10x barcode for this read:

	if (ten_x_barcode != (unsigned long)-1){
	    ten_x_barcode = ten_x_barcode | ((unsigned long)libInfo->libId << (sizeof(unsigned long)-1)*8);
	}

This is done to ensure that no two libraries share a barcode (sharing barcode tag does not mean anything useful in terms of 10x barcode generation).
