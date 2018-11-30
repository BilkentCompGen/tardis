tardis
======

Toolkit for Automated and Rapid DIscovery of Structural variants

Soylev, A., Kockan, C., Hormozdiari, F., & Alkan, C. (2017). Toolkit for automated and rapid discovery of structural variants. Methods, 129, 3-7. https://doi.org/10.1016/j.ymeth.2017.05.030

Soylev, A., Le, T., Amini, H., Alkan, C., & Hormozdiari, F. (2018). Discovery of tandem and interspersed segmental duplications using high throughput sequencing. bioRxiv, posted August 16, 2018. https://doi.org/10.1101/393694

TARDIS is developed and tested using Linux operating system (mostly Ubuntu flavors), and gcc versions 5.4 and 7.3. If it does not run as intended in other operating systems, we recommend to use the Docker container available at https://hub.docker.com/r/alkanlab/tardis/. We provide a sample script for Docker usage in this repository.

Requirements
============

 * gcc (version > 5)
 * zlib   (http://www.zlib.net)
 * htslib (included as submodule; http://htslib.org/)
 * sonic  (included as submodule; https://github.com/calkan/sonic)

htslib also requires:

 * libbz2
 * liblzma

Optional dependencies
=====================
 * mrFAST (https://github.com/BilkentCompGen/mrfast) - sensitive mode only.
 * gnuplot http://www.gnuplot.info) - for plotting fragment size distributions.

Fetching TARDIS
===============

	git clone https://github.com/BilkentCompGen/tardis.git --recursive

Compilation
===========

Type:

	make libs
	make
	cp tardis /path/to/your/favorite/binaries

Testing installation
====================

We provide a simple shell script in this repository (test_tardis.sh) to demonstrate a basic use case of TARDIS. This script downloads a single chromosome BAM file from the 1000 Genomes Project, and its associated reference genome and SONIC files. Note that these downloads may take some time to complete. The script then runs TARDIS with basic settings on this BAM file.

*Warning*: the downloaded files may be corrupted during transfer. This script does not check if the files are downloaded without problems.

To test, simply type:

	cd tardis
	sh ./test_tardis.sh


SONIC file (annotations container)
==================================

SONIC files for some human genome reference versions are available at external repo: https://github.com/BilkentCompGen/sonic-prebuilt

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


Running TARDIS - SENSITIVE mode (mrFAST Mappings)
=================================================

Sensitive mode uses mrFAST mappings (all possible mappings) with read-pair and read-depth signatures. 

	tardis -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--sensitive --out myoutput

This command first runs mrFAST and creates the DIVET file that contains all possible mappings of the reads in your BAM file. However, if you already have the DIVET files (there should be as many DIVET files as there are libraries in your BAM files), you can use --skip-mrfast

DIVET files should be inside the TARDIS directory or under the divet/ folder


Running TARDIS - Multiple BAM/CRAM files
========================================

You can run TARDIS on multiple input files to generate calls from a collection of samples. Note that we tested TARDIS with up to 25 input files, although hard limit is set to 100. 

There are three different ways of passing multiple input files to TARDIS:

* Use --bamlist with a text file. The text file should include paths to input BAM/CRAM files; one file per line:
	tardis --bamlist myinputs.txt --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic
		--out multiplesamples

* Use -i (or --input) multiple times:
        tardis -i myinput1.bam -i myinput2.bam -i myinput3.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic
                --out multiplesamples

* Use -i (or --input) with comma-separated BAM/CRAM list:
        tardis -i myinput1.bam,myinput2.bam,myinput3.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic
                --out multiplesamples


All parameters
==============

        Basic Parameters:
	--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.
	--input/-i [BAM files]     : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.
	--out   [output prefix]    : Prefix for the output file names.
	--ref   [reference genome] : Reference genome in FASTA format.
	--sonic [sonic file]       : SONIC file that contains assembly annotations.
	--hist-only                : Generate fragment size histograms only, then quit.

        Advanced Parameters:
	--read-cluster [int]       : # of clusters that a specific read can be involved in (Default is 10).
        --rp   [int]               : Minimum number of supporting read pairs in initial clustering (Default is 5).
	--mei   ["Alu:L1:SVA"]     : List of mobile element names.
	--no-soft-clip             : Skip soft clip remapping.
	--no-interdup              : Skip interspersed duplication clustering.
	--resolved                 : Output sequence resolved vcf calls.
	--xa                       : Look for the alternative mapping locations in BWA.
	--first-chr [chr_index]	   : Start running from a specific chromosome [0-based index in reference file]
	--last-chr  [chr_index]	   : Run up to a specific chromosome [0-based index in reference file]

	Additional parameters for sensitive mode:

	--sensitive                : Sensitive mode that uses all map locations. Requires mrFAST remapping.
	--skip-mrfast              : Skip mrFAST mapping. Use this only if you already have the correct divet file. Sensitive mode only
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


Docker Usage
============

To build a valor Docker image


	cd docker
	docker build . -t tardis:latest

Your image named "tardis" should be ready. You can run tardis using this image by

	docker run --user=$UID -v /path/to/inputs:/input -v /path/to/outputdir:/output tardis [args]

Alternatively, you can pull from Docker Hub:

	docker pull alkanlab/tardis

We also provide a test shell script for Docker usage based on DockerHub pull. This test is similar to the command line sample provided above, it first downloads test files.

	sh ./test_tardis_docker.sh

- ```[args]``` are usual arguments you would pass to valor executable. Be careful about mapping. You need to specify folders respective to container directory structure.
- You need to map host machine input and output directory to responding volume directories inside the container. These options are specified by '-v' argment.
- Docker works with root user by default. "--user" option saves your outputs.

Sample run:

	docker run --user=$UID -v /home/projects/tardis:/input -v /home/projects/tardis:/output tardis -i /input/my.input.grch38.bam --sonic /input/GRCh38.sonic --ref /input/human_v38.fasta --out /output/mydockertest

A pre-built docker image is also available at [Docker Hub](https://hub.docker.com/r/alkanlab/tardis/).

