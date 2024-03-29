tardis
======

Toolkit for Automated and Rapid DIscovery of Structural variants

Soylev, A., Kockan, C., Hormozdiari, F., & Alkan, C. (2017). Toolkit for automated and rapid discovery of structural variants. Methods, 129, 3-7. https://doi.org/10.1016/j.ymeth.2017.05.030

Soylev, A., Le, T., Amini, H., Alkan, C., & Hormozdiari, F. (2018). Discovery of tandem and interspersed segmental duplications using high throughput sequencing. Bioinformatics, Oct 15; 35(20): 3923–3930, 2019.

TARDIS is developed and tested using Linux operating system (mostly Ubuntu flavors), and gcc versions 5.4 and 7.3. If it does not run as intended in other operating systems, we recommend to use the Docker container available at https://hub.docker.com/r/alkanlab/tardis/. We provide a sample script for Docker usage in this repository.

Requirements
============

 * gcc (version > 5)
 * zlib   (development libraries; http://www.zlib.net)
	* Ubuntu example:
	sudo apt-get install zlib1g-dev

 * htslib (included as submodule; http://htslib.org/)
 * sonic  (included as submodule; https://github.com/calkan/sonic)

htslib also requires:

 * libbz2 (development libraries; http://www.bzip.org/)
	* Ubuntu example:
	sudo apt-get install libbz2-dev

 * liblzma (development libraries; https://tukaani.org/xz/)
	* Ubuntu example:
	sudo apt-get install liblzma-dev

Note that your BAM/CRAM files should be PCR free, otherwise you might get inaccurate results. Please eliminate PCR duplicates before running TARDIS with either Sambamba (https://lomereiter.github.io/sambamba/) or Picard (https://broadinstitute.github.io/picard/).

Optional dependencies
=====================
 * gnuplot (http://www.gnuplot.info) - for plotting fragment size distributions.

Fetching TARDIS
===============

	git clone https://github.com/BilkentCompGen/tardis.git --recursive

Installing development libraries in Ubuntu (requires sudo access)
======================================================

	sudo apt-get install zlib1g-dev liblzma-dev libbz2-dev

Compilation
===========

Type:

	make libs
	make
	cp tardis /path/to/your/favorite/binaries


Compiling without sudo access
=============================

If you do not have root access to install liblzma and/or libbz2, you can compile htslib without CRAM support. Note that lzma and libbz2 are prerequisites for htslib. The libz library is still required, talk to your admin if it is not available on your system.

Type:

	make nocram
	cp tardis-nocram /path/to/your/favorite/binaries

Note that this will disable CRAM support and you will be able to run TARDIS only with BAM files, and the name for the executable becomes tardis-nocram.


Testing installation
====================

We provide a simple shell script in this repository (test_tardis.sh) to demonstrate a basic use case of TARDIS. This script downloads a single chromosome BAM file from the 1000 Genomes Project, and its associated reference genome and SONIC files. Note that these downloads may take some time to complete. The script then runs TARDIS with basic settings on this BAM file.

*Warning*: the downloaded files may be corrupted during transfer. This script does not check if the files are downloaded without problems.

To test, simply type:

	cd tardis
	sh ./test_tardis.sh


Running TARDIS
===========================

	tardis -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--out myoutput


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


SONIC file (annotations container)
==================================

SONIC files for some human and mouse genome reference versions are available at external repo: https://github.com/BilkentCompGen/sonic-prebuilt

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
	* Also download the reference genome at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz. 
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * ucsc_hg38.sonic: SONIC file for the human reference genome build 38.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * mm9.sonic: SONIC file for the mouse reference genome version mm9.
	* Also download the reference genome at: http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * mm10.sonic: SONIC file for the mouse reference genome version mm10.
	* Also download the reference genome at: http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.


Building the SONIC file
=======================

Please refer to the SONIC development repository: https://github.com/calkan/sonic/

The README.md file includes documentation on how to obtain the necessary files for different genomes from the UCSC Genome Browser.

*Please note that you can also generate the SONIC file while running TARDIS:*
	
 	tardis -i myinput.bam --ref GRCh38.fa --gaps hg38.gap.bed --reps hg38.repeats.out --dups hg38.dups.bed \
	--make-sonic my_sonic.sonic --out myoutput


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
 
	--interdup                 : Run interspersed duplication clustering.
	--read-cluster [int]       : # of clusters that a specific read can be involved in (Default is 20).
	--rp   [int]               : Minimum number of supporting read pairs in initial clustering (Default is 5).
	--mei   [string	]          : List of mobile element names separated by colon (Default is ["Alu:L1:SVA"])
	--no-soft-clip             : Skip soft clip remapping.
	--no-mei                   : Skip mobile element insertion (MEI) clustering.
	--resolved                 : Output sequence resolved vcf calls.
	--first-chr [int]	   : Start running from a specific chromosome [0-based index in reference file]
	--last-chr  [int]	   : Run up to a specific chromosome [0-based index in reference file]

	Additional parameters to build SONIC file within TARDIS:

	--make-sonic [sonic file]  : SONIC file that will contain the assembly annotations.
	--sonic-info [string]      : SONIC information string to be used as the reference genome name, e.g., hg19.
	--gaps  [gaps file]        : Assembly gap coordinates in BED3 format.
	--dups  [dups file]        : Segmental duplication coordinates in BED3 format.
	--reps  [reps file]        : RepeatMasker annotation coordinates in RepeatMasker format. See manual for details.
	
	Information:
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.


Converting output VCF file to BED
==============
awk '! /\#/' tardis_sim_30.vcf | grep "<INV>" | awk '{print $1"\t"($2-1)"\t"(substr($8,match($8,/END=[0-9]+/)+length("END="),RLENGTH-length("END=")))}'
	awk '! /\#/' out.vcf |\
		awk '{print $1"\t"($2-1)"\t"(substr($8,match($8,/END=[0-9]+/)+length("END="),RLENGTH-length("END=")), $5)}' > out.bed

Alternatively, use VCFlib: https://github.com/vcflib/vcflib


Docker Usage
============

To build a TARDIS Docker image


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

