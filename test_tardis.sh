#!/bin/bash

# This script acts as a demonstration of TARDIS usage in a simple case
# It downloads a single-chromosome BAM file from the 1000 Genomes Project, and associated reference genome and SONIC files
# Downloading these files will take a while. If they are pre-downloaded, just comment out the WGET lines below
# It then runs TARDIS on this data set. Since it is a single-chrmosome file, TARDIS will print information that shows
# no data is found in all other chromosomes
# This script also assumes TARDIS is downloaded, compiled, and available in the PATH.
# Prerequsites:
#	TARDIS
#	wget
# 	gzip

TARDIS_EXE=`which tardis`

if [ -z ${TARDIS_EXE} ]
then
	echo "TARDIS binary is not in your PATH. Either move the TARDIS binary to your PATH, or add its location to the PATH environment variable."
else
	echo "Found TARDIS at " ${TARDIS_EXE}.
fi


WGET_EXE=`which wget`

if [ -z ${WGET_EXE} ]
then
	echo "WGET binary is not in your PATH. Either move the WGET binary to your PATH, or add its location to the PATH environment variable."
else
	echo "Found WGET at " ${WGET_EXE}.
fi

SAMTOOLS_EXE=`which samtools`

if [ -z ${SAMTOOLS_EXE} ]
then
	echo "SAMTOOLS binary is not in your PATH. Either move the SAMTOOLS binary to your PATH, or add its location to the PATH environment variable."
else
	echo "Found SAMTOOLS at " ${SAMTOOLS_EXE}.
fi

# Download NA11930 chromosome 20 BAM
${WGET_EXE} -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA11930/alignment/NA11930.chrom20.ILLUMINA.bwa.CEU.low_coverage.20130415.bam
# Download NA11930 chromosome 20 BAI
${WGET_EXE} -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA11930/alignment/NA11930.chrom20.ILLUMINA.bwa.CEU.low_coverage.20130415.bam.bai

# Download reference
${WGET_EXE} -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
${WGET_EXE} -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

# unzip reference

gunzip human_g1k_v37.fasta.gz
if [ $? -eq 0 ] 
then
	echo "Reference genome seems NOT to be downloaded correctly. Please rerun this script.";
	exit 1;
else
	echo "Reference genome seems to be downloaded correctly.";
fi


# Download SONIC file for the same reference
${WGET_EXE} -c https://github.com/BilkentCompGen/sonic-prebuilt/raw/master/human_g1k_v37.sonic

# these are the parameters. Commend out the WGET lines aboove if they are predownloaded
INPUTBAM=NA11930.chrom20.ILLUMINA.bwa.CEU.low_coverage.20130415.bam
REF=human_g1k_v37.fasta
SONIC=human_g1k_v37.sonic

${SAMTOOLS_EXE} quickcheck $INPUTBAM
if [ $? -eq 0 ] 
then
	echo $INPUTBAM " seems to be downloaded correctly according to samtools quickcheck";
else
	echo $INPUTBAM " seems NOT to be downloaded correctly according to samtools quickcheck. Please rerun this script.";
	exit 1;
fi

# ready to call TARDIS 

${TARDIS_EXE} -i $INPUTBAM --ref $REF --sonic $SONIC --out NA11930-chrom20-demo


