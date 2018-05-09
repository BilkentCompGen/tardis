#ifndef PROCESSREFGEN_H_
#define PROCESSREFGEN_H_

#include "common.h"
#include "processbam.h"


#define WINDOWSIZE 100
#define WINDOWSLIDE 100
#define STDEVLIMIT 5.0

/* Number of chromosomoes used in read depth analysis */
#define RDCHRCOUNT 24
#define CHRX 22
#define CHRY 23

int run_rd( bam_info** in_bam, parameters* param);
void init_rd( bam_info* in_bam, parameters* param);
void init_rd_per_chr( bam_info* in_bam, parameters* param, int chr_index);
void calc_mean( parameters *, bam_info *);
void calc_mean_per_chr( parameters *, bam_info *, int);

#endif /* PROCESSREFGEN_H_ */
