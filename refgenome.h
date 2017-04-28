#ifndef REFGENOME_H_
#define REFGENOME_H_

#include "common.h"
#include "processbam.h"


#define WINDOWSIZE 100
#define WINDOWSLIDE 100
#define STDEVLIMIT 5.0

/* Number of chromosomoes used in read depth analysis */
#define RDCHRCOUNT 24
#define CHRX 22
#define CHRY 23

int run_rd( bam_info** in_bam, parameters* param, ref_genome* ref);
int load_refgen(ref_genome** ref, parameters* params);
void calc_gc(ref_genome** ref);
void calc_gc_per_chr( ref_genome** ref, int chr_index);
void init_rd( bam_info* in_bam, parameters* param, ref_genome* ref);
void init_rd_per_chr( bam_info* in_bam, parameters* param, ref_genome* ref, int chr_index);
void calc_mean(ref_genome*, bam_info*);
void calc_mean_per_chr( ref_genome* ref, bam_info* in_bam, int chr_index);

#endif /* REFGENOME_H_ */
