#ifndef MAIN_H__
#define MAIN_H__

#include "../refgenome.h"
#include "vh_common.h"


void vh_clustering (bam_info** in_bams, ref_genome* ref, parameters *params, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit);
void vh_clustering_bam(bam_info**, ref_genome*, parameters*, double, char*, char*, int);
void vh_quitProgram (int);
int run_vh( ref_genome* ref, parameters *params, bam_info ** in_bams);

#endif
