#ifndef MAIN_H__
#define MAIN_H__

#include "../processrefgen.h"
#include "vh_common.h"


void vh_clustering (bam_info** in_bams, parameters *params, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit);
void vh_quitProgram (int);
int run_vh( parameters *params, bam_info ** in_bams);

#endif
