#ifndef __EXTERNAL
#define __EXTERNAL

#include "common.h"
#include "config.h"
#include "processbam.h"

/* Function Prototypes */
int remap_mrfast( parameters *params, bam_info ** in_bams, configuration *cfg);
void plot_histogram  ( parameters *params, configuration *cfg, char *sample_name, char *libname, int sample_size, int *fragment_lengths, float average, float stdev);


#endif
