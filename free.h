/*
 * free.h
 *
 *  Created on: Sep 26, 2017
 *      Author: tardis
 */

#ifndef FREE_H_
#define FREE_H_

#include "bamonly.h"

void free_alignments( bam_alignment_region** bam_align);
void free_alignments2( bam_alignment_region** bam_align);
void free_mappings( bam_info** in_bams, parameters* params);
void free_the_rest( bam_info** in_bams, parameters* params);
void free_sensitive(bam_info** in_bams, parameters *params);
void free_quick(bam_info** in_bams, parameters *params);
void free_libraries();

#endif /* FREE_H_ */
