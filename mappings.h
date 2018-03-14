/*
 * altmapping.h
 *
 *  Created on: Dec 25, 2017
 *      Author: tardis
 */

#ifndef MAPPINGS_H_
#define MAPPINGS_H_

#include "processfq.h"

void find_alt_mappings(  bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment,
		char* xa_string, int32_t *bamToRefIndex);
int primRightAltLeftMappings( bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment, int32_t *bamToRefIndex);
int primary_mapping( bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment, int32_t *bamToRefIndex);



#endif /* MAPPINGS_H_ */
