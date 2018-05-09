/*
 * altmapping.h
 *
 *  Created on: Dec 25, 2017
 *      Author: tardis
 */

#ifndef MAPPINGS_H_
#define MAPPINGS_H_

#include "processfq.h"

extern int altLeftPrimRight, primRightAltLeft, altLeftAltRight, primLeftAltRight;

void find_alt_mappings(  bam_info* in_bam, parameters* params, int lib_index, bam1_t* bam_alignment,
		char* xa_string);
int primRightAltLeftMappings( library_properties* library, bam1_t* bam_alignment);
int altLeftAltRightMappings( library_properties *library, bam1_t* bam_alignment, char altmap[4][1024]);
int primLeftAltRightMappings( library_properties *library, bam1_t* bam_alignment, char altmap[4][1024]);
int primary_mapping( bam_info* in_bam, parameters* params, int lib_index, bam1_t* bam_alignment, int32_t *bamToRefIndex);



#endif /* MAPPINGS_H_ */
