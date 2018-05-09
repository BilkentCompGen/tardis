/*
 * splitread.h
 *
 *  Created on: Feb 2, 2017
 *      Author: tardis
 */

#ifndef SPLITREAD_H_
#define SPLITREAD_H_

#include "processbam.h"
#include "processfq.h"
#include "bamonly.h"

#define HASHKMERLEN 10

typedef struct lociInRef
{
	int pos; /* pos in the referenceSeqInterest array (keeps the reference genome for each chromosome) */
	struct lociInRef *next;
}lociInRef;

extern char *ref_seq_per_chr;

void addSoftClip( library_properties * library, bam_alignment_region* bam_align, bam1_t* bam_alignment, char* chromosome_nameD);
void create_HashIndex( parameters* params, int chr_index);
void free_HashIndex();
void countNumSoftClipInCluster( parameters *params, bam_info* in_bam, int chr_index);
void mapSoftClipToRef( bam_info* in_bam,  parameters* params, int chr_index);
void readReferenceSeq( parameters *params, int chr_index);

#endif /* SPLITREAD_H_ */
