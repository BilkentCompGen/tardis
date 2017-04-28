/*
 * splitread.h
 *
 *  Created on: Feb 2, 2017
 *      Author: tardis
 */

#ifndef SPLITREAD_H_
#define SPLITREAD_H_

#include "processbam.h"

#define SR_HASH_SIZE 1000000

typedef struct lociInRef
{
	int pos; /* pos in the referenceSeqInterest array (keeps the reference genome for each chromosome) */
	struct lociInRef *next;
}lociInRef;

extern char *ref_seq_per_chr;

void addSoftClip( ref_genome* ref, bam_info * in_bam, bam1_t* bam_alignment, int library_index, int flag, int *op, int *opl, int chrID);
void create_10bp_HashIndex( ref_genome* ref, int chr_index);
void free_10bp_HashIndex();
void countNumSoftClipInCluster( parameters *params, ref_genome* ref, bam_info* in_bam, int chr_index);
void mapSoftClipToRef( bam_info* in_bam,  parameters* params, ref_genome* ref, int chr_index);
void readReferenceSeq(ref_genome *ref, parameters *params, int chr_index);

#endif /* SPLITREAD_H_ */
