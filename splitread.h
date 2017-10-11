/*
 * splitread.h
 *
 *  Created on: Feb 2, 2017
 *      Author: tardis
 */

#ifndef SPLITREAD_H_
#define SPLITREAD_H_

#include "processbam.h"
#include "bamonly.h"

#define HASHKMERLEN 10

typedef struct lociInRef
{
	int pos; /* pos in the referenceSeqInterest array (keeps the reference genome for each chromosome) */
	struct lociInRef *next;
}lociInRef;

extern char *ref_seq_per_chr;

void addSoftClip( ref_genome* ref, library_properties * library, bam_alignment_region* bam_align, bam1_t* bam_alignment, int chrID);
void create_HashIndex( ref_genome* ref, int chr_index);
void free_HashIndex();
void countNumSoftClipInCluster( parameters *params, ref_genome* ref, bam_info* in_bam, int chr_index);
void mapSoftClipToRef( bam_info* in_bam,  parameters* params, ref_genome* ref, int chr_index);
void readReferenceSeq(ref_genome *ref, parameters *params, int chr_index);

#endif /* SPLITREAD_H_ */
