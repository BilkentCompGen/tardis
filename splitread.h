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
#define MAX_SR_HIT 50000

typedef struct lociInRef
{
	int pos; /* pos in the referenceSeqInterest array (keeps the reference genome for each chromosome) */
	struct lociInRef *next;
}lociInRef;

extern char *ref_seq_per_chr;

void addSoftClip( library_properties * library, bam_alignment_region* bam_align, bam1_t* bam_alignment, char* chromosome_nameD);
void init_hash_count(void);
void build_hash_table( const char *ref, int len, int mode);
void create_hash_table( char *ref, int len);
void free_hash_table( void);
void init_hash_table( void);
void create_HashIndex( parameters* params, int chr_index);
void free_HashIndex(void);
void countNumSoftClipInCluster( parameters *params, bam_info* in_bam, int chr_index);
void mapSoftClipToRef( bam_info* in_bam,  parameters* params, int chr_index);
void readReferenceSeq( parameters *params, int chr_index);
unsigned int hash_function_next( unsigned int prev_hash, unsigned int mask, const char next_char);
unsigned int hash_function_ref( char *str);
int is_kmer_valid ( char *str);
  
#endif /* SPLITREAD_H_ */
