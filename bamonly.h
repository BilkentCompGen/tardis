/*
 * bo_bamhandler.h
 *
 *  Created on: Aug 23, 2016
 *      Author: tardis
 */

#ifndef BAMONLY_BO_BAMHANDLER_H_
#define BAMONLY_BO_BAMHANDLER_H_

#include "processbam.h"
#include "processfq.h"

#define MIN_INS_LEN 100
#define MAX_INV_LEN 100000000
#define MIN_SOFTCLIP_LEN 10
#define MAX_NUM_CIGAR 3
#define SOFTCLIP_WRONGMAP_WINDOW 20
#define MAX_WIN 10000000
#define MAX_READNAME_LEN 200
#define MIN_MEI_DISTANCE 1000000
#define NOTMEI -1
#define NOTNUMT -1
#define NUMTF 0
#define NUMTR 1

extern long sr_cnt_bam;
extern long alt_cnt_bam;

#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))

int bamonly_run( parameters *params, bam_info ** in_bams);
int read_mapping( library_properties *library, parameters* params, bam1_t* bam_alignment, int32_t *bamToRefIndex, bam_alignment_region* bam_align);

#endif /* BAMONLY_BO_BAMHANDLER_H_ */
