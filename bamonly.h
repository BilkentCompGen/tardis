/*
 * bo_bamhandler.h
 *
 *  Created on: Aug 23, 2016
 *      Author: tardis
 */

#ifndef BAMONLY_BO_BAMHANDLER_H_
#define BAMONLY_BO_BAMHANDLER_H_

#include "processbam.h"

#define MIN_INS_LEN 100
#define MAX_INV_LEN 100000000
#define MIN_SOFTCLIP_LEN 10
#define SOFTCLIP_WRONGMAP_WINDOW 20
#define MAX_WIN 10000000
#define MAX_READNAME_LEN 200
#define MIN_MEI_DISTANCE 1000000
#define NOTMEI -1

#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))

extern long del_cnt_div;
extern long ins_cnt_div;
extern long inv_cnt_div;
extern long mei_cnt_div;
extern long tandup_cnt_div;
extern long sr_cnt_div;
extern long sr_cnt_bam;

int bamonly_run(ref_genome* ref, parameters *params, bam_info ** in_bams);

#endif /* BAMONLY_BO_BAMHANDLER_H_ */
