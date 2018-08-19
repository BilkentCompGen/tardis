/*
 * vh_createMaxClusterNUMT.h
 *
 *  Created on: Jul 31, 2018
 *      Author: tardis
 */

#ifndef VH_VH_CREATEMAXCLUSTERNUMT_H_
#define VH_VH_CREATEMAXCLUSTERNUMT_H_

#include "../processbam.h"

typedef struct numt_Reads{
	char *readName;
	int pos;
	int pos_End;
	char orient;
	int mQual;
	int NUMT_Type; // 0: Forward, 1: Reverse
	char *libName;
	char *indName;
	int libId;
	unsigned long int ten_x_barcode;

	struct numt_Reads *next;
} numt_Reads;

void vh_finalizeReadMapping_NUMT( int chroSize);
void initializeReadMapping_NUMT( bam_info** in_bams, parameters *params, int chr_index);
void NUMTCluster_Region( parameters* params, int chr_index);

#endif /* VH_VH_CREATEMAXCLUSTERNUMT_H_ */
