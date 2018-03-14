/*
 * vh_createMaxClusterMEI.h
 *
 *  Created on: Sep 20, 2016
 *      Author: tardis
 */

#ifndef VH_VH_CREATEMAXCLUSTERMEI_H_
#define VH_VH_CREATEMAXCLUSTERMEI_H_

#include "../processbam.h"

#define lenSplitReadBrakWindow 20

typedef struct mei_Reads{
	char *readName;
	int pos;
	int pos_End;
	char orient;
	int mQual;
	int MEI_Type; // 0 : Alu +; 1: Alu -; 2: L1 +; 3; L1 -; 4: SVA +; 5; SVA -
	char* MEI_Class;
	char* MEI_Subclass;
	char *libName;
	char *indName;
	int libId;
	int readTypeSupport;// 0: ReadPair, 1: SplitRead (Soft clipping of the beginig of the read, 2: SplitRead (Soft clipping at the end of the read)

	struct mei_Reads *next;
} mei_Reads;

void initializeReadMapping_MEI( bam_info** in_bams, parameters *params, char *chromosome_name, int chroSize);
void MEICluster_Region( parameters* params, char* chromosome_name, int chroSize);
void vh_finalizeReadMapping_Mei( int chroSize);


#endif /* VH_VH_CREATEMAXCLUSTERMEI_H_ */
