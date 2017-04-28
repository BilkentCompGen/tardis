/*
 * vh_createMaxClusterMEI.h
 *
 *  Created on: Sep 20, 2016
 *      Author: tardis
 */

#ifndef VH_VH_CREATEMAXCLUSTERMEI_H_
#define VH_VH_CREATEMAXCLUSTERMEI_H_

#include "vh_maximalCluster.h"
#define lenSplitReadBrakWindow 20

typedef struct heapEl2{
	int priorityValue;
	int startPos;
	int stopPos;
	char cloneName[200];
	int editDist;
	int clonHitId;
	float seqSimProb;
	char orient1;
	char orient2;
	int indexOfMobileEl; // The Index of Mobile El
	int breakPointStart;//The starting position of breakpoint that the paired-end alignment will be valid cluster
	int breakPointEnd;//The end position of breakpoint that the paired-end alignment will not be valid cluster
	int libId;
	int mQual1;
	int mQual2;
	char libName[200];
	char indName[200];
	char* mei_subclass;
}HeapEl2;

typedef struct Heap2{
	HeapEl2 heapArray[MAX_CLUSTER_SIZE];
	int heapSize;
} Heap2;

typedef struct mei_Reads{
	char *readName;
	int pos;
	int pos_End;
	char orient;
	int mQual;
	int MEI_Type; // 0 : Alu +; 1: Alu -; 2: L1 +; 3; L1 -; 4: SVA +; 5; SVA -
	char* MEI_Subclass;
	char *libName;
	char *indName;
	int libId;
	int readTypeSupport;// 0: ReadPair, 1: SplitRead (Soft clipping of the beginig of the read, 2: SplitRead (Soft clipping at the end of the read)
	struct mei_Reads *next;
} mei_Reads;

void vh_initializeReadMapping_MEI (bam_info**, parameters *, char *, int);
void initializeReadMapping_MEI (bam_info** in_bams, ref_genome* ref, parameters *params, char *chroName, int chroSize);
void MEICluster_Region(bam_info**, char*, int);
void vh_finalizeReadMapping_Mei(int);

#endif /* VH_VH_CREATEMAXCLUSTERMEI_H_ */
