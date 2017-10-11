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

void initializeReadMapping_MEI( bam_info** in_bams, parameters *params, char *chromosome_name, int chroSize);
void MEICluster_Region(parameters*, char*, int);
void vh_finalizeReadMapping_Mei(int);

#endif /* VH_VH_CREATEMAXCLUSTERMEI_H_ */
