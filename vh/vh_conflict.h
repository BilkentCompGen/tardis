/*
 * vh_conflict.h
 *
 *  Created on: Apr 4, 2017
 *      Author: tardis
 */

#ifndef VH_VH_CONFLICT_H_
#define VH_VH_CONFLICT_H_

int conflictsAny(int i, int *supInd);
void addToListOfConflicts(int i, int j, int *countReads);
void addToConflict(int maxWeightSet, int *countReads);


#endif /* VH_VH_CONFLICT_H_ */
