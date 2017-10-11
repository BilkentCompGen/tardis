/*
 * vh_conflict.c
 *
 *  Created on: Apr 4, 2017
 *      Author: tardis
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vh_conflict.h"
#include "vh_setcover.h"

int numSV = 0; // total number of distinct SVs picked

void removeInd( int ind, int j)
{
	listClusterEl[j].indIdCount[ind] = -1;
	listClusterEl[j].oldBestIsGood = 0;
}

int conflictsBetweenTwoSV_Cord( int startCoord1, int stopCoord1,char SVtype1, int startCoord2, int stopCoord2, char SVtype2)
{
	if( ( startCoord1 > stopCoord2 || stopCoord1 < startCoord2))
		return 0;
	if( SVtype1 == MEIFORWARD || SVtype1 == MEIREVERSE || SVtype2 == MEIFORWARD || SVtype2 == MEIREVERSE)
		return 0;
	if( SVtype1 == INVERSION && SVtype2 == INVERSION)
	{
		if( ( ( startCoord1 > startCoord2) && ( stopCoord1 < stopCoord2)) || ( ( startCoord1 < startCoord2) && ( stopCoord1 > stopCoord2)))
			return 0;
	}
	else if( SVtype1 == INVERSION && SVtype2 == DELETION)
	{
		if( startCoord1 < startCoord2 && stopCoord1 > stopCoord2)
			return 0;
	}
	else if( SVtype1 == DELETION && SVtype2 == INVERSION)
	{
		if( startCoord2 < startCoord1 && stopCoord2 > stopCoord1)
			return 0;
	}
	return 1;
}

int conflictsAny(int i, int *supInd) // return the individual that SV_i is in conflict with any of the selected SVs
{
	SV_selected *ptrSV;
	int conflict=0;
	int count;

	for (count=0; count<numSV; count++)
	{
		if (strcmp(listSelectedSV[count].chromosome_name, listClusterEl[i].chromosome_name)==0)
		{
			if(conflictsBetweenTwoSV_Cord(listSelectedSV[count].posStart_SV, listSelectedSV[count].posEnd_SV, listSelectedSV[count].SVtype, listClusterEl[i].posStartSV, listClusterEl[i].posEndSV, listClusterEl[i].SVtype))
			{
				ptrSV=listSelectedSV[count].conflict_Next;
				while (ptrSV!=NULL)
				{
					if (conflictsBetweenTwoSV_Cord(ptrSV->posStart_SV, ptrSV->posEnd_SV, ptrSV->SVtype, listClusterEl[i].posStartSV, listClusterEl[i].posEndSV, listClusterEl[i].SVtype))
					{
						int countInd=0;
						for (countInd=0; countInd<multiIndCount; countInd++)
						{
							if (ptrSV->sup[countInd] > 0 && listClusterEl[i].indIdCount[countInd] <= 0 && listSelectedSV[count].sup[countInd] > 0 && supInd[countInd]>0)
							{
								removeInd(countInd, i);
								conflict=1;
							}
						}
					}
					ptrSV=ptrSV->conflict_Next;
				}
			}
		}
	}
	if (conflict==1) return 1;
	else return -1;
}

int checkConflictNewSelected(SV_selected *ptr1, SV_selected *ptr2, int NotSelectedClus)
{
	int stillCanBeSelected=0;
	int count=0;

	for (count=0; count<multiIndCount; count++)
	{
		if (listClusterEl[NotSelectedClus].indIdCount[count]==0)
			stillCanBeSelected=1;
	}
	if (stillCanBeSelected==0)
		return 0;

	if (conflictsBetweenTwoSV_Cord(ptr2->posStart_SV, ptr2->posEnd_SV, ptr2->SVtype,  listClusterEl[NotSelectedClus].posStartSV, listClusterEl[NotSelectedClus].posEndSV, listClusterEl[NotSelectedClus].SVtype))
	{
		int countInd;
		for (countInd=0; countInd<multiIndCount; countInd++)
		{
			if (listClusterEl[NotSelectedClus].indIdCount[countInd]==0 && ptr1->sup[countInd]>0 && ptr2->sup[countInd]>0)
			{
				listClusterEl[NotSelectedClus].indIdCount[countInd]=-1;
				listClusterEl[NotSelectedClus].oldBestIsGood=0;
			}
		}
	}
	return 0;
}

int conflictsTwoWay(int i, int j) //returns the individual that SV_i and SV_j are in conflict among the slected SVs
{
	int count=0;

	if (strcmp(listSelectedSV[i].chromosome_name, listSelectedSV[j].chromosome_name)!=0)
		return -1;
	if (conflictsBetweenTwoSV_Cord(listSelectedSV[i].posStart_SV, listSelectedSV[i].posEnd_SV, listSelectedSV[i].SVtype, listSelectedSV[j].posStart_SV, listSelectedSV[j].posEnd_SV, listSelectedSV[j].SVtype)==0)
		return -1;

	for (count=0; count<multiIndCount; count++)
	{
		if (listSelectedSV[i].sup[count]>0 && listSelectedSV[j].sup[count]>0)
			return count;
	}
	return -1;
}

int wasNotInConliftNowIsConflict(int i, int j, int *countReads)
// add do the two SV selected i and j not in conflict previously
//however with the additional supported given to cluster "j" they will be in conflict
{
	int count;

	if (strcmp(listSelectedSV[i].chromosome_name, listSelectedSV[j].chromosome_name)!=0)
		return -1;
	if (conflictsBetweenTwoSV_Cord(listSelectedSV[i].posStart_SV, listSelectedSV[i].posEnd_SV, listSelectedSV[i].SVtype, listSelectedSV[j].posStart_SV, listSelectedSV[j].posEnd_SV, listSelectedSV[j].SVtype)==0)
		return -1;

	for (count=0; count<multiIndCount; count++)
	{
		if (listSelectedSV[i].sup[count]>0 && listSelectedSV[j].sup[count]>0)
			return -1;
	}
	for (count=0; count<multiIndCount; count++)
	{
		if (listSelectedSV[i].sup[count]>0 && listSelectedSV[j].sup[count]==0 && countReads[count]>0)
			return count;
	}
	return -1;
}

void addToListOfConflicts(int i, int j, int *countReads) // adds SV i to SV j's conflict list and SV j to SV i's conflict list
{
	int count;
	SV_selected *newSV;
	newSV = (SV_selected *) getMem( sizeof( SV_selected));
	newSV->posStart_SV = listSelectedSV[i].posStart_SV;
	newSV->posEnd_SV = listSelectedSV[i].posEnd_SV;
	newSV->clusterId = listSelectedSV[i].clusterId;
	newSV->SVtype = listSelectedSV[i].SVtype;

	newSV->chromosome_name = NULL;
	set_str( &(newSV->chromosome_name), listSelectedSV[i].chromosome_name);

	for( count = 0; count < multiIndCount; count++)
	{
		newSV->sup[count] = listSelectedSV[i].sup[count];
	}
	newSV->conflict_Next = listSelectedSV[j].conflict_Next;
	listSelectedSV[j].conflict_Next=newSV;
	newSV=listSelectedSV[j].conflict_Next;
}


void addToConflict(int maxWeightSet, int *countReads)// adds the SV numSV to the conflicts tables of previous SVs
{
	int i;
	int idSV2Add = -1;
	SV_selected *ptrTmp;

	for( i = 0; i < numSV; i++)
	{
		if( listSelectedSV[i].clusterId == maxWeightSet)
			idSV2Add = i;
	}
	if( idSV2Add == -1)
	{
		listSelectedSV[numSV].chromosome_name = NULL;
		set_str( &(listSelectedSV[numSV].chromosome_name), listClusterEl[maxWeightSet].chromosome_name);

		listSelectedSV[numSV].posStart_SV = listClusterEl[maxWeightSet].posStartSV;
		listSelectedSV[numSV].clusterId = maxWeightSet;
		listSelectedSV[numSV].posEnd_SV = listClusterEl[maxWeightSet].posEndSV;
		listSelectedSV[numSV].conflict_Next = NULL;
		listSelectedSV[numSV].SVtype = listClusterEl[maxWeightSet].SVtype;

		for( i = 0; i < multiIndCount; i++)
			listSelectedSV[numSV].sup[i] = listClusterEl[maxWeightSet].indIdCount[i];

		for( i = 0; i < numSV; i++)
		{
			if( conflictsTwoWay( numSV, i) >= 0)
			{
				addToListOfConflicts( i, numSV, countReads);
				addToListOfConflicts( numSV, i, countReads);
			}
		}

		for( i = 0; i < sizeListClusterEl; i++)
		{
			ptrTmp = listSelectedSV[numSV].conflict_Next;
			if( conflictsBetweenTwoSV_Cord( listSelectedSV[numSV].posStart_SV, listSelectedSV[numSV].posEnd_SV,
					listSelectedSV[numSV].SVtype, listClusterEl[i].posStartSV, listClusterEl[i].posEndSV,
					listClusterEl[i].SVtype) && (strcmp(listSelectedSV[numSV].chromosome_name, listClusterEl[i].chromosome_name) == 0))
			{
				while( ptrTmp != NULL)
				{
					checkConflictNewSelected( &listSelectedSV[numSV], ptrTmp, i);
					ptrTmp = ptrTmp->conflict_Next;
				}
			}
		}
		numSV++;
	}
	else
	{
		for( i = 0; i < numSV; i++)
		{
			if( i != idSV2Add && wasNotInConliftNowIsConflict( i, idSV2Add, countReads) > 0)
			{
				addToListOfConflicts( i, idSV2Add, countReads);
				addToListOfConflicts( idSV2Add, i, countReads);
			}
		}
		for( i = 0; i < multiIndCount; i++)
			listSelectedSV[idSV2Add].sup[i] = listSelectedSV[idSV2Add].sup[i] + countReads[i];

		for( i = 0; i < sizeListClusterEl; i++)
		{
			ptrTmp = listSelectedSV[idSV2Add].conflict_Next;
			if( conflictsBetweenTwoSV_Cord( listSelectedSV[idSV2Add].posStart_SV, listSelectedSV[idSV2Add].posEnd_SV,
					listSelectedSV[idSV2Add].SVtype, listClusterEl[i].posStartSV, listClusterEl[i].posEndSV,
					listClusterEl[i].SVtype) && (strcmp(listSelectedSV[idSV2Add].chromosome_name, listClusterEl[i].chromosome_name)==0))
			{
				while( ptrTmp != NULL)
				{
					checkConflictNewSelected( &listSelectedSV[numSV], ptrTmp, i);
					ptrTmp = ptrTmp->conflict_Next;
				}
			}
		}
	}
}


