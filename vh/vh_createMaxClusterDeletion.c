#include "vh_common.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_heap.h"
#include "vh_maximalCluster.h"

int maxDelta;
int minDelta;

void vh_addToGenomeIndex_Deletion( sonic *this_sonic, int chr_index)
{
	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	MappingOnGenome *newEl;
	libInfo = g_libInfo;
	while (libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while (divetReadMappingPtr != NULL)
		{
			if (strcmp ( divetReadMappingPtr->chromosome_name, this_sonic->chromosome_names[chr_index]) == 0
					&& !sonic_is_gap ( this_sonic, this_sonic->chromosome_names[chr_index], divetReadMappingPtr->locMapLeftEnd, divetReadMappingPtr->locMapRightStart)
					&& ( divetReadMappingPtr->locMapRightStart - divetReadMappingPtr->locMapLeftEnd < maxDeletionLen)
					&& divetReadMappingPtr->locMapRightEnd < this_sonic->chromosome_lengths[chr_index] && divetReadMappingPtr->locMapLeftStart > 0)
			{
				if( divetReadMappingPtr->svType == DELETION)
				{
					newEl = (MappingOnGenome *) getMem (sizeof (MappingOnGenome));
					newEl->readMappingPtr = divetReadMappingPtr;
					newEl->next = g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd];
					g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd] = newEl;
				}
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
	}
}

void vh_initializeReadMapping_Deletion( sonic *this_sonic, int chr_index)
{
	//Gap info is global
	LibraryInfo *libInfoPtr = g_libInfo;
	int genomeIndexId;
	int i;

	g_genomeIndexStart = (MappingOnGenome **) getMem (this_sonic->chromosome_lengths[chr_index] * sizeof (MappingOnGenome *));

	if (g_genomeIndexStart == NULL)
		vh_logWarning ("Memory Problem in vh_createMaxClusterDeletion.c:48");

	for (genomeIndexId = 0; genomeIndexId < this_sonic->chromosome_lengths[chr_index]; genomeIndexId++)
		g_genomeIndexStart[genomeIndexId] = NULL;


	//Initing the List of begin and end of right side break point ranges - make sure to free
	g_listRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
	g_tempListRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
	for (i = 0; i < g_maxListBrkPointIntr; i++)
	{
		g_listRightBrkPointIntr[i].readMappingPtr = NULL;
		g_tempListRightBrkPointIntr[i].readMappingPtr = NULL;
	}
	g_listRightBrkPointIntrCount = 0;
	vh_addToGenomeIndex_Deletion( this_sonic, chr_index);


	g_intersectInterval = vh_newHeap(MAX_CLUSTER_SIZE);


	while (libInfoPtr != NULL)
	{
		g_maxDeltaAmongLibs = vh_max (g_maxDeltaAmongLibs, libInfoPtr->maxDelta);
		libInfoPtr = libInfoPtr->next;
	}
}


void vh_reevaluate_Deletion (int id, int brkPointLeft)
{
	int brkRightTemp, brkLeftTemp;
	brkRightTemp = vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
			vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->minDelta - (brkPointLeft -
					g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd), 0),
					g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd);
	brkLeftTemp = vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
			vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->maxDelta - (brkPointLeft -
					g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd), 0),
					g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd);

	if (g_tempListRightBrkPointIntr[id].keyLorR == LEFT)
	{
		g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].key = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;
	}
	else if (g_tempListRightBrkPointIntr[id].keyLorR == RIGHT)
	{
		g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].key = brkRightTemp;
		g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;
	}
}


int vh_createBreakPointIntervals_Deletion (int brkPointLeft)
{
	int newElAdded = 0;
	int maxDeltaTemp, minDeltaTemp;
	int locBrkPointLeftTemp, locBrkPointRightTemp;	// These the left and right loci of the interval of the breakpoint on the right
	int count = 0, tempListRightBrkPointIntrId = 0;
	RightBrkPointInterval *temp;
	MappingOnGenome *ptrMappingOnGenome;

	//qsort (g_listRightBrkPointIntr, g_listRightBrkPointIntrCount, sizeof (struct RightBrkPointInterval), vh_compare);
	if (g_listRightBrkPointIntrCount > 0)
	{
		//TODO: Can this be made more efficient using a Heap?
		while (count < g_listRightBrkPointIntrCount)
		{
			//fprintf(stderr, "%d - %d\n", (g_listRightBrkPointIntr[count].readMappingPtr->locMapLeftEnd + g_listRightBrkPointIntr[count].readMappingPtr->libInfo->maxDelta), brkPointLeft);
			if (g_listRightBrkPointIntr[count].readMappingPtr->locMapLeftEnd +
					g_listRightBrkPointIntr[count].readMappingPtr->libInfo->maxDelta == brkPointLeft)
			{
				count++;
			}
			else
			{
				vh_copyElBrkPointIntr( tempListRightBrkPointIntrId, count); // increaseByOneRightBrkPointIntr(tempListRightBrkPointIntrId);
				vh_reevaluate_Deletion( tempListRightBrkPointIntrId, brkPointLeft);

				count++;
				tempListRightBrkPointIntrId++;
			}
		}
		//fprintf(stderr, "\n");
	}

	ptrMappingOnGenome = g_genomeIndexStart[brkPointLeft];
	while( ptrMappingOnGenome != NULL)
	{
		newElAdded = 1;
		locBrkPointRightTemp = ptrMappingOnGenome->readMappingPtr->locMapRightStart -
				vh_max ((ptrMappingOnGenome->readMappingPtr->libInfo->minDelta - (brkPointLeft - ptrMappingOnGenome->readMappingPtr->locMapLeftEnd)), 0);
		locBrkPointLeftTemp = ptrMappingOnGenome->readMappingPtr->locMapRightStart -
				vh_max ((ptrMappingOnGenome->readMappingPtr->libInfo->maxDelta - (brkPointLeft - ptrMappingOnGenome->readMappingPtr->locMapLeftEnd)), 0);

		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointRight = locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft = locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key = locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr = ptrMappingOnGenome->readMappingPtr;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = LEFT;
		tempListRightBrkPointIntrId++;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointRight = locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft = locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key = locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr = ptrMappingOnGenome->readMappingPtr;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = RIGHT;
		tempListRightBrkPointIntrId++;
		ptrMappingOnGenome = ptrMappingOnGenome->next;
	}
	temp = g_listRightBrkPointIntr;
	g_listRightBrkPointIntr = g_tempListRightBrkPointIntr;
	g_tempListRightBrkPointIntr = temp;
	g_listRightBrkPointIntrCount = tempListRightBrkPointIntrId;

	return (newElAdded);
}


// This is the head function called for creating clusters. It starts the loop which iterates the left breakpoint
void vh_createDeletionClusters (int chroSize)
{
	int leftBreakPoint;		// The value of the left breakpoint considered
	g_listRightBrkPointIntrCount = 0;
	int newElAdded = 0;

	// The loop which iterates on the chromosome representing the left breakpoint
	for (leftBreakPoint = 1; leftBreakPoint < chroSize; leftBreakPoint++)
	{
		newElAdded = 0;
		newElAdded = vh_createBreakPointIntervals_Deletion( leftBreakPoint);
		if( newElAdded)		// For deletion only when we have added new element we need to check
			vh_createIntersectingIntervals( leftBreakPoint, DELETION);
	}
	vh_flushOut ( leftBreakPoint, DELETION);
}
