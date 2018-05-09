#include "vh_common.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_heap.h"
#include "vh_maximalCluster.h"

void vh_addToGenomeIndex_Insertion (sonic *this_sonic, int chr_index)
{
	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	MappingOnGenome *newEl, *newEl2;

	libInfo = g_libInfo;
	while (libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while (divetReadMappingPtr != NULL)
		{
			if (strcmp (divetReadMappingPtr->chromosome_name, this_sonic->chromosome_names[chr_index]) == 0 && divetReadMappingPtr->svType == INSERTION
			    && !sonic_is_gap (this_sonic, this_sonic->chromosome_names[chr_index], divetReadMappingPtr->locMapLeftEnd, divetReadMappingPtr->locMapRightStart)
					&& (divetReadMappingPtr->locMapRightStart > divetReadMappingPtr->locMapLeftEnd)
					&& divetReadMappingPtr->locMapRightEnd < this_sonic->chromosome_lengths[chr_index] && divetReadMappingPtr->locMapLeftStart > 0)
			{
				newEl = (MappingOnGenome *) getMem (sizeof (MappingOnGenome));
				newEl->readMappingPtr = divetReadMappingPtr;
				newEl->next = g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd];
				g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd] = newEl;
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
	}
}

void vh_initializeReadMapping_Insertion (sonic *this_sonic, int chr_index)
{
	//Gap info is global
	LibraryInfo *libInfoPtr = g_libInfo;
	int genomeIndexId, i;

	//Initing the Genome Array
	g_genomeIndexStart = (MappingOnGenome **) getMem (this_sonic->chromosome_lengths[chr_index] * sizeof (MappingOnGenome *));

	if (g_genomeIndexStart == NULL)
		vh_logWarning ("Memory Problem");
	for (genomeIndexId = 0; genomeIndexId < this_sonic->chromosome_lengths[chr_index]; genomeIndexId++)
	{
		g_genomeIndexStart[genomeIndexId] = NULL;
	}

	//Initing the List of begin and end of right side break point ranges - make sure to free
	g_listRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
	g_tempListRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
	for (i = 0; i < g_maxListBrkPointIntr; i++)
	{
		g_listRightBrkPointIntr[i].readMappingPtr = NULL;
		g_tempListRightBrkPointIntr[i].readMappingPtr = NULL;
	}
	g_listRightBrkPointIntrCount = 0;
	vh_addToGenomeIndex_Insertion (this_sonic, chr_index);

	/////Malocing the intersectingInterval (intersectinterval) heap
	g_intersectInterval = (Heap *) getMem (sizeof (Heap));
	g_intersectInterval->heapSize = 0;

	while (libInfoPtr != NULL)
	{
		g_maxDeltaAmongLibs = vh_max (g_maxDeltaAmongLibs, libInfoPtr->maxDelta);
		libInfoPtr = libInfoPtr->next;
	}
}

int vh_createBreakPointIntervals_Insertion (int brkPointLeft)
{
	MappingOnGenome *ptrMappingOnGenome;
	int newElAdded = 0;
	int genomeId;
	int maxDeltaTemp, minDeltaTemp;
	int locBrkPointLeftTemp, locBrkPointRightTemp;	// These are the min max value of the length of the insertion
	RightBrkPointInterval *temp;
	int listRightBrkPointIntrId = 0, tempListRightBrkPointIntrId = 0;

	if (g_listRightBrkPointIntrCount > 0)
	{
		//TODO: Can this be made more efficient using a Heap?
		while (listRightBrkPointIntrId < g_listRightBrkPointIntrCount)
		{
			if (g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapRightStart == brkPointLeft)
			{
				listRightBrkPointIntrId++;
			}
			else
			{
				vh_copyElBrkPointIntr (tempListRightBrkPointIntrId, listRightBrkPointIntrId);
				listRightBrkPointIntrId++;
				tempListRightBrkPointIntrId++;
			}
		}
	}

	ptrMappingOnGenome = g_genomeIndexStart[brkPointLeft];
	while (ptrMappingOnGenome != NULL)
	{
		newElAdded = 1;
		locBrkPointRightTemp = ptrMappingOnGenome->readMappingPtr->libInfo->maxDelta - (ptrMappingOnGenome->readMappingPtr->locMapRightStart -
				ptrMappingOnGenome->readMappingPtr->locMapLeftEnd);
		locBrkPointLeftTemp = ptrMappingOnGenome->readMappingPtr->libInfo->minDelta - (ptrMappingOnGenome->readMappingPtr->locMapRightStart -
				ptrMappingOnGenome->readMappingPtr->locMapLeftEnd);
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


// This is the head function called for creating clusters. It startes the loop which iterates the left breakpoint
void vh_createInsertionClusters (int chroSize)
{
	int leftBreakPoint;		// The value of the left breakpoint considered
	// The genome starts from 0
	int startLeftWindow, endLeftWindow;
	g_listRightBrkPointIntrCount = 0;
	int newElAdded = 0;
	// The loop which iterates on the chromosome repersenting the left breakpoint
	for (leftBreakPoint = 1; leftBreakPoint < chroSize; leftBreakPoint++)
	{
		newElAdded = 0;
		newElAdded = vh_createBreakPointIntervals_Insertion (leftBreakPoint);
		if (newElAdded)
			vh_createIntersectingIntervals (leftBreakPoint, INSERTION);
	}
	vh_flushOut ( leftBreakPoint, INSERTION);
}
