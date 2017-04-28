#include "vh_common.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_heap.h"
#include "vh_maximalCluster.h"

//#define NAME_STR_LEN 100



void vh_addToGenomeIndex_Inversion (char *chroName)
{
	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	MappingOnGenome *newEl, *newEl2;
	int leftWindowEnd, leftWindowStart;
	libInfo = g_libInfo;
	while (libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while (divetReadMappingPtr != NULL)
		{
			if (strcmp (divetReadMappingPtr->chroName, chroName) == 0 && divetReadMappingPtr->svType == 'V'
					&& (divetReadMappingPtr->locMapRightStart - divetReadMappingPtr->locMapLeftEnd < maxInversionLen))
			{
				newEl = (MappingOnGenome *) getMem (sizeof (MappingOnGenome));
				newEl2 = (MappingOnGenome *) getMem (sizeof (MappingOnGenome));
				newEl->readMappingPtr = divetReadMappingPtr;
				newEl2->readMappingPtr = divetReadMappingPtr;
				if (divetReadMappingPtr->orientationLeft == 'F' && divetReadMappingPtr->orientationRight == 'F')
				{
					newEl->next = g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd];
					leftWindowEnd = divetReadMappingPtr->locMapLeftEnd + libInfo->maxDelta;
					newEl2->next = g_genomeIndexEnd[leftWindowEnd];
					g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd] = newEl;
					g_genomeIndexEnd[leftWindowEnd] = newEl2;
				}
				else if (divetReadMappingPtr->orientationLeft == 'R' && divetReadMappingPtr->orientationRight == 'R')
				{
					leftWindowStart = vh_max (divetReadMappingPtr->locMapLeftStart - libInfo->maxDelta, 0);
					newEl->next = g_genomeIndexStart[leftWindowStart];
					newEl2->next = g_genomeIndexEnd[divetReadMappingPtr->locMapLeftStart];
					g_genomeIndexStart[leftWindowStart] = newEl;
					g_genomeIndexEnd[divetReadMappingPtr->locMapLeftStart] = newEl2;
				}
				/*
		newEl->next=g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd];
		leftWindowEnd=divetReadMappingPtr->locMapLeftEnd+libInfo->maxDelta;
		newEl2->next=g_genomeIndexEnd[leftWindowEnd];
		g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd]=newEl;
		g_genomeIndexEnd[leftWindowEnd]=newEl2;*/
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
	}

}

void vh_initializeReadMapping_Inversion (char *chroName, int chroSize)
{
	//Gap info is global
	LibraryInfo *libInfoPtr = g_libInfo;
	int genomeIndexId;
	int i;

	//Initing the Genome Array
	g_genomeIndexStart = (MappingOnGenome **) getMem (chroSize * sizeof (MappingOnGenome *));
	g_genomeIndexEnd = (MappingOnGenome **) getMem (chroSize * sizeof (MappingOnGenome *));

	if (g_genomeIndexStart == NULL || g_genomeIndexEnd == NULL)
		vh_logWarning ("Memory Problem");
	for (genomeIndexId = 0; genomeIndexId < chroSize; genomeIndexId++)
	{
		g_genomeIndexStart[genomeIndexId] = NULL;
		g_genomeIndexEnd[genomeIndexId] = NULL;
	}
	//Initing the List of begin and end of right side break point ranges
	g_listRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
	g_tempListRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
	for (i = 0; i < g_maxListBrkPointIntr; i++)
	{
		g_listRightBrkPointIntr[i].readMappingPtr = NULL;
		g_tempListRightBrkPointIntr[i].readMappingPtr = NULL;
	}
	g_listRightBrkPointIntrCount = 0;
	vh_addToGenomeIndex_Inversion (chroName);

	/////Malocing the intersectingInterval (intersectinterval) heap
	g_intersectInterval = (Heap *) getMem (sizeof (Heap));
	g_intersectInterval->heapSize = 0;

	while (libInfoPtr != NULL)
	{
		g_maxDeltaAmongLibs = vh_max (g_maxDeltaAmongLibs, libInfoPtr->maxDelta);
		libInfoPtr = libInfoPtr->next;
	}
}

///////////////////////////////////The left end side and right end side of the Right breakpoint//////////////
void vh_reevaluate_rightBrkPoint_FF (int id, int brkPointLeft)
{
	int brkRightTemp, brkLeftTemp;
	// we need to take the max with locMapRightStart, because it might be possible that (brkPointLeft-g_tempListRightPointIntr[id].readMappingPtr->locMapLeftEnd) be bigger than minDelta.
	brkRightTemp =
			vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd +
					g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->maxDelta -
					(brkPointLeft -
							g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd),
							g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd);

	brkLeftTemp =
			vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd +
					g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->minDelta -
					(brkPointLeft -
							g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd),
							g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd);

	if (g_tempListRightBrkPointIntr[id].keyLorR == 'L')
	{
		g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].key = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;
	}
	else if (g_tempListRightBrkPointIntr[id].keyLorR == 'R')
	{
		g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].key = brkRightTemp;
		g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;

	}


}

void vh_reevaluate_rightBrkPoint_RR (int id, int brkPointLeft)
{
	int brkRightTemp, brkLeftTemp;
	brkRightTemp =
			vh_min (vh_max	 (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -(g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->minDelta -
					(g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd -  brkPointLeft)),
					g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftStart),
					g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart);

	brkLeftTemp =
			vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
					(g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->maxDelta -
							(g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftStart -
									brkPointLeft)),
									g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftStart);

	if (g_tempListRightBrkPointIntr[id].keyLorR == 'L')
	{
		g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].key = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;
	}
	else if (g_tempListRightBrkPointIntr[id].keyLorR == 'R')
	{
		g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
		g_tempListRightBrkPointIntr[id].key = brkRightTemp;
		g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;

	}


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int vh_createBreakPointIntervals_Inversion (int brkPointLeft)
{
	MappingOnGenome *ptrMappingOnGenome;
	int newElAdded = 0;
	int genomeId;
	int maxDeltaTemp, minDeltaTemp;
	int locBrkPointLeftTemp, locBrkPointRightTemp;	// These the left and right loci of the interval of the breakpoint on the right
	RightBrkPointInterval *temp;
	int listRightBrkPointIntrId = 0, tempListRightBrkPointIntrId = 0;

	if (g_listRightBrkPointIntrCount > 0)
	{
		ptrMappingOnGenome = g_genomeIndexEnd[brkPointLeft];
		//TODO: Can this be made more efficient using a Heap?
		while (listRightBrkPointIntrId < g_listRightBrkPointIntrCount)
		{
			if (g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationLeft == 'F'
					&& g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationRight =='F')
			{
				if ((g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapLeftEnd +
						g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->libInfo->maxDelta == brkPointLeft)
						|| (brkPointLeft == g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapRightEnd))
				{
					listRightBrkPointIntrId++;
				}
				else
				{
					vh_copyElBrkPointIntr (tempListRightBrkPointIntrId, listRightBrkPointIntrId);
					//increaseByOneRightBrkPointIntr(tempListRightBrkPointIntrId);
					//decreaseByOneRightBrkPointIntr_FF(tempListRightBrkPointIntrId);
					vh_reevaluate_rightBrkPoint_FF (tempListRightBrkPointIntrId, brkPointLeft);
					listRightBrkPointIntrId++;
					tempListRightBrkPointIntrId++;
				}
			}
			else
				if (g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationLeft =='R'
						&& g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationRight =='R')
				{
					if ((g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapLeftStart == brkPointLeft)
							||(g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapRightStart - brkPointLeft <
									g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->libInfo->minDelta))
					{
						listRightBrkPointIntrId++;
					}
					else
					{
						vh_copyElBrkPointIntr (tempListRightBrkPointIntrId,listRightBrkPointIntrId);
						//decreaseByOneRightBrkPointIntr_RR(tempListRightBrkPointIntrId);
						vh_reevaluate_rightBrkPoint_RR (tempListRightBrkPointIntrId, brkPointLeft);
						listRightBrkPointIntrId++;
						tempListRightBrkPointIntrId++;
					}
				}
		}
	}


	ptrMappingOnGenome = g_genomeIndexStart[brkPointLeft];
	while (ptrMappingOnGenome != NULL)
	{
		newElAdded = 1;

		if (ptrMappingOnGenome->readMappingPtr->orientationLeft == 'F' && ptrMappingOnGenome->readMappingPtr->orientationRight == 'F')
		{
			locBrkPointRightTemp = vh_max (ptrMappingOnGenome->readMappingPtr->libInfo->maxDelta -
					(brkPointLeft - ptrMappingOnGenome->readMappingPtr->locMapLeftEnd), 0) +
							ptrMappingOnGenome->readMappingPtr->locMapRightEnd;
			locBrkPointLeftTemp = vh_max (ptrMappingOnGenome->readMappingPtr->libInfo->minDelta -
					(brkPointLeft - ptrMappingOnGenome->readMappingPtr->locMapLeftEnd), 0) +
							ptrMappingOnGenome->readMappingPtr->locMapRightEnd;
		}
		else if (ptrMappingOnGenome->readMappingPtr->orientationLeft == 'R' && ptrMappingOnGenome->readMappingPtr->orientationRight == 'R')
		{
			locBrkPointRightTemp = vh_max (ptrMappingOnGenome->readMappingPtr->locMapRightStart -
					vh_max (ptrMappingOnGenome->readMappingPtr->libInfo->minDelta -
							(ptrMappingOnGenome->readMappingPtr->locMapLeftStart -
									brkPointLeft), 0), ptrMappingOnGenome->readMappingPtr->locMapLeftStart);
			locBrkPointLeftTemp = vh_max (ptrMappingOnGenome->readMappingPtr->locMapRightStart -
					vh_max (ptrMappingOnGenome->readMappingPtr->libInfo->maxDelta -
							(ptrMappingOnGenome->readMappingPtr->locMapLeftStart - brkPointLeft), 0), ptrMappingOnGenome->readMappingPtr->locMapLeftStart);
		}

		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointRight = locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft = locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key = locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = 'L';
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr = ptrMappingOnGenome->readMappingPtr;
		tempListRightBrkPointIntrId++;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointRight = locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft = locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key = locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = 'R';
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr = ptrMappingOnGenome->readMappingPtr;
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
void vh_createInversionClusters (int chroSize)
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
		newElAdded = vh_createBreakPointIntervals_Inversion (leftBreakPoint);
		//if (newElAdded) For Inversions we need to check for every new BrkPoint
		vh_createIntersectingIntervals (leftBreakPoint, 3);
	}
	vh_flushOut (g_listPotClusterFound, leftBreakPoint, 3);
}
