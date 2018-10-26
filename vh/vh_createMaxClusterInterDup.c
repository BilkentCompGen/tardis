/*
 * vh_createMaxClusterInterDup.c
 *
 *  Created on: Mar 13, 2018
 *      Author: tardis
 */

#include "vh_common.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_heap.h"
#include "vh_maximalCluster.h"
#include "vh_setcover.h"

#define BRKPOINTCLUSTERS 100000
MappingOnGenome **g_genomeIndexStartFR;
MappingOnGenome **g_genomeIndexStartRF;

PEAlistEls **leftSidePEAPicked;//This is an array the size of chrN which keep all the cluster created for a particular breakpoint when the paired-end read is FR
int sizeOfLeftSidePEAPicked;

int *posOfBrkPointOfFR_Cluster;//[100000];
int countPosBrkPointFR_Clusters = 0;
int size_of_brkpoint_clusters;

void vh_finalizeReadMapping_InterDup (char *chromosome_name, int chroSize)
{
	int i;

	ClustersFound *tempPtrOldCluster;
	ClustersFound *ptrToOldClusterList = g_listPotClusterFound;
	while (ptrToOldClusterList != NULL)
	{
		tempPtrOldCluster = ptrToOldClusterList->next;
		free (ptrToOldClusterList->readMappingPtrArray);
		free (ptrToOldClusterList->readMappingIdArray);
		free (ptrToOldClusterList);
		ptrToOldClusterList = tempPtrOldCluster;
	}
	g_listPotClusterFound = NULL;

	for (i = 0; i < max_chromosome_size; i++)
	{
		vh_freeLinkedList (g_genomeIndexStartFR[i]);
		g_genomeIndexStartFR[i] = NULL;
	}
	for (i = 0; i < max_chromosome_size + 1000; i++)
	{
		if (leftSidePEAPicked[i] != NULL)
			free(leftSidePEAPicked[i]);
	}
	free(leftSidePEAPicked);

	free( g_genomeIndexStartFR);
	g_genomeIndexStartFR = NULL;

	for (i = 0; i < max_chromosome_size; i++)
	{
		vh_freeLinkedList (g_genomeIndexStartRF[i]);
		g_genomeIndexStartRF[i] = NULL;
	}
	free( g_genomeIndexStartRF);
	g_genomeIndexStartRF = NULL;

	//Free g_listRightBrkPointIntrCount and set the g_listRightBrkPointIntrCount
	free (g_listRightBrkPointIntr);
	g_listRightBrkPointIntr = NULL;
	g_listRightBrkPointIntrCount = 0;

	//Free g_intersectInterval -> Heap AND set the heapsize
	vh_free_heap (g_intersectInterval);
	g_intersectInterval = NULL;

	freeMem( posOfBrkPointOfFR_Cluster, (size_of_brkpoint_clusters * sizeof(int)));
}


int addToPotentialOutputInterDup( PEAlistEls *leftList, Heap* rightList, int leftBreakPoint, char SVtype)
{
	int count, leftListSize = 0;
	ClustersFound *newCluster;
	ClustersFound *ptrToOldClusterList = NULL;
	ClustersFound *tempPtrOldCluster;
	PEAlistEls *tempPtr;
	int newClusterIsMaximal = 1;
	int oldClusterIsMaximal = 1;

	tempPtr = leftList;
	while(tempPtr != NULL)
	{
		leftListSize++;
		tempPtr = tempPtr->next;
	}

	newCluster = (ClustersFound *) getMem (sizeof (ClustersFound));
	newCluster->next = NULL;
	newCluster->isMaximalCluster = 1;
	newCluster->leftBrkPoint = leftBreakPoint;
	newCluster->clusterSize = leftListSize + rightList->heapSize;
	newCluster->readMappingPtrArray = (DivetRow **) getMem (newCluster->clusterSize * sizeof (DivetRow *));
	newCluster->readMappingIdArray = (int *) getMem (newCluster->clusterSize * sizeof (int));

	tempPtr = leftList;
	for (count = 0; count < leftListSize; count++)
	{
		newCluster->readMappingIdArray[count] = tempPtr->readMappingPtr->divetRowId;
		newCluster->readMappingPtrArray[count] = tempPtr->readMappingPtr;

		tempPtr = tempPtr->next;
	}

	for (count = 0; count < rightList->heapSize; count++)
	{
		newCluster->readMappingIdArray[count + leftListSize] = rightList->heapArray[count].readMappingPtr->divetRowId;
		newCluster->readMappingPtrArray[count + leftListSize] = rightList->heapArray[count].readMappingPtr;
	}

	qsort (newCluster->readMappingIdArray, newCluster->clusterSize, sizeof( int), vh_compareInt);
	ptrToOldClusterList = g_listPotClusterFound;

	while (ptrToOldClusterList != NULL && newClusterIsMaximal == 1)
	{
		if (newClusterIsMaximal == 1)
			newClusterIsMaximal = (!vh_isItSubset(newCluster->readMappingIdArray, newCluster->clusterSize,
					ptrToOldClusterList->readMappingIdArray, ptrToOldClusterList->clusterSize));
		if (newClusterIsMaximal == 0)
			ptrToOldClusterList->leftBrkPoint = leftBreakPoint;

		oldClusterIsMaximal = (!vh_isItSubset( ptrToOldClusterList->readMappingIdArray,  ptrToOldClusterList->clusterSize,
				newCluster->readMappingIdArray,  newCluster->clusterSize));

		if (oldClusterIsMaximal == 0 && newClusterIsMaximal == 1)
			ptrToOldClusterList->isMaximalCluster = 0;

		ptrToOldClusterList = ptrToOldClusterList->next;
	}

	if (newClusterIsMaximal == 1)
	{
		newCluster->next = g_listPotClusterFound;
		g_listPotClusterFound = newCluster;
		ptrToOldClusterList = g_listPotClusterFound;
	}

	if (newClusterIsMaximal == 0)
	{
		free (newCluster->readMappingIdArray);
		free (newCluster->readMappingPtrArray);
		free (newCluster);
		newCluster = NULL;
	}
	vh_flushOut ( leftBreakPoint, SVtype);
}

/* To merge FR and RF into one cluster */
int addToOneSideInitOutputInterDup( Heap *heapName, int breakPointEnd, int brkPoint, char SVTYPE)
{
	int count;

	if (heapName->heapArray[0].readMappingPtr->orientationLeft == FORWARD
			&& heapName->heapArray[0].readMappingPtr->orientationRight == REVERSE)
	{
		PEAlistEls *newEl;

		if (leftSidePEAPicked[breakPointEnd] != NULL)
			vh_freeLinkList2( leftSidePEAPicked[breakPointEnd]);

		leftSidePEAPicked[breakPointEnd] = NULL;

		for (count = 0; count < heapName->heapSize; count++)
		{
			newEl = (PEAlistEls *) getMem(sizeof(PEAlistEls));
			newEl->brkPoint = brkPoint;
			newEl->readMappingPtr = heapName->heapArray[count].readMappingPtr;

			newEl->next = leftSidePEAPicked[breakPointEnd];
			leftSidePEAPicked[breakPointEnd] = newEl;
		}

		posOfBrkPointOfFR_Cluster[countPosBrkPointFR_Clusters] = breakPointEnd;
		countPosBrkPointFR_Clusters++;
		if (countPosBrkPointFR_Clusters > size_of_brkpoint_clusters){
		  posOfBrkPointOfFR_Cluster = resize_int_array ( posOfBrkPointOfFR_Cluster, size_of_brkpoint_clusters, size_of_brkpoint_clusters+BRKPOINTCLUSTERS);
		  size_of_brkpoint_clusters += BRKPOINTCLUSTERS;
		}
	}
	else if (heapName->heapArray[0].readMappingPtr->orientationLeft == REVERSE
			&& heapName->heapArray[0].readMappingPtr->orientationRight == FORWARD)
	{
		//fprintf(stderr, "%d ", countPosOBBrkPointFR_Clusters);
		if (countPosBrkPointFR_Clusters > 0)
		{
			for (count = 0; count < countPosBrkPointFR_Clusters; count++)
			{
				if (leftSidePEAPicked[posOfBrkPointOfFR_Cluster[count]]->brkPoint == brkPoint)
					addToPotentialOutputInterDup( leftSidePEAPicked[posOfBrkPointOfFR_Cluster[count]], heapName, brkPoint, SVTYPE);
			}
		}
	}
}

void vh_createBreakPointIntervals_InterDup( int brkPoint, Heap *h, char SVTYPE)
{
	int count, remSegLenMin, remSegLenMax;
	bool newElAdded;
	bool check = true;
	int locBrkPointLeftTemp, locBrkPointRightTemp; /* the min, max breakpoint positions for right side and left side */
	DivetRow *divet_row;

	g_listRightBrkPointIntrCount = 0;

	for( count = 0; count < h->heapSize; count++)
	{
		divet_row = h->heapArray[count].readMappingPtr;

		remSegLenMin = divet_row->libInfo->minDelta - ( brkPoint - divet_row->startPosition);
		remSegLenMax = divet_row->libInfo->maxDelta - ( brkPoint - divet_row->startPosition);

		if( remSegLenMin < 0)
			remSegLenMin = 0;
		if( remSegLenMax < 0)
		{
			/* Skip */
			check = false;
		}

		locBrkPointLeftTemp = divet_row->endPosition - remSegLenMax;
		locBrkPointRightTemp = divet_row->endPosition - remSegLenMin;

		if( check == true)
		{
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].readMappingPtr = divet_row;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].locBrkPointLeft = locBrkPointLeftTemp;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].locBrkPointRight = locBrkPointRightTemp;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].key = locBrkPointLeftTemp;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].keyLorR = LEFT;
			g_listRightBrkPointIntrCount++;

			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].readMappingPtr = divet_row;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].locBrkPointLeft = locBrkPointLeftTemp;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].locBrkPointRight = locBrkPointRightTemp;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].key = locBrkPointRightTemp;
			g_listRightBrkPointIntr[g_listRightBrkPointIntrCount].keyLorR = RIGHT;
			g_listRightBrkPointIntrCount++;
		}
	}

	newElAdded = false;
	if( g_listRightBrkPointIntrCount > 0)
		qsort( g_listRightBrkPointIntr, g_listRightBrkPointIntrCount , sizeof (struct RightBrkPointInterval), vh_compare);

	for( count = 0; count < g_listRightBrkPointIntrCount; count++)
	{
		if( g_listRightBrkPointIntr[count].keyLorR == LEFT)
		{
			vh_addToHeap( g_listRightBrkPointIntr[count].readMappingPtr, g_listRightBrkPointIntr[count].locBrkPointRight, g_intersectInterval);
			newElAdded = true;
		}
		else if( g_listRightBrkPointIntr[count].keyLorR == RIGHT)
		{
			if( newElAdded == true)
				addToOneSideInitOutputInterDup( g_intersectInterval, g_listRightBrkPointIntr[count].locBrkPointRight, brkPoint, SVTYPE);

			newElAdded = false;
			vh_heap_remove_top( g_intersectInterval);
		}
	}
}

void vh_createInterDupClusters (int chroSize, int inter_dup_location)
{
	int leftBreakPoint;
	int window_start_index, window_end_index, priority;
	int limit, i, j;
	g_listRightBrkPointIntrCount = 0;
	MappingOnGenome *tmp_mapping;

	//	Heap *FRHeap = ( Heap *) getMem ( sizeof ( Heap));
	//	Heap *RFHeap = ( Heap *) getMem ( sizeof ( Heap));
	Heap *FRHeap = vh_newHeap();
	Heap *RFHeap = vh_newHeap();


	g_intersectInterval->heapSize = 0;

	if( inter_dup_location == LEFTSIDE)
	{
		leftBreakPoint = 1;
		limit = chroSize;
	}
	else if( inter_dup_location == RIGHTSIDE)
	{
		leftBreakPoint = 0 - g_maxDeltaAmongLibs;
		limit = max_chromosome_size;
	}

	for( ; leftBreakPoint < limit; leftBreakPoint++)
	{
		countPosBrkPointFR_Clusters = 0;
		window_start_index = leftBreakPoint - g_maxDeltaAmongLibs;
		if( window_start_index < 0)
			window_start_index = 0;

		window_end_index = leftBreakPoint + g_maxDeltaAmongLibs;
		if( window_end_index >= chroSize)
			window_end_index = chroSize - 1;

		if( ( leftBreakPoint > 0) && g_genomeIndexStartFR[leftBreakPoint - 1] != NULL)
		{
			tmp_mapping = g_genomeIndexStartFR[leftBreakPoint - 1];
			while( tmp_mapping != NULL)
			{
				priority = tmp_mapping->readMappingPtr->startPosition;
				vh_addToHeap( tmp_mapping->readMappingPtr, priority, FRHeap);
				tmp_mapping = tmp_mapping->next;
			}
		}

		if( ( leftBreakPoint > 0) && g_genomeIndexStartRF[leftBreakPoint - 1] != NULL)
		{
			tmp_mapping =  g_genomeIndexStartRF[leftBreakPoint - 1];
			while( tmp_mapping != NULL)
			{
				priority = tmp_mapping->readMappingPtr->startPosition;
				vh_addToHeap( tmp_mapping->readMappingPtr, priority, RFHeap);
				tmp_mapping = tmp_mapping->next;
			}
		}

		if( ( window_start_index > 0) && ( g_genomeIndexStartFR[window_start_index - 1] != NULL))
			remove_Heap( FRHeap, window_start_index - 1);

		if( ( window_start_index > 0) && ( g_genomeIndexStartRF[window_start_index - 1] != NULL))
			remove_Heap( RFHeap, window_start_index - 1);

		//fprintf(stderr,"%d, %d, %d; Added %d, Removed %d; %d - %d\n", window_start_index, leftBreakPoint, window_end_index, i, j, FRHeap->heapSize, RFHeap->heapSize);

		if( FRHeap->heapSize > 0)
		{
			if( inter_dup_location == LEFTSIDE)
				vh_createBreakPointIntervals_InterDup( leftBreakPoint, FRHeap, INTERDUPLEFT);
			else if( inter_dup_location == RIGHTSIDE)
				vh_createBreakPointIntervals_InterDup( leftBreakPoint, FRHeap, INTERDUPRIGHT);
		}

		if( RFHeap->heapSize > 0)
		{
			if( inter_dup_location == LEFTSIDE)
				vh_createBreakPointIntervals_InterDup( leftBreakPoint, RFHeap, INTERDUPLEFT);
			else if( inter_dup_location == RIGHTSIDE)
				vh_createBreakPointIntervals_InterDup( leftBreakPoint, RFHeap, INTERDUPRIGHT);
		}
	}

	if( inter_dup_location == LEFTSIDE)
		vh_flushOut( leftBreakPoint - 1, INTERDUPLEFT);
	else if( inter_dup_location == RIGHTSIDE)
		vh_flushOut( leftBreakPoint - 1, INTERDUPRIGHT);

	//	free( FRHeap)
	vh_free_heap( FRHeap);
	FRHeap = NULL;

	//	free( RFHeap);
	vh_free_heap( RFHeap);       
	RFHeap = NULL;
}

void vh_addToGenomeIndex_InterDup (char *chromosome_name, sonic *this_sonic, int chroSize, int interdup_location)
{
	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	MappingOnGenome *newEl;
	int leftWindowStart;
	libInfo = g_libInfo;
	int count = 0;
	while (libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while (divetReadMappingPtr != NULL)
		{
			if (strcmp (divetReadMappingPtr->chromosome_name, chromosome_name) == 0
					&& !sonic_is_gap (this_sonic, chromosome_name, divetReadMappingPtr->locMapLeftEnd, divetReadMappingPtr->locMapRightStart)
					&& ( divetReadMappingPtr->locMapRightEnd < chroSize) && ( divetReadMappingPtr->locMapLeftStart > 0)
					&& ( divetReadMappingPtr->locMapRightStart - divetReadMappingPtr->locMapLeftEnd < maxDuplicationLen))
			{
				if( divetReadMappingPtr->svType == DELETION || divetReadMappingPtr->svType == INSERTION
						|| divetReadMappingPtr->svType == TANDEMDUP)
				{
					newEl = (MappingOnGenome *) getMem (sizeof (MappingOnGenome));
					newEl->readMappingPtr = divetReadMappingPtr;
					newEl->next = NULL;

					if( interdup_location == LEFTSIDE)
					{
						if( divetReadMappingPtr->orientationLeft == FORWARD && divetReadMappingPtr->orientationRight == REVERSE)
						{
							newEl->readMappingPtr->startPosition = divetReadMappingPtr->locMapLeftEnd;
							newEl->readMappingPtr->endPosition = divetReadMappingPtr->locMapRightStart;
							newEl->readMappingPtr->divetRowId = count++;

							leftWindowStart = divetReadMappingPtr->locMapLeftEnd;
							newEl->next = g_genomeIndexStartFR[leftWindowStart];
							g_genomeIndexStartFR[leftWindowStart] = newEl;
						}
						else if( divetReadMappingPtr->orientationLeft == REVERSE && divetReadMappingPtr->orientationRight == FORWARD)
						{
							newEl->readMappingPtr->startPosition = divetReadMappingPtr->locMapLeftStart;
							newEl->readMappingPtr->endPosition = divetReadMappingPtr->locMapRightEnd;
							newEl->readMappingPtr->divetRowId = count++;

							leftWindowStart = divetReadMappingPtr->locMapLeftStart;
							newEl->next = g_genomeIndexStartRF[leftWindowStart];
							g_genomeIndexStartRF[leftWindowStart] = newEl;
						}
					}
					else if( interdup_location == RIGHTSIDE)
					{
						if( divetReadMappingPtr->orientationLeft == FORWARD && divetReadMappingPtr->orientationRight == REVERSE)
						{
							newEl->readMappingPtr->startPosition = max_chromosome_size - divetReadMappingPtr->locMapRightStart;
							newEl->readMappingPtr->endPosition = max_chromosome_size - divetReadMappingPtr->locMapLeftEnd;
							newEl->readMappingPtr->divetRowId = count++;

							newEl->readMappingPtr->orientationLeft = REVERSE;
							newEl->readMappingPtr->orientationRight = FORWARD;

							leftWindowStart = max_chromosome_size - divetReadMappingPtr->locMapRightStart;
							newEl->next = g_genomeIndexStartRF[leftWindowStart];
							g_genomeIndexStartRF[leftWindowStart] = newEl;
						}
						else if( divetReadMappingPtr->orientationLeft == REVERSE && divetReadMappingPtr->orientationRight == FORWARD)
						{
							newEl->readMappingPtr->startPosition = max_chromosome_size - divetReadMappingPtr->locMapRightEnd;
							newEl->readMappingPtr->endPosition = max_chromosome_size - divetReadMappingPtr->locMapLeftStart;
							newEl->readMappingPtr->divetRowId = count++;

							newEl->readMappingPtr->orientationLeft = FORWARD;
							newEl->readMappingPtr->orientationRight = REVERSE;

							leftWindowStart = max_chromosome_size - divetReadMappingPtr->locMapRightEnd;
							newEl->next = g_genomeIndexStartFR[leftWindowStart];
							g_genomeIndexStartFR[leftWindowStart] = newEl;
						}
					}
				}
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
	}
	//fprintf(stderr,"Count %d\n", count);
}

void vh_initializeReadMapping_InterDup (sonic *this_sonic, int chr_index, int interdup_location)
{
	int i;
	LibraryInfo *libInfoPtr = g_libInfo;

	while (libInfoPtr != NULL)
	{
		g_maxDeltaAmongLibs = vh_max (g_maxDeltaAmongLibs, libInfoPtr->maxDelta);
		libInfoPtr = libInfoPtr->next;
	}

	g_listPotClusterFound = NULL;
	max_chromosome_size = this_sonic->chromosome_lengths[chr_index] + g_maxDeltaAmongLibs;

	/* Initializing the Genome Array */
	g_genomeIndexStartFR = (MappingOnGenome **) getMem ( max_chromosome_size * sizeof (MappingOnGenome *));
	g_genomeIndexStartRF = (MappingOnGenome **) getMem ( max_chromosome_size * sizeof (MappingOnGenome *));

	leftSidePEAPicked = (PEAlistEls **) getMem( ( max_chromosome_size + 1000) * sizeof(PEAlistEls **));

	size_of_brkpoint_clusters = BRKPOINTCLUSTERS;
	posOfBrkPointOfFR_Cluster = (int *) getMem(size_of_brkpoint_clusters * sizeof(int));
	
	countPosBrkPointFR_Clusters = 0;

	for (i = 0; i < max_chromosome_size + 1000; i++)
		leftSidePEAPicked[i] = NULL;

	if( g_genomeIndexStartFR == NULL || g_genomeIndexStartRF == NULL)
		vh_logWarning ("Memory Problem in vh_initializeReadMapping_InterDup()");

	for( i = 0; i < max_chromosome_size; i++)
	{
		g_genomeIndexStartFR[i] = NULL;
		g_genomeIndexStartRF[i] = NULL;
	}

	vh_addToGenomeIndex_InterDup( this_sonic->chromosome_names[chr_index], this_sonic, this_sonic->chromosome_lengths[chr_index], interdup_location);

	/* Initializing the list of begin and end of right side break point ranges */
	g_listRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));

	for( i = 0; i < g_maxListBrkPointIntr; i++)
		g_listRightBrkPointIntr[i].readMappingPtr = NULL;

	g_intersectInterval = vh_newHeap();

}
