#include "vh_maximalCluster.h"
#include "vh_common.h"
#include "vh_setcover.h"

#define MAXCLUSTERLIST 100000

MappingOnGenome **g_genomeIndexStart;
RightBrkPointInterval *g_listRightBrkPointIntr;
RightBrkPointInterval *g_tempListRightBrkPointIntr;
int g_maxListBrkPointIntr;
int g_listRightBrkPointIntrCount;
int g_maxDeltaAmongLibs = 0;
int max_chromosome_size;
ClustersFound *g_listPotClusterFound;

int test;
int vh_min (int x, int y)
{
	if (x < y)
		return x;
	else
		return y;
}


int vh_max (int x, int y)
{
	if (x < y)
		return y;
	else
		return x;
}

int vh_compare2 (const void *a, const void *b)
{

	struct RightBrkPointInterval *arg1 = (struct RightBrkPointInterval *) a;
	struct RightBrkPointInterval *arg2 = (struct RightBrkPointInterval *) b;
	if ((*arg1).readMappingPtr->locMapLeftEnd > (*arg2).readMappingPtr->locMapLeftEnd)
		return 1;
	if ((*arg1).readMappingPtr->locMapLeftEnd < (*arg2).readMappingPtr->locMapLeftEnd)
		return -1;
}


int vh_compare (const void *a, const void *b)
{

	struct RightBrkPointInterval *arg1 = (struct RightBrkPointInterval *) a;
	struct RightBrkPointInterval *arg2 = (struct RightBrkPointInterval *) b;
	if ((*arg1).key > (*arg2).key)
		return 1;
	if ((*arg1).key < (*arg2).key)
		return -1;
	/*	if ((*arg1).key==(*arg1).locBrkPointLeft)
	{
	return 1;
	} else {
	return -1;
	}
	 */

	if ((*arg1).key == (*arg2).key)
	{
		if ((*arg1).keyLorR == LEFT)
			return -1;
		else
			return 1;
	}
}

int vh_compareReadName (const void *a, const void *b)
{
	return strcmp ((*(DivetRow **) a)->readName->readName, (*(DivetRow **) b)->readName->readName);
}

int vh_compareInt (const void *a, const void *b)
{
	return *(int *) a - *(int *) b;
}

int vh_comparePtr (const void *a, const void *b)
{

	if (strcmp((*(DivetRow **) a)->readName->readName, (*(DivetRow **) b)->readName->readName) == 0)
	{
		if ((*(DivetRow **) a)->locMapLeftEnd == (*(DivetRow **) b)->locMapLeftEnd)
		{
			return ((*(DivetRow **) a)->locMapRightStart - (*(DivetRow **) b)->locMapRightStart);
		}
		else
			return (*(DivetRow **) a)->locMapLeftEnd - (*(DivetRow **) b)->locMapLeftEnd;
	}
	else
		return strcmp ((*(DivetRow **) a)->readName->readName, (*(DivetRow **) b)->readName->readName);
}


void vh_freeLinkedList (MappingOnGenome * cur)
{
	MappingOnGenome *next;
	while (cur != NULL)
	{
		next = cur->next;
		free (cur);
		cur = next;
	}
}
int vh_freeLinkList2(PEAlistEls *ptr)
{
	PEAlistEls* ptr2;

	while(ptr!=NULL)
	{
		ptr2=ptr->next;
		free(ptr);
		ptr=ptr2;
	}
}

void vh_finalizeReadMapping (char *chromosome_name, int chroSize)
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

	for (i = 0; i < chroSize; i++)
	{
		vh_freeLinkedList (g_genomeIndexStart[i]);
		g_genomeIndexStart[i] = NULL;
	}
	free (g_genomeIndexStart);
	g_genomeIndexStart = NULL;

	//Free g_listRightBrkPointIntrCount and set the g_listRightBrkPointIntrCount
	free (g_listRightBrkPointIntr);
	g_listRightBrkPointIntr = NULL;
	g_listRightBrkPointIntrCount = 0;

	//Free g_tempListRightBrkPointIntr
	free (g_tempListRightBrkPointIntr);
	g_tempListRightBrkPointIntr = NULL;

	//Free g_intersectInterval -> Heap AND set the heapsize
	free (g_intersectInterval);
	g_intersectInterval = NULL;

}

void vh_copyElBrkPointIntr (int dest, int src)
{
	g_tempListRightBrkPointIntr[dest].locBrkPointLeft = g_listRightBrkPointIntr[src].locBrkPointLeft;
	g_tempListRightBrkPointIntr[dest].locBrkPointRight = g_listRightBrkPointIntr[src].locBrkPointRight;
	g_tempListRightBrkPointIntr[dest].key = g_listRightBrkPointIntr[src].key;
	g_tempListRightBrkPointIntr[dest].readMappingPtr = g_listRightBrkPointIntr[src].readMappingPtr;
	g_tempListRightBrkPointIntr[dest].keyLorR = g_listRightBrkPointIntr[src].keyLorR;
}

/* Return 0 if it is not subset */
int vh_isItSubset (int *querySet, int querySetSize, int *patternSet, int patternSetSize)
{
	int idPatternSet = 0;
	int idQuerySet = 0;
	if (patternSetSize == 0)
		return 0;
	while (idQuerySet < querySetSize && idPatternSet < patternSetSize)
	{
		if (querySet[idQuerySet] == patternSet[idPatternSet])
		{
			idPatternSet++;
			idQuerySet++;
		}
		else if (querySet[idQuerySet] > patternSet[idPatternSet])
		{
			idPatternSet++;
		}
		else if (querySet[idQuerySet] < patternSet[idPatternSet])
			return 0;
	}

	if (idPatternSet == patternSetSize && idQuerySet < querySetSize)
		return 0;
	else
		return 1;
}

void vh_flushOut (int leftBreakPoint, char SVtype)
{
	ClustersFound *ptrToOldClusterList;
	ClustersFound *tempPtrOldCluster;
	ptrToOldClusterList = g_listPotClusterFound;

	while (ptrToOldClusterList != NULL)
	{
		if (ptrToOldClusterList->next != NULL && ptrToOldClusterList->next->isMaximalCluster == 0)
		{
			tempPtrOldCluster = ptrToOldClusterList->next->next;
			free (ptrToOldClusterList->next->readMappingPtrArray);
			free (ptrToOldClusterList->next->readMappingIdArray);
			free (ptrToOldClusterList->next);
			ptrToOldClusterList->next = tempPtrOldCluster;
		}
		else if (ptrToOldClusterList->next != NULL && leftBreakPoint > ptrToOldClusterList->next->leftBrkPoint + g_maxDeltaAmongLibs)
		{
			vh_outputCluster( ptrToOldClusterList->next, SVtype);
			tempPtrOldCluster = ptrToOldClusterList->next->next;
			free (ptrToOldClusterList->next->readMappingPtrArray);
			free (ptrToOldClusterList->next->readMappingIdArray);
			free (ptrToOldClusterList->next);
			ptrToOldClusterList->next = tempPtrOldCluster;
		}
		else
			ptrToOldClusterList = ptrToOldClusterList->next;
	}
	if (g_listPotClusterFound != NULL && g_listPotClusterFound->isMaximalCluster == 0)
	{
		tempPtrOldCluster = g_listPotClusterFound->next;
		free (g_listPotClusterFound->readMappingIdArray);
		free (g_listPotClusterFound->readMappingPtrArray);
		free (g_listPotClusterFound);
		g_listPotClusterFound = tempPtrOldCluster;
	}
	else if (g_listPotClusterFound != NULL && leftBreakPoint > g_listPotClusterFound->leftBrkPoint + g_maxDeltaAmongLibs)
	{
		vh_outputCluster (g_listPotClusterFound, SVtype);
		tempPtrOldCluster = g_listPotClusterFound->next;
		free (g_listPotClusterFound->readMappingPtrArray);
		free (g_listPotClusterFound->readMappingIdArray);
		free (g_listPotClusterFound);
		g_listPotClusterFound = tempPtrOldCluster;
	}
}

void vh_addToPotentialOutput (int leftBreakPoint, Heap *heapName, char SVtype)
{
	int count;
	ClustersFound *newCluster;
	ClustersFound *ptrToOldClusterList;
	ClustersFound *tempPtrOldCluster;
	int newClusterIsMaximal = 1;
	int oldClusterIsMaximal = 1;

	newCluster = (ClustersFound *) getMem (sizeof (ClustersFound));
	newCluster->next = NULL;
	newCluster->isMaximalCluster = 1;
	newCluster->leftBrkPoint = leftBreakPoint;
	newCluster->clusterSize = heapName->heapSize;
	newCluster->readMappingPtrArray = (DivetRow **) getMem (newCluster->clusterSize * sizeof (DivetRow *));
	newCluster->readMappingIdArray = (int *) getMem (newCluster->clusterSize * sizeof (int));
	for (count = 0; count < heapName->heapSize; count++)
	{
		newCluster->readMappingPtrArray[count] = heapName->heapArray[count].readMappingPtr;
		newCluster->readMappingIdArray[count] = heapName->heapArray[count].readMappingPtr->divetRowId;
	}

	qsort (newCluster->readMappingIdArray, newCluster->clusterSize, sizeof (int), vh_compareInt);
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

int **create_list_of_reads(int list_size)
{
  int **ret;
  int i;

  ret = (int **) getMem(list_size * sizeof (int *));

  for (i=0; i<list_size; i++){
    ret[i] = (int *) getMem(3 * sizeof (int));
  }

  return ret;
  
}

int **recreate_list_of_reads(int **old_list, int list_size, int new_list_size)
{
  int **ret;
  int i;

  ret = (int **) getMem(new_list_size * sizeof (int *));

  for (i=0; i<list_size; i++){
    ret[i] = (int *) getMem(3 * sizeof (int));
    ret[i][0] = old_list[i][0];
    ret[i][1] = old_list[i][1];
    ret[i][2] = old_list[i][2];
    freeMem (old_list[i], (3 * sizeof(int *)));
  }
  
  for (i=list_size; i<new_list_size; i++){
    ret[i] = (int *) getMem(3 * sizeof (int));
  }
  

  freeMem (old_list, list_size * sizeof (int *));
  
  return ret;
  
}

int vh_outputCluster (ClustersFound * cluster, char SVtype)
{
        int **list_of_written_reads;//[MAXCLUSTERLIST][3];
	// [0]:locMapLeftEnd, [1]: locMapRightStart, [2]: edit distance. This is used to remove the duplicated reads (clonal) as a result of PCR duplication from each cluster.
	int totalAddedToList = 0;
	int clonalRead = 0;
	int readMapCount;
	int countListOutputed;
	int start, end;
	int written = 0;
	clusters_final *cluster_new, *tmp, *prev;

	int cnt = 0;
	int list_size = MAXCLUSTERLIST;
	
	
	/* Count number of read-pairs inside the clusters excluding the ones >= cluster_of_reads */
	for (readMapCount = 0; readMapCount < cluster->clusterSize; readMapCount++)
	{
		if(cluster->readMappingPtrArray[readMapCount]->in_cluster_count < cluster_of_reads)
			cnt++;
	}

	if (cnt < 2 || (SVtype == INVERSION && vh_notBothDirections (cluster)))
		return 0;

	//if (cluster->clusterSize < 2 || (SVtype == INVERSION && vh_notBothDirections (cluster)))
	//return 0;

	list_of_written_reads = create_list_of_reads(list_size);

	clusters_all[cluster_count] = NULL;

	qsort (cluster->readMappingPtrArray, cluster->clusterSize, sizeof (DivetRow **), vh_compareReadName);

	for (readMapCount = 0; readMapCount < cluster->clusterSize; readMapCount++)
	{
		if(cluster->readMappingPtrArray[readMapCount]->in_cluster_count >= cluster_of_reads)
		{
			//fprintf(stderr,"REMOVED %s - %d\n", cluster->readMappingPtrArray[readMapCount]->readName->readName, cluster->readMappingPtrArray[readMapCount]->cluster_count);
			continue;
		}

		clonalRead = 0;
		if (readMapCount == 0 || strcmp (cluster->readMappingPtrArray[readMapCount]->readName->readName,
				cluster->readMappingPtrArray[readMapCount -1]->readName->readName) != 0)
		{
			for (countListOutputed = 0; countListOutputed < totalAddedToList; countListOutputed++)
			{
				if ((cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd == list_of_written_reads[countListOutputed][0])
						&& (cluster->readMappingPtrArray[readMapCount]->locMapRightStart == list_of_written_reads[countListOutputed][1])
						&& (cluster->readMappingPtrArray[readMapCount]->editDistance == list_of_written_reads[countListOutputed][2]))
					clonalRead = 1;
			}
			if (clonalRead == 0)
			{
				if( SVtype == INVDUPLEFT || SVtype == INTERDUPLEFT)
				{
					start = cluster->readMappingPtrArray[readMapCount]->startPosition;
					end = cluster->readMappingPtrArray[readMapCount]->endPosition;
				}
				else if( SVtype == INVDUPRIGHT || SVtype == INTERDUPRIGHT)
				{
					start = max_chromosome_size - cluster->readMappingPtrArray[readMapCount]->endPosition;
					end = max_chromosome_size - cluster->readMappingPtrArray[readMapCount]->startPosition;
				}
				else
				{
					start = cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd;
					end = cluster->readMappingPtrArray[readMapCount]->locMapRightStart;
				}

				if( debug_mode)
				{
					if (ten_x_flag != 1 && output_hs_flag != 1)
					{
						/* Output to .clusters file */
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->readName->readName);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->chromosome_name);
						fprintf (fileOutput, "%i ", start);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->chromosome_name);
						fprintf (fileOutput, "%i ", end);
						fprintf (fileOutput, "%c ", SVtype);
						fprintf (fileOutput, "%g ", cluster->readMappingPtrArray[readMapCount]->phredScore);
						fprintf (fileOutput, "%g ", cluster->readMappingPtrArray[readMapCount]->editDistance);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->libInfo->libName);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->libInfo->indName);
						fprintf (fileOutput, "%c ", cluster->readMappingPtrArray[readMapCount]->orientationLeft);
						fprintf (fileOutput, "%c ", cluster->readMappingPtrArray[readMapCount]->orientationRight);
						fprintf (fileOutput, "%d ", cluster->readMappingPtrArray[readMapCount]->mQual1);
						fprintf (fileOutput, "%d ", cluster->readMappingPtrArray[readMapCount]->mQual2);
					}
					else
					{
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->readName->readName);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->chromosome_name);
						fprintf (fileOutput, "%i ", start);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->chromosome_name);
						fprintf (fileOutput, "%i ", end);
						fprintf (fileOutput, "%c ", SVtype);
						fprintf (fileOutput, "%g ", cluster->readMappingPtrArray[readMapCount]->phredScore);
						fprintf (fileOutput, "%g ", cluster->readMappingPtrArray[readMapCount]->editDistance);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->libInfo->libName);
						fprintf (fileOutput, "%s ", cluster->readMappingPtrArray[readMapCount]->libInfo->indName);
						fprintf (fileOutput, "%c ", cluster->readMappingPtrArray[readMapCount]->orientationLeft);
						fprintf (fileOutput, "%c ", cluster->readMappingPtrArray[readMapCount]->orientationRight);
						fprintf (fileOutput, "%d ", cluster->readMappingPtrArray[readMapCount]->mQual1);
						fprintf (fileOutput, "%d ", cluster->readMappingPtrArray[readMapCount]->mQual2);
						fprintf (fileOutput, "%lu ", cluster->readMappingPtrArray[readMapCount]->ten_x_barcode);
					}
				}
				/* Fill the clusters struct */
				cluster_new = ( clusters_final *) getMem( sizeof( clusters_final));

				cluster_new->read_name = NULL;
				set_str( &cluster_new->read_name, cluster->readMappingPtrArray[readMapCount]->readName->readName);

				cluster_new->chromosome_name1 = NULL;
				set_str( &cluster_new->chromosome_name1, cluster->readMappingPtrArray[readMapCount]->chromosome_name);

				cluster_new->start_position = start;

				cluster_new->chromosome_name2 = NULL;
				set_str( &cluster_new->chromosome_name2, cluster->readMappingPtrArray[readMapCount]->chromosome_name);

				cluster_new->end_position = end;
				cluster_new->SV_type = SVtype;
				cluster_new->phred_score = cluster->readMappingPtrArray[readMapCount]->phredScore;
				cluster_new->edit_distance = cluster->readMappingPtrArray[readMapCount]->editDistance;

				cluster_new->library_name = NULL;
				set_str( &cluster_new->library_name, cluster->readMappingPtrArray[readMapCount]->libInfo->libName);

				cluster_new->individual_name = NULL;
				set_str( &cluster_new->individual_name, cluster->readMappingPtrArray[readMapCount]->libInfo->indName);

				cluster_new->orientation_left = cluster->readMappingPtrArray[readMapCount]->orientationLeft;
				cluster_new->orientation_right = cluster->readMappingPtrArray[readMapCount]->orientationRight;
				cluster_new->mapping_quality_left = cluster->readMappingPtrArray[readMapCount]->mQual1;
				cluster_new->mapping_quality_right = cluster->readMappingPtrArray[readMapCount]->mQual2;

				if( ten_x_flag != 1 && output_hs_flag != 1)
					cluster_new->ten_x_barcode = -1;
				else
					cluster_new->ten_x_barcode = cluster->readMappingPtrArray[readMapCount]->ten_x_barcode;

				cluster_new->mei_subclass = NULL;
				cluster_new->mei_type = NULL;

				if( ( clusters_all[cluster_count] == NULL) || ( clusters_all[cluster_count]->edit_distance >= cluster_new->edit_distance))
				{
					cluster_new->next = clusters_all[cluster_count];
					clusters_all[cluster_count] = cluster_new;
				}
				else
				{
					tmp = clusters_all[cluster_count];
					while( ( tmp->next != NULL) && ( tmp->next->edit_distance < cluster_new->edit_distance))
						tmp = tmp->next;

					cluster_new->next = tmp->next;
					tmp->next = cluster_new;
				}
				written = 1;

				/* Keep the number of clusters that a read is involved in */
				cluster->readMappingPtrArray[readMapCount]->in_cluster_count++;

				list_of_written_reads[totalAddedToList][0] = cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd;
				list_of_written_reads[totalAddedToList][1] = cluster->readMappingPtrArray[readMapCount]->locMapRightStart;
				list_of_written_reads[totalAddedToList][2] = (int) cluster->readMappingPtrArray[readMapCount]->editDistance;
				totalAddedToList++;
				if (totalAddedToList >= list_size) {
				  fprintf( stderr, "totaladdedto list: %d\n", totalAddedToList);
				  list_of_written_reads = recreate_list_of_reads(list_of_written_reads, list_size, list_size+MAXCLUSTERLIST);
				  list_size+=MAXCLUSTERLIST;
				}
				/* ARDA: there is no control here to check if totalAddedToList < MAXCLUSTERLIST. What happens if it exceeds? */
			}
		}
	}
	if(written)
	{
		if( debug_mode)
			fprintf (fileOutput, "END\n");
		cluster_count++;
	}

	for (cnt=0; cnt<list_size; cnt++){
	  freeMem (list_of_written_reads[cnt], (3 * sizeof(int *)));
	}
	
	freeMem (list_of_written_reads, list_size * sizeof (int *));
	
}

void vh_createIntersectingIntervals (int leftBreakPoint, char SVtype)
{
	int newElAdded = 0, countIntrEndPoints;

	qsort (g_listRightBrkPointIntr, g_listRightBrkPointIntrCount, sizeof (struct RightBrkPointInterval), vh_compare);

	for (countIntrEndPoints = 0; countIntrEndPoints < g_listRightBrkPointIntrCount; countIntrEndPoints++)
	{
		if (g_listRightBrkPointIntr[countIntrEndPoints].keyLorR == LEFT)
		{
			newElAdded = 1;
			vh_addToHeap (g_listRightBrkPointIntr[countIntrEndPoints].readMappingPtr, g_listRightBrkPointIntr[countIntrEndPoints].locBrkPointRight, g_intersectInterval);
		}
		else if (g_listRightBrkPointIntr[countIntrEndPoints].keyLorR == RIGHT)
		{
			if (vh_minValue_heap (g_intersectInterval) == g_listRightBrkPointIntr[countIntrEndPoints].key)
			{
				if (newElAdded == 1)
				{
					vh_addToPotentialOutput (leftBreakPoint, g_intersectInterval, SVtype);
					newElAdded = 0;
				}
				vh_heap_remove_top (g_intersectInterval);
			}/*
			else
			{
				printf ("An Error Occured\n");
				vh_writeHeap (g_intersectInterval);
			}*/
		}
	}
	g_intersectInterval->heapSize = 0;

}

int vh_notBothDirections (ClustersFound * cluster)
{
	int FF, RR, readMapCount;
	FF = 0;
	RR = 0;
	for (readMapCount = 0; readMapCount < cluster->clusterSize; readMapCount++)
	{
		if (cluster->readMappingPtrArray[readMapCount]->orientationLeft == FORWARD
				&& cluster->readMappingPtrArray[readMapCount]->orientationRight == FORWARD)
			FF = 1;

		if (cluster->readMappingPtrArray[readMapCount]->orientationLeft == REVERSE
				&& cluster->readMappingPtrArray[readMapCount]->orientationRight == REVERSE)
			RR = 1;

		if (FF && RR)
			return 0;
	}
	return 1;
}
