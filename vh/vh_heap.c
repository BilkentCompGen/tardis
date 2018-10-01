#include "vh_heap.h"
#include "vh_divethandler.h"

Heap *g_intersectInterval;

Heap *vh_newHeap(void)
{
       Heap *tmp_heap;
       tmp_heap = (Heap *) getMem (sizeof (Heap));
       tmp_heap->heapArray = (struct HeapEl *) getMem(sizeof(struct HeapEl) * MAX_CLUSTER_SIZE); // I really don't like this max_cluster_size thing
       return tmp_heap;
}

HeapMEI *vh_newHeapMEI(void)
{
       HeapMEI *tmp_heap;
       tmp_heap = (HeapMEI *) getMem (sizeof (HeapMEI));
       tmp_heap->heapArray = (struct HeapElMEI *) getMem(sizeof(struct HeapElMEI) * MAX_CLUSTER_SIZE); // I really don't like this max_cluster_size thing
       return tmp_heap;
}

HeapNUMT *vh_newHeapNUMT(void)
{
       HeapNUMT *tmp_heap;
       tmp_heap = (HeapNUMT *) getMem (sizeof (HeapNUMT));
       tmp_heap->heapArray = (struct HeapElNUMT *) getMem(sizeof(struct HeapElNUMT) * MAX_CLUSTER_SIZE); // I really don't like this max_cluster_size thing
       return tmp_heap;
}

void vh_free_heap(Heap *h)
{
  free(h->heapArray);
  free(h);
}

void vh_free_heap_mei(HeapMEI *h)
{
  free(h->heapArray);
  free(h);
}

void vh_free_heap_numt(HeapNUMT *h)
{
  free(h->heapArray);
  free(h);
}




void vh_writeHeap (Heap * heapName)
{
	int count;
	for (count = 0; count < heapName->heapSize; count++)
	{
		printf ("%i %i %i *", heapName->heapArray[count].priorityValue,
				heapName->heapArray[count].readMappingPtr->locMapLeftEnd,
				heapName->heapArray[count].readMappingPtr->locMapRightStart);
	}
	printf ("END\n");
}


void vh_copyHeapEl (HeapEl * dest, HeapEl * src)
{
	dest->readMappingPtr = src->readMappingPtr;
	dest->priorityValue = src->priorityValue;
}


void vh_push_heap (Heap * heapName, HeapEl * newEl)
{
	int heapIndex = heapName->heapSize;
	HeapEl tempEl;
	heapName->heapSize++;
	vh_copyHeapEl (&(heapName->heapArray[heapIndex]), newEl);
	heapIndex = heapName->heapSize - 1;;
	while (heapIndex > 0 && heapName->heapArray[heapIndex].priorityValue < heapName->heapArray[(heapIndex + 1) / 2 - 1].priorityValue)
	{
		vh_copyHeapEl (&tempEl, &heapName->heapArray[heapIndex]);
		vh_copyHeapEl (&(heapName->heapArray[heapIndex]), &(heapName->heapArray[(heapIndex + 1) / 2 - 1]));
		vh_copyHeapEl (&(heapName->heapArray[(heapIndex + 1) / 2 - 1]), &tempEl);
		heapIndex = (heapIndex + 1) / 2 - 1;
	}
}

int vh_minValue_heap (Heap * heapName)
{
	if (heapName->heapSize > 0)
		return (heapName->heapArray[0].priorityValue);
	else
		return -1000;
}

void vh_heapBubleDown (Heap * heapName)
{
	int heapIndex = 0;
	HeapEl tempEl;
	while ((heapIndex * 2 + 1 < heapName->heapSize && heapName->heapArray[heapIndex].priorityValue >
	heapName->heapArray[heapIndex * 2 + 1].priorityValue) || (heapIndex * 2 + 2 < heapName->heapSize
			&& heapName->heapArray[heapIndex].priorityValue >heapName->heapArray[heapIndex * 2 + 2].priorityValue))
	{
		if (heapIndex * 2 + 2 < heapName->heapSize)
		{
			if (heapName->heapArray[heapIndex * 2 + 1].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue)
			{
				vh_copyHeapEl (&tempEl, &(heapName->heapArray[heapIndex]));
				vh_copyHeapEl (&(heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 2]));
				vh_copyHeapEl (&(heapName->heapArray[heapIndex * 2 + 2]), &tempEl);
				heapIndex = heapIndex * 2 + 2;
			}
			else
			{
				vh_copyHeapEl (&tempEl, &(heapName->heapArray[heapIndex]));
				vh_copyHeapEl (&(heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
				vh_copyHeapEl (&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
				heapIndex = heapIndex * 2 + 1;
			}
		}
		else
		{
			vh_copyHeapEl (&tempEl, &(heapName->heapArray[heapIndex]));
			vh_copyHeapEl (&(heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
			vh_copyHeapEl (&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
			heapIndex = heapIndex * 2 + 1;
		}
	}
}

int remove_Heap(Heap * heapName, int position)
{
	while( vh_minValue_heap( heapName) == position)
		vh_heap_remove_top( heapName);
}

void vh_heap_remove_top (Heap * heapName)
{
	if (heapName->heapSize > 0)
	{
		vh_copyHeapEl (&(heapName->heapArray[0]), &(heapName->heapArray[heapName->heapSize - 1]));
		vh_heapBubleDown (heapName);
		heapName->heapSize--;
	}
}

void vh_addToHeap (DivetRow * readMappingPtr, int priorityValue, Heap * heapName)
{
	HeapEl *newEl;
	newEl = (HeapEl *) getMem (sizeof (HeapEl));
	newEl->readMappingPtr = readMappingPtr;
	newEl->priorityValue = priorityValue;

	vh_push_heap (heapName, newEl);
	free (newEl);
}

int minValue_heapMEI(HeapMEI *heapName)
{
	if (heapName->heapSize > 0)
		return (heapName->heapArray[0].priorityValue);
	else
		return -1000;
}

int minValue_heapNUMT(HeapNUMT *heapName)
{
	if (heapName->heapSize > 0)
		return (heapName->heapArray[0].priorityValue);
	else
		return -1000;
}

void copyHeapElMEI(HeapElMEI *dest, HeapElMEI *src)
{
	dest->mei_ReadsPtr = src->mei_ReadsPtr;
	dest->priorityValue = src->priorityValue;
}

void copyHeapElNUMT(HeapElNUMT *dest, HeapElNUMT *src)
{
	dest->numt_ReadsPtr = src->numt_ReadsPtr;
	dest->priorityValue = src->priorityValue;
}

void heapBubleDownMEI(HeapMEI *heapName)
{
	int heapIndex = 0;
	HeapElMEI tempEl;
	while((heapIndex * 2 + 1<heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 1].priorityValue)||(heapIndex * 2 + 2 < heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue))
	{
		if (heapIndex * 2 + 2<heapName->heapSize)
		{
			if (heapName->heapArray[heapIndex * 2 + 1].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue)
			{
				copyHeapElMEI(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapElMEI(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 2]));
				copyHeapElMEI(& (heapName->heapArray[heapIndex * 2 + 2]), &tempEl);
				heapIndex = heapIndex * 2 + 2;
			}
			else
			{
				copyHeapElMEI(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapElMEI(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
				copyHeapElMEI(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
				heapIndex = heapIndex * 2 + 1;
			}
		}
		else
		{
			copyHeapElMEI(&tempEl, & (heapName->heapArray[heapIndex]));
			copyHeapElMEI(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
			copyHeapElMEI(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
			heapIndex = heapIndex * 2 + 1;
		}
	}
}


void heapBubleDownNUMT(HeapNUMT *heapName)
{
	int heapIndex = 0;
	HeapElNUMT tempEl;
	while((heapIndex * 2 + 1<heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 1].priorityValue)||(heapIndex * 2 + 2 < heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue))
	{
		if (heapIndex * 2 + 2<heapName->heapSize)
		{
			if (heapName->heapArray[heapIndex * 2 + 1].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue)
			{
				copyHeapElNUMT(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapElNUMT(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 2]));
				copyHeapElNUMT(& (heapName->heapArray[heapIndex * 2 + 2]), &tempEl);
				heapIndex = heapIndex * 2 + 2;
			}
			else
			{
				copyHeapElNUMT(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapElNUMT(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
				copyHeapElNUMT(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
				heapIndex = heapIndex * 2 + 1;
			}
		}
		else
		{
			copyHeapElNUMT(&tempEl, & (heapName->heapArray[heapIndex]));
			copyHeapElNUMT(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
			copyHeapElNUMT(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
			heapIndex = heapIndex * 2 + 1;
		}
	}
}


void heap_remove_topMEI( HeapMEI *heapName)
{
	if( heapName->heapSize > 0)
	{
		copyHeapElMEI(&(heapName->heapArray[0]), &(heapName->heapArray[heapName->heapSize - 1]));
		heapBubleDownMEI( heapName);
		heapName->heapSize--;
	}
}

void heap_remove_topNUMT( HeapNUMT *heapName)
{
	if( heapName->heapSize > 0)
	{
		copyHeapElNUMT(&(heapName->heapArray[0]), &(heapName->heapArray[heapName->heapSize - 1]));
		heapBubleDownNUMT( heapName);
		heapName->heapSize--;
	}
}

void push_heap_mei(HeapMEI *heapName, HeapElMEI *newEl)
{
	int heapIndex = heapName->heapSize;
	HeapElMEI tempEl;
	heapName->heapSize++;
	copyHeapElMEI(&(heapName->heapArray[heapIndex]), newEl);
	heapIndex = heapName->heapSize - 1;

	while(heapIndex > 0 && heapName->heapArray[heapIndex].priorityValue < heapName->heapArray[(heapIndex + 1) / 2 -1].priorityValue)
	{
		copyHeapElMEI(&tempEl, &heapName->heapArray[heapIndex]);
		copyHeapElMEI(&(heapName->heapArray[heapIndex]),  &(heapName->heapArray[(heapIndex + 1) / 2 -1]));
		copyHeapElMEI(&(heapName->heapArray[(heapIndex + 1) / 2 -1]), &tempEl);
		heapIndex = (heapIndex + 1) / 2 - 1;
	}
}

void push_heap_numt(HeapNUMT *heapName, HeapElNUMT *newEl)
{
	int heapIndex = heapName->heapSize;
	HeapElNUMT tempEl;
	heapName->heapSize++;
	copyHeapElNUMT(&(heapName->heapArray[heapIndex]), newEl);
	heapIndex = heapName->heapSize - 1;

	while(heapIndex > 0 && heapName->heapArray[heapIndex].priorityValue < heapName->heapArray[(heapIndex + 1) / 2 -1].priorityValue)
	{
		copyHeapElNUMT(&tempEl, &heapName->heapArray[heapIndex]);
		copyHeapElNUMT(&(heapName->heapArray[heapIndex]),  &(heapName->heapArray[(heapIndex + 1) / 2 -1]));
		copyHeapElNUMT(&(heapName->heapArray[(heapIndex + 1) / 2 -1]), &tempEl);
		heapIndex = (heapIndex + 1) / 2 - 1;
	}
}
