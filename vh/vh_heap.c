#include "vh_heap.h"
#include "vh_divethandler.h"

int resize_cnt = 0;

Heap *g_intersectInterval;

Heap *vh_newHeap(int size)
{
       Heap *tmp_heap;
       tmp_heap = (Heap *) getMem (sizeof (Heap));
       tmp_heap->heapArray = (struct HeapEl *) getMem(sizeof(struct HeapEl) * size); 
       tmp_heap->heapSize = 0;
       tmp_heap->current_max_size = MAX_CLUSTER_SIZE;
       return tmp_heap;
}

HeapMEI *vh_newHeapMEI(int size)
{
       HeapMEI *tmp_heap;
       tmp_heap = (HeapMEI *) getMem (sizeof (HeapMEI));
       tmp_heap->heapArray = (struct HeapElMEI *) getMem(sizeof(struct HeapElMEI) * size); 
       tmp_heap->heapSize = 0;
       tmp_heap->current_max_size = size;
       return tmp_heap;
}

HeapNUMT *vh_newHeapNUMT(int size)
{
       HeapNUMT *tmp_heap;
       tmp_heap = (HeapNUMT *) getMem (sizeof (HeapNUMT));
       tmp_heap->heapArray = (struct HeapElNUMT *) getMem(sizeof(struct HeapElNUMT) * size); 
       tmp_heap->heapSize = 0;
       tmp_heap->current_max_size = MAX_CLUSTER_SIZE;
       return tmp_heap;
}


/*
  DEPRECATED FEB 25, 2019. WILL BE REMOVED FROM CODE BASE IN THE FUTURE

Heap *vh_renewHeap(Heap *h)
{
       Heap *tmp_heap;
       int new_size = h->current_max_size + MAX_CLUSTER_SIZE;
       tmp_heap = (Heap *) getMem (sizeof (Heap));
       tmp_heap->heapArray = (struct HeapEl *) getMem(sizeof(struct HeapEl) * new_size);
       tmp_heap->heapSize = h->heapSize;
       tmp_heap->current_max_size = new_size;
       memcpy(tmp_heap->heapArray, h->heapArray, (h->heapSize * sizeof(struct HeapEl)));
       vh_free_heap(h);
       return tmp_heap;
}

HeapMEI *vh_renewHeapMEI(HeapMEI *h)
{
       HeapMEI *tmp_heap;
       int i;
       int new_size = h->current_max_size + MAX_CLUSTER_SIZE;
       printf ("resizing heapMEI from %d to %d.\n", h->current_max_size, new_size);
       tmp_heap = vh_newHeapMEI(new_size);
       //	 (HeapMEI *) getMem (sizeof (HeapMEI));
       //       tmp_heap->heapArray = (struct HeapElMEI *) getMem(sizeof(struct HeapElMEI) * new_size);
       if (tmp_heap->heapArray == NULL)
	 fprintf(stderr, "cannot allocate memory for heapArray, size: %d\n", sizeof(struct HeapElMEI) * new_size);
       tmp_heap->heapSize = h->heapSize;
       //       tmp_heap->current_max_size = new_size;
       /*
	 memcpy(tmp_heap->heapArray, h->heapArray, (h->heapSize * sizeof(struct HeapElMEI)));
       
       for (i = 0; i < h->current_max_size; i++){
	 tmp_heap->heapArray[i].priorityValue = h->heapArray[i].priorityValue;
	 tmp_heap->heapArray[i].mei_ReadsPtr = h->heapArray[i].mei_ReadsPtr;
       }
       //vh_free_heap_mei(h);
       fprintf(stderr, "resize %d done. Current max size: %d\n", resize_cnt, tmp_heap->current_max_size);
       resize_cnt++;
       return tmp_heap;
}

HeapNUMT *vh_renewHeapNUMT(HeapNUMT *h)
{
       HeapNUMT *tmp_heap;
       int new_size = h->current_max_size + MAX_CLUSTER_SIZE;
       tmp_heap = (HeapNUMT *) getMem (sizeof (HeapNUMT));
       tmp_heap->heapArray = (struct HeapElNUMT *) getMem(sizeof(struct HeapElNUMT) * new_size); 
       tmp_heap->heapSize = h->heapSize;
       tmp_heap->current_max_size = new_size;
       memcpy(tmp_heap->heapArray, h->heapArray, (h->heapSize * sizeof(struct HeapElNUMT)));
       vh_free_heap_numt(h);
       return tmp_heap;
}

*/


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
	int size;
	
	/* HEAP_MEMORY_REALLOCATION */
	if (heapName->heapSize == heapName->current_max_size){
	  size = heapName->current_max_size + MAX_CLUSTER_SIZE;
	  heapName->heapArray = (struct HeapEl *) reallocMem(heapName->heapArray, sizeof(struct HeapEl) * heapName->heapSize, sizeof(struct HeapEl) * size);
	  heapName->current_max_size = size;
	}

	

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

	HeapElMEI tempEl;
        int heapIndex = heapName->heapSize;
	int size;
	
	/* HEAP_MEMORY_REALLOCATION  */
	if (heapName->heapSize == heapName->current_max_size){
	  size = heapName->current_max_size + MAX_CLUSTER_SIZE;
	  heapName->heapArray = (struct HeapElMEI *) reallocMem(heapName->heapArray, sizeof(struct HeapElMEI) * heapName->heapSize, sizeof(struct HeapElMEI) * size);
	  heapName->current_max_size = size;
	}
	
	
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
	int size;
	
	/* HEAP_MEMORY_REALLOCATION */
	if (heapName->heapSize == heapName->current_max_size){
	  size = heapName->current_max_size + MAX_CLUSTER_SIZE;
	  heapName->heapArray = (struct HeapElNUMT *) reallocMem(heapName->heapArray, sizeof(struct HeapElNUMT) * heapName->heapSize, sizeof(struct HeapElNUMT) * size);
	  heapName->current_max_size = size;
	}

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

int *resize_int_array (int *old_array, int old_size, int new_size){
  int *ret;
  ret = (int *) reallocMem(old_array, old_size * sizeof(int), new_size * sizeof(int));
  /*
  memcpy(ret, old_array, old_size * sizeof(int));
  freeMem (old_array, old_size * sizeof(int));
  memset(ret + (old_size*sizeof(int)), 0, (new_size-old_size) * sizeof(int));
  */
  return ret;
}
