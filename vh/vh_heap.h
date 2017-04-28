#ifndef _VH_HEAP__H_
#define _VH_HEAP__H_

#define MAX_CLUSTER_SIZE 10000000	//TODO: Ye karish bokon

#include "vh_common.h"

typedef struct HeapEl
{
  //int readMappingId;
  struct DivetRow *readMappingPtr;
  int priorityValue;
} HeapEl;

typedef struct Heap
{
  struct HeapEl heapArray[MAX_CLUSTER_SIZE];
  int heapSize;
} Heap;

extern Heap *g_intersectInterval;

void vh_copyHeapEl (struct HeapEl *dest, struct HeapEl *src);
void vh_push_heap (struct Heap *heapName, struct HeapEl *newEl);
int vh_minValue_heap (struct Heap *heapName);
void vh_heapBubleDown (struct Heap *heapName);
void vh_heap_remove_top (struct Heap *heapName);
void vh_addToHeap (struct DivetRow *readMappingPtr, int priorityValue,struct Heap *heapName);
void vh_writeHeap (struct Heap *heapName);

#endif
