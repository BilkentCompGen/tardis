#ifndef _VH_HEAP__H_
#define _VH_HEAP__H_

#define MAX_CLUSTER_SIZE 10000000	//TODO: Ye karish bokon

#include "vh_common.h"

typedef struct HeapEl
{
  struct DivetRow *readMappingPtr;
  int priorityValue;
} HeapEl;

typedef struct Heap
{
  struct HeapEl heapArray[MAX_CLUSTER_SIZE];
  int heapSize;
} Heap;


typedef struct mei_Reads{
	char *readName;
	int pos;
	int pos_End;
	char orient;
	int mQual;
	int MEI_Type; // 0 : Alu +; 1: Alu -; 2: L1 +; 3; L1 -; 4: SVA +; 5; SVA -
	char* MEI_Class;
	char* MEI_Subclass;
	char *libName;
	char *indName;
	int libId;
	int readTypeSupport;// 0: ReadPair, 1: SplitRead (Soft clipping of the beginig of the read, 2: SplitRead (Soft clipping at the end of the read)

	struct mei_Reads *next;
} mei_Reads;


typedef struct HeapElMEI{
	int priorityValue;
	mei_Reads *mei_ReadsPtr;
}HeapElMEI;


typedef struct HeapMEI{
	HeapElMEI heapArray[MAX_CLUSTER_SIZE];
	int heapSize;
} HeapMEI;


extern Heap *g_intersectInterval;

void vh_copyHeapEl (struct HeapEl *dest, struct HeapEl *src);
void vh_push_heap( struct Heap *heapName, struct HeapEl *newEl);
int vh_minValue_heap( struct Heap *heapName);
void vh_heapBubleDown( struct Heap *heapName);
void vh_heap_remove_top( struct Heap *heapName);
void vh_addToHeap( struct DivetRow *readMappingPtr, int priorityValue, struct Heap *heapName);
void vh_writeHeap( struct Heap *heapName);
void heapBubleDownMEI( HeapMEI *heapName);
void heap_remove_topMEI( HeapMEI *heapName);
void push_heap_mei( HeapMEI *heapName, HeapElMEI *newEl);
void freeHeapMEI( HeapMEI *heapName);
int minValue_heapMEI(HeapMEI *heapName);

#endif
