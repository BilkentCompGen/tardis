#ifndef _VH_HEAP__H_
#define _VH_HEAP__H_

#define MAX_CLUSTER_SIZE 50000


#include "vh_common.h"
#include "vh_createMaxClusterMEI.h"
#include "vh_createMaxClusterNUMT.h"


typedef struct PEAlistEl {
        int brkPoint;
        //int priorityValue;
        struct DivetRow *readMappingPtr;
        struct PEAlistEl *next;
        struct PEAlistEl *prev;
} PEAlistEls;

typedef struct HeapEl
{
  struct DivetRow *readMappingPtr;
  int priorityValue;
} HeapEl;

typedef struct Heap
{
        struct HeapEl *heapArray;
        int heapSize;
        int current_max_size;
} Heap;

typedef struct HeapElMEI{
	int priorityValue;
	mei_Reads *mei_ReadsPtr;
}HeapElMEI;


typedef struct HeapMEI{
        HeapElMEI *heapArray;
	int heapSize;
        int current_max_size;
} HeapMEI;

typedef struct HeapElNUMT{
	int priorityValue;
	numt_Reads *numt_ReadsPtr;
}HeapElNUMT;


typedef struct HeapNUMT{
        HeapElNUMT *heapArray;
	int heapSize;
        int current_max_size;
} HeapNUMT;


extern Heap *g_intersectInterval;

Heap *vh_newHeap(void);
HeapMEI *vh_newHeapMEI(void);
HeapNUMT *vh_newHeapNUMT(void);
void vh_renewHeap(Heap *);
void vh_renewHeapMEI(HeapMEI *);
void vh_renewHeapNUMT(HeapNUMT *);
  
void vh_free_heap(Heap *);
void vh_free_heap_mei(HeapMEI *);
void vh_free_heap_numt(HeapNUMT *);
void vh_copyHeapEl (struct HeapEl *dest, struct HeapEl *src);
void vh_push_heap( struct Heap *heapName, struct HeapEl *newEl);
int vh_minValue_heap( struct Heap *heapName);
void vh_heapBubleDown( struct Heap *heapName);
void vh_heap_remove_top( struct Heap *heapName);
void vh_addToHeap( struct DivetRow *readMappingPtr, int priorityValue, struct Heap *heapName);
void vh_writeHeap( struct Heap *heapName);
void heapBubleDownMEI( HeapMEI *heapName);
void heapBubleDownNUMT(HeapNUMT *heapName);
void heap_remove_topMEI( HeapMEI *heapName);
void heap_remove_topNUMT( HeapNUMT *heapName);
void push_heap_mei( HeapMEI *heapName, HeapElMEI *newEl);
void push_heap_numt(HeapNUMT *heapName, HeapElNUMT *newEl);
void freeHeapMEI( HeapMEI *heapName);
int minValue_heapMEI(HeapMEI *heapName);
int minValue_heapNUMT(HeapNUMT *heapName);
int remove_Heap(Heap * heapName, int position);

#endif
