#include "vh_common.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_heap.h"
#include "vh_maximalCluster.h"


void vh_addToGenomeIndex_TDup (char *chromosome_name, sonic *this_sonic)
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
	  if (strcmp (divetReadMappingPtr->chromosome_name, chromosome_name) == 0
	      && divetReadMappingPtr->svType == 'E'
	      && !sonic_is_gap (this_sonic, chromosome_name, divetReadMappingPtr->locMapLeftEnd,
				divetReadMappingPtr->locMapRightStart)
	      && (divetReadMappingPtr->locMapRightStart -
		  divetReadMappingPtr->locMapLeftEnd < maxDuplicationLen))
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

void vh_initializeReadMapping_TDup (char *chromosome_name, int chroSize, sonic *this_sonic)
{
  int i;
  LibraryInfo *libInfoPtr = g_libInfo;

  /* Initializing the Genome Array */
  g_genomeIndexStart = (MappingOnGenome **) getMem (chroSize * sizeof (MappingOnGenome *));

  if (g_genomeIndexStart == NULL)
    vh_logWarning ("Memory Problem in vh_createMaxClusterDeletion.cpp:53");
  for (i = 0; i < chroSize; i++)
    {
      g_genomeIndexStart[i] = NULL;
    }

  /* Initializing the list of begin and end of right side break point ranges */
  g_listRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
  g_tempListRightBrkPointIntr = (RightBrkPointInterval *) getMem (g_maxListBrkPointIntr * sizeof (RightBrkPointInterval));
  for (i = 0; i < g_maxListBrkPointIntr; i++)
    {
      g_listRightBrkPointIntr[i].readMappingPtr = NULL;
      g_tempListRightBrkPointIntr[i].readMappingPtr = NULL;
    }
  g_listRightBrkPointIntrCount = 0;
  vh_addToGenomeIndex_TDup (chromosome_name, this_sonic);


  g_intersectInterval = (Heap *) getMem (sizeof (Heap));
  g_intersectInterval->heapSize = 0;

  while (libInfoPtr != NULL)
    {
      g_maxDeltaAmongLibs = vh_max (g_maxDeltaAmongLibs, libInfoPtr->maxDelta);
      libInfoPtr = libInfoPtr->next;
    }

}

void vh_reevaluate_TDup (int id, int brkPointLeft)
{
  int brkRightTemp, brkLeftTemp;
  brkRightTemp = vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
			 vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->minDelta -
				 (brkPointLeft -g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd), 0),
			 g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd);
  brkLeftTemp = vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
			vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->
				libInfo->maxDelta - (brkPointLeft - g_tempListRightBrkPointIntr[id]. readMappingPtr->locMapLeftEnd), 0),
			g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd);

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

int vh_createBreakPointIntervals_TDup (int brkPointLeft)
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
      //TODO: Can this be made more efficient using a Heap?
      while (listRightBrkPointIntrId < g_listRightBrkPointIntrCount)
	{
	  if (g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapLeftEnd +
	      g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->libInfo->maxDelta == brkPointLeft)
	    {
	      listRightBrkPointIntrId++;
	    }
	  else
	    {
	      vh_copyElBrkPointIntr (tempListRightBrkPointIntrId,listRightBrkPointIntrId); // increaseByOneRightBrkPointIntr(tempListRightBrkPointIntrId);
	      vh_reevaluate_TDup (tempListRightBrkPointIntrId, brkPointLeft);
	      listRightBrkPointIntrId++;
	      tempListRightBrkPointIntrId++;
	    }

	}
    }

  ptrMappingOnGenome = g_genomeIndexStart[brkPointLeft];
  while (ptrMappingOnGenome != NULL)
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
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = 'L';
      tempListRightBrkPointIntrId++;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointRight = locBrkPointRightTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft= locBrkPointLeftTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key = locBrkPointRightTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr= ptrMappingOnGenome->readMappingPtr;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = 'R';
      tempListRightBrkPointIntrId++;
      ptrMappingOnGenome = ptrMappingOnGenome->next;
    }
  temp = g_listRightBrkPointIntr;
  g_listRightBrkPointIntr = g_tempListRightBrkPointIntr;
  g_tempListRightBrkPointIntr = temp;
  g_listRightBrkPointIntrCount = tempListRightBrkPointIntrId;
  return (newElAdded);
}

// This is the head function called for creating clusters.
void vh_createTDupClusters (int chroSize)
{
  int leftBreakPoint;
  int startLeftWindow, endLeftWindow;
  int newElAdded = 0;
  g_listRightBrkPointIntrCount = 0;

  for (leftBreakPoint = 1; leftBreakPoint < chroSize; leftBreakPoint++)
    {
      newElAdded = 0;
      newElAdded = vh_createBreakPointIntervals_TDup (leftBreakPoint);
      if (newElAdded)		// For deletion only when we have added new element we need to check
	vh_createIntersectingIntervals (leftBreakPoint, 4);
    }
  vh_flushOut (g_listPotClusterFound, leftBreakPoint,4);
}

