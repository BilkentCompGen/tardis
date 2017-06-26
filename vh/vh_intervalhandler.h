/*
 * vh_intervalhandler.h
 *
 *  Created on: Nov 1, 2016
 *      Author: tardis
 */

#ifndef VH_VH_INTERVALHANDLER_H_
#define VH_VH_INTERVALHANDLER_H_

#include "vh_divethandler.h"

typedef struct interval{
	char *chromosome_name;
	int start;
	int end;
	char strand[2];
	//char *type; // AluY, L1Hs for MEIs, or "satellite" for satellites, nothing for gaps and dups, etc.
	int size; //only for g_chroTable
}Interval;

typedef struct mei
{
  char *chromosome_name;
  int start;
  int end;
  char strand[2];
  char *subclass;
  char *superclass;
} Mei;

typedef struct meiIndex
{
  char *chromosome_name;
  int start;
  int end;
} MeiIndex;

#define STRMAX 50

extern Interval *dups;
extern int dup_count;

/*
extern Interval *gaps;
extern int gap_count;
*/

extern Interval *repeats;
extern int repeat_count;
extern int g_maxRepeatLength;

extern Mei *g_meiTable;
extern int meiTypeSize;
extern int selectedMeiTypeSize;

extern MeiIndex *mInd;


int mei_filtering(ref_genome* ref, parameters *params);

#endif /* VH_VH_INTERVALHANDLER_H_ */
