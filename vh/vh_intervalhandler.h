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
	char *chroName;
	int start;
	int end;
	char strand[2];
	char *type; // AluY, L1Hs for MEIs, or "satellite" for satellites, nothing for gaps and dups, etc.
	int size; //only for g_chroTable
}Interval;

typedef struct mei
{
  char *chroName;
  int start;
  int end;
  char strand[2];
  char *subclass;
  char *superclass;
} Mei;

typedef struct meiIndex
{
  char *chroName;
  int start;
  int end;
} MeiIndex;

#define STRMAX 50

extern Interval *dups;
extern int dup_count;

extern Interval *gaps;
extern int gap_count;

extern Interval *repeats;
extern int repeat_count;
extern int g_maxRepeatLength;

extern Mei *g_meiTable;
extern int meiTypeSize;
extern int selectedMeiTypeSize;

extern MeiIndex *mInd;

int meiIntervalSearch(ref_genome*, char*, int, int);
char* meiIntervalSearch2( ref_genome* ref, char *chroName, int pos);
int meiIntervalSearch_Span( ref_genome* ref, char* chr_name, int read_start,int read_end);
int mei_filtering(ref_genome* ref);
void vh_writeMeiIndices(ref_genome* ref);
int vh_notInRepeat (DivetRow *);
int notInRepeat( char* chroName, int pos);
int notInRepeat2 (char* chroName, int locMapLeftEnd, int locMapRightEnd);
int minim (int x, int y);
void vh_readRepeatTable (char *, char *);
void vh_readGapTable (char *);
void vh_readDupsTable(char *);
int dups_filtering(char *, int, int, float);

#endif /* VH_VH_INTERVALHANDLER_H_ */
