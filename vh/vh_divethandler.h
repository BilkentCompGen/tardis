#ifndef VH_DIVET_HANDLER__
#define VH_DIVET_HANDLER__

#include <stdlib.h>
#include <string.h>
#include "vh_common.h"
#include "vh_hash.h"
#include "../processfq.h"
#include "../processbam.h"

#define DIVET_ROW_DELIMITERS " \t\r\n"
#define BUFFER_SIZE 20000

#define MAX_LIBNAME_LEN  512
#define MAX_FILENAME_LEN 512

typedef struct LibraryInfo
{
	char libName[MAX_LIBNAME_LEN];
	char indName[MAX_LIBNAME_LEN];
	char libFileAdrs[MAX_FILENAME_LEN];
	int libId;
	int minDelta;
	int maxDelta;
	int readLen;
	int size;

	struct ReadName **hash;
	struct DivetRow *head;
	struct DivetRow *tail;

	struct LibraryInfo *next;
} LibraryInfo;


typedef struct DivetRow
{
	struct ReadName *readName;

	char *chromosome_name;
	int locMapLeftEnd;
	int locMapLeftStart;
	char orientationLeft;
	int locMapRightStart;
	int locMapRightEnd;
	int startPosition;
	int endPosition;
	char orientationRight;
	char svType;
	char splitOrientation;
	int mQual1;
	int mQual2;
	double avgQual;
	double editDistance;
	double phredScore;
	int divetRowId;
	char* meiType;
	int mei_code; //Only for mobile elements 0 : Alu +; 1: Alu -; 2: L1 +; 3; L1 -; 4: SVA +; 5; SVA -
	char* mei_subclass;
	unsigned long ten_x_barcode; // Only for 10x genomics data
	int in_cluster_count;

	struct LibraryInfo *libInfo;

	struct DivetRow *next;
} DivetRow;


struct DivetRow *createDivetRow (struct ReadName *hash[],
		char *readName,
		char *chromosome_name,
		char *locMapLeftStart,
		char *locMapLeftEnd,
		char *orientatinoLeft,
		char *locMapRightStart,
		char *locMapRightEnd,
		char *orientationRight,
		char *svType,
		char *editDistance,
		char *avgQual,
		char *phredScore,
		unsigned long ten_x_barcode,
		struct LibraryInfo *libInfo,
		int id);

typedef struct discordantMapping // RR, FF, RF or FR (Ins>\delta)
{
	char *chromosome_name;
	char *readName;
	int pos1;
	int pos1_End;// In most cases pos1_End=pos1+readLen (unless there is soft Clip);
	int pos2;
	int pos2_End;// In most cases pos2_End=pos2+readLen (unless there is soft Clip);
	char orient1;
	char orient2;
	int mQual1; // mapping qual1
	int mQual2;
	int editDistance_left;
	int editDistance;
	char svType;
	uint32_t* cigar;
	uint32_t n_cigar;
	uint16_t flag;
	unsigned long ten_x_barcode; // Only for 10x genomics data

	struct discordantMapping *next;
}discordantMapping;

typedef struct alternativeMapping // RR, FF, RF or FR (Ins>\delta)
{
	char *chromosome_name;
	char *readName;
	int pos1;
	int pos1_End;// In most cases pos1_End=pos1+readLen (unless there is soft Clip);
	int pos2;
	int pos2_End;// In most cases pos2_End=pos2+readLen (unless there is soft Clip);
	int mQual1; // mapping qual1
	int mQual2;
	int editDistance_left;
	int editDistance;
	int mapp;
	uint16_t flag;
	unsigned long ten_x_barcode; // Only for 10x genomics data

	struct alternativeMapping *next;
}alternativeMapping;


typedef struct discordantMappingNUMT
{
	char *chromosome_name;
	char *readName;
	int pos;
	int pos_End;
	int qual;
	char orientation; /*orientation of the pair (the one not inside the MT) */
	char NUMT_Type; /* numt in forward direction is 0, reverse direction is 1 */
	unsigned long ten_x_barcode; // Only for 10x genomics data
	struct discordantMappingNUMT *next;
}discordantMappingNUMT
;

typedef struct discordantMappingMEI
{
	char *chromosome_name; // chromosome_name
	char *readName;
	int pos; // the high quality end mapped
	int pos_End;
	int qual;
	char orient; //the high quality end mapped
	int MEI_Type; // 0: Alu +, 1: Alu -, 2: L1 +, 3: L1 -, 4: SVA +, 5: SVA -,
	char *MEI_class;
	char* MEI_subclass;
	unsigned long ten_x_barcode; // Only for 10x genomics data
	struct discordantMappingMEI *next;
}discordantMappingMEI;


typedef struct posMapSoftClip{
	int posMap;
	char orient;
	int mapq; // the mapq calculate based on number of mappings for this softclip + the distance

	struct posMapSoftClip *next;
}posMapSoftClip;


typedef struct softClip
{
	char *readName;
	char *chromosome_name;
	int pos;
	char orient;
	int qual;
	int op[5];
	int opl[5];
	int opCount;
	char *softClipString;
	int avgPhredQualSoftClip; // the average phred quality of the basepairs which have be cut (the soft clip part).
	int numSoftClipInConcordance; // The total number of reads which have a soft clip +-10 from this read. Too high or too low indicates something is wrong and not interesting.

	struct softClip *next;
	posMapSoftClip *ptrPosMapSoftClip; // a linked list of all the positions that the soft clipped part of the read maps in the reference genome (chromosome_name:windowStart-windowEnd)
}softClip;

void vh_pruneAndNormalizeDivets( struct LibraryInfo *lib, double preProsPrune, int overMapLimit);
int load_Divet_bam( bam_info** in_bams, parameters *params, int chr_index);
void vh_freeDivets ();
struct DivetRow *vh_loadDivetFile (struct LibraryInfo *, sonic *);
void vh_printDivet (struct DivetRow *);

#endif
