
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

	char *chroName;
	char *rName;
	int locMapLeftEnd;
	int locMapLeftStart;
	char orientationLeft;
	int locMapRightStart;
	int locMapRightEnd;
	char orientationRight;
	char svType;
	int mQual1;
	int mQual2;
	double avgQual;
	double editDistance;
	double phredScore;
	int divetRowId;
	char* meiType; //Only for mobile elements 0 : Alu +; 1: Alu -; 2: L1 +; 3; L1 -; 4: SVA +; 5; SVA -
	char* mei_subclass;
	unsigned long ten_x_barcode; // Only for 10x genomics data

	struct LibraryInfo *libInfo;

	struct DivetRow *next;
} DivetRow;


struct DivetRow *createDivetRow (struct ReadName *hash[],
		char *readName,
		char *chroName,
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

struct DivetRow *vh_loadDivetRowFromString (struct ReadName *hash[], char *line, struct LibraryInfo *libInfo, int id);
void vh_freeDivets ();
struct DivetRow *vh_loadDivetFile (struct LibraryInfo *);
int read_Divet_bam( discordantReadMapping_Info *discordantReadPtr, parameters *params, ref_genome* ref, LibraryInfo * libInfo, int chr_index, int counterDivetRow);
int read_Divet_bam_softClip( softClip *ptrSoftClip, parameters *params, ref_genome* ref, LibraryInfo * libInfo, int chr_index, int read_len, int divet_row_count);
void vh_printDivet (struct DivetRow *);

#endif
