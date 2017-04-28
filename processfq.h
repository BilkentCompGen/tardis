#ifndef __PROCESSFQ
#define __PROCESSFQ

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>
#include <stdbool.h>
#include "common.h"

#define MEMUSE 2047483648
#define MEMSCALE 1.5

/* Maximum sequence/quality length */
#define MAX_SEQ 1000

typedef struct read
{
	char* qname;
	char* seq;
	char* qual;
	char empty;
} _read;


typedef struct _discordantReadMapping_Info // RR, FF, RF or FR (Ins>\delta)
{
	char *chroName;
	char *readName;
	int pos1;
	int pos1_End;// In most cases pos1_End=pos1+readLen (unless there is soft Clip);
	int pos2;
	int pos2_End;// In most cases pos2_End=pos2+readLen (unless there is soft Clip);
	char orient1;
	char orient2;
	int mQual1; // mapping qual1
	int mQual2;
	int editDistance;
	char svType;
	unsigned long ten_x_barcode; // Only for 10x genomics data
	struct _discordantReadMapping_Info *next;
} discordantReadMapping_Info;


typedef struct _discordantReadMappingMEI_Info
{
	char *chroName; // chroName
	char *readName;
	char *subclass;
	int pos; // the high quality end mapped
	int pos_End;
	int qual;
	char orient; //the high quality end mapped
	int MEI_Type; // 0: Alu +, 1: Alu -, 2: L1 +, 3: L1 -, 4: SVA +, 5: SVA -,
	char* MEI_subclass;
	unsigned long ten_x_barcode; // Only for 10x genomics data
	struct _discordantReadMappingMEI_Info *next;
}discordantReadMappingMEI_Info;

typedef struct posMapSoftClip{
	int posMap;
	char orient;
	int mapq; // the mapq calculate based on number of mappings for this softclip + the distance
	struct posMapSoftClip *next;
}posMapSoftClip;

typedef struct softClip
{
	char *readName;
	char *chroName;
	int pos;
	char orient;
	int qual;
	int op[5];
	int opl[5];
	int opCount;
	char *softClipString;
	struct softClip *next;

	int avgPhredQualSoftClip; // the average phred quality of the basepairs which have be cut (the soft clip part).
	int numSoftClipInConcordance; // The total number of reads which have a soft clip +-10 from this read. Too high or too low indicates something is wrong and not interesting.

	posMapSoftClip *ptrPosMapSoftClip; // a linked list of all the positions that the soft clipped part of the read maps in the reference genome (chroName:windowStart-windowEnd)
}softClip;

struct library_properties
{
	char* libname; /* id/name of the library */
	float frag_avg; /* average fragment size */
	float frag_std; /* fragment size standard deviation */
	int frag_med; /* median of the fragment sizes */
	int conc_min; /* min cutoff for concordants */
	int conc_max; /* max cutoff for concordants */
	char* fastq1; /* file name for the FASTQ file of the /1 reads */
	char* fastq2; /* file name for the FASTQ file of the /2 reads */
	char *divet; /* file name for the DIVET file post mrfast mapping */
	int read_length; /* length of reads for this library */
	int num_sequences; /* number of paired-end sequences for this library */

	softClip *listSoftClip;
	discordantReadMapping_Info *listRR_FF_Mapping;
	discordantReadMapping_Info *listFR_Mapping;
	discordantReadMapping_Info *listRF_Mapping;
	discordantReadMappingMEI_Info *listMEI_Mapping;
};


/* Function Prototypes */
void fastq_match( char*, char*, int, int);
int load_reads( gzFile, struct read**, int);
static int fastq_qname_comp( const void*, const void*);
void alloc_reads( struct read***, int);
void realloc_reads( struct read***, int, int);
void free_reads( struct read***, int);
void create_fastq_library( struct library_properties* in_lib, char* sample_name, char* bam_path, parameters* params);

#endif
