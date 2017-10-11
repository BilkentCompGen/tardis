#ifndef __PROCESSFQ
#define __PROCESSFQ

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>
#include <stdbool.h>
#include "common.h"
#include "vh/vh_divethandler.h"

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

typedef struct library_properties
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

	struct softClip *listSoftClip;
	struct discordantMapping *listRR_FF_Mapping;
	struct discordantMapping *listFR_Mapping;
	struct discordantMapping *listRF_Mapping;
	struct discordantMappingMEI *listMEI_Mapping;
} library_properties;


/* Function Prototypes */
void fastq_match( char*, char*, int, int);
int load_reads( gzFile, struct read**, int);
static int fastq_qname_comp( const void*, const void*);
void alloc_reads( struct read***, int);
void realloc_reads( struct read***, int, int);
void free_reads( struct read***, int);
void create_fastq_library( struct library_properties* in_lib, char* sample_name, char* bam_path, parameters* params);

#endif
