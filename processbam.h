#ifndef __PROCESSBAM
#define __PROCESSBAM

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "processfq.h"

/* Sample this many fragments to calculate avg/median/std per library */
#define SAMPLEFRAG 1000000 

/* Maximum sequence/quality length */
#define MAX_SEQ 1000

/* Gender of the bam file */
enum gender{ MALE, FEMALE};

typedef struct _bam_info
{
	int read_count; /* total number of reads in this library */
	short** read_depth; /* read depth */
	short* read_depth_per_chr; /* read depth */
	float mean;
	float mean_x;
	float mean_y;
	float mean_rd_per_gc[101]; /* GC percentages, i.e., GC[13]=323 means 343 windows have GC of 13% */

	htsFile* bam_file; /* file pointer to the BAM file */
	hts_idx_t* bam_file_index;
	hts_itr_t *iter;
	bam_hdr_t* bam_header;

	enum gender sample_gender; /* gender of the sample */
	char* sample_name; /* name of the sample, parsed from SM in the BAM header */
	bool contribution; /* does this individual contribute to the particular SV? - checked in setcover.c */
	int num_libraries; /* number of libraries, counted from the RG tags in the BAM header */
	struct library_properties** libraries; /* each library_properties struct holds statistical/other info */
} bam_info;

/* Function Prototypes */
void load_bam( bam_info* in_bam, char* path);
void create_fastq( bam_info* in_bam, char* bam_path, parameters *params);
void print_bam( bam_info* in_bam);
void print_libs( bam_info* in_bam);
int find_library_index( bam_info* in_bam, char* library_name);
int sufficient_fragments_sampled( int* fragments_sampled, int num_libraries);
void set_library_min_max( struct library_properties* in_lib);
void freeDiscordantReads(bam_info* in_bam, char* bam_path, parameters *params);

/* BAM Utility functions */
void get_sample_name( bam_info* in_bam, char* header_text);
void get_library_count( bam_info* in_bam, char* header_text);
void get_library_names( bam_info* in_bam, char* header_text);
char* convertUCSCtoGRC( char* inputchr);


#endif
