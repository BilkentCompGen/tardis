#ifndef __COMMON
#define __COMMON

#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>
#include <stdbool.h>
#include "sonic/sonic.h"

/* Exit Codes */
#define EXIT_SUCCESS 0
#define EXIT_COMMON 1
#define EXIT_MAXBAMS 2
#define EXIT_PARAM_ERROR 3
#define EXIT_EXTERNAL_PROG_ERROR 4
#define EXIT_FILE_OPEN_ERROR 5
#define EXIT_READGROUP 6
#define EXIT_SONIC 7

/* Return Codes */
#define RETURN_SUCCESS 1
#define RETURN_ERROR 0

#define MAX_BAMS 256

/* Maximum filename length */
#define MAX_LENGTH 1024

/* MAPPING INFO */
#define RPUNMAPPED 0
#define RPCONC 1
#define RPDEL 2
#define RPINV 3
#define RPINS 4
#define RPTDUP 5
#define RPMEI 6
#define RPINTERCHR 7;

/* Running mode of tardis */
#define QUICK 0
#define SENSITIVE 1
extern int running_mode;
extern int ten_x_flag;
extern int output_hs_flag;


// Track memory usage
extern long long memUsage;
extern FILE *logFile; //Defined in tardis.c

typedef struct _params
{
	char* ref_genome; /* path to reference genome - fasta */
	char* reps; /* path to repeatmasker file - *rm.out */
	char* dups; /* path to segmental duplications file - bed */
	char* bam_files; /* paths to comma separated input BAM files as a single string before being tokenized */
	char* bam_list_path; /* path to a file that lists BAM file paths in advance */
	char** bam_file_list; /* the actual list that holds all bam file paths after tokenization */
	char* gaps; /* path to assembly gaps file - bed */
	char* mei;  /* regular expression-like MEI list */
	char* outprefix; /* prefix for the output files */
	int  force_read_length; /* force read length to a certain value, discard those that are shorter. Hidden feature due to GIAB */
	char run_vh; /* boolean stand-in to run VariationHunter */
	char run_rd; /* boolean stand-in to run Read Depth */
	char run_ns; /* boolean stand-in to run NovelSeq */
	char run_sr; /* boolean stand-in to run SPLITREAD */
	char no_soft_clip; /* boolean stand-in to skip soft clip */
	char skip_fastq; /* boolean stand-in to skip FASTQ dump */
	char skip_sort; /* boolean stand-in to skip FASTQ sort */
	char skip_remap; /* boolean stand-in to skip FASTQ remap */
	char skip_vhcluster; /* boolean stand-in to skip VH clustering */
	int  threads; /* number of threads to use for parallel mrFAST, and maybe future parallelization of TARDIS */
	int num_bams; /* number of input BAM files */
	int quick; /* boolean stand-in to work in bam-only mode (no divet) */
	int sensitive; /* boolean stand-in to work in sensitive mode (divet) */
	int ten_x; /*boolean for whether we're using 10x data*/
	int output_hs; /*boolean for whether to record the homogeneity score (HS) in VCF regardless whether HS is used in set cover or not*/
	int make_sonic; /*make SONIC file and exit*/
	int load_sonic; /*load SONIC file*/
        char *sonic_info; /* SONIC reference information string for building */
	int first_chrom; /*the first chromosome as indexed in the ref to be computer for. 0 by default*/
	int last_chrom; /*the last chromosome as indexed in the ref to be computer for. ref->chrom_count by default*/
	int rd_threshold; /* Threshold is used in RD filtering, calWeight() in vh_setcover.c */
	int mq_threshold; /* Minimum mapping quality */
	int rp_threshold; /* Minimum read pair support */
        char *sonic_file; /* SONIC file name */
	sonic *this_sonic; /* SONIC */
} parameters;


typedef struct _ref_genome
{
	char* ref_name; /* name of the chromosome */
        bool* in_bam; /* if the chromosome is available in the bam file */
  /* the following don't seem to be necessary any more. load_refgen now just copies the pointer from sonic to avoid substantial refactoring */
        int chrom_count; /* number of chromosomes */
        int* chrom_lengths; /* lengths of the chromosomes */
        char** chrom_names; /* names of the chromosomes */
}ref_genome;

/* Parameter related TARDIS functions */
void init_params( parameters**);
void print_params( parameters*);

/* FILE opening and error printing functions. For opening regular and BAM/SAM
 files safely */
void print_error( char*);
FILE* safe_fopen( char* path, char* mode);
gzFile safe_fopen_gz( char* path, char* mode);
htsFile* safe_hts_open( char* path, char* mode);

/* General BAM processing functions */
int is_proper( int flag);
int is_discordant( bam1_core_t bam_alignment_core, int min, int max);
int is_concordant( bam1_core_t bam_alignment_core, int min, int max);
char base_as_char( int base_as_int);
char complement_char( char base);
void qual_to_ascii( char* qual);

/* String functions */
void set_str( char **target, char *source); /* Even safer than strncpy */
void reverse_string( char* str);

/* Misc. Utility */
int compare_size_int( const void* p, const void* q);
void print_quote( void);
int find_chr_index_bam(ref_genome*, char*, bam_hdr_t*);
int max( int x, int y);
int min( int x, int y);
int hammingDistance( char *str1, char *str2, int len);

// Memory allocation/tracking functions
void* getMem( size_t size);
void freeMem( void* ptr, size_t size);
double getMemUsage();


unsigned long encode_ten_x_barcode(char* source);

#endif
