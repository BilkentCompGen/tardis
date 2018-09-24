#ifndef __COMMON
#define __COMMON

#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>
#include <stdbool.h>
#include "sonic/sonic.h"

//#define MAIN_DELETION_CLUSTER
#define INVERSION 'V'
#define INSERTION 'I'
#define DELETION 'D'
#define TANDEMDUP 'E'

#define NUMTFORWARD 'N'
#define NUMTREVERSE 'O'

#define INVDUP 'W'
#define INVDUPLEFT 'S'
#define INVDUPRIGHT 'T'

#define INTERDUP 'X'
#define INTERDUPLEFT 'H'
#define INTERDUPRIGHT 'M'

#define MEIFORWARD 'A'
#define MEIREVERSE 'B'

/* For split reads FFAB means reads are in Forward-Forward orientation and
 * position of A is smaller than position of B */
#define FFAB 'C'
#define FFBA 'G'
#define FRAB 'J'
#define FRBA 'K'
#define RFAB 'L'
#define RFBA 'N'


#define LEFT 'L'
#define RIGHT 'R'
#define NONE 'N'
#define FORWARD 'F'
#define REVERSE 'R'

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
#define MAXLISTBRKPOINTINTR 10000000;

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
#define RPNUMT 8;

/* Running mode of tardis */
#define QUICK 0
#define SENSITIVE 1
extern int running_mode;
extern int ten_x_flag;
extern int output_hs_flag;
extern int debug_mode; /* boolean stand-in to work in debug mode - .name and .clusters are created */
extern int cluster_of_reads;

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
	char* outdir; /* output directory */
	int force_read_length; /* force read length to a certain value, discard those that are shorter. Hidden feature due to GIAB */
	char no_soft_clip; /* boolean stand-in to skip soft clip */
	int alt_mapping; /* check the alternative mapping locations from the xa field in bwa */
	int seq_resolved; /* whether to output sequence resolved calls */
	int no_interdup; /* whether to cluster interspersed duplications */
	char skip_mrfast; /* boolean stand-in to skip mrFast mapping (If you already have the divet file created) */
	int threads; /* number of threads to use for parallel mrFAST, and maybe future parallelization of TARDIS */
	int num_bams; /* number of input BAM files */
	int quick; /* boolean stand-in to work in bam-only mode (no divet) */
	int sensitive; /* boolean stand-in to work in sensitive mode (divet) */
	int ten_x; /*boolean for whether we're using 10x data*/
	int output_hs; /*boolean for whether to record the homogeneity score (HS) in VCF regardless whether HS is used in set cover or not*/
	int make_sonic; /*make SONIC file and exit*/
	int load_sonic; /*load SONIC file*/
	char *sonic_info; /* SONIC reference information string for building */
	int first_chr; /*the first chromosome as indexed in sonic. 0 by default*/
	int last_chr; /*the last chromosome as indexed in sonic. chrom_count by default*/
	int mq_threshold; /* Minimum mapping quality */
	int rp_threshold; /* Minimum read pair support */
	int number_of_different_mei_types; /* Number of distinct MEI types e.g. ALU:L1:SVA has three different types */
	int cluster_of_read; /* Number of clusters that a read can be involved. 10 by default */
	char *sonic_file; /* SONIC file name */
	sonic *this_sonic; /* SONIC */
} parameters;

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
int find_chr_index_bam( char* chromosome_name, bam_hdr_t* bam_header);
int max( int x, int y);
int min( int x, int y);
int hammingDistance( char *str1, char *str2, int len);
int vh_cmprReadNameStr (const void *a, const void *b);
int32_t calculateInsertSize( int32_t pos_left, int32_t pos_right,uint16_t flag, int read_length);

// Memory allocation/tracking functions
void* getMem( size_t size);
void freeMem( void* ptr, size_t size);
double getMemUsage();


unsigned long encode_ten_x_barcode(char* source);

void get_working_directory(parameters *);
void clean_up_temp_files(parameters *);

#endif
