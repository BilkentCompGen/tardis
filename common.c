#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "common.h"
#include "quotes.h"

// Track memory usage
long long memUsage = 0;
int running_mode = QUICK;
int ten_x_flag = 0;
int output_hs_flag = 0;
int debug_mode = 0;

void init_params( parameters** params)
{
	int i;

	/* initialize parameters */
	*params = ( parameters*) getMem( sizeof( parameters));
	( *params)->ref_genome = NULL;
	( *params)->sonic_file = NULL;
	( *params)->this_sonic = NULL;
	( *params)->reps = NULL;
	( *params)->dups = NULL;
	( *params)->bam_files = NULL;
	( *params)->bam_list_path = NULL;
	( *params)->outprefix = NULL;
	( *params)->bam_file_list = ( char**) getMem( sizeof( char*) * MAX_BAMS);
	( *params)->gaps = NULL;
	( *params)->mei = NULL;
	( *params)->force_read_length = 0;
	( *params)->run_vh = 0; 
	( *params)->run_rd = 0;
	( *params)->run_ns = 0;
	( *params)->run_sr = 0;
	( *params)->threads = 1;
	( *params)->num_bams = 0;
	( *params)->skip_fastq = 0;
	( *params)->skip_sort = 0;
	( *params)->skip_remap = 0;
	( *params)->ten_x = 0;
	( *params)->output_hs = 0;
	( *params)->first_chrom = 0;
	( *params)->last_chrom = -1;
	( *params)->make_sonic = 0;
	( *params)->load_sonic = 0;
	( *params)->sonic_info = NULL;
	( *params)->number_of_different_mei_types = 0;
	( *params)->quick = 1;

	for( i = 0; i < MAX_BAMS; i++)
	{
		(*params)->bam_file_list[i] = NULL;
	}
}


void print_params( parameters* params)
{
	int i;
	printf("\n");
	for( i = 0; i < params->num_bams; i++)
	{
		printf( "%-30s%s\n","BAM input:",params->bam_file_list[i]);
		fprintf( logFile,"%-30s%s\n","BAM input:",params->bam_file_list[i]);
	}
	printf( "%-30s%s\n","Reference genome:", params->ref_genome);
	printf( "%-30s%s\n","SONIC file:", params->sonic_file);
	printf( "%-30s%s\n","Mobile Elements:", params->mei);

	fprintf( logFile, "%-30s%s\n","Reference genome:", params->ref_genome);
	fprintf( logFile, "%-30s%s\n","SONIC file:", params->sonic_file);
	fprintf( logFile, "%-30s%s\n","Mobile Elements:", params->mei);
	fprintf( logFile, "%-30s%d\n","10x Tag:", params->ten_x);
	fprintf( logFile, "%-30s%d\n","First chrom:", params->first_chrom);
	fprintf( logFile, "%-30s%d\n","Last chrom:", params->last_chrom);
}

void print_error( char* msg)
{
	/* print error message and exit */
	fprintf( stderr, "\n%s\n", msg);
	fprintf( stderr, "Invoke parameter -h for help.\n");
	exit( EXIT_COMMON);
}

void print_quote( void)
{
	/* print a quote by the Doctor */

	int quotenum;

	srand(time(NULL));
	quotenum = rand() % NUM_QUOTES;
	fprintf( stderr, "\n\t%s\n\n", quotes[quotenum]);
}

FILE* safe_fopen( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];

	file = fopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);

	}
	return file;
}

gzFile safe_fopen_gz( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	gzFile file;
	char err[500];

	file = gzopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);		
	}
	return file;
}

htsFile* safe_hts_open( char* path, char* mode)
{
	htsFile* bam_file;
	char err[500];

	bam_file = hts_open( path, mode);
	if( !bam_file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}

	return bam_file;
}


int is_proper( int flag)
{
	if ( (flag & BAM_FPAIRED) != 0 && (flag & BAM_FSECONDARY) == 0 && (flag & BAM_FSUPPLEMENTARY) == 0 && (flag & BAM_FDUP) == 0 && (flag & BAM_FQCFAIL) == 0)
		return 1;

	return 0;
}


int is_concordant( bam1_core_t bam_alignment_core, int min, int max)
{
	int flag = bam_alignment_core.flag;

	if( ( flag & BAM_FPAIRED) == 0)
	{
		/* Read is single-end. Skip this by calling it concordant */
		return RPCONC;
	}
	/*
	if( ( flag & BAM_FPROPER_PAIR) == 0)
	{
		/* Not proper pair
		return RPUNMAPPED;
	}*/

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
	{
		/* -- orientation = inversion */
		return RPINV;
	}

	if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
	{
		/* ++ orientation = inversion */
		return RPINV;
	}

	if( bam_alignment_core.tid != bam_alignment_core.mtid)
	{
		/* On different chromosomes */
		return RPINTERCHR;
	}

	if( bam_alignment_core.pos <= bam_alignment_core.mpos) // c.a.
	{
		/* Read is placed BEFORE its mate */
		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			/* -+ orientation = tandem duplication */
			return RPTDUP; //was 0 before
		}
	}
	else
	{
		/* Read is placed AFTER its mate */
		if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			/* +- orientation = tandem duplication */
			return RPTDUP; //was 0 before
		}
	}

	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs(bam_alignment_core.isize) < min) // c.a.
	{
		/* Deletion or Insertion */
		return RPINS;
	}
	else if(abs(bam_alignment_core.isize) > max)
	{
		return RPDEL;
	}

	/* All passed. Read is concordant */
	return RPCONC;
}





/* Decode 4-bit encoded bases to their corresponding characters */
char base_as_char( int base_as_int)
{
	if( base_as_int == 1)
	{
		return 'A';
	}
	else if( base_as_int == 2)
	{
		return 'C';
	}
	else if( base_as_int == 4)
	{
		return 'G';
	}
	else if( base_as_int == 8)
	{
		return 'T';
	}
	else if( base_as_int == 15)
	{
		return 'N';
	}
}

/* Return the complement of a base */
char complement_char( char base)
{
	switch( base)
	{
	case 'A':
		return 'T';
		break;
	case 'C':
		return 'G';
		break;
	case 'G':
		return 'C';
		break;
	case 'T':
		return 'A';
		break;
	default:
		return 'N';
		break;
	}
	return 'X';
}

/* Add 33 to the interger value of the qual characters to convert them to ASCII */
void qual_to_ascii( char* qual)
{
	int i;
	for( i = 0; i < strlen( qual); i++)
	{
		qual[i] = qual[i] + 33;
	}
}

/* Even safer than strncpy as it dynamically allocates space for the string if
 there hasn't been already */
void set_str( char** target, char* source)
{
	if( *target != NULL)
	{
		free( ( *target));
	}

	if (source != NULL)
	{
		( *target) = ( char*) getMem( sizeof( char) * ( strlen( source) + 1));
		strncpy( ( *target), source, ( strlen( source) + 1));
	}
	else
	{
		( *target) = NULL;
	}
}


/* Reverse a given string */
void reverse_string( char* str)
{
	int i;
	char swap;
	int len = strlen( str);

	for( i = 0; i < len / 2; i++)
	{
		swap = str[i];
		str[i] = str[len - i - 1];
		str[len - i - 1] = swap;
	}
}

int compare_size_int( const void* p, const void* q)
{
	int i = *( const int*) p;
	int j = *( const int*) q;

	if( i < j)
	{
		return -1;
	}
	else if( i == j)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void* getMem( size_t size)
{
	void* ret;

	ret = malloc( size);
	if( ret == NULL)
	{
		fprintf( stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory.\n", getMemUsage(), ( float) ( size / 1048576.0));
		exit( 0);
	}

	memUsage = memUsage + size;
	return ret;
}

void freeMem( void* ptr, size_t size)
{
	memUsage = memUsage - size;
	free( ptr);
}

double getMemUsage()
{
	return memUsage / 1048576.0;
}

int max( int x, int y)
{
	if( x < y)
		return y;
	else
		return x;
}

int min( int x, int y)
{
	if( x < y)
		return x;
	else
		return y;
}

int hammingDistance( char *str1, char *str2, int len)
{
	int dist = 0, i;
	for( i = 0; i < len; i++)
	{
		if( str1[i] != str2[i])
			dist++;
	}
	return dist;
}

int find_chr_index_bam( ref_genome* ref, char* chromosome_name, bam_hdr_t* bam_header)
{
	int i, len;
	char *tmp;

	for(i = 0; i < bam_header->n_targets; i++)
	{
		if( strcmp( chromosome_name, bam_header->target_name[i]) == 0)
			return i;

		/* Check if the chromosome name contains chr at the beginning */
		len = strlen( chromosome_name) + 4;
		tmp = ( char*) getMem( sizeof( char) * len);
		strcpy( tmp, "chr");
		strcat( tmp, chromosome_name);

		if( strcmp( tmp, bam_header->target_name[i]) == 0)
		{
			free( tmp);
			return i;
		}
		free( tmp);
	}
	return -1;
}

int vh_cmprReadNameStr (const void *a, const void *b)
{
	return strcmp (*(char **) a, *(char **) b);
}

unsigned long encode_ten_x_barcode(char* source){
	int i, len;
	unsigned long next_digit, result;

	if (source == NULL){
		return -1;
	}

	result = 0;
	next_digit = 0;
	len = strlen(source);
	for (i = 0; i < len; i++){
		switch(source[i]){
		case 'A':
			next_digit = 0;
			result = (result << 2)|next_digit;
			break;
		case 'C':
			next_digit = 1;
			result = (result << 2)|next_digit;
			break;
		case 'G' :
			next_digit = 2;
			result = (result << 2)|next_digit;
			break;
		case 'T':
			next_digit = 3;
			result = (result << 2)|next_digit;
			break;
		default:
			break;
		}
	}
	return result;
}
