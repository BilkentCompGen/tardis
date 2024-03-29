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
int cluster_of_reads;

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
	( *params)->outdir = NULL;
	( *params)->bam_file_list = ( char**) getMem( sizeof( char*) * MAX_BAMS);
	( *params)->gaps = NULL;
	( *params)->mei = NULL;
	( *params)->force_read_length = 0;
	//( *params)->threads = 1;
	( *params)->num_bams = 0;
	( *params)->num_samples = 0;
	( *params)->sample_names = NULL;
	( *params)->size_samples_array = MAX_SAMPLES;
	//( *params)->skip_mrfast = 0;
	( *params)->ten_x = 0;
	( *params)->output_hs = 0;
	( *params)->first_chr = 0;
	( *params)->last_chr = -1;
	( *params)->make_sonic = 0;
	( *params)->load_sonic = 0;
	( *params)->sonic_info = NULL;
	( *params)->number_of_different_mei_types = 0;
	( *params)->quick = 1;
	( *params)->ref_seq = NULL;
	( *params)->hash_size = 0;
	( *params)->histogram_only = 0;

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
	//fprintf( logFile, "%-30s%d\n","10x Tag:", params->ten_x);
	fprintf( logFile, "%-30s%d\n","First chrom:", params->first_chr);
	fprintf( logFile, "%-30s%d\n","Last chrom:", params->last_chr);
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

/* check if ACGT */
int is_dna_letter( char base){
  if (base == 'A')
    return 1;
  if (base == 'C')
    return 1;
  if (base == 'G')
    return 1;
  if (base == 'T')
    return 1;
  return 0;
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

void* reallocMem( void* ptr, size_t old_size, size_t new_size)
{
        ptr = realloc( ptr, new_size);
	memUsage = memUsage - old_size + new_size;
	return ptr;
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

int find_chr_index_bam( char* chromosome_name, bam_hdr_t* bam_header)
{
	int i, len;
	char *tmp;

	for(i = 0; i < bam_header->n_targets; i++)
	{
		if( strcmp( chromosome_name, bam_header->target_name[i]) == 0)
			return i;
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


int32_t calculateInsertSize( int32_t pos_left, int32_t pos_right, uint16_t flag, int read_length)
{
	int32_t isize = 0;

	/*
	if( chrID_left != chrID_right)
		isize = 0;
		*/

	/* Calculate the insert size */
	if( ( flag & BAM_FREVERSE) == 0) /* if left pair is forward */
	{
		/*If right pair is reverse */
		if( ( flag & BAM_FMREVERSE) != 0)
			isize = ( pos_right + read_length) - ( pos_left);
		else
			isize = pos_right - pos_left;
	}
	else if( ( flag & BAM_FREVERSE) != 0) /* if left pair is reverse */
	{
		/*If right pair is reverse */
		if( ( flag & BAM_FMREVERSE) != 0)
			isize = ( pos_right + read_length) - ( pos_left + read_length);
		else
			isize = ( pos_right) - ( pos_left + read_length);
	}
	return isize;
}

void get_working_directory(parameters *params){
  char *directory;
  char *prefix;
  int i;

  
  directory = strrchr(params->outprefix, '/');
  prefix = NULL;
  if (directory != NULL)
    set_str( &prefix, directory + 1);
  if ( debug_mode && prefix != NULL)
    fprintf (stderr, "prefix: %s\n", prefix);
  
  if (directory == NULL){
    //set_str(&(params->outdir), "./");
    set_str(&(params->outdir), "");
    return;
  }
    
  i = 0;
  while (params->outprefix+i != directory)
    i++;

  params->outdir = (char *) getMem((i+2) * sizeof(char));
  memcpy(params->outdir, params->outprefix, i*sizeof(char));
  params->outdir[i]='/';
  params->outdir[i+1] = 0;

  if ( debug_mode && prefix != NULL)
    fprintf (stderr, "prefix2: %s\n", prefix);

  if (prefix!=NULL){
    free( params->outprefix);
    params->outprefix = NULL;
    set_str( &(params->outprefix), prefix);
    free (prefix);
  }
}

char *get_file_name(char *path)
{
  char *ret;
  ret = strrchr (path, '/');
  if (ret != NULL)
    return ret + 1;
  return path;
}
    

void clean_up_temp_files(parameters *params){
  char *divetfile_path;

  divetfile_path = (char *) getMem(sizeof(char) * (3+strlen("divv.vh")+strlen(params->outprefix)+strlen(params->outdir)));
  sprintf(divetfile_path, "%s%s-%s", params->outdir, params->outprefix, "divv.vh");
  remove(divetfile_path);
  free(divetfile_path);
 
  char *outputscore_path;
  outputscore_path = (char *) getMem(sizeof(char) * (3+strlen("output.score")+strlen(params->outprefix)+strlen(params->outdir)));
  sprintf(outputscore_path, "%s%s-%s", params->outdir, params->outprefix, "output.score");
  remove(outputscore_path);
  free(outputscore_path);
  
  char *debugsr_path;
  debugsr_path = (char *) getMem(sizeof(char) * (3+strlen("debug.sr")+strlen(params->outprefix)+strlen(params->outdir)));
  sprintf(debugsr_path, "%s%s-%s", params->outdir, params->outprefix, "debug.sr");
  remove(debugsr_path);
  free(debugsr_path);
}


int compare_sonic_ref(parameters *params){

  faidx_t* ref_fai;
  int is_present = 0;
  int i;
  
  ref_fai = fai_load( params->ref_genome);

  for( i = params->first_chr; i <= params->last_chr; i++)
    {
      is_present = faidx_has_seq(ref_fai, params->this_sonic->chromosome_names[i]);
      if ( is_present)
	break;
    }
  
  fai_destroy( ref_fai);

  return is_present;
}
