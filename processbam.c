#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <inttypes.h>

/* tardis headers */

#include "common.h"
#include "processbam.h"
#include "processfq.h"

bool isRGAvailable = true;

void load_bam( parameters *params, bam_info* in_bam, char* path, int alternative, int bam_cnt, char* ref_genome_path)
{
	/* Variables */
	htsFile* bam_file;
	bam_hdr_t* bam_header;
	bam1_core_t bam_alignment_core;
	bam1_t*	bam_alignment;
	int** fragment_size;
	int** second_pass_fragments;
	int* fragments_sampled;
	int* second_test_pass;
	int* fragment_size_total;
	float* variance;
	char* library_name = NULL;
	int i, j, library_index, diff, return_value;
	char fai_file[MAX_SEQ];

	fprintf( stderr, "Processing BAM file %s.\n", path);


	/* Open the BAM file for reading. htslib automatically detects the format
		of the file, so appending "b" after "r" in mode is redundant. */
	bam_file = safe_hts_open( path, "r");

	/* Read in BAM header information */
	//bam_header = bam_hdr_read( ( bam_file->fp).bgzf);
	bam_header = sam_hdr_read(bam_file);

	/* Extract the Sample Name from the header text */
	get_sample_name( in_bam, bam_header->text);

	/* Extract the number of libraries within the BAM file */
	get_library_count( in_bam, bam_header->text);

	if ( in_bam->num_libraries != 0)
	{
		fprintf( stderr, "Total %d %s found in %s. Sample name is %s.\n", in_bam->num_libraries, in_bam->num_libraries > 1 ? "libraries":"library", path, in_bam->sample_name);
		fprintf( logFile, "\nSample name = %s, total library count = %d\n", in_bam->sample_name, in_bam->num_libraries);
	}
	else
	{
		fprintf( stderr, "[BAM FILE ERROR] No libraries (RG) found in the BAM file. Setting default library as Lib1.\n");
		in_bam->num_libraries = 1;
		isRGAvailable = false;
		sprintf(in_bam->sample_name, "Sample%d", bam_cnt);
		/*
		fprintf( stderr, "[BAM FILE ERROR] No libraries found in the BAM file. Do you have read groups set in the file? Use picard-tools (AddOrRemoveReadGroups) to fix your BAM file.\nExiting.\n");
		exit (EXIT_READGROUP);
		 */
	}

	/* Initialize the structures for library properties */
	in_bam->libraries = ( struct library_properties**) getMem( (in_bam->num_libraries + 1) * sizeof( struct library_properties*));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		( in_bam->libraries)[i] = ( struct library_properties*) getMem( sizeof( struct library_properties));
		( in_bam->libraries)[i]->libname = NULL;
		( in_bam->libraries)[i]->fastq1 = NULL;
		( in_bam->libraries)[i]->fastq2 = NULL;
		( in_bam->libraries)[i]->divet = NULL;
		( in_bam->libraries)[i]->mappings_discordant = ( discordantMapping **) getMem( (NHASH + 1) * sizeof( discordantMapping *));
		for( j = 0; j < NHASH; j++)
			( in_bam->libraries)[i]->mappings_discordant[j] = NULL;

		if( alternative != 0)
		{
			( in_bam->libraries)[i]->mappings_alternative = ( alternativeMapping **) getMem( (NHASH + 1) * sizeof( alternativeMapping *));
			for( j = 0; j < NHASH; j++)
				( in_bam->libraries)[i]->mappings_alternative[j] = NULL;
		}

		( in_bam->libraries)[i]->listMEI_Mapping = NULL;
		( in_bam->libraries)[i]->listNUMT_Mapping = NULL;
		( in_bam->libraries)[i]->listSoftClip = NULL;
	}

	/* Extract the ids/names for the libraries. A single Sample with multiple
	 possible libraries are assumed for each BAM file */
	if( isRGAvailable == true)
		get_library_names( in_bam, bam_header->text);
	else
		in_bam->libraries[0]->libname = "Lib1";

	/* For SAMPLEFRAG number of alignments, store the template length field.
	 Performed for each different library */
	fragment_size = ( int**) getMem( in_bam->num_libraries * sizeof( int*));
	
	for( i = 0; i < in_bam->num_libraries; i++)
	{
	        fragment_size[i] = NULL;
		fragment_size[i] = ( int*) getMem( SAMPLEFRAG * sizeof( int));
	}

	sprintf( fai_file,"%s.fai", ref_genome_path);
	hts_set_fai_filename( bam_file, fai_file);

	/* Initial read */
	bam_alignment = bam_init1();
	//return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
	return_value = sam_read1( bam_file, bam_header, bam_alignment);

	/* The remaining reads */
	fragments_sampled = ( int*) getMem( in_bam->num_libraries * sizeof( int));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		fragments_sampled[i] = 0;
	}

	fprintf( stderr, "Sampling reads from libraries to infer fragment sizes.\n");



	while( return_value != -1 && !sufficient_fragments_sampled( fragments_sampled, in_bam->num_libraries))
	{
		bam_alignment_core = bam_alignment->core;

		if( bam_alignment_core.isize > 0 && !bam_is_rev( bam_alignment) && bam_is_mrev( bam_alignment))
		{
			if( isRGAvailable == true)
			{
				set_str( &library_name, bam_aux_get( bam_alignment, "RG"));

				/* get rid of the leading 'Z' in front of the library_name */
				library_index = find_library_index( in_bam, library_name + 1);
			}
			else
				library_index = 0;

			/* Sample SAMPLEFRAG number of alignments for each library at most */
			if( library_index != -1 && fragments_sampled[library_index] < SAMPLEFRAG)
			{
				fragment_size[library_index][fragments_sampled[library_index]] = bam_alignment_core.isize;
				fragments_sampled[library_index] = fragments_sampled[library_index] + 1;
			}
		}

		/* Read next alignment */
		//return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
		return_value = sam_read1( bam_file, bam_header, bam_alignment);
	}

	fprintf( stderr, "Sampling finished. Now calculating library statistics.\n");

	/* Now we have SAMPLEFRAG number of fragment sizes which are positive and pass the flag conditions */
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		/* Sort the fragment sizes */
		qsort( fragment_size[i], SAMPLEFRAG, sizeof( int), compare_size_int);

		/* Get the medians */
		( in_bam->libraries)[i]->frag_med = fragment_size[i][( SAMPLEFRAG / 2) - 1];
	}

	/* Find the fragment sizes which pass the second test, and will contribute to the avg and std */
	second_pass_fragments = ( int**) getMem( in_bam->num_libraries * sizeof( int*));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		second_pass_fragments[i] = ( int*) getMem( SAMPLEFRAG * sizeof( int));
	}

	second_test_pass = ( int*) getMem( in_bam->num_libraries * sizeof( int));
	fragment_size_total = ( int*) getMem( in_bam->num_libraries * sizeof( int));

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		second_test_pass[i] = 0;
		fragment_size_total[i] = 0;
	}


	for( i = 0; i < in_bam->num_libraries; i++)
	{
		for( j = 0; j < SAMPLEFRAG; j++)
		{
			if( fragment_size[i][j] <= 2 * ( in_bam->libraries)[i]->frag_med)
			{
				fragment_size_total[i] = fragment_size_total[i] + fragment_size[i][j];
				second_pass_fragments[i][j] = fragment_size[i][j];
				second_test_pass[i] = second_test_pass[i] + 1;
			}
		}

	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		/* Compute the averages */
		( in_bam->libraries)[i]->frag_avg = ( float) fragment_size_total[i] / ( float) second_test_pass[i];

		/* Set the read length for the current library */
		( in_bam->libraries)[i]->read_length = bam_alignment_core.l_qseq;
	}

	/* Compute the variance and std */
	variance = ( float*) getMem( in_bam->num_libraries * sizeof( float));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		variance[i] = 0;
	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		for( j = 0; j < second_test_pass[i]; j++)
		{
			diff = ( second_pass_fragments[i][j] - ( in_bam->libraries)[i]->frag_avg);
			variance[i] = variance[i] +  diff * diff;
		}
	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		variance[i] = ( float) variance[i] / ( float) second_test_pass[i];
		( in_bam->libraries)[i]->frag_std = sqrt( variance[i]);
		fprintf( stderr, "\nLibrary %s\n\tMean: %f\n\tStdev: %f\n\n", ( in_bam->libraries)[i]->libname, ( in_bam->libraries)[i]->frag_avg, ( in_bam->libraries)[i]->frag_std);
		fprintf( logFile, "\nLibrary %s\n\tMean: %f\n\tStdev: %f\n\tRead Length: %d\n\n", ( in_bam->libraries)[i]->libname, ( in_bam->libraries)[i]->frag_avg, ( in_bam->libraries)[i]->frag_std, ( in_bam->libraries)[i]->read_length);
		set_library_min_max( ( in_bam->libraries)[i]);
	}

	/* Close the BAM file */
	return_value = hts_close( bam_file);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing BAM file\n");
		exit( EXIT_BAM_CLOSE);
	}

	/* Free Memory */
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		free( fragment_size[i]);
		free( second_pass_fragments[i]);
	}
	free( fragment_size);
	free( second_pass_fragments);
	free( fragments_sampled);
	free( second_test_pass);
	free( fragment_size_total);
	free( variance);
	free( library_name);
}


void get_sample_name(bam_info* in_bam, char* header_text)
{
	/* Delimit the BAM header text with tabs and newlines */

	char *tmp_header = NULL;
	set_str( &( tmp_header), header_text);
	char* p = strtok( tmp_header, "\t\n");
	char sample_name_buffer[1024];

	while( p != NULL)
	{
		/* If the current token has "SM" as the first two characters,
			we have found our Sample Name */
		if( p[0] == 'S' && p[1] == 'M')
		{
			/* Get the Sample Name */
			strncpy( sample_name_buffer, p + 3, strlen( p) - 3);

			/* Add the NULL terminator */
			sample_name_buffer[strlen( p) - 3] = '\0';

			/* Exit loop */
			break;
		}
		p = strtok( NULL, "\t\n");
	}

	set_str( &( in_bam->sample_name), sample_name_buffer);
	free( tmp_header);
}

void get_library_count( bam_info* in_bam, char* header_text)
{
	int number_of_libraries = 0;

	/* Delimit the BAM header text with newlines */
	char *tmp_header = NULL;
	set_str( &( tmp_header), header_text);
	char* p = strtok( tmp_header, "\n");

	while( p != NULL)
	{
		/* If the current token (which is a line) has "RG" as the second and third characters,
		 we have found a new library */
		if( p[1] == 'R' && p[2] == 'G')
		{
			number_of_libraries = number_of_libraries + 1;
		}
		p = strtok( NULL, "\n");
	}

	in_bam->num_libraries = number_of_libraries;
	free( tmp_header);
}

void get_library_names( bam_info* in_bam, char* header_text)
{
	char line_buffer[1024];
	char library_name_buffer[1024];
	int i;
	char *saveptr_p, *saveptr_q;

	/* Delimit the BAM header text with newlines */
	i = 0;

	char *tmp_header = NULL;
	set_str( &( tmp_header), header_text);

	char* p = strtok_r( tmp_header, "\n", &saveptr_p);

	while( p != NULL)
	{
		/* If the current token (which is a line) has "RG" as the second and third characters,
			we copy the current line to a new char* and tokenize it to get the id/name */
		if( p[1] == 'R' && p[2] == 'G')
		{
			/* Copy the current line to the line buffer */
			strncpy( line_buffer, p, strlen( p));
			line_buffer[strlen( p)] = '\0';

			/* Tokenize the line buffer by tabs */
			char* q = strtok_r( line_buffer, "\t", &saveptr_q);
			while( q != NULL)
			{
				if( q[0] == 'I' && q[1] == 'D')
				{
					/* Get the Library Name */
					strncpy( library_name_buffer, q + 3, strlen( q) - 3);

					/* Add the NULL terminator */
					library_name_buffer[strlen( q) - 3] = '\0';

					/* Exit loop */
					break;
				}

				q = strtok_r( NULL, "\t", &saveptr_q);
			}
			set_str( &( ( in_bam->libraries)[i]->libname), library_name_buffer);
			i = i + 1;
		}
		p = strtok_r( NULL, "\n", &saveptr_p);
	}
	free( tmp_header);
}

int find_library_index( bam_info* in_bam, char* library_name)
{
	int library_index;
	int i;

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		if( strcmp( ( in_bam->libraries)[i]->libname, library_name) == 0)
		{
			return i;
		}
	}

	return -1;
}

int sufficient_fragments_sampled( int* fragments_sampled, int num_libraries)
{
	int i;
	for( i = 0; i < num_libraries; i++)
	{
		if( fragments_sampled[i] != SAMPLEFRAG)
		{
			return 0;
		}
	}

	return 1;
}


void create_fastq( bam_info* in_bam, char* bam_path, parameters* params)
{
	int i;
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		fprintf( stderr, "Creating FASTQ files for the library: %s.\n", ( in_bam->libraries)[i]->libname);
		create_fastq_library( ( in_bam->libraries)[i], in_bam->sample_name, bam_path, params);
	}
}

void set_library_min_max( struct library_properties* in_lib)
{
	in_lib->conc_min = in_lib->frag_avg - ( 4 * in_lib->frag_std);
	if ( in_lib->conc_min < 0)
	{
		in_lib->conc_min = 0;
	}
	in_lib->conc_max = in_lib->frag_avg + ( 4 * in_lib->frag_std);
}


int is_concordant_bamonly( int pos_left, int pos_right, uint16_t flag, int32_t isize, int min, int max)
{
	if( ( flag & BAM_FPAIRED) == 0)
	{
		/* Read is single-end. Skip this by calling it concordant */
		return RPCONC;
	}

	/*if( ( flag & BAM_FPROPER_PAIR) == 0)
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
	/*
	if( chrID_left != chrID_right)
	{
		/* On different chromosomes
		return RPINTERCHR;
	}
	 */
	if( pos_left <= pos_right) // c.a.
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
	if( abs( isize) < min) // c.a.
	{
		/* Deletion or Insertion */
		return RPINS;
	}
	else if( abs( isize) > max)
	{
		return RPDEL;
	}

	/* All passed. Read is concordant */
	return RPCONC;
}

