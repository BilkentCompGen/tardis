#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <htslib/faidx.h>
#include <stdbool.h>
#include "refgenome.h"
#include "vh/vh_logger.h"

#define STR_SIZE 255


void calc_gc( ref_genome** ref)
{
	int i, j, k = 0, loclength, gc_count, length;
	long min, max;
	float avg;
	char *ref_seq;
	faidx_t* ref_fai;

	/* Initialize the gc profile */
	( *ref)->gc = ( float**) getMem( sizeof( float*) * ( *ref)->chrom_count);

	for( i = 0; i < ( *ref)->chrom_count; i++)
	{
		length = ( (*ref)->chrom_lengths[i]);
		( *ref)->gc[i] = ( float*) getMem( sizeof( float) * length / ( WINDOWSLIDE - 1));
	}

	ref_fai = fai_load( (*ref)->ref_name);
	for( i = 0; i < (*ref)->chrom_count; i++)
	{
		min = 0;
		max = WINDOWSIZE;

		/* Retrieve the sequence of each chromosome to ref_seq array in each iteration */
		ref_seq = faidx_fetch_seq( ref_fai, (*ref)->chrom_names[i], min, (*ref)->chrom_lengths[i]-1, &loclength);
		k = 0;
		do
		{
			j = min;
			gc_count = 0;

			/* Count the number of G and C characters in the sequence */
			do
			{
				if( ref_seq[j] == 'G' || ref_seq[j] == 'C')
					gc_count++;
				j++;
			}while( j < max);

			/* Find the average G and C over total A, C, G, T and N(machine cannot read) count in each window */
			avg = ( float)gc_count / ( max - min);
			(*ref)->gc[i][k] = avg;
			(*ref)->window_count[i]++;
			k++;
			min += WINDOWSLIDE;
			if( max + WINDOWSLIDE > loclength)
				max += ( loclength - max);
			else
				max += WINDOWSLIDE;
		}while( min < loclength);

		free( ref_seq);

		//(*ref)->window_count[i]--;
	}
	fai_destroy( ref_fai);
}

void calc_gc_per_chr( ref_genome** ref, int chr_index)
{
	int j, k = 0, loclength, gc_count, length;
	long min, max;
	float avg;
	char *ref_seq;
	faidx_t* ref_fai;


	length = ( (*ref)->chrom_lengths[chr_index]);
	( *ref)->gc_per_chr = ( float*) getMem( sizeof( float) * length / ( WINDOWSLIDE - 1));

	ref_fai = fai_load( (*ref)->ref_name);
	min = 0;
	max = WINDOWSIZE;

	/* Retrieve the sequence of each chromosome to ref_seq array in each iteration */
	ref_seq = faidx_fetch_seq( ref_fai, (*ref)->chrom_names[chr_index], min, (*ref)->chrom_lengths[chr_index]-1, &loclength);
	k = 0;
	do
	{
		j = min;
		gc_count = 0;

		/* Count the number of G and C characters in the sequence */
		do
		{
			if( ref_seq[j] == 'G' || ref_seq[j] == 'C')
				gc_count++;
			j++;
		}while( j < max);

		/* Find the average G and C over total A, C, G, T and N(machine cannot read) count in each window */
		avg = ( float)gc_count / ( max - min);
		(*ref)->gc_per_chr[k] = avg;
		(*ref)->window_count[chr_index]++;
		k++;
		min += WINDOWSLIDE;
		if( max + WINDOWSLIDE > loclength)
			max += ( loclength - max);
		else
			max += WINDOWSLIDE;
	}while( min < loclength);

	free( ref_seq);

	fai_destroy( ref_fai);
}


int load_refgen( ref_genome** ref, parameters* params)
{
	FILE* fai_file;
	int ln_count = 0, i, c, length;
	long total_length;
	char filename[STR_SIZE], first_arg[STR_SIZE], sec_arg[STR_SIZE];
	int retval;

	retval = RETURN_SUCCESS;

	(*ref) = ( ref_genome*) getMem( sizeof( ref_genome));
	(*ref)->ref_name = params->ref_genome;

	sprintf( filename, "%s.fai", (*ref)->ref_name);
	fai_file = safe_fopen( filename, "r");

	/* Count the number of chromosomes by counting the non-empty lines in the .fai file */
	do{
		c = fgetc(fai_file);
		if( c == '\n') ln_count++;

	}while( c != EOF);
	(*ref)->chrom_count = ln_count;

	/* Reset the file pointer to the start of the file */
	rewind( fai_file);

	(*ref)->chrom_lengths = ( int*) getMem( (*ref)->chrom_count * sizeof(int));
	(*ref)->chrom_names = ( char**) getMem( (*ref)->chrom_count * sizeof(char*));

	for( i = 0; i < (*ref)->chrom_count; i++)
	{
		retval = fscanf( fai_file, "%s%s%*s%*s%*d", first_arg, sec_arg);
		if ( retval <= 0)
			;
		(*ref)->chrom_names[i] = NULL;
		set_str( &((*ref)->chrom_names[i]), first_arg);
		(*ref)->chrom_lengths[i] = atoi( sec_arg);
	}
	fclose( fai_file);

	/* Calculate total length of the genome */
	total_length = 0;
	for( i = 0; i < (*ref)->chrom_count; i++)
	{
		total_length += ( *ref)->chrom_lengths[i];
	}
	( *ref)->gen_length = total_length;
	( *ref)->in_bam = ( bool*) getMem( sizeof( bool) * ( *ref)->chrom_count);

	/* Initialize window count array */
	( *ref)->window_count = ( int*) getMem( sizeof( int) * ( *ref)->chrom_count);

	for( i = 0; i < ( *ref)->chrom_count; i++)
	{
		( *ref)->window_count[i] = 0;
		( *ref)->in_bam[i] = false;
	}

	retval = RETURN_SUCCESS;
	/* TODO
	if ( params->load_sonic)
	{
		retval = load_sonic(params);
		} */

	return retval;
}

void init_rd( bam_info* in_bam, parameters* param, ref_genome* ref)
{
	int i, chr, lib_index;

	in_bam->read_depth = ( short**) getMem( RDCHRCOUNT * sizeof( int*));
	in_bam->read_count = 0;

	for( chr = 0; chr < RDCHRCOUNT; chr++)
	{
		in_bam->read_depth[chr] = ( short*) getMem( sizeof( short) * (ref->chrom_lengths[chr]));

		for( i = 0; i < ref->chrom_lengths[chr]; i++)
			in_bam->read_depth[chr][i] = 0;
	}
}


void init_rd_per_chr( bam_info* in_bam, parameters* param, ref_genome* ref, int chr_index)
{
	int i, lib_index;

	in_bam->read_count = 0;

	in_bam->read_depth_per_chr = ( short*) getMem( sizeof( short) * (ref->chrom_lengths[chr_index]));

	for( i = 0; i < ref->chrom_lengths[chr_index]; i++)
		in_bam->read_depth_per_chr[i] = 0;
}


void calc_mu_per_chr( ref_genome* ref, bam_info* in_bam, int chr_index)
{
	int i;
	long total_length = 0, rd_cnt = 0, window_total = 0;
	double cov;

	for( i = 0; i < ref->chrom_lengths[chr_index]; i++)
	{
		rd_cnt += (long) in_bam->read_depth_per_chr[i];
		window_total++;
	}
	/* Calculate mu values */
	in_bam->mean = ( double)(rd_cnt) / ( double)window_total;

	fprintf( logFile, "Read Count:%li  Window count:%li mean=%f\n", rd_cnt, window_total, in_bam->mean);

	//cov = ( double) rd_cnt * (in_bam->libraries[lib_index]->read_length) / ( double) ( window_total * WINDOWSIZE );
	//fprintf( logFile,"Coverage is %.0lfx\n", round(cov));

}


void calc_mu( ref_genome* ref, bam_info* in_bam)
{
	int chr,i;
	long total_length = 0, rd_cnt = 0, window_total = 0;
	double cov;

	/* Count the reads for each chromosome */
	for( chr = 0; chr < ref->chrom_count; chr++)
	{
		if( ref->in_bam[chr] == false)
			continue;

		for( i = 0; i < ref->chrom_lengths[chr]; i++)
		{
			rd_cnt += (long)in_bam->read_depth[chr][i];
			window_total++;
		}
	}

	/* Calculate mu values */
	in_bam->mean = ( double)(rd_cnt) / ( double)window_total;

	//fprintf( stderr, "RDX= %li RDY=%li length_x=%li length_y=%li\n", rd_cnt_x, rd_cnt_y, length_x, length_y);

	fprintf( logFile, "Read Count:%li  Total Read Length:%li mean=%f\n",rd_cnt, total_length, in_bam->mean);

	//cov = ( double) ( rd_cnt_x + rd_cnt_y + rd_cnt) * (in_bam->libraries[lib_index]->read_length) / ( double) ( total_length + length_x + length_y);
	fprintf( logFile,"Coverage is %.0lfx\n", round(in_bam->mean));
}


void calc_mean( ref_genome* ref, bam_info* in_bam)
{
	int lib_index, chr, i, gc_val, gc_window, window_per_gc[101];
	long rd_per_gc[101];

	calc_mu( ref, in_bam);

	/* Calculate mu_GC values */
	for( i = 0; i < 101; i++)
	{
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
	}
	for( chr = 0; chr < ref->chrom_count; chr++)
	{
		if( ref->in_bam[chr] == false)
			continue;

		for( i = 0; i < ref->chrom_lengths[chr]; i++)
		{
			gc_window = ( int) floorf( ( float) i / WINDOWSLIDE);
			gc_val = ( int) round( ( float)ref->gc[chr][gc_window] * 100);

			rd_per_gc[gc_val] += ( long) in_bam->read_depth[chr][i];
			window_per_gc[gc_val]++;
		}
	}
	in_bam->mean_rd_per_gc[0] = 0.0;
	for( i = 1; i < 101; i++)
	{
		in_bam->mean_rd_per_gc[i] = ( float)rd_per_gc[i] / ( window_per_gc[i]);
		if( isnanf( in_bam->mean_rd_per_gc[i]) || isinff( ( in_bam->mean_rd_per_gc[i])) == -1
				|| isinff( ( in_bam->mean_rd_per_gc[i])) == 1 )
			in_bam->mean_rd_per_gc[i] = 0;
	}
}

void calc_mean_per_chr( ref_genome* ref, bam_info* in_bam, int chr_index)
{
	int lib_index, i, gc_val, window_per_gc[101], gc_window;
	long rd_per_gc[101];

	calc_mu_per_chr( ref, in_bam, chr_index);

	/* Calculate mu_GC values */
	for( i = 0; i < 101; i++)
	{
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
	}

	for( i = 0; i < ref->chrom_lengths[chr_index]; i++)
	{
		gc_window = ( int) floorf( ( float) i / WINDOWSLIDE);
		gc_val = ( int) round( ( float)ref->gc_per_chr[gc_window] * 100);

		//fprintf(stderr,"i=%d gc_wind=%d gc_val=%d depht=%d\n",i,gc_window,gc_val, in_bam->read_depth[0][i]);
		rd_per_gc[gc_val] += ( long) in_bam->read_depth_per_chr[i];
		window_per_gc[gc_val]++;
	}

	in_bam->mean_rd_per_gc[0] = 0.0;
	for( i = 1; i < 101; i++)
	{
		in_bam->mean_rd_per_gc[i] = ( float)rd_per_gc[i] / ( window_per_gc[i]);
		if( isnanf( in_bam->mean_rd_per_gc[i]) || isinff( ( in_bam->mean_rd_per_gc[i])) == -1
				|| isinff( ( in_bam->mean_rd_per_gc[i])) == 1 )
			in_bam->mean_rd_per_gc[i] = 0;
		//fprintf(stderr, "\ni = %d rd_per_gc=%li - mean=%f",i, rd_per_gc[i], in_bam->mean_rd_per_gc[i] );
	}
}

int run_rd( bam_info** in_bam, parameters* param, ref_genome* ref)
{
	/* Variables */
	int return_value, lib_index, chr_index, bam_index, posLoc, i;
	int chrid_left;
	char* library_name = NULL, *chr_name;
	int *bamToRefIndex;
	fprintf( stderr,"10x flag: %d\n", ten_x_flag);

	htsFile* bam_file;
	bam_hdr_t* bam_header;
	bam1_t*	bam_alignment;
	bam1_core_t bam_alignment_core;
	hts_itr_t *iter;
	hts_idx_t *idx;

	for ( bam_index = 0; bam_index < param->num_bams; bam_index++)
	{
		/* Initialize the read depth(rd) array of each library for each bam file */
		init_rd( in_bam[bam_index], param, ref);

		/* HTS implementation */
		bam_file = safe_hts_open( param->bam_file_list[bam_index], "r");

		/* Read in BAM header information */
		bam_header = bam_hdr_read( ( bam_file->fp).bgzf);
		bam_alignment = bam_init1();

		/* Load the bam index file */
		idx = sam_index_load( bam_file, param->bam_file_list[bam_index]);
		if( idx == NULL)
			fprintf( stderr, "Error: IDX NULL");

		/* The array is used to map the chromosome indices in bam file to the ones in reference genome in case there is a difference */
		bamToRefIndex = ( int *) getMem( bam_header->n_targets * sizeof( int));
		for( i = 0; i < bam_header->n_targets; i++){
			bamToRefIndex[i] = find_chr_index_ref( ref, bam_header->target_name[i]);
			// if(bamToRefIndex[i] < 0){
			//     fprintf( stderr, "ReadDepth: Couldn't find chrm %s index in the reference :(\n", bam_header->target_name[i]);
			// }
		}
		for ( chr_index = 0; chr_index < bam_header->n_targets; chr_index++)
		{
			chr_name = bam_header->target_name[chr_index];

			/* If a chromosome in bam is missing in the reference, skip it */
			if( bamToRefIndex[chr_index] == -1)
				continue;

			fprintf( stderr, "\r                                                        ");
			fflush( stderr);
			fprintf( stderr, "\rReading the bam file of %s - chromosome %s", in_bam[bam_index]->sample_name, chr_name);
			fflush( stderr);

			iter = bam_itr_queryi( idx, chr_index, 0, bam_header->target_len[chr_index]);
			if( iter == NULL)
				fprintf( stderr, "Error: Iter NULL");

			ref->in_bam[bamToRefIndex[chr_index]] = true;

			while( bam_itr_next( bam_file, iter, bam_alignment) > 0)
			{
				bam_alignment_core = bam_alignment->core;
				chrid_left = bamToRefIndex[bam_alignment_core.tid];

				/* Get library index */
				set_str( &library_name, bam_aux_get( bam_alignment, "RG"));
				lib_index = find_library_index( in_bam[bam_index], library_name + 1);

				/* Increase the read count and depth */
				in_bam[bam_index]->read_count++;
				in_bam[bam_index]->read_depth[chrid_left][bam_alignment_core.pos]++;
			}
		}

		fflush( stderr);
		fprintf( stderr, "\n");
		vh_logInfo( "Finished reading the bam file");

		vh_logInfo ( "Calculating GC profile");
		calc_gc( &ref);

		/* Mean value (mu) calculation */
		calc_mean( ref, in_bam[bam_index]);

		/* Close the BAM file */
		return_value = hts_close( bam_file);
		if( return_value != 0)
		{
			fprintf( stderr, "Error closing BAM file\n");
			exit( 1);
		}
		free(bamToRefIndex);

		/* Free the bam related files */
		bam_hdr_destroy( bam_header);
	}
	return RETURN_SUCCESS;
}
