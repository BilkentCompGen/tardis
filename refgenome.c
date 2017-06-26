#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <htslib/faidx.h>
#include <stdbool.h>
#include "refgenome.h"
#include "vh/vh_logger.h"

#define STR_SIZE 255

int load_refgen( ref_genome** ref, parameters* params)
{
  //  FILE* fai_file;
        int i;
  //	int ln_count = 0, i, c, length;
  //	char filename[STR_SIZE], first_arg[STR_SIZE], sec_arg[STR_SIZE];	

	(*ref) = ( ref_genome*) getMem( sizeof( ref_genome));
	(*ref)->ref_name = params->ref_genome;
	(*ref)->chrom_count = params->this_sonic->number_of_chromosomes;
	(*ref)->chrom_names = params->this_sonic->chromosome_names;
	(*ref)->chrom_lengths = params->this_sonic->chromosome_lengths;
	( *ref)->in_bam = ( bool*) getMem( sizeof( bool) * ( *ref)->chrom_count);

	for( i = 0; i < ( *ref)->chrom_count; i++)
	{
		( *ref)->in_bam[i] = false;
	}
	
	return RETURN_SUCCESS;
}

void init_rd( bam_info* in_bam, parameters* param, ref_genome* ref)
{
	int i, chr, lib_index;


	in_bam->read_depth = ( short**) getMem( ref->chrom_count * sizeof( short*));
	in_bam->read_count = 0;

	for( chr = 0; chr < ref->chrom_count; chr++)
	{
		in_bam->read_depth[chr] = ( short*) getMem( sizeof( short) * (ref->chrom_lengths[chr]));

		memset (in_bam->read_depth[chr], 0, (ref->chrom_lengths[chr] * sizeof(short)));
	}
}


void init_rd_per_chr( bam_info* in_bam, parameters* param, ref_genome* ref, int chr_index)
{
	int i, lib_index;

	in_bam->read_count = 0;
	in_bam->read_depth_per_chr = ( short*) getMem( sizeof( short) * ( ref->chrom_lengths[chr_index]));

	memset (in_bam->read_depth_per_chr, 0, (ref->chrom_lengths[chr_index] * sizeof(short)));
}


void calc_mu_per_chr( ref_genome* ref, bam_info* in_bam, int chr_index)
{
	int i;
	long rd_cnt = 0, window_total = 0;
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
	long rd_cnt = 0, window_total = 0;
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

	fprintf( logFile, "Read Count:%li  mean=%f\n",rd_cnt, in_bam->mean);

	fprintf( logFile,"Coverage is %.0lfx\n", round(in_bam->mean));
}


void calc_mean( ref_genome* ref, parameters *params, bam_info* in_bam)
{
	int lib_index, chr, i, gc_val, window_per_gc[101];
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
		        gc_val = (int) round ( sonic_get_gc_content(params->this_sonic, ref->chrom_names[chr], i, i+WINDOWSLIDE));			
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

void calc_mean_per_chr( ref_genome* ref, parameters *params, bam_info* in_bam, int chr_index)
{
	int lib_index, i, gc_val, window_per_gc[101];
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
		gc_val = (int) round ( sonic_get_gc_content(params->this_sonic, ref->chrom_names[chr_index], i, i+WINDOWSLIDE));
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
	}
}

int run_rd( bam_info** in_bam, parameters* params, ref_genome* ref)
{
	/* Variables */
	int return_value, lib_index, chr_index, bam_index, posLoc, i;
	int chrid_left;
	char* library_name = NULL, *chr_name;
	int *bamToRefIndex;

	htsFile* bam_file;
	bam_hdr_t* bam_header;
	bam1_t*	bam_alignment;
	bam1_core_t bam_alignment_core;
	hts_itr_t *iter;
	hts_idx_t *idx;

	fprintf( stderr,"10x flag: %d\n", ten_x_flag);

	for ( bam_index = 0; bam_index < params->num_bams; bam_index++)
	{
		/* Initialize the read depth(rd) array of each library for each bam file */
		init_rd( in_bam[bam_index], params, ref);

		/* HTS implementation */
		bam_file = safe_hts_open( params->bam_file_list[bam_index], "r");

		/* Read in BAM header information */
		bam_header = bam_hdr_read( ( bam_file->fp).bgzf);
		bam_alignment = bam_init1();

		/* Load the bam index file */
		idx = sam_index_load( bam_file, params->bam_file_list[bam_index]);
		if( idx == NULL)
			fprintf( stderr, "Error: IDX NULL");

		/* The array is used to map the chromosome indices in bam file to the ones in reference genome in case there is a difference */
		bamToRefIndex = ( int *) getMem( bam_header->n_targets * sizeof( int));
		for( i = 0; i < bam_header->n_targets; i++){
		        bamToRefIndex[i] = sonic_refind_chromosome_index( params->this_sonic, bam_header->target_name[i]);
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

		/* Mean value (mu) calculation */
		calc_mean( ref, params, in_bam[bam_index]);

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
