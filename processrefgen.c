#include "processrefgen.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <htslib/faidx.h>
#include <stdbool.h>
#include "vh/vh_logger.h"
#include "processfq.h"

#define STR_SIZE 255

/*
int load_refgen( ref_genome** ref, parameters* params)
{
	int i;

	(*ref) = ( ref_genome*) getMem( sizeof( ref_genome));
	(*ref)->ref_name = params->ref_genome;
	(*ref)->chrom_count = params->this_sonic->number_of_chromosomes;
	(*ref)->chrom_names = params->this_sonic->chromosome_names;
	(*ref)->chrom_lengths = params->this_sonic->chromosome_lengths;
	( *ref)->in_bam = ( bool*) getMem( sizeof( bool) * ( *ref)->chrom_count);

	for( i = 0; i < ( *ref)->chrom_count; i++)
		( *ref)->in_bam[i] = false;

	return RETURN_SUCCESS;
}*/

void init_rd( bam_info* in_bam, parameters* param)
{
	int i, chr, lib_index;

	in_bam->read_depth = ( short**) getMem( param->this_sonic->number_of_chromosomes * sizeof( short*));
	in_bam->read_count = 0;

	for( chr = param->first_chr; chr <= param->last_chr; chr++)
	{
		in_bam->read_depth[chr] = ( short*) getMem( sizeof( short) * ( param->this_sonic->chromosome_lengths[chr]));
		memset (in_bam->read_depth[chr], 0, (param->this_sonic->chromosome_lengths[chr] * sizeof(short)));

		/* Initialize all the depths to zero */
		for( i = 0; i < param->this_sonic->chromosome_lengths[chr]; i++)
		{
			in_bam->read_depth[chr][i] = 0;
		}
	}
}


void init_rd_per_chr( bam_info* in_bam, parameters* param, int chr_index)
{
	int i, lib_index;

	in_bam->read_count = 0;
	in_bam->read_depth_per_chr = ( short*) getMem( sizeof( short) * ( param->this_sonic->chromosome_lengths[chr_index]));

	memset (in_bam->read_depth_per_chr, 0, (param->this_sonic->chromosome_lengths[chr_index] * sizeof(short)));
}


void calc_mu_per_chr( bam_info* in_bam, int chromosome_length)
{
	int i;
	long rd_cnt = 0, window_total = 0;
	double cov;

	for( i = 0; i < chromosome_length; i++)
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


void calc_mu( bam_info* in_bam, parameters *params)
{
	int chr,i;
	long rd_cnt = 0, window_total = 0;
	double cov;

	/* Count the reads for each chromosome */
	for( chr = params->first_chr; chr <= params->last_chr; chr++)
	{
		for( i = 0; i <= params->this_sonic->chromosome_lengths[chr]; i++)
		{
			rd_cnt += (long)in_bam->read_depth[chr][i];
			window_total++;
		}
	}

	/* Calculate mu values */
	in_bam->mean = ( double)(rd_cnt) / ( double)window_total;

	fprintf( logFile, "Total Read Count= %li\tMean= %f\tCoverage= %.0lfx\n",rd_cnt, in_bam->mean,
			round( in_bam->mean * in_bam->libraries[0]->read_length));
}


void calc_mean( parameters *params, bam_info* in_bam)
{
	int chr, i, gc_val, window_per_gc[101];
	long rd_per_gc[101];

	calc_mu( in_bam, params);

	/* Calculate mu_GC values */
	for( i = 0; i < 101; i++)
	{
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
	}
	for( chr = params->first_chr; chr <= params->last_chr; chr++)
	{
		for( i = 0; i < params->this_sonic->chromosome_lengths[chr]; i++)
		{
			gc_val = (int) round ( sonic_get_gc_content(params->this_sonic, params->this_sonic->chromosome_names[chr], i, ( i + WINDOWSLIDE)));
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

void calc_mean_per_chr( parameters *params, bam_info* in_bam, int chr_index)
{
	int lib_index, i, gc_val, window_per_gc[101];
	long rd_per_gc[101];

	calc_mu_per_chr( in_bam, params->this_sonic->chromosome_lengths[chr_index]);

	/* Calculate mu_GC values */
	for( i = 0; i < 101; i++)
	{
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
	}

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index]; i++)
	{
		gc_val = (int) round ( sonic_get_gc_content( params->this_sonic, params->this_sonic->chromosome_names[chr_index], i, i + WINDOWSLIDE));
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

int run_rd( bam_info** in_bam, parameters* params)
{
	/* Variables */
	int return_value, chr_index, bam_index, i;
	int chr_index_bam;

	bam1_t*	bam_alignment;
	bam1_core_t bam_alignment_core;

	if( ten_x_flag)
		fprintf( stderr,"10x flag: %d\n", ten_x_flag);

	for ( bam_index = 0; bam_index < params->num_bams; bam_index++)
	{
		/* Initialize the read depth(rd) array of each library for each bam file */
		init_rd( in_bam[bam_index], params);

		/* HTS implementation */
		in_bam[bam_index]->bam_file = safe_hts_open( params->bam_file_list[bam_index], "r");

		/* Read in BAM header information */
		in_bam[bam_index]->bam_header = bam_hdr_read( ( in_bam[bam_index]->bam_file->fp).bgzf);
		bam_alignment = bam_init1();

		/* Load the bam index file */
		in_bam[bam_index]->bam_file_index = sam_index_load( in_bam[bam_index]->bam_file, params->bam_file_list[bam_index]);
		if( in_bam[bam_index]->bam_file_index == NULL)
		{
			fprintf( stderr, "Error: Sam Index cannot be loaded (sam_index_load)\n");
			exit( 1);
		}

		for ( chr_index = params->first_chr; chr_index <= params->last_chr; chr_index++)
		{
			chr_index_bam = find_chr_index_bam( params->this_sonic->chromosome_names[chr_index], in_bam[bam_index]->bam_header);
			if( chr_index_bam == -1)
			{
				fprintf( stderr, "\nCannot find chromosome name %s in bam %s", params->this_sonic->chromosome_names[chr_index], in_bam[bam_index]->sample_name);
				continue;
			}

			fprintf( stderr, "\r                                                        ");
			fflush( stderr);
			fprintf( stderr, "\rReading BAM [%s] - Chromosome %s", in_bam[bam_index]->sample_name,
					in_bam[bam_index]->bam_header->target_name[chr_index_bam]);
			fflush( stderr);

			in_bam[bam_index]->iter = bam_itr_queryi( in_bam[bam_index]->bam_file_index, chr_index_bam, 0, in_bam[bam_index]->bam_header->target_len[chr_index_bam]);
			if( in_bam[bam_index]->iter == NULL)
			{
				fprintf( stderr, "Error: Iterator cannot be loaded (bam_itr_queryi)\n");
				exit( 1);
			}

			while( bam_itr_next( in_bam[bam_index]->bam_file, in_bam[bam_index]->iter, bam_alignment) > 0)
			{
				bam_alignment_core = bam_alignment->core;

				/* Increase the read count and depth */
				in_bam[bam_index]->read_count++;
				in_bam[bam_index]->read_depth[chr_index][bam_alignment_core.pos]++;
			}
		}

		fflush( stderr);
		fprintf( stderr, "\n");
		vh_logInfo( "Finished reading the bam file");

		/* Mean value (mu) calculation */
		calc_mean( params, in_bam[bam_index]);

		/* Close the BAM file */
		return_value = hts_close( in_bam[bam_index]->bam_file);
		if( return_value != 0)
		{
			fprintf( stderr, "Error closing BAM file\n");
			exit( 1);
		}
		/* Free the bam related files */
		bam_itr_destroy( in_bam[bam_index]->iter);
		bam_hdr_destroy( in_bam[bam_index]->bam_header);
		hts_idx_destroy( in_bam[bam_index]->bam_file_index);
	}
	return RETURN_SUCCESS;
}
