#include <stdio.h>
#include <time.h>
#include "common.h"
#include "cmdline.h"
#include "processbam.h"
#include "processfq.h"
#include "config.h"
#include "external.h"
#include "tardis.h"
#include "bamonly.h"
#include "variants.h"
#include "splitread.h"
#include "vh_main.h"
#include "sonic/sonic.h"
#include "processrefgen.h"
#include "free.h"

FILE *logFile = NULL;

int main( int argc, char** argv)
{
	bam_info** in_bams;
	parameters* params;
	configuration* cfg;
	int return_value;
	char* username;
	int i;
	time_t rawtime;
	struct tm * timeinfo;
	char *log_file_path;
	
	time ( &rawtime);
	timeinfo = localtime( &rawtime);

	print_quote();


	/* Set program parameters */
	init_params( &params);
	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params);


	if( return_value == 0){
	  /*        fprintf(stderr, "parse_command_line returned 0.\n"); */
		exit(EXIT_SUCCESS);
	}
	else if( return_value != 1){
	  /*      fprintf(stderr, "parse_command_line returned %d.\n", return_value); */
		exit( return_value);
	}

	/* Keeping simple logs in tardis.log file */
	log_file_path = (char *) getMem(sizeof(char) * (12+strlen(params->outprefix)));
	sprintf(log_file_path, "%s-%s", params->outprefix, "tardis.log");
	logFile = safe_fopen (log_file_path, "w");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday);
	free(log_file_path);

	/* Load configuration file (created if it does not exist) */
	cfg = ( configuration*) getMem( sizeof( configuration));
	load_config( cfg, params);

	/* make_sonic is standalone. Execute and return.  */
	if ( params->make_sonic)
		return sonic_build(params->ref_genome, params->gaps, params->reps, params->dups, params->sonic_info, params->sonic_file);

	/* Load SONIC */
	params->this_sonic = sonic_load(params->sonic_file);

	if (params->last_chr < params->first_chr)
		params->last_chr = params->this_sonic->number_of_chromosomes - 1;

	if ( TARDIS_DEBUG)
		print_params( params);

	print_quote();

	/* Read BAM files and calculate the median/avg/std of fragment sizes per library */
	in_bams = ( bam_info**) getMem( (params->num_bams) * sizeof( bam_info*));
	for( i = 0; i < params->num_bams; i++)
	{
		in_bams[i] = ( bam_info*) getMem( sizeof( bam_info));
		in_bams[i]->sample_name = NULL;
		load_bam( params, in_bams[i], params->bam_file_list[i], params->alt_mapping, i, params->ref_genome);
	}

	/* Passing the flags to VH */
	ten_x_flag = params->ten_x;
	output_hs_flag = params->output_hs;

	print_quote();
	fprintf( stderr, "\n\tRun. Run, you clever boy... And remember.\n");

	if ( !params->no_soft_clip)
	  init_hash_count( params);
	
	/* Sensitive Mode */
	if( running_mode == SENSITIVE)
	{
	  fprintf( stderr, "\nTARDIS (v%s) is running in Sensitive Mode - using mrFAST...\n\n", TARDIS_VERSION);
		fprintf( logFile, "(Running in sensitive mode - using mrFAST)\n\n");

		/* If you already have the correct divet file */
		if( params->skip_mrfast == 0)
		{
			/* Extract FASTQs of discordants, OEAs, and orphans */
			for( i = 0; i < params->num_bams; i++)
				create_fastq( in_bams[i], params->bam_file_list[i], params);

			fprintf( stderr, "All FASTQ files ready for remapping.\n");

			/* Remap with mrFAST */
			return_value = remap_mrfast( params, in_bams, cfg);
			if( return_value != RETURN_SUCCESS)
				return EXIT_EXTERNAL_PROG_ERROR;
		}
		/* Read depth calculation */
		return_value = run_rd( in_bams, params);
		if( return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;

		/* Read pair calculation */
		return_value = run_vh( params, in_bams);
		if( return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;

		clean_up_temp_files(params);
		free_sensitive( in_bams, params);
	}
	/* Quick Mode */
	else
	{
	        fprintf( stderr, "\nTARDIS (v%s) is running in Quick Mode\n\n", TARDIS_VERSION);
		fprintf( logFile, "(Running in quick mode)\n\n");
		return_value = bamonly_run( params, in_bams);
		if ( return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;

		clean_up_temp_files(params);
		free_quick( in_bams, params);
	}

	username = ( char*) getMem( MAX_SEQ * sizeof( char));
	getlogin_r( username, (MAX_SEQ - 1));
	fprintf( stderr, "\n%s, before I go, I just want to tell you: you were fantastic. Absolutely fantastic. And you know what? So was I.\n", username);

	fclose( logFile);

	free( username);
	return EXIT_SUCCESS;
}
