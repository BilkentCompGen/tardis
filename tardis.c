
#include <stdio.h>
#include <time.h>
#include "common.h"
#include "cmdline.h"
#include "processbam.h"
#include "processfq.h"
#include "config.h"
#include "external.h"
#include "vh.h"
#include "tardis.h"
#include "bamonly.h"


FILE *logFile = NULL;

int main( int argc, char** argv)
{
	bam_info** in_bams;
	parameters* params;
	configuration* cfg;
	ref_genome* refgen;
	int return_value;
	char username[MAX_SEQ];
	int i;
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime);
	timeinfo = localtime( &rawtime);

	print_quote();

	/* Keeping simple logs in tardis.log file */
	logFile = safe_fopen ("tardis.log", "w");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n", timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday);

	/* Load configuration file (created if it does not exist) */
	cfg = ( configuration*) getMem( sizeof( configuration));
	load_config( cfg);

	/* Set program parameters */
	init_params( &params);

	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params);

	if( return_value == 0)
	{
		exit(EXIT_SUCCESS);
	}
	else if( return_value != 1)
	{
		exit( return_value);
	}

	if ( TARDIS_DEBUG)
	{
		print_params( params);
	}

	print_quote();


	/* make_sonic is standalone. Execute and return.  TODO
	if ( params->make_sonic)
	{
		return make_sonic(params);
		} */


	/* Load reference genome properties*/
	return_value = load_refgen( &refgen, params);

	if (return_value != RETURN_SUCCESS){
		exit (EXIT_SONIC);
	}

	if (params->last_chrom < params->first_chrom){
		params->last_chrom = refgen->chrom_count - 1;
	}

	/* Read BAM files and calculate the median/avg/std of fragment sizes per library */
	in_bams = ( bam_info**) getMem( sizeof( bam_info*));
	for( i = 0; i < params->num_bams; i++)
	{
		in_bams[i] = ( bam_info*) getMem( sizeof( bam_info));
		in_bams[i]->sample_name = NULL;
		load_bam( in_bams[i], params->bam_file_list[i]);
	}

	//passing the flags to VH
	ten_x_flag = params->ten_x;
	output_hs_flag = params->output_hs;

	print_quote();
	fprintf( stderr,"\nStructural Variation Detection\n");
	fprintf( stderr, "\n\tRun. Run, you clever boy... And remember.\n");

	/* Sensitive Mode */
	if( running_mode == SENSITIVE)
	{
		fprintf( stderr, "\nTardis is running in Sensitive Mode - using mrFast...\n\n");
		fprintf( logFile, "(Running in sensitive mode - using mrFast)\n\n");

		/* Extract FASTQs of discordants, OEAs, and orphans */
		for( i = 0; i < params->num_bams; i++){
			create_fastq( in_bams[i], params->bam_file_list[i], params);
		}
		fprintf( stderr, "All FASTQ files ready for remapping.\n");


		/* Remap with mrFAST */
		if( params->skip_remap == 0)
		{
			return_value = remap_mrfast( params, in_bams, cfg);

			if( return_value != RETURN_SUCCESS)
				return EXIT_EXTERNAL_PROG_ERROR;
		}
		else
			/* TODO: check if the remapping output indeed exists, so it is safe to skip */
			fprintf( stderr, "Skipping remapping step.\n");

		/* Read depth calculation */
		fprintf( logFile,"\n--> Read Depth Method\n\n");
		fprintf( stderr, "\nRunning Read Depth Method...\n");
		return_value = run_rd( in_bams, params, refgen);
		if( return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;

		/* Read pair calculation */
		fprintf( logFile,"--> Read Pair Method\n\n");
		fprintf( stderr, "\nRunning Read Pair Method (VariationHunter/CommonLAW)...\n");
		return_value = run_vh( refgen, params, in_bams);
		if( return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;

		free_all( in_bams, params, refgen);
	}
	/* Quick Mode */
	else
	{
		fprintf( stderr, "\nTardis is running in Quick Mode - using BWA...\n\n");
		fprintf( logFile, "(Running in quick mode - using BWA)\n\n");
		return_value = bamonly_run( refgen, params, in_bams);
		if ( return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;
	}

	getlogin_r( username, MAX_SEQ);
	fprintf( stderr, "\n%s, before I go, I just want to tell you: you were fantastic. Absolutely fantastic. And you know what? So was I.\n", username);

	fclose( logFile);

	return EXIT_SUCCESS;
}
