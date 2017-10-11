#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tardis.h"
#include "cmdline.h"
#define FIRST_CHROM 10001
#define LAST_CHROM 10002

int count_mei_columns( char* mei_string)
{
	char *tok;
	char str[1024];
	int mei_code;

	strcpy(str, mei_string);

	mei_code = 0;
	tok = strtok(str, ":");

	while (tok != NULL)
	{
		mei_code++;
		tok = strtok(NULL, ":");
	}

	return mei_code;
}

int parse_command_line( int argc, char** argv, parameters* params)
{
	int index;
	int o;
	static int run_vh = 0, run_ns = 0, run_rd = 0, run_sr = 0, run_all = 1, sensitive = 0, no_soft_clip = 0, debug = 0;
	static int skip_fastq = 0, skip_sort = 0, skip_remap = 0, skip_cluster = 0, quick = 0, ten_x = 0, output_hs = 0;
	static int make_sonic = 0;
	static int load_sonic = 0;
	static int do_remap = 0;
	char *threshold = NULL, *mapping_qual = NULL, *rp_support = NULL, *alt_mapping = NULL;

	static struct option long_options[] = 
	{
			{"input"  , required_argument,   0, 'i'},
			{"ref"    , required_argument,   0, 'f'},
			{"gaps"   , required_argument,   0, 'g'},
			{"dups"   , required_argument,   0, 'd'},
			{"reps"   , required_argument,   0, 'r'},
			{"mei"    , required_argument,   0, 'm'},
			{"threads", required_argument,   0, 't'},
			{"help"   , no_argument,         0, 'h'},
			{"version", no_argument,         0, 'v'},
			{"bamlist",               no_argument,	 0, 'b'},
			{"force-read-length"    , required_argument,	 0, 'l'},
			{"out"    , required_argument,	 0, 'o'},
			{"make-sonic"    , required_argument,	 0, 'c'},
			{"sonic"    , required_argument,	 0, 's'},
			{"sonic-info"    , required_argument,	 0, 'n'},
			{"first-chrom", required_argument, 0, FIRST_CHROM},
			{"last-chrom", required_argument, 0, LAST_CHROM},
			{"rd-ratio", required_argument, 0, 'a'},
			{"mq", required_argument, 0, 'e'},
			{"rp", required_argument, 0, 'j'},
			{"vh"     , no_argument, &run_vh,     1 },
			{"no-soft-clip" , no_argument, &no_soft_clip, 1 },
			{"debug" , no_argument, &debug, 1 },
			{"xa" , required_argument, 0, 'k' },
			/*
			{"rd"     , no_argument, &run_rd,     1 },
			{"ns"     , no_argument, &run_ns,     1 },
			{"sr"     , no_argument, &run_sr,     1 },
			{"all"    , no_argument, &run_all,    1 },
			 */
			{"sensitive"    , no_argument, &sensitive,    1 },
			{"skip-fastq", no_argument, &skip_fastq,  1 },
			{"skip-sort" , no_argument, &skip_sort,  1 },
			{"skip-remap" , no_argument, &skip_remap,  1 },
			{"skip-cluster" , no_argument, &skip_cluster,  1 },
			{"quick" , no_argument, &quick,  1 },
			{"remap" , no_argument, &do_remap,  1 },
			{"10x", no_argument, &ten_x, 1},
			{"output-hs", no_argument, &output_hs, 1},
			{0        , 0,                   0,  0 }
	};

	if( argc == 1)
	{
		print_help();
		return 0;
	}

	while( ( o = getopt_long( argc, argv, "hvb:i:f:g:d:r:o:m:c:s:a:e:n:j:k", long_options, &index)) != -1)
	{
		switch( o)
		{
		case 'b':
			set_str( &( params->bam_list_path), optarg);
			parse_bam_list( &params);
			break;

		case 'i':
			if ( params->num_bams == MAX_BAMS)
			{
				fprintf( stderr, "Number of input BAMs exceeded the maximum value (%d). Exiting.\n", MAX_BAMS);
				exit( EXIT_MAXBAMS);
			}
			set_str( &( params->bam_file_list[( params->num_bams)++]), optarg);
			break;

		case 'f':
			set_str( &( params->ref_genome), optarg);
			break;

		case 'g':
			set_str( &( params->gaps), optarg);
			break;

		case 's':
			set_str( &( params->sonic_file), optarg);
			load_sonic = 1;
			break;

		case 'c':
			set_str( &( params->sonic_file), optarg);
			make_sonic = 1;
			break;

		case 'd':
			set_str( &( params->dups), optarg);
			break;

		case 'r':
			set_str( &( params->reps), optarg);
			break;

		case 'm':
			set_str( &( params->mei), optarg);
			break;

		case 'n':
			set_str( &( params->sonic_info), optarg);
			break;

		case 'o':
			set_str( &( params->outprefix), optarg);
			break;

		case 't':
			params->threads = atoi( optarg);
			break;

		case 'l':
			params->force_read_length = atoi( optarg);
			break;

		case 'h':
			print_help();
			return 0;
			break;

		case FIRST_CHROM:
			params->first_chrom = atoi(optarg);
			break;

		case LAST_CHROM:
			params->last_chrom = atoi(optarg);
			break;

		case 'a':
			set_str( &( threshold), optarg);
			break;

		case 'e':
			set_str( &( mapping_qual), optarg);
			break;

		case 'j':
			set_str( &( rp_support), optarg);
			break;

		case 'k':
			set_str( &( alt_mapping), optarg);
			break;

		case 'v':
			fprintf( stderr, "\nTARDIS: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
			fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", TARDIS_VERSION, TARDIS_UPDATE, BUILD_DATE);
			fprintf( stderr, "It is bigger on the inside!\n\n");
			return 0;
			break;
		}
	}

	/* check quick vs remap mode */
	if (quick && do_remap){
		fprintf(stderr, "Cannot run both in quick and remap mode. Resetting to default (quick mode).\n");
		do_remap = 0;
	}
	else if (do_remap){
		quick = 0;
	}



	/* TODO: check parameter validity */

	/* check algorithms to run; run_all is the default */	

	/*
	if( !run_vh && !run_rd && !run_sr && !run_ns)
	{
		run_all = 1;
	}

	if( run_all)
	{
		run_vh = 1; run_rd=1; run_sr = 1; run_ns = 1;
		} */


	/* check if outprefix is given */
	if( params->outprefix == NULL && !make_sonic)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the output file name prefix using the --out option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --num-bams > 0 */
	if( params->num_bams <= 0 && params->bam_list_path == NULL && !make_sonic)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Invalid number of input BAM files was entered (%d).\n", params->num_bams);
		return EXIT_PARAM_ERROR;
	}

	/* check if --input is invoked */
	if( params->bam_file_list[0] == NULL && !make_sonic)
	{
		if( params->bam_list_path == NULL)
		{
			fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter list of input BAM files using the --input or --bamlist option.\n");
			return EXIT_PARAM_ERROR;
		}
	}

	/* check if --ref   is invoked */
	if( params->ref_genome == NULL)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter reference genome file (FASTA) using the --ref option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --gaps  is invoked */
	if( params->gaps == NULL && !load_sonic)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the assembly gaps file (BED) using the --gaps option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --reps  is invoked */
	if( params->reps == NULL && !load_sonic)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the repeats file (RepeaMasker) using the --reps option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --dups  is invoked */
	if( params->dups == NULL && !load_sonic)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the segmental duplications file (BED) using the --dups option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --sonic  is invoked */
	if( params->sonic_file == NULL && load_sonic)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the SONIC file (BED) using the --sonic option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --mei is invoked. If not, set Alu:L1HS:SVA as default */
	if( params->mei == NULL)
	{
		set_str( &( params->mei), "Alu:L1:SVA");
	}

	/* check if threads>0 */
	if( params->threads <= 0)
	{
		fprintf( stderr, "[TARDIS CMDLINE WARNING] Invalid number of threads was entered (%d). Resetted to 1.\n", params->threads);
		params->threads = 1;
	}

	/* check if both --input and --bamlist are invoked */
	if ( params->bam_list_path != NULL && params->bam_file_list[0] != NULL && !make_sonic)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please use either --input or --bamlist parameter. Not both!\n");
		return EXIT_PARAM_ERROR;	        
	}

	/* check forced read length to be a positive integer or zero */
	if( params->force_read_length < 0)
	{
		fprintf( stderr, "[TARDIS CMDLINE WARNING] Invalid forced read length (%d). Resetted to 0 (disabled).\n", params->force_read_length);
		params->force_read_length = 0;
	}

	if( threshold == NULL)
		params->rd_threshold = 2;
	else
	{
		params->rd_threshold = atoi(threshold);
		free( threshold);
	}

	if( mapping_qual == NULL)
		params->mq_threshold = 5;
	else
	{
		params->mq_threshold = atoi(mapping_qual);
		free( mapping_qual);
	}

	if( rp_support == NULL)
		params->rp_threshold = 5;
	else
	{
		params->rp_threshold = atoi(rp_support);
		free( rp_support);
	}

	if( alt_mapping == NULL)
		params->alt_mapping = 0;
	else
	{
		params->alt_mapping = atoi(alt_mapping);
		free( alt_mapping);
	}

	/* set flags */
	params->run_vh = run_vh;
	params->run_rd = run_rd;
	params->run_sr = run_sr;
	params->run_ns = run_ns;
	params->no_soft_clip = no_soft_clip;
	params->skip_fastq = skip_fastq;
	params->skip_sort = skip_sort;
	params->skip_remap = skip_remap;
	params->skip_vhcluster = skip_cluster;
	params->quick = quick; 
	params->ten_x = ten_x;
	params->output_hs = output_hs | ten_x;
	params->make_sonic = make_sonic;
	params->sensitive = sensitive;
	params->number_of_different_mei_types = count_mei_columns( params->mei);
	debug_mode = debug;

	if( debug_mode)
		fprintf(stderr, "\n\n*** DEBUG MODE is on, you can check .NAME and .CLUSTER files ***\n\n");

	if( params->sensitive == 0)
		params->quick = 1;

	if( params->quick)
		running_mode = QUICK;
	else
		running_mode = SENSITIVE;

	if (load_sonic)
		params->load_sonic = load_sonic;

	if ( params->sonic_info == NULL)
		set_str( &(params->sonic_info), params->ref_genome);


	return RETURN_SUCCESS;

}

void print_help( void)
{  
	fprintf( stdout, "\nTARDIS: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", TARDIS_VERSION, TARDIS_UPDATE, BUILD_DATE);
	fprintf( stdout, "\tParameters:\n\n");
	fprintf( stdout, "\t--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.\n");
	fprintf( stdout, "\t--input [BAM files]        : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.\n");
	fprintf( stdout, "\t--out   [output prefix]    : Prefix for the output file names.\n");
	fprintf( stdout, "\t--ref   [reference genome] : Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--sonic [sonic file]       : SONIC file that contains assembly annotations.\n");	
	fprintf( stdout, "\t--mei   [\"Alu:L1:SVA\"]   : List of mobile element names.\n");
	fprintf( stdout, "\t--no-soft-clip             : Skip soft clip remapping.\n");

	fprintf( stdout, "\n\tAdditional parameters for sensitive mode:\n\n");
	fprintf( stdout, "\t--sensitive                : Sensitive mode that uses all map locations. Requires mrFAST remapping.\n");
	fprintf( stdout, "\t--skip-fastq               : Skip FASTQ dump for discordants. Use this only if you are regenerating the calls. Sensitive mode only.\n");
	fprintf( stdout, "\t--skip-sort                : Skip FASTQ sort for discordants. Use this only if you are regenerating the calls. Sensitive mode only. \n");
	fprintf( stdout, "\t--skip-remap               : Skip FASTQ remapping for discordants. Use this only if you are regenerating the calls. Sensitive mode only\n");
	fprintf( stdout, "\t--threads                  : Number of threads for mrFAST to remap discordant reads.\n");

	fprintf( stdout, "\n\tAdditional parameters to build SONIC file within TARDIS:\n\n");
	fprintf( stdout, "\t--make-sonic [sonic file]  : SONIC file that will contain the assembly annotations.\n");
	fprintf( stdout, "\t--sonic-info [\"string\"]  : SONIC information string to be used as the reference genome name.\n");
	fprintf( stdout, "\t--gaps  [gaps file]        : Assembly gap coordinates in BED3 format.\n");
	fprintf( stdout, "\t--dups  [dups file]        : Segmental duplication coordinates in BED3 format.\n"); 
	fprintf( stdout, "\t--reps  [reps file]        : RepeatMasker annotation coordinates in RepeatMasker format. See manual for details.\n");


	fprintf( stdout, "\n\n\tInformation:\n");
	fprintf( stdout, "\t--version                  : Print version and exit.\n");
	fprintf( stdout, "\t--help                     : Print this help screen and exit.\n\n");
	fprintf( stdout, "It is bigger on the inside!\n\n");
}

void parse_bam_list( parameters** params)
{
	FILE* bam_list;
	char next_path[1024];
	int i;

	bam_list = safe_fopen( ( *params)->bam_list_path, "r");

	i = 0;
	while( fscanf( bam_list, "%s\n", next_path) == 1)
	{
		set_str( &( ( *params)->bam_file_list)[i], next_path);
		i = i + 1;
	}

	fclose( bam_list);
}
