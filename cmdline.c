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
	static int sensitive = 0, no_soft_clip = 0, debug = 0, no_interdup = 0;
	static int skip_mrfast = 0, quick = 0, ten_x = 0, output_hs = 0, alt_mapping = 0, resolved = 0;
	static int make_sonic = 0;
	static int load_sonic = 0;
	static int do_remap = 0;
	static int histogram_only = 0;
	char *mapping_qual = NULL, *rp_support = NULL, *cluster_of_read = NULL;

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
			{"hist-only"   , no_argument,    &histogram_only, 'p'},
			{"version", no_argument,         0, 'v'},
			{"bamlist",               required_argument,	 0, 'b'},
			{"force-read-length"    , required_argument,	 0, 'l'},
			{"out"    , required_argument,	 0, 'o'},
			{"make-sonic"    , required_argument,	 0, 'c'},
			{"sonic"    , required_argument,	 0, 's'},
			{"sonic-info"    , required_argument,	 0, 'n'},
			{"first-chr", required_argument, 0, FIRST_CHROM},
			{"last-chr", required_argument, 0, LAST_CHROM},
			{"rd-ratio", required_argument, 0, 'a'},
			{"mq", required_argument, 0, 'e'},
			{"rp", required_argument, 0, 'j'},
			{"read-cluster", required_argument, 0, 'k'},
			{"no-soft-clip", no_argument, &no_soft_clip, 1 },
			{"no-interdup", no_argument, &no_interdup, 1 },
			{"debug", no_argument, &debug, 1 },
			{"xa", no_argument, &alt_mapping, 1 },
			{"resolved", no_argument, &resolved, 1 },
			{"sensitive", no_argument, &sensitive,    1 },
			{"skip-mrfast", no_argument, &skip_mrfast,  1 },
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
			//parse_bam_list( &params);
			break;

		case 'i':
			if ( params->num_bams == MAX_BAMS)
			{
				fprintf( stderr, "Number of input BAMs exceeded the maximum value (%d). Exiting.\n", MAX_BAMS);
				exit( EXIT_MAXBAMS);
			}
			params->num_bams = parse_bam_array(params, optarg);
			//			set_str( &( params->bam_file_list[( params->num_bams)++]), optarg);
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
			params->first_chr = atoi(optarg);
			break;

		case LAST_CHROM:
			params->last_chr = atoi(optarg);
			break;

		case 'e':
			set_str( &( mapping_qual), optarg);
			break;

		case 'j':
			set_str( &( rp_support), optarg);
			break;

		case 'k':
			set_str( &( cluster_of_read), optarg);
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
	if ( quick && do_remap){
		fprintf(stderr, "Cannot run both in quick and remap mode. Resetting to default (quick mode).\n");
		do_remap = 0;
	}
	else if ( do_remap){
		quick = 0;
	}

	/* histogram only mode */
	if ( histogram_only)
	  params->histogram_only = 1;
	  
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

	if ( params->bam_list_path != NULL)
	  params->num_bams = parse_bam_list( &params);

	/* check if outprefix is given */
	if( params->outprefix == NULL && !make_sonic)
	  {
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the output file name prefix using the --out option.\n");
		char *tmp_output_prefix = (char *) malloc(sizeof (char) * (strlen(params->bam_file_list[0]) + strlen("-output") + 2));
		sprintf( tmp_output_prefix, "%s-output", params->bam_file_list[0]);
		set_str( &( params->outprefix), tmp_output_prefix);
		free( tmp_output_prefix);
		//return EXIT_PARAM_ERROR;
	}



	/* check forced read length to be a positive integer or zero */
	if( params->force_read_length < 0)
	{
		fprintf( stderr, "[TARDIS CMDLINE WARNING] Invalid forced read length (%d). Resetted to 0 (disabled).\n", params->force_read_length);
		params->force_read_length = 0;
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

	if( cluster_of_read == NULL)
		params->cluster_of_read = 10;
	else
	{
		params->cluster_of_read = atoi(cluster_of_read);
		free( cluster_of_read);
	}
	cluster_of_reads = params->cluster_of_read;

	/* set flags */
	params->no_soft_clip = no_soft_clip;
	params->skip_mrfast = skip_mrfast;
	params->quick = quick; 
	params->ten_x = ten_x;
	params->output_hs = output_hs | ten_x;
	params->make_sonic = make_sonic;
	params->sensitive = sensitive;
	params->number_of_different_mei_types = count_mei_columns( params->mei);
	params->alt_mapping = alt_mapping;
	params->no_interdup = no_interdup;
	params->seq_resolved = resolved;

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


	get_working_directory(params);

	fprintf(stderr, "[TARDIS INFO] Working directory: %s\n", params->outdir); 
	
	return RETURN_SUCCESS;

}

void print_help( void)
{  
	fprintf( stdout, "\nTARDIS: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", TARDIS_VERSION, TARDIS_UPDATE, BUILD_DATE);
	fprintf( stdout, "\tBasic Parameters:\n\n");
	fprintf( stdout, "\t--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.\n");
	fprintf( stdout, "\t--input/-i [BAM files]     : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.\n");
	fprintf( stdout, "\t--out   [output prefix]    : Prefix for the output file names.\n");
	fprintf( stdout, "\t--ref   [reference genome] : Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--sonic [sonic file]       : SONIC file that contains assembly annotations.\n");
	fprintf( stdout, "\t--hist-only                : Generate fragment size histograms only, then quit.\n\n");
	fprintf( stdout, "\tAdvanced Parameters:\n\n");
	fprintf( stdout, "\t--read-count [int]         : # of clusters that a specific read can be involved in (Default is 10).\n");
	fprintf( stdout, "\t--rp   [int]               : Minimum number of supporting read pairs in initial clustering (Default is 5).\n");
	fprintf( stdout, "\t--mei   [\"Alu:L1:SVA\"]     : List of mobile element names.\n");
	fprintf( stdout, "\t--no-soft-clip             : Skip soft clip remapping.\n");
	fprintf( stdout, "\t--no-interdup              : Skip interspersed duplication clustering.\n");
	fprintf( stdout, "\t--resolved                 : Output sequence resolved vcf calls.\n");
	fprintf( stdout, "\t--xa                       : Look for the alternative mapping locations in BWA.\n");
	fprintf( stdout, "\t--first-chr [chr_index]	   : Start running from a specific chromosome [0-based index in reference file].\n");
	fprintf( stdout, "\t--last-chr [chr_index]	   : Run up to a specific chromosome [0-based index in reference file].\n");

	fprintf( stdout, "\n\tAdditional parameters for sensitive mode:\n\n");
	fprintf( stdout, "\t--sensitive                : Sensitive mode that uses all map locations. Requires mrFAST remapping.\n");
	fprintf( stdout, "\t--skip-mrfast              : Skip mrFAST mapping. Use this only if you already have the correct divet file. Sensitive mode only\n");
	fprintf( stdout, "\t--threads                  : Number of threads for mrFAST to remap discordant reads.\n");

	fprintf( stdout, "\n\tAdditional parameters to build SONIC file within TARDIS:\n\n");
	fprintf( stdout, "\t--make-sonic [sonic file]  : SONIC file that will contain the assembly annotations.\n");
	fprintf( stdout, "\t--sonic-info [\"string\"]    : SONIC information string to be used as the reference genome name.\n");
	fprintf( stdout, "\t--gaps  [gaps file]        : Assembly gap coordinates in BED3 format.\n");
	fprintf( stdout, "\t--dups  [dups file]        : Segmental duplication coordinates in BED3 format.\n"); 
	fprintf( stdout, "\t--reps  [reps file]        : RepeatMasker annotation coordinates in RepeatMasker format. See manual for details.\n");

	/* hidden
	fprintf( stdout, "\n\tAdditional parameters for 10X Genomics Linked Reads (under development):\n\n");
	fprintf( stdout, "\t--10x                      : Enable 10X Genomics Linked Reads mode.\n");
	fprintf( stdout, "\t--output-hs                : Output the selected clusters homogeneity scores to the VCF file.");
	*/
	
	fprintf( stdout, "\n\n\tInformation:\n");
	fprintf( stdout, "\t--version                  : Print version and exit.\n");
	fprintf( stdout, "\t--help                     : Print this help screen and exit.\n\n");
	fprintf( stdout, "It is bigger on the inside!\n\n");
}

int parse_bam_list( parameters** params)
{
	FILE* bam_list;
	char next_path[1024];
	int i;

	fprintf( stderr, "Parsing BAM/CRAM list file: %s ", ( *params)->bam_list_path);
	fflush( stderr);
	bam_list = safe_fopen( ( *params)->bam_list_path, "r");

	i = 0;
	while( fscanf( bam_list, "%s\n", next_path) == 1)
	{
		if ( i == MAX_BAMS){		  
		  fprintf( stderr, "Number of input BAMs exceeded the maximum value (%d). Exiting.\n", MAX_BAMS);
		  exit( EXIT_MAXBAMS);
		}
		set_str( &( ( *params)->bam_file_list)[i], next_path);
		i = i + 1;
	}

	fclose( bam_list);
	fprintf( stderr, " %d input files found.\n", i);
	return i;

}

int parse_bam_array( parameters* params, char *optarg)
{
  /* when the BAM/CRAM list is comma separated */
	char *next_bam;
	int i;

	i = params->num_bams;

	next_bam = strtok(optarg, ",");
	while( next_bam != NULL)
	{
	  printf ("Input file %s\n", next_bam);
	        if ( i == MAX_BAMS){		  
		  fprintf( stderr, "Number of input BAMs exceeded the maximum value (%d). Exiting.\n", MAX_BAMS);
		  exit( EXIT_MAXBAMS);
		}
		set_str( &( params->bam_file_list[i]), next_bam);
		i = i + 1;
		next_bam = strtok(NULL, ",");
	}

	fprintf( stderr, " %d input files found.\n", i);
	return i;
}
