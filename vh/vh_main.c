#include "vh_main.h"
#include <string.h>
#include <math.h>
#include "vh_createMaxClusterMEI.h"
#include "vh_divethandler.h"
#include "vh_maximalCluster.h"
#include "vh_setcover.h"
#include "../free.h"

FILE *fileOutput = NULL;

int vh_isValid ()
{
	//TODO: To be completed. Check for files to exists, etc.
	return 1;
}

void vh_quitProgram (int exitCode)
{
	if (exitCode != 0)
	{
		if (!strcmp (g_error_message, ""))
			strcpy (g_error_message, "Unexpected Error Found!\n");

		char errMsg[150];
		sprintf (errMsg, "Program Exited with Errors:%d", exitCode);
		vh_logError (errMsg);
		switch (exitCode)
		{
		case EXIT_CODE_ARG_ERROR:
			sprintf (errMsg, "Error in options: %s", g_error_message);
			vh_logError (errMsg);
			break;
		case EXIT_CODE_DIVET_ERROR:
			sprintf (errMsg, "Error in DIVET file: %s", g_error_message);
			vh_logError (errMsg);
			break;
		case EXIT_CODE_MEMORY_ERROR:
			sprintf (errMsg, "Memory Problem Occured: %s", g_error_message);
			vh_logError (errMsg);
			break;
		default:
			sprintf (errMsg, "Uncategorized Error Found: %s", g_error_message);
			vh_logError (errMsg);
			break;
		}
		exit (exitCode);
	}
	else
	{
		vh_logInfo ("Program Exited Successfully.\r\n");
		exit (0);
	}
}



void vh_init(bam_info** in_bams, parameters *params, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit)
{
	int totalNumUniqueReads = 0;
	int indexStart = 0;
	int count;
	int i,j,k;
	int total_lib_count = 0;
	struct LibraryInfo *newLibInfo, *cursor, *t;
	char **allReadNameList;
	FILE *fileOutputReadName;

	int num_bams = params->num_bams;

	g_maxListBrkPointIntr = MAXLISTBRKPOINTINTR;

	/*Hold the libraries in a linked list*/
	g_libInfo = NULL;


	
	for (k = 0; k < num_bams; k++)
	{
		for (i = 0; i < in_bams[k]->num_libraries; i++)
		{
			newLibInfo = (struct LibraryInfo *) getMem (sizeof (struct LibraryInfo));
			strcpy (newLibInfo->libName, in_bams[k]->libraries[i]->libname);
			strcpy (newLibInfo->indName, in_bams[k]->sample_name);
			strcpy (newLibInfo->libFileAdrs, in_bams[k]->libraries[i]->divet);
			newLibInfo->libId = total_lib_count++;
			newLibInfo->minDelta = in_bams[k]->libraries[i]->conc_min - 2 * in_bams[k]->libraries[i]->read_length;
			newLibInfo->maxDelta = in_bams[k]->libraries[i]->conc_max - 2 * in_bams[k]->libraries[i]->read_length;
			if (newLibInfo->minDelta < 0)
				newLibInfo->minDelta = 0;
			if (newLibInfo->maxDelta < 0)
				newLibInfo->maxDelta = 0;
			newLibInfo->readLen = in_bams[k]->libraries[i]->read_length;

			fprintf( logFile,"Min and Max Delta for %s = %d - %d\n", newLibInfo->libName, newLibInfo->minDelta, newLibInfo->maxDelta);

			/* We store the reads in hash[] based on the hash values of read names */
			newLibInfo->hash =(struct ReadName **) getMem (NHASH * sizeof (struct ReadName *));
			for (j = 0; j < NHASH; j++)
				newLibInfo->hash[j] = NULL;

			newLibInfo->head = NULL;
			newLibInfo->tail = NULL;
			newLibInfo->size = 0;
			newLibInfo->next = NULL;

			if (g_libInfo == NULL)
			{
				g_libInfo = newLibInfo;
			}
			else                  //add to the end of the linked list
			{
				for (t = g_libInfo; t->next != NULL; t = t->next)
					;        //Skip till end of LinkedList
				t->next = newLibInfo;
			}
		}
	}
	cursor = g_libInfo;

	if( debug_mode)
	{
		fileOutput = safe_fopen (outputFile, "w");
		fileOutputReadName = safe_fopen (outputRead, "w");
	}

	for (; cursor; cursor = cursor->next)
	{
		vh_logInfo( "Reading Divet Files ...");
		vh_loadDivetFile( cursor, params->this_sonic);
		sprintf( g_loggerMsgBuffer, "%d rows loaded successfully.",cursor->size);
		vh_logInfo( g_loggerMsgBuffer);

		vh_pruneAndNormalizeDivets( cursor, preProsPrune, overMapLimit);
		sprintf( g_loggerMsgBuffer, "%d rows after pruning.", cursor->size);
		vh_logInfo( g_loggerMsgBuffer);
	}

	vh_logInfo ("Writing ReadName Sorted");
	cursor = g_libInfo;

	while (cursor != NULL)
	{
		totalNumUniqueReads = totalNumUniqueReads + vh_countNumReads (cursor->hash);
		cursor = cursor->next;
	}
	cursor = g_libInfo;
	allReadNameList = (char **) getMem (totalNumUniqueReads * sizeof (char *));

	/* Put the names of all the reads in allReadNameList and write them to .name file in sorted order */
	for (; cursor; cursor = cursor->next)
	{
		indexStart = vh_exportToArray(cursor->hash, allReadNameList, indexStart);
	}

	qsort (allReadNameList, totalNumUniqueReads, sizeof (char *),vh_cmprReadNameStr);
	if( debug_mode)
		fprintf (fileOutputReadName, "%i\n", totalNumUniqueReads);

	fprintf( logFile, "Unique Read Count = %d\n", totalNumUniqueReads);

	/* Write the read names to read_names structure for set_cover */
	read_names = (readEl *) getMem( ( totalNumUniqueReads + 1) * sizeof( readEl));

	for (count = 0; count < totalNumUniqueReads; count++)
	{
		if( debug_mode)
			fprintf( fileOutputReadName, "%s\n", allReadNameList[count]);

		read_names[count].readName = NULL;
		set_str( &(read_names[count].readName), allReadNameList[count]);
		read_names[count].readCovered = 0;
		read_names[count].readId = count;
		read_names[count].libId = -1;
		read_names[count].indId = -1;
		read_names[count].next = NULL;
	}
	read_names_count = totalNumUniqueReads;

	if( debug_mode)
		fclose (fileOutputReadName);
}

void vh_clustering (bam_info** in_bams, parameters *params, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit)
{
	int chr_index, return_value, invdup_location, interdup_location;
	int mei_count;
	struct LibraryInfo *cursor, *t;

	/* Initialization function */
	vh_init( in_bams, params, preProsPrune, outputFile, outputRead, overMapLimit);

	for( chr_index = params->first_chr; chr_index <= params->last_chr; chr_index++)
	{

		fprintf(stderr, "\n                                                        ");
		fflush(stderr);
		fprintf(stderr, "\nProcessing chromosome %s\n", params->this_sonic->chromosome_names[chr_index]);
		fflush(stderr);

		/* Deletion */
		fprintf( stderr, "\nPreparing Deletion clusters");
		vh_initializeReadMapping_Deletion( params->this_sonic, chr_index);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_createDeletionClusters( params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);

		/* Inversion */
		fprintf( stderr, "\nPreparing Inversion clusters");
		vh_initializeReadMapping_Inversion( params->this_sonic, chr_index);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_createInversionClusters( params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);

		/* Insertion */
		fprintf( stderr, "\nPreparing Insertion clusters");
		vh_initializeReadMapping_Insertion( params->this_sonic, chr_index);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_createInsertionClusters( params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);

		/* Tandem Duplication */
		fprintf( stderr, "\nPreparing Tandem Duplication clusters");
		vh_initializeReadMapping_TDup( params->this_sonic, chr_index);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_createTDupClusters( params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);

		/* Interspersed Direct Duplication
		fprintf( stderr, "\nPreparing Interspersed Duplication clusters");
		for( interdup_location = 0; interdup_location <= RIGHTSIDE; interdup_location++)
		{
			vh_initializeReadMapping_InterDup( params->this_sonic, chr_index, interdup_location);
			fprintf( stderr, ".");
			fflush( stderr);
			vh_createInterDupClusters( params->this_sonic->chromosome_lengths[chr_index], interdup_location);
			fprintf( stderr, ".");
			fflush( stderr);
			vh_finalizeReadMapping_InterDup( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
			fprintf( stderr, ".");
			fflush( stderr);
		}

		/* Interspersed Inverted Duplication
		fprintf( stderr, "\nPreparing Interspersed Duplication (Inverted) clusters");
		for( invdup_location = 0; invdup_location <= RIGHTSIDE; invdup_location++)
		{
			vh_initializeReadMapping_InvDup( params->this_sonic, chr_index, invdup_location);
			fprintf( stderr, ".");
			fflush( stderr);
			vh_createInvDupClusters( params->this_sonic->chromosome_lengths[chr_index], invdup_location);
			fprintf( stderr, ".");
			fflush( stderr);
			vh_finalizeReadMapping_InvDup( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
			fprintf( stderr, ".");
			fflush( stderr);
		}*/

		/* Mei */
		fprintf( stderr, "\nPreparing MEI clusters");
		initializeReadMapping_MEI( in_bams, params, chr_index);
		fprintf(stderr, ".");
		fflush(stderr);
		MEICluster_Region( params, chr_index);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping_Mei( params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(stderr, ".");
		fflush(stderr);
	}
	fprintf(stderr, "\n");
	if( debug_mode)
		fclose (fileOutput);

	free_libraries();
}

int run_vh( parameters *params, bam_info ** in_bams)
{
	int i, j;
	double preProsPrune = 0.001;
	int overMapLimit = 500;
	char divetfile[MAX_SEQ];
	char outputfile[MAX_SEQ];
	char outputread[MAX_SEQ];
	char svfile[MAX_SEQ] = {};

	/* Print all structural variations in .vcf format */
	sprintf( svfile, "%s.vcf", params->outprefix);
	FILE* fpVcf = safe_fopen( svfile, "w");

	if( debug_mode)
	{
		sprintf( outputfile,"%s.clusters", params->outprefix);
		sprintf( outputread,"%s.name", params->outprefix);
	}

	for( i = 0; i < params->num_bams; i++)
	{
		for ( j = 0; j < in_bams[i]->num_libraries; j++)
		{
		        sprintf( divetfile, "%s/%s-%s.sam_DIVET.vh", params->outdir, in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			set_str( &( in_bams[i]->libraries[j]->divet), divetfile);
		}
	}
	print_vcf_header( fpVcf, in_bams, params);

	vh_logInfo( "Calculating maximal clusters.");
	vh_clustering( in_bams, params, preProsPrune, outputfile, outputread, overMapLimit);

	vh_logInfo( "Applying SET-COVER approximation to find putative structural variation.");
	vh_setcover( in_bams, params, fpVcf);

	print_sv_stats();

	fprintf( stderr, "\nTARDIS is complete. Found %d SVs. Results are in the %s file.", sv_count, svfile);

	return RETURN_SUCCESS;
}

