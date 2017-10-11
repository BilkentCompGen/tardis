#include "vh_main.h"
#include <string.h>
#include <math.h>
#include "vh_createMaxClusterMEI.h"
#include "vh_divethandler.h"
#include "vh_maximalCluster.h"
#include "vh_intervalhandler.h"
#include "vh_setcover.h"

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


void vh_pruneAndNormalizeDivets (struct LibraryInfo *lib, double preProsPrune, int overMapLimit)
{
	struct DivetRow *cursor = lib->head->next;
	if (cursor == NULL)
		return;

	//We are not checking the first and last one
	while (cursor != NULL && cursor->next != NULL)
	{
		if (((cursor->next->phredScore / cursor->next->readName->sumPhredValue) < preProsPrune)
				|| (cursor->next->readName->occurrences > overMapLimit))
		{
			struct DivetRow *toBeDeleted = cursor->next;
			cursor->next = cursor->next->next;
			lib->size--;
			free (toBeDeleted->chromosome_name);
			free (toBeDeleted);
		}
		else
		{
			cursor->next->phredScore = cursor->next->phredScore / cursor->next->readName->sumPhredValue;
			cursor = cursor->next;
		}
	}
}

void vh_init(bam_info** in_bams, ref_genome* ref, parameters *params, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit)
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

			fprintf(logFile,"Min Delta for %s = %d\nMax Delta for %s = %d\n",newLibInfo->libName, newLibInfo->minDelta, newLibInfo->libName, newLibInfo->maxDelta);

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
		vh_logInfo ("Reading Divet Files ...");
		vh_loadDivetFile (cursor, params->this_sonic);
		sprintf (g_loggerMsgBuffer, "%d rows loaded successfully.",cursor->size);
		vh_logInfo (g_loggerMsgBuffer);
		fprintf(logFile,"There are %d divet rows - ",cursor->size);

		vh_pruneAndNormalizeDivets (cursor, preProsPrune, overMapLimit);
		sprintf (g_loggerMsgBuffer, "%d rows after pruning.", cursor->size);
		vh_logInfo (g_loggerMsgBuffer);
		fprintf(logFile,"%d left after pruning\n",cursor->size);
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

	fprintf (logFile, "There are %d unique reads \n", totalNumUniqueReads);

	/* Write the read names to read_names structure for set_cover */
	read_names = (readEl *) getMem( ( totalNumUniqueReads + 1) * sizeof( readEl));

	for (count = 0; count < totalNumUniqueReads; count++)
	{
		if( debug_mode)
			fprintf (fileOutputReadName, "%s\n", allReadNameList[count]);

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

	//TODO: Someone should free the memory allocated for the divets and for the libraryinfos
	//free(newLibInfo);
}

void vh_clustering (bam_info** in_bams, ref_genome* ref, parameters *params, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit)
{
	int i, return_value;
	int mei_count;
	struct LibraryInfo *cursor, *t;

	/* Initialization function */
	vh_init(in_bams, ref, params, preProsPrune, outputFile, outputRead, overMapLimit);

	/* MEI Filtering */
	mei_count = mei_filtering( ref, params);
	fprintf(logFile,"%d mobile elements filtered\n", mei_count);
	sprintf (g_loggerMsgBuffer, "%d mobile elements filtered.", mei_count);
	vh_logInfo (g_loggerMsgBuffer);

	for (i = 0; i < ref->chrom_count; i++)
	{
		if(ref->in_bam[i] == false)
			continue;

		fprintf(stderr, "\r                                                        ");
		fflush(stderr);
		fprintf(stderr, "\rProcessing chromosome %s", ref->chrom_names[i]);
		fflush(stderr);

		/* Deletion */
		vh_initializeReadMapping_Deletion (ref->chrom_names[i], ref->chrom_lengths[i], params->this_sonic);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_createDeletionClusters (ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping (ref->chrom_names[i], ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);

		/* Inversion */
		vh_initializeReadMapping_Inversion (ref->chrom_names[i], ref->chrom_lengths[i], params->this_sonic);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_createInversionClusters (ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping (ref->chrom_names[i], ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);

		/* Insertion */
		vh_initializeReadMapping_Insertion (ref->chrom_names[i], ref->chrom_lengths[i], params->this_sonic);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_createInsertionClusters (ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping (ref->chrom_names[i], ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);

		/* Tandem Duplication 
		vh_initializeReadMapping_TDup (ref->chrom_names[i], ref->chrom_lengths[i], params->this_sonic);

		fprintf(stderr, ".");
		fflush(stderr);
		vh_createTDupClusters (ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping (ref->chrom_names[i], ref->chrom_lengths[i]);
		fprintf(stderr, ".");

		fflush(stderr); */

		/* Mei */
		initializeReadMapping_MEI (in_bams, params, ref->chrom_names[i], ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);
		MEICluster_Region(params, ref->chrom_names[i], ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);
		vh_finalizeReadMapping_Mei(ref->chrom_lengths[i]);
		fprintf(stderr, ".");
		fflush(stderr);
	}
	fprintf(stderr, "\n");
	if( debug_mode)
		fclose (fileOutput);

	/* Free g_libInfo and hash */
	cursor = g_libInfo;
	t = cursor;
	while (t != NULL)
	{
		cursor = cursor->next;
		for (i = 0; i < NHASH; i++)
			if (t->hash != NULL)
				if (t->hash[i] != NULL)
					free(t->hash[i]);
		if (t->hash != NULL)
			free(t->hash);
		if (t != NULL)
			free(t);
		t = cursor;
	}
}

int run_vh( ref_genome* ref, parameters *params, bam_info ** in_bams)
{
	int i, j;
	double preProsPrune= 0.001;
	int overMapLimit=500;
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
			sprintf( divetfile, "%s-%s.sam_DIVET.vh",  in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			set_str( &( in_bams[i]->libraries[j]->divet), divetfile);
		}
	}
	print_vcf_header(fpVcf, in_bams, params);

	vh_logInfo( "Calculating maximal clusters.");
	if ( !params->skip_vhcluster) // this parameter is only intended for debugging purposes. End users shouldn't use this
		vh_clustering( in_bams, ref, params, preProsPrune, outputfile, outputread, overMapLimit);

	vh_logInfo( "Applying SET-COVER approximation to find putative structural variation.");
	vh_setcover( in_bams, params, ref, fpVcf);

	if( debug_mode)
		fprintf( stderr, "\nTARDIS is complete. Found %d SVs and %d LowQual SVs. Results are in the %s file.", sv_count, sv_lowqual_count, svfile);
	else
		fprintf( stderr, "\nTARDIS is complete. Found %d SVs. Results are in the %s file.", sv_count, svfile);



	return RETURN_SUCCESS;
}

