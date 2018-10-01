#include <time.h>
#include <math.h>
#include <ctype.h>		//For isspace
#include "vh_divethandler.h"
#include "vh_main.h"
#include "vh_common.h"
#include "../bamonly.h"

#define LINE_READING_LENGTH 512
struct LibraryInfo *g_libInfo = NULL;
FILE *divet_file;

int scffab = 0, scffba = 0, scfrab = 0, scfrba = 0, scrfab = 0, scrfba = 0;
long del_cnt_div = 0, ins_cnt_div = 0, inv_cnt_div = 0, tandup_cnt_div = 0, sr_cnt_div = 0, alt_cnt_div = 0;

void vh_pruneAndNormalizeDivets( struct LibraryInfo *lib, double preProsPrune, int overMapLimit)
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

			switch( toBeDeleted->svType)
			{
			case TANDEMDUP:
				tandup_cnt_div--;
				break;
			case DELETION:
				del_cnt_div--;
				break;
			case INSERTION:
				ins_cnt_div--;
				break;
			case INVERSION:
				inv_cnt_div--;
				break;
			}

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
	fprintf( logFile, "After Pruning = %d; %li Deletion, %li Inversion, %li Insertion, %li Tandem Duplication divets\n", lib->size, del_cnt_div, inv_cnt_div, ins_cnt_div, tandup_cnt_div);
}

void vh_printDivet (DivetRow * divetRow)
{
	char msg[LINE_READING_LENGTH];
	sprintf (msg, "%s %s %d %d %c %d %d %c %f %f %g %lu\n",
			divetRow->readName->readName,
			divetRow->chromosome_name,
			divetRow->locMapLeftEnd,
			divetRow->locMapLeftStart,
			divetRow->orientationLeft,
			divetRow->locMapRightStart,
			divetRow->locMapRightEnd,
			divetRow->orientationRight,
			divetRow->editDistance,
			divetRow->avgQual,
			divetRow->phredScore,
			divetRow->ten_x_barcode);

	vh_logOutput (msg);
}

struct DivetRow *createDivetRow (
		struct ReadName *hash[],
		char *readName,
		char *chromosome_name,
		char *locMapLeftStart,
		char *locMapLeftEnd,
		char *orientationLeft,
		char *locMapRightStart,
		char *locMapRightEnd,
		char *orientationRight,
		char *svType,
		char *editDistance,
		char *avgQual,	//skip
		char *phredScore,
		unsigned long ten_x_barcode,
		struct LibraryInfo *libInfo,
		int id)
{
	struct DivetRow *newRow = (DivetRow *) getMem (sizeof (DivetRow));
	if (newRow == NULL)
	{
		strcpy (g_error_message, "Memory Problem in loading divet file!\n");
		vh_quitProgram (EXIT_CODE_MEMORY_ERROR);
	}

	newRow->readName = NULL;
	newRow->next = NULL;

	newRow->chromosome_name = NULL;

	set_str ( &(newRow->chromosome_name), chromosome_name);


	int mls = atoi (locMapLeftStart);
	int mrs = atoi (locMapRightStart);

	if (mls > mrs)		//swap
	{
		newRow->locMapLeftEnd = atoi (locMapRightEnd);
		newRow->locMapRightStart = atoi (locMapLeftStart);
		newRow->locMapLeftStart = atoi (locMapRightStart);
		newRow->locMapRightEnd = atoi (locMapLeftEnd);
		newRow->orientationLeft = orientationRight[0];
		newRow->orientationRight = orientationLeft[0];
	}
	else
	{
		newRow->locMapLeftEnd = atoi (locMapLeftEnd);
		newRow->locMapRightStart = atoi (locMapRightStart);
		newRow->locMapLeftStart = atoi (locMapLeftStart);
		newRow->locMapRightEnd = atoi (locMapRightEnd);
		newRow->orientationLeft = orientationLeft[0];
		newRow->orientationRight = orientationRight[0];
	}

	newRow->ten_x_barcode = ten_x_barcode;
	newRow->avgQual = atof (avgQual);
	newRow->editDistance = atof (editDistance);
	newRow->phredScore = atof (phredScore);
	newRow->svType = svType[0];

	newRow->libInfo = libInfo;
	newRow->in_cluster_count = 0;

	struct ReadName *r = vh_addReadName (hash, readName, newRow->editDistance, newRow->phredScore);
	newRow->readName = r;
	newRow->divetRowId = id;

	return newRow;

}

struct DivetRow *vh_loadDivetRowFromString (struct ReadName *hash[], char *line, struct LibraryInfo *libInfo, int id)
{
	char *readName = strtok (line, DIVET_ROW_DELIMITERS);
	unsigned long ten_x_barcode = -1;


	if (ten_x_flag || output_hs_flag){
		sscanf(readName + strlen(readName) - 20, "%20lu%0", &ten_x_barcode); // 20 is the number of digits in the largest unsigned long value. All read names have ten-x barcode info as a 20 digit, 0-padded suffix
		if (ten_x_barcode != (unsigned long)-1){
			ten_x_barcode = ten_x_barcode | ((unsigned long)libInfo->libId << (sizeof(unsigned long)-1)*8); //use the last byte of the ten x barcode to make sure no two libraries share a barcode
		}
	}

	char *chromosome_name = strtok (NULL, DIVET_ROW_DELIMITERS);

	char *locMapLeftStart = strtok (NULL, DIVET_ROW_DELIMITERS);	//skip
	char *locMapLeftEnd = strtok (NULL, DIVET_ROW_DELIMITERS);
	char *orientationLeft = strtok (NULL, DIVET_ROW_DELIMITERS);

	char *chromosome_name2 = strtok (NULL, DIVET_ROW_DELIMITERS);	// This is the same as chromosome_name for IDVE (indicated by =) and for T it is another chr name (We are not deadling with it here)

	char *locMapRightStart = strtok (NULL, DIVET_ROW_DELIMITERS);
	char *locMapRightEnd = strtok (NULL, DIVET_ROW_DELIMITERS);	//skip
	char *orientationRight = strtok (NULL, DIVET_ROW_DELIMITERS);

	char *svType = strtok (NULL, DIVET_ROW_DELIMITERS);	// sv_type

	char *editDistance = strtok (NULL, DIVET_ROW_DELIMITERS);
	char *avgQual = strtok (NULL, DIVET_ROW_DELIMITERS);
	char *phredScore = strtok (NULL, DIVET_ROW_DELIMITERS);


	struct DivetRow *newRow = createDivetRow (
			hash,
			readName,
			chromosome_name,
			locMapLeftStart, 
			locMapLeftEnd,
			orientationLeft,
			locMapRightStart, 
			locMapRightEnd,
			orientationRight,
			svType,
			editDistance, 
			avgQual,
			phredScore,
			ten_x_barcode,
			libInfo, 
			id);


	return newRow;
}

DivetRow *createDivetRowNew( ReadName *hash[], unsigned long ten_x_barcode, char *readName, char *chromosome_name, int locMapLeftStart, int locMapLeftEnd, char orientationLeft,
		int locMapRightStart, int locMapRightEnd, char orientationRight, char svType, char splitOrientation, int editDistance, int mQual1, int mQual2, float
		avgQual, float phredScore,  LibraryInfo *libInfo, int id)
{
	DivetRow *newRow = ( DivetRow*) getMem( sizeof( DivetRow));

	newRow->readName = NULL;
	newRow->next = NULL;

	newRow->chromosome_name = NULL;
	set_str( &(newRow->chromosome_name), chromosome_name);

	newRow->ten_x_barcode = ten_x_barcode;

	newRow->locMapLeftEnd = locMapLeftEnd;
	newRow->locMapLeftStart = locMapLeftStart;
	newRow->locMapRightStart = locMapRightStart;
	newRow->locMapRightEnd = locMapRightEnd;

	newRow->orientationLeft = orientationLeft;
	newRow->orientationRight = orientationRight;

	newRow->editDistance = ( double) editDistance;
	newRow->avgQual = 1;
	if(mQual1 > 40)
		mQual1 = 40;
	newRow->phredScore = mQual1;

	newRow->mQual1 = mQual1;
	newRow->mQual2 = mQual2;

	newRow->svType = svType;
	newRow->splitOrientation = splitOrientation;

	newRow->libInfo = libInfo;
	newRow->in_cluster_count = 0;

	ReadName *r = vh_addReadName( hash, readName, newRow->editDistance, newRow->phredScore);
	newRow->readName = r;

	newRow->divetRowId = id;

	return newRow;
}


void vh_freeDivets (struct LibraryInfo *libInfo)
{
	//TODO: To be implemented - FreeDivets
	struct DivetRow *temp = libInfo->head;
	while (temp != NULL)
	{
		//Don't free the readName, as it is shared between more than one divet row
		free (temp->chromosome_name);
		libInfo->head = temp->next;
		free (temp);
		temp = libInfo->head;
	}
	libInfo->head = NULL;
	libInfo->tail = NULL;
	libInfo->size = 0;
}



DivetRow *vh_loadDivetRowFromBam(discordantMapping *discordantReadPtr, LibraryInfo *libInfo, int counterDivetRow)
{
	if( discordantReadPtr->svType == TANDEMDUP || discordantReadPtr->svType == INSERTION
			|| discordantReadPtr->svType == DELETION || discordantReadPtr->svType == INVERSION)
	{
		DivetRow *newRow = createDivetRowNew (libInfo->hash, discordantReadPtr->ten_x_barcode, discordantReadPtr->readName, discordantReadPtr->chromosome_name,
				discordantReadPtr->pos1, discordantReadPtr->pos1_End, discordantReadPtr->orient1,
				discordantReadPtr->pos2, discordantReadPtr->pos2_End, discordantReadPtr->orient2,
				discordantReadPtr->svType, NONE, discordantReadPtr->editDistance, discordantReadPtr->mQual1, discordantReadPtr->mQual2, 0, 0, libInfo, counterDivetRow);
		return newRow;
	}
	return NULL;
}

DivetRow *vh_loadDivetRowFromBamAlternative( alternativeMapping *discordantReadPtr, LibraryInfo *libInfo, int conc_min,
		int conc_max, int counterDivetRow)
{
	int pos1_1, pos1_2, pos2_1, pos2_2, svType;
	int32_t isize;
	char orientLeft, orientRight, SV;
	int tmp;
	uint16_t tmp_flag;


	if( discordantReadPtr->pos1 > discordantReadPtr->pos2)
	{
		/*
		tmp = discordantReadPtr->pos1;
		discordantReadPtr->pos1 = discordantReadPtr->pos2;
		discordantReadPtr->pos2 = tmp;

		tmp = discordantReadPtr->pos1_End;
		discordantReadPtr->pos1_End = discordantReadPtr->pos2_End;
		discordantReadPtr->pos2_End = tmp;

		tmp = discordantReadPtr->mQual1;
		discordantReadPtr->mQual1 = discordantReadPtr->mQual2;
		discordantReadPtr->mQual2 = tmp;

		tmp_flag = discordantReadPtr->flag;

		if( ( tmp_flag & BAM_FREVERSE) != 0 )
			orientLeft = REVERSE;
		else
			orientLeft = FORWARD;


		if( ( tmp_flag & BAM_FMREVERSE) != 0 )
			orientRight = REVERSE;
		else
			orientRight = FORWARD;

		if( orientLeft == FORWARD)
			discordantReadPtr->flag = tmp_flag & 0xFFDF;
		else if( orientLeft == REVERSE)
			discordantReadPtr->flag = tmp_flag | BAM_FMREVERSE;

		tmp_flag = discordantReadPtr->flag;
		if( orientRight == FORWARD)
			discordantReadPtr->flag = tmp_flag & 0xFFEF;
		else if( orientRight == REVERSE)
			discordantReadPtr->flag = tmp_flag | BAM_FREVERSE;
		 */
		return NULL;
	}

	if( !is_proper( discordantReadPtr->flag))
		return NULL;

	isize = calculateInsertSize( discordantReadPtr->pos1, discordantReadPtr->pos2, discordantReadPtr->flag, libInfo->readLen);
	svType = is_concordant_bamonly( discordantReadPtr->pos1, discordantReadPtr->pos2, discordantReadPtr->flag, isize,
			conc_min, conc_max);

	switch( svType)
	{
	case RPTDUP:
		SV = TANDEMDUP;
		break;
	case RPDEL:
		SV = DELETION;
		break;
	case RPINS:
		SV = INSERTION;
		break;
	case RPINV:
		SV = INVERSION;
		break;
	default:
		return NULL;
	}

	if( ( discordantReadPtr->flag & BAM_FREVERSE) != 0)
		orientLeft = REVERSE;
	else
		orientLeft = FORWARD;

	if( ( discordantReadPtr->flag & BAM_FMREVERSE) != 0)
		orientRight = REVERSE;
	else
		orientRight = FORWARD;


	fprintf( divet_file, "%s\t%s\t%d\t%d\t%c\t=\t%d\t%d\t%c\t%c\t%d\t%d\t%d\n", discordantReadPtr->readName, discordantReadPtr->chromosome_name,
			discordantReadPtr->pos1, discordantReadPtr->pos1_End, orientLeft, discordantReadPtr->pos2,
			discordantReadPtr->pos2_End, orientRight, SV, discordantReadPtr->editDistance,
			discordantReadPtr->mQual1, discordantReadPtr->mQual2);

	DivetRow *newRow = createDivetRowNew (libInfo->hash, discordantReadPtr->ten_x_barcode,
			discordantReadPtr->readName, discordantReadPtr->chromosome_name, discordantReadPtr->pos1,
			discordantReadPtr->pos1_End, orientLeft, discordantReadPtr->pos2,
			discordantReadPtr->pos2_End, orientRight, SV, NONE, discordantReadPtr->editDistance,
			discordantReadPtr->mQual1, discordantReadPtr->mQual2, 0, 0, libInfo, counterDivetRow);

	return newRow;
}

DivetRow *vh_loadDivetRowFromBamSoftClip( softClip *ptrSoftClip, posMapSoftClip *ptrPosMapSoftClip, LibraryInfo *libInfo, int readLen, int counterDivetRow)
{
	int pos1_1, pos1_2, pos2_1, pos2_2;

	/* Length of read A(left) and read B(right) */
	int lengthRead, lengthSoft;

	/* If soft clip is at the beginning */
	if( ptrSoftClip->op[0] == BAM_CSOFT_CLIP)
	{
		lengthSoft = ptrSoftClip->opl[0];
		lengthRead = readLen - lengthSoft;

		if( ptrSoftClip->pos < ptrPosMapSoftClip->posMap)
		{
			pos1_1 = ptrSoftClip->pos;
			pos1_2 = ptrSoftClip->pos + lengthRead;
			pos2_1 = ptrPosMapSoftClip->posMap;
			pos2_2 = ptrPosMapSoftClip->posMap + lengthSoft;
		}
		else if( ptrPosMapSoftClip->posMap < ptrSoftClip->pos)
		{
			pos1_1 = ptrPosMapSoftClip->posMap;
			pos1_2 = ptrPosMapSoftClip->posMap + lengthSoft;
			pos2_1 = ptrSoftClip->pos;
			pos2_2 = ptrSoftClip->pos + lengthRead;
		}
	}
	/* If soft clip is at the end */
	else if( ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
	{
		lengthSoft = ptrSoftClip->opl[ptrSoftClip->opCount - 1];
		lengthRead = readLen - lengthSoft;

		if( ptrSoftClip->pos < ptrPosMapSoftClip->posMap)
		{
			pos1_1 = ptrSoftClip->pos;
			pos1_2 = ptrSoftClip->pos + lengthRead;
			pos2_1 = ptrPosMapSoftClip->posMap;
			pos2_2 = ptrPosMapSoftClip->posMap + lengthSoft;
		}
		else if( ptrPosMapSoftClip->posMap < ptrSoftClip->pos)
		{
			pos1_1 = ptrPosMapSoftClip->posMap;
			pos1_2 = ptrPosMapSoftClip->posMap + lengthSoft;
			pos2_1 = ptrSoftClip->pos;
			pos2_2 = ptrSoftClip->pos + lengthRead;
		}
	}
	if( pos1_2 >= pos2_1)
		return NULL;

	if( ptrSoftClip->orient == FORWARD && ptrPosMapSoftClip->orient == FORWARD)
	{
		if( ( ptrSoftClip->pos < ptrPosMapSoftClip->posMap && ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
				|| ( ptrPosMapSoftClip->posMap < ptrSoftClip->pos && ptrSoftClip->op[0] == BAM_CSOFT_CLIP))
		{

			DivetRow * newRow = createDivetRowNew (libInfo->hash, 0, ptrSoftClip->readName, ptrSoftClip->chromosome_name,
					pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, FORWARD, NONE, FFAB, 0, ptrSoftClip->qual,
					ptrPosMapSoftClip->mapq, 0, 0, libInfo, counterDivetRow);

			scffab++;
			return newRow;
		}
		else if( ( ptrPosMapSoftClip->posMap < ptrSoftClip->pos && ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
				|| ( ptrSoftClip->pos < ptrPosMapSoftClip->posMap && ptrSoftClip->op[0] == BAM_CSOFT_CLIP))
		{

			DivetRow * newRow = createDivetRowNew (libInfo->hash, 0, ptrSoftClip->readName, ptrSoftClip->chromosome_name,
					pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, FORWARD, NONE, FFBA, 0, ptrSoftClip->qual,
					ptrPosMapSoftClip->mapq, 0, 0, libInfo, counterDivetRow);
			//fprintf(stderr,"%d - %d - %d - %d\n", pos1_1, pos1_2, pos2_1, pos2_2);
			scffba++;
			return newRow;
		}
	}
	else if( ptrSoftClip->orient == FORWARD && ptrPosMapSoftClip->orient == REVERSE)
	{
		if( ( ptrSoftClip->pos < ptrPosMapSoftClip->posMap && ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP))
		{

			DivetRow * newRow = createDivetRowNew (libInfo->hash, 0, ptrSoftClip->readName, ptrSoftClip->chromosome_name,
					pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, REVERSE, NONE, FRAB, 0, ptrSoftClip->qual,
					ptrPosMapSoftClip->mapq, 0, 0, libInfo, counterDivetRow);
			scfrab++;
			return newRow;
		}
		else if( ptrPosMapSoftClip->posMap < ptrSoftClip->pos && ptrSoftClip->op[0] == BAM_CSOFT_CLIP)
		{
			DivetRow * newRow = createDivetRowNew (libInfo->hash, 0, ptrSoftClip->readName, ptrSoftClip->chromosome_name,
					pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, REVERSE, NONE, RFAB, 0, ptrSoftClip->qual,
					ptrPosMapSoftClip->mapq, 0, 0, libInfo, counterDivetRow);
			scrfab++;
			return newRow;
		}
		else if( ( ptrPosMapSoftClip->posMap < ptrSoftClip->pos && ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP))
		{
			DivetRow * newRow = createDivetRowNew (libInfo->hash, 0, ptrSoftClip->readName, ptrSoftClip->chromosome_name,
					pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, REVERSE, NONE, RFBA, 0, ptrSoftClip->qual,
					ptrPosMapSoftClip->mapq, 0, 0, libInfo, counterDivetRow);
			scrfba++;
			return newRow;
		}
		else if(ptrSoftClip->pos < ptrPosMapSoftClip->posMap && ptrSoftClip->op[0] == BAM_CSOFT_CLIP)
		{
			DivetRow * newRow = createDivetRowNew (libInfo->hash, 0, ptrSoftClip->readName, ptrSoftClip->chromosome_name,
					pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, REVERSE, NONE, FRBA, 0, ptrSoftClip->qual,
					ptrPosMapSoftClip->mapq, 0, 0, libInfo, counterDivetRow);
			//fprintf(stderr,"%d - %d - %d - %d\n", pos1_1, pos1_2, pos2_1, pos2_2);
			scfrba++;
			return newRow;
		}
	}
	return NULL;
}

int read_Divet_bam( discordantMapping **mapping, parameters *params, LibraryInfo * libInfo, int chr_index, int counterDivetRow)
{
	int is_satellite, i;

	struct DivetRow *newRow = NULL;
	discordantMapping *discordantReadPtr;

	for( i = 0; i < NHASH; i++)
	{
		discordantReadPtr = mapping[i];
		while( discordantReadPtr != NULL)
		{
			is_satellite = sonic_is_satellite( params->this_sonic, discordantReadPtr->chromosome_name , discordantReadPtr->pos1, discordantReadPtr->pos1_End)																																																																																															&& sonic_is_satellite( params->this_sonic, discordantReadPtr->chromosome_name, discordantReadPtr->pos2, discordantReadPtr->pos2_End)
									&& sonic_is_satellite( params->this_sonic, discordantReadPtr->chromosome_name , discordantReadPtr->pos2, discordantReadPtr->pos2_End);

			if ( is_satellite == 0 && discordantReadPtr->mQual1 > params->mq_threshold && discordantReadPtr->mQual2 > params->mq_threshold
					&& !sonic_is_gap( params->this_sonic, params->this_sonic->chromosome_names[chr_index], discordantReadPtr->pos1, discordantReadPtr->pos2)
					&& strcmp( discordantReadPtr->chromosome_name, params->this_sonic->chromosome_names[chr_index]) == 0
					&& discordantReadPtr->pos1 > 0 && discordantReadPtr->pos1_End < params->this_sonic->chromosome_lengths[chr_index] &&
					discordantReadPtr->pos2 > 0 && discordantReadPtr->pos2_End < params->this_sonic->chromosome_lengths[chr_index])
			{
				newRow = vh_loadDivetRowFromBam( discordantReadPtr, libInfo, counterDivetRow);
				if( newRow == NULL)
					fprintf( stderr, "ERROR loading divet from bam\n");

				switch( discordantReadPtr->svType)
				{
				case TANDEMDUP:
					tandup_cnt_div++;
					break;
				case DELETION:
					del_cnt_div++;
					break;
				case INSERTION:
					ins_cnt_div++;
					break;
				case INVERSION:
					inv_cnt_div++;
					break;
				}

				if( libInfo->head == NULL || libInfo->tail == NULL)
				{
					libInfo->head = newRow;
					libInfo->tail = newRow;
				}
				else
				{
					libInfo->tail->next = newRow;
					libInfo->tail = newRow;
				}
				libInfo->size++;
				counterDivetRow++;
			}
			discordantReadPtr = discordantReadPtr->next;
		}
	}
	return counterDivetRow;
}

int read_Divet_bam_alternative( alternativeMapping **mapping, parameters *params,
		LibraryInfo *libInfo, int chr_index, int conc_min, int conc_max, int counterDivetRow)
{
	int is_satellite, i, insLen;

	struct DivetRow *newRow = NULL;
	alternativeMapping *discordantReadPtr;

	long alt_cnt_div = 0;

	for( i = 0; i < NHASH; i++)
	{
		discordantReadPtr = mapping[i];
		while( discordantReadPtr != NULL)
		{
			/* Since MT is circular, we need to eliminate the read-pairs at both ends of the chromosome */
			insLen = abs( discordantReadPtr->pos1 - discordantReadPtr->pos2);
			is_satellite = sonic_is_satellite( params->this_sonic, discordantReadPtr->chromosome_name , discordantReadPtr->pos1, discordantReadPtr->pos1_End)																																																																																															&& sonic_is_satellite( params->this_sonic, discordantReadPtr->chromosome_name, discordantReadPtr->pos2, discordantReadPtr->pos2_End) &&
					sonic_is_satellite( params->this_sonic, discordantReadPtr->chromosome_name , discordantReadPtr->pos2, discordantReadPtr->pos2_End)																																																																																															&& sonic_is_satellite( params->this_sonic, discordantReadPtr->chromosome_name, discordantReadPtr->pos2, discordantReadPtr->pos2_End);

			if ( is_satellite == 0 && ( discordantReadPtr->mQual1 > params->mq_threshold) &&
					( !sonic_is_gap( params->this_sonic, params->this_sonic->chromosome_names[chr_index], discordantReadPtr->pos1, discordantReadPtr->pos2)) &&
					( discordantReadPtr->mQual2 > params->mq_threshold) &&
					( strcmp( discordantReadPtr->chromosome_name, params->this_sonic->chromosome_names[chr_index]) == 0) &&
					( discordantReadPtr->pos1 > 0) && ( discordantReadPtr->pos1_End < params->this_sonic->chromosome_lengths[chr_index]) &&
					( discordantReadPtr->pos2 > 0) && ( discordantReadPtr->pos2_End < params->this_sonic->chromosome_lengths[chr_index]) &&
					( insLen < (params->this_sonic->chromosome_lengths[chr_index]) - (2 * conc_max)))
			{
				newRow = vh_loadDivetRowFromBamAlternative( discordantReadPtr, libInfo, conc_min, conc_max, counterDivetRow);
				if( newRow != NULL)
				{
					if( libInfo->head == NULL || libInfo->tail == NULL)
					{
						libInfo->head = newRow;
						libInfo->tail = newRow;
					}
					else
					{
						libInfo->tail->next = newRow;
						libInfo->tail = newRow;
					}
					libInfo->size++;
					counterDivetRow++;
					alt_cnt_div++;
				}
			}
			discordantReadPtr = discordantReadPtr->next;
		}
	}

	fprintf( stderr,"There are %d alternative mapping rows\n", libInfo->size);
	vh_pruneAndNormalizeDivets( libInfo, 0.4, 100);
	sprintf( g_loggerMsgBuffer, "%d rows after pruning.", libInfo->size);
	vh_logInfo( g_loggerMsgBuffer);
	fprintf( logFile,"%d left after pruning\n", libInfo->size);

	return counterDivetRow;
}



int read_Divet_bam_softClip( softClip *ptrSoftClip, parameters *params, LibraryInfo * libInfo, int chr_index, int read_len, int divet_row_count)
{
	int i, is_satellite;
	struct DivetRow *newRow = NULL;
	posMapSoftClip *ptrPosMapSoftClip;

	while( ptrSoftClip != NULL)
	{
		ptrPosMapSoftClip = ptrSoftClip->ptrPosMapSoftClip;
		while( ptrPosMapSoftClip != NULL)
		{
			is_satellite = sonic_is_satellite( params->this_sonic, ptrSoftClip->chromosome_name, ptrSoftClip->pos, ptrSoftClip->pos + 1 ) &&
					sonic_is_satellite( params->this_sonic, ptrSoftClip->chromosome_name, ptrPosMapSoftClip->posMap, ptrPosMapSoftClip->posMap + 1);

			if ( is_satellite == 0 && ptrSoftClip->qual > params->mq_threshold && strcmp(ptrSoftClip->chromosome_name, params->this_sonic->chromosome_names[chr_index]) == 0
					&& ptrPosMapSoftClip->mapq > 0 && ptrSoftClip->pos > 0 && ptrPosMapSoftClip->posMap > 0
					&& ptrSoftClip->pos < params->this_sonic->chromosome_lengths[chr_index] && ptrPosMapSoftClip->posMap < params->this_sonic->chromosome_lengths[chr_index])
			{
				for( i = 0; i < 2; i++)
				{
					newRow = vh_loadDivetRowFromBamSoftClip( ptrSoftClip, ptrPosMapSoftClip, libInfo, read_len, divet_row_count);
					if( newRow != NULL)
					{
						/* For Inversion and Inverted Duplication */
						if( i == 0 && ( newRow->splitOrientation == FRAB || newRow->splitOrientation == FRBA
								|| newRow->splitOrientation == RFAB || newRow->splitOrientation == RFBA))
						{
							newRow->svType = INVERSION;
							newRow->orientationLeft = FORWARD;
							newRow->orientationRight = FORWARD;
							newRow->locMapLeftStart -= SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapLeftEnd -= SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightStart -= SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightEnd -= SOFTCLIP_WRONGMAP_WINDOW;
						}
						else if( i == 1 && ( newRow->splitOrientation == FRAB || newRow->splitOrientation == FRBA
								|| newRow->splitOrientation == RFAB || newRow->splitOrientation == RFBA))
						{
							newRow->svType = INVERSION;
							newRow->orientationLeft = REVERSE;
							newRow->orientationRight = REVERSE;
							newRow->locMapLeftStart += SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapLeftEnd += SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightStart += SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightEnd += SOFTCLIP_WRONGMAP_WINDOW;
						}
						/* For Deletion */
						else if( i == 0 && ( newRow->splitOrientation == FFAB))
						{
							newRow->svType = DELETION;
							newRow->orientationLeft = FORWARD;
							newRow->orientationRight = REVERSE;
							newRow->locMapLeftStart -= SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapLeftEnd -= SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightStart += SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightEnd += SOFTCLIP_WRONGMAP_WINDOW;
						}
						/* For Tandem Duplication */
						else if( i == 0 && ( newRow->splitOrientation == FFBA))
						{
							newRow->svType = TANDEMDUP;
							newRow->orientationLeft = REVERSE;
							newRow->orientationRight = FORWARD;
							newRow->locMapLeftStart -= SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapLeftEnd -= SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightStart += SOFTCLIP_WRONGMAP_WINDOW;
							newRow->locMapRightEnd += SOFTCLIP_WRONGMAP_WINDOW;
						}
						else if( i == 1)
							newRow = NULL;
					}

					if( newRow == NULL)
						;//fprintf( stderr, "ERROR loading divet from bam (soft clip)\n");
					else
					{
						if( libInfo->head == NULL || libInfo->tail == NULL)
						{
							libInfo->head = newRow;
							libInfo->tail = newRow;
						}
						else
						{
							libInfo->tail->next = newRow;
							libInfo->tail = newRow;
						}
						libInfo->size++;
						divet_row_count++;
						sr_cnt_div++;
					}
				}
			}
			ptrPosMapSoftClip = ptrPosMapSoftClip->next;
		}
		ptrSoftClip = ptrSoftClip->next;
	}
	//fprintf(stderr, "\nThere are %li sc in divets;%d FFAB, %d FFBA, %d FRAB, %d FRBA, %d RFAB, %d RFBA\n",
	//	sr_cnt_div, scffab, scffba, scfrab,scfrba, scrfab, scrfba);
	return divet_row_count;
}


int load_Divet_bam( bam_info** in_bams, parameters *params, int chr_index)
{
	int lib_cnt, read_cnt;
	int j, bam_index, chr, divet_row_count = 0;
	struct LibraryInfo *newLibInfo, *cursor, *t;
	char *divetfile_path;
	
	g_maxListBrkPointIntr = MAXLISTBRKPOINTINTR;
	g_libInfo = NULL;

	
	divetfile_path = (char *) getMem(sizeof(char) * (1+strlen("divv.vh")+strlen(params->outdir)));
	sprintf(divetfile_path, "%s%s", params->outdir, "divv.vh");

	divet_file = safe_fopen (divetfile_path, "w");
	free(divetfile_path);
	del_cnt_div = 0, ins_cnt_div = 0, inv_cnt_div = 0, tandup_cnt_div = 0, sr_cnt_div = 0, alt_cnt_div = 0;

	for( bam_index = 0; bam_index < params->num_bams; bam_index++)
	{
		for( lib_cnt = 0; lib_cnt < in_bams[bam_index]->num_libraries; lib_cnt++)
		{
			newLibInfo = ( struct LibraryInfo *) getMem( sizeof( struct LibraryInfo));
			strcpy( newLibInfo->libName, in_bams[bam_index]->libraries[lib_cnt]->libname);
			strcpy( newLibInfo->indName, in_bams[bam_index]->sample_name);
			newLibInfo->minDelta = in_bams[bam_index]->libraries[lib_cnt]->conc_min - 2 * in_bams[bam_index]->libraries[lib_cnt]->read_length;
			newLibInfo->maxDelta = in_bams[bam_index]->libraries[lib_cnt]->conc_max - 2 * in_bams[bam_index]->libraries[lib_cnt]->read_length;
			if( newLibInfo->minDelta < 0)
				newLibInfo->minDelta = 0;
			if( newLibInfo->maxDelta < 0)
				newLibInfo->maxDelta = 0;
			newLibInfo->readLen = in_bams[bam_index]->libraries[lib_cnt]->read_length;

			/* We store the reads in hash[] based on the hash values of read names */
			newLibInfo->hash = ( struct ReadName **) getMem( NHASH * sizeof( struct ReadName *));
			for( j = 0; j < NHASH; j++)
				newLibInfo->hash[j] = NULL;

			newLibInfo->head = NULL;
			newLibInfo->tail = NULL;
			newLibInfo->size = 0;
			newLibInfo->next = NULL;

			divet_row_count = read_Divet_bam( in_bams[bam_index]->libraries[lib_cnt]->mappings_discordant,
					params, newLibInfo, chr_index, divet_row_count);

			if( g_libInfo == NULL)
				g_libInfo = newLibInfo;
			else                  //add to the end of the linked list
			{
				for( t = g_libInfo; t->next != NULL; t = t->next)
					;        //Skip till end of LinkedList
				t->next = newLibInfo;
			}
		}
	}
	fprintf( logFile, "CHR %d - %li Deletion, %li Inversion, %li Insertion, %li Tandem Duplication, ", chr_index, del_cnt_div, inv_cnt_div, ins_cnt_div, tandup_cnt_div);

	if( params->alt_mapping != 0)
	{
		for( bam_index = 0; bam_index < params->num_bams; bam_index++)
		{
			for( lib_cnt = 0; lib_cnt < in_bams[bam_index]->num_libraries; lib_cnt++)
			{
				newLibInfo = ( struct LibraryInfo *) getMem( sizeof( struct LibraryInfo));
				strcpy( newLibInfo->libName, "Alternative\0");
				strcpy( newLibInfo->indName, in_bams[bam_index]->sample_name);
				newLibInfo->minDelta = in_bams[bam_index]->libraries[lib_cnt]->conc_min - 2 * in_bams[bam_index]->libraries[lib_cnt]->read_length;
				newLibInfo->maxDelta = in_bams[bam_index]->libraries[lib_cnt]->conc_max - 2 * in_bams[bam_index]->libraries[lib_cnt]->read_length;
				if( newLibInfo->minDelta < 0)
					newLibInfo->minDelta = 0;
				if( newLibInfo->maxDelta < 0)
					newLibInfo->maxDelta = 0;
				newLibInfo->readLen = in_bams[bam_index]->libraries[lib_cnt]->read_length;

				/* We store the reads in hash[] based on the hash values of read names */
				newLibInfo->hash = ( struct ReadName **) getMem( NHASH * sizeof( struct ReadName *));
				for( j = 0; j < NHASH; j++)
					newLibInfo->hash[j] = NULL;

				newLibInfo->head = NULL;
				newLibInfo->tail = NULL;
				newLibInfo->size = 0;
				newLibInfo->next = NULL;

				divet_row_count = read_Divet_bam_alternative( in_bams[bam_index]->libraries[lib_cnt]->mappings_alternative,
						params, newLibInfo, chr_index, in_bams[bam_index]->libraries[lib_cnt]->conc_min,
						in_bams[bam_index]->libraries[lib_cnt]->conc_max, divet_row_count);

				if( g_libInfo == NULL)
					g_libInfo = newLibInfo;
				else                  //add to the end of the linked list
				{
					for( t = g_libInfo; t->next != NULL; t = t->next)
						;        //Skip till end of LinkedList
					t->next = newLibInfo;
				}
			}
		}
	}
	fprintf( logFile, "%li Alternative Mapping, ", alt_cnt_div);
	if( !params->no_soft_clip)
	{
		for( bam_index = 0; bam_index < params->num_bams; bam_index++)
		{
			for( lib_cnt = 0; lib_cnt < in_bams[bam_index]->num_libraries; lib_cnt++)
			{
				newLibInfo = ( struct LibraryInfo *) getMem( sizeof( struct LibraryInfo));
				strcpy( newLibInfo->libName, "SplitRead\0");
				strcpy( newLibInfo->indName, in_bams[bam_index]->sample_name);
				newLibInfo->minDelta = 0;
				newLibInfo->maxDelta = 5 * SOFTCLIP_WRONGMAP_WINDOW;
				newLibInfo->readLen = in_bams[bam_index]->libraries[lib_cnt]->read_length;

				newLibInfo->hash = ( struct ReadName **) getMem( NHASH * sizeof( struct ReadName *));
				for( j = 0; j < NHASH; j++)
					newLibInfo->hash[j] = NULL;

				newLibInfo->head = NULL;
				newLibInfo->tail = NULL;
				newLibInfo->size = 0;
				newLibInfo->next = NULL;

				divet_row_count = read_Divet_bam_softClip( in_bams[bam_index]->libraries[lib_cnt]->listSoftClip, params,
						newLibInfo, chr_index, newLibInfo->readLen, divet_row_count);

				if( g_libInfo == NULL)
					g_libInfo = newLibInfo;
				else
				{
					for( t = g_libInfo; t->next != NULL; t = t->next)
						;        //Skip till end of LinkedList
					t->next = newLibInfo;
				}
			}
		}
	}
	fprintf( logFile, "%li Split Read Divets\n", sr_cnt_div);
	fclose( divet_file);
	return divet_row_count;
}

/**
 * Divet file is read and loaded into the linked lists
 * After running this method, the values of g_headDivet, g_tailDivet and g_divetRowLinkedListSize are set
 */
DivetRow *vh_loadDivetFile (LibraryInfo * libInfo, sonic *this_sonic)
{
	int i;
	char line[LINE_READING_LENGTH]; //shouldn't this be longer?
	char div_location[MAX_SEQ];
	int counter = 0;
	char *token;
	char *return_value = 0;

	del_cnt_div = 0, ins_cnt_div = 0, inv_cnt_div = 0, tandup_cnt_div = 0;
	FILE *divetFile = fopen (libInfo->libFileAdrs, "r");
	if (!divetFile)
	{
		sprintf( div_location, "divet/%s", libInfo->libFileAdrs);
		divetFile = safe_fopen( div_location, "r");
		if (divetFile == NULL)
		{
			sprintf (g_error_message, "Divet file '%s' could not be opened!", libInfo->libFileAdrs);
			vh_quitProgram (EXIT_CODE_DIVET_ERROR);
		}
	}

	//Initializing the linked list of divet rows
	libInfo->head = NULL;
	libInfo->tail = NULL;
	libInfo->size = 0;

	int counterDivetRow = 0;
	while (!feof (divetFile))
	{
		return_value = fgets (line, LINE_READING_LENGTH, divetFile);

		if (line == NULL)
			continue;
		// If the read line is empty
		int isEmpty = 1;
		int len = strlen (line);
		for (i = 0; i < len; i++)
		{
			if (!isspace (line[i]))
			{
				isEmpty = 0;
				break;
			}
		}
		if (isEmpty)
			continue;

		DivetRow *newRow = vh_loadDivetRowFromString (libInfo->hash, line, libInfo, counterDivetRow);

		switch( newRow->svType)
		{
		case TANDEMDUP:
			tandup_cnt_div++;
			break;
		case DELETION:
			del_cnt_div++;
			break;
		case INSERTION:
			ins_cnt_div++;
			break;
		case INVERSION:
			inv_cnt_div++;
			break;
		}

		if (!sonic_is_satellite( this_sonic, newRow->chromosome_name, newRow->locMapLeftStart, newRow->locMapLeftEnd) &&
				!sonic_is_satellite( this_sonic, newRow->chromosome_name, newRow->locMapRightStart, newRow->locMapRightEnd) &&
				!sonic_is_gap( this_sonic, newRow->chromosome_name, newRow->locMapLeftStart, newRow->locMapRightEnd))
		{
			if (libInfo->head == NULL || libInfo->tail == NULL)
			{
				libInfo->head = newRow;
				libInfo->tail = newRow;
			}
			else
			{
				libInfo->tail->next = newRow;
				libInfo->tail = newRow;
			}
			libInfo->size++;
			counterDivetRow++;
		}
		else
			free (newRow);

		line[0] = '\0';
	}

	vh_logInfo ("Finished reading the file");
	fclose (divetFile);
	fprintf( logFile, "Total = %d; %li Deletion, %li Inversion, %li Insertion, %li Tandem Duplication divets\n", libInfo->size, del_cnt_div, inv_cnt_div, ins_cnt_div, tandup_cnt_div);

	return libInfo->head;
}
