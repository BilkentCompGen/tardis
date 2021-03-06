#include "vh_createMaxClusterMEI.h"
#include "vh_common.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_setcover.h"
#include "vh_heap.h"
#include "../bamonly.h"

mei_Reads **mReads = NULL;
HeapMEI *H_F; /* The forward mapping heap from left of breakpoint */
HeapMEI *H_R; /* the reverse mapping heap from right of breakpoint */
HeapMEI *H_S; /* The soft clipping read */

int addToGenomeIndex_MEI (bam_info** in_bams, parameters *params, char *chromosome_name, int chroSize)
{
	mei_Reads *tmpMEI_Reads;
	discordantMappingMEI *discordantReadPtr;
	softClip *softClipPtr;
	int numSample, count, mei_count = 0, tmp;
	float is_satellite;

	sonic_repeat *mei;
	int libId = 0;

	for( numSample = 0; numSample < params->num_bams; numSample++)
	{
		for( count = 0; count < in_bams[numSample]->num_libraries; count++)
		{
			discordantReadPtr = in_bams[numSample]->libraries[count]->listMEI_Mapping;
			softClipPtr = in_bams[numSample]->libraries[count]->listSoftClip;
			while( discordantReadPtr != NULL)
			{
				is_satellite = sonic_is_satellite (params->this_sonic, discordantReadPtr->chromosome_name, discordantReadPtr->pos, discordantReadPtr->pos_End);
				if( is_satellite == 0 && strcmp( discordantReadPtr->chromosome_name, chromosome_name) == 0
						&& ( discordantReadPtr->pos > 0) && ( discordantReadPtr->pos_End < chroSize))
				{
					if( discordantReadPtr->MEI_Type > -1)
					{
						tmpMEI_Reads = (mei_Reads*) getMem( sizeof( mei_Reads));

						tmpMEI_Reads->readName = NULL;
						set_str( &(tmpMEI_Reads->readName), discordantReadPtr->readName);

						tmpMEI_Reads->MEI_Subclass = NULL;
						set_str( &(tmpMEI_Reads->MEI_Subclass), discordantReadPtr->MEI_subclass);

						tmpMEI_Reads->MEI_Class = NULL;
						set_str( &(tmpMEI_Reads->MEI_Class), discordantReadPtr->MEI_class);

						tmpMEI_Reads->pos = discordantReadPtr->pos;
						tmpMEI_Reads->pos_End = discordantReadPtr->pos_End;
						tmpMEI_Reads->orient = discordantReadPtr->orient;
						tmpMEI_Reads->mQual = discordantReadPtr->qual;
						tmpMEI_Reads->MEI_Type = discordantReadPtr->MEI_Type;
						tmpMEI_Reads->readTypeSupport = 0; // The support is type paired-end read mapping
						tmpMEI_Reads->libId = libId;
						tmpMEI_Reads->ten_x_barcode = discordantReadPtr->ten_x_barcode;

						tmpMEI_Reads->libName = NULL;
						set_str( &(tmpMEI_Reads->libName), in_bams[numSample]->libraries[count]->libname);

						tmpMEI_Reads->indName = NULL;
						set_str( &tmpMEI_Reads->indName, in_bams[numSample]->sample_name);

						if (discordantReadPtr->orient == FORWARD)
						{
							tmpMEI_Reads->next = mReads[discordantReadPtr->pos_End];
							mReads[discordantReadPtr->pos_End] = tmpMEI_Reads;
						}
						if (discordantReadPtr->orient == REVERSE)
						{
							tmp = max( 0, discordantReadPtr->pos - in_bams[numSample]->libraries[count]->conc_max);
							tmpMEI_Reads->next = mReads[tmp];
							mReads[tmp] = tmpMEI_Reads;
						}
						mei_count++;
					}
				}
				discordantReadPtr = discordantReadPtr->next;
			}/*
			while( softClipPtr != NULL)
			{
				is_satellite = sonic_is_satellite (params->this_sonic, softClipPtr->chromosome_name, softClipPtr->pos, softClipPtr->pos+1);
				if( is_satellite == 0 && strcmp( softClipPtr->chromosome_name, chromosome_name) == 0 && ( softClipPtr->pos > 0))
				{
					posMapSoftClip *ptrPosMapSoftClip = softClipPtr->ptrPosMapSoftClip;
					if( ptrPosMapSoftClip != NULL)
					{
						mei = sonic_is_mobile_element( params->this_sonic, chromosome_name, ptrPosMapSoftClip->posMap, ptrPosMapSoftClip->posMap+1, params->mei);
						if( mei != NULL && ( ptrPosMapSoftClip->posMap < chroSize))
						{
							tmpMEI_Reads = (mei_Reads*) getMem(sizeof(mei_Reads));

							tmpMEI_Reads->readName = NULL;
							set_str( &(tmpMEI_Reads->readName), softClipPtr->readName);

							tmpMEI_Reads->MEI_Subclass = NULL;
							set_str( &(tmpMEI_Reads->MEI_Subclass), mei->repeat_type);

							tmpMEI_Reads->MEI_Class = NULL;
							set_str( &(tmpMEI_Reads->MEI_Class), mei->repeat_class);

							tmpMEI_Reads->pos = softClipPtr->pos;
							tmpMEI_Reads->orient = softClipPtr->orient;
							tmpMEI_Reads->mQual = softClipPtr->qual;
							tmpMEI_Reads->libId = libId;

							tmpMEI_Reads->libName = NULL;
							set_str( &(tmpMEI_Reads->libName), in_bams[numSample]->libraries[count]->libname);

							tmpMEI_Reads->indName = NULL;
							set_str( &(tmpMEI_Reads->indName), in_bams[numSample]->sample_name);

							if( softClipPtr->op[0] == 4) //soft clipped at the beginning of the read. Horrible. Replace with #define
							{
								tmpMEI_Reads->readTypeSupport = 1;
								tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos-lenSplitReadBrakWindow];
								mReads[tmpMEI_Reads->pos - lenSplitReadBrakWindow] = tmpMEI_Reads;
							}
							else if( softClipPtr->op[softClipPtr->opCount - 1] == 4) // soft cliped at the end of the read. Horrible. Replace with #define
							{
								tmpMEI_Reads->readTypeSupport = 2;
								tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos + softClipPtr->opl[0] - lenSplitReadBrakWindow];
								mReads[tmpMEI_Reads->pos + softClipPtr->opl[0] - lenSplitReadBrakWindow] = tmpMEI_Reads;
							}
							mei_count++;
						}
					}
				}
				softClipPtr = softClipPtr->next;
			}*/
			libId++;
		}
	}
	return mei_count;
}


int addToGenomeIndex_MEI_sensitive( parameters *params, char *chromosome_name, int chroSize)
{
	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	mei_Reads *tmpMEI_Reads;
	int libId = 0, mei_count = 0;

	libInfo = g_libInfo;
	while (libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while (divetReadMappingPtr != NULL)
		{
			if (strcmp (divetReadMappingPtr->chromosome_name, chromosome_name) == 0
					&& ((divetReadMappingPtr->svType == MEIRIGHTPAIR || divetReadMappingPtr->svType == MEILEFTPAIR))
					&& ( divetReadMappingPtr->locMapLeftStart > 0) && ( divetReadMappingPtr->locMapRightEnd < chroSize))
			{
				tmpMEI_Reads = (mei_Reads *)getMem(sizeof(mei_Reads));
				tmpMEI_Reads->readName = NULL;
				set_str( &(tmpMEI_Reads->readName), divetReadMappingPtr->readName->readName);

				tmpMEI_Reads->mQual = divetReadMappingPtr->phredScore;
				tmpMEI_Reads->readTypeSupport = 0; // The support is type paired-end read mapping

				tmpMEI_Reads->libName = NULL;
				set_str( &(tmpMEI_Reads->libName), libInfo->libName);

				tmpMEI_Reads->indName = NULL;
				set_str( &(tmpMEI_Reads->indName), libInfo->indName);

				tmpMEI_Reads->MEI_Subclass = NULL; 
				set_str( &(tmpMEI_Reads->MEI_Subclass), divetReadMappingPtr->mei_subclass);

				tmpMEI_Reads->MEI_Class = NULL;
				set_str( &(tmpMEI_Reads->MEI_Class), divetReadMappingPtr->meiType);

				tmpMEI_Reads->libId = libId;
				tmpMEI_Reads->ten_x_barcode = divetReadMappingPtr->ten_x_barcode;

				tmpMEI_Reads->MEI_Type = divetReadMappingPtr->mei_code;

				if (divetReadMappingPtr->svType == MEIRIGHTPAIR)
				{
					if(divetReadMappingPtr->orientationLeft == FORWARD)
					{
						tmpMEI_Reads->pos = divetReadMappingPtr->locMapLeftStart;
						tmpMEI_Reads->pos_End = divetReadMappingPtr->locMapLeftEnd;
						tmpMEI_Reads->orient = FORWARD;
						tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos_End];
						mReads[tmpMEI_Reads->pos_End] = tmpMEI_Reads;
					}
					else if (divetReadMappingPtr->orientationLeft == REVERSE)
					{
						tmpMEI_Reads->pos = divetReadMappingPtr->locMapLeftStart;
						tmpMEI_Reads->pos_End = divetReadMappingPtr->locMapLeftEnd;
						tmpMEI_Reads->orient = REVERSE;
						tmpMEI_Reads->next = mReads[max(0,tmpMEI_Reads->pos-libInfo->maxDelta)];
						mReads[max( 0, tmpMEI_Reads->pos-libInfo->maxDelta)] = tmpMEI_Reads;
					}
				}
				else if (divetReadMappingPtr->svType == MEILEFTPAIR)
				{
					if(divetReadMappingPtr->orientationRight == FORWARD)
					{
						tmpMEI_Reads->pos = divetReadMappingPtr->locMapRightStart;
						tmpMEI_Reads->pos_End = divetReadMappingPtr->locMapRightEnd;
						tmpMEI_Reads->orient = FORWARD;
						tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos_End];
						mReads[tmpMEI_Reads->pos_End] = tmpMEI_Reads;
					}
					else if (divetReadMappingPtr->orientationRight == REVERSE)
					{
						tmpMEI_Reads->pos = divetReadMappingPtr->locMapRightStart;
						tmpMEI_Reads->pos_End = divetReadMappingPtr->locMapRightEnd;
						tmpMEI_Reads->orient = REVERSE;
						tmpMEI_Reads->next = mReads[max(0,tmpMEI_Reads->pos-libInfo->maxDelta)];
						mReads[max( 0, tmpMEI_Reads->pos-libInfo->maxDelta)] = tmpMEI_Reads;
					}
				}
				mei_count++;
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
		libId++;
	}
	return mei_count;
}


void add_R_Heap( int pos)
{
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapElMEI *newEl;
	while(mei_ReadsPtr!=NULL)
	{
		if (mei_ReadsPtr->orient == REVERSE && mei_ReadsPtr->readTypeSupport == 0)
		{
		        newEl = (HeapElMEI *)getMem(sizeof(HeapElMEI));
			newEl->mei_ReadsPtr = mei_ReadsPtr;
			newEl->priorityValue = mei_ReadsPtr->pos;
			push_heap_mei(H_R, newEl);
			free(newEl);
		}
		mei_ReadsPtr = mei_ReadsPtr->next;
	}
}

void add_S_Heap( int pos)
{
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapElMEI *newEl;
	while(mei_ReadsPtr!=NULL)
	{
		if ( mei_ReadsPtr->readTypeSupport == 1 || mei_ReadsPtr->readTypeSupport == 2)
		{
		        newEl = ( HeapElMEI *) getMem( sizeof( HeapElMEI));
			newEl->mei_ReadsPtr = mei_ReadsPtr;
			newEl->priorityValue = pos + ( 2 * lenSplitReadBrakWindow);
			push_heap_mei( H_S, newEl);
			free( newEl);
		}
		mei_ReadsPtr = mei_ReadsPtr->next;
	}
}

void add_F_Heap( int pos)
{
	LibraryInfo *libInfo;
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapElMEI *newEl;
	int max_delta_val;

	while( mei_ReadsPtr != NULL)
	{
		if( mei_ReadsPtr->orient == FORWARD && mei_ReadsPtr->readTypeSupport == 0)
		{
		        newEl = ( HeapElMEI *) getMem( sizeof( HeapElMEI));		 
			newEl->mei_ReadsPtr = mei_ReadsPtr;
			libInfo = g_libInfo;
			max_delta_val = 0;
			while( libInfo != NULL)
			{
				if( strcmp( mei_ReadsPtr->libName, libInfo->libName) == 0 && strcmp( mei_ReadsPtr->indName, libInfo->indName) == 0)
				{
					max_delta_val = libInfo->maxDelta;
					newEl->priorityValue = mei_ReadsPtr->pos + libInfo->maxDelta;
					break;
				}
				libInfo = libInfo->next;
			}
			if( max_delta_val == 0)
			{
				max_delta_val = 5 * SOFTCLIP_WRONGMAP_WINDOW;
				newEl->priorityValue = mei_ReadsPtr->pos + max_delta_val;
			}
			push_heap_mei( H_F, newEl);
			free( newEl);
		}
		mei_ReadsPtr = mei_ReadsPtr->next;
	}
}

void outputMEIClusters( parameters* params, char* chromosome_name)
{
	int MEIType;
	int MEIType2;
	int R_count, F_count;
	int count;
	char orientMEI;
	int written; /* Work around. Arda should fix this. Some clusters come out empty */
	clusters_final *cluster_new, *tmp;

	if( H_F->heapSize > 0 && H_R->heapSize > 0)
	{
		/* Iterate through all MEI types including the reverse orientations of each*/
		for( MEIType = 0; MEIType < (params->number_of_different_mei_types * 2); MEIType++)
		{
			clusters_all[cluster_count] = NULL;

			R_count = 0;
			F_count = 0;

			/* If the MEI is in forward direction, MEIType2 is in Reverse */
			if( ( MEIType % 2) == 0)
			{
				MEIType2 = MEIType + 1;
				orientMEI = MEIFORWARD;
			}
			else /* If the MEI is in reverse direction, MEIType2 is in Forward */
			{
				MEIType2 = MEIType - 1;
				orientMEI = MEIREVERSE;
			}

			/* Count the number of MEIs in H_F that matches the MEIType */
			for( count = 0; count < H_F->heapSize; count++)
			{
				if( H_F->heapArray[count].mei_ReadsPtr->MEI_Type == MEIType)
					F_count++;
			}
			/* Count the number of MEIs in H_R that matches the MEIType2 */
			for( count = 0; count < H_R->heapSize; count++)
			{
				if( H_R->heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2)
					R_count++;
			}

			if( F_count > 2 && R_count > 2)
			{
				written = 0;
				for( count = 0; count < H_F->heapSize; count++)
				{
					if( H_F->heapArray[count].mei_ReadsPtr->MEI_Type == MEIType)
					{
						if( debug_mode)
							fprintf( fileOutput, "%s %s %i %s %i %c %i %i %s %s F F %i %i ",
									H_F->heapArray[count].mei_ReadsPtr->readName,
									chromosome_name,
									H_F->heapArray[count].mei_ReadsPtr->pos,
									H_F->heapArray[count].mei_ReadsPtr->MEI_Subclass,
									H_F->heapArray[count].mei_ReadsPtr->pos,
									orientMEI,
									H_F->heapArray[count].mei_ReadsPtr->mQual,
									H_F->heapArray[count].mei_ReadsPtr->mQual,
									H_F->heapArray[count].mei_ReadsPtr->libName,
									H_F->heapArray[count].mei_ReadsPtr->indName,
									H_F->heapArray[count].mei_ReadsPtr->mQual,
									H_F->heapArray[count].mei_ReadsPtr->mQual);

						/* Fill the clusters struct */
						cluster_new = ( clusters_final *) getMem( sizeof( clusters_final));

						cluster_new->id = cluster_count;

						cluster_new->read_name = NULL;
						set_str( &cluster_new->read_name, H_F->heapArray[count].mei_ReadsPtr->readName);

						cluster_new->chromosome_name1 = NULL;
						set_str( &cluster_new->chromosome_name1, chromosome_name);

						cluster_new->start_position = H_F->heapArray[count].mei_ReadsPtr->pos;
						cluster_new->chromosome_name2 = NULL;

						cluster_new->end_position = H_F->heapArray[count].mei_ReadsPtr->pos;
						cluster_new->SV_type = orientMEI;
						cluster_new->phred_score = H_F->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->edit_distance = H_F->heapArray[count].mei_ReadsPtr->mQual;

						cluster_new->library_name = NULL;
						set_str( &cluster_new->library_name, H_F->heapArray[count].mei_ReadsPtr->libName);

						cluster_new->individual_name = NULL;
						set_str( &cluster_new->individual_name, H_F->heapArray[count].mei_ReadsPtr->indName);

						cluster_new->orientation_left = FORWARD;
						cluster_new->orientation_right = FORWARD;
						cluster_new->mapping_quality_left = H_F->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->mapping_quality_right = H_F->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->ten_x_barcode = H_F->heapArray[count].mei_ReadsPtr->ten_x_barcode;

						cluster_new->mei_type = NULL;
						set_str( &cluster_new->mei_type, H_F->heapArray[count].mei_ReadsPtr->MEI_Class);

						cluster_new->mei_subclass = NULL;
						set_str( &cluster_new->mei_subclass, H_F->heapArray[count].mei_ReadsPtr->MEI_Subclass);

						if( ( clusters_all[cluster_count] == NULL) || ( clusters_all[cluster_count]->edit_distance >= cluster_new->edit_distance))
						{
							cluster_new->next = clusters_all[cluster_count];
							clusters_all[cluster_count] = cluster_new;
						}
						else
						{
							tmp = clusters_all[cluster_count];
							while( ( tmp->next != NULL) && ( tmp->next->edit_distance < cluster_new->edit_distance))
								tmp = tmp->next;

							cluster_new->next = tmp->next;
							tmp->next = cluster_new;
						}
						written = 1;
					}
				}
				for( count = 0; count < H_R->heapSize; count++)
				{
					if( H_R->heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2)
					{
						if( debug_mode)
							fprintf( fileOutput, "%s %s %i %s %i %c %i %i %s %s R R %i %i ",
									H_R->heapArray[count].mei_ReadsPtr->readName,
									chromosome_name,
									H_R->heapArray[count].mei_ReadsPtr->pos,
									H_R->heapArray[count].mei_ReadsPtr->MEI_Subclass,
									H_R->heapArray[count].mei_ReadsPtr->pos,
									orientMEI,
									H_R->heapArray[count].mei_ReadsPtr->mQual,
									H_R->heapArray[count].mei_ReadsPtr->mQual,
									H_R->heapArray[count].mei_ReadsPtr->libName,
									H_R->heapArray[count].mei_ReadsPtr->indName,
									H_R->heapArray[count].mei_ReadsPtr->mQual,
									H_R->heapArray[count].mei_ReadsPtr->mQual);

						/* Fill the clusters struct */
						cluster_new = ( clusters_final *) getMem( sizeof( clusters_final));

						cluster_new->id = cluster_count;

						cluster_new->read_name = NULL;
						set_str( &cluster_new->read_name, H_R->heapArray[count].mei_ReadsPtr->readName);

						cluster_new->chromosome_name1 = NULL;
						set_str( &cluster_new->chromosome_name1, chromosome_name);

						cluster_new->start_position = H_R->heapArray[count].mei_ReadsPtr->pos;
						cluster_new->chromosome_name2 = NULL;

						cluster_new->end_position = H_R->heapArray[count].mei_ReadsPtr->pos;
						cluster_new->SV_type = orientMEI;
						cluster_new->phred_score = H_R->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->edit_distance = H_R->heapArray[count].mei_ReadsPtr->mQual;

						cluster_new->library_name = NULL;
						set_str( &cluster_new->library_name, H_R->heapArray[count].mei_ReadsPtr->libName);

						cluster_new->individual_name = NULL;
						set_str( &cluster_new->individual_name, H_R->heapArray[count].mei_ReadsPtr->indName);

						cluster_new->orientation_left = REVERSE;
						cluster_new->orientation_right = REVERSE;
						cluster_new->mapping_quality_left = H_R->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->mapping_quality_right = H_R->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->ten_x_barcode = H_R->heapArray[count].mei_ReadsPtr->ten_x_barcode;

						cluster_new->mei_type = NULL;
						set_str( &cluster_new->mei_type, H_R->heapArray[count].mei_ReadsPtr->MEI_Class);

						cluster_new->mei_subclass = NULL;
						set_str( &cluster_new->mei_subclass, H_R->heapArray[count].mei_ReadsPtr->MEI_Subclass);

						if( ( clusters_all[cluster_count] == NULL) || ( clusters_all[cluster_count]->edit_distance >= cluster_new->edit_distance))
						{
							cluster_new->next = clusters_all[cluster_count];
							clusters_all[cluster_count] = cluster_new;
						}
						else
						{
							tmp = clusters_all[cluster_count];
							while( ( tmp->next != NULL) && ( tmp->next->edit_distance < cluster_new->edit_distance))
								tmp = tmp->next;

							cluster_new->next = tmp->next;
							tmp->next = cluster_new;
						}
						written = 1;
					}
				}

				for( count = 0; count < H_S->heapSize; count++)
				{
					if( H_S->heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2 &&
							H_S->heapArray[count].mei_ReadsPtr->readName != NULL)
					{
						if( debug_mode)
						{
							fprintf(stderr,"%d MEITYPE%d MEITYPE2%d -", count, MEIType, MEIType2);
							fprintf( stderr, "%s ", H_S->heapArray[count].mei_ReadsPtr->readName);
							fprintf( stderr, "%s ", chromosome_name);
							fprintf( stderr, "%d ", H_S->heapArray[count].mei_ReadsPtr->pos);
							fprintf( stderr, "%s ", H_S->heapArray[count].mei_ReadsPtr->MEI_Subclass);
							fprintf( stderr, "%s ", H_S->heapArray[count].mei_ReadsPtr->indName);
							fprintf( stderr, "%d ", H_S->heapArray[count].mei_ReadsPtr->mQual);
							fprintf( stderr, "\n");

							fprintf( fileOutput, "%s %s %i %s %i %c %i %i SplitRead %s S S %i %i ",
									H_S->heapArray[count].mei_ReadsPtr->readName,
									chromosome_name,
									H_S->heapArray[count].mei_ReadsPtr->pos,
									H_S->heapArray[count].mei_ReadsPtr->MEI_Subclass,
									H_S->heapArray[count].mei_ReadsPtr->pos,
									orientMEI,
									H_S->heapArray[count].mei_ReadsPtr->mQual,
									H_S->heapArray[count].mei_ReadsPtr->mQual,
									H_S->heapArray[count].mei_ReadsPtr->indName,
									H_S->heapArray[count].mei_ReadsPtr->mQual,
									H_S->heapArray[count].mei_ReadsPtr->mQual);
						}
						/* Fill the clusters struct */
						cluster_new = ( clusters_final *) getMem( sizeof( clusters_final));

						cluster_new->read_name = NULL;
						set_str( &cluster_new->read_name, H_S->heapArray[count].mei_ReadsPtr->readName);

						cluster_new->chromosome_name1 = NULL;
						set_str( &cluster_new->chromosome_name1, chromosome_name);

						cluster_new->start_position = H_S->heapArray[count].mei_ReadsPtr->pos;
						cluster_new->chromosome_name2 = NULL;

						cluster_new->end_position = H_S->heapArray[count].mei_ReadsPtr->pos;
						cluster_new->SV_type = orientMEI;
						cluster_new->phred_score = H_S->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->edit_distance = H_S->heapArray[count].mei_ReadsPtr->mQual;

						cluster_new->library_name = NULL;
						set_str( &cluster_new->library_name, "SplitRead");

						cluster_new->individual_name = NULL;
						set_str( &cluster_new->individual_name, H_S->heapArray[count].mei_ReadsPtr->indName);

						cluster_new->orientation_left = 'S';
						cluster_new->orientation_right = 'S';
						cluster_new->mapping_quality_left = H_S->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->mapping_quality_right = H_S->heapArray[count].mei_ReadsPtr->mQual;
						cluster_new->ten_x_barcode = H_S->heapArray[count].mei_ReadsPtr->ten_x_barcode;

						cluster_new->mei_type = NULL;
						set_str( &cluster_new->mei_type, H_S->heapArray[count].mei_ReadsPtr->MEI_Class);

						cluster_new->mei_subclass = NULL;
						set_str( &cluster_new->mei_subclass, H_S->heapArray[count].mei_ReadsPtr->MEI_Subclass);

						if( ( clusters_all[cluster_count] == NULL) || ( clusters_all[cluster_count]->edit_distance >= cluster_new->edit_distance))
						{
							cluster_new->next = clusters_all[cluster_count];
							clusters_all[cluster_count] = cluster_new;
						}
						else
						{
							tmp = clusters_all[cluster_count];
							while( ( tmp->next != NULL) && ( tmp->next->edit_distance < cluster_new->edit_distance))
								tmp = tmp->next;

							cluster_new->next = tmp->next;
							tmp->next = cluster_new;
						}
						written = 1;
					}
				}
				if (written)
				{
					if( debug_mode)
						fprintf( fileOutput, "END\n");
					cluster_count++;
				}
			}
		}
	}
}

void MEICluster_Region( parameters* params, int chr_index)
{
	int leftBreakPoint;
	int boolMEITypeNewAdded = 0; // 0 or 1 indicates an new insertion of one of the 6 different types of MEI
	int boolMEITypeNewRemoved = 0; // 0 or 1 indicates a new deletion of one of the 6 different types of MEI

	for( leftBreakPoint = 0; leftBreakPoint < params->this_sonic->chromosome_lengths[chr_index]; leftBreakPoint++)
	{
		if( mReads[leftBreakPoint] != NULL)
		{
			add_F_Heap( leftBreakPoint);
			add_R_Heap( leftBreakPoint);
			add_S_Heap( leftBreakPoint);
			boolMEITypeNewAdded = 1;
		}

		if( ( ( H_R->heapSize > 0 && minValue_heapMEI( H_R) == leftBreakPoint) ||
				( H_F->heapSize > 0 && minValue_heapMEI( H_F) == leftBreakPoint) ||
				( H_S->heapSize > 0 && minValue_heapMEI( H_S) == leftBreakPoint)) && boolMEITypeNewAdded)
		{
			if( ( H_R->heapSize + H_F->heapSize + H_S->heapSize) > 0)
				outputMEIClusters( params, params->this_sonic->chromosome_names[chr_index]);

			boolMEITypeNewAdded = 0;
		}
		while( H_R->heapSize > 0 && minValue_heapMEI( H_R) == leftBreakPoint)
			heap_remove_topMEI( H_R);

		while( H_F->heapSize > 0 && minValue_heapMEI( H_F) == leftBreakPoint)
			heap_remove_topMEI( H_F);

		while( H_S->heapSize > 0 && minValue_heapMEI( H_S) == leftBreakPoint)
			heap_remove_topMEI( H_S);
	}
}

void vh_freeLinkedListMEI( mei_Reads * cur)
{
	mei_Reads *next;
	while (cur != NULL)
	{
		next = cur->next;
		free( cur->indName);
		free( cur->libName);
		free( cur->readName);
		free( cur->MEI_Class);
		free( cur->MEI_Subclass);
		free( cur);
		cur = next;
	}
}

void vh_finalizeReadMapping_Mei( int chroSize)
{
	int i;
	mei_Reads *ptr, *ptrNext;
	for( i = 0; i < chroSize; i++)
	{
		vh_freeLinkedListMEI(mReads[i]);
		mReads[i] = NULL;
	}

	if( mReads != NULL)
	{
		free( mReads);
		mReads = NULL;
	}
	vh_free_heap_mei(H_F);
	vh_free_heap_mei(H_R);
	vh_free_heap_mei(H_S);
}

int mei_regions( parameters *params, char* chromosome_name)
{
	int i, mei_count = 0;
	char *svtype;

	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;

	sonic_repeat *mei_left, *mei_right;

	libInfo = g_libInfo;

	while ( libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while ( ( divetReadMappingPtr != NULL) && ( strcmp( divetReadMappingPtr->chromosome_name, chromosome_name) == 0))
		{
			/* Check the mei table whether the right or left read is inside a mobile element */
			mei_left  = sonic_is_mobile_element( params->this_sonic, divetReadMappingPtr->chromosome_name, divetReadMappingPtr->locMapLeftStart, divetReadMappingPtr->locMapLeftEnd, params->mei);
			mei_right = sonic_is_mobile_element( params->this_sonic, divetReadMappingPtr->chromosome_name, divetReadMappingPtr->locMapRightStart, divetReadMappingPtr->locMapRightEnd, params->mei);

			if( mei_left != NULL)
			{
				divetReadMappingPtr->svType = MEILEFTPAIR;
				divetReadMappingPtr->mei_subclass = NULL;
				set_str( &(divetReadMappingPtr->mei_subclass), mei_left->repeat_type);

				divetReadMappingPtr->meiType = NULL;
				set_str( &(divetReadMappingPtr->meiType), mei_left->repeat_class);

				/*IF THE MEI INSERT IS + STRAND UPPERCASE IF - STRAND LOWER CASE*/
				if ( ( divetReadMappingPtr->orientationLeft == FORWARD && mei_left->strand == SONIC_STRAND_REV)
						|| ( divetReadMappingPtr->orientationLeft == REVERSE && mei_left->strand == SONIC_STRAND_FWD  ))
					divetReadMappingPtr->mei_code = (mei_left->mei_code * 2) + 1;
				else
					divetReadMappingPtr->mei_code = mei_left->mei_code * 2;
				mei_count++;
			}
			if( mei_right != NULL)
			{
				divetReadMappingPtr->svType = MEIRIGHTPAIR;
				divetReadMappingPtr->mei_subclass = NULL;
				set_str( &(divetReadMappingPtr->mei_subclass), mei_right->repeat_type);

				divetReadMappingPtr->meiType = NULL;
				set_str( &(divetReadMappingPtr->meiType), mei_right->repeat_class);

				/*IF THE MEI INSERT IS + STRAND UPPERCASE IF - STRAND LOWER CASE*/
				if ( ( divetReadMappingPtr->orientationRight == FORWARD && mei_right->strand == SONIC_STRAND_REV)
						|| ( divetReadMappingPtr->orientationRight == REVERSE && mei_right->strand == SONIC_STRAND_FWD ))
					divetReadMappingPtr->mei_code = (mei_right->mei_code * 2) + 1;
				else
					divetReadMappingPtr->mei_code = mei_right->mei_code * 2;
				mei_count++;
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
	}
	return mei_count;
}


void initializeReadMapping_MEI( bam_info** in_bams, parameters *params, int chr_index)
{
	int i, mei_count;

	mReads = (mei_Reads **) getMem( params->this_sonic->chromosome_lengths[chr_index] * sizeof( mei_Reads *));

	if( mReads == NULL)
		vh_logWarning ("Memory Problem in vh_createMaxClusterMEI");

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index]; i++)
		mReads[i] = NULL;

	H_F = vh_newHeapMEI(MAX_CLUSTER_SIZE);
	H_R = vh_newHeapMEI(MAX_CLUSTER_SIZE);
	H_S = vh_newHeapMEI(MAX_CLUSTER_SIZE);

	if( running_mode == QUICK)
		mei_count = addToGenomeIndex_MEI( in_bams, params, params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
	else
	{
		/* find the MEI sites */
		mei_count = mei_regions( params, params->this_sonic->chromosome_names[chr_index]);
		fprintf(logFile,"MEI regions count= %d\t", mei_count);
		mei_count = addToGenomeIndex_MEI_sensitive( params, params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		fprintf(logFile,"%d MEI in clusters\n", mei_count);
	}
}
