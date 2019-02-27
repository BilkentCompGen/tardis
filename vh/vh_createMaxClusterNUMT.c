/*
 * vh_createMaxClusterNUMT.c
 *
 *  Created on: Jul 31, 2018
 *      Author: tardis
 */

#include "vh_createMaxClusterNUMT.h"
#include "vh_common.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_setcover.h"
#include "vh_heap.h"
#include "../bamonly.h"

numt_Reads **numtReads = NULL;

HeapNUMT *Heap_F; /* The forward mapping heap from left of breakpoint */
HeapNUMT *Heap_R; /* the reverse mapping heap from right of breakpoint */

void vh_freeLinkedListNUMT( numt_Reads * cur)
{
	numt_Reads *next;
	while (cur != NULL)
	{
		next = cur->next;
		free( cur->indName);
		free( cur->libName);
		free( cur->readName);
		free( cur);
		cur = next;
	}
}

void vh_finalizeReadMapping_NUMT( int chroSize)
{
	int i;
	numt_Reads *ptr, *ptrNext;
	for( i = 0; i < chroSize; i++)
	{
		vh_freeLinkedListNUMT(numtReads[i]);
		numtReads[i] = NULL;
	}

	if( numtReads != NULL)
	{
		free( numtReads);
		numtReads = NULL;
	}
	vh_free_heap_numt(Heap_F);
	vh_free_heap_numt(Heap_R);
}

void outputNUMTClusters( parameters* params, char* chromosome_name)
{
	int count, numt_type, numt_type_alt;
	int R_count, F_count;
	int written;
	char orientation_NUMT;
	clusters_final *cluster_new, *tmp;

	//fprintf(stderr,"%d-%d, ", Heap_F->heapSize, Heap_R->heapSize);
	if( Heap_F->heapSize > 0 && Heap_R->heapSize > 0)
	{
		for( numt_type = 0; numt_type < 2; numt_type++)
		{
			clusters_all[cluster_count] = NULL;

			R_count = 0;
			F_count = 0;

			if( numt_type == NUMTF)
			{
				numt_type_alt = NUMTR;
				orientation_NUMT = NUMTFORWARD;
			}
			else
			{
				numt_type_alt = NUMTF;
				orientation_NUMT = NUMTREVERSE;
			}

			for( count = 0; count < Heap_F->heapSize; count++)
			{
				if( Heap_F->heapArray[count].numt_ReadsPtr->NUMT_Type == numt_type)
					F_count++;
			}
			for( count = 0; count < Heap_R->heapSize; count++)
			{
				if( Heap_R->heapArray[count].numt_ReadsPtr->NUMT_Type == numt_type_alt)
					R_count++;
			}
			if( F_count > 0 && R_count > 0)
			{
				written = 0;
				for( count = 0; count < Heap_F->heapSize; count++)
				{
					if( Heap_F->heapArray[count].numt_ReadsPtr->NUMT_Type == numt_type)
					{
						if( debug_mode)
							fprintf( fileOutput, "%s %s %i %c %i %i %s %s F F %i %i ",
									Heap_F->heapArray[count].numt_ReadsPtr->readName,
									chromosome_name,
									Heap_F->heapArray[count].numt_ReadsPtr->pos,
									orientation_NUMT,
									Heap_F->heapArray[count].numt_ReadsPtr->mQual,
									Heap_F->heapArray[count].numt_ReadsPtr->mQual,
									Heap_F->heapArray[count].numt_ReadsPtr->libName,
									Heap_F->heapArray[count].numt_ReadsPtr->indName,
									Heap_F->heapArray[count].numt_ReadsPtr->mQual,
									Heap_F->heapArray[count].numt_ReadsPtr->mQual);

						/* Fill the clusters struct */
						cluster_new = ( clusters_final *) getMem( sizeof( clusters_final));

						cluster_new->id = cluster_count;

						cluster_new->read_name = NULL;
						set_str( &cluster_new->read_name, Heap_F->heapArray[count].numt_ReadsPtr->readName);

						cluster_new->chromosome_name1 = NULL;
						set_str( &cluster_new->chromosome_name1, chromosome_name);

						cluster_new->start_position = Heap_F->heapArray[count].numt_ReadsPtr->pos;
						cluster_new->chromosome_name2 = NULL;

						cluster_new->end_position = Heap_F->heapArray[count].numt_ReadsPtr->pos;
						cluster_new->SV_type = orientation_NUMT;
						cluster_new->phred_score = Heap_F->heapArray[count].numt_ReadsPtr->mQual;
						cluster_new->edit_distance = Heap_F->heapArray[count].numt_ReadsPtr->mQual;

						cluster_new->library_name = NULL;
						set_str( &cluster_new->library_name, Heap_F->heapArray[count].numt_ReadsPtr->libName);

						cluster_new->individual_name = NULL;
						set_str( &cluster_new->individual_name, Heap_F->heapArray[count].numt_ReadsPtr->indName);

						cluster_new->orientation_left = FORWARD;
						cluster_new->orientation_right = FORWARD;
						cluster_new->mapping_quality_left = Heap_F->heapArray[count].numt_ReadsPtr->mQual;
						cluster_new->mapping_quality_right = Heap_F->heapArray[count].numt_ReadsPtr->mQual;
						cluster_new->ten_x_barcode = Heap_F->heapArray[count].numt_ReadsPtr->ten_x_barcode;

						cluster_new->mei_type = NULL;
						cluster_new->mei_subclass = NULL;

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
				for( count = 0; count < Heap_R->heapSize; count++)
				{
					if( Heap_R->heapArray[count].numt_ReadsPtr->NUMT_Type == numt_type_alt)
					{
						if( debug_mode)
							fprintf( fileOutput, "%s %s %i %c %i %i %s %s R R %i %i ",
									Heap_R->heapArray[count].numt_ReadsPtr->readName,
									chromosome_name,
									Heap_R->heapArray[count].numt_ReadsPtr->pos,
									orientation_NUMT,
									Heap_R->heapArray[count].numt_ReadsPtr->mQual,
									Heap_R->heapArray[count].numt_ReadsPtr->mQual,
									Heap_R->heapArray[count].numt_ReadsPtr->libName,
									Heap_R->heapArray[count].numt_ReadsPtr->indName,
									Heap_R->heapArray[count].numt_ReadsPtr->mQual,
									Heap_R->heapArray[count].numt_ReadsPtr->mQual);

						/* Fill the clusters struct */
						cluster_new = ( clusters_final *) getMem( sizeof( clusters_final));

						cluster_new->id = cluster_count;

						cluster_new->read_name = NULL;
						set_str( &cluster_new->read_name, Heap_R->heapArray[count].numt_ReadsPtr->readName);

						cluster_new->chromosome_name1 = NULL;
						set_str( &cluster_new->chromosome_name1, chromosome_name);

						cluster_new->chromosome_name2 = NULL;

						cluster_new->start_position = Heap_R->heapArray[count].numt_ReadsPtr->pos;

						cluster_new->end_position = Heap_R->heapArray[count].numt_ReadsPtr->pos;
						cluster_new->SV_type = orientation_NUMT;
						cluster_new->phred_score = Heap_R->heapArray[count].numt_ReadsPtr->mQual;
						cluster_new->edit_distance = Heap_R->heapArray[count].numt_ReadsPtr->mQual;

						cluster_new->library_name = NULL;
						set_str( &cluster_new->library_name, Heap_R->heapArray[count].numt_ReadsPtr->libName);

						cluster_new->individual_name = NULL;
						set_str( &cluster_new->individual_name, Heap_R->heapArray[count].numt_ReadsPtr->indName);

						cluster_new->orientation_left = REVERSE;
						cluster_new->orientation_right = REVERSE;
						cluster_new->mapping_quality_left = Heap_R->heapArray[count].numt_ReadsPtr->mQual;
						cluster_new->mapping_quality_right = Heap_R->heapArray[count].numt_ReadsPtr->mQual;
						cluster_new->ten_x_barcode = Heap_R->heapArray[count].numt_ReadsPtr->ten_x_barcode;

						cluster_new->mei_type = NULL;
						cluster_new->mei_subclass = NULL;

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


void add_R_Heap_Numt( int pos)
{
	numt_Reads *numt_ReadsPtr;
	numt_ReadsPtr = numtReads[pos];
	HeapElNUMT *newEl;

	while(numt_ReadsPtr != NULL)
	{
		if (numt_ReadsPtr->orient == REVERSE)
		{
		  newEl = (HeapElNUMT *)getMem( sizeof( HeapElNUMT));

			newEl->numt_ReadsPtr = numt_ReadsPtr;
			newEl->priorityValue = numt_ReadsPtr->pos;
			//fprintf(stderr,"ADDING REVERSE - Pos %d, PRIORITY %d, breakpoint %d\n", numt_ReadsPtr->pos, newEl->priorityValue, pos);
			push_heap_numt(Heap_R, newEl);
			free(newEl);
		}
		numt_ReadsPtr = numt_ReadsPtr->next;
	}
}

void add_F_Heap_Numt( int pos)
{
	LibraryInfo *libInfo;
	numt_Reads *numt_ReadsPtr;
	numt_ReadsPtr = numtReads[pos];
	HeapElNUMT *newEl;
	int max_delta_val;

	while( numt_ReadsPtr != NULL)
	{
		if( numt_ReadsPtr->orient == FORWARD)
		{
		       newEl = ( HeapElNUMT *) getMem( sizeof( HeapElNUMT));
			newEl->numt_ReadsPtr = numt_ReadsPtr;

			libInfo = g_libInfo;
			max_delta_val = 0;
			while( libInfo != NULL)
			{
				if( strcmp( numt_ReadsPtr->libName, libInfo->libName) == 0 && strcmp( numt_ReadsPtr->indName, libInfo->indName) == 0)
				{
					max_delta_val = libInfo->maxDelta;
					newEl->priorityValue = numt_ReadsPtr->pos + max_delta_val + ( 2 * libInfo->readLen);
					break;
				}
				libInfo = libInfo->next;
			}
			if( max_delta_val == 0)
			{
				max_delta_val = 5 * SOFTCLIP_WRONGMAP_WINDOW;
				newEl->priorityValue = numt_ReadsPtr->pos + max_delta_val;
			}
			//fprintf(stderr,"ADDING FORWARD - Pos %d, PRIORITY %d, breakpoint %d\n", numt_ReadsPtr->pos, newEl->priorityValue, pos);
			push_heap_numt( Heap_F, newEl);
			free( newEl);
		}
		numt_ReadsPtr = numt_ReadsPtr->next;
	}
}


void NUMTCluster_Region( parameters* params, int chr_index)
{
	int leftBreakPoint;
	int boolNUMTTypeNewAdded = 0; // 0 or 1 indicates a new insertion NUMT

	for( leftBreakPoint = 0; leftBreakPoint < params->this_sonic->chromosome_lengths[chr_index]; leftBreakPoint++)
	{
		if( numtReads[leftBreakPoint] != NULL)
		{
			add_F_Heap_Numt( leftBreakPoint);
			add_R_Heap_Numt( leftBreakPoint);

			boolNUMTTypeNewAdded = 1;
		}

		if( ( ( Heap_R->heapSize > 0 && minValue_heapNUMT( Heap_R) == leftBreakPoint) ||
				( Heap_F->heapSize > 0 && minValue_heapNUMT( Heap_F) == leftBreakPoint)) && boolNUMTTypeNewAdded)
		{
			if( ( Heap_R->heapSize + Heap_F->heapSize) > 0)
				outputNUMTClusters( params, params->this_sonic->chromosome_names[chr_index]);

			boolNUMTTypeNewAdded = 0;
		}
		while( Heap_R->heapSize > 0 && minValue_heapNUMT( Heap_R) == leftBreakPoint)
			heap_remove_topNUMT( Heap_R);

		while( Heap_F->heapSize > 0 && minValue_heapNUMT( Heap_F) == leftBreakPoint)
			heap_remove_topNUMT( Heap_F);
	}
}

int addToGenomeIndex_NUMT (bam_info** in_bams, parameters *params, int chr_index)
{
	LibraryInfo *libInfo;
	numt_Reads *tmpNUMT_Reads;
	discordantMappingNUMT *discordantReadPtr;
	int numSample, count, numt_count = 0, tmp;
	float is_satellite;

	int libId = 0;

	for( numSample = 0; numSample < params->num_bams; numSample++)
	{
		for( count = 0; count < in_bams[numSample]->num_libraries; count++)
		{
			discordantReadPtr = in_bams[numSample]->libraries[count]->listNUMT_Mapping;

			while( discordantReadPtr != NULL)
			{
				is_satellite = sonic_is_satellite (params->this_sonic, discordantReadPtr->chromosome_name, discordantReadPtr->pos, discordantReadPtr->pos_End);
				if( is_satellite == 0 && strcmp( discordantReadPtr->chromosome_name, params->this_sonic->chromosome_names[chr_index]) == 0
						&& ( discordantReadPtr->pos > 0) && ( discordantReadPtr->pos_End < params->this_sonic->chromosome_lengths[chr_index]))
				{
					if( discordantReadPtr->NUMT_Type > -1)
					{
						//fprintf(stderr, "%s %d-%d %c %d %d\n",discordantReadPtr->readName, discordantReadPtr->pos, discordantReadPtr->pos_End, discordantReadPtr->orientation,
						//	discordantReadPtr->qual, discordantReadPtr->NUMT_Type);
						tmpNUMT_Reads = (numt_Reads*) getMem( sizeof( numt_Reads));

						tmpNUMT_Reads->readName = NULL;
						set_str( &(tmpNUMT_Reads->readName), discordantReadPtr->readName);

						tmpNUMT_Reads->pos = discordantReadPtr->pos;
						tmpNUMT_Reads->pos_End = discordantReadPtr->pos_End;

						tmpNUMT_Reads->orient = discordantReadPtr->orientation;
						tmpNUMT_Reads->mQual = discordantReadPtr->qual;
						tmpNUMT_Reads->NUMT_Type = discordantReadPtr->NUMT_Type;
						tmpNUMT_Reads->libId = libId;
						tmpNUMT_Reads->ten_x_barcode = discordantReadPtr->ten_x_barcode;

						tmpNUMT_Reads->libName = NULL;
						set_str( &(tmpNUMT_Reads->libName), in_bams[numSample]->libraries[count]->libname);

						tmpNUMT_Reads->indName = NULL;
						set_str( &tmpNUMT_Reads->indName, in_bams[numSample]->sample_name);

						if (discordantReadPtr->orientation == FORWARD)
						{
							tmpNUMT_Reads->next = numtReads[discordantReadPtr->pos_End];
							numtReads[discordantReadPtr->pos_End] = tmpNUMT_Reads;
						}
						if (discordantReadPtr->orientation == REVERSE)
						{
							//tmp = max( 0, discordantReadPtr->pos);  -
							//(2 * in_bams[numSample]->libraries[count]->read_length))
							tmp = max( 0, discordantReadPtr->pos - (in_bams[numSample]->libraries[count]->conc_max));
							tmpNUMT_Reads->next = numtReads[tmp];
							numtReads[tmp] = tmpNUMT_Reads;
						}
						numt_count++;
					}
				}
				discordantReadPtr = discordantReadPtr->next;
			}
			libId++;
		}
	}
	return numt_count;
}

void initializeReadMapping_NUMT( bam_info** in_bams, parameters *params, int chr_index)
{
	int i, numt_count;

	numtReads = (numt_Reads **) getMem( params->this_sonic->chromosome_lengths[chr_index] * sizeof( numt_Reads *));

	if( numtReads == NULL)
		vh_logWarning ("Memory Problem in vh_createMaxClusterNUMT");

	
	
	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index]; i++)
		numtReads[i] = NULL;

	Heap_F = vh_newHeapNUMT(MAX_CLUSTER_SIZE);
	Heap_R = vh_newHeapNUMT(MAX_CLUSTER_SIZE);

	numt_count = addToGenomeIndex_NUMT( in_bams, params, chr_index);
}
