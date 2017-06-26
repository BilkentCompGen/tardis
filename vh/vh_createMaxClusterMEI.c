#include "vh_createMaxClusterMEI.h"
#include "vh_common.h"
#include "vh_intervalhandler.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_setcover.h"
#include "vh_heap.h"

mei_Reads **mReads;
HeapMEI H_F; /* The forward mapping heap from left of breakpoint */
HeapMEI H_R; /* the reverse mapping heap from right of breakpoint */
HeapMEI H_S; /* The soft clipping read */

void addToGenomeIndex_MEI (bam_info** in_bam, ref_genome* ref, parameters *params, char *chromosome_name)
{
	mei_Reads *tmpMEI_Reads;
	discordantMappingMEI *discordantReadPtr;
	softClip *softClipPtr;
	int numSample, count;
	int is_satellite;
	int countRR = 0, countFF = 0, countTotal = 0;

	sonic_repeat *mei;
	int libId = 0;

	for( numSample = 0; numSample < params->num_bams; numSample++)
	{
		for( count = 0; count < in_bam[numSample]->num_libraries; count++)
		{
			discordantReadPtr = in_bam[numSample]->libraries[count]->listMEI_Mapping;
			softClipPtr = in_bam[numSample]->libraries[count]->listSoftClip;

			while( discordantReadPtr != NULL)
			{
			  is_satellite = sonic_is_satellite (params->this_sonic, discordantReadPtr->chromosome_name, discordantReadPtr->pos, discordantReadPtr->pos_End);
			  if( is_satellite == 0 && strcmp(discordantReadPtr->chromosome_name, chromosome_name) == 0)
				{
					if( discordantReadPtr->MEI_Type > -1)
					{
						tmpMEI_Reads = (mei_Reads*) getMem( sizeof( mei_Reads));

						tmpMEI_Reads->readName = NULL;
						set_str( &(tmpMEI_Reads->readName), discordantReadPtr->readName);

						tmpMEI_Reads->MEI_Subclass = NULL;
						set_str( &(tmpMEI_Reads->MEI_Subclass), discordantReadPtr->MEI_subclass);


						tmpMEI_Reads->pos = discordantReadPtr->pos;
						tmpMEI_Reads->pos_End = discordantReadPtr->pos_End;
						tmpMEI_Reads->orient = discordantReadPtr->orient;
						tmpMEI_Reads->mQual = discordantReadPtr->qual;
						tmpMEI_Reads->MEI_Type = discordantReadPtr->MEI_Type;
						tmpMEI_Reads->readTypeSupport = 0; // The support is type paired-end read mapping
						tmpMEI_Reads->libId = libId;
	    
						tmpMEI_Reads->libName = NULL;
						set_str( &(tmpMEI_Reads->libName), in_bam[numSample]->libraries[count]->libname);

						tmpMEI_Reads->indName = NULL;
						set_str( &tmpMEI_Reads->indName, in_bam[numSample]->sample_name);

						if (discordantReadPtr->orient=='F')
						{
							countFF++;
							tmpMEI_Reads->next = mReads[discordantReadPtr->pos_End];
							mReads[discordantReadPtr->pos_End] = tmpMEI_Reads;
						}
						else if (discordantReadPtr->orient == 'R')
						{
							countRR++;
							tmpMEI_Reads->next = mReads[max(0, discordantReadPtr->pos - in_bam[numSample]->libraries[count]->conc_max)];
							mReads[max(0, discordantReadPtr->pos - in_bam[numSample]->libraries[count]->conc_max)] = tmpMEI_Reads;
						}
					}
				}
				discordantReadPtr = discordantReadPtr->next;
			}
			while( softClipPtr != NULL)
			{
			  is_satellite = sonic_is_satellite (params->this_sonic, softClipPtr->chromosome_name, softClipPtr->pos, softClipPtr->pos+1);			  
			  if( is_satellite == 0 && strcmp( softClipPtr->chromosome_name, chromosome_name) == 0)
			    {
					posMapSoftClip *ptrPosMapSoftClip = softClipPtr->ptrPosMapSoftClip;
					if( ptrPosMapSoftClip != NULL)
					{

					  mei = sonic_is_mobile_element( params->this_sonic, chromosome_name, ptrPosMapSoftClip->posMap, ptrPosMapSoftClip->posMap+1, params->mei);
						if( mei != NULL)
						{
							tmpMEI_Reads = (mei_Reads*) getMem(sizeof(mei_Reads));

							tmpMEI_Reads->readName = NULL;
							set_str( &(tmpMEI_Reads->readName), softClipPtr->readName);
							
							tmpMEI_Reads->MEI_Subclass = NULL;
							set_str( &(tmpMEI_Reads->MEI_Subclass), mei->repeat_type);


							tmpMEI_Reads->pos = softClipPtr->pos;
							tmpMEI_Reads->orient = softClipPtr->orient;
							tmpMEI_Reads->mQual = softClipPtr->qual;
							tmpMEI_Reads->libId = libId;

							tmpMEI_Reads->libName = NULL;
							set_str( &(tmpMEI_Reads->libName), in_bam[numSample]->libraries[count]->libname);

							tmpMEI_Reads->indName = NULL;
							set_str( &(tmpMEI_Reads->indName), in_bam[numSample]->sample_name);


							if( softClipPtr->op[0] == 4) //soft clipped at the beginning of the read. Horrible. Replace with #define
							{
								tmpMEI_Reads->readTypeSupport = 1;
								tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos-lenSplitReadBrakWindow];
								mReads[tmpMEI_Reads->pos - lenSplitReadBrakWindow] = tmpMEI_Reads;
							}
							else if( softClipPtr->op[softClipPtr->opCount - 1] == 4) // soft cliped at the end of the read. Horrible. Replace with #define
							{
								tmpMEI_Reads->readTypeSupport=2;
								tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos + softClipPtr->opl[0] - lenSplitReadBrakWindow];
								mReads[tmpMEI_Reads->pos + softClipPtr->opl[0] - lenSplitReadBrakWindow] = tmpMEI_Reads;
							}
						}
					}
				}
				softClipPtr = softClipPtr->next;
			}
			libId++;
		}
	}
}



void vh_addToGenomeIndex_MEI (bam_info** in_bams, parameters *params, char *chromosome_name)
{
	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	mei_Reads *tmpMEI_Reads;
	int countRR = 0, countFF = 0, countTotal = 0, libId = 0;

	libInfo = g_libInfo;
	while (libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while (divetReadMappingPtr != NULL)
		{
			if (strcmp (divetReadMappingPtr->chromosome_name, chromosome_name) == 0 && ((divetReadMappingPtr->svType == 'M' || divetReadMappingPtr->svType == 'X')))
			{
				tmpMEI_Reads=(mei_Reads *)getMem(sizeof(mei_Reads));
				tmpMEI_Reads->readName = NULL;
				set_str( &(tmpMEI_Reads->readName), divetReadMappingPtr->readName->readName);

				tmpMEI_Reads->mQual=divetReadMappingPtr->phredScore;
				tmpMEI_Reads->readTypeSupport=0; // The support is type paired-end read mapping

				tmpMEI_Reads->libName = NULL;
				set_str( &(tmpMEI_Reads->libName), libInfo->libName);

				tmpMEI_Reads->indName = NULL;
				set_str( &(tmpMEI_Reads->indName), libInfo->indName);

				/* Modify to get subclass name instead */
				tmpMEI_Reads->MEI_Subclass = NULL; 
				set_str( &(tmpMEI_Reads->MEI_Subclass), divetReadMappingPtr->mei_subclass);


				tmpMEI_Reads->libId=libId;

				/* ARDA: Horrible coding. Fix. */
				/* sonic_repeat now has mei_code that points to the index of the MEI type in the --mei string */
				/* just use it, instead of hardcoding meiType[]. For reverse complements, add n = number_of_mei_types */
				if(divetReadMappingPtr->meiType[0]=='A')
					tmpMEI_Reads->MEI_Type=0;
				else if (divetReadMappingPtr->meiType[0]=='a')
					tmpMEI_Reads->MEI_Type=1;
				else if(divetReadMappingPtr->meiType[0]=='L')
					tmpMEI_Reads->MEI_Type=2;
				else if (divetReadMappingPtr->meiType[0]=='l')
					tmpMEI_Reads->MEI_Type=3;
				else if(divetReadMappingPtr->meiType[0]=='S')
					tmpMEI_Reads->MEI_Type=4;
				else if (divetReadMappingPtr->meiType[0]=='s')
					tmpMEI_Reads->MEI_Type=5;

				if (divetReadMappingPtr->svType == 'M')
				{
					if(divetReadMappingPtr->orientationLeft=='F')
					{
						tmpMEI_Reads->pos=divetReadMappingPtr->locMapLeftStart;
						tmpMEI_Reads->pos_End=divetReadMappingPtr->locMapLeftEnd;
						tmpMEI_Reads->orient='F';
						countFF++;
						tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos_End];
						mReads[tmpMEI_Reads->pos_End]=tmpMEI_Reads;
					}
					else if (divetReadMappingPtr->orientationLeft=='R')
					{
						tmpMEI_Reads->pos=divetReadMappingPtr->locMapLeftStart;
						tmpMEI_Reads->pos_End=divetReadMappingPtr->locMapLeftEnd;
						tmpMEI_Reads->orient='R';
						countRR++;
						tmpMEI_Reads->next=mReads[max(0,tmpMEI_Reads->pos-libInfo->maxDelta)];
						mReads[max(0, tmpMEI_Reads->pos-libInfo->maxDelta)]=tmpMEI_Reads;
					}
				}
				else if (divetReadMappingPtr->svType == 'X')
				{
					if(divetReadMappingPtr->orientationRight=='F')
					{
						tmpMEI_Reads->pos=divetReadMappingPtr->locMapRightStart;
						tmpMEI_Reads->pos_End=divetReadMappingPtr->locMapRightEnd;
						tmpMEI_Reads->orient='F';
						countFF++;
						tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos_End];
						mReads[tmpMEI_Reads->pos_End]=tmpMEI_Reads;
					}
					else if (divetReadMappingPtr->orientationRight=='R')
					{
						tmpMEI_Reads->pos=divetReadMappingPtr->locMapRightStart;
						tmpMEI_Reads->pos_End=divetReadMappingPtr->locMapRightEnd;
						tmpMEI_Reads->orient='R';
						countRR++;
						tmpMEI_Reads->next = mReads[max(0,tmpMEI_Reads->pos-libInfo->maxDelta)];
						mReads[max(0, tmpMEI_Reads->pos-libInfo->maxDelta)]=tmpMEI_Reads;
					}
				}
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
		libId++;
	}
}




void add_R_Heap(int pos)
{
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapElMEI *newEl;
	while(mei_ReadsPtr!=NULL)
	{
		if (mei_ReadsPtr->orient=='R' && mei_ReadsPtr->readTypeSupport==0)
		{
			newEl=(HeapElMEI *)getMem(sizeof(HeapElMEI));
			newEl->mei_ReadsPtr = mei_ReadsPtr;
			newEl->priorityValue = mei_ReadsPtr->pos;
			push_heap_mei(&H_R, newEl);
			free(newEl);
		}
		mei_ReadsPtr=mei_ReadsPtr->next;
	}
}

void add_S_Heap(int pos)
{
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapElMEI *newEl;
	while(mei_ReadsPtr!=NULL)
	{
		if ( mei_ReadsPtr->readTypeSupport==1 || mei_ReadsPtr->readTypeSupport==2)
		{
			newEl=(HeapElMEI *)getMem(sizeof(HeapElMEI));
			newEl->mei_ReadsPtr = mei_ReadsPtr;
			newEl->priorityValue = pos+2*lenSplitReadBrakWindow;
			push_heap_mei(& H_S, newEl);
			free(newEl);
		}
		mei_ReadsPtr=mei_ReadsPtr->next;
	}
}

void add_F_Heap(bam_info** in_bams, int pos)
{
	LibraryInfo *libInfo;
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapElMEI *newEl;
	int max_delta_val;

	while(mei_ReadsPtr!=NULL)
	{
		if (mei_ReadsPtr->orient=='F' && mei_ReadsPtr->readTypeSupport==0)
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
			push_heap_mei( & H_F, newEl);
			free( newEl);
		}
		mei_ReadsPtr = mei_ReadsPtr->next;
	}
}

void outputMEIClusters(char* chromosome_name)
{
	int MEIType;
	int MEIType2;
	int R_count, F_count;
	int count;
	int orientMEI;
	int written; /* Work around. Arda should fix this. Some clusters come out empty */
	
	if( H_F.heapSize > 0  && H_R.heapSize > 0)
	{
		for( MEIType = 0; MEIType < 6; MEIType++)
		{
			R_count = 0;
			F_count = 0;

			/* match reads with the same type of MEI (Alu, etc) and the same orientation of mapping
			 *  so if it is H_R assuming the read is in R orientation then mapping
			 * inside a Alu with R orientation then the insertion is type F*/

			/* Horrible. Impossible to understand what is going on here */
			if( ( MEIType % 2) == 0)
				MEIType2 = MEIType + 1;
			else
				MEIType2 = MEIType - 1;
			if( ( MEIType % 2) == 0)
				orientMEI = 10;
			else
			  orientMEI = 12;

			for( count = 0; count < H_F.heapSize; count++)
			{
				if( H_F.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType)
					F_count++;
			}
			for( count = 0; count < H_R.heapSize; count++)
			{
				if( H_R.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2)
					R_count++;
			}

			if( F_count > 2 && R_count > 2)
			{
			        written = 0;
  			        for( count = 0; count < H_F.heapSize; count++)
				  {
					if( H_F.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType && MEIType > -1 &&  MEIType < 2)
					{
						fprintf( fileOutput, "%s %s %i %s %i %i %i %i %s %s F F %i %i ",
								H_F.heapArray[count].mei_ReadsPtr->readName,
								chromosome_name,
								H_F.heapArray[count].mei_ReadsPtr->pos,
								H_F.heapArray[count].mei_ReadsPtr->MEI_Subclass,
								H_F.heapArray[count].mei_ReadsPtr->pos,
								orientMEI,
								H_F.heapArray[count].mei_ReadsPtr->mQual,
								H_F.heapArray[count].mei_ReadsPtr->mQual,
								H_F.heapArray[count].mei_ReadsPtr->libName,
								H_F.heapArray[count].mei_ReadsPtr->indName,
								H_F.heapArray[count].mei_ReadsPtr->mQual,
								H_F.heapArray[count].mei_ReadsPtr->mQual);
						written = 1;
					}
					else if( H_F.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType && MEIType > 1 &&  MEIType < 4)
					{
						fprintf( fileOutput, "%s %s %i %s %i %i %i %i %s %s F F %i %i ",
								H_F.heapArray[count].mei_ReadsPtr->readName,
								chromosome_name,
								H_F.heapArray[count].mei_ReadsPtr->pos,
								H_F.heapArray[count].mei_ReadsPtr->MEI_Subclass,
								H_F.heapArray[count].mei_ReadsPtr->pos,
								orientMEI,
								H_F.heapArray[count].mei_ReadsPtr->mQual,
								H_F.heapArray[count].mei_ReadsPtr->mQual,
								H_F.heapArray[count].mei_ReadsPtr->libName,
								H_F.heapArray[count].mei_ReadsPtr->indName,
								H_F.heapArray[count].mei_ReadsPtr->mQual,
								H_F.heapArray[count].mei_ReadsPtr->mQual);
						written = 1;
					}
				}
				for( count = 0; count < H_R.heapSize; count++)
				{
					if( H_R.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2 && MEIType > -1 &&  MEIType < 2)
					{
						fprintf( fileOutput, "%s %s %i %s %i %i %i %i %s %s R R %i %i ",
								H_R.heapArray[count].mei_ReadsPtr->readName,
								chromosome_name,
								H_R.heapArray[count].mei_ReadsPtr->pos,
								H_R.heapArray[count].mei_ReadsPtr->MEI_Subclass,
								H_R.heapArray[count].mei_ReadsPtr->pos,
								orientMEI,
								H_R.heapArray[count].mei_ReadsPtr->mQual,
								H_R.heapArray[count].mei_ReadsPtr->mQual,
								H_R.heapArray[count].mei_ReadsPtr->libName,
								H_R.heapArray[count].mei_ReadsPtr->indName,
								H_R.heapArray[count].mei_ReadsPtr->mQual,
								H_R.heapArray[count].mei_ReadsPtr->mQual);
						written = 1;
					}
					else if( H_R.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2 && MEIType > 1 &&  MEIType < 4)
					{
						fprintf( fileOutput, "%s %s %i %s %i %i %i %i %s %s R R %i %i ",
								H_R.heapArray[count].mei_ReadsPtr->readName,
								chromosome_name,
								H_R.heapArray[count].mei_ReadsPtr->pos,
								H_R.heapArray[count].mei_ReadsPtr->MEI_Subclass,
								H_R.heapArray[count].mei_ReadsPtr->pos,
								orientMEI,
								H_R.heapArray[count].mei_ReadsPtr->mQual,
								H_R.heapArray[count].mei_ReadsPtr->mQual,
								H_R.heapArray[count].mei_ReadsPtr->libName,
								H_R.heapArray[count].mei_ReadsPtr->indName,
								H_R.heapArray[count].mei_ReadsPtr->mQual,
								H_R.heapArray[count].mei_ReadsPtr->mQual);
						written = 1;
					}
				}
				for( count = 0; count < H_S.heapSize; count++)
				{
					if( H_S.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2 && MEIType > -1 &&  MEIType < 2)
					  {
						fprintf( fileOutput, "%s %s %i %s %i %i %i %i SplitRead %s S S %i %i ",
								H_S.heapArray[count].mei_ReadsPtr->readName,
								chromosome_name,
								H_S.heapArray[count].mei_ReadsPtr->pos,
								H_S.heapArray[count].mei_ReadsPtr->MEI_Subclass,
								H_S.heapArray[count].mei_ReadsPtr->pos,
								orientMEI,
								H_S.heapArray[count].mei_ReadsPtr->mQual,
								H_S.heapArray[count].mei_ReadsPtr->mQual,
								H_S.heapArray[count].mei_ReadsPtr->indName,
								H_S.heapArray[count].mei_ReadsPtr->mQual,
								H_S.heapArray[count].mei_ReadsPtr->mQual);
						written = 1;
					  }
					else if( H_S.heapArray[count].mei_ReadsPtr->MEI_Type == MEIType2 && MEIType > 1 &&  MEIType < 4)
					  {
					        fprintf( fileOutput, "%s %s %i %s %i %i %i %i SplitRead %s S S %i %i ",
								H_S.heapArray[count].mei_ReadsPtr->readName,
								chromosome_name,
								H_S.heapArray[count].mei_ReadsPtr->pos,
								H_S.heapArray[count].mei_ReadsPtr->MEI_Subclass,
								H_S.heapArray[count].mei_ReadsPtr->pos,
								orientMEI,
								H_S.heapArray[count].mei_ReadsPtr->mQual,
								H_S.heapArray[count].mei_ReadsPtr->mQual,
								H_S.heapArray[count].mei_ReadsPtr->indName,
								H_S.heapArray[count].mei_ReadsPtr->mQual,
								H_S.heapArray[count].mei_ReadsPtr->mQual);
						written = 1;
					  }
				}
				if (written){
				  fprintf( fileOutput, "END\n");
				  cluster_count++;
				}
			}
		}
	}
}

void MEICluster_Region(bam_info** in_bams, char* chromosome_name, int chroSize)
{
	int brkPointIndex=0;
	int boolMEITypeNewAdded=0; // 0 or 1 indicates an new insertion of one of the 6 different types of MEI
	int boolMEITypeNewRemoved=0; // 0 or 1 indicates a new deletion of one of the 6 different types of MEI

	for (brkPointIndex=0; brkPointIndex < chroSize; brkPointIndex++)
	{
		if (mReads[brkPointIndex]!=NULL)
		{
			add_F_Heap(in_bams, brkPointIndex);
			add_R_Heap(brkPointIndex);
			add_S_Heap(brkPointIndex);
			boolMEITypeNewAdded=1;
		}
		if (((H_R.heapSize>0 && minValue_heapMEI(&H_R) == brkPointIndex) || (minValue_heapMEI(&H_F) == brkPointIndex && H_F.heapSize>0) || (minValue_heapMEI(&H_S)==brkPointIndex && H_S.heapSize>0)) && boolMEITypeNewAdded)
		{
			if ((H_R.heapSize+H_F.heapSize + H_S.heapSize)>0)
			{
				outputMEIClusters(chromosome_name);
			}
			boolMEITypeNewAdded=0;
		}
		while(H_R.heapSize>0 && minValue_heapMEI(&H_R)==brkPointIndex)
		{
			heap_remove_topMEI(&H_R);
		}
		while(H_F.heapSize>0 && minValue_heapMEI(&H_F)==brkPointIndex)
		{
			heap_remove_topMEI(&H_F);
		}
		while(H_S.heapSize>0 && minValue_heapMEI(&H_S)==brkPointIndex)
		{
			heap_remove_topMEI(&H_S);
		}
	}
}

void vh_freeLinkedListMEI (mei_Reads * cur)
{
	mei_Reads *next;
	while (cur != NULL)
	{
		next = cur->next;
		free( cur->indName);
		free( cur->libName);
		free( cur->readName);
		free( cur->MEI_Subclass);
		free( cur);
		cur = next;
	}
}

void vh_finalizeReadMapping_Mei( int chroSize)
{
	int count;
	mei_Reads *ptr, *ptrNext;
	for( count = 0; count <= chroSize; count++)
	{
		vh_freeLinkedListMEI( mReads[count]);
		mReads[count] = NULL;
	}
	if( mReads != NULL)
	{
		free( mReads);
		mReads = NULL;
	}
	freeHeapMEI( &H_S);
	freeHeapMEI( &H_R);
	freeHeapMEI( &H_F);
}


void initializeReadMapping_MEI( bam_info** in_bams, ref_genome* ref, parameters *params, char *chromosome_name, int chroSize)
{
	int count;
	mReads = (mei_Reads **) getMem( (chroSize + 1) * sizeof( mei_Reads *));

	if( mReads == NULL)
		vh_logWarning ("Memory Problem in vh_createMaxClusterMEI");

	for( count = 0; count <= chroSize; count++)
		mReads[count] = NULL;
	addToGenomeIndex_MEI( in_bams, ref, params, chromosome_name);
}


void vh_initializeReadMapping_MEI (bam_info** in_bams, parameters *params, char *chromosome_name, int chroSize)
{
	int count;
	mReads = (mei_Reads **) getMem( ( chroSize + 1) * sizeof( mei_Reads *));
	if( mReads == NULL)
		vh_logWarning ("Memory Problem in vh_createMaxClusterMEI");

	for( count = 0; count <= chroSize; count++)
		mReads[count] = NULL;

	vh_addToGenomeIndex_MEI (in_bams, params, chromosome_name);
}
