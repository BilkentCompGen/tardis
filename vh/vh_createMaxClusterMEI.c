#include "vh_createMaxClusterMEI.h"
#include "vh_common.h"
#include "vh_intervalhandler.h"
#include "vh_divethandler.h"
#include "vh_hash.h"
#include "vh_setcover.h"

mei_Reads **mReads;

Heap2 H_F; // The forward mapping heap from left of breakpoint
Heap2 H_R; // the reverse mapping heap from right of breakpoint
Heap2 H_S; // The soft clipping read

int minValue_heap2(Heap2 *heapName)
{
	if (heapName->heapSize > 0)
		return (heapName->heapArray[0].priorityValue);
	else
		return -1000;
}

void copyHeapEl2(HeapEl2 *dest, HeapEl2 *src)
{
	int len;

	dest->startPos = src->startPos;
	dest->stopPos = src->stopPos;
	strcpy(dest->cloneName, src->cloneName);
	dest->seqSimProb = src->seqSimProb;
	dest->orient1 = src->orient1;
	dest->orient2 = src->orient2;
	dest->indexOfMobileEl = src->indexOfMobileEl;
	dest->breakPointStart = src->breakPointStart;
	dest->mQual1=src->mQual1;
	dest->mQual2=src->mQual2;
	dest->editDist=src->editDist;
	dest->clonHitId=src->clonHitId;
	dest->libId=src->libId;
	strcpy(dest->libName, src->libName);
	strcpy(dest->indName, src->indName);

	len = strlen( src->mei_subclass);
	dest->mei_subclass = (char *) getMem ( sizeof (char) * (len + 1));
	strcpy(dest->mei_subclass, src->mei_subclass);

	dest->breakPointEnd = src->breakPointEnd;
	dest->priorityValue = src->priorityValue;
}

void heapBubleDown2(Heap2 *heapName)
{
	int heapIndex = 0;
	HeapEl2 tempEl;
	while((heapIndex * 2 + 1<heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 1].priorityValue)||(heapIndex * 2 + 2 < heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue))
	{
		if (heapIndex * 2 + 2<heapName->heapSize)
		{
			if (heapName->heapArray[heapIndex * 2 + 1].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue)
			{
				copyHeapEl2(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapEl2(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 2]));
				copyHeapEl2(& (heapName->heapArray[heapIndex * 2 + 2]), &tempEl);
				heapIndex = heapIndex * 2 + 2;
			}
			else
			{
				copyHeapEl2(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapEl2(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
				copyHeapEl2(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
				heapIndex = heapIndex * 2 + 1;
			}
		}
		else
		{
			copyHeapEl2(&tempEl, & (heapName->heapArray[heapIndex]));
			copyHeapEl2(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
			copyHeapEl2(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
			heapIndex = heapIndex * 2 + 1;
		}
	}
}

void heap_remove_top2(Heap2 *heapName)
{
	if (heapName->heapSize > 0)
	{
		copyHeapEl2(&(heapName->heapArray[0]), &(heapName->heapArray[heapName->heapSize - 1]));
		heapBubleDown2(heapName);
		heapName->heapSize--;
	}
}

void addToGenomeIndex_MEI (bam_info** in_bam, ref_genome* ref, parameters *params, char *chroName)
{
	mei_Reads *tmpMEI_Reads;
	discordantReadMappingMEI_Info *discordantReadPtr;
	softClip *softClipPtr;
	int numSample, count;
	int notRepeat;
	int countRR = 0, countFF = 0, countTotal = 0;
	char* MEI_Type;
	int libId = 0, len;

	for (numSample = 0; numSample < params->num_bams; numSample++)
	{
		for (count = 0; count < in_bam[numSample]->num_libraries; count++)
		{
			discordantReadPtr = in_bam[numSample]->libraries[count]->listMEI_Mapping;
			softClipPtr = in_bam[numSample]->libraries[count]->listSoftClip;

			while(discordantReadPtr != NULL)
			{
				notRepeat = notInRepeat(discordantReadPtr->chroName, discordantReadPtr->pos);
				if (notRepeat==1 && strcmp(discordantReadPtr->chroName, chroName) == 0)
				{
					if (discordantReadPtr->MEI_Type > -1)
					{
						tmpMEI_Reads = (mei_Reads*) getMem( sizeof( mei_Reads));

						len = strlen( discordantReadPtr->readName);
						tmpMEI_Reads->readName = (char *) getMem ( sizeof (char) * (len + 1));
						strcpy(tmpMEI_Reads->readName, discordantReadPtr->readName);

						len = strlen( discordantReadPtr->MEI_subclass);
						tmpMEI_Reads->MEI_Subclass = (char *) getMem ( sizeof (char) * (len + 1));
						set_str( &tmpMEI_Reads->MEI_Subclass, discordantReadPtr->MEI_subclass);

						tmpMEI_Reads->pos = discordantReadPtr->pos;
						tmpMEI_Reads->pos_End = discordantReadPtr->pos_End;
						tmpMEI_Reads->orient = discordantReadPtr->orient;
						tmpMEI_Reads->mQual = discordantReadPtr->qual;
						tmpMEI_Reads->MEI_Type = discordantReadPtr->MEI_Type;
						tmpMEI_Reads->readTypeSupport = 0; // The support is type paired-end read mapping
						tmpMEI_Reads->libId = libId;

						len = strlen(in_bam[numSample]->libraries[count]->libname);
						tmpMEI_Reads->libName = (char *) getMem (sizeof (char) * (len + 1));
						strcpy(tmpMEI_Reads->libName, in_bam[numSample]->libraries[count]->libname); //setstr

						len = strlen(in_bam[numSample]->sample_name);
						tmpMEI_Reads->indName = (char *) getMem (sizeof (char) * (len + 1));
						strcpy(tmpMEI_Reads->indName, in_bam[numSample]->sample_name);

						if (discordantReadPtr->orient=='F')
						{
							countFF++;
							tmpMEI_Reads->next = mReads[discordantReadPtr->pos_End];
							mReads[discordantReadPtr->pos_End] = tmpMEI_Reads;
						}
						else if (discordantReadPtr->orient == 'R')
						{
							countRR++;
							tmpMEI_Reads->next = mReads[max(0, discordantReadPtr->pos- in_bam[numSample]->libraries[count]->conc_max)];
							mReads[max(0, discordantReadPtr->pos - in_bam[numSample]->libraries[count]->conc_max)] = tmpMEI_Reads;
						}
					}
				}
				discordantReadPtr = discordantReadPtr->next;
			}
			while( softClipPtr != NULL)
			{
				notRepeat = notInRepeat(softClipPtr->chroName, softClipPtr->pos);
				if( notRepeat == 1 && strcmp( softClipPtr->chroName, chroName) == 0)
				{
					posMapSoftClip *ptrPosMapSoftClip = softClipPtr->ptrPosMapSoftClip;
					if( ptrPosMapSoftClip != NULL)
					{
						MEI_Type = meiIntervalSearch2( ref, chroName, ptrPosMapSoftClip->posMap);
						if( MEI_Type != NULL)
						{
							tmpMEI_Reads = (mei_Reads*) getMem(sizeof(mei_Reads));

							len = strlen(softClipPtr->readName);
							tmpMEI_Reads->readName = (char *) getMem (sizeof (char) * len + 1);
							strcpy(tmpMEI_Reads->readName, softClipPtr->readName);

							len = strlen( MEI_Type);
							tmpMEI_Reads->MEI_Subclass = (char *) getMem ( sizeof (char) * (len + 1));
							set_str( &tmpMEI_Reads->MEI_Subclass, MEI_Type);

							tmpMEI_Reads->pos = softClipPtr->pos;
							tmpMEI_Reads->orient = softClipPtr->orient;
							tmpMEI_Reads->mQual = softClipPtr->qual;
							tmpMEI_Reads->libId = libId;

							len = strlen(in_bam[numSample]->libraries[count]->libname);
							tmpMEI_Reads->libName = (char *) getMem (sizeof (char) * len + 1);
							strcpy(tmpMEI_Reads->libName, in_bam[numSample]->libraries[count]->libname);

							len = strlen(in_bam[numSample]->sample_name);
							tmpMEI_Reads->indName = (char *) getMem (sizeof (char) * len + 1);
							strcpy(tmpMEI_Reads->indName, in_bam[numSample]->sample_name);

							if( softClipPtr->op[0] == 4) //soft clipped at the beginning of the read
							{
								tmpMEI_Reads->readTypeSupport = 1;
								tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos-lenSplitReadBrakWindow];
								mReads[tmpMEI_Reads->pos - lenSplitReadBrakWindow] = tmpMEI_Reads;
							}
							else if( softClipPtr->op[softClipPtr->opCount- 1] == 4) // soft cliped at the end of the read
							{
								//printf("Soft Clipping at end\n");
								tmpMEI_Reads->readTypeSupport=2;
								tmpMEI_Reads->next = mReads[tmpMEI_Reads->pos + softClipPtr->opl[0] - lenSplitReadBrakWindow];
								mReads[tmpMEI_Reads->pos + softClipPtr->opl[0] - lenSplitReadBrakWindow] = tmpMEI_Reads;
							}
						}
					}
				}
				softClipPtr=softClipPtr->next;
			}
			libId++;
		}
	}
}


void vh_addToGenomeIndex_MEI (bam_info** in_bams, char *chroName)
{
	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	mei_Reads *tmpMEI_Reads;
	int countRR=0, countFF=0, countTotal=0, libId=0, len;

	libInfo = g_libInfo;
	while (libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while (divetReadMappingPtr != NULL)
		{
			if (strcmp (divetReadMappingPtr->chroName, chroName) == 0 && ((divetReadMappingPtr->svType == 'M' || divetReadMappingPtr->svType == 'X')))
			{
				tmpMEI_Reads=(mei_Reads *)getMem(sizeof(mei_Reads));
				len=strlen(divetReadMappingPtr->readName->readName);
				tmpMEI_Reads->readName = (char *) getMem (sizeof (char) * len + 1);
				strcpy(tmpMEI_Reads->readName, divetReadMappingPtr->readName->readName);

				tmpMEI_Reads->mQual=divetReadMappingPtr->phredScore;
				tmpMEI_Reads->readTypeSupport=0; // The support is type paired-end read mapping

				len=strlen(libInfo->libName);
				tmpMEI_Reads->libName=(char *) getMem (sizeof (char) * len + 1);
				strcpy(tmpMEI_Reads->libName, libInfo->libName);

				len=strlen(divetReadMappingPtr->libInfo->indName);
				tmpMEI_Reads->indName=(char *) getMem (sizeof (char) * len + 1);
				strcpy(tmpMEI_Reads->indName,libInfo->indName);

				/* Modify to get subclass name instead */
				len = strlen( divetReadMappingPtr->mei_subclass);
				tmpMEI_Reads->MEI_Subclass = (char *) getMem ( sizeof (char) * (len + 1));
				set_str( &tmpMEI_Reads->MEI_Subclass, divetReadMappingPtr->mei_subclass);

				tmpMEI_Reads->libId=libId;

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

void push_heap2(Heap2 *heapName, HeapEl2 *newEl)
{
	int heapIndex = heapName->heapSize;
	HeapEl2 tempEl;
	heapName->heapSize++;
	copyHeapEl2(&(heapName->heapArray[heapIndex]), newEl);
	heapIndex = heapName->heapSize - 1;

	while(heapIndex > 0 && heapName->heapArray[heapIndex].priorityValue < heapName->heapArray[(heapIndex + 1) / 2 -1].priorityValue)
	{
		copyHeapEl2(&tempEl, &heapName->heapArray[heapIndex]);
		copyHeapEl2(&(heapName->heapArray[heapIndex]),  &(heapName->heapArray[(heapIndex + 1) / 2 -1]));
		copyHeapEl2(&(heapName->heapArray[(heapIndex + 1) / 2 -1]), &tempEl);
		heapIndex = (heapIndex + 1) / 2 - 1;
	}
}

void add_R_Heap(int pos)
{
	int len;
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapEl2 *newEl;
	while(mei_ReadsPtr!=NULL)
	{
		if (mei_ReadsPtr->orient=='R' && mei_ReadsPtr->readTypeSupport==0)
		{
			newEl=(HeapEl2 *)getMem(sizeof(HeapEl2));
			strcpy(newEl->cloneName, mei_ReadsPtr->readName);
			newEl->startPos=mei_ReadsPtr->pos;
			newEl->orient1=mei_ReadsPtr->orient;
			newEl->mQual1=mei_ReadsPtr->mQual;
			newEl->indexOfMobileEl=mei_ReadsPtr->MEI_Type;
			newEl->libId=mei_ReadsPtr->libId;
			strcpy(newEl->libName, mei_ReadsPtr->libName);
			strcpy(newEl->indName, mei_ReadsPtr->indName);

			len = strlen( mei_ReadsPtr->MEI_Subclass);
			newEl->mei_subclass = (char *) getMem ( sizeof (char) * (len + 1));
			set_str( &newEl->mei_subclass, mei_ReadsPtr->MEI_Subclass);

			newEl->priorityValue = mei_ReadsPtr->pos;
			push_heap2(&H_R, newEl);
			free(newEl);
		}
		mei_ReadsPtr=mei_ReadsPtr->next;
	}
}

void add_S_Heap(int pos)
{
	int len;
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapEl2 *newEl;
	while(mei_ReadsPtr!=NULL)
	{
		if ( mei_ReadsPtr->readTypeSupport==1 || mei_ReadsPtr->readTypeSupport==2)
		{
			newEl=(HeapEl2 *)getMem(sizeof(HeapEl2));
			strcpy(newEl->cloneName, mei_ReadsPtr->readName);
			newEl->startPos=mei_ReadsPtr->pos;
			newEl->orient1=mei_ReadsPtr->orient;
			newEl->mQual1=mei_ReadsPtr->mQual;
			newEl->indexOfMobileEl=mei_ReadsPtr->MEI_Type;
			newEl->libId=mei_ReadsPtr->libId;
			strcpy(newEl->libName, mei_ReadsPtr->libName);
			strcpy(newEl->indName, mei_ReadsPtr->indName);

			len = strlen( mei_ReadsPtr->MEI_Subclass);
			newEl->mei_subclass = (char *) getMem ( sizeof (char) * (len + 1));
			set_str( &newEl->mei_subclass, mei_ReadsPtr->MEI_Subclass);

			newEl->priorityValue = pos+2*lenSplitReadBrakWindow;
			push_heap2(& H_S, newEl);
			free(newEl);
		}
		mei_ReadsPtr=mei_ReadsPtr->next;
	}
}

void add_F_Heap(bam_info** in_bams, int pos)
{
	int len;
	LibraryInfo *libInfo;
	mei_Reads *mei_ReadsPtr;
	mei_ReadsPtr = mReads[pos];
	HeapEl2 *newEl;
	int max_delta_val;

	while(mei_ReadsPtr!=NULL)
	{
		if (mei_ReadsPtr->orient=='F' && mei_ReadsPtr->readTypeSupport==0)
		{
			newEl=(HeapEl2 *)getMem(sizeof(HeapEl2));
			strcpy(newEl->cloneName, mei_ReadsPtr->readName);
			newEl->orient1=mei_ReadsPtr->orient;
			newEl->mQual1=mei_ReadsPtr->mQual;
			newEl->indexOfMobileEl=mei_ReadsPtr->MEI_Type;
			newEl->libId=mei_ReadsPtr->libId;
			newEl->startPos=mei_ReadsPtr->pos_End;
			strcpy(newEl->libName, mei_ReadsPtr->libName);
			strcpy(newEl->indName, mei_ReadsPtr->indName);

			len = strlen( mei_ReadsPtr->MEI_Subclass);
			newEl->mei_subclass = (char *) getMem ( sizeof (char) * (len + 1));
			set_str( &newEl->mei_subclass, mei_ReadsPtr->MEI_Subclass);

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

			//fprintf(stderr, " %s delta=%d - priority=%d\n",mei_ReadsPtr->libName, max_delta_val, newEl->priorityValue);
			push_heap2(& H_F, newEl);
			free(newEl);
		}
		mei_ReadsPtr=mei_ReadsPtr->next;
	}
}

void outputClusters(char* chroName)
{
	int MEIType;
	int MEIType2;
	int R_count=0, F_count=0;
	int count;
	int orientMEI;

	if( H_F.heapSize > 0  && H_R.heapSize > 0)
	{
		for( MEIType = 0; MEIType < 6; MEIType++)
		{
			R_count=0;
			F_count=0;

			/* match reads with the same type of MEI (Alu, etc) and the same orientation of mapping
			 *  so if it is H_R assuming the read is in R orientation then mapping
			 * inside a Alu with R orientation then the insertion is type F*/

			if ((MEIType%2)==0)
				MEIType2=MEIType+1;
			else
				MEIType2=MEIType-1;
			if (MEIType%2==0)
				orientMEI=10;
			else
				orientMEI=12;

			for (count=0; count<H_F.heapSize; count++)
			{
				if (H_F.heapArray[count].indexOfMobileEl==MEIType)
					F_count++;
			}
			for (count=0; count<H_R.heapSize; count++)
			{
				if (H_R.heapArray[count].indexOfMobileEl==MEIType2)
					R_count++;
			}

			if (F_count>2 && R_count>2)
			{
				for (count=0; count<H_F.heapSize; count++)
				{
					if (H_F.heapArray[count].indexOfMobileEl==MEIType && MEIType>-1 &&  MEIType<2)
					{
						fprintf(fileOutput, "%s %s %i %s %i %i %i %i %s %s F F %i %i ",
								H_F.heapArray[count].cloneName,
								chroName,
								H_F.heapArray[count].startPos,
								H_F.heapArray[count].mei_subclass,
								H_F.heapArray[count].startPos,
								orientMEI,
								H_F.heapArray[count].mQual1,
								H_F.heapArray[count].mQual1,
								H_F.heapArray[count].libName,
								H_F.heapArray[count].indName,
								H_F.heapArray[count].mQual1,
								H_F.heapArray[count].mQual1);
					}
					else if (H_F.heapArray[count].indexOfMobileEl==MEIType && MEIType>1 &&  MEIType<4)
					{
						fprintf(fileOutput, "%s %s %i %s %i %i %i %i %s %s F F %i %i ",
								H_F.heapArray[count].cloneName,
								chroName,
								H_F.heapArray[count].startPos,
								H_F.heapArray[count].mei_subclass,
								H_F.heapArray[count].startPos,
								orientMEI,
								H_F.heapArray[count].mQual1,
								H_F.heapArray[count].mQual1,
								H_F.heapArray[count].libName,
								H_F.heapArray[count].indName,
								H_F.heapArray[count].mQual1,
								H_F.heapArray[count].mQual1);
					}
				}
				for (count=0; count<H_R.heapSize; count++)
				{
					if (H_R.heapArray[count].indexOfMobileEl==MEIType2 && MEIType>-1 &&  MEIType<2)
					{
						fprintf(fileOutput, "%s %s %i %s %i %i %i %i %s %s R R %i %i ",
								H_R.heapArray[count].cloneName,
								chroName,
								H_R.heapArray[count].startPos,
								H_R.heapArray[count].mei_subclass,
								H_R.heapArray[count].startPos,
								orientMEI,
								H_R.heapArray[count].mQual1,
								H_R.heapArray[count].mQual1,
								H_R.heapArray[count].libName,
								H_R.heapArray[count].indName,
								H_R.heapArray[count].mQual1,
								H_R.heapArray[count].mQual1);
					}
					else if (H_R.heapArray[count].indexOfMobileEl==MEIType2 && MEIType>1 &&  MEIType<4)
					{
						fprintf(fileOutput, "%s %s %i %s %i %i %i %i %s %s R R %i %i ",
								H_R.heapArray[count].cloneName,
								chroName,
								H_R.heapArray[count].startPos,
								H_R.heapArray[count].mei_subclass,
								H_R.heapArray[count].startPos,
								orientMEI,
								H_R.heapArray[count].mQual1,
								H_R.heapArray[count].mQual1,
								H_R.heapArray[count].libName,
								H_R.heapArray[count].indName,
								H_R.heapArray[count].mQual1,
								H_R.heapArray[count].mQual1);
					}
				}
				for (count=0; count<H_S.heapSize; count++)
				{
					if (H_S.heapArray[count].indexOfMobileEl==MEIType2 && MEIType>-1 &&  MEIType<2)
						fprintf(fileOutput, "%s %s %i %s %i %i %i %i SplitRead %s S S %i %i ",
								H_S.heapArray[count].cloneName,
								chroName,
								H_S.heapArray[count].startPos,
								H_S.heapArray[count].mei_subclass,
								H_S.heapArray[count].startPos,
								orientMEI,
								H_S.heapArray[count].mQual1,
								H_S.heapArray[count].mQual1,
								H_S.heapArray[count].indName,
								H_S.heapArray[count].mQual1,
								H_S.heapArray[count].mQual1);
					else if (H_S.heapArray[count].indexOfMobileEl==MEIType2 && MEIType>1 &&  MEIType<4)
						fprintf(fileOutput, "%s %s %i %s %i %i %i %i SplitRead %s S S %i %i ",
								H_S.heapArray[count].cloneName,
								chroName,
								H_S.heapArray[count].startPos,
								H_S.heapArray[count].mei_subclass,
								H_S.heapArray[count].startPos,
								orientMEI,
								H_S.heapArray[count].mQual1,
								H_S.heapArray[count].mQual1,
								H_S.heapArray[count].indName,
								H_S.heapArray[count].mQual1,
								H_S.heapArray[count].mQual1);
				}
				fprintf(fileOutput, "END\n");
				cluster_count++;
			}
		}
	}
}

void MEICluster_Region(bam_info** in_bams, char* chroName, int chroSize)
{
	int indexWindowStart=0, indexWindowEnd=0, i=0;
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
		if (((H_R.heapSize>0 && minValue_heap2(&H_R) == brkPointIndex) || (minValue_heap2(&H_F) == brkPointIndex && H_F.heapSize>0) || (minValue_heap2(&H_S)==brkPointIndex && H_S.heapSize>0)) && boolMEITypeNewAdded)
		{
			if ((H_R.heapSize+H_F.heapSize + H_S.heapSize)>0)
			{
				outputClusters(chroName);
			}
			boolMEITypeNewAdded=0;
		}
		while(H_R.heapSize>0 && minValue_heap2(&H_R)==brkPointIndex)
		{
			heap_remove_top2(&H_R);
		}
		while(H_F.heapSize>0 && minValue_heap2(&H_F)==brkPointIndex)
		{
			heap_remove_top2(&H_F);
		}
		while(H_S.heapSize>0 && minValue_heap2(&H_S)==brkPointIndex)
		{
			heap_remove_top2(&H_S);
		}
	}
}

void freeHeap2(Heap2 *heapName)
{
	heapName->heapSize=0;
}

void vh_finalizeReadMapping_Mei( int chroSize)
{
	int count;
	mei_Reads *ptr, *ptrNext;
	for( count = 0; count < chroSize; count++)
	{
		ptr = mReads[count];
		while( ptr!=NULL)
		{
			ptrNext = ptr->next;
			free( ptr->indName);
			free( ptr->libName);
			free( ptr->readName);
			free( ptr->MEI_Subclass);
			free( ptr);
			ptr = ptrNext;
		}
	}
	if( mReads != NULL)
	{
		free( mReads);
	}
	freeHeap2( &H_S);
	freeHeap2( &H_R);
	freeHeap2( &H_F);
}


void initializeReadMapping_MEI( bam_info** in_bams, ref_genome* ref, parameters *params, char *chroName, int chroSize)
{
	int count;
	mReads = (mei_Reads **) getMem( chroSize * sizeof( mei_Reads *));
	for( count = 0; count < chroSize; count++)
	{
		mReads[count] = NULL;
	}
	addToGenomeIndex_MEI( in_bams, ref, params, chroName);
}


void vh_initializeReadMapping_MEI (bam_info** in_bams, parameters *params, char *chroName, int chroSize)
{
	int count;
	mReads = (mei_Reads **) getMem(chroSize * sizeof(mei_Reads *));
	for (count=0; count<chroSize; count++)
	{
		mReads[count]=NULL;
	}
	vh_addToGenomeIndex_MEI (in_bams, chroName);
}
