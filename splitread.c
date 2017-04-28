/*
 * splitread.c


 *
 *  Created on: Feb 2, 2017
 *      Author: tardis
 */

#include "splitread.h"
#include "bamonly.h"
#include <ctype.h>
#include <math.h>
#include "vh/vh_common.h"
#include <htslib/faidx.h>

lociInRef *hash_table_SR[SR_HASH_SIZE];
//lociInRef **hash_table_SR;
char *ref_seq_per_chr;

/* Read the reference genome */
void readReferenceSeq( ref_genome *ref, parameters *params, int chr_index)
{
	int i = 0, min, max, loc_length;
	char *ref_seq;
	long bp_cnt = 0;
	faidx_t* ref_fai;

	min = 0, max = 999;
	ref_fai = fai_load( ref->ref_name);
	ref_seq_per_chr = ( char *) getMem( (ref->chrom_lengths[chr_index]+1) * sizeof( char));

	while ( max < ref->chrom_lengths[chr_index])
	{
		ref_seq = faidx_fetch_seq( ref_fai, ref->chrom_names[chr_index], min, max, &loc_length);

		for( i = 0; i < loc_length; i++)
		{
			if( bp_cnt < ref->chrom_lengths[chr_index])
				ref_seq_per_chr[bp_cnt] = toupper(ref_seq[i]);
			bp_cnt++;
		}
		if( bp_cnt >= ref->chrom_lengths[chr_index])
			break;

		min += loc_length;
		max += loc_length;
	}
	fai_destroy(ref_fai);
	free(ref_seq);

	ref_seq_per_chr[bp_cnt] = '\0';
	create_10bp_HashIndex( ref, chr_index);
}

int hash_function_ref( ref_genome *ref, char *str)
{
	int i, count=0;

	for( i = 0; i < strlen( str); i++)
		count = (2 * count) + (int)str[i];

	return count % 1000000;
}

void create_10bp_HashIndex( ref_genome* ref, int chr_index)
{
	int i, index;
	char str[11];
	lociInRef *ptr;

	//hash_table_SR = ( lociInRef **) getMem( SR_HASH_SIZE * sizeof( lociInRef*));

	for( i = 0; i < SR_HASH_SIZE; i++)
		hash_table_SR[i] = NULL;

	for( i = 0; i < ref->chrom_lengths[chr_index] - 10; i++)
	{
		strncpy( str, &(ref_seq_per_chr[i]), 10);
		str[10] = '\0';

		if( strcmp( str, "NNNNNNNNNN\0") != 0)
		{
			index = hash_function_ref( ref, str);
			ptr = ( lociInRef *) getMem( sizeof( lociInRef));
			ptr->pos = i;
			ptr->next = hash_table_SR[index];
			hash_table_SR[index] = ptr;
		}
	}
}

void free_10bp_HashIndex()
{
	int i;
	lociInRef *ptr, *ptrTemp;

	for( i = 0; i < SR_HASH_SIZE; i++)
	{
		ptr = hash_table_SR[i];
		while( ptr != NULL)
		{
			ptrTemp = ptr->next;
			free( ptr);
			ptr = ptrTemp;
		}
	}
	//if( hash_table_SR != NULL)
		//free( hash_table_SR);

	free( ref_seq_per_chr);
}

posMapSoftClip *almostPerfect_match_seq_ref(ref_genome* ref, int chr_index, char *str, int pos)
{
	int i, index, posMapSize, posMap[10000], hammingDisMap[10000];
	char orient[10000];// orient of the mapping
	int dist, reverseMatch;
	lociInRef *ptr;
	char seed[11];
	char *strRev;
	posMapSoftClip *tmpSoftClipMap, *returnPtr;
	returnPtr = NULL;

	strncpy( seed, str, 10);
	seed[10] = '\0';
	index = hash_function_ref( ref, seed);
	ptr = hash_table_SR[index];

	posMapSize = 0;

	while( ptr != NULL)
	{
		if( abs( ptr->pos - pos) < 100000)
		{
			dist = hammingDistance( &( ref_seq_per_chr[ptr->pos]), str, strlen( str));
			if( dist <= ( 0.05 * strlen( str)))
			{
				posMap[posMapSize] = ptr->pos;
				orient[posMapSize] = 'F';
				hammingDisMap[posMapSize] = dist;
				posMapSize = posMapSize + 1;
			}
		}
		ptr = ptr->next;
		if( posMapSize > 10)
			ptr = NULL;
	}

	if( posMapSize < 11)
	{
		strRev = ( char *)getMem( ( strlen( str) + 1) * sizeof( char));
		for( i = 0; i < strlen( str); i++)
		{
			if( str[i] == 'A')
				strRev[strlen( str) - i - 1] = 'T';
			if( str[i] == 'T')
				strRev[strlen( str) - i - 1] = 'A';
			if( str[i] == 'G')
				strRev[strlen( str) - i - 1] = 'C';
			if( str[i] == 'C')
				strRev[strlen( str) - i - 1] = 'G';
			if( str[i] == 'N')
				strRev[strlen( str) - i - 1] = 'N';
		}
		strRev[strlen( str)] = '\0';
		strncpy( seed, strRev, 10);
		seed[10] = '\0';
		index = hash_function_ref( ref, seed);
		ptr = hash_table_SR[index];

		reverseMatch = 0;
		while( ptr != NULL)
		{
			if( abs( ptr->pos - pos) < 100000)
			{
				dist = hammingDistance( &( ref_seq_per_chr[ptr->pos]), strRev, strlen( strRev));
				if( dist <= ( 0.05 * strlen( strRev)))
				{
					posMap[posMapSize] = ptr->pos;
					orient[posMapSize] = 'R';
					hammingDisMap[posMapSize] = dist;
					posMapSize = posMapSize + 1;
					reverseMatch = 1;
				}
			}
			ptr = ptr->next;
			if( posMapSize > 10)
				ptr = NULL;
		}

		free( strRev);
	}

	if( posMapSize < 12 && posMapSize > 0)
	{
		for( i = 0; i < posMapSize; i++)
		{
			tmpSoftClipMap = ( posMapSoftClip *) getMem( sizeof( posMapSoftClip));
			tmpSoftClipMap->posMap = posMap[i];
			tmpSoftClipMap->orient = orient[i];

			if( posMapSize < 11)
				tmpSoftClipMap->mapq = 60 / posMapSize;
			else
				posMapSize = 0;
			tmpSoftClipMap->next = returnPtr;
			returnPtr = tmpSoftClipMap;
		}
	}
	return returnPtr;
}


void countNumSoftClipInCluster( parameters *params, ref_genome* ref, bam_info* in_bam, int chr_index)
{
	int i, j, countNumSR, posStartSF;
	int* countSoftClip = ( int *)getMem( ( ref->chrom_lengths[chr_index] + 1000) * sizeof( int));
	softClip *ptrSoftClip;

	for( i = 0; i < ref->chrom_lengths[chr_index] + 1000; i++)
		countSoftClip[i] = 0;

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		ptrSoftClip = in_bam->libraries[i]->listSoftClip;
		while( ptrSoftClip != NULL)
		{
			if( ptrSoftClip->qual > params->mq_threshold && ptrSoftClip->avgPhredQualSoftClip > params->mq_threshold)
			{
				if( ptrSoftClip->op[0] == BAM_CSOFT_CLIP)
					countSoftClip[max( 0, ptrSoftClip->pos + ptrSoftClip->opl[0])]++;
				else if( ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
					countSoftClip[max( 0, ptrSoftClip->pos - ptrSoftClip->opl[ptrSoftClip->opCount - 1] + in_bam->libraries[i]->read_length)]++;
			}
			ptrSoftClip = ptrSoftClip->next;
		}
	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		ptrSoftClip = in_bam->libraries[i]->listSoftClip;
		while( ptrSoftClip != NULL)
		{
			countNumSR = 0;
			if( ptrSoftClip->op[0] == BAM_CSOFT_CLIP)
			{
				posStartSF = ptrSoftClip->pos + ptrSoftClip->opl[0];
				for( j = max( 0, posStartSF - MIN_SOFTCLIP_LEN); j < min( ref->chrom_lengths[chr_index], posStartSF + MIN_SOFTCLIP_LEN); j++)
					countNumSR = countNumSR + countSoftClip[j];
			}
			else if( ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
			{
				posStartSF = ptrSoftClip->pos - ptrSoftClip->opl[ptrSoftClip->opCount-1] + in_bam->libraries[i]->read_length;

				for(j = max(0, posStartSF - MIN_SOFTCLIP_LEN); j < min(ref->chrom_lengths[chr_index], posStartSF + MIN_SOFTCLIP_LEN); j++)
					countNumSR = countNumSR + countSoftClip[j];
			}

			ptrSoftClip->numSoftClipInConcordance = countNumSR;
			ptrSoftClip = ptrSoftClip->next;
		}
	}
	if( countSoftClip != NULL)
		free( countSoftClip);
}


void mapSoftClipToRef( bam_info* in_bam,  parameters* params, ref_genome* ref, int chr_index)
{
	int i;
	int mappingPos[100000];
	int numMap = 0;
	char *str;
	softClip *ptrSoftClip;

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		str = ( char *)malloc( sizeof( char) * in_bam->libraries[i]->read_length);

		ptrSoftClip = in_bam->libraries[i]->listSoftClip;
		while( ptrSoftClip != NULL)
		{
			if( ptrSoftClip->avgPhredQualSoftClip > params->mq_threshold && ptrSoftClip->numSoftClipInConcordance > 1
					&& ptrSoftClip->numSoftClipInConcordance < 100)
			{
				if( ptrSoftClip->op[0] == BAM_CSOFT_CLIP)
				{
					strncpy( str, ptrSoftClip->softClipString, ptrSoftClip->opl[0]);
					str[ptrSoftClip->opl[0]] = '\0';
				}
				else if( ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
				{
					strncpy( str, &( ptrSoftClip->softClipString[in_bam->libraries[i]->read_length - ptrSoftClip->opl[ptrSoftClip->opCount - 1]]), in_bam->libraries[i]->read_length - ptrSoftClip->opl[ptrSoftClip->opCount - 1]);
					str[in_bam->libraries[i]->read_length - ptrSoftClip->opl[ptrSoftClip->opCount - 1]] = '\0';
				}
				ptrSoftClip->ptrPosMapSoftClip = almostPerfect_match_seq_ref( ref, chr_index, str, ptrSoftClip->pos);
			}
			ptrSoftClip = ptrSoftClip->next;
		}
		if( str != NULL)
			free( str);
	}
}

void addSoftClip( ref_genome* ref, bam_info * in_bam, bam1_t* bam_alignment, int library_index, int flag, int *op, int *opl, int chrID)
{
	int len, i;
	float avgPhredQual = 0;
	bam1_core_t bam_alignment_core;
	bam_alignment_core = bam_alignment->core;

	softClip *newEl = ( softClip *) malloc( sizeof( softClip));
	newEl->readName = ( char*)malloc( ( strlen( bam_alignment->data) + 7) * sizeof( char));
	strncpy(newEl->readName, bam_alignment->data, strlen( bam_alignment->data));

	newEl->readName[strlen( bam_alignment->data)] = '_';
	newEl->readName[strlen( bam_alignment->data) + 1] = 'S';
	newEl->readName[strlen( bam_alignment->data) + 2] = 'p';
	newEl->readName[strlen( bam_alignment->data) + 3] = 'l';
	newEl->readName[strlen( bam_alignment->data) + 4] = 'i';
	newEl->readName[strlen( bam_alignment->data) + 5] = 't';
	newEl->readName[strlen( bam_alignment->data) + 6] = '\0';

	/* Get the name of the chromosome */
	len = strlen( ref->chrom_names[chrID]);
	newEl->chroName = ( char *) getMem( sizeof( char) * len + 1);
	strcpy( newEl->chroName, ref->chrom_names[chrID]);

	newEl->pos = bam_alignment_core.pos;
	newEl->qual = bam_alignment_core.qual;
	newEl->softClipString = ( char *)malloc( ( bam_alignment_core.l_qseq + 1) * sizeof( char));

	if( ( flag & BAM_FREVERSE) != 0)
		newEl->orient='R';
	else
		newEl->orient='F';

	uint8_t *a_qual = bam_get_qual(bam_alignment);
	for( i = 0; i < bam_alignment->core.n_cigar; i++)
	{
		newEl->op[i] = op[i];
		newEl->opl[i] = opl[i];
	}
	newEl->opCount = bam_alignment->core.n_cigar;

	for( i = 0; i < bam_alignment_core.l_qseq; i++)
	{
		if( bam_seqi( bam_get_seq( bam_alignment), i) == 1)
			newEl->softClipString[i] = 'A';
		else if( bam_seqi( bam_get_seq( bam_alignment), i) == 2)
			newEl->softClipString[i] = 'C';
		else if(bam_seqi( bam_get_seq( bam_alignment), i) == 4)
			newEl->softClipString[i] = 'G';
		else if( bam_seqi( bam_get_seq( bam_alignment), i) == 8)
			newEl->softClipString[i] = 'T';
		else if( bam_seqi( bam_get_seq( bam_alignment), i) == 15)
			newEl->softClipString[i] = 'N';
	}

	if( op[0] == BAM_CSOFT_CLIP && opl[0] > MIN_SOFTCLIP_LEN)
	{
		for( i = 0; i < opl[0]; i++)
			avgPhredQual = avgPhredQual + a_qual[i];

		avgPhredQual = ( float)avgPhredQual / ( float)opl[0];
	}
	else if( op[newEl->opCount - 1] == BAM_CSOFT_CLIP && opl[newEl->opCount - 1] > MIN_SOFTCLIP_LEN)
	{
		for( i = in_bam->libraries[library_index]->read_length - opl[newEl->opCount - 1] - 1; i < in_bam->libraries[library_index]->read_length; i++)
			avgPhredQual = avgPhredQual + a_qual[i];

		avgPhredQual = ( float)avgPhredQual / opl[newEl->opCount - 1];
	}

	newEl->softClipString[bam_alignment_core.l_qseq] = '\0';
	newEl->avgPhredQualSoftClip = ( int)floorf( avgPhredQual);
	newEl->ptrPosMapSoftClip = NULL;
	newEl->next = in_bam->libraries[library_index]->listSoftClip;
	in_bam->libraries[library_index]->listSoftClip = newEl;

	sr_cnt_bam++;
}
