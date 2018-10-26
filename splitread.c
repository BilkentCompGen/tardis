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
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#define HASH_COUNT 0
#define HASH_BUILD 1

#define SR_LOOKAHEAD 100000


int **hash_table_array;
//char *ref_seq_per_chr = NULL;
int *hash_table_count;
int *hash_table_iter;
//int hash_size;
float total_hash_build_time = 0;

/* Read the reference genome */
void readReferenceSeq( parameters *params, int chr_index)
{
	int i, min, max, loc_length;
	char *ref_seq;
	long bp_cnt = 0;
	faidx_t* ref_fai;

	if (params->ref_seq != NULL)
	  {
	    fprintf (stderr, ". Reference genome is already loaded.\n");
	    return;
	  }
	
	min = 0, max = 999;
	ref_fai = fai_load( params->ref_genome);

	params->ref_seq = ( char *) getMem( (params->this_sonic->chromosome_lengths[chr_index] + 1) * sizeof( char));

	while ( max < params->this_sonic->chromosome_lengths[chr_index])
	{
		ref_seq = faidx_fetch_seq( ref_fai, params->this_sonic->chromosome_names[chr_index], min, max, &loc_length);

		for( i = 0; i < loc_length; i++)
		{
		  /* can we do this faster with memcpy? */
			if( bp_cnt < params->this_sonic->chromosome_lengths[chr_index])
				params->ref_seq[bp_cnt] = toupper( ref_seq[i]);
			bp_cnt++;
		}
		if( bp_cnt >= params->this_sonic->chromosome_lengths[chr_index])
			break;

		min += loc_length;
		max += loc_length;
		free( ref_seq);
	}
	fai_destroy( ref_fai);

	params->ref_seq[bp_cnt] = '\0';

	build_hash_table(params->ref_seq, bp_cnt, params->hash_size, HASH_COUNT);

	create_hash_table(params, bp_cnt);
}

void init_hash_count(parameters *params){

        params->hash_size = pow (4, HASHKMERLEN);
	hash_table_count = (int *) getMem( params->hash_size * sizeof (int));
	hash_table_iter = (int *) getMem( params->hash_size * sizeof (int));
	params->ref_seq = NULL;
}

void init_hash_table(parameters *params){
  int i;
  int hash_size = params->hash_size;
  hash_table_array = (int **) getMem (hash_size * sizeof (int *));

  for (i=0; i<hash_size; i++){
    if ( hash_table_count[i] != 0 && hash_table_count[i] < MAX_SR_HIT)
      hash_table_array[i] = (int *) getMem (hash_table_count[i] * sizeof (int));
    else{
      hash_table_count[i] = 0;
      hash_table_array[i] = NULL;
    }
  }
}

void free_hash_table(parameters *params){
  int i;
  int hash_size = params->hash_size;
  
  for (i=0; i<hash_size; i++){
    if ( hash_table_count[i] != 0)
      freeMem (hash_table_array[i], (hash_table_count[i] * sizeof (int)));
  }
  freeMem (hash_table_array, ( hash_size * sizeof (int)));

  freeMem( params->ref_seq, strlen(params->ref_seq));
  params->ref_seq = NULL;
}

void build_hash_table(const char *ref, int len, int hash_size, int mode){
  int i = 0, j = 0;
  char seed[HASHKMERLEN + 1];
  int hash_val;
  int mask;
  int len_limit;
  
  mask = (hash_size - 1) >> 2;

  if (mode == HASH_COUNT)
    memset (hash_table_count, 0, hash_size * sizeof (int));
  else
    memset (hash_table_iter, 0, hash_size * sizeof (int));


  len_limit = len - HASHKMERLEN + 1;

  /* get first kmer */
  while (!is_dna_letter(ref[i])) i++;

  strncpy (seed, ref+i, HASHKMERLEN);
  seed[HASHKMERLEN] = 0;
  
  while (!is_kmer_valid (seed)){
    i++;
    strncpy (seed, ref+i, HASHKMERLEN);
    seed[HASHKMERLEN] = 0;
  } 
  
  hash_val = hash_function_ref (seed);
  if (mode == HASH_COUNT)
    ( hash_table_count[hash_val])++;
  else if (hash_table_count[hash_val] != 0) {
    hash_table_array[hash_val][hash_table_iter[hash_val]] = i;
    ( hash_table_iter[hash_val])++;
  }

  j = i + HASHKMERLEN;
  
  while (i < len_limit){
    i++; //shift 
    if (is_dna_letter(ref[j])){
      hash_val = hash_function_next( hash_val, mask, ref[j++]);
      if (mode == HASH_COUNT)
	( hash_table_count[hash_val])++;
      else if (hash_table_count[hash_val] != 0) {
	hash_table_array[hash_val][hash_table_iter[hash_val]] = i;
	( hash_table_iter[hash_val])++;
      }
    }
    else{
      /* recover from non-ACGT */
      while (!is_dna_letter(ref[j]) && i < len_limit) { i++; j++; }
      if (i >= len_limit)
	break;

      strncpy (seed, ref+i, HASHKMERLEN);
      seed[HASHKMERLEN] = 0;

      while (!is_kmer_valid (seed)){
	i++;
	strncpy (seed, ref+i, HASHKMERLEN);
	seed[HASHKMERLEN] = 0;
      } 

      j = i + HASHKMERLEN;
      hash_val = hash_function_ref (seed);

      if (mode == HASH_COUNT)
	( hash_table_count[hash_val])++;
      else if (hash_table_count[hash_val] != 0) {
	hash_table_array[hash_val][hash_table_iter[hash_val]] = i;
	( hash_table_iter[hash_val])++;
      }
    }
  }

}

unsigned int hash_function_ref( char *str)
{
	/* this strictly assumes HASHKMERLEN < 16 and is_kmer_valid is already called and retuned TRUE */

	int i = 0;
	unsigned int val = 0; unsigned int numericVal = 0;
	while(i < HASHKMERLEN)
	{
	        numericVal = (str[i++] & 0x6) >> 1;
		val = (val << 2) | numericVal;
	}
	return val;
}

unsigned int hash_function_next( unsigned int prev_hash, unsigned int mask, const char next_char)
{

	/* this strictly assumes HASHKMERLEN < 16 */
        return (((prev_hash & mask) << 2) |  ((next_char & 0x6) >> 1));

}

int is_kmer_valid (char *str){

        int i, l;
	l = strlen(str);

	if (l < HASHKMERLEN)
		return 0;

	for (i=0; i<HASHKMERLEN; i++)
	  if (str[i] != 'A' && str[i] != 'C' && str[i] != 'G' && str[i] != 'T')
	    {
	      return 0;
	    }
	
	return 1;
}

void create_hash_table( parameters *params, int len){

	
        init_hash_table( params);
	build_hash_table( params->ref_seq, len, params->hash_size, HASH_BUILD);

}

/* deprecated 
void create_HashIndex( parameters* params, int chr_index)
{
	int i,k, ind;
	char str[HASHKMERLEN + 1];
	lociInRef *ptr;

	hash_table_SR = ( lociInRef **) getMem( ( hash_size + 1) * sizeof( lociInRef*));


	for( i = 0; i < hash_size; i++)
		hash_table_SR[i] = NULL;

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index] - HASHKMERLEN; i++)
	{
		if( ref_seq_per_chr[i] == '\0')
			break;

		strncpy( str, &( ref_seq_per_chr[i]), HASHKMERLEN);
		str[HASHKMERLEN] = '\0';
		
		if( is_kmer_valid(str))
		{
			ind = hash_function_ref( str);
			/*
			if( ind < 0)
			continue;

			ptr = ( lociInRef *) getMem( sizeof( lociInRef));
			ptr->pos = i;
			ptr->next = hash_table_SR[ind];
			hash_table_SR[ind] = ptr;
		}
	}

}


void free_HashIndex(void)
{
	int i;
	lociInRef *ptr, *ptrNext;


	for( i = 0; i < hash_size; i++)
	{
		ptr = hash_table_SR[i];
		while( ptr != NULL)
		{
			ptrNext = ptr->next;
			free( ptr);
			ptr = ptrNext;
		}
	}

	if( hash_table_SR != NULL)
		free( hash_table_SR);

	freeMem( ref_seq_per_chr, strlen(ref_seq_per_chr));
}

*/

posMapSoftClip *almostPerfect_match_seq_ref( parameters *params, int chr_index, char *str, int pos)
{
	int i, index, posMapSize, posMap[10000], hammingDisMap[10000];
	char orient[10000];// orient of the mapping
	int dist, reverseMatch;
	lociInRef *ptr;
	char seed[HASHKMERLEN + 1];
	char *strRev;
	posMapSoftClip *tmpSoftClipMap, *returnPtr;

	int *hash_ptr;
	int num_hits;
	int cnt_hits;
	int str_length;
	int dist_max;
	
	returnPtr = NULL;

	str_length = strlen(str);
	dist_max = 0.05 * str_length;
	
	strncpy( seed, str, HASHKMERLEN);
	seed[HASHKMERLEN] = '\0';
	
	if ( is_kmer_valid (seed) )
	  {	
	    index = hash_function_ref( seed);
	    num_hits = hash_table_count[index];
	    hash_ptr = hash_table_array[index];
	  }
	else{
	  num_hits = 0;
	  hash_ptr = NULL;
	}
	
	posMapSize = 0;

	for (cnt_hits = 0; cnt_hits < num_hits; cnt_hits++)
	  {
	    if ( abs (hash_ptr[cnt_hits] - pos) < SR_LOOKAHEAD)
	      {
	                dist = hammingDistance( &( params->ref_seq[hash_ptr[cnt_hits]]), str, str_length);
			if( dist <= dist_max)
			  {
 			        posMap[posMapSize] = hash_ptr[cnt_hits];
				orient[posMapSize] = FORWARD;
				hammingDisMap[posMapSize] = dist;
				posMapSize = posMapSize + 1;
			  }
	      }
	  }
	

	/*  deprecated
	while( ptr != NULL)
	{
		if( abs( ptr->pos - pos) < 100000)
		{
			dist = hammingDistance( &( params->ref_seq[ptr->pos]), str, strlen( str));
			if( dist <= ( 0.05 * strlen( str)))
			{
				posMap[posMapSize] = ptr->pos;
				orient[posMapSize] = FORWARD;
				hammingDisMap[posMapSize] = dist;
				posMapSize = posMapSize + 1;
			}
		}
		ptr = ptr->next;
		if( posMapSize > 10)
			ptr = NULL;
			}*/

	if( posMapSize < 11)
	  {
	         strRev = ( char *)getMem( (str_length + 1) * sizeof( char));
		 for( i = 0; i < str_length; i++)
		   {
			if( str[i] == 'A')
				strRev[str_length - i - 1] = 'T';
			else if( str[i] == 'T')
				strRev[str_length - i - 1] = 'A';
			else if( str[i] == 'G')
				strRev[str_length - i - 1] = 'C';
			else if( str[i] == 'C')
				strRev[str_length - i - 1] = 'G';
			else if( str[i] == 'N')
				strRev[str_length - i - 1] = 'N';
		   }
		
		strRev[str_length] = '\0';
		strncpy( seed, strRev, HASHKMERLEN);
		seed[HASHKMERLEN] = '\0';
		
		if (is_kmer_valid( seed))
		  {
		    index = hash_function_ref( seed);
		    num_hits = hash_table_count[index];
		    hash_ptr = hash_table_array[index];
		  }
		else {
		  hash_ptr = NULL;
		  num_hits = 0;
		}

		reverseMatch = 0;
		for (cnt_hits = 0; cnt_hits < num_hits; cnt_hits++)
		  {
		        if ( abs (hash_ptr[cnt_hits] - pos) < SR_LOOKAHEAD)
			  {
				dist = hammingDistance( &( params->ref_seq[hash_ptr[cnt_hits]]), strRev, str_length);
				if( dist <= dist_max)
				  {
				        posMap[posMapSize] = hash_ptr[cnt_hits];
					orient[posMapSize] = REVERSE;
					hammingDisMap[posMapSize] = dist;
					posMapSize = posMapSize + 1;
					reverseMatch = 1;
				  }
			  }
			if( posMapSize > 10)
			  break;
		  }

		freeMem( strRev, str_length+1);
	  }
		
		
/*  deprecated
	if( posMapSize < 11)
	{
		strRev = ( char *)getMem( str_length + 1) * sizeof( char));
		for( i = 0; i < strlen( str); i++)
		{
			if( str[i] == 'A')
				strRev[str_length - i - 1] = 'T';
			else if( str[i] == 'T')
				strRev[str_length - i - 1] = 'A';
			else if( str[i] == 'G')
				strRev[str_length - i - 1] = 'C';
			else if( str[i] == 'C')
				strRev[str_length - i - 1] = 'G';
			else if( str[i] == 'N')
				strRev[str_length - i - 1] = 'N';
		}
		strRev[str_length] = '\0';
		strncpy( seed, strRev, HASHKMERLEN);
		seed[HASHKMERLEN] = '\0';
		
		if (is_kmer_valid( seed)){
		  index = hash_function_ref( seed);
		  ptr = hash_table_SR[index];
		}
		else
		  ptr = NULL;

		reverseMatch = 0;
		while( ptr != NULL)
		{
			if( abs( ptr->pos - pos) < 100000)
			{
				dist = hammingDistance( &( ref_seq_per_chr[ptr->pos]), strRev, str_length);
				if( dist <= ( 0.05 * strlen( strRev)))
				{
					posMap[posMapSize] = ptr->pos;
					orient[posMapSize] = REVERSE;
					hammingDisMap[posMapSize] = dist;
					posMapSize = posMapSize + 1;
					reverseMatch = 1;
				}
			}
			ptr = ptr->next;
			if( posMapSize > 10)
				ptr = NULL;
		}

		freeMem( strRev, str_length+1);
	}
*/

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


void countNumSoftClipInCluster( parameters *params, bam_info* in_bam, int chr_index)
{
	int i, j, countNumSR, posStartSF;
	int *countSoftClip = ( int *)getMem( ( params->this_sonic->chromosome_lengths[chr_index] + 1000) * sizeof( int));
	softClip *ptrSoftClip;

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index] + 1000; i++)
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
				for( j = max( 0, posStartSF - MIN_SOFTCLIP_LEN); j < min( params->this_sonic->chromosome_lengths[chr_index], posStartSF + MIN_SOFTCLIP_LEN); j++)
					countNumSR = countNumSR + countSoftClip[j];
			}
			else if( ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
			{
				posStartSF = ptrSoftClip->pos - ptrSoftClip->opl[ptrSoftClip->opCount-1] + in_bam->libraries[i]->read_length;
				for( j = max( 0, posStartSF - MIN_SOFTCLIP_LEN); j < min( params->this_sonic->chromosome_lengths[chr_index], posStartSF + MIN_SOFTCLIP_LEN); j++)
					countNumSR = countNumSR + countSoftClip[j];
			}

			ptrSoftClip->numSoftClipInConcordance = countNumSR;
			ptrSoftClip = ptrSoftClip->next;
		}
	}
	if( countSoftClip != NULL)
		free( countSoftClip);
}


void mapSoftClipToRef( bam_info* in_bam, parameters* params, int chr_index)
{
	int i;
	char *str;
	softClip *ptrSoftClip;
	int tmp;

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		str = ( char *) getMem( sizeof( char) * in_bam->libraries[i]->read_length);

		ptrSoftClip = in_bam->libraries[i]->listSoftClip;
		while( ptrSoftClip != NULL)
		{
			if( ptrSoftClip->avgPhredQualSoftClip > params->mq_threshold && ptrSoftClip->numSoftClipInConcordance > 1
			    && ptrSoftClip->numSoftClipInConcordance < 100 && ptrSoftClip->softClipString != NULL)
			{
				if( ptrSoftClip->op[0] == BAM_CSOFT_CLIP)
				{
					strncpy( str, ptrSoftClip->softClipString, ptrSoftClip->opl[0]);
					str[ptrSoftClip->opl[0]] = '\0';
				}
				else if( ptrSoftClip->op[ptrSoftClip->opCount - 1] == BAM_CSOFT_CLIP)
				  {
				        tmp = in_bam->libraries[i]->read_length - ptrSoftClip->opl[ptrSoftClip->opCount - 1];
					if (tmp <= 0)
					  tmp = in_bam->libraries[i]->read_length;
					strncpy( str, &( ptrSoftClip->softClipString[tmp]), tmp);
					str[in_bam->libraries[i]->read_length - ptrSoftClip->opl[ptrSoftClip->opCount - 1]] = '\0';
				}
				ptrSoftClip->ptrPosMapSoftClip = almostPerfect_match_seq_ref( params, chr_index, str, ptrSoftClip->pos);
			}
			ptrSoftClip = ptrSoftClip->next;
		}
		if( str != NULL)
			free( str);
	}
}

void addSoftClip( library_properties * library, bam_alignment_region* bam_align, bam1_t* bam_alignment, char* chromosome_nameD)
{
	int i;
	float avgPhredQual = 0;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	softClip *newEl = ( softClip *) getMem( sizeof( softClip));
	newEl->readName = ( char*) getMem( ( strlen( bam_alignment->data) + 7) * sizeof( char));
	strncpy( newEl->readName, bam_alignment->data, strlen( bam_alignment->data));

	newEl->readName[strlen( bam_alignment->data)] = '_';
	newEl->readName[strlen( bam_alignment->data) + 1] = 'S';
	newEl->readName[strlen( bam_alignment->data) + 2] = 'p';
	newEl->readName[strlen( bam_alignment->data) + 3] = 'l';
	newEl->readName[strlen( bam_alignment->data) + 4] = 'i';
	newEl->readName[strlen( bam_alignment->data) + 5] = 't';
	newEl->readName[strlen( bam_alignment->data) + 6] = '\0';

	/* Get the name of the chromosome */
	newEl->chromosome_name = NULL;
	set_str( &(newEl->chromosome_name), chromosome_nameD);

	newEl->pos = bam_align->pos_left;
	newEl->qual = bam_alignment_core.qual;
	newEl->softClipString = ( char *)getMem( ( bam_alignment_core.l_qseq + 1) * sizeof( char));

	/* when the read is in reverse strand, bwa does not map it in reverse order
	 * (important for split reads, we don't use paired-ends) */
	newEl->orient = FORWARD;

	uint8_t *a_qual = bam_get_qual( bam_alignment);

	for( i = 0; i < bam_align->n_cigar; i++)
	{
		newEl->op[i] = bam_cigar_op( bam_align->cigar[i]);
		newEl->opl[i] = bam_cigar_oplen( bam_align->cigar[i]);
	}
	newEl->opCount = bam_align->n_cigar;

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

	if( newEl->op[0] == BAM_CSOFT_CLIP && newEl->opl[0] > MIN_SOFTCLIP_LEN)
	{
		for( i = 0; i < newEl->opl[0]; i++)
			avgPhredQual = avgPhredQual + a_qual[i];

		avgPhredQual = ( float)avgPhredQual / ( float)newEl->opl[0];
	}
	else if( newEl->op[newEl->opCount - 1] == BAM_CSOFT_CLIP && newEl->opl[newEl->opCount - 1] > MIN_SOFTCLIP_LEN)
	{
		for( i = library->read_length - newEl->opl[newEl->opCount - 1] - 1; i < library->read_length; i++)
			avgPhredQual = avgPhredQual + a_qual[i];

		avgPhredQual = ( float)avgPhredQual / newEl->opl[newEl->opCount - 1];
	}

	newEl->softClipString[bam_alignment_core.l_qseq] = '\0';
	newEl->avgPhredQualSoftClip = ( int)floorf( avgPhredQual);
	newEl->ptrPosMapSoftClip = NULL;
	newEl->next = library->listSoftClip;
	library->listSoftClip = newEl;

	sr_cnt_bam++;
}
