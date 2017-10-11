/*
 * bo_bamhandler.c
 *
 *  Created on: Aug 23, 2016
 *      Author: tardis
 */

#include "bamonly.h"
#include "common.h"
#include <math.h>
#include "variants.h"
#include "vh/vh_common.h"
#include "vh/vh_maximalCluster.h"
#include "vh/vh_setcover.h"
#include "vh_createMaxClusterMEI.h"
#include "splitread.h"
#include "processfq.h"

long del_cnt_bam = 0;
long ins_cnt_bam = 0;
long inv_cnt_bam = 0;
long mei_cnt_bam = 0;
long tandup_cnt_bam = 0;
long sr_cnt_bam = 0;
long total_read_count = 0;

long del_cnt_div = 0;
long ins_cnt_div = 0;
long inv_cnt_div = 0;
long tandup_cnt_div = 0;
long sr_cnt_div = 0;

char **allReadNameList;

int is_concordant_quick( bam_alignment_region* bam_align, int min, int max)
{
	int flag = bam_align->flag;

	if( ( flag & BAM_FPAIRED) == 0)
	{
		/* Read is single-end. Skip this by calling it concordant */
		return RPCONC;
	}
	/*
	if( ( flag & BAM_FPROPER_PAIR) == 0)
	{
		/* Not proper pair
		return RPUNMAPPED;
	}*/

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
	{
		/* -- orientation = inversion */
		return RPINV;
	}

	if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
	{
		/* ++ orientation = inversion */
		return RPINV;
	}

	if( bam_align->chrID_left != bam_align->chrID_right)
	{
		/* On different chromosomes */
		return RPINTERCHR;
	}

	if( bam_align->pos_left <= bam_align->pos_right) // c.a.
	{
		/* Read is placed BEFORE its mate */
		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			/* -+ orientation = tandem duplication */
			return RPTDUP; //was 0 before
		}
	}
	else
	{
		/* Read is placed AFTER its mate */
		if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			/* +- orientation = tandem duplication */
			return RPTDUP; //was 0 before
		}
	}

	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs( bam_align->isize) < min) // c.a.
	{
		/* Deletion or Insertion */
		return RPINS;
	}
	else if( abs( bam_align->isize) > max)
	{
		return RPDEL;
	}

	/* All passed. Read is concordant */
	return RPCONC;
}


void free_alignments( bam_alignment_region* bam_align)
{
	int i;
	bam_alignment_region* bam_align_next;

	while( ( bam_align) != NULL)
	{
		bam_align_next = (bam_align)->next;
		free( bam_align);
		bam_align = bam_align_next;
	}
	bam_align = NULL;
}


/* Get the alternative mapping locations from the XA field of bam file */
void get_alt_mappings( ref_genome* ref, library_properties * library, bam1_core_t bam_alignment_core, bam_hdr_t* bam_header,
		bam_alignment_region** alt_map, char* xa_string, int map_limit)
{
	int i, j, k, cigar_count, mapping_count = 0;
	uint32_t cigar_op, cigar_opl, cigar_opl_shifted, cigar_final;
	char *tok, *tok2, a[4][200];
	char str[1024], num[3];
	bam_alignment_region* new_mapping;

	strcpy( str, xa_string);
	//fprintf(stderr,"\n\n\n\n\n\nMAIN= %s\n", str);

	tok = strchr( str, ';');
	while (tok != NULL)
	{
		*tok++ = '\0';
		i = 0;
		tok2 = strtok(str, ",");

		while (tok2 != NULL)
		{
			strcpy( a[i], tok2);
			i++;
			tok2 = strtok( NULL, ",");
		}
		new_mapping = ( bam_alignment_region*) getMem( sizeof( bam_alignment_region));
		new_mapping->chrID_left = find_chr_index_bam( ref, a[0], bam_header);
		new_mapping->pos_left = atoi( a[1] + 1);
		new_mapping->chrID_right = bam_alignment_core.mtid;
		new_mapping->pos_right = bam_alignment_core.mpos;
		new_mapping->flag = bam_alignment_core.flag;

		//fprintf(stderr,"%d - ",new_mapping->flag);
		if( a[1][0] == '-')
			new_mapping->flag = new_mapping->flag | BAM_FREVERSE;
		else if( a[1][0] == '+')
			new_mapping->flag = new_mapping->flag & 0xFFEF;

		if( new_mapping->chrID_left != new_mapping->chrID_right)
			new_mapping->isize = 0;

		/*If left pair is forward */
		else if(  new_mapping->flag & BAM_FREVERSE == 0)
		{
			/*If right pair is reverse */
			if( bam_alignment_core.flag & BAM_FMREVERSE != 0)
				new_mapping->isize = (new_mapping->pos_right + library->read_length)  - (new_mapping->pos_left);
			else
				new_mapping->isize = new_mapping->pos_right - new_mapping->pos_left;
		}
		/*If left pair is reverse */
		else if(  new_mapping->flag & BAM_FREVERSE != 0)
		{
			/*If right pair is reverse */
			if(  bam_alignment_core.flag & BAM_FMREVERSE != 0)
				new_mapping->isize = (new_mapping->pos_right + library->read_length) - (new_mapping->pos_left + library->read_length);
			else
				new_mapping->isize = (new_mapping->pos_right) - (new_mapping->pos_left + library->read_length);
		}

		/* Cigar */
		j = 0;
		cigar_count = 0;
		new_mapping->cigar = (uint32_t*) getMem( 20 * sizeof( uint32_t));

		while( a[2][j] != '\0')
		{
			k = 0;
			while( isdigit( a[2][j]))
			{
				num[k] = a[2][j];
				k++;
				j++;
			}
			num[k] = '\0';
			cigar_opl = atoi( num);
			cigar_opl_shifted = cigar_opl << 4;

			switch( a[2][j])
			{
			case 'M': cigar_op = 0; break;
			case 'I': cigar_op = 1; break;
			case 'D': cigar_op = 2;	break;
			case 'N': cigar_op = 4; break;
			}

			cigar_final = cigar_op | cigar_opl_shifted;
			//fprintf(stderr, "%s - op=0x%.8X opl=0x%.8X opl_shifted=0x%.8X final=0x%.8X\n\n",a[2], cigar_op, cigar_opl, cigar_opl_shifted, cigar_final );

			j++;
			new_mapping->cigar[cigar_count] = cigar_final;
			cigar_count++;
		}

		new_mapping->n_cigar = cigar_count;
		new_mapping->edit_distance = atoi( a[3]);
		new_mapping->next = *alt_map;
		*alt_map = new_mapping;

		//fprintf(stderr,"INNN %d %d %d\n", new_mapping->chrID_left, new_mapping->pos_left, new_mapping->n_cigar);

		strcpy( str, tok);
		tok = strchr( str, ';');

		mapping_count++;
		if( mapping_count == map_limit)
			break;
	}
}


void findUniqueReads( bam_info** in_bam, parameters *params, ref_genome* ref, char *outputread)
{
	int i, j;
	int totalCountRead = 0;
	long read_name_count;

	softClip *softClipPtr;
	discordantMappingMEI *discordantReadPtrMEI;
	discordantMapping *discordantReadPtr;

	FILE *fileOutputReadName;

	if( debug_mode)
		fileOutputReadName = safe_fopen ( outputread, "w");

	allReadNameList = ( char **) getMem ( ( total_read_count + 1) * sizeof ( char *));
	if( allReadNameList == NULL)
		fprintf( stderr,"Memory problem\n");

	read_name_count = 0;

	/* Put the names of all the reads in allReadNameList and write them to .name file in sorted order */
	for( i = 0; i < params->num_bams; i++)
	{
		for( j = 0; j < in_bam[i]->num_libraries; j++)
		{
			discordantReadPtr = in_bam[i]->libraries[j]->listFR_Mapping;
			while( discordantReadPtr != NULL)
			{
				allReadNameList[read_name_count] = NULL;
				set_str( &(allReadNameList[read_name_count]), discordantReadPtr->readName);
				read_name_count++;
				discordantReadPtr = discordantReadPtr->next;
			}

			discordantReadPtr = in_bam[i]->libraries[j]->listRR_FF_Mapping;
			while( discordantReadPtr != NULL)
			{
				allReadNameList[read_name_count] = NULL;
				set_str( &(allReadNameList[read_name_count]), discordantReadPtr->readName);
				read_name_count++;
				discordantReadPtr = discordantReadPtr->next;
			}

			discordantReadPtrMEI = in_bam[i]->libraries[j]->listMEI_Mapping;
			while( discordantReadPtrMEI != NULL)
			{
				allReadNameList[read_name_count] = NULL;
				set_str( &(allReadNameList[read_name_count]), discordantReadPtrMEI->readName);
				read_name_count++;
				discordantReadPtrMEI = discordantReadPtrMEI->next;
			}

			discordantReadPtr = in_bam[i]->libraries[j]->listRF_Mapping;
			while( discordantReadPtr != NULL)
			{
				allReadNameList[read_name_count] = NULL;
				set_str( &(allReadNameList[read_name_count]), discordantReadPtr->readName);
				read_name_count++;
				discordantReadPtr = discordantReadPtr->next;
			}
			if( !params->no_soft_clip)
			{
				softClipPtr = in_bam[i]->libraries[j]->listSoftClip;
				while(softClipPtr!=NULL)
				{
					allReadNameList[read_name_count] = NULL;
					set_str( &(allReadNameList[read_name_count]), softClipPtr->readName);
					read_name_count++;
					softClipPtr = softClipPtr->next;
				}
			}
		}
	}
	qsort( allReadNameList, read_name_count, sizeof( char *), vh_cmprReadNameStr);

	for( i = 0; i < read_name_count; i++)
	{
		if ( i == 0 || strcmp( allReadNameList[i], allReadNameList[i-1]) != 0)
			totalCountRead++;
	}
	if( debug_mode)
		fprintf( fileOutputReadName, "%i\n", totalCountRead);

	/* Write the read names to read_names structure for set_cover */
	read_names = (readEl *) getMem( ( read_name_count + 1) * sizeof( readEl));
	j = 0;

	for( i = 0; i < read_name_count; i++)
	{
		if( i == 0 || strcmp( allReadNameList[i], allReadNameList[i-1]) != 0)
		{
			if( debug_mode)
				fprintf( fileOutputReadName, "%s\n", allReadNameList[i]);

			read_names[j].readName = NULL;
			set_str( &(read_names[j].readName), allReadNameList[i]);
			read_names[j].readCovered = 0;
			read_names[j].readId = j;
			read_names[j].libId = -1;
			read_names[j].indId = -1;
			read_names[j].next = NULL;
			j++;
		}
	}
	read_names_count = j;

	if( debug_mode)
		fclose( fileOutputReadName);

	for( i = 0; i < total_read_count; i++)
	{
		if( allReadNameList[i] != NULL)
			free( allReadNameList[i]);
	}
	if( allReadNameList != NULL)
		free( allReadNameList);

	total_read_count = 0;
}


int load_Divet_bam( bam_info** in_bams, ref_genome* ref, parameters *params, int chr_index)
{
	int lib_cnt, read_cnt;
	int j, bam_index, chr, divet_row_count = 0;
	struct LibraryInfo *newLibInfo, *cursor, *t;

	g_maxListBrkPointIntr = MAXLISTBRKPOINTINTR;
	g_libInfo = NULL;

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

			divet_row_count = read_Divet_bam( in_bams[bam_index]->libraries[lib_cnt]->listFR_Mapping, params, ref, newLibInfo, chr_index, divet_row_count);
			divet_row_count = read_Divet_bam( in_bams[bam_index]->libraries[lib_cnt]->listRF_Mapping, params, ref, newLibInfo, chr_index, divet_row_count);
			divet_row_count = read_Divet_bam( in_bams[bam_index]->libraries[lib_cnt]->listRR_FF_Mapping, params, ref, newLibInfo, chr_index, divet_row_count);
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
	fprintf( logFile, "%li del, %li inv, %li ins, %li tandup divet rows loaded sucessfully\n", del_cnt_div, inv_cnt_div, ins_cnt_div, tandup_cnt_div);

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

				divet_row_count = read_Divet_bam_softClip( in_bams[bam_index]->libraries[lib_cnt]->listSoftClip, params, ref, newLibInfo, chr_index, newLibInfo->readLen, divet_row_count);

				if( g_libInfo == NULL)
				{
					g_libInfo = newLibInfo;
				}
				else                  //add to the end of the linked list
				{
					for( t = g_libInfo; t->next != NULL; t = t->next)
						;        //Skip till end of LinkedList
					t->next = newLibInfo;
				}
			}
		}
	}
	return divet_row_count;
}

void free_mappings( bam_info** in_bams, ref_genome* ref, parameters* params)
{
	int i, j, bam_index;
	struct LibraryInfo *cursor, *t;

	softClip *sfcPtr, *sfcPtrNext;
	discordantMapping *ptrDisMap, *ptrDisMapNext;
	discordantMappingMEI *ptrMEIMap, *ptrMEIMapNext;

	for( bam_index = 0; bam_index < params->num_bams; bam_index++)
	{
		for( i = 0; i < in_bams[bam_index]->num_libraries; i++)
		{
			ptrDisMap = in_bams[bam_index]->libraries[i]->listRR_FF_Mapping;
			while( ptrDisMap != NULL)
			{
				ptrDisMapNext = ptrDisMap->next;
				if( ptrDisMap->readName != NULL)
					free( ptrDisMap->readName);
				if( ptrDisMap->chromosome_name!=NULL)
					free( ptrDisMap->chromosome_name);
				if( ptrDisMap != NULL)
					free( ptrDisMap);
				ptrDisMap = ptrDisMapNext;
			}
			in_bams[bam_index]->libraries[i]->listRR_FF_Mapping = NULL;

			ptrDisMap = in_bams[bam_index]->libraries[i]->listRF_Mapping;
			while( ptrDisMap != NULL)
			{
				ptrDisMapNext = ptrDisMap->next;
				if( ptrDisMap->readName != NULL)
					free( ptrDisMap->readName);
				if( ptrDisMap->chromosome_name!=NULL)
					free( ptrDisMap->chromosome_name);
				if( ptrDisMap != NULL)
					free( ptrDisMap);
				ptrDisMap = ptrDisMapNext;
			}

			in_bams[bam_index]->libraries[i]->listRF_Mapping = NULL;

			ptrDisMap = in_bams[bam_index]->libraries[i]->listFR_Mapping;
			while( ptrDisMap != NULL)
			{
				ptrDisMapNext = ptrDisMap->next;
				if(ptrDisMap->chromosome_name != NULL)
					free( ptrDisMap->chromosome_name);
				if( ptrDisMap->readName != NULL)
					free( ptrDisMap->readName);
				if( ptrDisMap != NULL)
					free( ptrDisMap);
				ptrDisMap = ptrDisMapNext;
			}
			in_bams[bam_index]->libraries[i]->listFR_Mapping = NULL;

			ptrMEIMap = in_bams[bam_index]->libraries[i]->listMEI_Mapping;
			while( ptrMEIMap != NULL)
			{
				ptrMEIMapNext = ptrMEIMap->next;
				if( ptrMEIMap->chromosome_name != NULL)
					free( ptrMEIMap->chromosome_name);
				if( ptrMEIMap->readName != NULL)
					free( ptrMEIMap->readName);
				if( ptrMEIMap->MEI_subclass != NULL)
					free( ptrMEIMap->MEI_subclass);
				if( ptrMEIMap != NULL)
					free( ptrMEIMap);
				ptrMEIMap = ptrMEIMapNext;
			}
			in_bams[bam_index]->libraries[i]->listMEI_Mapping = NULL;

			if( !params->no_soft_clip)
			{
				sfcPtr = in_bams[bam_index]->libraries[i]->listSoftClip;
				while( sfcPtr != NULL)
				{
					sfcPtrNext = sfcPtr->next;
					if( sfcPtr->readName != NULL)
						free( sfcPtr->readName);
					if( sfcPtr->chromosome_name != NULL)
						free( sfcPtr->chromosome_name);
					if( sfcPtr->softClipString != NULL)
						free( sfcPtr->softClipString);
					free( sfcPtr);
					sfcPtr = sfcPtrNext;
				}
				in_bams[bam_index]->libraries[i]->listSoftClip = NULL;
			}
		}
	}
}

void freeAll( bam_info** in_bams, ref_genome* ref, parameters* params)
{
	int i, j, bam_index;
	struct LibraryInfo *cursor, *t;
	LibraryInfo *libInfo, *libInfoNext;
	DivetRow *tmp, *tmp_next;

	for( bam_index = 0; bam_index < params->num_bams; bam_index++)
	{
		/* Free the read depth array*/
		free( in_bams[bam_index]->read_depth_per_chr);
	}
	//free( ref->gc_per_chr);

	del_cnt_div = 0;
	ins_cnt_div = 0;
	inv_cnt_div = 0;
	tandup_cnt_div = 0;
	sr_cnt_div = 0;

	libInfo = g_libInfo;
	while( libInfo != NULL)
	{
		libInfoNext = libInfo->next;
		for( i = 0; i < NHASH; i++)
		{
			if( libInfo->hash != NULL)
				if( libInfo->hash[i] != NULL)
					free( libInfo->hash[i]);
		}
		if( libInfo->hash != NULL)
			free( libInfo->hash);

		tmp = libInfo->head;
		while( tmp != NULL)
		{
			tmp_next = tmp->next;
			free( tmp->chromosome_name);
			free( tmp);
			tmp = tmp_next;
		}

		libInfo->head = NULL;
		libInfo->tail = NULL;
		free( libInfo);
		libInfo = libInfoNext;
	}
	g_libInfo = NULL;
}

void add_discordant_RR_FF( ref_genome* ref, library_properties *library, parameters* params, bam1_t* bam_alignment, bam_alignment_region* bam_align, int chrID)
{
	int len;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	discordantMapping *newEl;
	newEl = ( discordantMapping *) getMem( sizeof( discordantMapping));

	/* Need to be put into into divet row */
	if ( params->ten_x)
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get( bam_alignment, "BX"));

	/* Get the read name */
	char* read_name = bam_get_qname( bam_alignment);

	if( bam_align->pos_left < bam_align->pos_right)
	{
		newEl->readName = NULL;
		set_str( &(newEl->readName), bam_get_qname(bam_alignment));

		/* Get the name of the chromosome */
		len = strlen( ref->chrom_names[chrID]);
		newEl->chromosome_name = ( char *) getMem( sizeof( char) * len + 1);
		strcpy( newEl->chromosome_name, ref->chrom_names[chrID]);

		newEl->pos1 = bam_align->pos_left;
		newEl->pos2 = bam_align->pos_right;
		newEl->mQual1 = bam_alignment_core.qual;
		newEl->mQual2 = 0;

		if( ( bam_align->flag & BAM_FREVERSE) != 0 && ( bam_align->flag & BAM_FMREVERSE) != 0)
		{
			newEl->orient1 = 'R';
			newEl->orient2 = 'R';
		}
		else if( ( bam_align->flag & BAM_FREVERSE) == 0 && ( bam_align->flag & BAM_FMREVERSE) == 0)
		{
			newEl->orient1 = 'F';
			newEl->orient2 = 'F';
		}

		inv_cnt_bam++;
		newEl->svType = INVERSION;
		newEl->editDistance = bam_align->edit_distance;

		newEl->pos1_End = newEl->pos1 + library->read_length;

		if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
			newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

		if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
			newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( bam_align->cigar[0]);

		newEl->pos1_End = newEl->pos1_End - 1;
		newEl->pos1 = newEl->pos1 + 1;

		newEl->next = library->listRR_FF_Mapping;
		library->listRR_FF_Mapping = newEl;
	}

	else if( bam_align->pos_left > bam_align->pos_right)
	{
		newEl = library->listRR_FF_Mapping;

		while( newEl != NULL && strcmp( newEl->readName, read_name) != 0)
			newEl = newEl->next;

		if( newEl != NULL)
		{
			newEl->editDistance += bam_align->edit_distance;

			newEl->pos2_End = newEl->pos2 + library->read_length;
			if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

			if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( bam_align->cigar[0]);


			newEl->pos2_End = newEl->pos2_End - 1;
			newEl->pos2 = newEl->pos2 + 1;
			newEl->mQual2 = bam_alignment_core.qual;
		}
	}
}


void add_discordant_RF( ref_genome* ref, library_properties *library, parameters* params, bam1_t* bam_alignment, bam_alignment_region* bam_align, int chrID)
{
	int len;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	discordantMapping *newEl;
	newEl = ( discordantMapping *) getMem( sizeof( discordantMapping));

	/* Need to be put into into divet row */
	if ( params->ten_x)
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get( bam_alignment, "BX"));

	if( bam_align->pos_left < bam_align->pos_right)
	{
		newEl->readName = NULL;
		set_str( &(newEl->readName), bam_get_qname(bam_alignment));

		/* Get the name of the chromosome */
		len = strlen( ref->chrom_names[chrID]);
		newEl->chromosome_name = ( char *) getMem( sizeof ( char) * len + 1);
		strcpy( newEl->chromosome_name, ref->chrom_names[chrID]);

		newEl->pos1 = bam_align->pos_left;
		newEl->pos2 = bam_align->pos_right;
		newEl->mQual1 = bam_alignment_core.qual;
		newEl->mQual2 = 0;

		if( ( bam_align->flag & BAM_FREVERSE) != 0 && ( bam_align->flag & BAM_FMREVERSE) == 0)
		{
			newEl->orient1 = 'R';
			newEl->orient2 = 'F';
		}
		else
			fprintf(stderr, "Error in add_discordant_RF()\n");

		newEl->editDistance = bam_align->edit_distance;

		tandup_cnt_bam++;
		newEl->svType = TANDEMDUP;


		newEl->pos1_End = newEl->pos1 + library->read_length;

		if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
			newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

		if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
			newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( bam_align->cigar[0]);

		newEl->pos1_End = newEl->pos1_End - 1;
		newEl->pos1 = newEl->pos1 + 1;

		newEl->next = library->listRF_Mapping;
		library->listRF_Mapping = newEl;
	}

	else if( bam_align->pos_left > bam_align->pos_right)
	{
		newEl = library->listRF_Mapping;
		while( newEl != NULL && strncmp( newEl->readName, bam_alignment->data, bam_alignment_core.l_qname) != 0)
			newEl = newEl->next;

		if( newEl != NULL)
		{
			newEl->editDistance += bam_align->edit_distance;

			newEl->pos2_End = newEl->pos2 + library->read_length;
			if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

			if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( bam_align->cigar[0]);

			newEl->pos2_End = newEl->pos2_End - 1;
			newEl->pos2 = newEl->pos2 + 1;
			newEl->mQual2 = bam_alignment_core.qual;
		}
	}
}

void add_discordant_FR( ref_genome* ref, library_properties *library, parameters* params, bam1_t* bam_alignment, bam_alignment_region* bam_align, int svType, int chrID)
{
	int len;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	discordantMapping *newEl;
	newEl = ( discordantMapping *) getMem( sizeof( discordantMapping));

	/* Need to be put into into divet row */
	if ( params->ten_x)
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get( bam_alignment, "BX"));

	if( bam_align->pos_left < bam_align->pos_right)
	{
		newEl->readName = NULL;
		set_str( &(newEl->readName), bam_get_qname( bam_alignment));

		/* Get the name of the chromosome */
		len = strlen( ref->chrom_names[chrID]);
		newEl->chromosome_name = ( char *) getMem( sizeof ( char) * len + 1);
		strcpy( newEl->chromosome_name, ref->chrom_names[chrID]);

		newEl->pos1 = bam_align->pos_left;
		newEl->pos2 = bam_align->pos_right;
		newEl->mQual1 = bam_alignment_core.qual;
		newEl->mQual2 = 0;

		if( ( bam_align->flag & BAM_FREVERSE) == 0 && ( bam_align->flag & BAM_FMREVERSE) != 0)
		{
			newEl->orient1 = 'F';
			newEl->orient2 = 'R';
		}
		else
		{
			fprintf(stderr, "Error in add_discordant_FR() - flag:%016x\n", bam_align->flag);
		}

		if( svType == RPDEL)
		{
			newEl->svType = DELETION;
			del_cnt_bam++;
		}
		if( svType == RPINS)
		{
			newEl->svType = INSERTION;
			ins_cnt_bam++;
		}

		newEl->editDistance = bam_align->edit_distance;

		newEl->pos1_End = newEl->pos1 + library->read_length;

		if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
			newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

		if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
			newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( bam_align->cigar[0]);

		newEl->pos1_End = newEl->pos1_End - 1;
		newEl->pos1 = newEl->pos1 + 1;

		newEl->pos2_End = newEl->pos2 + library->read_length;

		newEl->next = library->listFR_Mapping;
		library->listFR_Mapping = newEl;

		//fprintf(stderr, "pos1=%d posend=%d ed=%d name=%s chr_name=%s lib=%d qual=%d\n",newEl->pos1, newEl->pos1_End, newEl->editDistance, newEl->readName, newEl->chromosome_name, bam_align->lib_index, newEl->mQual1);
	}

	else if( bam_align->pos_left > bam_align->pos_right)
	{
		newEl = library->listFR_Mapping;
		while( newEl != NULL && strncmp( newEl->readName, bam_alignment->data, bam_alignment_core.l_qname) != 0)
			newEl = newEl->next;

		if( newEl != NULL)
		{
			newEl->editDistance += bam_align->edit_distance;

			newEl->pos2_End = newEl->pos2 + library->read_length;
			if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

			if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( bam_align->cigar[0]);

			newEl->pos2_End = newEl->pos2_End - 1;
			newEl->pos2 = newEl->pos2 + 1;
			newEl->mQual2 = bam_alignment_core.qual;
		}
	}

}

void add_discordant_MEI( ref_genome* ref, library_properties *library, parameters* params, bam1_t* bam_alignment, bam_alignment_region* bam_align, char* mei_subclass, char* mei_class, int MEI_Type, int chrID)
{
	int len;
	int countOp;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	discordantMappingMEI *newEl;
	newEl = ( discordantMappingMEI *) getMem( sizeof( discordantMappingMEI));


	/* Need to be put into into divet row */        
	if ( params->ten_x){
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get( bam_alignment, "BX"));
	}

	newEl->readName = NULL;
	set_str( &(newEl->readName), bam_get_qname(bam_alignment));

	/* Get the name of the chromosome */
	len = strlen( ref->chrom_names[chrID]);
	newEl->chromosome_name = NULL;
	set_str( &(newEl->chromosome_name), ref->chrom_names[chrID]);

	/* Get the name of the mei subclass */
	len = strlen( mei_subclass);
	newEl->MEI_subclass = NULL;

	set_str( &(newEl->MEI_subclass), mei_subclass);

	/* Get the name of the mei class */
	len = strlen( mei_class);
	newEl->MEI_class = NULL;

	set_str( &(newEl->MEI_class), mei_class);

	newEl->pos = bam_align->pos_left;
	newEl->qual = bam_alignment_core.qual;

	if( ( bam_align->flag & BAM_FREVERSE) != 0)
		newEl->orient = 'R';
	else
		newEl->orient = 'F';

	mei_cnt_bam++;
	newEl->MEI_Type = MEI_Type;

	newEl->pos_End = newEl->pos;
	for( countOp = 0; countOp < bam_align->n_cigar; countOp++)
	{
		if( bam_cigar_opchr( bam_align->cigar[countOp]) == 'S')
			newEl->pos_End = newEl->pos_End + bam_cigar_oplen( bam_align->cigar[countOp]);
	}

	newEl->next = library->listMEI_Mapping;
	library->listMEI_Mapping = newEl;
}


int find_mei_bam( ref_genome* ref, parameters *params, char *chromosome_name, char** mei_subclass, char** mei_class, int start, int end, int flag)
{
	int ind, len;
	int return_type = NOTMEI;

	sonic_repeat *repeat_item;

	/* Check if the right end is inside the annotated transposon */
	repeat_item = sonic_is_mobile_element( params->this_sonic, chromosome_name, start, end, params->mei );
	if( repeat_item == NULL)
		return NOTMEI;

	(*mei_subclass) = NULL;
	set_str( mei_subclass, repeat_item->repeat_type);

	(*mei_class) = NULL;
	set_str( mei_class, repeat_item->repeat_class);

	/* NOTE: SONIC keeps repeat_class as SINE/Alu, LINE/L1, etc. */
	if( ( ( flag & BAM_FMREVERSE) == 0 && repeat_item->strand == SONIC_STRAND_REV)
			|| ( ( flag & BAM_FMREVERSE) != 0 && repeat_item->strand == SONIC_STRAND_FWD))
	{
		return_type = repeat_item->mei_code * 2;
		/*
		if( repeat_item->repeat_type[0] == 'A')
			return_type = 0;
		else if( repeat_item->repeat_type[0] == 'L')
			return_type = 2;
		else if( repeat_item->repeat_type[0] == 'S')
			return_type = 4;
		 */
	}
	else
	{
		return_type = (repeat_item->mei_code * 2) + 1;
		/*
		if( repeat_item->repeat_type[0] == 'A')
			return_type = 1;
		else if( repeat_item->repeat_type[0] == 'L')
			return_type = 3;
		else if( repeat_item->repeat_type[0] == 'S')
			return_type = 5;
		 */
	}
	return return_type;
}

int find_SV_regions( library_properties *library, parameters* params, ref_genome* ref, bam1_t* bam_alignment, int32_t *bamToRefIndex, bam_alignment_region* bam_align)
{
	int svType, meiType = NOTMEI, left_end_id, right_end_id, i;
	char* mei_subclass, *mei_class;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	/* Make sure the chromosome name is within the correct range */
	if( bam_align->chrID_left >= 0 && bam_align->chrID_left < params->this_sonic->number_of_chromosomes &&
			bam_align->chrID_right >= 0 && bam_align->chrID_right < params->this_sonic->number_of_chromosomes)
	{
		left_end_id = bamToRefIndex[bam_align->chrID_left];
		right_end_id = bamToRefIndex[bam_align->chrID_right];

		if( right_end_id >= params->this_sonic->number_of_chromosomes || right_end_id < 0 ||
				left_end_id >= params->this_sonic->number_of_chromosomes || left_end_id < 0)
			return -1;
	}
	else
		return -1;

	/* Find the SVs */
	if( !sonic_is_satellite( params->this_sonic, ref->chrom_names[left_end_id], bam_align->pos_left, bam_align->pos_left + bam_alignment_core.l_qseq)
			&& !sonic_is_satellite( params->this_sonic, ref->chrom_names[right_end_id], bam_align->pos_right, bam_align->pos_right + bam_alignment_core.l_qseq))
	{
		int insLen = abs( bam_align->pos_left - bam_align->pos_right);

		/* Remove the paired-ends that both end overlap each other */
		if( is_proper( bam_align->flag) && !( ( left_end_id == right_end_id) && ( insLen < 100)))
		{
			svType = is_concordant_quick( bam_align, library->conc_min, library->conc_max);
			meiType = find_mei_bam( ref, params, ref->chrom_names[right_end_id], &mei_subclass, &mei_class, bam_align->pos_right,
					bam_align->pos_right + library->read_length, bam_align->flag);

			if( svType != RPCONC && bam_alignment_core.qual > params->mq_threshold)
			{
				/* MEI */
				if( meiType != NOTMEI && ( left_end_id != right_end_id ||
						abs( bam_align->pos_left - bam_align->pos_right) > MIN_MEI_DISTANCE))
					add_discordant_MEI( ref, library, params, bam_alignment, bam_align, mei_subclass, mei_class, meiType, left_end_id);
				/* Deletion or Insertion */
				else if( ( svType == RPDEL || svType == RPINS) && left_end_id == right_end_id)
					add_discordant_FR( ref, library, params, bam_alignment, bam_align, svType, left_end_id);
				/* Tandem Duplication */
				else if( svType == RPTDUP && left_end_id == right_end_id)
				{
					/* Since MT is circular, we need to eliminate the read-pairs at both ends of the chromosome */
					if( insLen < (ref->chrom_lengths[left_end_id]) - (2 * library->conc_max))
						add_discordant_RF( ref, library, params, bam_alignment, bam_align, left_end_id);
				}
			}
			/* Inversion */
			if( svType == RPINV && meiType == NOTMEI && left_end_id == right_end_id && abs( bam_align->pos_left - bam_align->pos_right) < MAX_INV_LEN)
				add_discordant_RR_FF( ref, library, params, bam_alignment, bam_align, left_end_id);

			/* Soft Clipping */
			if( !params->no_soft_clip)
			{
				if( bam_alignment_core.qual > params->mq_threshold && bam_align->n_cigar < MAX_NUM_CIGAR)
				{
					/* We need to have soft clip at the beginning, end, or both with length >MIN_SOFTCLIP_LEN */
					if( ( bam_cigar_opchr( bam_align->cigar[0]) == 'S' && bam_cigar_oplen( bam_align->cigar[0]) > MIN_SOFTCLIP_LEN) ||
							( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S' &&
									bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]) > MIN_SOFTCLIP_LEN))
						addSoftClip( ref, library, bam_align, bam_alignment, left_end_id);
				}
			}
		}
	}
	return 1;
}


void read_bam( bam_info* in_bam, parameters* params, ref_genome* ref)
{
	/* Variables */
	int i, chr_index_bam, return_type, ed, len, lib_index;
	char *xa_str = NULL, *library_name = NULL;
	int32_t bamToRefIndex[in_bam->bam_header->n_targets];
	uint8_t *tmp;
	bam1_core_t bam_alignment_core;
	bam_alignment_region* bam_align = NULL, *bam_align_tmp;

	bam1_t* bam_alignment = bam_init1();

	/* The array is used to map the chromosome indices in bam file to the ones in reference genome */
	for( i = 0; i < in_bam->bam_header->n_targets; i++)
		bamToRefIndex[i] = sonic_refind_chromosome_index( params->this_sonic, in_bam->bam_header->target_name[i]);

	while( bam_itr_next( in_bam->bam_file, in_bam->iter, bam_alignment) > 0)
	{
		bam_alignment_core = bam_alignment->core;

		if( bam_align != NULL)
			free_alignments( bam_align);

		/* Get library index */
		set_str( &library_name, bam_aux_get( bam_alignment, "RG"));
		lib_index = find_library_index( in_bam, library_name + 1);

		/* Put the alignment data to bam_alignment_region data structure */
		bam_align = ( bam_alignment_region*) getMem( sizeof( bam_alignment_region));
		bam_align->chrID_left = bam_alignment_core.tid;
		bam_align->chrID_right = bam_alignment_core.mtid;
		bam_align->pos_left = bam_alignment_core.pos;
		bam_align->pos_right = bam_alignment_core.mpos;
		bam_align->flag = bam_alignment_core.flag;
		bam_align->n_cigar = bam_alignment_core.n_cigar;
		bam_align->cigar = bam_get_cigar( bam_alignment);
		bam_align->isize = bam_alignment_core.isize;
		bam_align->next = NULL;

		tmp = bam_aux_get( bam_alignment, "NM");
		bam_align->edit_distance = 0;
		if( tmp != 0)
			bam_align->edit_distance = bam_aux2i( tmp);

		/* Increase the read depth and read count for RD filtering */
		in_bam->read_depth_per_chr[bam_align->pos_left]++;
		in_bam->read_count++;

		return_type = find_SV_regions( in_bam->libraries[lib_index], params, ref, bam_alignment, bamToRefIndex, bam_align);
		if( return_type == -1)
			continue;

		/* Get alternative mapping locations - XA field */
		if( params->alt_mapping != 0)
		{
			xa_str = NULL;
			set_str( &xa_str, bam_aux_get( bam_alignment, "XA"));

			if( xa_str != NULL)
			{
				get_alt_mappings( ref, in_bam->libraries[lib_index], bam_alignment_core, in_bam->bam_header, &bam_align, ( xa_str + 1), params->alt_mapping);

				bam_align_tmp = bam_align;
				while(	bam_align_tmp != NULL)
				{
					//fprintf(stderr,"%d %d %d\n", bam_align_tmp->chrID_left, bam_align_tmp->flag, bam_align_tmp->edit_distance);
					return_type = find_SV_regions( in_bam->libraries[lib_index], params, ref, bam_alignment, bamToRefIndex, bam_align_tmp);
					bam_align_tmp = bam_align_tmp->next;
				}
				//fprintf(stderr,"\n\n\n");
			}
		}
	}
	bam_destroy1( bam_alignment);
	fprintf( stderr, "\n%li DEL, %li INV, %li INS, %li TANDUP, %li MEI clusters and %li split reads found in BAM.\n", del_cnt_bam, inv_cnt_bam, ins_cnt_bam, tandup_cnt_bam, mei_cnt_bam, sr_cnt_bam);
	fprintf( logFile, "\n%li DEL, %li INV, %li INS, %li TANDUP, %li MEI clusters and %li split reads found in BAM.\n", del_cnt_bam, inv_cnt_bam, ins_cnt_bam, tandup_cnt_bam, mei_cnt_bam, sr_cnt_bam);
}


void bamonly_vh_clustering( bam_info** in_bams, ref_genome* ref, parameters *params)
{
	int i, bam_index, chr_index, chr_index_bam, return_value, not_in_bam = 0;
	int total_sv = 0, total_sv_lowqual = 0, divet_row_count;
	char outputfile[MAX_SEQ];
	char outputread[MAX_SEQ];
	char svfile[MAX_SEQ];
	FILE *fpVcf = NULL;

	sprintf( outputread, "%s.name", params->outprefix);
	sprintf( outputfile, "%s.clusters", params->outprefix);

	/* Print all structural variations in .vcf format */
	sprintf( svfile, "%s.vcf", params->outprefix);
	fpVcf = safe_fopen( svfile,"w");

	print_vcf_header( fpVcf, in_bams, params);

	for( chr_index = 0; chr_index < params->this_sonic->number_of_chromosomes; chr_index++)
	{
		if (chr_index < params->first_chrom)
			chr_index = params->first_chrom;

		if (chr_index > params->last_chrom)
		{
			chr_index = params->this_sonic->number_of_chromosomes;
			continue;
		}

		if( !params->no_soft_clip)
		{
			fprintf( stderr, "\nReading reference genome");
			readReferenceSeq( ref, params, chr_index);
		}

		for ( bam_index = 0; bam_index < params->num_bams; bam_index++)
		{
			/* HTS implementation */
			in_bams[bam_index]->bam_file = safe_hts_open( params->bam_file_list[bam_index], "r");

			/* Read in BAM header information */
			in_bams[bam_index]->bam_header = bam_hdr_read( ( in_bams[bam_index]->bam_file->fp).bgzf);

			/* Load the bam index file */
			in_bams[bam_index]->bam_file_index = sam_index_load( in_bams[bam_index]->bam_file, params->bam_file_list[bam_index]);
			if( in_bams[bam_index]->bam_file_index == NULL)
			{
				fprintf( stderr, "Error: Sam Index cannot be loaded (sam_index_load)\n");
				exit( 1);
			}

			chr_index_bam = find_chr_index_bam( ref, ref->chrom_names[chr_index], in_bams[bam_index]->bam_header);
			not_in_bam = 0;
			if( chr_index_bam == -1)
			{
				fprintf( stderr, "\nCannot find chromosome name %s in bam %s", ref->chrom_names[chr_index], in_bams[bam_index]->sample_name);
				not_in_bam = 1;
				continue;
			}

			in_bams[bam_index]->iter = bam_itr_queryi( in_bams[bam_index]->bam_file_index, chr_index_bam, 0, ref->chrom_lengths[chr_index]);
			if( in_bams[bam_index]->iter == NULL)
			{
				fprintf( stderr, "Error: Iterator cannot be loaded (bam_itr_queryi)\n");
				exit( 1);
			}
			ref->in_bam[chr_index] = true;

			fprintf( stderr, "\n                                                        ");
			fflush( stderr);
			fprintf( stderr, "\nReading BAM [%s] - Chromosome: %s", in_bams[bam_index]->sample_name, in_bams[bam_index]->bam_header->target_name[chr_index_bam]);
			fflush( stderr);

			/* Initialize the read depth and read count */
			init_rd_per_chr( in_bams[bam_index], params, ref, chr_index);

			/* Read bam file for this chromosome */
			read_bam( in_bams[bam_index], params, ref);

			if( !params->no_soft_clip)
			{
				/* Count the number of softclip reads which are clustering for each read */
				fprintf( stderr, "\nCollecting soft clipped read information");
				countNumSoftClipInCluster( params, ref, in_bams[bam_index], chr_index);
				fprintf( stderr, "\nRemapping soft clipped reads");
				mapSoftClipToRef( in_bams[bam_index], params, ref, chr_index);
			}

			/* Mean value (mu) calculation */
			calc_mean_per_chr( ref, params, in_bams[bam_index], chr_index);

			/* Close the BAM file */
			return_value = hts_close( in_bams[bam_index]->bam_file);
			if( return_value != 0)
			{
				fprintf( stderr, "Error closing BAM file\n");
				exit( 1);
			}
			/* Free the bam related files */
			bam_itr_destroy( in_bams[bam_index]->iter);
			bam_hdr_destroy( in_bams[bam_index]->bam_header);
			hts_idx_destroy(in_bams[bam_index]->bam_file_index);

			/* Increase the total read count for this chromosome and make the number of SVs 0 for this bam */
			total_read_count += del_cnt_bam + inv_cnt_bam + ins_cnt_bam + tandup_cnt_bam + mei_cnt_bam + sr_cnt_bam;
			del_cnt_bam = 0, ins_cnt_bam = 0, inv_cnt_bam = 0, mei_cnt_bam = 0, tandup_cnt_bam = 0, sr_cnt_bam = 0;
		}
		if( not_in_bam == 1)
			continue;
		if( !params->no_soft_clip)
		{
			/* Free the hash */
			free_HashIndex();
		}

		divet_row_count = load_Divet_bam( in_bams, ref, params, chr_index);

		/* Open the .clusters file if running in debug mode */
		if( debug_mode)
			fileOutput = safe_fopen( outputfile, "w");

		/* Deletion */
		fprintf( stderr, "\nPreparing Deletion clusters");
		vh_initializeReadMapping_Deletion( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index], params->this_sonic);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_createDeletionClusters( ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_finalizeReadMapping( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);

		/* Inversion */
		fprintf( stderr, "\nPreparing Inversion clusters");
		vh_initializeReadMapping_Inversion( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index], params->this_sonic);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_createInversionClusters( ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_finalizeReadMapping( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);

		/* Insertion */
		fprintf( stderr, "\nPreparing Insertion clusters");
		vh_initializeReadMapping_Insertion( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index], params->this_sonic);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_createInsertionClusters( ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_finalizeReadMapping( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);

		/* Tandem Duplication */
		fprintf( stderr, "\nPreparing Tandem Dup clusters");
		vh_initializeReadMapping_TDup( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index], params->this_sonic);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_createTDupClusters( ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_finalizeReadMapping( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);


		/* Mei filtering */
		fprintf( stderr, "\nPreparing MEI clusters");
		initializeReadMapping_MEI( in_bams, params, ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		MEICluster_Region( params, ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_finalizeReadMapping_Mei( ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);

		fprintf( stderr, "\n");

		if( debug_mode)
			fclose( fileOutput);

		findUniqueReads( in_bams, params, ref, outputread);

		/* Free the mappings in bam structure */
		free_mappings( in_bams, ref, params);

		fprintf( stderr, "\nApplying set cover\n");
		vh_setcover( in_bams, params, ref, fpVcf);
		total_sv += sv_count;
		total_sv_lowqual += sv_lowqual_count;

		freeAll( in_bams, ref, params);
	}
	fprintf( stderr, "\n");
	fclose( fpVcf);

	fflush( stderr);
	fprintf( stderr, "\n");

	if( debug_mode)
		fprintf( stderr, "TARDIS is complete. Found %d SVs and %d LowQual.", total_sv, total_sv_lowqual);
	else
		fprintf( stderr, "TARDIS is complete. Found %d SVs", total_sv);

	print_sv_stats();
}

int bamonly_run( ref_genome* ref, parameters *params, bam_info ** in_bams)
{
	int rd_del_filtered, bam_index;
	int sv_total, i, len;

	/* Initialize and read bam file */
	fprintf( stderr, "Processing bam file for read pair and read depth filtering\n"
			"(RD Threshold: %d; Mapping Quality Threshold: %d; RP Support Threshold: %d)\n\n"
			, params->rd_threshold, params->mq_threshold, params->rp_threshold);

	fprintf( logFile,"\n--> Processing bam file for read pair and read depth filtering\n"
			"(RD Threshold: %d; Mapping Quality Threshold: %d; RP Support Threshold: %d)\n\n"
			, params->rd_threshold, params->mq_threshold, params->rp_threshold);

	bamonly_vh_clustering( in_bams, ref, params);

	return RETURN_SUCCESS;
}
