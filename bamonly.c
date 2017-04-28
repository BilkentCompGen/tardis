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
#include "vh/vh_divethandler.h"
#include "vh/vh_main.h"
#include "vh/vh_maximalCluster.h"
#include "vh/vh_setcover.h"
#include "vh/vh_intervalhandler.h"
#include "vh_createMaxClusterMEI.h"
#include "splitread.h"

int  STR_SIZE = 100; //maximum size of a constant string, such as 10x barcode length

long del_cnt_bam = 0;
long ins_cnt_bam = 0;
long inv_cnt_bam = 0;
long mei_cnt_bam = 0;
long tandup_cnt_bam = 0;
long sr_cnt_bam = 0;

long del_cnt_div = 0;
long ins_cnt_div = 0;
long inv_cnt_div = 0;
long tandup_cnt_div = 0;
long sr_cnt_div = 0;

char **allReadNameList;


void add_discordant_RR_FF( ref_genome* ref, bam_info *in_bam, parameters* params, bam1_t* bam_alignment, int library_index, int chrID)
{
	//Assuming the BAM file is sorted by position only add the mapping to linked list if pos<mpos
	int len, ed, flag;
	char *rname;
	int cigar_op[100];
	int cigar_opl[100];
	int countOp;

	bam1_core_t bam_alignment_core;
	bam_alignment_core = bam_alignment->core;

	discordantReadMapping_Info *newEl;
	newEl = ( discordantReadMapping_Info *) getMem( sizeof( discordantReadMapping_Info));

	/* Get the read name */
	rname = bam_get_qname( bam_alignment);

	/* Need to be put into into divet row */
	if ( params->ten_x){
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get(bam_alignment, "BX"));
	}

	flag = bam_alignment_core.flag;

	if( bam_alignment_core.pos < bam_alignment_core.mpos)
	{
		newEl->readName = ( char*) getMem( ( strlen(rname) + 1) * sizeof( char));
		strcpy( newEl->readName, rname);

		/* Get the name of the chromosome */
		len = strlen( ref->chrom_names[chrID]);
		newEl->chroName = ( char *) getMem( sizeof( char) * len + 1);
		strcpy( newEl->chroName, ref->chrom_names[chrID]);

		newEl->pos1 = bam_alignment_core.pos;
		newEl->pos2 = bam_alignment_core.mpos;
		newEl->mQual1 = bam_alignment_core.qual;
		newEl->mQual2 = 0;

		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			newEl->orient1 = 'R';
			newEl->orient2 = 'R';
		}
		else if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			newEl->orient1 = 'F';
			newEl->orient2 = 'F';
		}

		inv_cnt_bam++;
		newEl->svType = 'V';
		ed = bam_aux2i(bam_aux_get( bam_alignment, "NM"));
		newEl->editDistance = ed;

		for(countOp=0; countOp<bam_alignment_core.n_cigar; countOp++)
		{
			cigar_op[countOp]=bam1_cigar(bam_alignment)[countOp] & BAM_CIGAR_MASK;
			cigar_opl[countOp]=bam1_cigar(bam_alignment)[countOp] >> BAM_CIGAR_SHIFT;
		}

		newEl->pos1_End = newEl->pos1 + in_bam->libraries[library_index]->read_length;

		if (cigar_op[bam_alignment_core.n_cigar-1]==BAM_CSOFT_CLIP)
			newEl->pos1_End = newEl->pos1_End - cigar_opl[bam_alignment_core.n_cigar-1];
		if (cigar_op[0]==BAM_CSOFT_CLIP)
			newEl->pos1_End = newEl->pos1_End - cigar_opl[0];

		newEl->pos1_End=newEl->pos1_End-1;
		newEl->pos1=newEl->pos1+1;

		newEl->next = in_bam->libraries[library_index]->listRR_FF_Mapping;
		in_bam->libraries[library_index]->listRR_FF_Mapping = newEl;
	}

	else if( bam_alignment_core.pos > bam_alignment_core.mpos)
	{
		newEl = in_bam->libraries[library_index]->listRR_FF_Mapping;

		while( newEl != NULL && strcmp( newEl->readName, rname) != 0)
			newEl = newEl->next;

		if( newEl != NULL)
		{
			ed = bam_aux2i(bam_aux_get( bam_alignment, "NM"));
			newEl->editDistance += ed;

			for(countOp=0; countOp<bam_alignment_core.n_cigar; countOp++)
			{
				cigar_op[countOp]=bam1_cigar(bam_alignment)[countOp] & BAM_CIGAR_MASK;
				cigar_opl[countOp]=bam1_cigar(bam_alignment)[countOp] >> BAM_CIGAR_SHIFT;
			}

			newEl->pos2_End=newEl->pos2 + in_bam->libraries[library_index]->read_length;
			if (cigar_op[bam_alignment_core.n_cigar-1]==BAM_CSOFT_CLIP)
				newEl->pos2_End=newEl->pos2_End-cigar_opl[bam_alignment_core.n_cigar-1];

			if (cigar_op[0]==BAM_CSOFT_CLIP)
				newEl->pos2_End=newEl->pos2_End-cigar_opl[0];


			newEl->pos2_End=newEl->pos2_End-1;
			newEl->pos2=newEl->pos2+1;
			newEl->mQual2 = bam_alignment_core.qual;
		}
	}
}


void add_discordant_RF( ref_genome* ref, bam_info *in_bam, parameters* params, bam1_t* bam_alignment, int library_index, int chrID)
{
	int len, ed, flag;
	char *rname;
	int cigar_op[100];
	int cigar_opl[100];
	int countOp;

	bam1_core_t bam_alignment_core;
	bam_alignment_core = bam_alignment->core;
	discordantReadMapping_Info *newEl;
	newEl = ( discordantReadMapping_Info *) malloc( sizeof( discordantReadMapping_Info));

	/* Get the read name */
	rname = bam_get_qname( bam_alignment);
	flag = bam_alignment_core.flag;

	/* Need to be put into into divet row */        
	if ( params->ten_x){
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get(bam_alignment, "BX"));
	}

	if( bam_alignment_core.pos <= bam_alignment_core.mpos)
	{
		newEl->readName = ( char*) getMem( ( strlen( rname) + 1) * sizeof( char));
		strcpy( newEl->readName, rname);

		/* Get the name of the chromosome */
		len = strlen( ref->chrom_names[chrID]);
		newEl->chroName = ( char *) getMem( sizeof ( char) * len + 1);
		strcpy( newEl->chroName, ref->chrom_names[chrID]);

		newEl->pos1 = bam_alignment_core.pos;
		newEl->pos2 = bam_alignment_core.mpos;
		newEl->mQual1 = bam_alignment_core.qual;
		newEl->mQual2 = 0;

		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			newEl->orient1 = 'R';
			newEl->orient2 = 'F';
		}
		else
			fprintf(stderr, "Error in add_discordant_RF()\n");

		ed = bam_aux2i(bam_aux_get( bam_alignment, "NM"));
		newEl->editDistance = ed;

		tandup_cnt_bam++;
		newEl->svType = 'E';


		for( countOp=0; countOp<bam_alignment_core.n_cigar; countOp++)
		{
			cigar_op[countOp]=bam1_cigar(bam_alignment)[countOp] & BAM_CIGAR_MASK;
			cigar_opl[countOp]=bam1_cigar(bam_alignment)[countOp] >> BAM_CIGAR_SHIFT;
		}

		newEl->pos1_End = newEl->pos1 + in_bam->libraries[library_index]->read_length;

		if (cigar_op[bam_alignment_core.n_cigar-1]==BAM_CSOFT_CLIP)
			newEl->pos1_End = newEl->pos1_End - cigar_opl[bam_alignment_core.n_cigar-1];
		if (cigar_op[0]==BAM_CSOFT_CLIP)
			newEl->pos1_End = newEl->pos1_End - cigar_opl[0];

		newEl->pos1_End=newEl->pos1_End-1;
		newEl->pos1=newEl->pos1+1;

		newEl->next = in_bam->libraries[library_index]->listRF_Mapping;
		in_bam->libraries[library_index]->listRF_Mapping = newEl;
	}

	else if( bam_alignment_core.pos > bam_alignment_core.mpos)
	{
		newEl = in_bam->libraries[library_index]->listRF_Mapping;
		while( newEl != NULL && strncmp( newEl->readName, bam_alignment->data, bam_alignment_core.l_qname) != 0)
			newEl = newEl->next;

		if( newEl != NULL)
		{
			ed = bam_aux2i(bam_aux_get( bam_alignment, "NM"));
			newEl->editDistance += ed;
			for(countOp=0; countOp<bam_alignment_core.n_cigar; countOp++)
			{
				cigar_op[countOp]=bam1_cigar(bam_alignment)[countOp] & BAM_CIGAR_MASK;
				cigar_opl[countOp]=bam1_cigar(bam_alignment)[countOp] >> BAM_CIGAR_SHIFT;
			}

			newEl->pos2_End=newEl->pos2 + in_bam->libraries[library_index]->read_length;
			if (cigar_op[bam_alignment_core.n_cigar-1]==BAM_CSOFT_CLIP)
				newEl->pos2_End=newEl->pos2_End-cigar_opl[bam_alignment_core.n_cigar-1];

			if (cigar_op[0]==BAM_CSOFT_CLIP)
				newEl->pos2_End=newEl->pos2_End-cigar_opl[0];


			newEl->pos2_End=newEl->pos2_End-1;
			newEl->pos2=newEl->pos2+1;
			newEl->mQual2 = bam_alignment_core.qual;
		}
	}
}

void add_discordant_FR( ref_genome* ref, bam_info *in_bam, parameters* params, bam1_t* bam_alignment, int library_index, int svType, int chrID)
{
	int len, ed, flag;
	char *rname;
	int cigar_op[100];
	int cigar_opl[100];
	int countOp;


	bam1_core_t bam_alignment_core;
	bam_alignment_core = bam_alignment->core;
	discordantReadMapping_Info *newEl;
	newEl = ( discordantReadMapping_Info *) getMem( sizeof( discordantReadMapping_Info));
	/* Get the read name */
	rname = bam_get_qname( bam_alignment);

	/* Need to be put into into divet row */
	if ( params->ten_x){
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get(bam_alignment, "BX"));
	}

	flag = bam_alignment_core.flag;

	if( bam_alignment_core.pos < bam_alignment_core.mpos)
	{
		newEl->readName = ( char*) getMem( ( strlen( rname) + 1) * sizeof( char));
		strcpy( newEl->readName, rname);

		/* Get the name of the chromosome */
		len = strlen( ref->chrom_names[chrID]);
		newEl->chroName = ( char *) getMem( sizeof ( char) * len + 1);
		strcpy( newEl->chroName, ref->chrom_names[chrID]);

		newEl->pos1 = bam_alignment_core.pos;
		newEl->pos2 = bam_alignment_core.mpos;
		newEl->mQual1 = bam_alignment_core.qual;
		newEl->mQual2 = 0;

		if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			newEl->orient1 = 'F';
			newEl->orient2 = 'R';
		}
		else
		{
			fprintf(stderr, "Error in add_discordant_FR()\n");
		}

		if(svType == RPDEL)
		{
			newEl->svType = 'D';
			del_cnt_bam++;
		}
		if(svType == RPINS)
		{
			newEl->svType = 'I';
			ins_cnt_bam++;
		}

		ed = bam_aux2i(bam_aux_get( bam_alignment, "NM"));
		newEl->editDistance = ed;

		for(countOp = 0; countOp < bam_alignment_core.n_cigar; countOp++)
		{
			cigar_op[countOp] = bam1_cigar( bam_alignment)[countOp] & BAM_CIGAR_MASK;
			cigar_opl[countOp] = bam1_cigar( bam_alignment)[countOp] >> BAM_CIGAR_SHIFT;
		}

		newEl->pos1_End = newEl->pos1 + in_bam->libraries[library_index]->read_length;

		if( cigar_op[bam_alignment_core.n_cigar-1] == BAM_CSOFT_CLIP)
			newEl->pos1_End = newEl->pos1_End - cigar_opl[bam_alignment_core.n_cigar - 1];

		if( cigar_op[0] == BAM_CSOFT_CLIP)
			newEl->pos1_End = newEl->pos1_End - cigar_opl[0];

		newEl->pos1_End = newEl->pos1_End - 1;
		newEl->pos1 = newEl->pos1 + 1;


		newEl->next = in_bam->libraries[library_index]->listFR_Mapping;
		in_bam->libraries[library_index]->listFR_Mapping = newEl;
	}
	else if( bam_alignment_core.pos > bam_alignment_core.mpos)
	{
		newEl = in_bam->libraries[library_index]->listFR_Mapping;
		while( newEl != NULL && strncmp( newEl->readName, bam_alignment->data, bam_alignment_core.l_qname) != 0)
			newEl = newEl->next;

		if( newEl != NULL)
		{
			ed = bam_aux2i(bam_aux_get( bam_alignment, "NM"));
			newEl->editDistance += ed;
			for( countOp = 0; countOp < bam_alignment_core.n_cigar; countOp++)
			{
				cigar_op[countOp] = bam1_cigar( bam_alignment)[countOp] & BAM_CIGAR_MASK;
				cigar_opl[countOp] = bam1_cigar( bam_alignment)[countOp] >> BAM_CIGAR_SHIFT;
			}

			newEl->pos2_End = newEl->pos2 + in_bam->libraries[library_index]->read_length;
			if( cigar_op[bam_alignment_core.n_cigar-1] == BAM_CSOFT_CLIP)
				newEl->pos2_End = newEl->pos2_End - cigar_opl[bam_alignment_core.n_cigar-1];

			if( cigar_op[0] == BAM_CSOFT_CLIP)
				newEl->pos2_End = newEl->pos2_End - cigar_opl[0];

			newEl->pos2_End = newEl->pos2_End - 1;
			newEl->pos2 = newEl->pos2 + 1;
			newEl->mQual2 = bam_alignment_core.qual;
		}
	}
}

void add_discordant_MEI( ref_genome* ref, bam_info * in_bam, parameters* params, bam1_t* bam_alignment, char* mei_subclass, int library_index, int MEI_Type, int chrID)
{
	int len, flag;
	char *rname;
	int countOp, cigar_op, cigar_opl;

	bam1_core_t bam_alignment_core;
	bam_alignment_core = bam_alignment->core;
	discordantReadMappingMEI_Info *newEl;
	newEl = ( discordantReadMappingMEI_Info *) getMem( sizeof( discordantReadMappingMEI_Info));

	/* Get the read name */
	rname = bam_get_qname(bam_alignment);

	/* Need to be put into into divet row */        
	if ( params->ten_x){
		newEl->ten_x_barcode = encode_ten_x_barcode(bam_aux_get(bam_alignment, "BX"));
	}

	flag = bam_alignment_core.flag;

	newEl->readName = ( char*) getMem( ( strlen( rname) + 1) * sizeof( char));
	strcpy( newEl->readName, rname);

	/* Get the name of the chromosome */
	len = strlen( ref->chrom_names[chrID]);
	newEl->chroName = ( char *) getMem( sizeof( char) * len + 1);
	strcpy( newEl->chroName, ref->chrom_names[chrID]);

	/* Get the name of the mei */
	len = strlen( mei_subclass);
	newEl->MEI_subclass = ( char *) getMem( sizeof (char) * len + 1);
	set_str( &newEl->MEI_subclass, mei_subclass);
	free( mei_subclass);

	newEl->pos = bam_alignment_core.pos;
	newEl->qual = bam_alignment_core.qual;

	if( ( flag & BAM_FREVERSE) != 0)
		newEl->orient = 'R';
	else
		newEl->orient='F';

	mei_cnt_bam++;
	newEl->MEI_Type = MEI_Type;

	newEl->pos_End = newEl->pos;
	for( countOp = 0; countOp < bam_alignment_core.n_cigar; countOp++)
	{
		cigar_op = bam1_cigar(bam_alignment)[countOp] & BAM_CIGAR_MASK;
		cigar_opl = bam1_cigar(bam_alignment)[countOp] >> BAM_CIGAR_SHIFT;
		if (cigar_op == BAM_CMATCH)
			newEl->pos_End = newEl->pos_End + cigar_opl;
	}

	newEl->next = in_bam->libraries[library_index]->listMEI_Mapping;
	in_bam->libraries[library_index]->listMEI_Mapping = newEl;
}


int find_mei_bam(ref_genome* ref, char *chroName, char** mei_subclass, int start, int end, int flag)
{
	int ind, len;
	int return_type = NOTMEI;
	char* MEI_Type;
	/* Check if the right end is inside the annotated transposon */
	ind = meiIntervalSearch( ref, chroName, start, end );
	if( ind == 0)
		return NOTMEI;

	len = strlen( g_meiTable[ind].subclass);
	(*mei_subclass) = ( char *) getMem( sizeof (char) * len + 1);
	set_str( mei_subclass, g_meiTable[ind].subclass);

	len = strlen( g_meiTable[ind].superclass);
	MEI_Type = ( char *) getMem( sizeof (char) * len + 1);
	set_str( &MEI_Type, g_meiTable[ind].superclass);

	/*IF THE MEI INSERT IS + STRAND UPPERCASE, IF - STRAND LOWER CASE*/
	if( ( ( flag & BAM_FMREVERSE) == 0 && (g_meiTable[ind].strand[0] == 'C' || g_meiTable[ind].strand[0] == '-'))
			|| ( ( flag & BAM_FMREVERSE) != 0 && g_meiTable[ind].strand[0] == '+'))
	{
		if( MEI_Type[0] == 'A')
			return_type = 0;
		else if( MEI_Type[0] == 'L')
			return_type = 2;
		else if( MEI_Type[0] == 'S')
			return_type = 4;
	}
	else
	{
		if( MEI_Type[0] == 'A')
			return_type = 1;
		else if( MEI_Type[0] == 'L')
			return_type = 3;
		else if( MEI_Type[0] == 'S')
			return_type = 5;
	}
	free( MEI_Type);
	return return_type;
}


void read_bam(bam_info* in_bam, parameters* params, ref_genome* ref, bam_hdr_t* bam_header)
{
	/* Variables */
	int lib_index, meiType, ind, i, window_start, window_end;
	int left_end_id, right_end_id, svType, flag;
	char* library_name = NULL, *chr_name;
	int *bamToRefIndex;
	int cigar_op[100], cigar_opl[100];
	char* mei_subclass;

	bam1_core_t bam_alignment_core;
	bam1_t* bam_alignment;

	bam_alignment = bam_init1();

	/* The array is used to map the chromosome indices in bam file to the ones in reference genome in case there is a difference */
	bamToRefIndex = ( int *) getMem( bam_header->n_targets * sizeof( int));
	for( i = 0; i < bam_header->n_targets; i++){
		bamToRefIndex[i] = find_chr_index_ref( ref, bam_header->target_name[i]);
		//if(bamToRefIndex[i] < 0){
		//            fprintf( stderr, "BAM Only: Read BAM: Couldn't find chrm %s index in the reference :(\n", bam_header->target_name[i]);
		//}
	}
	while( bam_itr_next( in_bam->bam_file, in_bam->iter, bam_alignment) > 0)
	{
		bam_alignment_core = bam_alignment->core;
		left_end_id = bamToRefIndex[bam_alignment_core.tid];
		right_end_id = bamToRefIndex[bam_alignment_core.mtid];

		/* Get library index */
		set_str( &library_name, bam_aux_get( bam_alignment, "RG"));
		lib_index = find_library_index( in_bam, library_name + 1);

		flag = bam_alignment_core.flag;

		in_bam->read_depth_per_chr[bam_alignment_core.pos]++;
		in_bam->read_count++;


		if( right_end_id >= ref->chrom_count || right_end_id < 0)
			continue;

		if (notInRepeat(ref->chrom_names[left_end_id], bam_alignment_core.pos ) && notInRepeat(ref->chrom_names[right_end_id], bam_alignment_core.mpos ))
		{
			// Remove the paired-ends that both end overlap each other
			if( is_proper( flag) && !( ( left_end_id == right_end_id) && ( abs( bam_alignment_core.pos - bam_alignment_core.mpos) < 100)))
			{
				svType = is_concordant( bam_alignment_core, in_bam->libraries[lib_index]->conc_min, in_bam->libraries[lib_index]->conc_max);

				/* Inversion */
				if( svType == RPINV && left_end_id == right_end_id && abs( bam_alignment_core.pos - bam_alignment_core.mpos) < MAX_INV_LEN)
					add_discordant_RR_FF( ref, in_bam, params, bam_alignment, lib_index, left_end_id);

				if( svType != RPCONC && bam_alignment_core.qual > params->mq_threshold)
				{
					if( ( svType == RPDEL || svType == RPINS) && left_end_id == right_end_id)	/* Deletion or Insertion */
						add_discordant_FR( ref, in_bam, params, bam_alignment, lib_index, svType, left_end_id);

					else if( svType == RPTDUP && left_end_id == right_end_id)	/* Tandem Duplication */
						add_discordant_RF( ref, in_bam, params, bam_alignment, lib_index, left_end_id);

					/* MEI */
					meiType = find_mei_bam( ref, ref->chrom_names[right_end_id], &mei_subclass, bam_alignment_core.mpos,
							bam_alignment_core.mpos + in_bam->libraries[lib_index]->read_length, flag);

					if (meiType != NOTMEI && ( left_end_id != right_end_id ||
							abs( bam_alignment_core.mpos - bam_alignment_core.pos) > MIN_MEI_DISTANCE))
						add_discordant_MEI( ref, in_bam, params, bam_alignment, mei_subclass, lib_index, meiType, left_end_id);
				}
				if( !params->no_soft_clip)
				{
					/* Soft Clipping */
					if( bam_alignment_core.qual > params->mq_threshold && bam_alignment_core.n_cigar < 3)
					{
						for( i = 0; i < bam_alignment_core.n_cigar; i++)
						{
							cigar_op[i] = bam1_cigar(bam_alignment)[i] & BAM_CIGAR_MASK;
							cigar_opl[i] = bam1_cigar(bam_alignment)[i] >> BAM_CIGAR_SHIFT;
						}
						if( ( cigar_op[0] == BAM_CSOFT_CLIP && cigar_opl[0] > MIN_SOFTCLIP_LEN) ||
								( cigar_op[bam_alignment_core.n_cigar - 1] == BAM_CSOFT_CLIP && cigar_opl[bam_alignment_core.n_cigar - 1] > MIN_SOFTCLIP_LEN))
							addSoftClip( ref, in_bam, bam_alignment, lib_index, flag, cigar_op, cigar_opl, left_end_id);
					}
				}
			}
		}
	}
	free(bamToRefIndex);
	fprintf( stderr, "\n%li DEL, %li INV, %li INS, %li TANDUP, %li MEI clusters and %li split reads found in BAM.\n", del_cnt_bam, inv_cnt_bam, ins_cnt_bam, tandup_cnt_bam, mei_cnt_bam, sr_cnt_bam);
	fprintf( logFile, "\n%li DEL, %li INV, %li INS, %li TANDUP, %li MEI clusters and %li split reads found in BAM.\n", del_cnt_bam, inv_cnt_bam, ins_cnt_bam, tandup_cnt_bam, mei_cnt_bam, sr_cnt_bam);
}

void findUniqueReads( bam_info** in_bam, parameters *params, ref_genome* ref, char *outputread)
{
	int i, j, len;
	int totalCountRead = 0;
	long read_name_count;
	long total_read_count;

	softClip *softClipPtr;
	discordantReadMappingMEI_Info *discordantReadPtrMEI;
	discordantReadMapping_Info *discordantReadPtr;

	FILE *fileOutputReadName;
	fileOutputReadName = safe_fopen ( outputread, "w");

	total_read_count = del_cnt_bam + inv_cnt_bam + ins_cnt_bam + tandup_cnt_bam + mei_cnt_bam + sr_cnt_bam;

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
				len = strlen(discordantReadPtr->readName);
				allReadNameList[read_name_count] = ( char *) malloc( (len + 1) * sizeof( char));
				strcpy( allReadNameList[read_name_count], discordantReadPtr->readName);
				read_name_count++;
				discordantReadPtr = discordantReadPtr->next;
			}

			discordantReadPtr = in_bam[i]->libraries[j]->listRR_FF_Mapping;
			while( discordantReadPtr != NULL)
			{
				len = strlen(discordantReadPtr->readName);
				allReadNameList[read_name_count] = ( char *) malloc( (len + 1) * sizeof( char));
				strcpy( allReadNameList[read_name_count], discordantReadPtr->readName);
				read_name_count++;
				discordantReadPtr = discordantReadPtr->next;
			}

			discordantReadPtrMEI = in_bam[i]->libraries[j]->listMEI_Mapping;
			while( discordantReadPtrMEI != NULL)
			{
				len = strlen(discordantReadPtrMEI->readName);
				allReadNameList[read_name_count] = ( char *) malloc( (len + 1) * sizeof( char));
				strcpy( allReadNameList[read_name_count], discordantReadPtrMEI->readName);
				read_name_count++;
				discordantReadPtrMEI = discordantReadPtrMEI->next;
			}

			discordantReadPtr = in_bam[i]->libraries[j]->listRF_Mapping;
			while( discordantReadPtr != NULL)
			{
				len = strlen(discordantReadPtr->readName);
				allReadNameList[read_name_count] = ( char *) malloc( (len + 1) * sizeof( char));
				strcpy( allReadNameList[read_name_count], discordantReadPtr->readName);
				read_name_count++;
				discordantReadPtr = discordantReadPtr->next;
			}
			if( !params->no_soft_clip)
			{
				softClipPtr = in_bam[i]->libraries[j]->listSoftClip;
				while(softClipPtr!=NULL)
				{
					len = strlen(softClipPtr->readName);
					allReadNameList[read_name_count] = ( char *) malloc( (len + 1) * sizeof( char));
					strcpy(allReadNameList[read_name_count], softClipPtr->readName);
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
	fprintf( fileOutputReadName, "%i\n", totalCountRead);

	for( i = 0; i < read_name_count; i++)
	{
		if( i == 0 || strcmp( allReadNameList[i], allReadNameList[i-1]) != 0)
			fprintf( fileOutputReadName, "%s\n", allReadNameList[i]);
	}

	fclose( fileOutputReadName);

	for( i = 0; i < total_read_count; i++)
	{
		if( allReadNameList[i] != NULL)
			free( allReadNameList[i]);
	}
	if( allReadNameList != NULL)
		free( allReadNameList);
}


void load_Divet_bam( bam_info** in_bams, ref_genome* ref, parameters *params, int chr_index)
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
			strcpy (newLibInfo->indName, in_bams[bam_index]->sample_name);
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

			divet_row_count += read_Divet_bam( in_bams[bam_index]->libraries[lib_cnt]->listFR_Mapping, params, ref, newLibInfo, chr_index, divet_row_count);
			divet_row_count += read_Divet_bam( in_bams[bam_index]->libraries[lib_cnt]->listRF_Mapping, params, ref, newLibInfo, chr_index, divet_row_count);
			divet_row_count += read_Divet_bam( in_bams[bam_index]->libraries[lib_cnt]->listRR_FF_Mapping, params, ref, newLibInfo, chr_index, divet_row_count);

			if(g_libInfo == NULL)
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

				read_Divet_bam_softClip( in_bams[bam_index]->libraries[lib_cnt]->listSoftClip, params, ref, newLibInfo, chr_index, newLibInfo->readLen, divet_row_count);

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
}

void free_mappings( bam_info** in_bams, ref_genome* ref, parameters* params)
{
	int i, j, bam_index;
	struct LibraryInfo *cursor, *t;

	softClip *sfcPtr, *sfcPtrNext;
	discordantReadMapping_Info *ptrDisMap, *ptrDisMapNext;
	discordantReadMappingMEI_Info *ptrMEIMap, *ptrMEIMapNext;

	for( bam_index = 0; bam_index < params->num_bams; bam_index++)
	{
		for( i = 0; i < in_bams[bam_index]->num_libraries; i++)
		{
			ptrDisMap = in_bams[bam_index]->libraries[i]->listRR_FF_Mapping;
			while( ptrDisMap != NULL)
			{
				ptrDisMapNext=ptrDisMap->next;
				if( ptrDisMap->readName != NULL)
					free( ptrDisMap->readName);
				if( ptrDisMap->chroName!=NULL)
					free( ptrDisMap->chroName);
				if( ptrDisMap!=NULL)
					free( ptrDisMap);
				ptrDisMap = ptrDisMapNext;
			}
			in_bams[bam_index]->libraries[i]->listRR_FF_Mapping=NULL;

			ptrDisMap = in_bams[bam_index]->libraries[i]->listRF_Mapping;
			while( ptrDisMap!=NULL)
			{
				ptrDisMapNext = ptrDisMap->next;
				if( ptrDisMap->readName != NULL)
					free( ptrDisMap->readName);
				if( ptrDisMap->chroName!=NULL)
					free( ptrDisMap->chroName);
				if( ptrDisMap!=NULL)
					free( ptrDisMap);
				ptrDisMap = ptrDisMapNext;
			}
			in_bams[bam_index]->libraries[i]->listRF_Mapping = NULL;

			ptrDisMap = in_bams[bam_index]->libraries[i]->listFR_Mapping;
			while( ptrDisMap != NULL)
			{
				ptrDisMapNext = ptrDisMap->next;
				if(ptrDisMap->chroName != NULL)
					free( ptrDisMap->chroName);
				if( ptrDisMap->readName != NULL)
					free( ptrDisMap->readName);
				if( ptrDisMap!=NULL)
					free( ptrDisMap);
				ptrDisMap = ptrDisMapNext;
			}
			in_bams[bam_index]->libraries[i]->listFR_Mapping = NULL;

			ptrMEIMap = in_bams[bam_index]->libraries[i]->listMEI_Mapping;
			while( ptrMEIMap != NULL)
			{
				ptrMEIMapNext = ptrMEIMap->next;
				if( ptrMEIMap->chroName != NULL)
					free( ptrMEIMap->chroName);
				if( ptrMEIMap->readName != NULL)
					free( ptrMEIMap->readName);
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
					if( sfcPtr->chroName != NULL)
						free( sfcPtr->chroName);
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

	for( bam_index = 0; bam_index < params->num_bams; bam_index++)
	{
		/* Free the read depth array*/
		free(in_bams[bam_index]->read_depth_per_chr);
	}
	free( ref->gc_per_chr);

	del_cnt_bam = 0;
	ins_cnt_bam = 0;
	inv_cnt_bam = 0;
	mei_cnt_bam = 0;
	tandup_cnt_bam = 0;
	sr_cnt_bam = 0;

	del_cnt_div = 0;
	ins_cnt_div = 0;
	inv_cnt_div = 0;
	tandup_cnt_div = 0;
	sr_cnt_div = 0;

	libInfo = g_libInfo;
	while( libInfo!= NULL)
	{
		libInfoNext = libInfo->next;
		for( i = 0; i < NHASH; i++)
			if( libInfo->hash != NULL)
				if( libInfo->hash[i] != NULL)
					free( libInfo->hash[i]);
		if( libInfo->hash != NULL)
			free( libInfo->hash);
		libInfo = libInfoNext;
	}
}


void bamonly_vh_clustering( bam_info** in_bams, ref_genome* ref, parameters *params)
{
	int i, bam_index, chr_index, chr_index_bam, return_value;
	int total_sv = 0, total_sv_lowqual = 0;
	char outputfile[MAX_SEQ];
	char outputread[MAX_SEQ];
	char svfile[MAX_SEQ];
	FILE *fpVcf = NULL;

	bam_hdr_t* bam_header;
	bam1_t*	bam_alignment;
	bam1_core_t bam_alignment_core;
	hts_idx_t *idx;

	sprintf( outputread, "%s.name", params->outprefix);
	sprintf( outputfile, "%s.clusters", params->outprefix);

	/* Print all structural variations in .vcf format */
	sprintf( svfile, "%s.vcf", params->outprefix);
	fpVcf = safe_fopen( svfile,"w");
	print_vcf_header(fpVcf, in_bams, params);


	for( chr_index = 0; chr_index < ref->chrom_count; chr_index++)
	{
		if (chr_index < params->first_chrom){
			chr_index = params->first_chrom;
		}
		if (chr_index > params->last_chrom){
			chr_index = ref->chrom_count;
			continue;
		}

		if( strstr( ( ref->chrom_names[chr_index]), "_") || strstr( ( ref->chrom_names[chr_index]), "M") || strstr( ( ref->chrom_names[chr_index]), "GL"))
			continue;

		fprintf( stderr, "\nCalculating GC profile");
		calc_gc_per_chr( &ref, chr_index);

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
			bam_header = bam_hdr_read( ( in_bams[bam_index]->bam_file->fp).bgzf);

			/* Load the bam index file */
			idx = sam_index_load( in_bams[bam_index]->bam_file, params->bam_file_list[bam_index]);
			if( idx == NULL)
			{
				fprintf( stderr, "Error: Sam Index cannot be loaded (sam_index_load)\n");
				exit( 1);
			}

			chr_index_bam = find_chr_index_bam(ref, ref->chrom_names[chr_index], bam_header);
			if( chr_index_bam == -1)
			{
				fprintf( stderr, "Error: Cannot find the chromosome name %s in bam\n", ref->chrom_names[chr_index]);
				exit( 1);
			}

			in_bams[bam_index]->iter = bam_itr_queryi( idx, chr_index_bam, 0, ref->chrom_lengths[chr_index]);
			if( in_bams[bam_index]->iter == NULL)
			{
				//continue;
				fprintf( stderr, "Error: Iterator cannot be loaded (bam_itr_queryi)\n");
				exit( 1);
			}
			ref->in_bam[chr_index] = true;

			fprintf( stderr, "\n                                                        ");
			fflush( stderr);
			fprintf( stderr, "\nReading BAM [%s] - chromosome %s", in_bams[bam_index]->sample_name, bam_header->target_name[chr_index_bam]);
			fflush( stderr);

			/* Initialize the read depth and read count */
			init_rd_per_chr( in_bams[bam_index], params, ref, chr_index);

			/* Read bam file for this chromosome */
			read_bam( in_bams[bam_index], params, ref, bam_header);

			if( !params->no_soft_clip)
			{
				/* Count the number of softclip reads which are clustering for each read */
				fprintf( stderr, "\nCollecting soft clipped read infromation");
				countNumSoftClipInCluster( params, ref, in_bams[bam_index], chr_index);
				fprintf( stderr, "\nRemapping soft clipped reads");
				mapSoftClipToRef( in_bams[bam_index], params, ref, chr_index);
			}

			/* Mean value (mu) calculation */
			calc_mean_per_chr( ref, in_bams[bam_index], chr_index);

			/* Close the BAM file */
			return_value = hts_close( in_bams[bam_index]->bam_file);
			if( return_value != 0)
			{
				fprintf( stderr, "Error closing BAM file\n");
				exit( 1);
			}
		}
		if( !params->no_soft_clip)
		{
			/* Free the hash */
			free_10bp_HashIndex();
		}

		/* Free the bam related files */
		bam_itr_destroy( in_bams[bam_index - 1]->iter);
		bam_hdr_destroy( bam_header);
		hts_idx_destroy(idx);

		load_Divet_bam( in_bams, ref, params, chr_index);

		/* Open the .clusters file */
		fileOutput = safe_fopen( outputfile, "w");

		/* Deletion */
		fprintf( stderr, "\nPreparing Deletion clusters");
		vh_initializeReadMapping_Deletion( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
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
		vh_initializeReadMapping_Inversion( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
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
		vh_initializeReadMapping_Insertion( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
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
		vh_initializeReadMapping_TDup( ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
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
		initializeReadMapping_MEI( in_bams, ref, params, ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		MEICluster_Region( in_bams, ref->chrom_names[chr_index], ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);
		vh_finalizeReadMapping_Mei( ref->chrom_lengths[chr_index]);
		fprintf( stderr, ".");
		fflush( stderr);

		fprintf( stderr, "\n");
		fclose( fileOutput);

		findUniqueReads( in_bams, params, ref, outputread);

		/* Free the mappings in bam structure */
		free_mappings( in_bams, ref, params);

		fprintf( stderr, "\nApplying set cover\n");
		vh_setcover( in_bams, params, ref, outputread, outputfile, fpVcf);
		total_sv += sv_count;
		total_sv_lowqual += sv_lowqual_count;

		freeAll( in_bams, ref, params);
	}
	fprintf( stderr, "\n");
	fclose( fpVcf);

	fflush( stderr);
	fprintf( stderr, "\n");

	if ( !TARDIS_DEBUG)
		remove( outputfile);
	fprintf( stderr, "TARDIS is complete. Found %d SVs and %d LowQual.", total_sv, total_sv_lowqual);
	print_sv_stats();
}

int bamonly_run( ref_genome* ref, parameters *params, bam_info ** in_bams)
{
	int rd_del_filtered, bam_index;
	int sv_total;

	char *gapFileName = params->gaps;
	char *repeatFileName = params->reps;

	/* Initialize and read bam file */
	fprintf( stderr, "Processing bam file for read pair and read depth filtering\n\n"
			"RD Threshold = %d\nMapping Quality Threshold = %d\nRP Support Threshold = %d\n\n"
			, params->rd_threshold, params->mq_threshold, params->rp_threshold);

	fprintf( logFile,"\n--> Processing bam file for read pair and read depth filtering\n\n"
			"RD Threshold = %d\nMapping Quality Threshold = %d\nRP Support Threshold = %d\n\n"
			, params->rd_threshold, params->mq_threshold, params->rp_threshold);

	/* Initialize the tables */
	vh_readGapTable( gapFileName);
	vh_readRepeatTable( repeatFileName, params->mei);
	vh_writeMeiIndices( ref);

	bamonly_vh_clustering( in_bams, ref, params);

	return RETURN_SUCCESS;
}
