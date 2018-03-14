/*
 * altmapping.c
 *
 *  Created on: Dec 25, 2017
 *      Author: tardis
 */

#include "mappings.h"
#include "common.h"
#include "bamonly.h"
#include "free.h"

int32_t calculateInsertSize( bam_alignment_region* new_mapping, int read_length)
{
	int32_t isize = 0;

	/* Calculate the insert size */
	if( new_mapping->chrID_left != new_mapping->chrID_right)
		isize = 0;
	else if( ( new_mapping->flag & BAM_FREVERSE) == 0) /* if left pair is forward */
	{
		/*If right pair is reverse */
		if( ( new_mapping->flag & BAM_FMREVERSE) != 0)
			isize = ( new_mapping->pos_right + read_length) - ( new_mapping->pos_left);
		else
			isize = new_mapping->pos_right - new_mapping->pos_left;
	}
	else if( ( new_mapping->flag & BAM_FREVERSE) != 0) /* if left pair is reverse */
	{
		/*If right pair is reverse */
		if( ( new_mapping->flag & BAM_FMREVERSE) != 0)
			isize = ( new_mapping->pos_right + read_length) - ( new_mapping->pos_left + read_length);
		else
			isize = ( new_mapping->pos_right) - ( new_mapping->pos_left + read_length);
	}
	return isize;
}

void calculateCigar( bam_alignment_region** new_mapping, char altmap[])
{
	int j = 0, k;
	int cigar_count = 0;
	char num[3];
	uint32_t cigar_op, cigar_opl, cigar_opl_shifted, cigar_final;
	(*new_mapping)->cigar = ( uint32_t*) getMem( 20 * sizeof( uint32_t));

	while( altmap[j] != '\0')
	{
		k = 0;
		while( isdigit( altmap[j]))
		{
			num[k] = altmap[j];
			k++;
			j++;
		}
		num[k] = '\0';
		cigar_opl = atoi( num);
		cigar_opl_shifted = cigar_opl << BAM_CIGAR_SHIFT;

		switch( altmap[j])
		{
		case 'M': cigar_op = 0; break;
		case 'I': cigar_op = 1; break;
		case 'D': cigar_op = 2;	break;
		case 'N': cigar_op = 3; break;
		case 'S': cigar_op = 4; break;
		case 'H': cigar_op = 5; break;
		case 'P': cigar_op = 6; break;
		case '=': cigar_op = 7; break;
		case 'X': cigar_op = 8; break;
		case 'B': cigar_op = 9; break;
		default: cigar_op = -1; break;
		}
		if(cigar_op == -1)
		{
			j++;
			continue;
		}

		cigar_final = cigar_op | cigar_opl_shifted;
		//fprintf(stderr, "%s - op=0x%.8X opl=0x%.8X opl_shifted=0x%.8X final=0x%.8X\n",altmap, cigar_op, cigar_opl, cigar_opl_shifted, cigar_final );

		j++;
		(*new_mapping)->cigar[cigar_count] = cigar_final;

		cigar_count++;
	}
	//fprintf(stderr,"%s - ", altmap);
	//for(j = 0; j<cigar_count;j++)
		//fprintf(stderr,"%d%c",bam_cigar_oplen( (*new_mapping)->cigar[j]), bam_cigar_opchr( (*new_mapping)->cigar[j]));
	//fprintf(stderr,"\n");

	(*new_mapping)->n_cigar = cigar_count;
}


/* Pair left's alternatives with right's alternatives */
int allLeftAltRightMappings( bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment,
		int32_t *bamToRefIndex, char altmap[4][1024])
{
	bam_alignment_region* new_mapping = NULL;
	discordantMapping* newEl;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	int h = vh_getHash( bam_alignment->data);
	newEl = in_bam->libraries[lib_index]->mappings_discordant[h];

	/* We need to create mappings with all the alternative mappings of the left pair
	 * For these mappings, newEl->side = LEFT */
	while( newEl != NULL)
	{
		if( strncmp( newEl->readName, bam_alignment->data, bam_alignment_core.l_qname) == 0 && ( newEl->side == LEFT))
		{
			/* Left-end is the first read */
			new_mapping = ( bam_alignment_region*) getMem( sizeof( bam_alignment_region));
			new_mapping->read_name = NULL;
			set_str( &(new_mapping->read_name), bam_get_qname( bam_alignment));

			new_mapping->chrID_left = find_chr_index_bam( ref, newEl->chromosome_name, in_bam->bam_header);
			new_mapping->pos_left = newEl->pos1;
			new_mapping->chrID_right = find_chr_index_bam( ref, altmap[0], in_bam->bam_header);
			new_mapping->pos_right = atoi( altmap[1] + 1);

			new_mapping->edit_distance = newEl->editDistance_left;
			new_mapping->qual = newEl->mQual1;
			new_mapping->flag = newEl->flag;
			new_mapping->xa = true;
			new_mapping->side = NONE;

			if( new_mapping->pos_left > new_mapping->pos_right)
			{
				free( new_mapping);
				return -1;
			}

			/* Check the orientation of the altmap and modify the flag */
			if( altmap[1][0] == '-')
				new_mapping->flag = new_mapping->flag | BAM_FMREVERSE;
			else if( altmap[1][0] == '+')
				new_mapping->flag = new_mapping->flag & 0xFFDF;

			/* Calculate the insert size */
			new_mapping->isize = calculateInsertSize( new_mapping, in_bam->libraries[lib_index]->read_length);

			/* If the mapping is unmapped, then continue */
			if( is_concordant_bamonly( new_mapping, in_bam->libraries[lib_index]->conc_min, in_bam->libraries[lib_index]->conc_max) == 0)
			{
				free( new_mapping);
				continue;
			}

			/* Get the Cigar of the left read */
			new_mapping->cigar = ( uint32_t*) getMem( ( newEl->n_cigar + 1) * sizeof( uint32_t));
			memcpy( new_mapping->cigar, newEl->cigar, sizeof( newEl->cigar));
			new_mapping->n_cigar = newEl->n_cigar;
			//fprintf(stderr,"CIGAR = %c %c %c\n", bam_cigar_opchr( new_mapping->cigar[0]), bam_cigar_opchr( new_mapping->cigar[1]), bam_cigar_opchr( new_mapping->cigar[2]));

			read_mapping( in_bam->libraries[lib_index], params, ref, bam_alignment, bamToRefIndex, new_mapping);

			if( new_mapping != NULL)
				free_alignments2( &new_mapping);
		}
		newEl = newEl->next;
	}
	return RETURN_SUCCESS;
}


int altRightAllLeftMappings(  bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment,
		int32_t *bamToRefIndex, char altmap[4][1024])
{
	bam_alignment_region* new_mapping = NULL;
	discordantMapping* newEl;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	int h = vh_getHash( bam_alignment->data);
	newEl = in_bam->libraries[lib_index]->mappings_discordant[h];

	while( newEl != NULL)
	{
		/* We then create the right pair */
		if( strncmp( newEl->readName, bam_alignment->data, bam_alignment_core.l_qname) == 0 &&
				( newEl->pos2 == atoi( altmap[1] + 1)))
		{
			/* Right-end of the xa reads */
			new_mapping = ( bam_alignment_region*) getMem( sizeof( bam_alignment_region));
			new_mapping->read_name = NULL;
			set_str( &(new_mapping->read_name), bam_get_qname( bam_alignment));

			new_mapping->chrID_left = find_chr_index_bam( ref, altmap[0], in_bam->bam_header);
			new_mapping->pos_left = atoi( altmap[1] + 1);
			new_mapping->chrID_right = find_chr_index_bam( ref, newEl->chromosome_name, in_bam->bam_header);
			new_mapping->pos_right = newEl->pos1;

			new_mapping->edit_distance = atoi( altmap[3]);
			new_mapping->flag = bam_alignment_core.flag;
			new_mapping->qual = bam_alignment_core.qual;
			new_mapping->isize = newEl->isize * -1;
			new_mapping->side = NONE;
			new_mapping->xa = true;

			/* Check the orientation of the altmap and modify the flag */
			if( altmap[1][0] == '-')
				new_mapping->flag = new_mapping->flag | BAM_FREVERSE;
			else if( altmap[1][0] == '+')
				new_mapping->flag = new_mapping->flag & 0xFFEF;

			if( newEl->orient1 == REVERSE)
				new_mapping->flag = new_mapping->flag | BAM_FMREVERSE;
			else if( newEl->orient1 == FORWARD)
				new_mapping->flag = new_mapping->flag & 0xFFDF;

			/* If the mapping is unmapped, then continue */
			if( is_concordant_bamonly( new_mapping, in_bam->libraries[lib_index]->conc_min, in_bam->libraries[lib_index]->conc_max) == 0)
			{
				free( new_mapping);
				continue;
			}

			/* Cigar */
			calculateCigar( &new_mapping, altmap[2]);

			read_mapping( in_bam->libraries[lib_index], params, ref, bam_alignment, bamToRefIndex, new_mapping);

			if( new_mapping != NULL)
				free_alignments2( &new_mapping);
		}
		newEl = newEl->next;
	}
	return RETURN_SUCCESS;
}

int altLeftPrimRightMappings(  bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment,
		int32_t *bamToRefIndex, char altmap[4][1024])
{
	bam_alignment_region* new_mapping = NULL;
	bam1_core_t bam_alignment_core = bam_alignment->core;

	new_mapping = ( bam_alignment_region*) getMem( sizeof( bam_alignment_region));
	new_mapping->read_name = NULL;
	set_str( &(new_mapping->read_name), bam_get_qname( bam_alignment));

	new_mapping->chrID_left = find_chr_index_bam( ref, altmap[0], in_bam->bam_header);
	new_mapping->pos_left = atoi( altmap[1] + 1);
	new_mapping->chrID_right = bam_alignment_core.mtid;
	new_mapping->pos_right = bam_alignment_core.mpos;

	new_mapping->edit_distance = atoi( altmap[3]);
	new_mapping->flag = bam_alignment_core.flag;
	new_mapping->qual = bam_alignment_core.qual;
	new_mapping->side = LEFT;
	new_mapping->xa = true;

	if( new_mapping->pos_left > new_mapping->pos_right)
	{
		free( new_mapping);
		return -1;
	}

	if( altmap[1][0] == '-')
		new_mapping->flag = new_mapping->flag | BAM_FREVERSE;
	else if( altmap[1][0] == '+')
		new_mapping->flag = new_mapping->flag & 0xFFEF;

	new_mapping->isize = calculateInsertSize( new_mapping, in_bam->libraries[lib_index]->read_length);

	/* If the mapping is unmapped, then continue */
	if( is_concordant_bamonly( new_mapping, in_bam->libraries[lib_index]->conc_min, in_bam->libraries[lib_index]->conc_max) == 0)
	{
		free( new_mapping);
		return -1;
	}

	/* Cigar */
	calculateCigar( &new_mapping, altmap[2]);
	//fprintf(stderr, "%d - %c%d\n ", new_mapping->n_cigar, bam_cigar_opchr( new_mapping->cigar[0]), bam_cigar_oplen( new_mapping->cigar[0]));


	read_mapping( in_bam->libraries[lib_index], params, ref, bam_alignment, bamToRefIndex, new_mapping);

	if( new_mapping != NULL)
		free_alignments2( &new_mapping);

	return RETURN_SUCCESS;
}

/* Get the alternative mapping locations from the XA field of bam file */
void find_alt_mappings( bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment,
		char* xa_string, int32_t *bamToRefIndex)
{
	int i, return_value;

	char *tok, *tok2, altmap[4][1024];
	char *str;
	char *end_str;
	bam1_core_t bam_alignment_core = bam_alignment->core;

	//strcpy( str, xa_string);
	str = NULL;
	set_str( &str, xa_string);
	//fprintf(stderr,"\nMAIN= %s\n", str);

	tok = strtok_r( str, ";", &end_str);

	/* Each mapping in XA field ends with ; and each mapping has
	 * position(starting with -/+ for orientation, cigar and edit distance */
	while (tok != NULL)
	{
		//*tok++ = '\0';
		i = 0;
		char *end_token;

		tok2 = strtok_r( tok, ",", &end_token);

		while (tok2 != NULL)
		{
			//fprintf(stderr,"%s ", tok2);
			strcpy( altmap[i], tok2);
			tok2 = strtok_r( NULL, ",", &end_token);
			i++;
		}
		//fprintf(stderr,"%s %s %s %s\n", altmap[0], altmap[1], altmap[2], altmap[3]);
		/* We only use the alternative mappings mapped to the same chromosome
		 * Also some alternative mappings have the same location with the primary mapping, these are eliminated*/
		if( find_chr_index_bam( ref, altmap[0], in_bam->bam_header) == bam_alignment_core.tid
				&& (atoi( altmap[1] + 1) != bam_alignment_core.pos) && (atoi( altmap[1] + 1) != bam_alignment_core.pos + 1))
		{
			/* If this is the left pair */
			if( bam_alignment_core.pos < bam_alignment_core.mpos)
			{
				return_value = altLeftPrimRightMappings( in_bam, ref, params, lib_index, bam_alignment, bamToRefIndex, altmap);
			}
			/* If this is the right pair of the mapping */
			else if( bam_alignment_core.pos > bam_alignment_core.mpos)
			{
				allLeftAltRightMappings( in_bam, ref, params, lib_index, bam_alignment, bamToRefIndex, altmap);
				return_value = altRightAllLeftMappings( in_bam, ref, params, lib_index, bam_alignment, bamToRefIndex, altmap);
			}
		}

		strcpy( str, tok);
		tok = strtok_r( NULL, ";", &end_str);
	}
}

/* Get the alternative mapping locations from the XA field of bam file */
int primary_mapping( bam_info* in_bam, ref_genome* ref, parameters* params, int lib_index, bam1_t* bam_alignment, int32_t *bamToRefIndex)
{
	uint8_t *tmp;
	int return_type;

	bam_alignment_region* bam_align = NULL;
	bam1_core_t bam_alignment_core = bam_alignment->core;

	/* Put the alignment data to bam_alignment_region data structure */
	bam_align = ( bam_alignment_region*) getMem( sizeof( bam_alignment_region));

	bam_align->read_name = NULL;
	set_str( &(bam_align->read_name), bam_get_qname( bam_alignment));

	bam_align->chrID_left = bam_alignment_core.tid;
	bam_align->chrID_right = bam_alignment_core.mtid;
	bam_align->pos_left = bam_alignment_core.pos;
	bam_align->pos_right = bam_alignment_core.mpos;
	bam_align->flag = bam_alignment_core.flag;
	bam_align->n_cigar = bam_alignment_core.n_cigar;
	bam_align->cigar = bam_get_cigar( bam_alignment);
	bam_align->isize = bam_alignment_core.isize;
	bam_align->qual = bam_alignment_core.qual;

	if( bam_alignment_core.pos < bam_alignment_core.mpos)
		bam_align->side = LEFT;
	else
		bam_align->side = RIGHT;

	bam_align->xa = false;

	if ( params->ten_x)
		bam_align->ten_x_barcode = encode_ten_x_barcode( bam_aux_get( bam_alignment, "BX"));

	tmp = bam_aux_get( bam_alignment, "NM");
	bam_align->edit_distance = 0;
	if( tmp != 0)
		bam_align->edit_distance = bam_aux2i( tmp);


	return_type = read_mapping( in_bam->libraries[lib_index], params, ref, bam_alignment, bamToRefIndex, bam_align);

	if( bam_align != NULL)
		free_alignments( &bam_align);

	return return_type;
}
