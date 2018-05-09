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

int altLeftPrimRight = 0, primRightAltLeft = 0, altLeftAltRight = 0, primLeftAltRight = 0;
int ll=0, rr=0;

int calculateCigar( uint32_t** cigar, char altmap[])
{
	int j = 0, k;
	int cigar_count = 0;
	char num[3];
	uint32_t cigar_op, cigar_opl, cigar_opl_shifted, cigar_final;
	(*cigar) = ( uint32_t*) getMem( 20 * sizeof( uint32_t));

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
		(*cigar)[cigar_count] = cigar_final;

		cigar_count++;
	}
	//fprintf(stderr,"%s - ", altmap);
	//for(j = 0; j<cigar_count;j++)
	//fprintf(stderr,"%d%c",bam_cigar_oplen( (*new_mapping)->cigar[j]), bam_cigar_opchr( (*new_mapping)->cigar[j]));
	//fprintf(stderr,"\n");

	return cigar_count;
}


/* Pair left's alternatives with right's alternatives */
int primLeftAltRightMappings( library_properties *library, bam1_t* bam_alignment, char altmap[4][1024])
{
	int cigar_cnt, h;
	uint32_t *cigar;
	char* str = NULL;

	alternativeMapping* newEl = NULL;
	bam1_core_t bam_alignment_core = bam_alignment->core;

	set_str( &(str), bam_get_qname( bam_alignment));
	h = vh_getHash( str);
	discordantMapping* primLeft = library->mappings_discordant[h];

	//fprintf(stderr,"\n");
	/* We need to create mappings with all the alternative mappings of the left pair
	 * For these mappings, newEl->side = LEFT */
	while( primLeft != NULL)
	{
		//fprintf(stderr,"%s - %s - %d....\n", primLeft->readName, bam_alignment->data, strlen(str));
		if( strncmp( primLeft->readName, bam_alignment->data, strlen(str)) == 0)
		{
			//fprintf(stderr,"INN\n");
			/* Left-end is the first read */
			newEl = ( alternativeMapping*) getMem( sizeof( alternativeMapping));

			newEl->mapp = 0;
			int len = strlen( primLeft->readName);
			newEl->readName = ( char*) getMem( ( len + 5) * sizeof( char));
			strncpy( newEl->readName, primLeft->readName, len);

			newEl->readName[len] = '_';
			newEl->readName[len + 1] = 'A';
			newEl->readName[len + 2] = 'l';
			newEl->readName[len + 3] = 't';
			newEl->readName[len + 4] = '\0';

			/* Get the name of the chromosome */
			len = strlen( primLeft->chromosome_name);
			newEl->chromosome_name = ( char *) getMem( sizeof ( char) * len + 1);
			strcpy( newEl->chromosome_name, primLeft->chromosome_name);

			newEl->pos1 = primLeft->pos1;
			newEl->pos2 = atoi( altmap[1] + 1);
			newEl->mQual1 = primLeft->mQual1;
			newEl->mQual2 = bam_alignment_core.qual;
			newEl->editDistance = primLeft->editDistance_left + atoi( altmap[3]);

			newEl->pos1_End = primLeft->pos1_End;
			newEl->pos2_End = newEl->pos2 + library->read_length;

			newEl->flag = primLeft->flag;

			/* Check the orientation of the altmap and modify the flag */
			if( altmap[1][0] == '-')
				newEl->flag = newEl->flag | BAM_FMREVERSE;
			else if( altmap[1][0] == '+')
				newEl->flag = newEl->flag & 0xFFDF;

			/* Cigar of right read */
			cigar_cnt = calculateCigar( &cigar, altmap[2]);
			if( bam_cigar_opchr( cigar[cigar_cnt - 1]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( cigar[cigar_cnt - 1]);

			if( bam_cigar_opchr( cigar[0]) == 'S')
				newEl->pos2 = newEl->pos2 + bam_cigar_oplen( cigar[0]);

			free( cigar);

			newEl->pos2_End--;

			newEl->next = library->mappings_alternative[h];
			library->mappings_alternative[h] = newEl;

			primLeftAltRight++;
			alt_cnt_bam++;
		}
		primLeft = primLeft->next;
	}
	free( str);
	return RETURN_SUCCESS;
}

/* Pair left's alternatives with right's alternatives */
int altLeftAltRightMappings( library_properties *library, bam1_t* bam_alignment, char altmap[4][1024])
{
	int cigar_cnt, h;
	uint32_t *cigar;
	char* str = NULL;

	alternativeMapping* newEl = NULL;
	bam1_core_t bam_alignment_core = bam_alignment->core;

	set_str( &(str), bam_get_qname( bam_alignment));

	h = vh_getHash( str);
	alternativeMapping* altLeft = library->mappings_alternative[h];
	//fprintf(stderr,"\n");
	/* We need to create mappings with all the alternative mappings of the left pair
	 * For these mappings, newEl->side = LEFT */
	while( altLeft != NULL)
	{
		//fprintf(stderr,"%s - %s - %d....\n", altLeft->readName, bam_alignment->data, strlen(str));
		if( strncmp( altLeft->readName, bam_alignment->data, strlen( str)) == 0 && ( altLeft->mapp == 1))
		{
			//fprintf(stderr,"INN\n");
			/* Left-end is the first read */
			newEl = ( alternativeMapping*) getMem( sizeof( alternativeMapping));

			newEl->mapp = 0;

			newEl->readName = NULL;
			set_str( &(newEl->readName), altLeft->readName);

			/* Get the name of the chromosome */
			int len = strlen( altLeft->chromosome_name);
			newEl->chromosome_name = ( char *) getMem( sizeof ( char) * len + 1);
			strcpy( newEl->chromosome_name, altLeft->chromosome_name);

			newEl->pos1 = altLeft->pos1;
			newEl->pos2 = atoi( altmap[1] + 1);
			newEl->mQual1 = altLeft->mQual1;
			newEl->mQual2 = bam_alignment_core.qual;
			newEl->editDistance = altLeft->editDistance_left + atoi( altmap[3]);

			newEl->pos1_End = altLeft->pos1_End;
			newEl->pos2_End = newEl->pos2 + library->read_length;

			newEl->flag = altLeft->flag;

			/* Check the orientation of the altmap and modify the flag */
			if( altmap[1][0] == '-')
				newEl->flag = newEl->flag | BAM_FMREVERSE;
			else if( altmap[1][0] == '+')
				newEl->flag = newEl->flag & 0xFFDF;

			/* Cigar of right read */
			cigar_cnt = calculateCigar( &cigar, altmap[2]);
			if( bam_cigar_opchr( cigar[cigar_cnt - 1]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( cigar[cigar_cnt - 1]);

			if( bam_cigar_opchr( cigar[0]) == 'S')
				newEl->pos2 = newEl->pos2 + bam_cigar_oplen( cigar[0]);

			free( cigar);

			newEl->pos2_End--;

			newEl->next = library->mappings_alternative[h];
			library->mappings_alternative[h] = newEl;

			altLeftAltRight++;
			alt_cnt_bam++;
		}
		altLeft = altLeft->next;
	}
	free( str);
	return RETURN_SUCCESS;
}

int primRightAltLeftMappings( library_properties* library, bam1_t* bam_alignment)
{
	char *str = NULL;
	bam1_core_t bam_alignment_core = bam_alignment->core;
	uint8_t *tmp;
	uint32_t* cigar;
	int ed;

	set_str( &(str), bam_get_qname( bam_alignment));
	int h = vh_getHash( str);
	alternativeMapping *newEl = library->mappings_alternative[h];

	//fprintf(stderr,"\n");
	/* We need to create mappings with all the alternative mappings of the left pair
	 * For these mappings, newEl->side = LEFT */
	//fprintf(stderr,"\n");
	while( newEl != NULL)
	{
		if( strncmp( newEl->readName, bam_alignment->data, strlen(str)) == 0 &&
				( newEl->pos2 == bam_alignment_core.pos) && newEl->mapp == 1)
		{
			//fprintf(stderr,"%s - %s - %d....%d - %d\n", newEl->readName, bam_alignment->data, strlen(str),
			//	newEl->pos1, newEl->pos2);
			//fprintf(stderr,"INN\n");
			tmp = bam_aux_get( bam_alignment, "NM");
			ed = 0;
			if( tmp != 0)
				ed = bam_aux2i( tmp);

			newEl->editDistance = newEl->editDistance_left + ed;
			cigar = bam_get_cigar( bam_alignment);

			if( bam_cigar_opchr( cigar[bam_alignment_core.n_cigar - 1]) == 'S')
				newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( cigar[bam_alignment_core.n_cigar - 1]);

			if( bam_cigar_opchr( cigar[0]) == 'S')
				newEl->pos2 = newEl->pos2 + bam_cigar_oplen( cigar[0]);

			newEl->pos2_End--;
			newEl->mQual2 = bam_alignment_core.qual;

			primRightAltLeft++;
		}
		newEl = newEl->next;
	}
	return RETURN_SUCCESS;
}


int altLeftPrimRightMappings( bam_info* in_bam, parameters* params, int lib_index, bam1_t* bam_alignment, char altmap[4][1024])
{
	int len, cigar_cnt, chrID_bam, chrID_ref;
	uint32_t* cigar;

	alternativeMapping *newEl;
	newEl = ( alternativeMapping *) getMem( sizeof( alternativeMapping));

	bam1_core_t bam_alignment_core = bam_alignment->core;

	len = strlen( bam_get_qname( bam_alignment));
	newEl->readName = ( char*) getMem( ( len + 5) * sizeof( char));
	strncpy( newEl->readName, bam_get_qname( bam_alignment), len);

	int h = vh_getHash( newEl->readName);

	newEl->readName[len] = '_';
	newEl->readName[len + 1] = 'A';
	newEl->readName[len + 2] = 'l';
	newEl->readName[len + 3] = 't';
	newEl->readName[len + 4] = '\0';

	newEl->mapp = 1;

	/* Get the name of the chromosome */
	chrID_bam = find_chr_index_bam( altmap[0], in_bam->bam_header);
	chrID_ref = sonic_refind_chromosome_index( params->this_sonic, in_bam->bam_header->target_name[chrID_bam]);

	len = strlen( params->this_sonic->chromosome_names[chrID_ref]);
	newEl->chromosome_name = ( char *) getMem( sizeof ( char) * len + 1);
	strcpy( newEl->chromosome_name, params->this_sonic->chromosome_names[chrID_ref]);

	newEl->pos1 = atoi( altmap[1] + 1);
	newEl->pos2 = bam_alignment_core.mpos;
	newEl->mQual1 = bam_alignment_core.qual;
	newEl->mQual2 = 0;

	newEl->editDistance = atoi( altmap[3]);
	newEl->editDistance_left = atoi( altmap[3]);

	newEl->pos1_End = newEl->pos1 + in_bam->libraries[lib_index]->read_length;

	cigar_cnt = calculateCigar( &cigar, altmap[2]);

	if( bam_cigar_opchr( cigar[cigar_cnt - 1]) == 'S')
		newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( cigar[cigar_cnt - 1]);

	if( bam_cigar_opchr( cigar[0]) == 'S')
		newEl->pos1 = newEl->pos1 + bam_cigar_oplen( cigar[0]);

	newEl->pos1_End--;
	newEl->pos2_End = newEl->pos2 + in_bam->libraries[lib_index]->read_length;

	newEl->flag = bam_alignment_core.flag;

	if( altmap[1][0] == '-')
		newEl->flag = newEl->flag | BAM_FREVERSE;
	else if( altmap[1][0] == '+')
		newEl->flag = newEl->flag & 0xFFEF;


	newEl->next = in_bam->libraries[lib_index]->mappings_alternative[h];
	in_bam->libraries[lib_index]->mappings_alternative[h] = newEl;

	//fprintf(stderr,"%s - %d - %d\n", newEl->readName, newEl->pos1, newEl->pos2);
	altLeftPrimRight++;
	alt_cnt_bam++;
	free( cigar);

	return RETURN_SUCCESS;
}

/* Get the alternative mapping locations from the XA field of bam file */
void find_alt_mappings( bam_info* in_bam, parameters* params, int lib_index, bam1_t* bam_alignment, char* xa_string)
{
	int i, return_value;

	char *tok, *tok2, altmap[4][1024];
	char *str;
	char *end_str;
	bam1_core_t bam_alignment_core = bam_alignment->core;

	str = NULL;
	set_str( &str, xa_string);
	//fprintf(stderr,"\nMAIN= %s\n", str);

	tok = strtok_r( str, ";", &end_str);

	/* Each mapping in XA field ends with ; and each mapping has
	 * position(starting with -/+ for orientation, cigar and edit distance */
	while( tok != NULL)
	{
		i = 0;
		char *end_token;

		tok2 = strtok_r( tok, ",", &end_token);

		while (tok2 != NULL)
		{
			strcpy( altmap[i], tok2);
			tok2 = strtok_r( NULL, ",", &end_token);
			i++;
		}

		/* We only use the alternative mappings mapped to the same chromosome
		 * Also some alternative mappings have the same location with the primary mapping, these are eliminated*/
		if( find_chr_index_bam( altmap[0], in_bam->bam_header) == bam_alignment_core.tid)
		{
			//fprintf(stderr,"%s %s %s %s\n", altmap[0], altmap[1], altmap[2], altmap[3]);
			//&& (atoi( altmap[1] + 1) != bam_alignment_core.pos) && (atoi( altmap[1] + 1) != bam_alignment_core.pos + 1
			/* If this is the left pair */
			if( bam_alignment_core.pos < bam_alignment_core.mpos)
			{
				altLeftPrimRightMappings( in_bam, params, lib_index, bam_alignment, altmap);
			}
			else if( bam_alignment_core.pos > bam_alignment_core.mpos)
			{
				/* For alt left - alt right mappings */
				altLeftAltRightMappings( in_bam->libraries[lib_index], bam_alignment, altmap);

				/* For primary left - alt right mappings */
				primLeftAltRightMappings( in_bam->libraries[lib_index], bam_alignment, altmap);
			}
		}

		strcpy( str, tok);
		tok = strtok_r( NULL, ";", &end_str);
	}
	if( bam_alignment_core.pos > bam_alignment_core.mpos)
	{
		/* Add the edit distance of right primary mapping to left alternative mappings' */
		primRightAltLeftMappings( in_bam->libraries[lib_index], bam_alignment);
	}
}

int primary_mapping( bam_info* in_bam, parameters* params, int lib_index, bam1_t* bam_alignment, int32_t *bamToRefIndex)
{
	uint8_t *tmp;
	int return_type;
	unsigned long tenx;

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
	bam_align->qual = bam_alignment_core.qual;
	bam_align->isize = bam_alignment_core.isize;


	tmp = bam_aux_get( bam_alignment, "NM");
	bam_align->edit_distance = 0;
	if( tmp != 0)
		bam_align->edit_distance = bam_aux2i( tmp);


	return_type = read_mapping( in_bam->libraries[lib_index], params, bam_alignment, bamToRefIndex, bam_align);

	if( bam_align != NULL)
		free_alignments( &bam_align);

	return return_type;
}
