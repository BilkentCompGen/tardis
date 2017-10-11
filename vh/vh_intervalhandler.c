/*
 * vh_intervalhandler.c
 *
 *  Created on: Nov 1, 2016
 *      Author: tardis
 */

#define _GNU_SOURCE
#include "vh_intervalhandler.h"
#include "stdio.h"
#include "string.h"
#include "vh_common.h"


int mei_filtering( ref_genome* ref, parameters *params)
{
	int i, mei_count = 0, chr = 0, indL, indR;
	char *svtype;
	char* chromosome_name = NULL;

	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;

	sonic_repeat *mei_left, *mei_right;

	libInfo = g_libInfo;

	while ( libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while ( divetReadMappingPtr != NULL)
		{
			/* Check the mei table whether the right or left read is inside a mobile element */
			mei_left  = sonic_is_mobile_element( params->this_sonic, divetReadMappingPtr->chromosome_name, divetReadMappingPtr->locMapLeftStart, divetReadMappingPtr->locMapLeftEnd, params->mei);
			mei_right = sonic_is_mobile_element( params->this_sonic, divetReadMappingPtr->chromosome_name, divetReadMappingPtr->locMapRightStart, divetReadMappingPtr->locMapRightEnd, params->mei);

			/* Horrible. Fix. There is no meaning of the --mei parameter with this type of bad coding. Everything is hard-coded, parameter becomes useless */
			/* At the very least, create int mei_code(char *class, int upper_lower_case) and use it. */
			/* I had to convert RepeatMasker classification, i.e. SINE/Alu to your style. Still bad. See above how to fix this. make meiType an integer and #define codes */

			if( mei_left != NULL)
			{
				/* Horrible. Replace svType with #define */
				divetReadMappingPtr->svType = MEILEFTPAIR;
				/* NOTE: repeat_class in SONIC is like SINE/Alu;  not Alu */
				//set_str( &(divetReadMappingPtr->meiType), mei_left->repeat_class);

				/* Horrible. See above
				if ( !strcmp(divetReadMappingPtr->meiType, "SINE/Alu")){
					free(divetReadMappingPtr->meiType);
					divetReadMappingPtr->meiType = NULL;
					set_str (&(divetReadMappingPtr->meiType), "Alu");
				}
				else if ( !strcmp(divetReadMappingPtr->meiType, "LINE/L1")){
					free(divetReadMappingPtr->meiType);
					divetReadMappingPtr->meiType = NULL;
					set_str (&(divetReadMappingPtr->meiType), "L1");
				}
				else if ( !strcmp(divetReadMappingPtr->meiType, "Other")){
					free(divetReadMappingPtr->meiType);
					divetReadMappingPtr->meiType = NULL;
					set_str (&(divetReadMappingPtr->meiType), "SVA");
				}
				 */

				divetReadMappingPtr->mei_subclass = NULL;
				set_str( &(divetReadMappingPtr->mei_subclass), mei_left->repeat_type);

				divetReadMappingPtr->meiType = NULL;
				set_str( &(divetReadMappingPtr->meiType), mei_left->repeat_class);

				/*IF THE MEI INSERT IS + STRAND UPPERCASE IF - STRAND LOWER CASE*/
				if ( ( divetReadMappingPtr->orientationLeft == 'F' && mei_left->strand == SONIC_STRAND_REV)
						|| ( divetReadMappingPtr->orientationLeft == 'R' && mei_left->strand == SONIC_STRAND_FWD  ))
				{
					/* Horrible. Change this to the mei_code as described above */
					//divetReadMappingPtr->meiType[0] = tolower(divetReadMappingPtr->meiType[0]);
					divetReadMappingPtr->mei_code = (mei_left->mei_code * 2) + 1;
				}
				else
					divetReadMappingPtr->mei_code = mei_left->mei_code * 2;
				mei_count++;
			}
			if( mei_right != NULL)
			{
				divetReadMappingPtr->svType = MEIRIGHTPAIR;
				/*
				divetReadMappingPtr->meiType = NULL;
				set_str( &(divetReadMappingPtr->meiType), mei_right->repeat_class);

				/* Horrible. See above
				if ( !strcmp(divetReadMappingPtr->meiType, "SINE/Alu")){
					free(divetReadMappingPtr->meiType);
					divetReadMappingPtr->meiType = NULL;
					set_str (&(divetReadMappingPtr->meiType), "Alu");
				}
				else if ( !strcmp(divetReadMappingPtr->meiType, "LINE/L1")){
					free(divetReadMappingPtr->meiType);
					divetReadMappingPtr->meiType = NULL;
					set_str (&(divetReadMappingPtr->meiType), "L1");
				}
				else if ( !strcmp(divetReadMappingPtr->meiType, "Other")){
					free(divetReadMappingPtr->meiType);
					divetReadMappingPtr->meiType = NULL;
					set_str (&(divetReadMappingPtr->meiType), "SVA");
				}
				 */

				divetReadMappingPtr->mei_subclass = NULL;
				set_str( &(divetReadMappingPtr->mei_subclass), mei_right->repeat_type);

				divetReadMappingPtr->meiType = NULL;
				set_str( &(divetReadMappingPtr->meiType), mei_right->repeat_class);

				/*IF THE MEI INSERT IS + STRAND UPPERCASE IF - STRAND LOWER CASE*/
				if ( ( divetReadMappingPtr->orientationRight == 'F' && mei_right->strand == SONIC_STRAND_REV)
						|| ( divetReadMappingPtr->orientationRight == 'R' && mei_right->strand == SONIC_STRAND_FWD ))
				{
					/* Horrible. Change this to the mei_code as described above */
					//divetReadMappingPtr->meiType[0] = tolower(divetReadMappingPtr->meiType[0]);
					divetReadMappingPtr->mei_code = (mei_right->mei_code * 2) + 1;
				}
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
