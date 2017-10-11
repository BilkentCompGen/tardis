#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vh_setcover.h"
#include "vh_buffer.h"
#include "vh_intervalhandler.h"
#include "vh_conflict.h"

#define hashtable_size 10000

int multiIndCount;//total number of individual <= totalNumInd
int sizeListClusterEl; // total number of clusters
int conflictResFlag = 1;
int weightedHeuristicFlag = 0;
int listMobileElSize;
int numCallsRequested; //maximum number of SVs the users wants to us to output
int multiLibsCount;
int cluster_count = 0;
int minimumSupNeeded;
int minimumSupNeededNoRD;
float mismatchPenalty = 0;
char multiInd[totalNumInd][strSize]; //Name of each individual
struct barcode_list_element* hashtable[hashtable_size];

struct clusters_final* clusters_all[MaxClusterCount];

struct readEl *read_names; // Array of all reads
long read_names_count;// number of total reads

//clusterElRead clusterElRead_Single; // the name of cluster which each new cluster is read to it

int free_readMappingEl( readMappingEl *ptr)
{
	if(ptr == NULL)
		return 0;
	else
	{
		free_readMappingEl( ptr->next);
		free( ptr->chromosome_name);
		free( ptr);
	}
	return 1;
}

int free_clusterIdEl( clusterIdEl *ptr)
{
	if( ptr == NULL)
		return 0;
	else
	{
		free_clusterIdEl( ptr->next);
		free( ptr);
	}
	return 1;
}

int compareReadMappingEl( const void *a, const void *b) // the compare function for qsort procudre, should sort in increasing order
{
	struct readMappingEl *arg1 = ( struct readMappingEl *)a;
	struct readMappingEl *arg2 = ( struct readMappingEl *)b;
	return ( ( *arg1).editDistance - ( *arg2).editDistance);
}


int findLibId( char *libName, char* indName)
{
	int count;
	for( count = 0; count < multiLibsCount; count++)
	{
		if( strcmp( libName, multiLibs[count].libName) == 0 && strcmp( indName, multiInd[multiLibs[count].indId]) == 0)
			return count;
	}
	fprintf(stderr, "ERROR: LIBNAME MISMATCH - libname = %s  indname = %s\n", libName, indName);
	exit( 1);
}

int findIndId( char *indName)
{
	int count;
	for( count = 0; count < multiIndCount; count++)
	{
		if( strcmp( multiInd[count], indName) == 0)
			return count;
	}
	fprintf(stderr, "ERROR - Could not find the individual id\n");
	exit( 1);
}

/* Add all the reads that support the SV */
readMappingEl *addToListOfReads(readMappingEl *linkList, readMappingEl element)
{

	readMappingEl *newEl;
	newEl = (readMappingEl *) getMem(sizeof(readMappingEl));
	newEl->readId = element.readId;

	newEl->chromosome_name = NULL;
	set_str( &(newEl->chromosome_name), element.chromosome_name);
	newEl->posMapLeft = element.posMapLeft;
	newEl->posMapRight = element.posMapRight;
	newEl->indId = element.indId;
	newEl->orient1 = element.orient1;
	newEl->orient2 = element.orient2;
	newEl->editDistance = element.editDistance;
	newEl->probEditDist = element.probEditDist;
	newEl->mapq1 = element.mapq1;
	newEl->mapq2 = element.mapq2;
	newEl->correctMappingQual = element.correctMappingQual;
	newEl->next = linkList;
	newEl->ten_x_barcode = element.ten_x_barcode;

	return newEl;
}

/* Add all the reads that support the SV */
readMappingEl *addToListOfReads2(readMappingEl *linkList, clusters_final *element)
{
	readMappingEl *newEl;
	newEl = (readMappingEl *) getMem(sizeof(readMappingEl));
	newEl->readId = element->id;

	newEl->chromosome_name = NULL;
	set_str( &(newEl->chromosome_name), element->chromosome_name1);
	newEl->posMapLeft = element->start_position;
	newEl->posMapRight = element->end_position;
	newEl->indId = findIndId(element->individual_name);
	newEl->orient1 = element->orientation_left;
	newEl->orient2 = element->orientation_right;
	newEl->editDistance = element->edit_distance;
	newEl->probEditDist = element->phred_score;
	newEl->mapq1 = element->mapping_quality_left;
	newEl->mapq2 = element->mapping_quality_right;
	//newEl->correctMappingQual = element.correctMappingQual;
	newEl->next = linkList;
	newEl->ten_x_barcode = element->ten_x_barcode;

	return newEl;
}

void calculateLikelihoodCNV(bam_info **in_bam, ref_genome* ref, parameters *params, int posS, int posE, int clusterId, double cnvProb[][10])
{
	int count, gc_val, i;
	float expectedReadCount, lambda;
	int totalReadCount, chr = -1;

	if( running_mode == SENSITIVE)
	{
		chr = sonic_refind_chromosome_index(params->this_sonic, listClusterEl[clusterId].chromosome_name);
		if( chr == -1)
		{
			fprintf(stderr, "ERROR - cannot find the chromosome in ref");
			exit( 1);
		}
	}

	for( count = 0; count < multiIndCount; count++)
	{
		expectedReadCount = 0;
		totalReadCount = 0;

		for( i = posS; i < posE; i++)
		{
			gc_val = (int) round ( sonic_get_gc_content(params->this_sonic, listClusterEl[clusterId].chromosome_name, i, i+WINDOWSLIDE));
			expectedReadCount += in_bam[count]->mean_rd_per_gc[gc_val];
			if( running_mode == QUICK)
			{
				totalReadCount += in_bam[count]->read_depth_per_chr[i];
			}
			else if( running_mode == SENSITIVE)
			{
				totalReadCount += in_bam[count]->read_depth[chr][i];

			}
		}
		listClusterEl[clusterId].readDepth[count] = ( long)totalReadCount;

		for( i = 0; i < 10; i++)
		{
			if( i > 0)
			{
				lambda = ( ( ( float)i / ( float)2 ) * expectedReadCount);
				cnvProb[count][i] =  exp( totalReadCount * log( lambda) - lambda - lgamma( totalReadCount + 1));
			}
			else if( i == 0)
			{
				if( totalReadCount < ( 0.2) * expectedReadCount)
					cnvProb[count][i] = 1;
				else
				{
					lambda = ( ( 0.01) * expectedReadCount);
					cnvProb[count][i] = exp( totalReadCount * log( lambda) - lambda - lgamma( totalReadCount + 1));
				}
			}
		}
	}
}

void calculateExpectedCN( bam_info **in_bam, ref_genome* ref, parameters *params, int posS, int posE, int clusterId, float *finalCNV)
{
	int count, pos, gc_val, i;
	float totalCNV, normalizedReadDepth;
	int chr = -1;

	if( running_mode == SENSITIVE)
	{
		chr = sonic_refind_chromosome_index(params->this_sonic, listClusterEl[clusterId].chromosome_name);
		if( chr == -1)
		{
			fprintf(stderr, "ERROR - cannot find the chromosome in ref");
			exit( 1);
		}
	}

	for( count = 0; count < multiIndCount; count++)
	{
		totalCNV = 0;
		normalizedReadDepth = 0;

		for( i = posS; i < posE; i++)
		{
			gc_val = (int) round ( sonic_get_gc_content(params->this_sonic, listClusterEl[clusterId].chromosome_name, i, i+WINDOWSLIDE));
			normalizedReadDepth += in_bam[count]->mean_rd_per_gc[gc_val];

			if( running_mode == QUICK)
			{
				totalCNV += ( float)in_bam[count]->read_depth_per_chr[i];
			}
			else if( running_mode == SENSITIVE)
			{
				totalCNV += ( float)in_bam[count]->read_depth[chr][i];
			}
		}
		finalCNV[count] = ( float)( 2 * totalCNV) / ( float)( normalizedReadDepth);
		//fprintf(stderr, "totalCNV=%f - normalizedRD=%f - final =%f - round=%f\n", totalCNV, normalizedReadDepth, finalCNV[count], round(finalCNV[count]));
	}
}



int findReadName( char *readName)
{
	int minId = 0;
	int maxId = read_names_count - 1;
	int middleId = read_names_count/2;
	if( strcmp( readName, read_names[minId].readName) == 0)
		return minId;
	if( strcmp( readName, read_names[maxId].readName) == 0)
		return maxId;

	while( strcmp( readName, read_names[middleId].readName) != 0)
	{
		if( strcmp( readName, read_names[middleId].readName) > 0)
		{
			minId = middleId;
			middleId = ( minId + maxId) / 2;
		}
		else if( strcmp( readName, read_names[middleId].readName) < 0)
		{
			maxId = middleId;
			middleId = ( minId + maxId) / 2;
		}
		if( strcmp( readName, read_names[minId].readName) == 0)
			return minId;
		if( strcmp( readName, read_names[maxId].readName) == 0)
			return maxId;

		if (minId >= maxId && strcmp( readName, read_names[minId].readName) != 0)
			return -1;
		if (minId == middleId && strcmp( readName, read_names[middleId].readName) != 0)
			return -1;
	}
	return middleId;
}


int addNewInd( char *indName)
{
	int count;

	for( count = 0; count < multiIndCount; count++)
	{
		if( strcmp( multiInd[count], indName) == 0)
			return count;
	}
	strcpy( multiInd[multiIndCount], indName);
	multiIndCount++;
	return multiIndCount - 1;
}


int freeLinkList( readMappingEl *ptr)
{
	if( ptr == NULL) return 0;
	else
	{
		freeLinkList( ptr->next);
		free( ptr);
	}
}


float scoreForEditDistanceFunc( int totalEditDist)
{
	return totalEditDist * mismatchPenalty;
}

//taken from http://www.cse.yorku.ca/~oz/hash.html
unsigned long hashing_function(unsigned char *str)
{
	unsigned long hash = 5381;
	int c;

	while (c = *str++)
		hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

	return hash;
}

/*Iterate over the cluster. Count total number of reads in cluster.
 *Count total number of reads per barcode using hashtable structure.
 *Iterate over the hashtable calculating sum for all over |B_i|^2, where |B_i| is the number of reads that has barcode B_i
 *Return 1/n^2 times the sum
 */

void init_barcode_homogeneity_hashtable(){
	int hashtable_index;
	//initialize the heads of the chains in the hashtable
	for (hashtable_index = 0; hashtable_index < hashtable_size; hashtable_index++){
		hashtable[hashtable_index] = (barcode_list_element*) malloc(sizeof(barcode_list_element));
		hashtable[hashtable_index]->ten_x_barcode = 0;
		hashtable[hashtable_index]->count = 0;
		hashtable[hashtable_index]->next = NULL;
	}
}

void free_barcode_homogeneity_hashtable(){
	int i;
	struct barcode_list_element *temp_ptr;
	for (i = 0; i < hashtable_size; i++){
		while(hashtable[i]->next){
			temp_ptr = hashtable[i]->next->next;
			free(hashtable[i]->next);
			hashtable[i]->next = temp_ptr;
		}

		free(hashtable[i]);
	}
}

double barcode_homogeneity_score(int clusterId){

	unsigned long total_read_count = 0;
	double sum = 0.0;
	double score = 1.0;
	int i;
	int hashtable_index;
	struct barcode_list_element *current_list_element;
	struct barcode_list_element *temp_ptr;
	readMappingEl *ptrReadMapping;
	int found;

	ptrReadMapping = listClusterEl[clusterId].next;
	while( ptrReadMapping != NULL){
		if( read_names[ptrReadMapping->readId].readCovered == 0
				&& listClusterEl[clusterId].indIdCount[ptrReadMapping->indId] > -1
				&& (signed long)ptrReadMapping->ten_x_barcode != -1){

			total_read_count++;

			hashtable_index = ptrReadMapping->ten_x_barcode % hashtable_size;
			current_list_element = hashtable[hashtable_index];

			found = 0;
			while (current_list_element->next){
				if (current_list_element->next->ten_x_barcode == ptrReadMapping->ten_x_barcode){
					current_list_element->next->count++;
					found = 1;
					break;
				}
				current_list_element = current_list_element->next;
			}

			if (found == 0){
				current_list_element->next = (barcode_list_element*)malloc(sizeof(barcode_list_element));
				current_list_element->next->ten_x_barcode = ptrReadMapping->ten_x_barcode;
				current_list_element->next->count = 1;
				current_list_element->next->next = NULL;
			}
		}
		ptrReadMapping=ptrReadMapping->next;
	}

	for (hashtable_index = 0; hashtable_index < hashtable_size; hashtable_index++){
		current_list_element = hashtable[hashtable_index];

		while(current_list_element->next){
			sum = sum + ((double)current_list_element->next->count*(double)current_list_element->next->count);
			current_list_element = current_list_element->next;
		}
	}

	for (i = 0; i < hashtable_size; i++){
		while(hashtable[i]->next){
			temp_ptr = hashtable[i]->next->next;
			free(hashtable[i]->next);
			hashtable[i]->next = temp_ptr;
		}
	}

	if (total_read_count == 0){
		listClusterEl[clusterId].homogeneity_score = score;
		return 1.0;
	}


	score = sum/(total_read_count*total_read_count);

	listClusterEl[clusterId].homogeneity_score = score;

	if (ten_x_flag != 1){
		return 1.0;
	}

	if (!(score > 1 || score < 0))
		return score;
	printf("Err: score isn't between 0-1: %f\n", score);

	return inf;
}


float avgEditDistReadSupportingCluster(readMappingEl *list)
{
	int countTotal = 0, editDistanceSum = 0;

	while( list != NULL)
	{
		editDistanceSum = editDistanceSum + list->editDistance;
		countTotal++;
		list=list->next;
	}
	return ( ( float)( editDistanceSum) / ( float)countTotal);
}


float avgEditDistIndReadSupportCluster(readMappingEl *list, int countInd)
{
	int countTotal = 0, editDistanceSum = 0;
	while( list != NULL)
	{
		if( list->indId == countInd)
		{
			editDistanceSum = editDistanceSum + list->editDistance;
			countTotal++;
		}
		list = list->next;
	}
	if( countTotal == 0)
		return 0;
	else
		return ( ( float)( editDistanceSum) / ( float)countTotal);
}


double* calculateTheHeuristicScore(readMappingEl *list)
{
	double *result;
	result = (double *) getMem(multiIndCount*sizeof(double));
	int count = 0;
	for( count = 0; count < multiIndCount; count++)
	{
		result[count]=0;
	}

	while( list != NULL)
	{
		result[read_names[list->readId].indId] = result[read_names[list->readId].indId] + list->probEditDist;
		list = list->next;
	}
	return result;

}


void markReadsCovered( int clusterId, int *countReads)
{
	int minDelLength = 1000000000, maxDelLength = 1000000000, totalReadRemove = 0, count;
	int readsToMark[totalNumInd];

	listClusterEl[clusterId].weight_without_homogeniety_score_at_read_covering = listClusterEl[clusterId].weight_without_homogeniety_score;
	readMappingEl *ptrRead;
	clusterIdEl * ptrClusterId;
	ptrRead=listClusterEl[clusterId].next;

	for( count = 0; count < multiIndCount; count++)
	{
		totalReadRemove = totalReadRemove + countReads[count];
		readsToMark[count] = countReads[count];
	}

	//This part will mark the reads covered and calculate the minimum and maximum length of deletion
	while( ptrRead != NULL && totalReadRemove > 0)
	{
		if( readsToMark[ptrRead->indId] > 0 && read_names[ptrRead->readId].readCovered == 0)
		{
			ptrClusterId = read_names[ptrRead->readId].next;
			listClusterEl[clusterId].indIdCount[ptrRead->indId]++;

			/* Count the number of split reads supporting this cluster */
			if( strcmp( multiLibs[read_names[ptrRead->readId].libId].libName, "SplitRead") == 0)
				listClusterEl[clusterId].sr_support[ptrRead->indId]++;

			listClusterEl[clusterId].readMappingSelected = addToListOfReads(listClusterEl[clusterId].readMappingSelected, *ptrRead);

			if( listClusterEl[clusterId].SVtype == DELETION && minDelLength > ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[read_names[ptrRead->readId].libId].maxInstSize)
				minDelLength=ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[read_names[ptrRead->readId].libId].maxInstSize;

			if( listClusterEl[clusterId].SVtype == DELETION && maxDelLength > ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[read_names[ptrRead->readId].libId].minInstSize)
				maxDelLength = ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[read_names[ptrRead->readId].libId].minInstSize;

			read_names[ptrRead->readId].readCovered = 1;
			readsToMark[ptrRead->indId]--;

			while( ptrClusterId != NULL)
			{
				listClusterEl[ptrClusterId->clusterId].oldBestIsGood = 0;
				ptrClusterId = ptrClusterId->next;
			}
			totalReadRemove--;
		}
		ptrRead = ptrRead->next;
	}

	if (listClusterEl[clusterId].SVtype == 'D')
	{
		listClusterEl[clusterId].minDelLength = minDelLength;
		listClusterEl[clusterId].maxDelLength = maxDelLength;
	}
}


void outputCluster( bam_info** in_bams, parameters* params, ref_genome* ref, int cluster_id, FILE *fpVcf)
{
	int i;
	char mobileName[strSize];

	int is_start_satellite, is_end_satellite, is_mid_satellite;
	sonic_repeat *mei_start, *mei_end, *mei_mid;
	bool MEI_Filter = false;
	struct strvar* var_example;

	for( i = 0; i < params->num_bams; i++)
		in_bams[i]->contribution = false;

	for( i = 0; i < params->num_bams; i++)
	{
		if( listClusterEl[cluster_id].indIdCount[i] > 0)
			in_bams[i]->contribution = true;
	}

	is_start_satellite = sonic_is_satellite(params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV, listClusterEl[cluster_id].posStartSV+1);
	is_end_satellite = sonic_is_satellite(params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posEndSV-1, listClusterEl[cluster_id].posEndSV);
	is_mid_satellite = sonic_is_satellite(params->this_sonic, listClusterEl[cluster_id].chromosome_name, (int)( listClusterEl[cluster_id].posStartSV + listClusterEl[cluster_id].posEndSV) / 2,
			(int)( listClusterEl[cluster_id].posStartSV + listClusterEl[cluster_id].posEndSV) / 2 + 1);

	mei_start = sonic_is_mobile_element( params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV, listClusterEl[cluster_id].posStartSV+1, params->mei);
	mei_end = sonic_is_mobile_element( params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posEndSV-1, listClusterEl[cluster_id].posEndSV, params->mei);
	mei_mid = sonic_is_mobile_element( params->this_sonic, listClusterEl[cluster_id].chromosome_name, (int)( listClusterEl[cluster_id].posStartSV + listClusterEl[cluster_id].posEndSV) / 2,
			(int)( listClusterEl[cluster_id].posStartSV + listClusterEl[cluster_id].posEndSV) / 2 + 1, params->mei);

	if( listClusterEl[cluster_id].SVtype == MEIFORWARD)
	{
		strcpy( mobileName, listClusterEl[cluster_id].mobileName);

		if( ( mei_start != NULL && strcmp( mei_start->repeat_type, mobileName) == 0 &&
				mei_start->strand == SONIC_STRAND_FWD) || ( is_start_satellite != 0))
			MEI_Filter = true;

		if( (mei_end != NULL && strcmp( mei_end->repeat_type, mobileName) == 0 &&
				mei_end->strand == SONIC_STRAND_FWD) || ( is_end_satellite != 0))
			MEI_Filter = true;

		if( ( mei_mid != NULL && strcmp( mei_mid->repeat_type, mobileName) == 0 &&
				mei_mid->strand == SONIC_STRAND_FWD == 0) || ( is_mid_satellite != 0))
			MEI_Filter = true;
	}
	else if( listClusterEl[cluster_id].SVtype == MEIREVERSE)
	{
		strcpy( mobileName, listClusterEl[cluster_id].mobileName);

		if( ( mei_start != NULL && strcmp( mei_start->repeat_type, mobileName) == 0 &&
				mei_start->strand == SONIC_STRAND_REV) || ( is_start_satellite != 0))
			MEI_Filter = true;

		if( ( mei_end != NULL && strcmp( mei_end->repeat_type, mobileName) == 0 &&
				mei_end->strand == SONIC_STRAND_REV) || ( is_end_satellite != 0))
			MEI_Filter = true;

		if( ( mei_mid != NULL && strcmp( mei_mid->repeat_type, mobileName) == 0 &&
				mei_mid->strand == SONIC_STRAND_REV) || ( is_mid_satellite != 0))
			MEI_Filter = true;
	}

	if( listClusterEl[cluster_id].SVtype == MEIFORWARD)
	{
		var_example = new_strvar( listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV_Outer, listClusterEl[cluster_id].posStartSV,
				listClusterEl[cluster_id].posEndSV_Outer, listClusterEl[cluster_id].posEndSV, listClusterEl[cluster_id].SVtype,
				avgEditDistReadSupportingCluster( listClusterEl[cluster_id].readMappingSelected), listClusterEl[cluster_id].minDelLength,
				listClusterEl[cluster_id].maxDelLength, NULL, 1.1, MEI_Filter, false, listClusterEl[cluster_id].mobileName,
				listClusterEl[cluster_id].readDepth, listClusterEl[cluster_id].CNV_Interest, listClusterEl[cluster_id].Del_Likelihood,
				listClusterEl[cluster_id].indIdCount, listClusterEl[cluster_id].sr_support, listClusterEl[cluster_id].homogeneity_score,
				listClusterEl[cluster_id].weight_without_homogeniety_score_at_read_covering);

		print_strvar( in_bams, params, var_example, fpVcf);
	}
	else if( listClusterEl[cluster_id].SVtype == MEIREVERSE)
	{
		var_example = new_strvar( listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV_Outer, listClusterEl[cluster_id].posStartSV,
				listClusterEl[cluster_id].posEndSV_Outer, listClusterEl[cluster_id].posEndSV, listClusterEl[cluster_id].SVtype,
				avgEditDistReadSupportingCluster( listClusterEl[cluster_id].readMappingSelected),listClusterEl[cluster_id].minDelLength,
				listClusterEl[cluster_id].maxDelLength, NULL, 1.1, MEI_Filter, false, listClusterEl[cluster_id].mobileName,
				listClusterEl[cluster_id].readDepth, listClusterEl[cluster_id].CNV_Interest, listClusterEl[cluster_id].Del_Likelihood,
				listClusterEl[cluster_id].indIdCount, listClusterEl[cluster_id].sr_support, listClusterEl[cluster_id].homogeneity_score,
				listClusterEl[cluster_id].weight_without_homogeniety_score_at_read_covering);

		print_strvar( in_bams, params, var_example, fpVcf);
	}
	else
	{
		var_example = new_strvar( listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV_Outer, listClusterEl[cluster_id].posStartSV,
				listClusterEl[cluster_id].posEndSV_Outer, listClusterEl[cluster_id].posEndSV, listClusterEl[cluster_id].SVtype,
				avgEditDistReadSupportingCluster( listClusterEl[cluster_id].readMappingSelected), listClusterEl[cluster_id].minDelLength,
				listClusterEl[cluster_id].maxDelLength, NULL, 1.1, listClusterEl[cluster_id].LowQual, listClusterEl[cluster_id].MEI_Del,
				listClusterEl[cluster_id].mobileName, listClusterEl[cluster_id].readDepth, listClusterEl[cluster_id].CNV_Interest,
				listClusterEl[cluster_id].Del_Likelihood, listClusterEl[cluster_id].indIdCount, listClusterEl[cluster_id].sr_support,
				listClusterEl[cluster_id].homogeneity_score, listClusterEl[cluster_id].weight_without_homogeniety_score_at_read_covering);
		print_strvar( in_bams, params, var_example, fpVcf);
	}
}

void outputPickedCluster( bam_info** in_bams, parameters* params, ref_genome* ref, FILE *fpVcf)
{
	int cluster_id, ind_id, total_rp_sup, total_sr_sup;

	for( cluster_id = 0; cluster_id < sizeListClusterEl; cluster_id++)
	{
		if (listClusterEl[cluster_id].indIdCount != NULL){

			total_rp_sup = 0;
			total_sr_sup = 0;
			for( ind_id = 0; ind_id < multiIndCount; ind_id++)
			{
				total_rp_sup = total_rp_sup + listClusterEl[cluster_id].indIdCount[ind_id];
				total_sr_sup = total_sr_sup + listClusterEl[cluster_id].sr_support[ind_id];
			}

			if( ( total_rp_sup > params->rp_threshold) || (listClusterEl[cluster_id].LowQual == true && debug_mode))
				outputCluster( in_bams, params, ref, cluster_id, fpVcf);
		}
	}
}


float calWeight( ref_genome* ref, parameters *params, int clusterId, int *countBestSetPicked)
{
	int idIndCanBePicked = 0; // the id of individuals which we are considering to have this SV in this round
	int editDistanceSum[totalNumInd]; // the total sum of edit distance for paired-end reads considered till now
	int supOfIndSeen[totalNumInd]; //the number of paired-end read supports for each individual seen in this round
	double weightedSupOfIndSeen[totalNumInd];// a normalized number of paired-end read supports for each individual seen
	int totalEditDist = 0; // total of edit distance of all paired-end reads for individuals which we are considering
	float bestScore, weightNew, normalizedWeightedSup = 0;
	int numReadCanBeCovered = 0, indCount, count, i;

	bool passMinSup;
	float bestScoreRD = 0;
	double maxRDlikelihoodDel;
	double sumLikelihoodDup, sumLikelihoodDel, sumLikelihoodNorm;
	double epsilon = 1e-5; /* was 1e-200 before */
	double scoreRDlikelihoodDel[MAX_BAMS];
	sonic_repeat *mei;

	if( listClusterEl[clusterId].oldBestIsGood == 1)
	{
		for( count = 0; count < totalNumInd; count++)
			countBestSetPicked[count] = listClusterEl[clusterId].bestReadToRemove[count];
		return listClusterEl[clusterId].oldBestScore;
	}

	for( indCount = 0; indCount < totalNumInd; indCount++)
	{
		supOfIndSeen[indCount] = 0;
		editDistanceSum[indCount] = 0;
		countBestSetPicked[indCount] = 0;
		weightedSupOfIndSeen[indCount] = 0;
	}

	bestScore = inf;
	readMappingEl *ptrReadMapping;
	ptrReadMapping = listClusterEl[clusterId].next;

	while( ptrReadMapping != NULL)
	{
		if( read_names[ptrReadMapping->readId].readCovered == 0 && listClusterEl[clusterId].indIdCount[ptrReadMapping->indId] > -1)
		{
			supOfIndSeen[ptrReadMapping->indId]++;

			//////////////////////////////////////THE HUERISTIC//////////////////////////
			//Instead of taking the total support as number of reads in each cluster, we use summation of
			//normalized number of reads. i.e each read contributes as 1/(total number of mappings) it has.

			weightedSupOfIndSeen[ptrReadMapping->indId] = weightedSupOfIndSeen[ptrReadMapping->indId] + ptrReadMapping->probEditDist;
			if( supOfIndSeen[ptrReadMapping->indId] == minimumSupNeeded)
			{
				idIndCanBePicked = idIndCanBePicked + ( int)pow( 2, ptrReadMapping->indId);
				numReadCanBeCovered = numReadCanBeCovered + supOfIndSeen[ptrReadMapping->indId];
				totalEditDist = totalEditDist + editDistanceSum[ptrReadMapping->indId];
			}
			else if( supOfIndSeen[ptrReadMapping->indId]>minimumSupNeeded)
			{
				numReadCanBeCovered++;
				totalEditDist = totalEditDist + ptrReadMapping->editDistance;
			}
		}
		ptrReadMapping = ptrReadMapping->next;
	}

	//weightNew = weightsForCombination[idIndCanBePicked]+scoreForEditDistanceFunc(totalEditDist);
	weightNew = 1;

	if( numReadCanBeCovered > 0)
	{
		for( count = 0; count < multiIndCount; count++)
		{
			if( supOfIndSeen[count] >= minimumSupNeeded)
			{
				normalizedWeightedSup = normalizedWeightedSup + supOfIndSeen[count];
			}
		}
		bestScore = ( ( float)weightNew / ( float)normalizedWeightedSup);

		mei = sonic_is_mobile_element(params->this_sonic, listClusterEl[clusterId].chromosome_name, listClusterEl[clusterId].posStartSV - 20,
				listClusterEl[clusterId].posEndSV + 20, params->mei);

		if( ( listClusterEl[clusterId].SVtype == DELETION) && mei == NULL)
		{
			listClusterEl[clusterId].MEI_Del = false;
			maxRDlikelihoodDel = 0;
			bestScoreRD = 0;
			for( i = 0; i < multiIndCount; i++)
			{
				sumLikelihoodDup = 0;
				sumLikelihoodDel = 0;
				sumLikelihoodNorm = 0;
				sumLikelihoodNorm = listClusterEl[clusterId].probabilityCNV[i][2];
				sumLikelihoodDel = listClusterEl[clusterId].probabilityCNV[i][0] + listClusterEl[clusterId].probabilityCNV[i][1];
				sumLikelihoodDup = listClusterEl[clusterId].probabilityCNV[i][3] + listClusterEl[clusterId].probabilityCNV[i][4] +
						listClusterEl[clusterId].probabilityCNV[i][5] + listClusterEl[clusterId].probabilityCNV[i][6];
				listClusterEl[clusterId].Del_Likelihood[i] = sumLikelihoodDel / (double)( sumLikelihoodDup + sumLikelihoodNorm + epsilon);

				if( supOfIndSeen[i] >= minimumSupNeeded)
				{
					scoreRDlikelihoodDel[i] = sumLikelihoodDel / ( sumLikelihoodDup + sumLikelihoodNorm + epsilon);

					if( scoreRDlikelihoodDel[i] > maxRDlikelihoodDel)
						maxRDlikelihoodDel = scoreRDlikelihoodDel[i];

					if( sumLikelihoodDel / ( sumLikelihoodDup + sumLikelihoodNorm + epsilon) > params->rd_threshold)
						bestScoreRD++;
					else
						bestScoreRD--;
				}
			}
			if( ( maxRDlikelihoodDel * multiIndCount > params->rd_threshold) && (listClusterEl[clusterId].posEndSV - listClusterEl[clusterId].posStartSV > 300))
			{
				if( bestScoreRD >= 0)
					bestScore = ( float)bestScore / ( float)( bestScoreRD + 1);
			}
			else if( maxRDlikelihoodDel * multiIndCount <= params->rd_threshold)
			{
				bestScore = inf;
				listClusterEl[clusterId].LowQual = true;
			}
		}
		else if( ( listClusterEl[clusterId].SVtype == DELETION) && mei != NULL)
		{
			listClusterEl[clusterId].mobileName = NULL;
			set_str( &(listClusterEl[clusterId].mobileName), mei->repeat_type);
			passMinSup = false;
			listClusterEl[clusterId].MEI_Del = true;
			for( i = 0; i < multiIndCount; i++)
			{
				listClusterEl[clusterId].Del_Likelihood[i] = -1;

				if( supOfIndSeen[i] > minimumSupNeededNoRD)
					passMinSup = true;
			}
			if( passMinSup == false)
				bestScore = inf;
		}
		else if( listClusterEl[clusterId].SVtype == MEIFORWARD || listClusterEl[clusterId].SVtype == MEIREVERSE)
		{
			float averageCNV_Temp = 0;
			for( i = 0; i < multiIndCount; i++)
			{
				averageCNV_Temp += listClusterEl[clusterId].CNV_Interest[i];
			}

			if( ( float)averageCNV_Temp / ( float)multiIndCount > 5)
				bestScore = inf;
		}
	}
	else
		bestScore = inf;

	for( count = 0; count < multiIndCount; count++)
	{
		if( supOfIndSeen[count] >= minimumSupNeeded)
		{
			countBestSetPicked[count] = supOfIndSeen[count];
		}
	}


	listClusterEl[clusterId].weight_without_homogeniety_score = bestScore;
	if( bestScore < inf && (ten_x_flag == 1 || output_hs_flag == 1)){
		bestScore = ( float) bestScore * barcode_homogeneity_score(clusterId);
	}

	listClusterEl[clusterId].oldBestIsGood = 1;
	listClusterEl[clusterId].oldBestScore = bestScore;

	for( count = 0; count < totalNumInd; count++)
		listClusterEl[clusterId].bestReadToRemove[count] = countBestSetPicked[count];

	return bestScore;
}


int pickSet( ref_genome* ref, parameters *params)
{
	float bestWeight;
	float newWeight;
	int countReads[totalNumInd];
	int bestReads[totalNumInd];
	int bestWeightSet;
	int count, clusterCounter;
	maxScoreInBuffer = inf - 1;

	if (ten_x_flag == 1 || output_hs_flag ==1)
		init_barcode_homogeneity_hashtable();

	while( numCallsRequested > 0)
	{
		bestWeight = inf;
		bestWeightSet = -1;

		if( !bufferIsUseful( ref, params))
		{
			emptyBuffer();
			for( clusterCounter = 0; clusterCounter < sizeListClusterEl; clusterCounter++)
			{
				newWeight = calWeight( ref, params, clusterCounter, countReads);
				addToBuffer( newWeight, clusterCounter);

				if( newWeight < bestWeight)
				{
					bestWeightSet = clusterCounter;

					for( count = 0; count < multiIndCount; count++)
						bestReads[count] = countReads[count];

					bestWeight = newWeight;
				}
			}
		}
		else
		{
			bestWeightSet = bestFromBuffer();
			bestWeight = listClusterEl[bestWeightSet].oldBestScore;

			for ( count = 0; count < multiIndCount; count++)
				bestReads[count] = listClusterEl[bestWeightSet].bestReadToRemove[count];
		}

		if( bestWeight < inf)
		{
			if( conflictResFlag == 0 || conflictsAny( bestWeightSet, bestReads) == -1)
			{
				markReadsCovered( bestWeightSet, bestReads);

				if( conflictResFlag == 1)
					addToConflict( bestWeightSet, bestReads);

				numCallsRequested--;
				listClusterEl[bestWeightSet].oldBestIsGood = 0;
			}
		}
		else
		{
			if (ten_x_flag == 1 || output_hs_flag == 1)
				free_barcode_homogeneity_hashtable();
			return 0;
		}
	}
	//freeConflict();
	if (ten_x_flag == 1 || output_hs_flag == 1)
		free_barcode_homogeneity_hashtable();
	return 1;
}

/* Each new cluster is copied into a new cell in array of clusters (listClusterEl) in index (clusterId) */
void processTheSV( bam_info **in_bams, ref_genome* ref, parameters *params, int listClusterId, int cluster_count)
{
	int posStartSV = -1, posEndSV = maxChroSize;
	int posStartSV_Outer = maxChroSize, posEndSV_Outer = -1;
	int count;
	clusters_final *tmp_cluster;

	//qsort( clusterElRead_Single.readMappingElArray, clusterElRead_Single.sizeOfCluster, sizeof( readMappingEl), compareReadMappingEl);
	tmp_cluster = clusters_all[cluster_count];

	listClusterEl[listClusterId].clusterId = listClusterId;
	listClusterEl[listClusterId].MEI_Del = false;
	listClusterEl[listClusterId].LowQual = false;
	listClusterEl[listClusterId].readMappingSelected = NULL;

	listClusterEl[listClusterId].chromosome_name = NULL;
	set_str( &(listClusterEl[listClusterId].chromosome_name), tmp_cluster->chromosome_name1);

	listClusterEl[listClusterId].mobileName = NULL;
	set_str( &(listClusterEl[listClusterId].mobileName), tmp_cluster->mei_type);

	listClusterEl[listClusterId].SVtype = tmp_cluster->SV_type;
	listClusterEl[listClusterId].next = NULL;

	listClusterEl[listClusterId].readDepth = (long *) getMem( multiIndCount * sizeof( long));
	listClusterEl[listClusterId].Del_Likelihood = (double *) getMem( multiIndCount * sizeof( double));
	listClusterEl[listClusterId].CNV_Interest = (float *) getMem( multiIndCount * sizeof( float));
	listClusterEl[listClusterId].indIdCount = (int *) getMem( multiIndCount * sizeof( int));
	listClusterEl[listClusterId].sr_support = (int *) getMem( multiIndCount * sizeof( int));

	for( count = 0; count < multiIndCount; count++)
	{
		listClusterEl[listClusterId].Del_Likelihood[count] = 0;
		listClusterEl[listClusterId].indIdCount[count] = 0;
		listClusterEl[listClusterId].sr_support[count] = 0;
		listClusterEl[listClusterId].readDepth[count] = 0;
		listClusterEl[listClusterId].CNV_Interest[count] = 0;
	}

	if( listClusterEl[listClusterId].SVtype == MEIFORWARD)
	{
		while( tmp_cluster != NULL)
		{
			if( tmp_cluster->orientation_left == 'F')
			{
				if( tmp_cluster->start_position > posStartSV)
					posStartSV = tmp_cluster->start_position;
				if( tmp_cluster->start_position < posStartSV_Outer)
					posStartSV_Outer = tmp_cluster->start_position;
			}
			if( tmp_cluster->orientation_left == 'R')
			{
				if( tmp_cluster->start_position < posEndSV)
					posEndSV = tmp_cluster->start_position;
				if( tmp_cluster->start_position > posEndSV_Outer)
					posEndSV_Outer = tmp_cluster->start_position;
			}
			listClusterEl[listClusterId].next = addToListOfReads2( listClusterEl[listClusterId].next, tmp_cluster);
			tmp_cluster = tmp_cluster->next;
		}
	}
	else if( listClusterEl[listClusterId].SVtype == MEIREVERSE)
	{
		while( tmp_cluster != NULL)
		{
			if (tmp_cluster->orientation_left == 'F')
			{
				if( tmp_cluster->start_position > posStartSV)
					posStartSV = tmp_cluster->start_position;
				if( tmp_cluster->start_position < posStartSV_Outer)
					posStartSV_Outer = tmp_cluster->start_position;
			}
			if( tmp_cluster->orientation_left == 'R')
			{
				if( tmp_cluster->start_position < posEndSV)
					posEndSV = tmp_cluster->start_position;
				if( tmp_cluster->start_position > posEndSV_Outer)
					posEndSV_Outer = tmp_cluster->start_position;
			}
			listClusterEl[listClusterId].next = addToListOfReads2( listClusterEl[listClusterId].next, tmp_cluster);
			tmp_cluster = tmp_cluster->next;
		}
	}
	else
	{
		while( tmp_cluster != NULL)
		{
			if( tmp_cluster->start_position > posStartSV)
				posStartSV = tmp_cluster->start_position;
			if( tmp_cluster->start_position < posStartSV_Outer)
				posStartSV_Outer = tmp_cluster->start_position;
			if( tmp_cluster->end_position < posEndSV)
				posEndSV = tmp_cluster->end_position;
			if( tmp_cluster->end_position > posEndSV_Outer)
				posEndSV_Outer = tmp_cluster->end_position;

			listClusterEl[listClusterId].next = addToListOfReads2( listClusterEl[listClusterId].next, tmp_cluster);
			tmp_cluster = tmp_cluster->next;
		}
	}
	listClusterEl[listClusterId].posStartSV = posStartSV;
	listClusterEl[listClusterId].posStartSV_Outer = posStartSV_Outer;
	listClusterEl[listClusterId].posEndSV = posEndSV;
	listClusterEl[listClusterId].posEndSV_Outer = posEndSV_Outer;

	if( listClusterEl[listClusterId].SVtype == DELETION)
	{
		calculateExpectedCN( in_bams, ref, params, listClusterEl[listClusterId].posStartSV, listClusterEl[listClusterId].posEndSV, listClusterId, listClusterEl[listClusterId].CNV_Interest);
		calculateLikelihoodCNV( in_bams, ref, params, listClusterEl[listClusterId].posStartSV, listClusterEl[listClusterId].posEndSV, listClusterId, listClusterEl[listClusterId].probabilityCNV);
	}
	else if( listClusterEl[listClusterId].SVtype == MEIFORWARD || listClusterEl[listClusterId].SVtype == MEIREVERSE)
	{
		calculateExpectedCN( in_bams, ref, params, listClusterEl[listClusterId].posStartSV_Outer, listClusterEl[listClusterId].posEndSV_Outer, listClusterId, listClusterEl[listClusterId].CNV_Interest);
		calculateLikelihoodCNV( in_bams, ref, params, listClusterEl[listClusterId].posStartSV, listClusterEl[listClusterId].posEndSV, listClusterId, listClusterEl[listClusterId].probabilityCNV);
	}
}


void init(bam_info **in_bams, parameters *params, ref_genome* ref)
{
	int i, j, count, multiIndId, return_value, listClusterElId, readId, cluster_size = 0, len;

	clusterIdEl *clusterIdElNew;
	readMappingEl *ptrMapping;
	clusters_final *tmp_cluster;

	multiIndCount = 0;

	for( i = 0; i < params->num_bams; i++)
		multiLibsCount += in_bams[i]->num_libraries;

	multiLibs = ( multiLib *) getMem( sizeof( multiLib) * ( multiLibsCount + params->num_bams + 1));
	multiLibsCount = 0;

	for( j = 0; j < params->num_bams; j++)
	{
		for( i = 0; i < in_bams[j]->num_libraries; i++)
		{
			multiLibs[multiLibsCount].libName = NULL;
			set_str( &(multiLibs[multiLibsCount].libName), in_bams[j]->libraries[i]->libname);
			multiLibs[multiLibsCount].indId = addNewInd( in_bams[j]->sample_name);
			multiLibs[multiLibsCount].maxInstSize = in_bams[j]->libraries[i]->conc_max - 2 * in_bams[j]->libraries[i]->read_length;
			multiLibs[multiLibsCount].minInstSize = in_bams[j]->libraries[i]->conc_min - 2 * in_bams[j]->libraries[i]->read_length;
			if( multiLibs[multiLibsCount].maxInstSize < 0)
				multiLibs[multiLibsCount].maxInstSize = 0;
			if( multiLibs[multiLibsCount].minInstSize < 0)
				multiLibs[multiLibsCount].minInstSize = 0;
			multiLibs[multiLibsCount].readLen = in_bams[j]->libraries[i]->read_length;
			multiLibsCount++;
		}
		if( running_mode == QUICK && !params->no_soft_clip)
		{
			multiLibs[multiLibsCount].libName = NULL;
			set_str( &multiLibs[multiLibsCount].libName, "SplitRead");
			multiLibs[multiLibsCount].indId = addNewInd( in_bams[j]->sample_name);
			multiLibs[multiLibsCount].maxInstSize =  0;
			multiLibs[multiLibsCount].minInstSize = 0;
			multiLibs[multiLibsCount].readLen = 0;
			multiLibsCount++;
		}
	}
	indCount = multiIndCount;

	listClusterEl = (clusterEl *) getMem( (cluster_count + 1) * sizeof( clusterEl));

	for( count = 0; count < cluster_count; count++)
	{
		listClusterEl[count].clusterId = 0;
		listClusterEl[count].oldBestIsGood = -1;
		listClusterEl[count].next = NULL;
		listClusterEl[count].homogeneity_score = 1.0;
	}
	listClusterElId = 0;

	for( count = 0; count < cluster_count; count++)
	{
		tmp_cluster = clusters_all[count];

		/* Process the SVs in the cluster */
		while( tmp_cluster != NULL)
		{
			/* Adding the cluster to the read */
			readId = findReadName( tmp_cluster->read_name);
			if (readId == -1)
				continue;
			tmp_cluster->id = readId;

			read_names[readId].indId = findIndId( tmp_cluster->individual_name);
			read_names[readId].libId = findLibId( tmp_cluster->library_name, tmp_cluster->individual_name);

			/* Creating a cluster identification element to be added to the list of clusters that this read is in */
			clusterIdElNew = (clusterIdEl *) getMem( sizeof( clusterIdEl));
			clusterIdElNew->clusterId = listClusterElId;
			clusterIdElNew->next = read_names[readId].next;
			read_names[readId].next = clusterIdElNew;

			if( weightedHeuristicFlag != 1)
				tmp_cluster->phred_score = 1;

			cluster_size++;
			tmp_cluster = tmp_cluster->next;
		}
		/* If all the SVs are processed in the cluster */
		if (cluster_size > 0)
			processTheSV( in_bams, ref, params, listClusterElId, count);

		len = listClusterEl[listClusterElId].posEndSV - listClusterEl[listClusterElId].posStartSV;

		/* Check the new cluster's SV length. If bigger than the threshold, prune out */
		if( ( len > maxLengthSV_Del && listClusterEl[listClusterElId].SVtype == DELETION )
				|| ( len > maxLengthSV_Inv && listClusterEl[listClusterElId].SVtype == INVERSION )
				|| ( len > maxLengthSV_TDup && listClusterEl[listClusterElId].SVtype == TANDEMDUP ))
		{
			ptrMapping = listClusterEl[listClusterElId].next;

			while( ptrMapping != NULL)
			{
				read_names[ptrMapping->readId].next = read_names[ptrMapping->readId].next->next;
				ptrMapping = ptrMapping->next;
			}
			ptrMapping = listClusterEl[listClusterElId].next;
			freeLinkList( ptrMapping);
		}
		else
			listClusterElId++;

		listClusterEl[listClusterElId].clusterId = listClusterElId;
		cluster_size = 0;
	}
	sizeListClusterEl = listClusterElId;
}

void free_clusters()
{
	int i;
	clusters_final* tmp_cluster, *tmp_cluster_next;

	for( i = 0; i < cluster_count; i++)
	{
		tmp_cluster = clusters_all[i];
		while( tmp_cluster != NULL)
		{
			tmp_cluster_next = tmp_cluster->next;

			free( tmp_cluster->chromosome_name1);
			free( tmp_cluster->chromosome_name2);
			free( tmp_cluster->individual_name);
			free( tmp_cluster->library_name);
			free( tmp_cluster->mei_subclass);
			free( tmp_cluster->mei_type);
			free( tmp_cluster->read_name);
			free( tmp_cluster);
			tmp_cluster = tmp_cluster_next;
		}
	}
	cluster_count = 0;
}

void vh_setcover( bam_info **in_bams, parameters *params, ref_genome* ref, FILE *fpVcf)
{
	int i;
	//fprintf( stderr,"Inside vh_setcover: 10x flag: %d\n", ten_x_flag);

	numCallsRequested = maxNumSV;

	sv_count = 0;
	sv_lowqual_count = 0;
	minimumSupNeeded = params->rp_threshold;
	minimumSupNeededNoRD = params->rp_threshold;

	init( in_bams, params, ref);
	free_clusters(); /* The clusters created in clustering step are freed */
	pickSet( ref, params);
	outputPickedCluster( in_bams, params, ref, fpVcf);

	if( running_mode == QUICK && debug_mode)
		fprintf( stderr, "There are %d SVs and %d LowQual\n", sv_count, sv_lowqual_count);
	else
		fprintf( stderr, "There are %d SVs\n", sv_count);

	/* Free the structures, etc */
	if( listClusterEl != NULL)
	{
		for( i = 0; i < sizeListClusterEl; i++)
		{
			free( listClusterEl[i].mobileName);
			free( listClusterEl[i].chromosome_name);
			free( listClusterEl[i].indIdCount);
			free( listClusterEl[i].readDepth);
			free( listClusterEl[i].CNV_Interest);
			free( listClusterEl[i].Del_Likelihood);
			free_readMappingEl( listClusterEl[i].next);
			free_readMappingEl( listClusterEl[i].readMappingSelected);
		}
		free( listClusterEl);
	}
	free( multiLibs->libName);
	free( multiLibs);
	multiLibs = NULL;

	emptyBuffer();

	for( i = 0; i < read_names_count; i++)
	{
		free( read_names[i].readName);
		free_clusterIdEl( read_names[i].next);
	}
	free( read_names);
	read_names = NULL;
}
