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

struct readEl *listReadEl; // Array of all reads
int sizeListReadEl;// number of total reads

clusterElRead clusterElRead_Single; // the name of cluster which each new cluster is read to it

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

/* Each new cluster read (clusterElRead_Single) is copied into a new cell in array of clusters (listClusterEl) in index (clusterId) */
void processTheSV( bam_info **in_bams, ref_genome* ref, parameters *params, int clusterId)
{
	int posStartSV = -1, posEndSV = maxChroSize;
	int posStartSV_Outer = maxChroSize, posEndSV_Outer = -1;
	int count;

	qsort( clusterElRead_Single.readMappingElArray, clusterElRead_Single.sizeOfCluster, sizeof( readMappingEl), compareReadMappingEl);

	listClusterEl[clusterId].clusterId = clusterId;
	listClusterEl[clusterId].MEI_Del = false;
	listClusterEl[clusterId].LowQual = false;
	listClusterEl[clusterId].readMappingSelected = NULL;

	listClusterEl[clusterId].chromosome_name = NULL;
	set_str( &(listClusterEl[clusterId].chromosome_name), clusterElRead_Single.readMappingElArray[0].chromosome_name);
	listClusterEl[clusterId].next = NULL;

	listClusterEl[clusterId].mobileName = NULL;

	listClusterEl[clusterId].readDepth = (long *) getMem( multiIndCount * sizeof( long));
	listClusterEl[clusterId].Del_Likelihood = (double *) getMem( multiIndCount * sizeof( double));
	listClusterEl[clusterId].CNV_Interest = (float *) getMem( multiIndCount * sizeof( float));
	listClusterEl[clusterId].indIdCount = (int *) getMem( multiIndCount * sizeof( int));
	listClusterEl[clusterId].sr_support = (int *) getMem( multiIndCount * sizeof( int));

	for( count = 0; count < multiIndCount; count++)
	{
		listClusterEl[clusterId].Del_Likelihood[count] = 0;
		listClusterEl[clusterId].indIdCount[count] = 0;
		listClusterEl[clusterId].sr_support[count] = 0;
		listClusterEl[clusterId].readDepth[count] = 0;
		listClusterEl[clusterId].CNV_Interest[count] = 0;
	}

	if( listClusterEl[clusterId].SVtype == 'A')
	{
		for( count = 0; count < clusterElRead_Single.sizeOfCluster; count++)
		{
			//Move each paired-end read in the cluster (clusterElRead_Single) into the array of clusters listClusterEl
			if( clusterElRead_Single.readMappingElArray[count].orient1 == 'F')
			{
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft > posStartSV)
					posStartSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft < posStartSV_Outer)
					posStartSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
			}
			if( clusterElRead_Single.readMappingElArray[count].orient1 == 'R')
			{
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft < posEndSV)
					posEndSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft > posEndSV_Outer)
					posEndSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
			}
			listClusterEl[clusterId].next = addToListOfReads( listClusterEl[clusterId].next, clusterElRead_Single.readMappingElArray[count]);
		}

		listClusterEl[clusterId].posStartSV = posStartSV;
		listClusterEl[clusterId].posStartSV_Outer = posStartSV_Outer;
		listClusterEl[clusterId].posEndSV = posEndSV;
		listClusterEl[clusterId].posEndSV_Outer = posEndSV_Outer;
	}
	else if( listClusterEl[clusterId].SVtype == 'B')
	{
		for( count = 0; count < clusterElRead_Single.sizeOfCluster; count++)
		{
			//Move each paired-end read in the cluster (clusterElRead_Single) into the array of clusters listClusterEl
			if (clusterElRead_Single.readMappingElArray[count].orient1 == 'F')
			{
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft > posStartSV)
					posStartSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft < posStartSV_Outer)
					posStartSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
			}
			if( clusterElRead_Single.readMappingElArray[count].orient1 == 'R')
			{
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft < posEndSV)
					posEndSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
				if( clusterElRead_Single.readMappingElArray[count].posMapLeft > posEndSV_Outer)
					posEndSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
			}
			listClusterEl[clusterId].next = addToListOfReads( listClusterEl[clusterId].next, clusterElRead_Single.readMappingElArray[count]);
		}
		listClusterEl[clusterId].posStartSV = posStartSV;
		listClusterEl[clusterId].posStartSV_Outer = posStartSV_Outer;
		listClusterEl[clusterId].posEndSV = posEndSV;
		listClusterEl[clusterId].posEndSV_Outer = posEndSV_Outer;
	}
	else
	{
		for( count = 0; count < clusterElRead_Single.sizeOfCluster; count++)
		{
			//Move each paired-end read in the cluster (clusterElRead_Single) into the array of clusters listClusterEl
			if( clusterElRead_Single.readMappingElArray[count].posMapLeft > posStartSV)
				posStartSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
			if( clusterElRead_Single.readMappingElArray[count].posMapLeft < posStartSV_Outer)
				posStartSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
			if( clusterElRead_Single.readMappingElArray[count].posMapRight < posEndSV)
				posEndSV = clusterElRead_Single.readMappingElArray[count].posMapRight;
			if( clusterElRead_Single.readMappingElArray[count].posMapRight > posEndSV_Outer)
				posEndSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapRight;
			listClusterEl[clusterId].next = addToListOfReads( listClusterEl[clusterId].next, clusterElRead_Single.readMappingElArray[count]);
		}
		listClusterEl[clusterId].posStartSV = posStartSV;
		listClusterEl[clusterId].posStartSV_Outer = posStartSV_Outer;
		listClusterEl[clusterId].posEndSV = posEndSV;
		listClusterEl[clusterId].posEndSV_Outer = posEndSV_Outer;
	}

	if( listClusterEl[clusterId].SVtype == 'D')
	{
		calculateExpectedCN( in_bams, ref, params, listClusterEl[clusterId].posStartSV, listClusterEl[clusterId].posEndSV, clusterId, listClusterEl[clusterId].CNV_Interest);
		calculateLikelihoodCNV( in_bams, ref, params, listClusterEl[clusterId].posStartSV, listClusterEl[clusterId].posEndSV, clusterId, listClusterEl[clusterId].probabilityCNV);
	}
	else if( listClusterEl[clusterId].SVtype == 'A' || listClusterEl[clusterId].SVtype == 'B')
	{
		calculateExpectedCN( in_bams, ref, params, listClusterEl[clusterId].posStartSV_Outer, listClusterEl[clusterId].posEndSV_Outer, clusterId, listClusterEl[clusterId].CNV_Interest);
		calculateLikelihoodCNV( in_bams, ref, params, listClusterEl[clusterId].posStartSV, listClusterEl[clusterId].posEndSV, clusterId, listClusterEl[clusterId].probabilityCNV);
	}
}

int findReadName( char *readName)
{
	int minId = 0;
	int maxId = sizeListReadEl - 1;
	int middleId = sizeListReadEl/2;
	if( strcmp( readName, listReadEl[minId].readName) == 0)
		return minId;
	if( strcmp( readName, listReadEl[maxId].readName) == 0)
		return maxId;

	while( strcmp( readName, listReadEl[middleId].readName) != 0)
	{
		if( strcmp( readName, listReadEl[middleId].readName) > 0)
		{
			minId = middleId;
			middleId = ( minId + maxId) / 2;
		}
		else if( strcmp( readName, listReadEl[middleId].readName) < 0)
		{
			maxId = middleId;
			middleId = ( minId + maxId) / 2;
		}
		if( strcmp( readName, listReadEl[minId].readName) == 0)
			return minId;
		if( strcmp( readName, listReadEl[maxId].readName) == 0)
			return maxId;

		if (minId >= maxId && strcmp( readName, listReadEl[minId].readName) != 0)
			return -1;
		if (minId == middleId && strcmp( readName, listReadEl[middleId].readName) != 0)
			return -1;
	}
	return middleId;
}


int findIndId( int multiLibsCount, char *libName)
{
	int count;
	for( count = 0; count < multiLibsCount; count++)
	{
		if( strcmp( multiLibs[count].libName, libName) == 0)
			return multiLibs[count].indId;
	}
}

int findIndId2( char *indName)
{
	int count;
	for( count = 0; count < multiIndCount; count++)
	{
		if( strcmp( multiInd[count], indName) == 0)
			return count;
	}
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
		if( listReadEl[ptrReadMapping->readId].readCovered == 0
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
	int countTotal=0, editDistanceSum=0;

	while(list!=NULL)
	{
		editDistanceSum = editDistanceSum + list->editDistance;
		countTotal++;
		list=list->next;
	}
	return ((float)(editDistanceSum)/(float)countTotal);
}


float avgEditDistIndReadSupportCluster(readMappingEl *list, int countInd)
{
	int countTotal=0, editDistanceSum=0;
	while(list!=NULL)
	{
		if (list->indId==countInd)
		{
			editDistanceSum = editDistanceSum + list->editDistance;
			countTotal++;
		}
		list=list->next;
	}
	if (countTotal==0)
		return 0; else
			return ((float)(editDistanceSum)/(float)countTotal);
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
		result[listReadEl[list->readId].indId] = result[listReadEl[list->readId].indId] + list->probEditDist;
		list = list->next;
	}
	return result;

}


void markReadsCovered(int clusterId, int *countReads)
{
	int minDelLength = 1000000000, maxDelLength = 1000000000;

	readMappingEl *ptrRead;
	clusterIdEl * ptrClusterId;
	ptrRead=listClusterEl[clusterId].next;
	int totalReadRemove=0;
	int readsToMark[totalNumInd];
	int count;
	for ( count = 0; count < multiIndCount; count++)
	{
		totalReadRemove=totalReadRemove+countReads[count];
		readsToMark[count]=countReads[count];
	}

	//This part will mark the reads covered and calculate the minimum and maximum length of deletion
	while(ptrRead!=NULL && totalReadRemove>0)
	{
		if (readsToMark[ptrRead->indId]>0 && listReadEl[ptrRead->readId].readCovered==0)
		{
			ptrClusterId=listReadEl[ptrRead->readId].next;
			listClusterEl[clusterId].indIdCount[ptrRead->indId]++;

			/* Count the number of split reads supporting this cluster */
			if( strcmp( multiLibs[listReadEl[ptrRead->readId].libId].libName, "SplitRead") == 0)
			{
				listClusterEl[clusterId].sr_support[ptrRead->indId]++;
			}

			listClusterEl[clusterId].readMappingSelected=addToListOfReads(listClusterEl[clusterId].readMappingSelected, *ptrRead);

			if (listClusterEl[clusterId].SVtype=='D' && minDelLength > ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].maxInstSize)
			{
				minDelLength=ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].maxInstSize;
			}

			if (listClusterEl[clusterId].SVtype=='D' && maxDelLength > ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].minInstSize)
			{
				maxDelLength=ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].minInstSize;
			}

			listReadEl[ptrRead->readId].readCovered=1;
			readsToMark[ptrRead->indId]--;
			while (ptrClusterId!=NULL)
			{
				listClusterEl[ptrClusterId->clusterId].oldBestIsGood=0;
				ptrClusterId=ptrClusterId->next;
			}
			totalReadRemove--;
		}
		ptrRead=ptrRead->next;
	}

	if (listClusterEl[clusterId].SVtype=='D')
	{
		listClusterEl[clusterId].minDelLength=minDelLength;
		listClusterEl[clusterId].maxDelLength=maxDelLength;
	}
}


void outputCluster( bam_info** in_bams, parameters* params, ref_genome* ref, int set, FILE *fpVcf)
{
	int i;
	char mobileName[strSize];

	int is_start_satellite, is_end_satellite, is_mid_satellite;
	sonic_repeat *mei_start, *mei_end, *mei_mid;
	bool MEI_Filter = false;
	struct strvar* var_example;

	for( i = 0; i < params->num_bams; i++)
	{
		in_bams[i]->contribution = false;
	}

	for( i = 0; i < params->num_bams; i++)
	{
		if( listClusterEl[set].indIdCount[i] > 0)
			in_bams[i]->contribution = true;
	}

	is_start_satellite = sonic_is_satellite(params->this_sonic, listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV, listClusterEl[set].posStartSV+1);
	is_end_satellite   = sonic_is_satellite(params->this_sonic, listClusterEl[set].chromosome_name, listClusterEl[set].posEndSV-1, listClusterEl[set].posEndSV);
	is_mid_satellite   = sonic_is_satellite(params->this_sonic, listClusterEl[set].chromosome_name, (int)( listClusterEl[set].posStartSV + listClusterEl[set].posEndSV) / 2,
			(int)( listClusterEl[set].posStartSV + listClusterEl[set].posEndSV) / 2 + 1);

	mei_start = sonic_is_mobile_element( params->this_sonic, listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV, listClusterEl[set].posStartSV+1, params->mei);
	mei_end   = sonic_is_mobile_element( params->this_sonic, listClusterEl[set].chromosome_name, listClusterEl[set].posEndSV-1, listClusterEl[set].posEndSV, params->mei);
	mei_mid   = sonic_is_mobile_element( params->this_sonic, listClusterEl[set].chromosome_name, (int)( listClusterEl[set].posStartSV + listClusterEl[set].posEndSV) / 2,
			(int)( listClusterEl[set].posStartSV + listClusterEl[set].posEndSV) / 2 + 1, params->mei);

	/* Horrible. Replace all SVtypes with #defines */
	if( listClusterEl[set].SVtype == 'A')
	{
		strcpy( mobileName, listClusterEl[set].mobileName);
		//idMobile = search_MobileName( ref, listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV);
		//idSimpleRep = search_SimpleRepeat( listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV);
		if( ( mei_start != NULL && strcmp( mei_start->repeat_type, mobileName) == 0 &&
				mei_start->strand == SONIC_STRAND_FWD) || ( is_start_satellite != 0))
		{
			MEI_Filter = true;
		}
		//idMobile = search_MobileName( ref, listClusterEl[set].chromosome_name, listClusterEl[set].posEndSV);
		//idSimpleRep = search_SimpleRepeat( listClusterEl[set].chromosome_name, listClusterEl[set].posEndSV);
		if( (mei_end != NULL && strcmp( mei_end->repeat_type, mobileName) == 0 &&
				mei_end->strand == SONIC_STRAND_FWD) || ( is_end_satellite != 0))
		{
			MEI_Filter = true;
		}
		//idMobile = search_MobileName( ref, listClusterEl[set].chromosome_name, ( int)( listClusterEl[set].posStartSV + listClusterEl[set].posEndSV) / 2);
		//idSimpleRep = search_SimpleRepeat( listClusterEl[set].chromosome_name, ( int)( listClusterEl[set].posStartSV + listClusterEl[set].posEndSV) / 2);
		if( ( mei_mid != NULL && strcmp( mei_mid->repeat_type, mobileName) == 0 &&
				mei_mid->strand == SONIC_STRAND_FWD == 0) || ( is_mid_satellite != 0))
		{
			MEI_Filter = true;
		}
	}
	else if( listClusterEl[set].SVtype == 'B')
	{
		strcpy( mobileName, listClusterEl[set].mobileName);
		//idMobile = search_MobileName( ref, listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV);
		//idSimpleRep = search_SimpleRepeat( listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV);
		if( ( mei_start != NULL && strcmp( mei_start->repeat_type, mobileName) == 0 &&
				mei_start->strand == SONIC_STRAND_REV) || ( is_start_satellite != 0))
		{
			MEI_Filter = true;
		}
		//idMobile = search_MobileName( ref, listClusterEl[set].chromosome_name, listClusterEl[set].posEndSV);
		//idSimpleRep = search_SimpleRepeat( listClusterEl[set].chromosome_name, listClusterEl[set].posEndSV);
		if( ( mei_end != NULL && strcmp( mei_end->repeat_type, mobileName) == 0 &&
				mei_end->strand == SONIC_STRAND_REV) || ( is_end_satellite != 0))
		{
			MEI_Filter = true;
		}
		//idMobile = search_MobileName( ref, listClusterEl[set].chromosome_name, ( int)( listClusterEl[set].posEndSV + listClusterEl[set].posStartSV) / 2);
		//idSimpleRep = search_SimpleRepeat( listClusterEl[set].chromosome_name, ( int)( listClusterEl[set].posEndSV+listClusterEl[set].posStartSV) / 2);
		if( ( mei_mid != NULL && strcmp( mei_mid->repeat_type, mobileName) == 0 &&
				mei_mid->strand == SONIC_STRAND_REV) || ( is_mid_satellite != 0))
		{
			MEI_Filter = true;
		}
	}

	if( listClusterEl[set].SVtype == 'A')
	{
		var_example = new_strvar( listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV_Outer, listClusterEl[set].posStartSV,
				listClusterEl[set].posEndSV_Outer, listClusterEl[set].posEndSV, getEnum( listClusterEl[set].SVtype),
				avgEditDistReadSupportingCluster( listClusterEl[set].readMappingSelected), listClusterEl[set].minDelLength,
				listClusterEl[set].maxDelLength, NULL, 1.1, MEI_Filter, false, listClusterEl[set].mobileName,
				listClusterEl[set].readDepth, listClusterEl[set].CNV_Interest, listClusterEl[set].Del_Likelihood,
				listClusterEl[set].indIdCount, listClusterEl[set].sr_support, listClusterEl[set].homogeneity_score);

		print_strvar( in_bams, params, var_example, fpVcf);
	}
	else if( listClusterEl[set].SVtype == 'B')
	{
		var_example = new_strvar( listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV_Outer, listClusterEl[set].posStartSV,
				listClusterEl[set].posEndSV_Outer, listClusterEl[set].posEndSV, getEnum( listClusterEl[set].SVtype),
				avgEditDistReadSupportingCluster( listClusterEl[set].readMappingSelected),listClusterEl[set].minDelLength,
				listClusterEl[set].maxDelLength, NULL, 1.1, MEI_Filter, false, listClusterEl[set].mobileName,
				listClusterEl[set].readDepth, listClusterEl[set].CNV_Interest, listClusterEl[set].Del_Likelihood,
				listClusterEl[set].indIdCount, listClusterEl[set].sr_support, listClusterEl[set].homogeneity_score);
		print_strvar( in_bams, params, var_example, fpVcf);
	}
	else
	{
		var_example = new_strvar( listClusterEl[set].chromosome_name, listClusterEl[set].posStartSV_Outer, listClusterEl[set].posStartSV,
				listClusterEl[set].posEndSV_Outer, listClusterEl[set].posEndSV, getEnum(listClusterEl[set].SVtype),
				avgEditDistReadSupportingCluster( listClusterEl[set].readMappingSelected), listClusterEl[set].minDelLength,
				listClusterEl[set].maxDelLength, NULL, 1.1, listClusterEl[set].LowQual, listClusterEl[set].MEI_Del,
				listClusterEl[set].mobileName, listClusterEl[set].readDepth, listClusterEl[set].CNV_Interest,
				listClusterEl[set].Del_Likelihood, listClusterEl[set].indIdCount, listClusterEl[set].sr_support, listClusterEl[set].homogeneity_score);
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
			//fprintf(stderr,"%d - %d\n", total_sr_sup, total_rp_sup);

			if( ( total_rp_sup > params->rp_threshold) || (listClusterEl[cluster_id].LowQual == true && running_mode == QUICK))
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
		if( listReadEl[ptrReadMapping->readId].readCovered == 0 && listClusterEl[clusterId].indIdCount[ptrReadMapping->indId] > -1)
		{
			supOfIndSeen[ptrReadMapping->indId]++;



			//////////////////////////////////////THE HUERISTIC//////////////////////////
			//Instead of taking the total support as number of reads in each cluster, we use summation of
			//normalized number of reads. i.e each read contributes as 1/(total number of mappings) it has.

			weightedSupOfIndSeen[ptrReadMapping->indId] = weightedSupOfIndSeen[ptrReadMapping->indId]+ptrReadMapping->probEditDist;
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

		mei =  sonic_is_mobile_element(params->this_sonic, listClusterEl[clusterId].chromosome_name, listClusterEl[clusterId].posStartSV - 20,
				listClusterEl[clusterId].posEndSV + 20, params->mei);

		if( ( listClusterEl[clusterId].SVtype == 'D') && mei == NULL)
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
		else if( ( listClusterEl[clusterId].SVtype == 'D') && mei != NULL)
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
		else if( listClusterEl[clusterId].SVtype == 'A' || listClusterEl[clusterId].SVtype == 'B')
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
			{
				bestReads[count] = listClusterEl[bestWeightSet].bestReadToRemove[count];
			}
		}

		if( bestWeight < inf)
		{
			if( conflictResFlag == 0 || conflictsAny( bestWeightSet, bestReads)==-1)
			{
				markReadsCovered( bestWeightSet, bestReads);

				if( conflictResFlag == 1)
					addToConflict( bestWeightSet, bestReads);

				numCallsRequested--;
				listClusterEl[bestWeightSet].oldBestIsGood = 0;
			}
		}
		else{
			if (ten_x_flag == 1 || output_hs_flag ==1)
				free_barcode_homogeneity_hashtable();
			return 0;
		}
	}
	if (ten_x_flag == 1 || output_hs_flag ==1)
		free_barcode_homogeneity_hashtable();
	return 1;
}


void init(bam_info **in_bams, parameters *params, ref_genome* ref, FILE *fpRead, FILE *fpCluster)
{
	char orient1, orient2, readNameStr[strSize], chromosome_name1[strSize], libName[strSize], indName[strSize], chromosome_name2[strSize], MEIName[strSize];
	unsigned long ten_x_barcode;
	int startPos, stopPos, SVtype=0, listReadElId = 0, listClusterElId = 0, readId, readLen, indId, editDist;
	int i, j, count, multiIndId, mQual1, mQual2;
	float editProbDist;
	clusterIdEl *clusterIdElNew;
	readMappingEl *ptrMapping;
	int return_value;

	multiIndCount = 0;

	return_value = fscanf( fpRead, "%i\n", &sizeListReadEl);
	listReadEl = (readEl *) getMem( ( sizeListReadEl) * sizeof( readEl));


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

	/* Put the names of all reads to listReadEl */
	while( fscanf( fpRead, "%s\n", readNameStr) != EOF)
	{
		listReadEl[listReadElId].readName = NULL;
		set_str( &(listReadEl[listReadElId].readName), readNameStr);
		listReadEl[listReadElId].readCovered = 0;
		listReadEl[listReadElId].readId = listReadElId;
		listReadEl[listReadElId].next = NULL;
		listReadElId++;
	}
	listClusterEl = (clusterEl *) getMem( (cluster_count + 1) * sizeof( clusterEl));

	for( count = 0; count < cluster_count; count++)
	{
		listClusterEl[count].clusterId = 0;
		listClusterEl[count].oldBestIsGood = -1;
		listClusterEl[count].next = NULL;
		listClusterEl[count].homogeneity_score = 1.0;
	}
	listClusterElId = 0;
	listClusterEl[listClusterElId].clusterId = listClusterElId;

	while( fscanf( fpCluster, "%s ", readNameStr) != EOF)
	{
		if( strcmp( readNameStr, "END") == 0)
		{
			if( SVtype == 4)
				listClusterEl[listClusterElId].SVtype = 'E';
			if( SVtype == 3)
				listClusterEl[listClusterElId].SVtype = 'V';
			if( SVtype == 2)
				listClusterEl[listClusterElId].SVtype = 'D';
			if(SVtype == 1)
				listClusterEl[listClusterElId].SVtype = 'I';
			if( SVtype == 10) // Insertion from chrN in Forward orientation
				listClusterEl[listClusterElId].SVtype = 'A';
			if( SVtype == 12) // Insertion from chrN in Reverse orientation
				listClusterEl[listClusterElId].SVtype = 'B';

			/* If a new cluster is read indicates an SV bigger than what is expected (user defined threshold) it is pruned out */
			if (clusterElRead_Single.sizeOfCluster > 0)
				processTheSV( in_bams, ref, params, listClusterElId);

			if( SVtype == 10 || SVtype == 12)
			{
				listClusterEl[listClusterElId].mobileName = NULL;
				set_str( &(listClusterEl[listClusterElId].mobileName), MEIName);
			}

			if( ( listClusterEl[listClusterElId].posEndSV - listClusterEl[listClusterElId].posStartSV > maxLengthSV_Del && listClusterEl[listClusterElId].SVtype == 'D' )
					|| ( listClusterEl[listClusterElId].posEndSV - listClusterEl[listClusterElId].posStartSV > maxLengthSV_Inv && listClusterEl[listClusterElId].SVtype == 'V' )
					|| ( listClusterEl[listClusterElId].posEndSV - listClusterEl[listClusterElId].posStartSV > maxLengthSV_TDup && listClusterEl[listClusterElId].SVtype == 'E' ))
			{
				//TROW THE LINK LIST OUT
				ptrMapping = listClusterEl[listClusterElId].next;

				while( ptrMapping != NULL)
				{
					listReadEl[ptrMapping->readId].next = listReadEl[ptrMapping->readId].next->next;
					ptrMapping = ptrMapping->next;
				}
				ptrMapping = listClusterEl[listClusterElId].next;
				freeLinkList( ptrMapping);
			}
			else
				listClusterElId++;

			listClusterEl[listClusterElId].clusterId = listClusterElId;
			clusterElRead_Single.sizeOfCluster = 0;
		}
		else
		{
			if (output_hs_flag != 1 && ten_x_flag != 1 ) {
				return_value = fscanf( fpCluster,"%s %i %s %i %i %g %i %s %s %c %c %d %d ",
						chromosome_name1, &startPos, chromosome_name2, &stopPos, &SVtype,
						&editProbDist, &editDist, libName, indName, &orient1, &orient2, &mQual1, &mQual2);
			} else {
				return_value = fscanf( fpCluster,"%s %i %s %i %i %g %i %s %s %c %c %d %d %lu ",
						chromosome_name1, &startPos, chromosome_name2, &stopPos, &SVtype,
						&editProbDist, &editDist, libName, indName, &orient1, &orient2, &mQual1, &mQual2, &ten_x_barcode);
			}

			/* For MEI (chrN), startPos is for the reference genome and stopPos is for chrN */

			/* Adding the cluster to the read */
			readId = findReadName( readNameStr);
			if (readId == -1)
				continue;
			listReadEl[readId].indId = findIndId2(indName);
			listReadEl[readId].libId = findLibId( libName, indName);

			//Creating a cluster identification element to be added to the list of clusters that this read is in
			clusterIdElNew = (clusterIdEl *) getMem( sizeof( clusterIdEl));
			clusterIdElNew->clusterId = listClusterElId;
			clusterIdElNew->next = listReadEl[readId].next;
			listReadEl[readId].next = clusterIdElNew;

			/* Adding the read to the cluster */
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].readId = readId;
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].indId = findIndId2( indName);
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].editDistance = editDist;
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].mapq1 = mQual1;
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].mapq2 = mQual2;
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].correctMappingQual =
					(1 - pow( 10, -1 * ( ( float)mQual1 / ( float) 10))) * ( 1 - pow( 10, -1 * ( ( float)mQual2 / ( float)10)));

			//fprintf(stderr, "Ind= %s Lib = %s MAPQ = %d - %d, Prob =%f\n", indName, libName, mQual1, mQual2, clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].correctMappingQual);
			if( weightedHeuristicFlag == 1)
				clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].probEditDist = editProbDist;
			else
				clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].probEditDist = 1; /* Since it is unweighted */

			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].posMapLeft = startPos;
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].posMapRight = stopPos;
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].orient1 = orient1;
			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].orient2 = orient2;

			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].chromosome_name = NULL;
			set_str( &clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].chromosome_name, chromosome_name1);

			clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].ten_x_barcode = ten_x_barcode;

			if( SVtype == 10 || SVtype == 12)
				strcpy( MEIName, chromosome_name2);
			clusterElRead_Single.sizeOfCluster++;
		}
	}
	sizeListClusterEl = listClusterElId;
}


void vh_setcover( bam_info **in_bams, parameters *params, ref_genome* ref, char* read_file, char* cluster_file, FILE *fpVcf)
{
	int i, count;
	FILE *readFp = NULL, *clusterFp = NULL;
	//fprintf( stderr,"Inside vh_setcover: 10x flag: %d\n", ten_x_flag);

	numCallsRequested = 10000000;
	sv_count = 0;
	sv_lowqual_count = 0;
	minimumSupNeeded = params->rp_threshold;
	minimumSupNeededNoRD = params->rp_threshold;

	readFp = safe_fopen( read_file,"r");
	clusterFp = safe_fopen( cluster_file,"r");

	init( in_bams, params, ref, readFp, clusterFp);
	pickSet( ref, params);
	outputPickedCluster(in_bams, params, ref, fpVcf);

	cluster_count = 0;

	if( running_mode == QUICK)
		fprintf(stderr, "There are %d SVs and %d LowQual\n", sv_count, sv_lowqual_count);


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
			free_readMappingEl(listClusterEl[i].next);
			free_readMappingEl(listClusterEl[i].readMappingSelected);
		}
		free( listClusterEl);
	}
	free( multiLibs->libName);
	free( multiLibs);
	multiLibs = NULL;

	emptyBuffer();

	for( i = 0; i < sizeListReadEl; i++)
	{
		free( listReadEl[i].readName);
		free_clusterIdEl( listReadEl[i].next);
	}
	free( listReadEl);
	listReadEl = NULL;

	fclose(clusterFp);
	fclose(readFp);
}
