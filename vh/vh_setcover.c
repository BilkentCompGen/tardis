#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vh_setcover.h"
#include "vh_buffer.h"
#include "vh_conflict.h"
#include "../processfq.h"

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
float mismatchPenalty = 0;
char multiInd[totalNumInd][strSize]; //Name of each individual
struct barcode_list_element* hashtable[hashtable_size];

struct clusters_final* clusters_all[MaxClusterCount];

struct readEl *read_names; // Array of all reads
long read_names_count;// number of total reads

FILE *cnvscoreFile = NULL;
FILE *debugSR = NULL;
//clusterElRead clusterElRead_Single; // the name of cluster which each new cluster is read to it

void free_clusters()
{
	int i;
	clusters_final* tmp_cluster, *tmp_cluster_next;

	for( i = 0; i < MaxClusterCount; i++)
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


int free_readMappingEl( readMappingEl *ptr)
{
	if(ptr == NULL)
		return 0;
	else
	{
		free_readMappingEl( ptr->next);
		if(ptr->chromosome_name != NULL)
			free( ptr->chromosome_name);
		free( ptr);
	}
	return 1;
}
/*
int free_readMappingEl( readMappingEl *ptr)
{
	readMappingEl *ptr_next;

	while( ptr != NULL)
	{
		ptr_next = ptr->next;
		if(ptr->chromosome_name != NULL)
			free( ptr->chromosome_name);
		free( ptr);
		ptr = ptr_next;
	}
	return 1;
}*/

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

/* Add all the read-pairs that support the SV */
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

/* Add all the read-pairs that support the SV */
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

double lpoisson(int d, double lambda, int type) {
	double penalty;
	if (type == 0)
		// read pair
		penalty = 0.001;
	else
		penalty = 0.01;
	if  (lambda == 0.0)
		return log(penalty) * d;
	return d*log(lambda) - lambda - lgamma(d+1);
}

// Calculate score for Deletion, Inversion, Inverted Dup, Inter Dup, Tandem
void calculateCNVScore(bam_info ** in_bam, parameters *params, int clusterId, char SV_type, bool isMEI, double *score) {
	int gc_val, i, count, gc_window, chr = -1;
	double observedReadDepth, expectedReadDepth, observedDiscordantRead, expectedDiscordantRead, mean, std;
	int delta_min, delta_max, read_length;
	int start = 0, end = 0, breakpoint;

	double lhomo, lhete, lnone;

	if( running_mode == SENSITIVE)
	{
		chr = sonic_refind_chromosome_index(params->this_sonic, listClusterEl[clusterId].chromosome_name);
		if( chr == -1)
		{
			fprintf(stderr, "ERROR - cannot find the chromosome in ref");
			exit( 1);
		}
	}

	if (SV_type == DELETION) {
		start = listClusterEl[clusterId].posStartSV;
		end = listClusterEl[clusterId].posEndSV;
		breakpoint = listClusterEl[clusterId].posStartSV_Outer;
	}
	else if (SV_type == INVERSION) {
		start = listClusterEl[clusterId].posStartSV_Outer;
		end = listClusterEl[clusterId].posStartSV;
		breakpoint = listClusterEl[clusterId].posStartSV_Outer;
	}
	else if (SV_type == INVDUPLEFT || SV_type == INTERDUPLEFT) {
		start = listClusterEl[clusterId].posEndSV;
		end = listClusterEl[clusterId].posEndSV_Outer;
		breakpoint = listClusterEl[clusterId].posStartSV_Outer;
	}
	else if (SV_type == INVDUPRIGHT || SV_type == INTERDUPRIGHT) {
		start = listClusterEl[clusterId].posStartSV_Outer;
		end = listClusterEl[clusterId].posStartSV;
		breakpoint = listClusterEl[clusterId].posEndSV_Outer;
	}
	else if (SV_type == TANDEMDUP) {
		start = listClusterEl[clusterId].posStartSV;
		end = listClusterEl[clusterId].posEndSV;
		breakpoint = listClusterEl[clusterId].posEndSV_Outer;
	}
	else if (SV_type == MEIFORWARD) {
		breakpoint = listClusterEl[clusterId].posStartSV_Outer;
	}
	else if (SV_type == MEIREVERSE) {
		breakpoint = listClusterEl[clusterId].posEndSV_Outer;
	}
	else if (SV_type == NUMTFORWARD) {
		breakpoint = listClusterEl[clusterId].posStartSV_Outer;
	}
	else if (SV_type == NUMTREVERSE) {
		breakpoint = listClusterEl[clusterId].posEndSV_Outer;
	}
	else if(SV_type == INSERTION) {
		breakpoint = listClusterEl[clusterId].posStartSV_Outer;
	}
	if (SV_type == DELETION) {
		start = listClusterEl[clusterId].posStartSV;
		end = listClusterEl[clusterId].posEndSV;
		breakpoint = listClusterEl[clusterId].posStartSV_Outer;
	}

	for( count = 0; count < multiIndCount; count++)
	{
		observedReadDepth = 0;
		expectedReadDepth = 0;

		for( i = start; i < end; i++)
		{
			gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, listClusterEl[clusterId].chromosome_name, i, i + WINDOWSLIDE));
			expectedReadDepth += in_bam[count]->mean_rd_per_gc[gc_val];

			if( running_mode == QUICK)
			{
				observedReadDepth += ( double)in_bam[count]->read_depth_per_chr[i];
			}
			else if( running_mode == SENSITIVE)
			{
				observedReadDepth += ( double)in_bam[count]->read_depth[chr][i];
			}
		}

		observedDiscordantRead = 0;
		expectedDiscordantRead = 0;
		readMappingEl *ptr = listClusterEl[clusterId].next;

		while (ptr != NULL) {
			if (read_names[ptr->readId].readCovered == 0){
				observedDiscordantRead++;
			}
			ptr = ptr->next;
		}

		delta_min = in_bam[count]->libraries[0]->conc_min;
		delta_max = in_bam[count]->libraries[0]->conc_max;
		mean = in_bam[count]->libraries[0]->frag_avg;
		std = in_bam[count]->libraries[0]->frag_std;
		read_length = in_bam[count]->libraries[0]->read_length;

		for( i = read_length; i <= delta_max; i++){
			expectedDiscordantRead += 1 - 0.5 * ( 1 + erf( ( i - mean) / ( std * sqrt( 2.))));
		}

		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, listClusterEl[clusterId].chromosome_name, breakpoint, breakpoint + WINDOWSLIDE));
		expectedDiscordantRead *= in_bam[count]->mean_rd_per_gc[gc_val];


		fprintf(cnvscoreFile,"chr= %s Id= %d observedRD= %lf expectedRD= %lf observedRP= %lf expectedRP= %lf\n", listClusterEl[clusterId].chromosome_name, clusterId, observedReadDepth, expectedReadDepth, observedDiscordantRead, expectedDiscordantRead);
		fprintf(cnvscoreFile, "Id= %d start= %d end= %d breakpoint= %d gc_val= %d meanrdpergc= %lf\n", clusterId, start, end, breakpoint, gc_val, in_bam[count]->mean_rd_per_gc[gc_val]);
		//fprintf(cnvscoreFile,"Id= %d delta_min= %d delta_max= %d mean= %lf std= %lf\n", clusterId, delta_min, delta_max, mean, std);

		int READPAIR = 0;
		int READDEPTH = 1;


		if (SV_type == DELETION && isMEI == true) {
			lhomo = lpoisson(observedDiscordantRead,  expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= DELETION MEI lhomo= %lf lhete= %lf lnone= %lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == DELETION && isMEI == false) {
			lhomo = lpoisson(observedReadDepth, 0.0, READDEPTH) + lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedReadDepth, 0.5 * expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedReadDepth, expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= DELETION lhomo= %lf lhete= %lf lnone= %lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == INVERSION) {
			lhomo = lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= INVERSION lhomo= %lf lhete= %lf lnone= %lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == INVDUPLEFT || SV_type == INVDUPRIGHT) {
			lhomo = lpoisson(observedReadDepth, 2 * expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedReadDepth, 1.5 * expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedReadDepth, expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= INVDUP lhomo= %lf lhete= %lf lnone=%lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == INTERDUPLEFT || SV_type == INTERDUPRIGHT) {
			lhomo = lpoisson(observedReadDepth, 2 * expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedReadDepth, 1.5 * expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedReadDepth, expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= INTERDUP lhomo= %lf lhete= %lf lnone=%lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == TANDEMDUP) {
			lhomo = lpoisson(observedReadDepth, 2 * expectedReadDepth, READDEPTH) +  lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedReadDepth, 1.5 * expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedReadDepth, expectedReadDepth, READDEPTH) + lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= TANDEMDUP lhomo= %lf lhete= %lf lnone= %lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == MEIFORWARD || SV_type == MEIREVERSE) {
			lhomo = lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= MEI lhomo= %lf lhete= %lf lnone= %lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == NUMTFORWARD || SV_type == NUMTREVERSE) {
			lhomo = lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= NUMT lhomo= %lf lhete= %lf lnone= %lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else if (SV_type == INSERTION) {
			lhomo = lpoisson(observedDiscordantRead, expectedDiscordantRead, READPAIR);
			lhete = lpoisson(observedDiscordantRead, 0.5 * expectedDiscordantRead, READPAIR);
			lnone = lpoisson(observedDiscordantRead, 0.0, READPAIR);
			score[count] = max(lhomo, lhete) / lnone;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= Insertion lhomo= %lf lhete= %lf lnone= %lf score= %lf\n", clusterId, count, lhomo, lhete, lnone, score[count]);
		}
		else {
			score[count] = inf;
			fprintf(cnvscoreFile,"Id= %d BAM= %d SV_type= %c lhomo= xxx lhete= xxx lnone= xxx score= inf\n", clusterId, count, SV_type);
		}

		if (observedDiscordantRead != 0){
			score[count] = score[count] / observedDiscordantRead;
			fprintf(cnvscoreFile,"Id= %d newscore(score/num_reads)= %lf\n", clusterId, score[count]);
		}
		else{
			score[count] = inf;
			fprintf(cnvscoreFile,"Id= %d newscore(score/num_reads)=inf\n", clusterId);
		}
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

void init_barcode_homogeneity_hashtable()
{
	int hashtable_index;

	//initialize the heads of the chains in the hashtable
	for( hashtable_index = 0; hashtable_index < hashtable_size; hashtable_index++)
	{
		hashtable[hashtable_index] = ( barcode_list_element*) getMem( sizeof( barcode_list_element));
		hashtable[hashtable_index]->ten_x_barcode = 0;
		hashtable[hashtable_index]->count = 0;
		hashtable[hashtable_index]->next = NULL;
	}
}

void free_barcode_homogeneity_hashtable()
{
	int i;
	struct barcode_list_element *temp_ptr;
	for( i = 0; i < hashtable_size; i++)
	{
		while(hashtable[i]->next)
		{
			temp_ptr = hashtable[i]->next->next;
			free( hashtable[i]->next);
			hashtable[i]->next = temp_ptr;
		}
		free( hashtable[i]);
	}
}

double barcode_homogeneity_score(int clusterId)
{
	unsigned long total_read_count = 0;
	double sum = 0.0, score = 1.0;
	int i, found, hashtable_index;

	struct barcode_list_element *current_list_element;
	struct barcode_list_element *temp_ptr;
	readMappingEl *ptrReadMapping;

	ptrReadMapping = listClusterEl[clusterId].next;
	while( ptrReadMapping != NULL)
	{
		if( ( read_names[ptrReadMapping->readId].readCovered == 0)
				&& ( listClusterEl[clusterId].indIdCount[ptrReadMapping->indId] > -1)
				&& ( signed long)ptrReadMapping->ten_x_barcode != -1)
		{
			total_read_count++;
			//fprintf(stderr,"%d - %lu\n", ptrReadMapping->readId, ptrReadMapping->ten_x_barcode);
			hashtable_index = ptrReadMapping->ten_x_barcode % hashtable_size;
			current_list_element = hashtable[hashtable_index];

			found = 0;
			while( current_list_element->next)
			{
				if( current_list_element->next->ten_x_barcode == ptrReadMapping->ten_x_barcode)
				{
					current_list_element->next->count++;
					found = 1;
					break;
				}
				current_list_element = current_list_element->next;
			}

			if( found == 0)
			{
				current_list_element->next = ( barcode_list_element*)getMem( sizeof( barcode_list_element));
				current_list_element->next->ten_x_barcode = ptrReadMapping->ten_x_barcode;
				current_list_element->next->count = 1;
				current_list_element->next->next = NULL;
			}
		}
		ptrReadMapping = ptrReadMapping->next;
	}

	for( hashtable_index = 0; hashtable_index < hashtable_size; hashtable_index++)
	{
		current_list_element = hashtable[hashtable_index];

		while( current_list_element->next)
		{
			sum = sum + ( ( double)current_list_element->next->count * ( double)current_list_element->next->count);
			current_list_element = current_list_element->next;
		}
	}

	for( i = 0; i < hashtable_size; i++)
	{
		while( hashtable[i]->next)
		{
			temp_ptr = hashtable[i]->next->next;
			free(hashtable[i]->next);
			hashtable[i]->next = temp_ptr;
		}
	}


	if( total_read_count == 0)
	{
		listClusterEl[clusterId].homogeneity_score = score;
		return 1.0;
	}

	score = sum / ( total_read_count * total_read_count);

	//fprintf(stderr, "score= %lf, sum=%.4lf - count=%li\n\n\n\n", score, sum, total_read_count);

	listClusterEl[clusterId].homogeneity_score = score;

	if( ten_x_flag != 1)
		return 1.0;

	if( !( score > 1 || score < 0))
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

void markReadsCovered( int clusterId, int *countReads)
{
	int minDelLength = 1000000000, maxDelLength = 1000000000, totalReadRemove = 0, count;
	int readsToMark[totalNumInd];

	listClusterEl[clusterId].weight_without_homogeniety_score_at_read_covering = listClusterEl[clusterId].weight_without_homogeniety_score;
	readMappingEl *ptrRead;
	clusterIdEl *ptrClusterId;
	ptrRead = listClusterEl[clusterId].next;

	for( count = 0; count < multiIndCount; count++)
	{
		totalReadRemove += countReads[count];
		readsToMark[count] = countReads[count];
	}

	//This part will mark the reads covered and calculate the minimum and maximum length of deletion
	while( ptrRead != NULL && totalReadRemove > 0)
	{
		if( readsToMark[ptrRead->indId] > 0 && read_names[ptrRead->readId].readCovered == 0)
		{
			ptrClusterId = read_names[ptrRead->readId].next;

			/* Count the number of read pairs supporting this cluster */
			listClusterEl[clusterId].indIdCount[ptrRead->indId]++;

			/* Count the number of split reads supporting this cluster */
			if( strcmp( multiLibs[read_names[ptrRead->readId].libId].libName, "SplitRead") == 0)
				listClusterEl[clusterId].sr_support[ptrRead->indId]++;

			listClusterEl[clusterId].readMappingSelected = addToListOfReads( listClusterEl[clusterId].readMappingSelected, *ptrRead);

			if( listClusterEl[clusterId].SVtype == DELETION && minDelLength > ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[read_names[ptrRead->readId].libId].maxInstSize)
				minDelLength = ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[read_names[ptrRead->readId].libId].maxInstSize;

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

	if (listClusterEl[clusterId].SVtype == DELETION)
	{
		listClusterEl[clusterId].minDelLength = minDelLength;
		listClusterEl[clusterId].maxDelLength = maxDelLength;
	}
}


void outputCluster( bam_info** in_bams, parameters* params, int cluster_id, FILE *fpVcf)
{
	int i;
	char mobileName[strSize];

	int is_start_satellite, is_end_satellite, is_mid_satellite;
	sonic_repeat *mei_start, *mei_end, *mei_mid;
	bool MEI_Filter = false;
	struct strvar* var_example = NULL;

	for( i = 0; i < params->num_bams; i++)
		in_bams[i]->contribution = false;

	for( i = 0; i < params->num_bams; i++)
	{
		if( listClusterEl[cluster_id].indIdCount[i] > 0)
			in_bams[i]->contribution = true;
	}

	is_start_satellite = sonic_is_satellite(params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV, listClusterEl[cluster_id].posStartSV + 1);
	is_end_satellite = sonic_is_satellite(params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posEndSV - 1, listClusterEl[cluster_id].posEndSV);
	is_mid_satellite = sonic_is_satellite(params->this_sonic, listClusterEl[cluster_id].chromosome_name, (int)( listClusterEl[cluster_id].posStartSV + listClusterEl[cluster_id].posEndSV) / 2,
			(int)( listClusterEl[cluster_id].posStartSV + listClusterEl[cluster_id].posEndSV) / 2 + 1);

	mei_start = sonic_is_mobile_element( params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV, listClusterEl[cluster_id].posStartSV + 1, params->mei);
	mei_end = sonic_is_mobile_element( params->this_sonic, listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posEndSV - 1, listClusterEl[cluster_id].posEndSV, params->mei);
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

	if( listClusterEl[cluster_id].SVtype == MEIFORWARD || listClusterEl[cluster_id].SVtype == NUMTFORWARD)
	{
		var_example = new_strvar( listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV_Outer, listClusterEl[cluster_id].posStartSV,
				listClusterEl[cluster_id].posEndSV_Outer, listClusterEl[cluster_id].posEndSV, listClusterEl[cluster_id].SVtype,
				MEI_Filter, false, listClusterEl[cluster_id].mobileName, listClusterEl[cluster_id].mei_type, listClusterEl[cluster_id].CNV_Score,
				listClusterEl[cluster_id].indIdCount, listClusterEl[cluster_id].sr_support, listClusterEl[cluster_id].homogeneity_score,
				listClusterEl[cluster_id].weight_without_homogeniety_score_at_read_covering);

		print_strvar( in_bams, params, var_example, fpVcf);
	}
	else if( listClusterEl[cluster_id].SVtype == MEIREVERSE || listClusterEl[cluster_id].SVtype == NUMTREVERSE)
	{
		var_example = new_strvar( listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV_Outer, listClusterEl[cluster_id].posStartSV,
				listClusterEl[cluster_id].posEndSV_Outer, listClusterEl[cluster_id].posEndSV, listClusterEl[cluster_id].SVtype,
				MEI_Filter, false, listClusterEl[cluster_id].mobileName, listClusterEl[cluster_id].mei_type, listClusterEl[cluster_id].CNV_Score,
				listClusterEl[cluster_id].indIdCount, listClusterEl[cluster_id].sr_support, listClusterEl[cluster_id].homogeneity_score,
				listClusterEl[cluster_id].weight_without_homogeniety_score_at_read_covering);

		print_strvar( in_bams, params, var_example, fpVcf);
	}
	else
	{
		var_example = new_strvar( listClusterEl[cluster_id].chromosome_name, listClusterEl[cluster_id].posStartSV_Outer, listClusterEl[cluster_id].posStartSV,
				listClusterEl[cluster_id].posEndSV_Outer, listClusterEl[cluster_id].posEndSV, listClusterEl[cluster_id].SVtype,
				listClusterEl[cluster_id].LowQual, listClusterEl[cluster_id].MEI_Del,
				listClusterEl[cluster_id].mobileName, listClusterEl[cluster_id].mei_type, listClusterEl[cluster_id].CNV_Score,
				listClusterEl[cluster_id].indIdCount, listClusterEl[cluster_id].sr_support,
				listClusterEl[cluster_id].homogeneity_score, listClusterEl[cluster_id].weight_without_homogeniety_score_at_read_covering);

		print_strvar( in_bams, params, var_example, fpVcf);
	}
	free( var_example);
}

void outputPickedCluster( bam_info** in_bams, parameters* params, FILE *fpVcf)
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
				outputCluster( in_bams, params, cluster_id, fpVcf);
		}
	}
}


float calWeight( bam_info **in_bams, parameters *params, int clusterId, int *countBestSetPicked)
{
	//int idIndCanBePicked = 0; // the id of individuals which we are considering to have this SV in this round
	//int editDistanceSum[totalNumInd]; // the total sum of edit distance for paired-end reads considered till now
	int supOfIndSeen[totalNumInd]; //the number of paired-end read supports for each individual seen in this round
	//double weightedSupOfIndSeen[totalNumInd];// a normalized number of paired-end read supports for each individual seen
	//int totalEditDist = 0; // total of edit distance of all paired-end reads for individuals which we are considering
	float bestScore, weightNew, normalizedWeightedSup = 0;
	int numReadCanBeCovered = 0, indCount, count, i, k;

	bool passMinSup;
	float bestScoreRD = 0;
	double maxRDlikelihoodDel, maxRDlikelihoodDup;
	double sumLikelihoodDup, sumLikelihoodDel, sumLikelihoodNorm;
	double epsilon = 1e-5; /* was 1e-200 before */
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
		//editDistanceSum[indCount] = 0;
		countBestSetPicked[indCount] = 0;
		//weightedSupOfIndSeen[indCount] = 0;
	}

	bestScore = inf;
	readMappingEl *ptrReadMapping;
	ptrReadMapping = listClusterEl[clusterId].next;

	while( ptrReadMapping != NULL)
	{
		if( read_names[ptrReadMapping->readId].readCovered == 0 && listClusterEl[clusterId].indIdCount[ptrReadMapping->indId] > -1)
		{
			supOfIndSeen[ptrReadMapping->indId]++;

			/* Count the number of split reads supporting this cluster */
			if( strcmp( multiLibs[read_names[ptrReadMapping->readId].libId].libName, "SplitRead") == 0)
				listClusterEl[clusterId].sr_support[ptrReadMapping->indId]++;

			//////////////////////////////////////THE HUERISTIC//////////////////////////
			//Instead of taking the total support as number of reads in each cluster, we use summation of
			//normalized number of reads. i.e each read contributes as 1/(total number of mappings) it has.

			//weightedSupOfIndSeen[ptrReadMapping->indId] = weightedSupOfIndSeen[ptrReadMapping->indId] + ptrReadMapping->probEditDist;
			if( supOfIndSeen[ptrReadMapping->indId] == minimumSupNeeded)
			{
				//idIndCanBePicked = idIndCanBePicked + ( int)pow( 2, ptrReadMapping->indId);
				numReadCanBeCovered = numReadCanBeCovered + supOfIndSeen[ptrReadMapping->indId];
				//totalEditDist = totalEditDist + editDistanceSum[ptrReadMapping->indId];
			}
			else if( supOfIndSeen[ptrReadMapping->indId] > minimumSupNeeded)
			{
				numReadCanBeCovered++;
				//totalEditDist = totalEditDist + ptrReadMapping->editDistance;
			}
		}
		ptrReadMapping = ptrReadMapping->next;
	}

	//weightNew = weightsForCombination[idIndCanBePicked]+scoreForEditDistanceFunc(totalEditDist);
	weightNew = 1;

	if( numReadCanBeCovered > 0)
	{
		mei = sonic_is_mobile_element(params->this_sonic, listClusterEl[clusterId].chromosome_name, listClusterEl[clusterId].posStartSV - 20,
				listClusterEl[clusterId].posEndSV + 20, params->mei);

		if( mei != NULL)
		{
			calculateCNVScore(in_bams, params, clusterId, listClusterEl[clusterId].SVtype, true, listClusterEl[clusterId].CNV_Score);
			bestScore = listClusterEl[clusterId].CNV_Score[0];
			if( listClusterEl[clusterId].SVtype == DELETION)
			{
				listClusterEl[clusterId].MEI_Del = true;
				listClusterEl[clusterId].mei_type = NULL;
				set_str( &(listClusterEl[clusterId].mei_type), mei->repeat_class);
			}
		}
		else
		{
			calculateCNVScore(in_bams, params, clusterId, listClusterEl[clusterId].SVtype, false, listClusterEl[clusterId].CNV_Score);
			bestScore = listClusterEl[clusterId].CNV_Score[0];
		}
	}
	else
		bestScore = inf;

	for( count = 0; count < multiIndCount; count++)
	{
		if( supOfIndSeen[count] >= minimumSupNeeded)
			countBestSetPicked[count] = supOfIndSeen[count];
	}
	double tmp = 0;
	listClusterEl[clusterId].weight_without_homogeniety_score = bestScore;
	if( bestScore < inf && (ten_x_flag == 1 || output_hs_flag == 1))
	{
		tmp = barcode_homogeneity_score(clusterId);
		bestScore = ( float) bestScore * tmp;
	}

	//fprintf(stderr, "%.4lf - %.4lf\n", tmp, listClusterEl[clusterId].weight_without_homogeniety_score);

	listClusterEl[clusterId].oldBestIsGood = 1;
	listClusterEl[clusterId].oldBestScore = bestScore;

	for( count = 0; count < totalNumInd; count++)
		listClusterEl[clusterId].bestReadToRemove[count] = countBestSetPicked[count];
	return bestScore;
}

void debug_split_reads(int clusterId, float bestWeight) {
	fprintf(debugSR,"Id= %d weight= %f StartOuter= %d Start= %d End= %d EndOuter= %d SVtype= %c sr_support[0]= %d readNamesNotCovered= \t",
			clusterId, bestWeight, listClusterEl[clusterId].posStartSV_Outer, listClusterEl[clusterId].posStartSV,
			listClusterEl[clusterId].posEndSV, listClusterEl[clusterId].posEndSV_Outer,
			listClusterEl[clusterId].SVtype, listClusterEl[clusterId].sr_support[0]);

	readMappingEl *ptrReadMapping;
	ptrReadMapping = listClusterEl[clusterId].next;

	int countSupport = 0;

	while( ptrReadMapping != NULL)
	{
		if( read_names[ptrReadMapping->readId].readCovered == 0)
		{
			fprintf(debugSR,"%s\t", read_names[ptrReadMapping->readId].readName);
			countSupport++;
		}
		ptrReadMapping = ptrReadMapping->next;
	}
	fprintf(debugSR,"CountSupport= %d\n", countSupport);
}

int pickSet(bam_info **in_bams, parameters *params)
{
	float bestWeight;
	float newWeight;
	int countReads[totalNumInd];
	int bestReads[totalNumInd];
	int bestWeightCluster;
	int count, clusterCounter;
	maxScoreInBuffer = inf - 1;
	//int isClusterSelected[sizeListClusterEl];

	//for( clusterCounter = 0; clusterCounter < sizeListClusterEl; clusterCounter++)
	//isClusterSelected[clusterCounter] = 0;

	if (ten_x_flag == 1 || output_hs_flag ==1)
		init_barcode_homogeneity_hashtable();

	while( numCallsRequested > 0)
	{
		bestWeight = inf;
		bestWeightCluster = -1;

		if( !bufferIsUseful( in_bams, params))
		{
			emptyBuffer();
			//fprintf(stderr, "EMPTY BUFFER\n");
			for( clusterCounter = 0; clusterCounter < sizeListClusterEl; clusterCounter++)
			{
				//if (isClusterSelected[clusterCounter]) continue;

				newWeight = calWeight(in_bams, params, clusterCounter, countReads);

				addToBuffer( newWeight, clusterCounter);

				if( newWeight < bestWeight)
				{
					bestWeightCluster = clusterCounter;

					for( count = 0; count < multiIndCount; count++)
						bestReads[count] = countReads[count];

					bestWeight = newWeight;
				}
			}
			//fprintf(stderr,"k=%d, weight = %f, support =%d, bw %f- bwclu %d\n",k++, newWeight, countReads[0], bestWeight, bestWeightCluster);
		}
		else
		{
			bestWeightCluster = bestFromBuffer();
			bestWeight = listClusterEl[bestWeightCluster].oldBestScore;

			for ( count = 0; count < multiIndCount; count++)
				bestReads[count] = listClusterEl[bestWeightCluster].bestReadToRemove[count];
		}

		if( bestWeight < inf)
		{
			if( conflictResFlag == 0 || conflictsAny( bestWeightCluster, bestReads) == -1)
			{
				//debug_split_reads(bestWeightCluster, bestWeight);
				markReadsCovered( bestWeightCluster, bestReads);
				//isClusterSelected[bestWeightCluster] = 1; // mark cluster as selected

				if( conflictResFlag == 1)
					addToConflict( bestWeightCluster, bestReads);

				numCallsRequested--;
				listClusterEl[bestWeightCluster].oldBestIsGood = 0;
			}
		}
		else
		{
			if (ten_x_flag == 1 || output_hs_flag == 1)
				free_barcode_homogeneity_hashtable();
			return 0;
		}
	}
	if (ten_x_flag == 1 || output_hs_flag == 1)
		free_barcode_homogeneity_hashtable();
	return 1;
}


/* Each new cluster is copied into a new cell in array of clusters (listClusterEl) in index (clusterId) */
void processTheSV( bam_info **in_bams, parameters *params, int listClusterId, int cluster_count)
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
	set_str( &(listClusterEl[listClusterId].mobileName), tmp_cluster->mei_subclass);

	listClusterEl[listClusterId].mei_type = NULL;
	set_str( &(listClusterEl[listClusterId].mei_type), tmp_cluster->mei_type);

	listClusterEl[listClusterId].SVtype = tmp_cluster->SV_type;
	listClusterEl[listClusterId].next = NULL;

	listClusterEl[listClusterId].CNV_Score = (double *) getMem( multiIndCount * sizeof( double));
	listClusterEl[listClusterId].indIdCount = (int *) getMem( multiIndCount * sizeof( int));
	listClusterEl[listClusterId].sr_support = (int *) getMem( multiIndCount * sizeof( int));

	for( count = 0; count < multiIndCount; count++)
	{
		listClusterEl[listClusterId].CNV_Score[count] = 0;
		listClusterEl[listClusterId].indIdCount[count] = 0;
		listClusterEl[listClusterId].sr_support[count] = 0;
	}

	if( listClusterEl[listClusterId].SVtype == MEIFORWARD || listClusterEl[listClusterId].SVtype == MEIREVERSE)
	{
		while( tmp_cluster != NULL)
		{
			if( tmp_cluster->orientation_left == FORWARD)
			{
				if( tmp_cluster->start_position > posStartSV)
					posStartSV = tmp_cluster->start_position;
				if( tmp_cluster->start_position < posStartSV_Outer)
					posStartSV_Outer = tmp_cluster->start_position;
			}
			if( tmp_cluster->orientation_left == REVERSE)
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
	else if( listClusterEl[listClusterId].SVtype == NUMTFORWARD || listClusterEl[listClusterId].SVtype == NUMTREVERSE)
	{
		while( tmp_cluster != NULL)
		{
			if (tmp_cluster->orientation_left == FORWARD)
			{
				if( tmp_cluster->start_position > posStartSV)
					posStartSV = tmp_cluster->start_position;
				if( tmp_cluster->start_position < posStartSV_Outer)
					posStartSV_Outer = tmp_cluster->start_position;
			}
			if( tmp_cluster->orientation_left == REVERSE)
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
	else if( listClusterEl[listClusterId].SVtype == INVDUPLEFT || listClusterEl[listClusterId].SVtype == INVDUPRIGHT
			|| listClusterEl[listClusterId].SVtype == INTERDUPLEFT || listClusterEl[listClusterId].SVtype == INTERDUPRIGHT)
	{
		while( tmp_cluster != NULL)
		{
			if( tmp_cluster->orientation_left == FORWARD)
			{
				if( tmp_cluster->start_position > posStartSV)
					posStartSV = tmp_cluster->start_position;

				if( tmp_cluster->start_position - 100 < posStartSV_Outer)
					posStartSV_Outer = max(0, tmp_cluster->start_position - 100);

				if( tmp_cluster->end_position - 100 < posEndSV)
					posEndSV = max(0, tmp_cluster->end_position - 100);

				if( tmp_cluster->end_position > posEndSV_Outer)
					posEndSV_Outer = tmp_cluster->end_position;

				listClusterEl[listClusterId].next = addToListOfReads2( listClusterEl[listClusterId].next, tmp_cluster);
			}
			else if( tmp_cluster->orientation_left == REVERSE)
			{
				if( tmp_cluster->start_position + 100 > posStartSV)
					posStartSV = tmp_cluster->start_position + 100;

				if( tmp_cluster->start_position < posStartSV_Outer)
					posStartSV_Outer = tmp_cluster->start_position;

				if( tmp_cluster->end_position < posEndSV)
					posEndSV = tmp_cluster->end_position;

				if( tmp_cluster->end_position + 100 > posEndSV_Outer)
					posEndSV_Outer = max(0, tmp_cluster->end_position + 100);

				listClusterEl[listClusterId].next = addToListOfReads2( listClusterEl[listClusterId].next, tmp_cluster);
			}
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
}


void init(bam_info **in_bams, parameters *params)
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
		if( running_mode == QUICK && params->alt_mapping != 0)
		{
			multiLibs[multiLibsCount].libName = NULL;
			set_str( &multiLibs[multiLibsCount].libName, "Alternative");
			multiLibs[multiLibsCount].indId = addNewInd( in_bams[j]->sample_name);
			multiLibs[multiLibsCount].maxInstSize =  0;
			multiLibs[multiLibsCount].minInstSize = 0;
			multiLibs[multiLibsCount].readLen = 0;
			multiLibsCount++;
		}
	}
	/* We use indCount in variants.c */
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
		if( cluster_size > 0)
			processTheSV( in_bams, params, listClusterElId, count);

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

	/* The clusters created in clustering step are freed */
	free_clusters();
}


void vh_setcover( bam_info **in_bams, parameters *params, FILE *fpVcf)
{
	int i;
	char *debugsr_path;
	char *outputscore_path;
	
	numCallsRequested = maxNumSV;

	/* create conflict array */
	listSelectedSV = (struct SV_selected *) getMem(sizeof(struct SV_selected) * maxNumSV);
	
	sv_count = 0;
	sv_lowqual_count = 0;
	minimumSupNeeded = params->rp_threshold;

	fprintf( stderr, "\nApplying set cover for %d clusters", cluster_count);

	if( running_mode == SENSITIVE)
		fprintf( logFile, "Total Cluster Count= %d\n", cluster_count);

	outputscore_path = (char *) getMem(sizeof(char) * (1+strlen("output.score")+strlen(params->outdir)));
	sprintf(outputscore_path, "%s%s", params->outdir, "output.score");
	cnvscoreFile = safe_fopen(outputscore_path, "w");
	free(outputscore_path);
	
	debugsr_path = (char *) getMem(sizeof(char) * (1+strlen("debug.sr")+strlen(params->outdir)));
	sprintf(debugsr_path, "%s%s", params->outdir, "debug.sr");
	debugSR = safe_fopen(debugsr_path, "w");
	free(debugsr_path);
	
	init( in_bams, params);
	fprintf( stderr, ".");
	fflush( stderr);
	pickSet( in_bams, params);
	fprintf( stderr, ".");
	fflush( stderr);
	outputPickedCluster( in_bams, params, fpVcf);
	fprintf( stderr, ".");
	fflush( stderr);
	fprintf( stderr, "\n\n");


	fprintf( stderr, "There are %d SVs\n", sv_count);

	fclose(cnvscoreFile);
	fclose(debugSR);


	/* Free the structures, etc */
	if( listClusterEl != NULL)
	{
		for( i = 0; i < sizeListClusterEl; i++)
		{
			free( listClusterEl[i].mobileName);
			free( listClusterEl[i].mei_type);
			free( listClusterEl[i].chromosome_name);
			free( listClusterEl[i].indIdCount);
			free( listClusterEl[i].CNV_Score);
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
	free( listSelectedSV);
}
