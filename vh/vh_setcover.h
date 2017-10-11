/*
 *  weighted_SetCover_MultiColor_W_ConflictResghted_SetCover_MultiColor_W_ConflictRes2_W_EditEdit.cpp2.h
 *  
 *
 *  Created by Fereydoun Hormozdiari on 24/10/11.
 *  Copyright 2011 Simon Fraser University. All rights reserved.
 *
 */
#ifndef __WEIGHTED_SETCOVER_MULTICOLOR__
#define __WEIGHTED_SETCOVER_MULTICOLOR__

#include <stdlib.h>
#include <stdio.h>
#include "../processbam.h"
#include "../common.h"
#include "../variants.h"

#define SOFTCLIP_WRONGMAP_WINDOW 20
#define maxChroSize 1000000000
#define MaxClusterCount 1000000
#define inf 10000000 // a constant to represent infinity
#define maxLengthSV_Del 500000 // maximum length of SV allowed
#define maxLengthSV_Inv 20000000
#define maxLengthSV_TDup 500000
#define maxListClusterSize 500000 // maximum size of a cluster
#define TRUE 1
#define FALSE 0
#define true 1
#define false 0

extern int cluster_count;
extern int multiIndCount;
extern int sizeListClusterEl;
extern long read_names_count;
extern struct readEl *read_names; // Array of all reads
extern struct clusters_final* clusters_all[MaxClusterCount];

typedef struct clusters_final
{
	int id;
	char* read_name;
	char* chromosome_name1;
	char* chromosome_name2;
	char* mei_type;
	char* mei_subclass;
	char* library_name;
	char* individual_name;
	int start_position;
	int end_position;
	char SV_type;
	char orientation_left;
	char orientation_right;
	float phred_score;
	int edit_distance;
	int mapping_quality_left;
	int mapping_quality_right;
	float correct_mapping_qual;
	bool split_read;
	bool ten_x_barcode;

	struct clusters_final *next;
} clusters_final;

typedef struct multiLib{
	char *libName;
	int indId; // index to the multiInd[X]
	int maxInstSize;
	int minInstSize;
	int readLen;
}multiLib;

multiLib *multiLibs;

////////////////////////////////////////////////THE DATA STRUCTURES TO HOLD INFORMATION FOR EACH MAPPING

/* Creates a linked list for cluster ids */
typedef struct clusterIdEl{
	int clusterId;
	struct clusterIdEl *next;
}clusterIdEl;

/* Holds the information for each read */
typedef struct readEl{
	char *readName;
	int readId;
	int indId;
	int libId;
	int readCovered;//0 not covered, 1 covered;
	struct clusterIdEl *next; // a link list of the clusters which have this read as support
} readEl;


///////////////////////////////////////////////THE DATA STRUCTURES FOR HOLDING INFORMATION FOR EACH CLUSTERS

/*  Linked list to hold information of each Read Mapping for each cluster */
typedef struct readMappingEl{
	int readId; //read id [0 to sizeListReadEl] is an unique identifier for each read (and index to array listReadEL).
	char *chromosome_name;
	char *mei_subclass;
	int posMapLeft;
	int posMapRight;
	int indId;
	int editDistance;
	char orient1;
	char orient2;
	int mapq1;
	int mapq2;
	float probEditDist;
	float correctMappingQual;
	unsigned long ten_x_barcode;
	struct readMappingEl *next; //ptr to the next readMappingEl in this cluster
}readMappingEl;

/* Holds information for each cluster of paired-end reads */
typedef struct clusterEl{
	int clusterId;
	char *chromosome_name;
	int posStartSV;
	int posEndSV;
	int posStartSV_Outer;
	int posEndSV_Outer;
	int minDelLength;// Only used for deletion. Represents the minimum size of deletion predicted
	int maxDelLength;//Only used for deletion. Represents the maximum size of deletion predicted
	char SVtype; //V: inversion, D: Deletion, I: insertion, E: tandem duplication
	int *indIdCount; //If this cluster is picked as one of the SVs it shows the number of supporting paired-end reads selected for each individual (0 : means that we have not picked any support for that inidividual for this SV. -1: means that for this individual this SV is in conflict with SVs picked before - NEVER PICK AN SV FOR AN IND WITH SUP -1 ).
	int *sr_support;
	int oldBestIsGood; // is the last best score computed for this SV still the best or things have changes
	float oldBestScore; // Whats is the best old score
	int bestReadToRemove[totalNumInd]; // the support for each individual for the best old score(if oldBestIsGood is true then this array is still the best selection from this cluster).
	char *mobileName;
	struct readMappingEl *next; //a array of number of individuals for each list of read mappings (i.e. for each individual we have a different list of read mappings and they should be sorted).
	struct readMappingEl *readMappingSelected; // The link list of all the read mappings which have been selected by set cover for this SV.

	float *CNV_Interest; // the CNV calculate for the region of interest
	double *Del_Likelihood;
	double probabilityCNV[totalNumInd][10];// Probability of CNV calculate for 0 to 9
	long *readDepth;

    double weight_without_homogeniety_score;
	double homogeneity_score;
	bool MEI_Del;
	bool LowQual;

    double weight_without_homogeniety_score_at_read_covering;
}clusterEl;

clusterEl *listClusterEl; // the array of all the cluster reads

/* A cluster which is being read and processed before putting all inside a large array (temp cluster) */
typedef struct clusterElRead{
	int clusterId;
	char chromosome_name[strSize];
	char SVtype;
	int sizeOfCluster;
	struct readMappingEl readMappingElArray[maxListClusterSize];
}clusterElRead;


//Used to calculate barcode_homogeneity_score. It makes tha hash table elements
typedef struct barcode_list_element{
	unsigned long ten_x_barcode;
	unsigned long count;
	struct barcode_list_element* next;
}barcode_list_element;

float calWeight( ref_genome* ref, parameters *params, int clusterId, int *countBestSetPicked);
void vh_setcover( bam_info **in_bams, parameters *params, ref_genome* ref, FILE *fpVcf);

#endif
