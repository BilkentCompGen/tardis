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

Interval *gaps;
Interval *repeats;
Interval *dups;

int repeat_count;
int dup_count;
int gap_count = 0;

Mei *g_meiTable;
int meiTypeSize; //Number of mobile element types that we want to check e.g. ALU, SVA
int selectedMeiTypeSize; // Total number of mobile elements of the selected mei types
int g_maxRepeatLength = 0;

MeiIndex *mInd;


int vh_compareIntInterval (const void *a, const void *b)
{
	return (((Interval *) a)->start - ((Interval *) b)->start);
}

int vh_compareCharInterval (const void *a, const void *b)
{
	return strcmp(((Interval *) a)->chroName, ((Interval *) b)->chroName);
}

int compareStartMeiTable (const void * a, const void * b)
{
	return (((Mei *) a)->start - ((Mei *) b)->start);
}

int compareChrNameMeiTable (const void * a, const void * b)
{
	return strcmp(((Mei *) a)->chroName, ((Mei *) b)->chroName);
}

int vh_noGap (char *chroName, int startMapping, int endMapping)
{
	int count;
	for (count = 0; count < gap_count; count++)
	{
		if (strcmp (chroName, gaps[count].chroName) == 0)
		{
			if (startMapping <= gaps[count].start && endMapping >= gaps[count].end)
				return 0;
			if (startMapping >= gaps[count].start && startMapping <= gaps[count].end)
				return 0;
			if (endMapping >= gaps[count].start && endMapping <= gaps[count].end)
				return 0;
		}
	}
	return 1;
}

/*Decomposes the input for MEI into tokens e.g Alu:L1Hs is an input */
char** vh_parseMeiString(char *meiString)
{
	char **mei_type;
	char a[50];
	char *tok;
	int i=0;
	char *str=meiString;
	meiTypeSize=1;

	/* Count number of MEIs. */
	while (*str!='\0')
	{
		if (':' == *str)
		{
			meiTypeSize++;
		}
		str++;
	}
	mei_type = (char **)getMem( ( meiTypeSize + 1 ) * sizeof( char* ) );
	str=meiString;

	tok = strtok(str, ":");
	while(tok!=NULL)
	{
		mei_type[i]=(char*)getMem(strlen(tok)+1 * sizeof(char));
		strcpy(mei_type[i], tok);

		i++;
		tok = strtok(NULL, ":");
	}

	return mei_type;
}


void vh_readGapTable (char *gapFileName)
{
	int c,ln_count=0;
	char chrom[STRMAX];
	FILE *gapFile = safe_fopen (gapFileName, "r");
	if (!gapFile);
	//TODO: handle me

	/* Count the number of lines in the gap file */
	do{
		c=fgetc(gapFile);
		if(c=='\n') ln_count++;

	}while (c!=EOF);

	/* Reset the file pointer to the start of the file */
	rewind(gapFile);

	gaps = (Interval*) getMem( ln_count * sizeof(Interval));

	int index = 0;
	while (!feof (gapFile))
	{
		int res = fscanf (gapFile, "%s%d%d\n", chrom, &(gaps[index].start), &(gaps[index].end));
		gaps[index].chroName=(char *) getMem(sizeof(char) * (strlen(chrom)+1));
		strcpy(gaps[index].chroName, chrom);

		if (res != 3)
			break;

		index++;
	}
	gap_count = index;

	fclose (gapFile);
}


void vh_readDupsTable( char *dupsFileName)
{
	char chrom[STRMAX];
	int start, end, res, i;
	FILE *dupsFile = safe_fopen (dupsFileName, "r");
	dup_count=0;

	while( !feof( dupsFile))
	{
		res = fscanf( dupsFile,"%s\t%d\t%d\n",chrom, &start, &end);
		if (res != 3)
		{
			// todo: global exit function
			fprintf(stderr, "The dups file %s is not in the correct BED6 format. Consult the manual on how to generate a valid file.\n", dupsFileName);
			exit(-1);
		}
		dup_count++;
	}
	if ( dup_count == 0)
		fprintf( stderr, "No duplications found\n");
	else
		fprintf( stderr, "Loaded %d duplication regions.\n", dup_count);

	dups = ( Interval*) getMem( sizeof( Interval) * dup_count);
	rewind( dupsFile);

	dup_count = 0;

	while( !feof( dupsFile))
	{
		res = fscanf( dupsFile,"%s\t%d\t%d\n",chrom, &start, &end);
		dups[dup_count].chroName = (char *) getMem( sizeof( char) * ( strlen( chrom) + 1));
		strcpy( dups[dup_count].chroName, chrom);
		dups[dup_count].start = start;
		dups[dup_count].end = end;
		dup_count++;
	}

	qsort( dups, dup_count, sizeof( Interval), vh_compareCharInterval);

	fclose( dupsFile);
}

void vh_writeMeiIndices( ref_genome* ref)
{
	int i, index, index_tmp, old_index;
	char* chrom;

	/* Make sure the mei table is sorted according to chromosome name */
	qsort( g_meiTable, selectedMeiTypeSize, sizeof(Mei), compareChrNameMeiTable);

	mInd = ( MeiIndex*) getMem( ref->chrom_count * sizeof( MeiIndex));
	for( i = 0; i < ref->chrom_count; i++)
	{
		mInd[i].chroName = NULL;
		mInd[i].start = -1;
		mInd[i].end = -1;
	}

	/* I suppose that the first one is in the reference, will fix it */
	index = find_chr_index_ref( ref, g_meiTable[0].chroName);

	mInd[index].start = 0;
	mInd[index].chroName = ( char*) getMem( sizeof( char) * strlen( g_meiTable[0].chroName));
	strcpy( mInd[index].chroName, g_meiTable[0].chroName);

	old_index = index;

	for( i = 1; i < selectedMeiTypeSize; i++)
	{
		if( strcmp( g_meiTable[i].chroName, mInd[index].chroName))
		{
			index_tmp = find_chr_index_ref( ref, g_meiTable[i].chroName);
			if( index_tmp == -1)
				continue;
			index = index_tmp;
			mInd[index].start = i;
			mInd[old_index].end = i - 1;
			mInd[index].chroName = ( char*) getMem( sizeof(char) * strlen( g_meiTable[i].chroName));
			strcpy( mInd[index].chroName, g_meiTable[i].chroName);
			old_index = index;
		}
	}
	mInd[old_index].end = i;

	//for(i=0; i<ref->chrom_count;i++)
		//fprintf(stderr,"%d %s %d %d\n", i, mInd[i].chroName, mInd[i].start, mInd[i].end);
}

void vh_readRepeatTable (char *repeatFileName, char *meiStr)
{
	int satellites, i,j;
	char chrom[STRMAX], strand[2], type[STRMAX], class[STRMAX];
	int start, end;
	int res, chr=0;
	char **selectedMeiType;
	FILE *repeatFile = safe_fopen (repeatFileName, "r");
	if (repeatFile == NULL)
	{
		//TODO: handle me
		fprintf (stderr, "Cannot open file %s\n", repeatFileName);
		exit (-1);
	}

	/* need to implement file format check */
	satellites = 0;
	selectedMeiTypeSize=0;
	selectedMeiType=vh_parseMeiString(meiStr);

	while (!feof (repeatFile))
	{
		res = fscanf(repeatFile, "%s\t%d\t%d\t%s\t%s\t%s\n", chrom, &start, &end, strand, type, class);

		if (feof (repeatFile))
			break;

		if (res != 6)
		{
			// todo: global exit function
			fprintf(stderr, "The repeats file %s is not in the correct BED6 format. Consult the manual on how to generate a valid file.\n", repeatFileName);
			exit(-1);
		}

		if (strstr(class, "Satellite"))
			satellites++;
		else
		{
			for(i=0;i<meiTypeSize;i++)
			{
				if(strcasestr(type,selectedMeiType[i]))
				{
					selectedMeiTypeSize++;
					break;
				}
			}
		}
	}
	if (satellites == 0)
	{
		fprintf(stderr, "No satellites found\n");
		fprintf(logFile, "No satellite regions found.\n");
	}
	else
	{
		fprintf(stderr, "Loaded %d satellite regions.\n", satellites);
		fprintf(logFile, "There are %d satellite regions.\n", satellites);
	}

	repeats = (Interval *) getMem(sizeof(Interval) * satellites);
	g_meiTable = (Mei *) getMem(sizeof(Mei) * selectedMeiTypeSize);
	rewind(repeatFile);

	satellites = 0;
	selectedMeiTypeSize=0;
	while (!feof (repeatFile))
	{
		res = fscanf(repeatFile, "%s\t%d\t%d\t%s\t%s\t%s\n", chrom, &start, &end, strand, type, class);
		if (strstr(class, "Satellite"))
		{
			repeats[satellites].chroName = (char *) getMem(sizeof(char) * (strlen(chrom)+1));
			strcpy(repeats[satellites].chroName, chrom);
			repeats[satellites].start = start;
			repeats[satellites].end = end;
			if (g_maxRepeatLength < repeats[satellites].end - repeats[satellites].start)
				g_maxRepeatLength =  repeats[satellites].end - repeats[satellites].start;
			satellites++;
		}
		else
		{
			for(i=0;i<meiTypeSize;i++)
			{
				if(strcasestr(type,selectedMeiType[i]))
				{
					g_meiTable[selectedMeiTypeSize].chroName= (char *) getMem(sizeof(char) * (strlen(chrom)+1));
					g_meiTable[selectedMeiTypeSize].subclass=(char *) getMem(sizeof(char) * (strlen(type)+1));
					g_meiTable[selectedMeiTypeSize].superclass=(char *) getMem(sizeof(char) * (strlen(class)+1));
					strcpy(g_meiTable[selectedMeiTypeSize].chroName, chrom);
					strcpy(g_meiTable[selectedMeiTypeSize].subclass, type);
					strcpy(g_meiTable[selectedMeiTypeSize].superclass, selectedMeiType[i]);
					strcpy(g_meiTable[selectedMeiTypeSize].strand,strand);
					g_meiTable[selectedMeiTypeSize].start = start;
					g_meiTable[selectedMeiTypeSize].end = end;
					//fprintf(stderr,"%s %s %d %d\n",g_meiTable[selectedMeiTypeSize].chroName, g_meiTable[selectedMeiTypeSize].superclass,g_meiTable[selectedMeiTypeSize].start, g_meiTable[selectedMeiTypeSize].end);
					selectedMeiTypeSize++;
					break;
				}
			}
		}

		if (feof (repeatFile))
			break;

		if (res != 6)
			break;
	}
	if (selectedMeiTypeSize == 0)
	{
		fprintf(stderr, "No mobile element found.\n");
		fprintf(logFile, "No mobile element regions found.\n");
	}
	else
	{
		fprintf(stderr, "Loaded %d mobile element regions.\n", selectedMeiTypeSize);
		fprintf(logFile, "There are %d mobile element regions.\n", selectedMeiTypeSize);
	}

	repeat_count = satellites;

	qsort (repeats, repeat_count, sizeof (Interval),vh_compareIntInterval);

	fclose (repeatFile);
}

int vh_binarySearchInterval (DivetRow * newRow, int pos)	// return 1 if one of the ends falls inside a repeat
{
	int minId = 0, maxId = repeat_count, midId, i;
	int stopId, startId;
	int posT = pos - g_maxRepeatLength;
	if (posT < 0)
		posT = 0;
	midId = (minId + maxId) / 2;
	while (!(minId > maxId - 3 || repeats[minId].start == pos|| repeats[maxId].start == pos || repeats[midId].start == pos))
	{
		midId = (minId + maxId) / 2;
		if (repeats[midId].start > pos)
		{
			maxId = midId;

		}
		else if (repeats[midId].start < pos)
		{
			minId = midId;
		}
	}
	stopId = maxId + 23;
	if (stopId > repeat_count)
		stopId = repeat_count;
	startId = 0;
	minId = 0;
	maxId = repeat_count;
	midId = (minId + maxId) / 2;
	while (!(minId > maxId - 3 || repeats[minId].end == posT || repeats[maxId].end == posT || repeats[midId].end == posT))
	{
		midId = (minId + maxId) / 2;
		if (repeats[midId].end > posT)
		{
			maxId = midId;
		}
		else if (repeats[midId].end < posT)
		{
			minId = midId;
		}
	}

	startId = minId - 23;
	if (startId < 0)
		startId = 0;
	for (i = startId; i < stopId; i++)
	{
		if (strcmp (repeats[i].chroName, newRow->chroName) == 0 && ((repeats[i].end > newRow->locMapLeftEnd
				&& repeats[i].start < newRow->locMapLeftEnd) || (repeats[i].end > newRow->locMapRightStart
						&& repeats[i].start < newRow->locMapRightStart)))
		{
			return 1;
		}
	}
	return 0;
}


int binarySearchInterval( char *chroName, int pos) // return true if one of the ends falls inside a repeat
{
	int minId = 0, maxId = repeat_count, midId, i;
	int stopId, startId;
	int posT = pos - g_maxRepeatLength;
	if (posT < 0)
		posT = 0;
	midId = (minId + maxId) / 2;
	while (!(minId > maxId - 3 || repeats[minId].start == pos|| repeats[maxId].start == pos || repeats[midId].start == pos))
	{
		midId = (minId + maxId) / 2;
		if (repeats[midId].start > pos)
		{
			maxId = midId;

		}
		else if (repeats[midId].start < pos)
		{
			minId = midId;
		}
	}
	stopId = maxId + 23;
	if (stopId > repeat_count)
		stopId = repeat_count;
	startId = 0;
	minId = 0;
	maxId = repeat_count;
	midId = (minId + maxId) / 2;
	while (!(minId > maxId - 3 || repeats[minId].end == posT || repeats[maxId].end == posT || repeats[midId].end == posT))
	{
		midId = (minId + maxId) / 2;
		if (repeats[midId].end > posT)
		{
			maxId = midId;
		}
		else if (repeats[midId].end < posT)
		{
			minId = midId;
		}
	}

	startId = minId - 23;
	if (startId < 0)
		startId = 0;

	for (i = startId; i < stopId; i++)
	{
		if (strcmp (repeats[i].chroName, chroName) == 0 && (repeats[i].end > pos && repeats[i].start < pos))
		{
			return 1;
		}
	}
	return 0;
}

int binarySearchInterval2( char *chroName, int pos, int leftEnd, int rightStart) // return true if one of the ends falls inside a repeat
{
	int minId = 0, maxId = repeat_count, midId, i;
	int stopId, startId;
	int posT = pos - g_maxRepeatLength;
	if (posT < 0)
		posT = 0;
	midId = (minId + maxId) / 2;
	while (!(minId > maxId - 3 || repeats[minId].start == pos|| repeats[maxId].start == pos || repeats[midId].start == pos))
	{
		midId = (minId + maxId) / 2;
		if (repeats[midId].start > pos)
		{
			maxId = midId;

		}
		else if (repeats[midId].start < pos)
		{
			minId = midId;
		}
	}
	stopId = maxId + 23;
	if (stopId > repeat_count)
		stopId = repeat_count;
	startId = 0;
	minId = 0;
	maxId = repeat_count;
	midId = (minId + maxId) / 2;
	while (!(minId > maxId - 3 || repeats[minId].end == posT || repeats[maxId].end == posT || repeats[midId].end == posT))
	{
		midId = (minId + maxId) / 2;
		if (repeats[midId].end > posT)
		{
			maxId = midId;
		}
		else if (repeats[midId].end < posT)
		{
			minId = midId;
		}
	}

	startId = minId - 23;
	if (startId < 0)
		startId = 0;
	for (i = startId; i < stopId; i++)
	{
		if (strcmp (repeats[i].chroName, chroName) == 0 && ((repeats[i].end > leftEnd
				&& repeats[i].start < leftEnd) || (repeats[i].end > rightStart
						&& repeats[i].start < rightStart)))
			return 1;
	}
	return 0;
}

int notInRepeat2 (char* chroName, int locMapLeftEnd, int locMapRightEnd)
{
	if (binarySearchInterval2 (chroName, locMapLeftEnd, locMapLeftEnd, locMapRightEnd) || binarySearchInterval2 (chroName, locMapRightEnd, locMapLeftEnd, locMapRightEnd))
		return 0;
	else
		return 1;
}

int notInRepeat( char* chroName, int pos)
{
	if( binarySearchInterval( chroName, pos))
		return 0;
	else
		return 1;
}


int vh_notInRepeat (DivetRow * newRow)
{
	if (vh_binarySearchInterval (newRow, newRow->locMapLeftEnd) || vh_binarySearchInterval (newRow, newRow->locMapRightStart))
		return 0;
	else
		return 1;
}


int meiIntervalSearch( ref_genome* ref, char* chr_name, int read_start, int read_end)
{
	int start, end, center, i, chr_index;

	chr_index = find_chr_index_ref( ref, chr_name);
	if( chr_index == -1 || mInd[chr_index].chroName == NULL)
		return 0;

	start = mInd[chr_index].start;
	end = mInd[chr_index].end;
	center = ( start + end) / 2;

	while( start <= end)
	{
		if( g_meiTable[center].start <= read_start && g_meiTable[center].end >= read_end)
			return center;
		else if( g_meiTable[center].start < read_start)
			start = center + 1;
		else
			end = center - 1;

		center = ( start + end) / 2;
	}
	if( start > 0)
		start--;
	end++;

	for( i = start; i < end; i++)
	{
		if( g_meiTable[i].start <= read_start && g_meiTable[i].end >= read_end)
			return i;
	}
	return 0;
}

char* meiIntervalSearch2( ref_genome* ref, char *chroName, int pos) // return true if one of the ends falls inside a mei
{
	int minId, maxId, midId, i, chr_index;

	chr_index = find_chr_index_ref( ref, chroName);
	if( chr_index == -1 || mInd[chr_index].chroName == NULL)
		return NULL;

	minId = mInd[chr_index].start;
	maxId = mInd[chr_index].end;
	midId = ( minId + maxId) / 2;

	while (!(minId > maxId - 3 || g_meiTable[minId].start == pos|| g_meiTable[maxId].start == pos || g_meiTable[midId].start == pos))
	{
		midId = (minId + maxId) / 2;
		//if( g_meiTable[midId].start <= pos && g_meiTable[midId].end >= pos)
		//return g_meiTable[midId].superclass;
		if( g_meiTable[midId].start < pos)
			minId = midId + 1;
		else
			maxId = midId - 1;

		midId = ( minId + maxId) / 2;
	}
	if( minId > 0)
		minId--;
	maxId++;

	for( i = minId; i < maxId; i++)
	{
		if( strcmp (g_meiTable[i].chroName, chroName) == 0 && (g_meiTable[i].end > pos && g_meiTable[i].start < pos))
			return g_meiTable[i].subclass;
	}
	return NULL;
}

int meiIntervalSearch_Span( ref_genome* ref, char* chr_name, int read_start,int read_end)
{
	int start, end, center, i, chr_index;
	chr_index = find_chr_index_ref( ref, chr_name);
	if( chr_index == -1 || mInd[chr_index].chroName == NULL)
		return 0;
	start = mInd[chr_index].start;
	end = mInd[chr_index].end;
	center = ( start + end) / 2;
	while( (end - start) > 3 && g_meiTable[center].start != read_start && g_meiTable[start].start != read_start && g_meiTable[end].start != read_start)
	{
		if( g_meiTable[center].start < read_start)
		{
			start = center;
			center = (start+end)/2;
		}
		else if( g_meiTable[center].start > read_start)
		{
			end = center;
			center = (start+end)/2;
		}
	}
	if( start > 0)
		start--;
	end++;
	for( i = start; i < end; i++)
	{
		if( ( read_start <= g_meiTable[i].start) && ( read_end >= g_meiTable[i].end) && ( ( read_end - read_start) <= ( g_meiTable[i].end - g_meiTable[i].start + 500)))
		{
			return i;
		}
	}
	return 0;
}


int mei_filtering( ref_genome* ref)
{
	int i, mei_count = 0, len, chr = 0, indL, indR;
	char *svtype;
	char* chroName = NULL;

	LibraryInfo *libInfo;
	DivetRow *divetReadMappingPtr;
	libInfo = g_libInfo;

	while ( libInfo != NULL)
	{
		divetReadMappingPtr = libInfo->head;
		while ( divetReadMappingPtr != NULL)
		{
			/* Check the mei table whether the right or left read is inside a mobile element */
			indL = meiIntervalSearch( ref, divetReadMappingPtr->chroName, divetReadMappingPtr->locMapLeftStart, divetReadMappingPtr->locMapLeftEnd);
			indR = meiIntervalSearch( ref, divetReadMappingPtr->chroName, divetReadMappingPtr->locMapRightStart, divetReadMappingPtr->locMapRightEnd);
			if( indL != 0)
			{
				divetReadMappingPtr->svType = 'X';
				len = strlen( g_meiTable[indL].superclass);
				divetReadMappingPtr->meiType = ( char *) getMem ( sizeof ( char) * len + 1);
				strcpy( divetReadMappingPtr->meiType, g_meiTable[indL].superclass);

				len = strlen( g_meiTable[indL].subclass);
				divetReadMappingPtr->mei_subclass = ( char *) getMem ( sizeof ( char) * len + 1);
				strcpy( divetReadMappingPtr->mei_subclass, g_meiTable[indL].subclass);

				/*IF THE MEI INSERT IS + STRAND UPPERCASE IF - STRAND LOWER CASE*/
				if ( ( divetReadMappingPtr->orientationLeft == 'F' && g_meiTable[indL].strand[0] == 'C')
						|| ( divetReadMappingPtr->orientationLeft == 'F' && g_meiTable[indL].strand[0] == '-')
						|| ( divetReadMappingPtr->orientationLeft == 'R' && g_meiTable[indL].strand[0] == '+'  ))
				{
					if ( divetReadMappingPtr->meiType[0] == 'A')
						divetReadMappingPtr->meiType[0] = 'a';
					if ( divetReadMappingPtr->meiType[0] == 'L')
						divetReadMappingPtr->meiType[0] = 'l';
					if ( divetReadMappingPtr->meiType[0] == 'S')
						divetReadMappingPtr->meiType[0] = 's';
				}
				mei_count++;
			}
			if( indR != 0)
			{
				divetReadMappingPtr->svType = 'M';
				len = strlen( g_meiTable[indR].superclass);
				divetReadMappingPtr->meiType = ( char *) getMem ( sizeof ( char) * len + 1);
				strcpy( divetReadMappingPtr->meiType, g_meiTable[indR].superclass);

				len = strlen( g_meiTable[indR].subclass);
				divetReadMappingPtr->mei_subclass = ( char *) getMem ( sizeof ( char) * len + 1);
				strcpy( divetReadMappingPtr->mei_subclass, g_meiTable[indL].subclass);

				/*IF THE MEI INSERT IS + STRAND UPPERCASE IF - STRAND LOWWER CASE*/
				if ( ( divetReadMappingPtr->orientationRight == 'F' && g_meiTable[indR].strand[0] == 'C')
						|| ( divetReadMappingPtr->orientationRight == 'F' && g_meiTable[indR].strand[0] == '-')
						|| ( divetReadMappingPtr->orientationRight == 'R' && g_meiTable[indR].strand[0] == '+'  ))
				{
					if ( divetReadMappingPtr->meiType[0] == 'A')
						divetReadMappingPtr->meiType[0] = 'a';
					if ( divetReadMappingPtr->meiType[0] == 'L')
						divetReadMappingPtr->meiType[0] = 'l';
					if ( divetReadMappingPtr->meiType[0] == 'S')
						divetReadMappingPtr->meiType[0] = 's';
				}
				mei_count++;
			}
			divetReadMappingPtr = divetReadMappingPtr->next;
		}
		libInfo = libInfo->next;
	}
	return mei_count;
}


int dups_filtering( char* chr, int start, int end, float sim)
{
	int length_dup, length_sv=0;
	int i = 0;
	float tmp;

	if( start > end)
		return 0;

	/* Find the chromosome in the table */
	while( strcmp( dups[i].chroName, chr))
		i++;

	while( !strcmp( dups[i].chroName, chr))
	{
		length_dup = dups[i].end - dups[i].start;

		if( start > dups[i].start && end < dups[i].end)
		{
			if( ( end - start) > length_sv)
				length_sv = end - start;
		}
		else if( start > dups[i].start && start < dups[i].end && end > dups[i].end)
		{
			if( ( dups[i].end - start) > length_sv)
				length_sv = dups[i].end - start;
		}
		else if( start < dups[i].start && end < dups[i].end && end > dups[i].start)
		{
			if( ( end - dups[i].start) > length_sv)
				length_sv = end-dups[i].start;
		}
		else if( start < dups[i].start && end > dups[i].end)
		{
			if( ( dups[i].end - dups[i].start) > length_sv)
				length_sv = dups[i].end - dups[i].start;
		}
		i++;

		tmp = length_dup * sim;
		if( tmp <= length_sv && length_sv!=0)
		{
			return 1;
		}
		if(i == dup_count)
			return 0;
	}
	return 0;
}

