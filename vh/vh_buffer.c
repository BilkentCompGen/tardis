#include <stdio.h>
#include <stdlib.h>
#include "vh_buffer.h"


clusterInBuffer listClusterInBuffer[maxSizeOfBuffer];
int countInBuffer;
float maxScoreInBuffer;


int bufferIsUseful( ref_genome* ref, parameters *params)
{
	float bestWeight = inf;
	int bestWeightId = -1;
	int count;

	for( count = 0; count < countInBuffer; count++)
	{
		if( listClusterEl[listClusterInBuffer[count].clusterId].oldBestIsGood == 0)
		{
			listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore = calWeight( ref, params, listClusterInBuffer[count].clusterId, listClusterEl[listClusterInBuffer[count].clusterId].bestReadToRemove );
			listClusterEl[listClusterInBuffer[count].clusterId].oldBestIsGood = 1;
			listClusterInBuffer[count].score=listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore;
		}
	}

	for( count = 0; count <  countInBuffer; count++)
	{
		if( bestWeight > listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore)
		{
			bestWeight = listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore;
			bestWeightId = count;
		}
	}

	if( bestWeight <= maxScoreInBuffer && bestWeight != inf)
		return 1;
	else
		return 0;
}

int bestFromBuffer()
{
	int bestSet;
	float bestSetScore=10000000;
	int count;
	for (count=0; count<countInBuffer; count++)
	{
		if (bestSetScore>listClusterInBuffer[count].score)
		{
			bestSetScore=listClusterInBuffer[count].score;
			bestSet=count;
		}
	}
	return listClusterInBuffer[bestSet].clusterId;
}


void emptyBuffer()
{
	countInBuffer = 0;
	maxScoreInBuffer = 0;
}

int addToBuffer(float score, int clusterId)
{
	if( countInBuffer < maxSizeOfBuffer)
	{
		listClusterInBuffer[countInBuffer].score = score;
		listClusterInBuffer[countInBuffer].clusterId = clusterId;
		countInBuffer++;
		if(score > maxScoreInBuffer)
			maxScoreInBuffer = score;
		return 0;
	}
	int count;
	int count2;
	if( score < maxScoreInBuffer)
	{
		for( count = 0; count < maxSizeOfBuffer; count++)
		{
			if( listClusterInBuffer[count].score == maxScoreInBuffer)
			{
				listClusterInBuffer[count].score = score;
				listClusterInBuffer[count].clusterId = clusterId;
				//listClusterInBuffer[count].valid=true;
				maxScoreInBuffer = 0;

				for( count2 = 0; count2 < maxSizeOfBuffer; count2++)
				{
					if( listClusterInBuffer[count2].score > maxScoreInBuffer)
						maxScoreInBuffer = listClusterInBuffer[count2].score;
				}
				return 0;
			}
		}
	}
}
