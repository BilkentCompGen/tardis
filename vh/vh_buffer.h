#ifndef __BUFFER__
#define __BUFFER__

#include "vh_setcover.h"

#define maxSizeOfBuffer 100

/*Buffer keeps a list of top score (minimim score) picked from previous calculation of
  scores for each cluster (in an array named listClusterInBuffer). 
  We also keep the largest score in the listClusterInBuffer as upperbound (in a variable names maxScoreInBuffer), so
  if at some iteration all of the scores in listClusterInBuffer become greater than maxScoreInBuffer
  we need to recalculate the top set of clusters (i.e. recreate the array listClusterInBuffer).
  The function BufferIsUsefull is responsible to tell us if the buffer still contains clusters
  which have score less than upperbound maxScoreInBuffer. 
  The big assumption is that in each iteration the score of all other clusters
  will increase (become worse).

  After selection of each cluster we need to update the score of all
  the clusters which are effected, includeing the ones in listClusterInBuffer.
  The function bestFromBuffer() gives the clusterId which is the best in listClusterInBuffer.
  Not if none of the clusters in listClusterInBuffer are less than maxScoreInBuffer the 
  function bestFromBuffer() should not be used, and we need to recalculate the top clusters.

*/

typedef struct clusterInBuffer{
  int clusterId;// id of the cluster;
  float score;//the score of picking this cluster
}clusterInBuffer;


extern struct clusterInBuffer listClusterInBuffer[maxSizeOfBuffer];
extern int countInBuffer;
extern float maxScoreInBuffer;
int addToBuffer(float, int);
int bestFromBuffer();
void emptyBuffer();
int bufferIsUseful( ref_genome* ref, parameters *params);
#endif
