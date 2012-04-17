/*
 * =====================================================================================
 *
 *       Filename:  TGM_Cluster.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/25/2011 14:25:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_CLUSTER_H
#define  TGM_CLUSTER_H

#include "TGM_ReadPairAttrbt.h"

typedef struct TGM_ClusterElmnt
{
    double mean[2];

    double std[2];

    double median[2];

    double min[2];

    double max[2];

    unsigned int numReadPair;

    unsigned int startIndex;

}TGM_ClusterElmnt;

typedef struct TGM_ClusterElmntArray
{
    TGM_ClusterElmnt* data;

    unsigned int size;

    unsigned int length;

    unsigned int capacity;

}TGM_ClusterElmntArray;

typedef struct TGM_Cluster
{
    const TGM_ReadPairAttrbtArray* pAttrbtArray;

    TGM_ClusterElmntArray* pElmntArray;

    int* pNext;

    int* pMap;
    
    int* pCount;

    unsigned int size;

    unsigned int capacity;

    int minReadPairNum;

    double minStd[2];

}TGM_Cluster;

TGM_Cluster* TGM_ClusterAlloc(int minReadPairNum);

void TGM_ClusterFree(TGM_Cluster* pCluster);

void TGM_ClusterSetMinStd(TGM_Cluster* pCluster, double minStd[2]);

void TGM_ClusterInit(TGM_Cluster* pCluster, const TGM_ReadPairAttrbtArray* pAttrbtArray);

void TGM_ClusterMake(TGM_Cluster* pCluster);

void TGM_ClusterBuild(TGM_Cluster* pCluster);

void TGM_ClusterConnect(TGM_Cluster* pCluster, unsigned int centerIndex, unsigned int memberIndex);

void TGM_ClusterFinalize(TGM_Cluster* pCluster);

int TGM_ClusterClean(TGM_Cluster* pCluster);

void TGM_ClusterPrint(const TGM_Cluster* pCluster, int32_t refID, const char* specialID);

int TGM_ClusterMerge(TGM_Cluster* pCluster, int dstIndex, int srcIndex);

#endif  /*TGM_CLUSTER_H*/
