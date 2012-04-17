/*
 * =====================================================================================
 *
 *       Filename:  TGM_Cluster.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  11/25/2011 23:01:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (),
 *        Company:
 *
 * =====================================================================================
 */

#include <math.h>

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Cluster.h"
#include "TGM_Utilities.h"

#define DEFAULT_CLUSTER_SIZE 20

KHASH_MAP_INIT_INT(clusterMap, int);


static void TGM_ClusterUpdateElmnt(TGM_Cluster* pCluster, int clusterID, int attrbtID)
{
    TGM_ClusterElmnt* pElmnt = TGM_ARRAY_GET_PT(pCluster->pElmntArray, clusterID);
    TGM_ReadPairAttrbt* pAttrbt = TGM_ARRAY_GET_PT(pCluster->pAttrbtArray, attrbtID);

    ++pElmnt->numReadPair;

    pElmnt->mean[0] += pAttrbt->firstAttribute;
    pElmnt->mean[1] += pAttrbt->secondAttribute;

    pElmnt->std[0] += pAttrbt->firstAttribute * pAttrbt->firstAttribute;
    pElmnt->std[1] += pAttrbt->secondAttribute * pAttrbt->secondAttribute;

    pElmnt->min[0] = (pAttrbt->firstAttribute < pElmnt->min[0] ? pAttrbt->firstAttribute : pElmnt->min[0]);
    pElmnt->min[1] = (pAttrbt->secondAttribute < pElmnt->min[1] ? pAttrbt->secondAttribute : pElmnt->min[1]);

    pElmnt->max[0] = (pAttrbt->firstAttribute > pElmnt->max[0] ? pAttrbt->firstAttribute : pElmnt->max[0]);
    pElmnt->max[1] = (pAttrbt->secondAttribute > pElmnt->max[1] ? pAttrbt->secondAttribute : pElmnt->max[1]);
}

static int TGM_ClusterAddElmnt(TGM_Cluster* pCluster, int attrbtID)
{
    int clusterID = pCluster->pElmntArray->size;

    if (TGM_ARRAY_IS_FULL(pCluster->pElmntArray))
    {
        TGM_ARRAY_RESIZE(pCluster->pElmntArray, pCluster->pElmntArray->capacity * 2, TGM_ClusterElmnt);
    }

    ++(pCluster->pElmntArray->size);

    TGM_ClusterElmnt* pElmnt = TGM_ARRAY_GET_PT(pCluster->pElmntArray, clusterID);
    TGM_ReadPairAttrbt* pAttrbt = TGM_ARRAY_GET_PT(pCluster->pAttrbtArray, attrbtID);

    pElmnt->numReadPair = 1;
    pElmnt->startIndex = attrbtID;

    pElmnt->mean[0] = pAttrbt->firstAttribute;
    pElmnt->mean[1] = pAttrbt->secondAttribute;

    pElmnt->std[0] = pAttrbt->firstAttribute * pAttrbt->firstAttribute;
    pElmnt->std[1] = pAttrbt->secondAttribute * pAttrbt->secondAttribute;

    pElmnt->min[0] = pAttrbt->firstAttribute;
    pElmnt->min[1] = pAttrbt->secondAttribute;

    pElmnt->max[0] = pAttrbt->firstAttribute;
    pElmnt->max[1] = pAttrbt->secondAttribute;

    return clusterID;
}

TGM_Cluster* TGM_ClusterAlloc(int minReadPairNum)
{
    TGM_Cluster* pCluster = (TGM_Cluster*) malloc(sizeof(TGM_Cluster));
    if (pCluster == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a cluster object.\n");

    TGM_ARRAY_ALLOC(pCluster->pElmntArray, DEFAULT_CLUSTER_SIZE, TGM_ClusterElmntArray, TGM_ClusterElmnt);

    pCluster->pAttrbtArray = NULL;
    pCluster->pCount = NULL;
    pCluster->pMap = NULL;
    pCluster->pNext = NULL;

    pCluster->size = 0;
    pCluster->capacity = 0;

    pCluster->minReadPairNum = minReadPairNum;
    pCluster->minStd[0] = -1;
    pCluster->minStd[1] = -1;

    return pCluster;
}

void TGM_ClusterFree(TGM_Cluster* pCluster)
{
    if (pCluster != NULL)
    {
        free(pCluster->pNext);
        free(pCluster->pMap);
        free(pCluster->pCount);

        TGM_ARRAY_FREE(pCluster->pElmntArray, TRUE);

        free(pCluster);
    }
}

void TGM_ClusterSetMinStd(TGM_Cluster* pCluster, double minStd[2])
{
    pCluster->minStd[0] = minStd[0];
    pCluster->minStd[1] = minStd[1];
}

void TGM_ClusterInit(TGM_Cluster* pCluster, const TGM_ReadPairAttrbtArray* pAttrbtArray)
{
    TGM_ARRAY_RESET(pCluster->pElmntArray);

    pCluster->pAttrbtArray = pAttrbtArray;
    pCluster->size = pAttrbtArray->size;

    if (pAttrbtArray->size > pCluster->capacity)
    {
        free(pCluster->pNext);
        free(pCluster->pMap);
        free(pCluster->pCount);

        pCluster->capacity = pAttrbtArray->size * 2;

        pCluster->pNext = (int*) malloc(pCluster->capacity * sizeof(int));
        if (pCluster->pNext == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of linked list array in a cluster object.\n");

        pCluster->pMap = (int*) malloc(pCluster->capacity * sizeof(int));
        if (pCluster->pMap == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the cluster map in a cluster object.\n");

        pCluster->pCount = (int*) malloc(pCluster->capacity * sizeof(int));
        if (pCluster->pCount == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the neighbour count in a cluster object.\n");
    }

    for (unsigned int i = 0; i != pCluster->size; ++i)
    {
        pCluster->pNext[i] = i;
        pCluster->pMap[i] = i;
        pCluster->pCount[i] = 0;
    }
}

void TGM_ClusterMake(TGM_Cluster* pCluster)
{
    unsigned int lastLowIndex = 0;
    for (unsigned int i = 0; i != pCluster->size; ++i)
    {
        unsigned int j = lastLowIndex;

        double firstBound = TGM_ReadPairAttrbtArrayGetFirstBound(pCluster->pAttrbtArray, i);
        double secondBound = TGM_ReadPairAttrbtArrayGetSecondBound(pCluster->pAttrbtArray, i);

        while ((j < pCluster->size) && (pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) > firstBound)
            ++j;

        lastLowIndex = j;

        while (fabs(pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) <= firstBound)
        {
            if (fabs(pCluster->pAttrbtArray->data[i].secondAttribute - pCluster->pAttrbtArray->data[j].secondAttribute) <= secondBound)
                ++(pCluster->pCount[i]);

            ++j;

            if (j == pCluster->size)
                break;
        }
    }
}

void TGM_ClusterBuild(TGM_Cluster* pCluster)
{
    unsigned int lastLowIndex = 0;
    for (unsigned int i = 0; i != pCluster->size; ++i)
    {
        unsigned int j = lastLowIndex;

        double firstBound = TGM_ReadPairAttrbtArrayGetFirstBound(pCluster->pAttrbtArray, i);
        double secondBound = TGM_ReadPairAttrbtArrayGetSecondBound(pCluster->pAttrbtArray, i);

        while ((j < pCluster->size) && (pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) > firstBound)
            ++j;

        lastLowIndex = j;

        unsigned int maxCount = 0;
        unsigned int centerIndex = 0;
        while (fabs(pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) <= firstBound)
        {
            if (fabs(pCluster->pAttrbtArray->data[i].secondAttribute - pCluster->pAttrbtArray->data[j].secondAttribute) <= secondBound)
            {
                if (pCluster->pCount[j] > maxCount)
                {
                    centerIndex = j;
                    maxCount = pCluster->pCount[j];
                }
            }

            ++j;

            if (j == pCluster->size)
                break;
        }

        TGM_ClusterConnect(pCluster, centerIndex, i);
    }
}

void TGM_ClusterConnect(TGM_Cluster* pCluster, unsigned int centerIndex, unsigned int memberIndex)
{
    if (centerIndex == memberIndex)
        return;

    int centerID = pCluster->pMap[centerIndex];
    int centerNext = pCluster->pNext[centerIndex];

    pCluster->pNext[centerIndex] = memberIndex;

    unsigned int i = memberIndex;
    for (; pCluster->pNext[i] != memberIndex; i = pCluster->pNext[i])
    {
        pCluster->pMap[i] = centerID;
    }

    pCluster->pNext[i] = centerNext;
    pCluster->pMap[i] = centerID;
}


void TGM_ClusterFinalize(TGM_Cluster* pCluster)
{
    // create a hash to transfer the cluster center ID to cluster element ID
    khash_t(clusterMap)* clusterHash = kh_init(clusterMap);
    kh_resize(clusterMap, clusterHash, DEFAULT_CLUSTER_SIZE);

    int khRet = 0;
    khiter_t khIter = 0;
    int clusterIndex = 0;

    for (unsigned int i = 0; i != pCluster->size; ++i)
    {
        // skip those very small clusters
        unsigned int clusterSize = pCluster->pCount[pCluster->pMap[i]];
        if (clusterSize < pCluster->minReadPairNum)
            continue;

        khIter = kh_put(clusterMap, clusterHash, pCluster->pMap[i], &khRet);
        if (khRet == 0)
        {
            clusterIndex = kh_value(clusterHash, khIter);
            TGM_ClusterUpdateElmnt(pCluster, clusterIndex, i);
        }
        else
        {
            clusterIndex = TGM_ClusterAddElmnt(pCluster, i);
            kh_value(clusterHash, khIter) = clusterIndex;
        }

        pCluster->pMap[i] = clusterIndex;
    }

    kh_destroy(clusterMap, clusterHash);

    for (unsigned int i = 0; i != pCluster->pElmntArray->size; ++i)
    {
        TGM_ClusterElmnt* pElmnt = TGM_ARRAY_GET_PT(pCluster->pElmntArray, i);

        for (unsigned int j = 0; j != 2; ++j)
        {
            pElmnt->mean[j] /= pElmnt->numReadPair;
            pElmnt->std[j] = sqrt(pElmnt->std[j] / pElmnt->numReadPair - (pElmnt->mean[j] * pElmnt->mean[j]));
        }
    }
}

int TGM_ClusterClean(TGM_Cluster* pCluster)
{
    unsigned int oldSize = pCluster->pElmntArray->size;
    pCluster->pElmntArray->length = pCluster->pElmntArray->size;

    // lazy deletion
    // the deleted elements' number of read pairs will be set to zero
    for (unsigned int i = 0; i != oldSize; ++i)
    {
        if (pCluster->pElmntArray->data[i].numReadPair < pCluster->minReadPairNum
            || pCluster->pElmntArray->data[i].std[0] <= pCluster->minStd[0]
            || pCluster->pElmntArray->data[i].std[1] <= pCluster->minStd[1]
            || pCluster->pElmntArray->data[i].std[0] != pCluster->pElmntArray->data[i].std[0]
            || pCluster->pElmntArray->data[i].std[1] != pCluster->pElmntArray->data[i].std[1])
        {
            pCluster->pElmntArray->data[i].numReadPair = 0;
            --(pCluster->pElmntArray->size);
        }
    }

    return oldSize;
}

void TGM_ClusterPrint(const TGM_Cluster* pCluster, int32_t refID, const char* specialID)
{
    const TGM_ClusterElmntArray* pArray = pCluster->pElmntArray;

    unsigned int i = 0;
    unsigned int count = 0;
    while (count != pArray->size)
    {
        int firstAttribute = pCluster->pAttrbtArray->data[pArray->data[i].startIndex].firstAttribute;
        printf("chr%d\t%d\t%d\t%d\t%f\n", refID + 1, firstAttribute, firstAttribute + 1, pArray->data[i].numReadPair, pArray->data[i].std[0]);

        if (pArray->data[i].numReadPair != 0)
            ++count;

        ++i;
    }

    printf("\n");
}

int TGM_ClusterMerge(TGM_Cluster* pCluster, int dstIndex, int srcIndex)
{
    TGM_ClusterElmnt* pSrcElmnt = pCluster->pElmntArray->data + srcIndex;
    TGM_ClusterElmnt* pDstElmnt = pCluster->pElmntArray->data + dstIndex;

    TGM_ClusterConnect(pCluster, pDstElmnt->startIndex, pSrcElmnt->startIndex);

    for (unsigned int i = 0; i != 2; ++i)
    {
        pDstElmnt->min[i] = (pDstElmnt->min[i] < pSrcElmnt->min[i] ? pDstElmnt->min[i] : pSrcElmnt->min[i]);
        pDstElmnt->max[i] = (pDstElmnt->max[i] > pSrcElmnt->max[i] ? pDstElmnt->max[i] : pSrcElmnt->max[i]);

        pDstElmnt->std[i] = (pDstElmnt->std[i] * pDstElmnt->std[i] + pDstElmnt->mean[i] + pDstElmnt->mean[i]) * pDstElmnt->numReadPair;
        pSrcElmnt->std[i] = (pSrcElmnt->std[i] * pSrcElmnt->std[i] + pSrcElmnt->mean[i] + pSrcElmnt->mean[i]) * pSrcElmnt->numReadPair;

        pDstElmnt->std[i] += pSrcElmnt->std[i];

        pDstElmnt->mean[i] = (pDstElmnt->mean[i] * pDstElmnt->numReadPair + pSrcElmnt->mean[i] * pSrcElmnt->numReadPair) / (pDstElmnt->numReadPair + pSrcElmnt->numReadPair);
    }

    pDstElmnt->numReadPair += pSrcElmnt->numReadPair;
    pSrcElmnt->numReadPair = 0;
    --(pCluster->size);

    pDstElmnt->std[0] = sqrt(pDstElmnt->std[0] / pDstElmnt->numReadPair - pDstElmnt->mean[0] * pDstElmnt->mean[0]);
    pDstElmnt->std[1] = sqrt(pDstElmnt->std[1] / pDstElmnt->numReadPair - pDstElmnt->mean[1] * pDstElmnt->mean[1]);

    return pDstElmnt->numReadPair;
}
