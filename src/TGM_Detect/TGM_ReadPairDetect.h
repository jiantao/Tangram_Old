/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairDetect.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/23/2012 01:38:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_READPAIRDETECT_H
#define  TGM_READPAIRDETECT_H

#include <pthread.h>

#include "TGM_Cluster.h"
#include "TGM_ReadPairBuild.h"

typedef struct TGM_ReadPairDetectPars
{
    char* workingDir;

    int numThreads;

    int workingRefID[2];

    int minNumClustered;

    int minEventLength;

}TGM_ReadPairDetectPars;

typedef struct TGM_SpecialID
{
    char (*names)[3];
    
    uint32_t size;

}TGM_SpecialID;

typedef struct TGM_DetectJob
{
    void* data;

    int32_t refID;

    SV_EventType eventType;

}TGM_DetectJob;

typedef struct TGM_DetectJobStack
{
    TGM_DetectJob* pJobs;

    const TGM_LibInfoTable* pLibTable;

    unsigned int top;

    unsigned int capacity;

    pthread_mutex_t mutex;

    unsigned int numThread;

    unsigned int loadingNum;

    unsigned int loadingLimit;

    int nextJobID;

}TGM_DetectJobStack;

typedef struct TGM_DelEvent
{
    int32_t clusterID;

    int32_t refID;
    uint32_t pos;
    uint32_t end;
    uint32_t length;

    int pos5[3];
    int pos3[3];

    int CIpos[2];
    int CIend[2];
    int CIlen[2];

    unsigned char quality;
    unsigned char mapQ5;
    unsigned char mapQ3;

}TGM_DelEvent;

typedef struct TGM_DelArray
{
    TGM_DelEvent* data;

    unsigned int size;

    unsigned int capacity;

}TGM_DelArray;

typedef struct TGM_SpecialEvent
{
    int32_t refID;
    uint32_t pos;
    uint32_t posUncertainty;
    uint32_t length;

    int pos5[2];
    int pos3[2];

    int clusterID[2];
    int numFrag[2];
    double sense[2];

}TGM_SpecialEvent;

typedef struct TGM_SpecialEventArray
{
    TGM_SpecialEvent* data;

    unsigned int size;

    unsigned int capacity;

}TGM_SpecialEventArray;

typedef struct SV_AssistArray
{
    int* pFragLenDiff;

    int* pMapQ5;

    int* pMapQ3;

    unsigned int size;

    unsigned int capacity;

}SV_AssistArray;

TGM_DetectJobStack* TGM_DetectJobStackAlloc(unsigned int numThread, unsigned int loadingLimit, const TGM_LibInfoTable* pLibTable);

void TGM_DetectJobStackFree(TGM_DetectJobStack* pJobStack);

SV_AssistArray* SV_AssistArrayAlloc(void);

void SV_AssistArrayFree(SV_AssistArray* pAssistArray);


void TGM_ReadPairDetect(const TGM_ReadPairDetectPars* pDetectPars);

int TGM_DetectJobStackLock(TGM_DetectJobStack* pJobStack);

int TGM_DetectJobStackUnlock(TGM_DetectJobStack* pJobStack);

static inline TGM_Bool TGM_DetectJobStackIsEmpty(const TGM_DetectJobStack* pJobStack)
{
    return (pJobStack->top == 0);
}

void TGM_DetectJobStackPush(TGM_DetectJobStack* pJobStack, const TGM_DetectJob* pJob);

TGM_Status TGM_DetectJobStackPop(TGM_DetectJob* pJob, TGM_DetectJobStack* pJobStack);

TGM_SpecialID* TGM_SpecialIDRead(FILE* libInput);

void TGM_SpecialIDFree(TGM_SpecialID* pSpecialID);

void TGM_LocalPairArrayRead(TGM_LocalPairArray* pLongPairArray, FILE* input);

void TGM_CrossPairArrayRead(TGM_CrossPairArray* pCrossPairArray, FILE* input);

void TGM_SpecialPairArrayRead(TGM_SpecialPairArray* pSpecialPairArray, FILE* input);

// void TGM_SpecialPairTableReadID(TGM_SpecialPairTable* pSpeicalPairTable, FILE* libInput);

void SV_AssistArrayResize(SV_AssistArray* pAssistArray, unsigned int newSize);

void TGM_DetectDel(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable);

void TGM_DetectTandemDup(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable);

void TGM_DetectInversion(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable);

void TGM_DetectSpecial(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable, const TGM_SpecialID* pSpecialID);

void TGM_SpecialEventMake(TGM_SpecialEvent* pSpecialEvent, const TGM_Cluster* pCluster, unsigned int index, 
                         const TGM_SpecialPairArray* pSpecialPairArray, const TGM_LibInfoTable* pLibTable);

void TGM_SpecialEventMerge(TGM_SpecialEvent* mergedEvent, const TGM_SpecialEvent* pHeadEvent, const TGM_SpecialEvent* pTailEvent, TGM_Cluster* pCluster3, TGM_Cluster* pCluster5);

void TGM_SpecialEventPrint(const TGM_SpecialEvent* pSpecialEvent);

void TGM_DetectTranslocation(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable);

void TGM_ReadPairFindDel(TGM_DelArray* pDelArray, SV_AssistArray* pAssistArray, const TGM_LocalPairArray* pLongPairArray,
                        TGM_Cluster* pDelCluster, const TGM_LibInfoTable* pLibTable, const TGM_ReadPairDetectPars* pPars);

/*  
void TGM_DelEventGenotype(TGM_DelArray* pDelArray, const TGM_LibInfoTable* pLibTable);
*/


#endif  /*TGM_READPAIRDETECT_H*/
