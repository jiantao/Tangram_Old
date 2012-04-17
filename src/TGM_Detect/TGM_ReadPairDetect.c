/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairDetect.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/23/2012 08:14:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <limits.h>

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Utilities.h"
#include "TGM_ReadPairDetect.h"
#include "TGM_ReadPairAttrbt.h"

#define DEFAULT_SV_CAPACITY 50

static const char* TGM_LibTableFileName = "lib_table.dat";

static const char* TGM_READ_PAIR_FILE_NAME_TEMPLATE[] = 
{
    "refXXX_long_pairs.dat",

    "refXXX_short_pairs.dat",

    "refXXX_reversed_pairs.dat",

    "refXXX_inverted_pairs.dat",

    "refXXX_special_pairs.dat",

    "refXXX_cross_pairs.dat",
};

KHASH_MAP_INIT_STR(name, uint32_t);

KHASH_SET_INIT_INT(readGrp);

static void TGM_DelEventMerge(TGM_DelArray* pDelArray, TGM_Cluster* pDelCluster)
{
    TGM_DelEvent* pLastEvent = pDelArray->data + (pDelArray->size - 1);
    TGM_DelEvent* pNewEvent = pDelArray->data + pDelArray->size;

    pLastEvent->pos5[0] = (pLastEvent->pos5[0] < pNewEvent->pos5[0] ? pLastEvent->pos5[0] : pNewEvent->pos5[0]);
    pLastEvent->pos5[1] = (pLastEvent->pos5[1] > pNewEvent->pos5[1] ? pLastEvent->pos5[1] : pNewEvent->pos5[1]);
    pLastEvent->pos5[2] = (pLastEvent->pos5[2] > pNewEvent->pos5[2] ? pLastEvent->pos5[2] : pNewEvent->pos5[2]);

    pLastEvent->pos3[0] = (pLastEvent->pos3[0] < pNewEvent->pos3[0] ? pLastEvent->pos3[0] : pNewEvent->pos3[0]);
    pLastEvent->pos3[1] = (pLastEvent->pos3[1] > pNewEvent->pos3[1] ? pLastEvent->pos3[1] : pNewEvent->pos3[1]);
    pLastEvent->pos3[2] = (pLastEvent->pos3[2] > pNewEvent->pos3[2] ? pLastEvent->pos3[2] : pNewEvent->pos3[2]);

    int numReadPair = TGM_ClusterMerge(pDelCluster, pLastEvent->clusterID, pNewEvent->clusterID);

    pLastEvent->pos = pLastEvent->pos5[1] + 1;
    pLastEvent->end = pLastEvent->pos3[0];
    pLastEvent->length = DoubleRoundToInt(pDelCluster->pElmntArray->data[pLastEvent->clusterID].mean[1]);
    pLastEvent->quality = (int) ((numReadPair * 100.0) / (numReadPair + 10.0));

    pLastEvent->CIpos[0] = -(pLastEvent->pos5[2] - pLastEvent->pos5[0]) / numReadPair;
    pLastEvent->CIpos[1] = (pLastEvent->pos3[2] - pLastEvent->pos3[0]) / numReadPair;

    pLastEvent->CIend[0] = -(pLastEvent->pos5[2] - pLastEvent->pos5[0]) / numReadPair;
    pLastEvent->CIend[1] = (pLastEvent->pos3[2] - pLastEvent->pos3[0]) / numReadPair;

    pLastEvent->CIlen[0] = 0;
    pLastEvent->CIlen[1] = 0;
}

static int CompareDelEvents(const void* pEvent1, const void* pEvent2)
{
    const TGM_DelEvent* pE1 = pEvent1;
    const TGM_DelEvent* pE2 = pEvent2;

    if (pE1->refID < pE2->refID)
        return -1;
    else if (pE1->refID > pE2->refID)
        return 1;
    else
    {
        if (pE1->pos < pE2->pos)
            return -1;
        else if (pE1->pos > pE2->pos)
            return 1;
        else
        {
            if (pE1->length < pE2->length)
                return -1;
            else if (pE1->length > pE2->length)
                return 1;
            else
            {
                if (pE1->quality < pE2->quality)
                    return -1;
                else if (pE1->quality > pE2->quality)
                    return 1;
            }
        }
    }

    return 0;
}

static int CompareSpecialEvents(const void* pEvent1, const void* pEvent2)
{
    const TGM_SpecialEvent* pE1 = pEvent1;
    const TGM_SpecialEvent* pE2 = pEvent2;

    if (pE1->refID < pE2->refID)
        return -1;
    else if (pE1->refID > pE2->refID)
        return 1;
    else
    {
        if (pE1->pos < pE2->pos)
            return -1;
        else if (pE1->pos > pE2->pos)
            return 1;
        else
        {
            if (pE1->length < pE2->length)
                return -1;
            else if (pE1->length > pE2->length)
                return 1;
            else
                return 0;
        }
    }

    return 0;
}

TGM_DetectJobStack* TGM_DetectJobStackAlloc(unsigned int numThread, unsigned int loadingLimit, const TGM_LibInfoTable* pLibTable)
{
    TGM_DetectJobStack* pJobStack = (TGM_DetectJobStack*) malloc(sizeof(TGM_DetectJobStack));
    if (pJobStack == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a detection job queue.\n");

    pJobStack->pJobs = (TGM_DetectJob*) calloc(sizeof(TGM_DetectJob), 2 * numThread);
    if (pJobStack->pJobs == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a detection job queue.\n");

    pJobStack->pLibTable = pLibTable;

    pJobStack->top = 0;
    pJobStack->capacity = 2 * numThread;

    int status = pthread_mutex_init(&(pJobStack->mutex), NULL);
    if (status != 0)
        TGM_ErrQuit("ERROR: Unable to initialize the mutex.\n");

    pJobStack->numThread = numThread;
    pJobStack->loadingNum = 0;
    pJobStack->loadingLimit = loadingLimit;

    pJobStack->nextJobID = 0;

    return pJobStack;
}

void TGM_DetectJobStackFree(TGM_DetectJobStack* pJobStack)
{
    if (pJobStack != NULL)
    {
        free(pJobStack->pJobs);
        pthread_mutex_destroy(&(pJobStack->mutex));

        free(pJobStack);
    }
}

int TGM_DetectJobStackLock(TGM_DetectJobStack* pJobStack)
{
    return (pthread_mutex_lock(&(pJobStack->mutex)));
}

int TGM_DetectJobStackUnlock(TGM_DetectJobStack* pJobStack)
{
    return (pthread_mutex_unlock(&(pJobStack->mutex)));
}

void TGM_DetectJobStackPush(TGM_DetectJobStack* pJobStack, const TGM_DetectJob* pJob)
{
    if (pJobStack->top == pJobStack->capacity)
    {
        pJobStack->capacity *= 2;
        pJobStack->pJobs = (TGM_DetectJob*) realloc(pJobStack->pJobs, sizeof(TGM_DetectJob) * pJobStack->capacity);
        if (pJobStack->pJobs == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the detection job queue.\n");
    }

    pJobStack->pJobs[pJobStack->top] = *pJob;
    ++(pJobStack->top);
}

TGM_Status TGM_DetectJobStackPop(TGM_DetectJob* pJob, TGM_DetectJobStack* pJobStack)
{
    if (TGM_DetectJobStackIsEmpty(pJobStack))
        return TGM_ERR;

    *pJob = pJobStack->pJobs[pJobStack->top - 1];
    --(pJobStack->top);

    return TGM_OK;
}

void TGM_ReadPairDetect(const TGM_ReadPairDetectPars* pDetectPars)
{
    // create a buffer to store the file name 
    int dirLen = strlen(pDetectPars->workingDir);
    char nameBuff[dirLen + 60];

    // create the library information file name
    strncpy(nameBuff, pDetectPars->workingDir, dirLen);
    strcpy(nameBuff + dirLen, TGM_LibTableFileName);

    // read the library information into memory
    FILE* pLibInput = fopen(nameBuff, "rb");
    if (pLibInput == NULL)
        TGM_ErrQuit("ERROR: Cannot open library information file: \"%s\".\n", nameBuff);

    TGM_LibInfoTable* pLibTable = TGM_LibInfoTableRead(pLibInput);


    uint32_t detectSet = 0;
    uint32_t readSize = fread(&detectSet, sizeof(uint32_t), 1, pLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read detect set.\n");

    for (unsigned int i = SV_DELETION; i != SV_INTER_CHR_TRNSLCTN; ++i)
    {
        switch(i)
        {
            case SV_DELETION:
                if ((detectSet & (1 << i)) != 0)
                    TGM_DetectDel(pDetectPars, pLibTable);
                break;
            case SV_TANDEM_DUP:
                if ((detectSet & (1 << i)) != 0)
                    TGM_DetectTandemDup(pDetectPars, pLibTable);
                break;
            case SV_INVERSION:
                if ((detectSet & (1 << i)) != 0)
                    TGM_DetectInversion(pDetectPars, pLibTable);
                break;
            case SV_SPECIAL:
                if ((detectSet & (1 << i)) != 0)
                {
                    TGM_SpecialID* pSpecialID = TGM_SpecialIDRead(pLibInput);
                    if (pSpecialID != NULL)
                        TGM_DetectSpecial(pDetectPars, pLibTable, pSpecialID);

                    TGM_SpecialIDFree(pSpecialID);
                }
                break;
            case SV_INTER_CHR_TRNSLCTN:
                if ((detectSet & (1 << i)) != 0)
                    TGM_DetectTranslocation(pDetectPars, pLibTable);
                break;
            default:
                break;
        }
    }

    fclose(pLibInput);
    TGM_LibInfoTableFree(pLibTable);
}

TGM_SpecialID* TGM_SpecialIDRead(FILE* libInput)
{
    uint32_t size = 0;
    unsigned int readSize = fread(&size, sizeof(uint32_t), 1, libInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of special references.\n");

    if (size == 0)
        return NULL;

    TGM_SpecialID* pSpecialID = (TGM_SpecialID*) malloc(sizeof(TGM_SpecialID));
    if (pSpecialID == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a special ID object.\n");

    pSpecialID->names = (char (*)[3]) calloc(sizeof(char), 3 * size);
    if (pSpecialID->names == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the special IDs.\n");

    for (unsigned int i = 0; i != size; ++i)
    {
        readSize = fread(pSpecialID->names[i], sizeof(char), 2, libInput);
        if (readSize != 2)
            TGM_ErrQuit("ERROR: Cannot read special ID from the library information file.\n");
    }

    pSpecialID->size = size;

    return pSpecialID;
}

void TGM_SpecialIDFree(TGM_SpecialID* pSpecialID)
{
    if (pSpecialID != NULL)
    {
        free(pSpecialID->names);
        free(pSpecialID);
    }
}

void TGM_LocalPairArrayRead(TGM_LocalPairArray* pLocalPairArray, FILE* input)
{
    int64_t size = 0;
    unsigned int readSize = fread(&(size), sizeof(int64_t), 1, input);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the size of local pair array.\n");

    if (size > pLocalPairArray->capacity)
        TGM_ARRAY_RESIZE_NO_COPY(pLocalPairArray, size, TGM_LocalPair);
    else if (size <= 0)
        return;

    pLocalPairArray->size = size;
    readSize = fread(pLocalPairArray->data, sizeof(TGM_LocalPair), size, input);
    if (readSize != size)
        TGM_ErrQuit("ERROR: Cannot read the local pair array.\n");
}

void TGM_CrossPairArrayRead(TGM_CrossPairArray* pCrossPairArray, FILE* input)
{
    int64_t size = 0;
    unsigned int readSize = fread(&(size), sizeof(int64_t), 1, input);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the size of cross pair array.\n");

    if (size > pCrossPairArray->capacity)
        TGM_ARRAY_RESIZE_NO_COPY(pCrossPairArray, size, TGM_CrossPair);
    else if (size <= 0)
        return;

    pCrossPairArray->size = size;
    readSize = fread(pCrossPairArray->data, sizeof(TGM_CrossPair), size, input);
    if (readSize != size)
        TGM_ErrQuit("ERROR: Cannot read the cross pair array.\n");
}

void TGM_SpecialPairArrayRead(TGM_SpecialPairArray* pSpecialPairArray, FILE* input)
{
    int64_t size = 0;
    unsigned int readSize = fread(&(size), sizeof(int64_t), 1, input);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the size of special pair array.\n");

    if (size > pSpecialPairArray->capacity)
        TGM_ARRAY_RESIZE_NO_COPY(pSpecialPairArray, size, TGM_SpecialPair);
    else if (size <= 0)
        return;

    pSpecialPairArray->size = size;
    readSize = fread(pSpecialPairArray->data, sizeof(TGM_SpecialPair), size, input);
    if (readSize != size)
        TGM_ErrQuit("ERROR: Cannot read the special pair array.\n");
}

/*  
void TGM_SpecialPairTableReadID(TGM_SpecialPairTable* pSpeicalPairTable, FILE* libInput)
{
    unsigned int readSize = fread(&(pSpeicalPairTable->size), sizeof(uint32_t), 1, libInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the size of the sp")

    if (pSpeicalPairTable->size > pSpeicalPairTable->capacity)
    {
        free(pSpeicalPairTable->names);
        pSpeicalPairTable->names = (char (*)[3]) malloc(sizeof(char) * 3 * pSpeicalPairTable->size);
        pSpeicalPairTable->capacity = pSpeicalPairTable->size;
    }
    else if (pSpeicalPairTable->size == 0)
        return;

    int ret = 0;
    khiter_t khIter = 0;
    for (unsigned int i = 0; i != pSpeicalPairTable->size; ++i)
    {
        fread(pSpeicalPairTable->names[i], sizeof(char), 2, libInput);
        pSpeicalPairTable->names[i][2] = '\0';

        khash_t(name)* pHash = pSpeicalPairTable->nameHash;
        khIter = kh_put(name, pHash, pSpeicalPairTable->names[i], &ret);
        kh_value(pHash, khIter) = i;
    }
}
*/

void SV_AssistArrayResize(SV_AssistArray* pAssistArray, unsigned int newSize)
{
    if (newSize > pAssistArray->capacity)
    {
        free(pAssistArray->pFragLenDiff);
        free(pAssistArray->pMapQ3);
        free(pAssistArray->pMapQ5);

        pAssistArray->capacity = newSize * 2;

        pAssistArray->pFragLenDiff = (int*) malloc(sizeof(int) * pAssistArray->capacity);
        pAssistArray->pMapQ3 = (int*) malloc(sizeof(int) * pAssistArray->capacity);
        pAssistArray->pMapQ5 = (int*) malloc(sizeof(int) * pAssistArray->capacity);
    }

    pAssistArray->size = newSize;
}

void TGM_DetectDel(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable)
{

}

void TGM_DetectTandemDup(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable)
{

}

void TGM_DetectInversion(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable)
{

}

void TGM_DetectSpecial(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable, const TGM_SpecialID* pSpecialID)
{
    char specialFileName[TGM_MAX_LINE];
    unsigned int workingDirLen = strlen(pDetectPars->workingDir);
    strcpy(specialFileName, pDetectPars->workingDir);
    strcpy(specialFileName + workingDirLen, TGM_READ_PAIR_FILE_NAME_TEMPLATE[4]);

    char* refPos = specialFileName + workingDirLen + 3;
    TGM_SpecialPairArray* pSpecialPairArray = NULL;
    TGM_ARRAY_ALLOC(pSpecialPairArray, 10, TGM_SpecialPairArray, TGM_SpecialPair);

    unsigned int numArray = pSpecialID->size * 2;
    TGM_ReadPairAttrbtArray** pAttrbtArrays = (TGM_ReadPairAttrbtArray**) malloc(sizeof(TGM_ReadPairAttrbtArray*) * numArray);
    if (pAttrbtArrays == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the read pair attribute arrays.\n");

    for (unsigned int i = 0; i != numArray; ++i)
        pAttrbtArrays[i] = TGM_ReadPairAttrbtArrayAlloc(pLibTable->size);

    TGM_Cluster* pCluster3 = TGM_ClusterAlloc(pDetectPars->minNumClustered);
    TGM_Cluster* pCluster5 = TGM_ClusterAlloc(pDetectPars->minNumClustered);

    TGM_SpecialEventArray* pSpecialEventArray = NULL;
    TGM_ARRAY_ALLOC(pSpecialEventArray, DEFAULT_SV_CAPACITY, TGM_SpecialEventArray, TGM_SpecialEvent);

    for (unsigned int k = pDetectPars->workingRefID[0]; k <= pDetectPars->workingRefID[1]; ++k)
    {
        sprintf(refPos, "%03u", k);
        // fix the null character introduced by sprintf
        *(refPos + 3) = '_';
        FILE* specialInput = fopen(specialFileName, "rb");
        if (specialInput == NULL)
            TGM_ErrQuit("ERROR: Cannot open special pair file.\n");

        TGM_SpecialPairArrayRead(pSpecialPairArray, specialInput);
        TGM_ReadPairMakeSpecial(pAttrbtArrays, numArray, pSpecialPairArray, pLibTable);

        for (unsigned int i = 0; i != numArray; i += 2)
        {
            TGM_ARRAY_RESET(pSpecialEventArray);
            TGM_ClusterInit(pCluster3, pAttrbtArrays[i]);
            TGM_ClusterInit(pCluster5, pAttrbtArrays[i + 1]);

            TGM_ClusterMake(pCluster3);
            TGM_ClusterMake(pCluster5);

            TGM_ClusterBuild(pCluster3);
            TGM_ClusterBuild(pCluster5);

            TGM_ClusterFinalize(pCluster3);
            TGM_ClusterFinalize(pCluster5);

            TGM_ClusterClean(pCluster3);
            TGM_ClusterClean(pCluster5);

            // TGM_ClusterPrint(pCluster3, k, pSpecialID->names[i / 2]);
            // TGM_ClusterPrint(pCluster5, k, pSpecialID->names[i / 2]);

            const TGM_Cluster* pCluster = pCluster3;
            unsigned int numClsElmnts = pCluster->pElmntArray->length;
            for (unsigned int j = 0; j != numClsElmnts; ++j)
            {
                if (pCluster->pElmntArray->data[j].numReadPair == 0)
                    continue;

                if (TGM_ARRAY_IS_FULL(pSpecialEventArray))
                    TGM_ARRAY_RESIZE(pSpecialEventArray, pSpecialEventArray->capacity * 2, TGM_SpecialEvent);

                TGM_SpecialEvent* pSpecialEvent = pSpecialEventArray->data + pSpecialEventArray->size;
                TGM_SpecialEventMake(pSpecialEvent, pCluster, j, pSpecialPairArray, pLibTable);
                ++(pSpecialEventArray->size);
            }

            pCluster = pCluster5;
            numClsElmnts = pCluster->pElmntArray->length;
            for (unsigned int j = 0; j != numClsElmnts; ++j)
            {
                if (pCluster->pElmntArray->data[j].numReadPair == 0)
                    continue;

                if (TGM_ARRAY_IS_FULL(pSpecialEventArray))
                    TGM_ARRAY_RESIZE(pSpecialEventArray, pSpecialEventArray->capacity * 2, TGM_SpecialEvent);

                TGM_SpecialEvent* pSpecialEvent = pSpecialEventArray->data + pSpecialEventArray->size;
                TGM_SpecialEventMake(pSpecialEvent, pCluster, j, pSpecialPairArray, pLibTable);
                ++(pSpecialEventArray->size);
            }

            qsort(pSpecialEventArray->data, pSpecialEventArray->size, sizeof(TGM_SpecialEvent), CompareSpecialEvents);

            unsigned int headIndex = 1;
            unsigned int tailIndex = 0;
            unsigned int newSize = pSpecialEventArray->size;
            while (headIndex < pSpecialEventArray->size)
            {
                TGM_SpecialEvent mergedEvent;
                if (abs(pSpecialEventArray->data[headIndex].pos - pSpecialEventArray->data[tailIndex].pos) < pLibTable->fragLenMax)
                {
                    TGM_SpecialEvent* pHeadEvent = pSpecialEventArray->data + headIndex;
                    TGM_SpecialEvent* pTailEvent = pSpecialEventArray->data + tailIndex;

                    TGM_SpecialEventMerge(&mergedEvent, pHeadEvent, pTailEvent, pCluster3, pCluster5);
                    memcpy(pTailEvent, &mergedEvent, sizeof(TGM_SpecialEvent));

                    --newSize;
                }
                else
                {
                    ++tailIndex;
                    if (headIndex != tailIndex)
                        memcpy(pSpecialEventArray->data + tailIndex, pSpecialEventArray->data + headIndex, sizeof(TGM_SpecialEvent));
                }

                ++headIndex;
            }

            pSpecialEventArray->size = newSize;
            qsort(pSpecialEventArray->data, pSpecialEventArray->size, sizeof(TGM_SpecialEvent), CompareSpecialEvents);

            for (unsigned int h = 0; h != pSpecialEventArray->size; ++h)
                TGM_SpecialEventPrint(pSpecialEventArray->data + h);
        }

        fclose(specialInput);
    }

    TGM_ClusterFree(pCluster3);
    TGM_ClusterFree(pCluster5);

    TGM_ARRAY_FREE(pSpecialPairArray, TRUE);
    TGM_ARRAY_FREE(pSpecialEventArray, TRUE);

    for (unsigned int i = 0; i != numArray; ++i)
        TGM_ReadPairAttrbtArrayFree(pAttrbtArrays[i]);

    free(pAttrbtArrays);
}

void TGM_SpecialEventMake(TGM_SpecialEvent* pSpecialEvent, const TGM_Cluster* pCluster, unsigned int index, 
                         const TGM_SpecialPairArray* pSpecialPairArray, const TGM_LibInfoTable* pLibTable)
{
    const TGM_ClusterElmnt* pClusterElmnt = pCluster->pElmntArray->data + index;
    unsigned int i = pClusterElmnt->startIndex;
    unsigned int numReadPair = pClusterElmnt->numReadPair;

    unsigned int posMin5 = UINT_MAX;
    unsigned int endMax5 = 0;

    unsigned int posMin3 = UINT_MAX;
    unsigned int endMax3 = 0;

    const TGM_SpecialPair* pSpecialPair = NULL;

    do
    {
        unsigned int origIndex = pCluster->pAttrbtArray->data[i].origIndex;
        pSpecialPair = pSpecialPairArray->data + origIndex;

        int fragLenMedian = pLibTable->pLibInfo[pSpecialPair->readGrpID].fragLenMedian;

        int posUniq = pSpecialPair->pos[0];
        int endUniq = pSpecialPair->end[0];

        if (pSpecialPair->readPairType == PT_SPECIAL3)
        {
            int posMultiple = posUniq + fragLenMedian - (pSpecialPair->end[1] - pSpecialPair->pos[1] + 1);
            int endMultiple = posUniq + fragLenMedian - 1;

            if (posUniq < posMin5)
                posMin5 = posUniq;

            if (endUniq > endMax5)
                endMax5 = endUniq;

            if (posMultiple < posMin3)
                posMin3 = posMultiple;

            if (endMultiple > endMax3)
                endMax3 = endMultiple;
        }
        else
        {
            int posMultiple = endUniq + 1 - fragLenMedian;
            int endMultiple = endUniq + pSpecialPair->end[1] + 1 - fragLenMedian - pSpecialPair->pos[1];

            if (posMultiple < posMin5)
                posMin5 = posMultiple;

            if (endMultiple > endMax5)
                endMax5 = endMultiple;

            if (posUniq < posMin3)
                posMin3 = posUniq;

            if (endUniq > endMax3)
                endMax3 = endUniq;
        }

        i = pCluster->pNext[i];

    }while(i != pClusterElmnt->startIndex);

    pSpecialEvent->refID = pSpecialPair->refID[0];
    if (pSpecialPair->readPairType == PT_SPECIAL3)
    {
        pSpecialEvent->pos = endMax5;
        pSpecialEvent->numFrag[0] = numReadPair;
        pSpecialEvent->numFrag[1] = 0;
        pSpecialEvent->clusterID[0] = pClusterElmnt->startIndex;
    }
    else
    {
        pSpecialEvent->pos = posMin3;
        pSpecialEvent->numFrag[0] = 0;
        pSpecialEvent->numFrag[1] = numReadPair;
        pSpecialEvent->clusterID[1] = pClusterElmnt->startIndex;
    }

    pSpecialEvent->length = 0;

    pSpecialEvent->pos5[0] = posMin5;
    pSpecialEvent->pos5[1] = endMax5;

    pSpecialEvent->pos3[0] = posMin3;
    pSpecialEvent->pos3[1] = endMax3;

    pSpecialEvent->posUncertainty = DoubleRoundToInt((double) ((endMax5 - posMin5) + (endMax3 - posMin3)) / (double) (2 * numReadPair));
}

void TGM_SpecialEventMerge(TGM_SpecialEvent* pMergedEvent, const TGM_SpecialEvent* pHeadEvent, const TGM_SpecialEvent* pTailEvent, TGM_Cluster* pCluster3, TGM_Cluster* pCluster5)
{
    if (pTailEvent->numFrag[0] == 0)
        TGM_SWAP(pHeadEvent, pTailEvent, const TGM_SpecialEvent*);

    memcpy(pMergedEvent, pTailEvent, sizeof(TGM_SpecialEvent));

    TGM_Bool bothForward = (pHeadEvent->numFrag[0] > 0) && (pTailEvent->numFrag[0] > 0);
    TGM_Bool bothReversed = (pHeadEvent->numFrag[1] > 0) && (pTailEvent->numFrag[1] > 0);

    if (bothForward)
    {
        pMergedEvent->pos5[0] = (pTailEvent->pos5[0] < pHeadEvent->pos5[0] ? pTailEvent->pos5[0] : pHeadEvent->pos5[0]);
        pMergedEvent->pos5[1] = (pTailEvent->pos5[1] > pHeadEvent->pos5[1] ? pTailEvent->pos5[1] : pHeadEvent->pos5[1]);
        TGM_ClusterConnect(pCluster3, pHeadEvent->clusterID[0], pTailEvent->clusterID[0]);
        pMergedEvent->clusterID[0] = pHeadEvent->clusterID[0];
    }
    else if (bothReversed)
    {
        pMergedEvent->pos3[0] = (pTailEvent->pos3[0] < pHeadEvent->pos3[0] ? pTailEvent->pos3[0] : pHeadEvent->pos3[0]);
        pMergedEvent->pos3[1] = (pTailEvent->pos3[1] > pHeadEvent->pos3[1] ? pTailEvent->pos3[1] : pHeadEvent->pos3[1]);
        TGM_ClusterConnect(pCluster5, pHeadEvent->clusterID[1], pTailEvent->clusterID[1]);
        pMergedEvent->clusterID[1] = pHeadEvent->clusterID[1];
    }
    else
    {
        pMergedEvent->pos3[0] = pHeadEvent->pos3[0];
        pMergedEvent->pos3[1] = pHeadEvent->pos3[1];
        pMergedEvent->clusterID[1] = pHeadEvent->clusterID[1];
    }

    int numFrag5 = pMergedEvent->numFrag[0] = pHeadEvent->numFrag[0] + pTailEvent->numFrag[0];
    pMergedEvent->numFrag[1] = pHeadEvent->numFrag[1] + pTailEvent->numFrag[1];

    int32_t posTail = (pTailEvent->numFrag[0] > 0 ? pTailEvent->pos5[1] : pTailEvent->pos3[0]);
    int32_t posHead = (pHeadEvent->numFrag[1] > 0 ? pHeadEvent->pos3[0] : pHeadEvent->pos5[1]);

    pMergedEvent->pos = (posHead < posTail ? posHead : posTail);
    pMergedEvent->posUncertainty = DoubleRoundToInt((double) (pMergedEvent->pos5[1] - pMergedEvent->pos5[0]) / numFrag5);
}

void TGM_SpecialEventPrint(const TGM_SpecialEvent* pSpecialEvent)
{
    printf("chr%d\t%d\t%d\t%d\t%d\n", pSpecialEvent->refID + 1, pSpecialEvent->pos, pSpecialEvent->pos + 1, pSpecialEvent->numFrag[0], pSpecialEvent->numFrag[1]);
}

void TGM_DetectTranslocation(const TGM_ReadPairDetectPars* pDetectPars, const TGM_LibInfoTable* pLibTable)
{

}

void TGM_ReadPairFindDel(TGM_DelArray* pDelArray, SV_AssistArray* pAssistArray, const TGM_LocalPairArray* pLongPairArray,
                        TGM_Cluster* pDelCluster, const TGM_LibInfoTable* pLibTable, const TGM_ReadPairDetectPars* pPars)
{
    TGM_ARRAY_RESIZE_NO_COPY(pDelArray, pDelCluster->pElmntArray->size, TGM_DelEvent);

    unsigned int i = 0;
    unsigned int counter = 0;
    unsigned int numClusters = pDelCluster->pElmntArray->size;

    while (counter != numClusters)
    {
        if (pDelCluster->pElmntArray->data[i].numReadPair == 0)
        {
            ++i;
            continue;
        }

        int posMin5 = INT_MAX;
        int posMax5 = 0;

        int posMin3 = INT_MAX;
        int posMax3 = 0;

        int endMax5 = 0;
        int endMax3 = 0;

        int fragLenDiffMin = INT_MAX;
        int fragLenDiffMax = 0;

        int fragLenMax = 0;

        int numReadPair = pDelCluster->pElmntArray->data[i].numReadPair;
        SV_AssistArrayResize(pAssistArray, numReadPair);

        unsigned int startIndex = pDelCluster->pElmntArray->data[i].startIndex;
        unsigned int j = startIndex;
        unsigned int assistIndex = 0;
        const TGM_LocalPair* pLongPair = NULL;

        do
        {
            unsigned int origIndex = pDelCluster->pAttrbtArray->data[j].origIndex;
            pLongPair = (pLongPairArray->data + origIndex);

            if (pLongPair->upPos < posMin5)
                posMin5 = pLongPair->upPos;

            if (pLongPair->upPos > posMax5)
                posMax5 = pLongPair->upPos;

            if (pLongPair->upEnd > endMax5)
                endMax5 = pLongPair->upEnd;

            if (pLongPair->downPos < posMin3)
                posMin3 = pLongPair->downPos;

            if (pLongPair->downPos > posMax3)
                posMax3 = pLongPair->downPos;

            unsigned int downEnd = pLongPair->upPos + pLongPair->fragLen;
            if (downEnd > endMax3)
                endMax3 = downEnd;

            unsigned int upperFragLen = pLibTable->pLibInfo[pLongPair->readGrpID].fragLenHigh;
            unsigned int medianFragLen = pLibTable->pLibInfo[pLongPair->readGrpID].fragLenMedian;

            if (upperFragLen > fragLenMax)
                fragLenMax = upperFragLen;

            pAssistArray->pFragLenDiff[assistIndex] = pLongPair->fragLen - medianFragLen;
            pAssistArray->pMapQ5[assistIndex] = pLongPair->upMapQ;
            pAssistArray->pMapQ3[assistIndex] = pLongPair->downMapQ;

            if (pAssistArray->pFragLenDiff[assistIndex] < fragLenDiffMin)
                fragLenDiffMin = pAssistArray->pFragLenDiff[assistIndex];

            if (pAssistArray->pFragLenDiff[assistIndex] > fragLenDiffMax)
                fragLenDiffMax = pAssistArray->pFragLenDiff[assistIndex];

            j = pDelCluster->pNext[j];
            ++assistIndex;

        }while(j != startIndex);

        int eventLength = FindMedianInt(pAssistArray->pFragLenDiff, pAssistArray->size);
        if (eventLength < 1)
            continue;

        int k = pDelArray->size;

        pDelArray->data[k].clusterID = i;

        pDelArray->data[k].refID = pLongPair->refID;
        pDelArray->data[k].pos = endMax5 + 1;
        pDelArray->data[k].end = posMin3;
        pDelArray->data[k].length = FindMedianInt(pAssistArray->pFragLenDiff, pAssistArray->size);
        pDelArray->data[k].quality = (int) ((numReadPair * 100.0) / (numReadPair + 10.0));

        pDelArray->data[k].pos5[0] = posMin5;
        pDelArray->data[k].pos5[1] = endMax5;
        pDelArray->data[k].pos5[2] = posMax5;

        pDelArray->data[k].pos3[0] = posMin3;
        pDelArray->data[k].pos3[1] = endMax3;
        pDelArray->data[k].pos3[2] = posMax3;

        pDelArray->data[k].CIpos[0] = -(posMax5 - posMin5) / numReadPair;
        pDelArray->data[k].CIpos[1] = (posMax3 - posMin3) / numReadPair;

        pDelArray->data[k].CIend[0] = -(posMax5 - posMin5) / numReadPair;
        pDelArray->data[k].CIend[1] = (posMax3 - posMin3) / numReadPair;

        pDelArray->data[k].CIlen[0] = -(eventLength - fragLenDiffMin) / numReadPair;
        pDelArray->data[k].CIlen[1] = (fragLenDiffMax - eventLength) / numReadPair;

        pDelArray->data[k].mapQ5 = FindMedianInt(pAssistArray->pMapQ5, pAssistArray->size);
        pDelArray->data[k].mapQ3 = FindMedianInt(pAssistArray->pMapQ3, pAssistArray->size);

        TGM_Bool isOverlap = FALSE;
        if (pDelArray->size > 0)
        {
            int lastIndex = pDelArray->size - 1;

            if (pDelArray->data[lastIndex].end > pDelArray->data[k].pos 
                && abs(pDelArray->data[lastIndex].pos - pDelArray->data[k].pos) < fragLenMax)
            {
                isOverlap = TRUE;
            }
        }

        // correction to pos and len for deletions bracketed by repeat region 
        int fragLenDiffError = eventLength - DoubleRoundToInt(pDelCluster->pElmntArray->data[i].mean[1]);
        if (fragLenDiffError > 4 * (pDelArray->data[k].CIpos[1]))
        {
            pDelArray->data[k].pos += (int) (fragLenDiffError / 2);
            pDelArray->data[k].CIpos[0] = (int) (-fragLenDiffError / 2);
            pDelArray->data[k].CIpos[1] = (int) (fragLenDiffError / 2);

            pDelArray->data[k].length = DoubleRoundToInt(pDelCluster->pElmntArray->data[i].mean[1]);
            pDelArray->data[k].CIlen[0] = DoubleRoundToInt(-pDelCluster->pElmntArray->data[i].std[1] / sqrt(numReadPair));
            pDelArray->data[k].CIlen[1] = DoubleRoundToInt(pDelCluster->pElmntArray->data[i].std[1] / sqrt(numReadPair));
        }

        if (isOverlap)
            TGM_DelEventMerge(pDelArray, pDelCluster);
        else
        {
            if (numReadPair >= pPars->minNumClustered && pDelArray->data[k].length >= pPars->minEventLength)
                ++(pDelArray->size);
        }

        ++i;
        ++counter;
    }

    qsort(pDelArray->data, pDelArray->size, sizeof(TGM_DelEvent), CompareDelEvents);

    //TGM_DelEventGenotype(pDelArray, pLibTable);
}
