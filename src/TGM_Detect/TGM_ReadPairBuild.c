/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairBuild.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/08/2011 11:45:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <math.h>
#include <unistd.h>
#include <assert.h>

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Utilities.h"
#include "TGM_BamPairAux.h"
#include "TGM_BamInStream.h"
#include "TGM_ReadPairBuild.h"

#define DEFAULT_RP_INFO_CAPACITY 50

#define TGM_NUM_RP_FILE 6

#define DEFAULT_RP_TABLE_CAPACITY 10

#define MAX_MQ_FOR_MULTIPLE 10

#define MIN_MQ_FOR_UNIQUE 20

static const char* TGM_LibTableFileName = "lib_table.dat";

static const char* TGM_HistFileName = "hist.dat";

static const char* TGM_READ_PAIR_FILE_NAME_TEMPLATE[] = 
{
    "refXXX_long_pairs.dat",

    "refXXX_short_pairs.dat",

    "refXXX_reversed_pairs.dat",

    "refXXX_inverted_pairs.dat",

    "refXXX_special_pairs.dat",

    "refXXX_cross_pairs.dat",
};

enum
{
    TGM_LONG_PAIR_FILE = 0,

    TGM_SHORT_PAIR_FILE = 1,

    TGM_REVERSED_PAIR_FILE = 2,

    TGM_INVERTED_PAIR_FILE = 3,

    TGM_SPEICAL_PAIR_FILE = 4,

    TGM_CROSS_PAIR_FILE = 5
};

typedef struct TGM_ReadPairOutStream
{
    FILE* output;

    int64_t numPairs;

    int64_t splitInfoPos;

}TGM_ReadPairOutStream;

KHASH_MAP_INIT_STR(name, uint32_t);

KHASH_MAP_INIT_INT(file, TGM_ReadPairOutStream);

// check the read pair type
static SV_ReadPairType TGM_CheckReadPairType(const TGM_ZAtag* pZAtag, const TGM_PairStats* pPairStats, const bam1_t* pDownAlgn, 
                                            const TGM_MateInfo* pMateInfo, uint8_t minMQ, uint8_t minSpMQ, const TGM_LibInfoTable* pLibTable)
{
    // check the za tag first if it is available
    if (pZAtag != NULL)
    {
        // special reference name starting with ' ' means empty
        if (pZAtag->spRef[0][0] != ' ' && pZAtag->spRef[1][0] != ' ')
            return PT_UNKNOWN;
        else if (pZAtag->spRef[0][0] == ' ' && pZAtag->spRef[1][0] != ' ')
        {
            if (pZAtag->bestMQ[0] >= minSpMQ)
            {
                switch(pPairStats->pairMode)
                {
                    case TGM_1F2F:
                    case TGM_1F2R:
                    case TGM_2F1F:
                    case TGM_2F1R:
                        return PT_SPECIAL3;
                        break;
                    case TGM_1R2F:
                    case TGM_1R2R:
                    case TGM_2R1F:
                    case TGM_2R1R:
                        return PT_SPECIAL5;
                        break;
                    default:
                        return PT_UNKNOWN;
                        break;
                }
            }
            else
                return PT_UNKNOWN;
        }
        else if (pZAtag->spRef[0][0] != ' ' && pZAtag->spRef[1][0] == ' ')
        {
            if (pZAtag->bestMQ[1] >= minSpMQ)
            {
                switch(pPairStats->pairMode)
                {
                    case TGM_2F1F:
                    case TGM_2R1F:
                    case TGM_1F2F:
                    case TGM_1R2F:
                        return PT_SPECIAL3;
                        break;
                    case TGM_2F1R:
                    case TGM_2R1R:
                    case TGM_1F2R:
                    case TGM_1R2R:
                        return PT_SPECIAL5;
                        break;
                    default:
                        return PT_UNKNOWN;
                        break;
                }
            }
            else
                return PT_UNKNOWN;
        }

        /*
        // only leave those unique-unique pairs for further judgement
        if (pZAtag->numMappings[0] != 1 || pZAtag->numMappings[1] != 1)
            return PT_UNKNOWN;
        */
    }

    // TODO: evaluate the possibilty to detect MEI without ZA
    /*
    if (pMateInfo != NULL)
    {
        if (pMateInfo->mapQ > MIN_MQ_FOR_UNIQUE && pDownAlgn->core.qual < MAX_MQ_FOR_MULTIPLE)
            return PT_SPECIAL3;
        else if (pMateInfo->mapQ < MAX_MQ_FOR_MULTIPLE && pDownAlgn->core.qual > MIN_MQ_FOR_UNIQUE)
            return PT_SPECIAL5;
    }
    */

    // two mates of a read pair aligned to different references
    if (pPairStats->fragLen == -1)
        return PT_CROSS;

    double fragLenHigh = pLibTable->pLibInfo[pPairStats->readGrpID].fragLenHigh;
    double fragLenLow = pLibTable->pLibInfo[pPairStats->readGrpID].fragLenLow;

    SV_ReadPairType type1 = SV_ReadPairTypeMap[0][pPairStats->pairMode];
    SV_ReadPairType type2 = SV_ReadPairTypeMap[1][pPairStats->pairMode];

    if (type1 == PT_NORMAL || type2 == PT_NORMAL)
    {
        if (pPairStats->fragLen > fragLenHigh)
            return PT_LONG;
        else if (pPairStats->fragLen < fragLenLow)
            return PT_SHORT;
        else
            return PT_NORMAL;
    }
    else if ((type1 == PT_UNKNOWN && type2 == PT_UNKNOWN) || (type1 != PT_UNKNOWN && type2 != PT_UNKNOWN))
    {
        return PT_UNKNOWN;
    }
    else
    {
        return (type1 == PT_UNKNOWN ? type2 : type1);
    }
}

// update the local pair array
static void TGM_LocalPairArrayUpdate(TGM_LocalPairArray* pLocalPairArray, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const TGM_ZAtag* pZAtag, const TGM_PairStats* pPairStats,
                                    const TGM_MateInfo* pMateInfo, const TGM_LibInfoTable* pLibTable, const TGM_FragLenHistArray* pHistArray,  SV_ReadPairType readPairType)
{
    if (pLocalPairArray->size == pLocalPairArray->capacity)
        TGM_ARRAY_RESIZE(pLocalPairArray, pLocalPairArray->capacity * 2, TGM_LocalPair);

    TGM_LocalPair* pLocalPair = pLocalPairArray->data + pLocalPairArray->size;

    if (pMateInfo != NULL)
    {
        pLocalPair->refID = pDownAlgn->core.tid;
        pLocalPair->upPos = pDownAlgn->core.mpos;
        pLocalPair->downPos = pDownAlgn->core.pos;
        pLocalPair->upEnd = pMateInfo->upEnd;

        pLocalPair->upNumMM = pMateInfo->numMM;
        pLocalPair->downNumMM = TGM_GetNumMismatchFromBam(pDownAlgn);

        pLocalPair->upMapQ = pMateInfo->mapQ;
        pLocalPair->downMapQ = pDownAlgn->core.qual;
    }
    else
    {
        pLocalPair->refID = pUpAlgn->core.tid;
        pLocalPair->upPos = pUpAlgn->core.pos;
        pLocalPair->downPos = pUpAlgn->core.mpos;
        pLocalPair->upEnd = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));

        if (pZAtag != NULL)
        {
            pLocalPair->upNumMM = pZAtag->numMM[0];
            pLocalPair->downNumMM = pZAtag->numMM[1];

            pLocalPair->upMapQ = pZAtag->bestMQ[0];
            pLocalPair->downMapQ = pZAtag->bestMQ[1];
        }
        else
        {
            pLocalPair->upNumMM = TGM_GetNumMismatchFromBam(pUpAlgn);
            pLocalPair->downNumMM = TGM_GetNumMismatchFromBam(pDownAlgn);

            pLocalPair->upMapQ = pUpAlgn->core.qual;
            pLocalPair->downMapQ = pDownAlgn->core.qual;
        }
    }

    pLocalPair->fragLen = pPairStats->fragLen;
    pLocalPair->readGrpID = pPairStats->readGrpID;
    pLocalPair->pairMode = pPairStats->pairMode;
    pLocalPair->readPairType = readPairType;

    pLocalPair->fragLenQual = TGM_FragLenHistArrayGetFragLenQual(pHistArray, pLibTable->size - pPairStats->readGrpID, pPairStats->fragLen);

    ++(pLocalPairArray->chrCount[pLocalPair->refID]);
    ++(pLocalPairArray->size);
}

// update the inverted array
static void TGM_InvertedPairArrayUpdate(TGM_LocalPairArray* pInvertedPairArray, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const TGM_ZAtag* pZAtag, 
                                       const TGM_PairStats* pPairStats, const TGM_MateInfo* pMateInfo, SV_ReadPairType readPairType)
{
    if (pInvertedPairArray->size == pInvertedPairArray->capacity)
        TGM_ARRAY_RESIZE(pInvertedPairArray, pInvertedPairArray->capacity * 2, TGM_LocalPair);

    TGM_LocalPair* pInvertedPair = pInvertedPairArray->data + pInvertedPairArray->size;

    if (pMateInfo != NULL)
    {
        pInvertedPair->refID = pDownAlgn->core.tid;
        pInvertedPair->upPos = pDownAlgn->core.mpos;
        pInvertedPair->downPos = pDownAlgn->core.pos;
        pInvertedPair->upEnd = pMateInfo->upEnd;

        pInvertedPair->upNumMM = pMateInfo->numMM;
        pInvertedPair->downNumMM = TGM_GetNumMismatchFromBam(pDownAlgn);

        pInvertedPair->upMapQ = pMateInfo->mapQ;
        pInvertedPair->downMapQ = pDownAlgn->core.qual;
    }
    else
    {
        pInvertedPair->refID = pUpAlgn->core.tid;
        pInvertedPair->upPos = pUpAlgn->core.pos;
        pInvertedPair->downPos = pUpAlgn->core.mpos;
        pInvertedPair->upEnd = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));

        if (pZAtag != NULL)
        {
            pInvertedPair->upNumMM = pZAtag->numMM[0];
            pInvertedPair->downNumMM = pZAtag->numMM[1];

            pInvertedPair->upMapQ = pZAtag->bestMQ[0];
            pInvertedPair->downMapQ = pZAtag->bestMQ[1];
        }
        else
        {
            pInvertedPair->upNumMM = TGM_GetNumMismatchFromBam(pUpAlgn);
            pInvertedPair->downNumMM = TGM_GetNumMismatchFromBam(pDownAlgn);

            pInvertedPair->upMapQ = pUpAlgn->core.qual;
            pInvertedPair->downMapQ = pDownAlgn->core.qual;
        }
    }

    pInvertedPair->fragLen = pPairStats->fragLen;
    pInvertedPair->readGrpID = pPairStats->readGrpID;
    pInvertedPair->pairMode = pPairStats->pairMode;
    pInvertedPair->readPairType = readPairType;

    pInvertedPair->fragLenQual = INVALID_FRAG_LEN_QUAL;

    ++(pInvertedPairArray->chrCount[pInvertedPair->refID]);
    ++(pInvertedPairArray->size);
}

// update the cross pair array
static void TGM_CrossPairArrayUpdate(TGM_CrossPairArray* pCrossPairArray, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const TGM_ZAtag* pZAtag, const TGM_PairStats*  pPairStats)
{
    if (pCrossPairArray->size == pCrossPairArray->capacity)
        TGM_ARRAY_RESIZE(pCrossPairArray, pCrossPairArray->capacity * 2, TGM_CrossPair);

    TGM_CrossPair* pCrossPair = pCrossPairArray->data + pCrossPairArray->size;

    pCrossPair->readGrpID = pPairStats->readGrpID;

    assert(pUpAlgn->core.tid >= 0);
    assert(pUpAlgn->core.mtid >= 0);

    pCrossPair->upRefID = pUpAlgn->core.tid;
    pCrossPair->downRefID = pUpAlgn->core.mtid;

    pCrossPair->upPos = pUpAlgn->core.pos;
    pCrossPair->downPos = pUpAlgn->core.mpos;

    if (pZAtag != NULL)
    {
        pCrossPair->upEnd = pZAtag->end[0];
        pCrossPair->downEnd = pZAtag->end[1];

        pCrossPair->upNumMM = pZAtag->numMM[0];
        pCrossPair->downNumMM = pZAtag->numMM[1];

        pCrossPair->upMapQ = pZAtag->bestMQ[0];
        pCrossPair->downMapQ = pZAtag->bestMQ[1];
    }
    else
    {
        pCrossPair->upEnd = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));
        pCrossPair->downEnd = bam_calend(&(pDownAlgn->core), bam1_cigar(pDownAlgn));

        pCrossPair->upNumMM = TGM_GetNumMismatchFromBam(pUpAlgn);
        pCrossPair->downNumMM = TGM_GetNumMismatchFromBam(pDownAlgn);

        pCrossPair->upMapQ = pUpAlgn->core.qual;
        pCrossPair->downMapQ = pDownAlgn->core.qual;
    }

    pCrossPair->pairMode = pPairStats->pairMode;
    pCrossPair->readPairType = PT_CROSS;

    pCrossPair->fragLenQual = INVALID_FRAG_LEN_QUAL;

    ++(pCrossPairArray->chrCount[pCrossPair->upRefID]);
    ++(pCrossPairArray->size);
}

// update the special pair table
static void TGM_SpecialPairTableUpdate(TGM_SpecialPairTable* pSpecialPairTable, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const TGM_LibInfoTable* pLibTable, 
        const TGM_FragLenHistArray* pHistArray, const TGM_PairStats* pPairStats, const TGM_ZAtag* pZAtag, SV_ReadPairType readPairType)
{
    if (pSpecialPairTable->size == pSpecialPairTable->capacity)
    {
        pSpecialPairTable->capacity *= 2;
        pSpecialPairTable->names = (char (*)[3]) realloc(pSpecialPairTable->names, sizeof(char) * pSpecialPairTable->capacity);
        if (pSpecialPairTable->names == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the special reference names.\n");

        for (unsigned int i = pSpecialPairTable->size; i != pSpecialPairTable->capacity; ++i)
            pSpecialPairTable->names[i][2] = '\0';
    }

    // check which mate is the anchor
    int anchorIndex = 0;
    int specialIndex = 0;

    if (pZAtag->spRef[0][0] == ' ')
    {
        anchorIndex = 0;
        specialIndex = 1;
    }
    else
    {
        anchorIndex = 1;
        specialIndex = 0;
    }

    // copy the special reference name
    pSpecialPairTable->names[pSpecialPairTable->size][0] = pZAtag->spRef[specialIndex][0];
    pSpecialPairTable->names[pSpecialPairTable->size][1] = pZAtag->spRef[specialIndex][1];

    //  get the special reference ID
    int ret = 0;
    uint16_t spRefID = 0;
    khiter_t khIter = kh_put(name, pSpecialPairTable->nameHash, pSpecialPairTable->names[pSpecialPairTable->size], &ret);

    if (ret != 0)
    {
        kh_value((khash_t(name)*) pSpecialPairTable->nameHash, khIter) = pSpecialPairTable->size;
        spRefID = pSpecialPairTable->size;
        ++(pSpecialPairTable->size);
    }
    else
        spRefID = kh_value((khash_t(name)*) pSpecialPairTable->nameHash, khIter);

    const int32_t pRefID[2] = {pUpAlgn->core.tid, pUpAlgn->core.mtid};
    const int32_t pPos[2] = {pUpAlgn->core.pos, pUpAlgn->core.mpos};
    TGM_SpecialPairArray* pSpecialPairArray = NULL;

    // to see which special pairs array should we use
    if (pRefID[anchorIndex] == pUpAlgn->core.tid)
        pSpecialPairArray = &(pSpecialPairTable->array);
    else
        pSpecialPairArray = &(pSpecialPairTable->crossArray);

    if (pSpecialPairArray->capacity == pSpecialPairArray->size)
        TGM_ARRAY_RESIZE(pSpecialPairArray, pSpecialPairArray->capacity * 2, TGM_SpecialPair);

    TGM_SpecialPair* pSpecialPair = pSpecialPairArray->data + pSpecialPairArray->size;

    if (pUpAlgn->core.tid == pUpAlgn->core.mtid)
    {
        int32_t fragLenHigh = pLibTable->pLibInfo[pPairStats->readGrpID].fragLenHigh;
        int32_t fragLenLow = pLibTable->pLibInfo[pPairStats->readGrpID].fragLenLow;
        int32_t fragLenMedian = pLibTable->pLibInfo[pPairStats->readGrpID].fragLenMedian;

        // widen the fragment length window to clean up the clusters
        int32_t spFragLenHigh = fragLenHigh + (fragLenHigh - fragLenMedian);
        int32_t spFragLenLow = fragLenLow - (fragLenMedian - fragLenLow);

        int32_t mapDist = 0;
        if (readPairType == PT_SPECIAL3)                                        // anchor is on the forward strand
            mapDist = pZAtag->end[specialIndex] - pPos[anchorIndex];
        else                                                                    // anchor is on the reversed strand
            mapDist = pZAtag->end[anchorIndex] - pPos[specialIndex];

        // we have to filter out those special pairs with normal fragment length
        if (mapDist >= spFragLenLow && mapDist <= spFragLenHigh)
            return;
    }

    pSpecialPair->readGrpID = pPairStats->readGrpID;
    pSpecialPair->refID[0] = pRefID[anchorIndex];
    pSpecialPair->refID[1] = pRefID[specialIndex];

    pSpecialPair->pos[0] = pPos[anchorIndex];
    pSpecialPair->pos[1] = pPos[specialIndex];

    pSpecialPair->end[0] = pZAtag->end[anchorIndex];
    pSpecialPair->end[1] = pZAtag->end[specialIndex];

    pSpecialPair->numMM[0] = pZAtag->numMM[anchorIndex];
    pSpecialPair->numMM[1] = pZAtag->numMM[specialIndex];

    pSpecialPair->bestMQ[0] = pZAtag->bestMQ[anchorIndex];
    pSpecialPair->bestMQ[1] = pZAtag->bestMQ[specialIndex];

    pSpecialPair->numSpeicalHits = pZAtag->numMappings[specialIndex];
    pSpecialPair->pairMode = pPairStats->pairMode;
    pSpecialPair->readPairType = readPairType;
    pSpecialPair->specialID = spRefID;

    if (pRefID[anchorIndex] == pUpAlgn->core.tid)
        ++(pSpecialPairArray->chrCount[pSpecialPair->refID[0]]);

    ++(pSpecialPairArray->size);
}

static void TGM_SpecialPairTableClear(TGM_SpecialPairTable* pSpecialPairTable, unsigned int numChr)
{
    pSpecialPairTable->array.size = 0;
    memset(pSpecialPairTable->array.chrCount, 0, sizeof(uint64_t) * numChr);

    pSpecialPairTable->crossArray.size = 0;
}

static char* TGM_SpecialPairTableGetID(const TGM_SpecialPairTable* pSpecialPairTable)
{
    if (pSpecialPairTable == NULL || pSpecialPairTable->size == 0)
        return NULL;

    char* pSpecialID = (char*) malloc((pSpecialPairTable->size * 2 + 1) * sizeof(char));
    if (pSpecialID == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the special reference ID.\n");

    for (unsigned int i = 0; i != pSpecialPairTable->size; ++i)
        strncpy(pSpecialID + i * 2, pSpecialPairTable->names[i], 2);

    pSpecialID[pSpecialPairTable->size * 2] = '\0';

    return pSpecialID;
}

static void TGM_SplitPairArrayUpdate(TGM_SplitPairArray* pSplitPairArray, const bam1_t* pFirstPartial, const bam1_t* pSecondPartial, int32_t readGrpID, const TGM_ZAtag* pZAtag)
{
    if (pSplitPairArray->size == pSplitPairArray->capacity)
        TGM_ARRAY_RESIZE(pSplitPairArray, pSplitPairArray->capacity * 2, TGM_SplitPair);

    TGM_SplitPair* pSplitPair = pSplitPairArray->data + pSplitPairArray->size;

    pSplitPair->readGrpID = readGrpID;

    pSplitPair->refID[0] = pFirstPartial->core.tid;
    pSplitPair->refID[1] = pSecondPartial->core.tid;

    pSplitPair->pos[0] = pFirstPartial->core.pos;
    pSplitPair->pos[1] = pSecondPartial->core.pos;

    pSplitPair->end[0] = bam_calend(&(pFirstPartial->core), bam1_cigar(pFirstPartial));
    pSplitPair->end[1] = bam_calend(&(pSecondPartial->core), bam1_cigar(pSecondPartial));

    pSplitPair->numMM[0] = TGM_GetNumMismatchFromBam(pFirstPartial);
    pSplitPair->numMM[1] = TGM_GetNumMismatchFromBam(pSecondPartial);

    pSplitPair->bestMQ[0] = pFirstPartial->core.qual;
    pSplitPair->bestMQ[1] = pSecondPartial->core.qual;

    // FIXME: ZA tag has not been added in the split bam
    // change this line back when ZA is available
    // pSplitPair->readPairType = pZAtag->readPairType;
    pSplitPair->readPairType = 0;

    ++(pSplitPairArray->chrCount[pSplitPair->refID[0]]);
    ++(pSplitPairArray->size);
}

//===============================
// Constructors and Destructors
//===============================

TGM_SpecialPairTable* TGM_SpecialPairTableAlloc(unsigned int numChr)
{
    TGM_SpecialPairTable* pSpecialPairTable = (TGM_SpecialPairTable*) malloc(sizeof(TGM_SpecialPairTable));
    if (pSpecialPairTable == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a special pair table object.\n");

    pSpecialPairTable->size = 0;
    pSpecialPairTable->capacity = DEFAULT_RP_TABLE_CAPACITY;

    pSpecialPairTable->names = (char (*)[3]) calloc(sizeof(char), DEFAULT_RP_TABLE_CAPACITY * 3);
    if (pSpecialPairTable->names == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the special reference namese.\n");

    pSpecialPairTable->nameHash = kh_init(name);

    TGM_SpecialPairArray* pArray = &(pSpecialPairTable->array);
    TGM_ARRAY_INIT(pArray, DEFAULT_RP_INFO_CAPACITY, TGM_SpecialPair);

    TGM_SpecialPairArray* pCrossArray = &(pSpecialPairTable->crossArray);
    TGM_ARRAY_INIT(pCrossArray, DEFAULT_RP_INFO_CAPACITY, TGM_SpecialPair);

    pSpecialPairTable->array.chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
    if (pSpecialPairTable->array.chrCount == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

    return pSpecialPairTable;
}

void TGM_SpecialPairTableFree(TGM_SpecialPairTable* pSpecialPairTable)
{
    if (pSpecialPairTable != NULL)
    {
        free(pSpecialPairTable->array.chrCount);
        TGM_ARRAY_FREE(&(pSpecialPairTable->array), FALSE);

        TGM_ARRAY_FREE(&(pSpecialPairTable->crossArray), FALSE);

        free(pSpecialPairTable->names);
        kh_destroy(name, pSpecialPairTable->nameHash);

        free(pSpecialPairTable);
    }
}

TGM_ReadPairTable* TGM_ReadPairTableAlloc(uint32_t numChr, uint32_t detectSet)
{
    TGM_ReadPairTable* pNewInfoTable = (TGM_ReadPairTable*) calloc(sizeof(TGM_ReadPairTable), 1);
    if (pNewInfoTable == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for read pair information table.\n");

    pNewInfoTable->detectSet = detectSet;
    pNewInfoTable->numChr = numChr;
    pNewInfoTable->numPairs = 0;

    if ((detectSet & (1 << SV_DELETION)) != 0)
    {
        TGM_ARRAY_ALLOC(pNewInfoTable->pLongPairArray, DEFAULT_RP_INFO_CAPACITY, TGM_LocalPairArray, TGM_LocalPair);
        pNewInfoTable->pLongPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pLongPairArray->chrCount == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_TANDEM_DUP)) != 0)
    {
        TGM_ARRAY_ALLOC(pNewInfoTable->pShortPairArray, DEFAULT_RP_INFO_CAPACITY, TGM_LocalPairArray, TGM_LocalPair);
        pNewInfoTable->pShortPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pShortPairArray->chrCount == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

        TGM_ARRAY_ALLOC(pNewInfoTable->pReversedPairArray, DEFAULT_RP_INFO_CAPACITY, TGM_LocalPairArray, TGM_LocalPair);
        pNewInfoTable->pReversedPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pReversedPairArray->chrCount == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_INVERSION)) != 0)
    {
        TGM_ARRAY_ALLOC(pNewInfoTable->pInvertedPairArray, DEFAULT_RP_INFO_CAPACITY, TGM_LocalPairArray, TGM_LocalPair);
        pNewInfoTable->pInvertedPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pInvertedPairArray->chrCount == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_INTER_CHR_TRNSLCTN)) != 0)
    {
        TGM_ARRAY_ALLOC(pNewInfoTable->pCrossPairArray, DEFAULT_RP_INFO_CAPACITY, TGM_CrossPairArray, TGM_CrossPair);
        pNewInfoTable->pCrossPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pCrossPairArray->chrCount == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_SPECIAL)) != 0)
        pNewInfoTable->pSpecialPairTable = TGM_SpecialPairTableAlloc(numChr);

    return pNewInfoTable;
}

void TGM_ReadPairTableFree(TGM_ReadPairTable* pReadPairTable)
{
    if (pReadPairTable != NULL)
    {
        if (pReadPairTable->pLongPairArray != NULL)
        {
            free(pReadPairTable->pLongPairArray->chrCount);
            TGM_ARRAY_FREE(pReadPairTable->pLongPairArray, TRUE);
        }

        if (pReadPairTable->pShortPairArray != NULL)
        {
            free(pReadPairTable->pShortPairArray->chrCount);
            TGM_ARRAY_FREE(pReadPairTable->pShortPairArray, TRUE);
        }

        if (pReadPairTable->pReversedPairArray != NULL)
        {
            free(pReadPairTable->pReversedPairArray->chrCount);
            TGM_ARRAY_FREE(pReadPairTable->pReversedPairArray, TRUE);
        }

        if (pReadPairTable->pInvertedPairArray != NULL)
        {
            free(pReadPairTable->pInvertedPairArray->chrCount);
            TGM_ARRAY_FREE(pReadPairTable->pInvertedPairArray, TRUE);
        }

        if (pReadPairTable->pCrossPairArray != NULL)
        {
            free(pReadPairTable->pCrossPairArray->chrCount);
            TGM_ARRAY_FREE(pReadPairTable->pCrossPairArray, TRUE);
        }

        TGM_SpecialPairTableFree(pReadPairTable->pSpecialPairTable);

        free(pReadPairTable);
    }
}

TGM_SplitPairTable* TGM_SplitPairTableAlloc(uint32_t numChr, uint32_t detectSet)
{
    TGM_SplitPairTable* pSplitPairTable = (TGM_SplitPairTable*) malloc (sizeof(TGM_SplitPairTable));
    if (pSplitPairTable == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the split pair table object.\n");

    pSplitPairTable->numChr = numChr;
    pSplitPairTable->detectSet = detectSet;

    TGM_ARRAY_ALLOC(pSplitPairTable->pSplitLongArray, DEFAULT_RP_INFO_CAPACITY, TGM_SplitPairArray, TGM_SplitPair);
    pSplitPairTable->pSplitLongArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
    if (pSplitPairTable->pSplitLongArray->chrCount == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

    TGM_ARRAY_ALLOC(pSplitPairTable->pSplitShortArray, DEFAULT_RP_INFO_CAPACITY, TGM_SplitPairArray, TGM_SplitPair);
    pSplitPairTable->pSplitShortArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
    if (pSplitPairTable->pSplitShortArray->chrCount == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

    TGM_ARRAY_ALLOC(pSplitPairTable->pSplitReversedArray, DEFAULT_RP_INFO_CAPACITY, TGM_SplitPairArray, TGM_SplitPair);
    pSplitPairTable->pSplitReversedArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
    if (pSplitPairTable->pSplitReversedArray->chrCount == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

    TGM_ARRAY_ALLOC(pSplitPairTable->pSplitInvertedArray, DEFAULT_RP_INFO_CAPACITY, TGM_SplitPairArray, TGM_SplitPair);
    pSplitPairTable->pSplitInvertedArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
    if (pSplitPairTable->pSplitInvertedArray->chrCount == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

    TGM_ARRAY_ALLOC(pSplitPairTable->pSplitSpecialArray, DEFAULT_RP_INFO_CAPACITY, TGM_SplitPairArray, TGM_SplitPair);
    pSplitPairTable->pSplitSpecialArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
    if (pSplitPairTable->pSplitSpecialArray->chrCount == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

    return pSplitPairTable;
}

void TGM_SplitPairTableFree(TGM_SplitPairTable* pSplitPairTable)
{
    if (pSplitPairTable != NULL)
    {
        if (pSplitPairTable->pSplitLongArray != NULL)
        {
            free(pSplitPairTable->pSplitLongArray->chrCount);
            TGM_ARRAY_FREE(pSplitPairTable->pSplitLongArray, TRUE);
        }

        if (pSplitPairTable->pSplitShortArray != NULL)
        {
            free(pSplitPairTable->pSplitShortArray->chrCount);
            TGM_ARRAY_FREE(pSplitPairTable->pSplitShortArray, TRUE);
        }

        if (pSplitPairTable->pSplitReversedArray != NULL)
        {
            free(pSplitPairTable->pSplitReversedArray->chrCount);
            TGM_ARRAY_FREE(pSplitPairTable->pSplitReversedArray, TRUE);
        }

        if (pSplitPairTable->pSplitInvertedArray != NULL)
        {
            free(pSplitPairTable->pSplitInvertedArray->chrCount);
            TGM_ARRAY_FREE(pSplitPairTable->pSplitInvertedArray, TRUE);
        }

        if (pSplitPairTable->pSplitSpecialArray != NULL)
        {
            free(pSplitPairTable->pSplitSpecialArray->chrCount);
            TGM_ARRAY_FREE(pSplitPairTable->pSplitSpecialArray, TRUE);
        }

        free(pSplitPairTable);
    }
}


//======================
// Interface functions
//======================

void TGM_ReadPairBuild(const TGM_ReadPairBuildPars* pBuildPars)
{
    // some default capacity of the containers
    unsigned int capAnchor = 150;
    unsigned int capSample = 20;
    unsigned int capReadGrp = 20;
    unsigned int capHist = 10;


    // initialize the library information table
    TGM_LibInfoTable* pLibTable = TGM_LibInfoTableAlloc(capAnchor, capSample, capReadGrp, pBuildPars->specialPrefix);
    TGM_LibInfoTableSetCutoff(pLibTable, pBuildPars->cutoff);
    TGM_LibInfoTableSetTrimRate(pLibTable, pBuildPars->trimRate);

    // get the library table output file name
    char* libTableOutputFile = TGM_CreateFileName(pBuildPars->workingDir, TGM_LibTableFileName);

    // write the library information table into the file
    FILE* libTableOutput = fopen(libTableOutputFile, "wb");
    if (libTableOutput == NULL)
        TGM_ErrQuit("ERROR: Cannot open the library output file: %s\n", libTableOutputFile);

    free(libTableOutputFile);

    // structure initialization
    TGM_BamInStreamLite* pBamInStreamLite = TGM_BamInStreamLiteAlloc();

    // boolean variable control if we want to
    // write the fragment length info into the library information file
    TGM_Bool writeLibInfo = FALSE;

    // read pair file hash
    khash_t(file)* pFileHash = NULL;


    char* pSpecialID = NULL;

    // buffer used to hold the bam file name
    char bamFileName[TGM_MAX_LINE];

    if (pBuildPars->fileListInput != NULL) // if we are going to use the read pair information
    {
        TGM_FragLenHistArray* pHistArray = TGM_FragLenHistArrayAlloc(capHist);
        TGM_BamHeader* pBamHeader = NULL;

        // a flag to indicate if we have created the read pair table
        TGM_Bool hasReadPairTable = FALSE;
        TGM_ReadPairTable* pReadPairTable = NULL;

        // open the fragment length histogram output file
        char* histOutputFile = TGM_CreateFileName(pBuildPars->workingDir, TGM_HistFileName);
        FILE* histOutput = fopen(histOutputFile, "w");
        if (histOutput == NULL)
            TGM_ErrQuit("ERROR: Cannot open fragment length histogram file: %s\n", histOutputFile);

        free(histOutputFile);

        while (TGM_GetNextLine(bamFileName, TGM_MAX_LINE, pBuildPars->fileListInput) == TGM_OK)
        {
            // open the bam file
            TGM_BamInStreamLiteOpen(pBamInStreamLite, bamFileName);

            // load the bam header before read any alignments
            pBamHeader = TGM_BamInStreamLiteLoadHeader(pBamInStreamLite);

            // get the file position of the bam aligments
            int64_t bamPos = TGM_BamInStreamLiteTell(pBamInStreamLite);

            // process the header information
            unsigned int oldSize = 0;
            if (TGM_LibInfoTableSetRG(pLibTable, &oldSize, pBamHeader) != TGM_OK)
                TGM_ErrQuit("ERROR: Found an error when loading the bam file.\n");

            // initialize the fragment length histogram array with the number of newly added libraries in the bam file
            TGM_FragLenHistArrayInit(pHistArray, pLibTable->size - oldSize);

            // do not load cross pairs when build the fragment length distribution
            TGM_Bool loadCross = FALSE;

            // filter data for the no-za sort mode
            TGM_FilterDataNoZA filterData = {pLibTable, TRUE};

            // get the sorting order from the bam header
            TGM_SortMode sortMode = TGM_BamHeaderGetSortMode(pBamHeader);
            if (sortMode == TGM_SORTED_COORDINATE_NO_ZA || sortMode == TGM_SORTED_NAME || sortMode == TGM_SORTED_SPLIT)
            {
                // test if the bam file has za tag or not
                TGM_Status status = TGM_BamInStreamLiteTestZA(pBamInStreamLite);

                if (status == TGM_OK)
                    sortMode = TGM_SORTED_COORDINATE_ZA;
                else if (status == TGM_ERR)
                    continue;

                TGM_BamInStreamLiteSetSortMode(pBamInStreamLite, sortMode);
            }
            else
                TGM_ErrQuit("ERROR: Invalid sorting order.\n");

            // set the sort order for the bam instream
            if (sortMode != TGM_SORTED_COORDINATE_NO_ZA)
            {
                TGM_BamInStreamLiteSetFilter(pBamInStreamLite, TGM_ReadPairFilter);
                TGM_BamInStreamLiteSetFilterData(pBamInStreamLite, &loadCross);
            }
            else
            {
                TGM_BamInStreamLiteSetFilter(pBamInStreamLite, TGM_ReadPairNoZAFilter);
                TGM_BamInStreamLiteSetFilterData(pBamInStreamLite, &filterData);
            }

            int retNum = 0;
            const bam1_t* pAlgns[3] = {NULL, NULL, NULL};

            TGM_Status bamStatus = TGM_OK;

            // read the primary bam the first time to build fragment length distribution
            do
            {
                int64_t index = -1;
                bamStatus = TGM_BamInStreamLiteRead(pAlgns, &retNum, &index, pBamInStreamLite);

                if (retNum > 0)
                {
                    // check if the incoming read pair is normal (unique-unique pair)
                    // if yes, then update the corresponding fragment length histogram
                    TGM_PairStats pairStats;
                    TGM_ZAtag zaTag;
                    const TGM_MateInfo* pMateInfo = TGM_BamInStreamLiteGetMateInfo(pBamInStreamLite, index);

                    unsigned int backHistIndex = 0;
                    TGM_Status zaStatus = TGM_ERR;
                    if (TGM_IsNormalPair(&pairStats, &zaTag, &zaStatus, pMateInfo, &backHistIndex, pAlgns, retNum, pLibTable, pBuildPars->minMQ))
                        TGM_FragLenHistArrayUpdate(pHistArray, backHistIndex, pairStats.fragLen);
                }

            }while(bamStatus == TGM_OK);

            // finish the process of the histogram and update the library information table
            TGM_FragLenHistArrayFinalize(pHistArray);
            TGM_LibInfoTableUpdate(pLibTable, pHistArray, oldSize);

            // write the fragment length histogram into the file
            TGM_FragLenHistArrayWrite(pHistArray, histOutput);


            // clear the bam in stream
            TGM_BamInStreamLiteClear(pBamInStreamLite);

            // rewind the bam file to the beginning of the alignments (right after the header)
            TGM_BamInStreamLiteSeek(pBamInStreamLite, bamPos, SEEK_SET);

            // we only have to create the read pair table and open the read pair files once
            if (!hasReadPairTable)
            {
                pReadPairTable = TGM_ReadPairTableAlloc(pLibTable->pAnchorInfo->size, pBuildPars->detectSet); 
                pFileHash = TGM_ReadPairFilesOpen(pLibTable, pBuildPars->detectSet, pBuildPars->workingDir);
                hasReadPairTable =TRUE;
            }

            // we need to load the cross pair if we want to detect inter-chromosome translocation
            if ((pBuildPars->detectSet & SV_INTER_CHR_TRNSLCTN) != 0)
                loadCross = TRUE;
            else
                loadCross = FALSE;

            // turn off the keep normal flag
            // we are going to select those SV candidates now
            filterData.keepNormal = FALSE;

            // set the filter data for different sorting order
            if (sortMode != TGM_SORTED_COORDINATE_NO_ZA)
                TGM_BamInStreamLiteSetFilterData(pBamInStreamLite, &loadCross);
            else
                TGM_BamInStreamLiteSetFilterData(pBamInStreamLite, &filterData);

            do
            {
                // load the next read pair
                int64_t index = -1;
                bamStatus = TGM_BamInStreamLiteRead(pAlgns, &retNum, &index, pBamInStreamLite);

                if (retNum > 0)
                {
                    TGM_ZAtag zaTag;
                    TGM_PairStats pairStats;
                    const TGM_MateInfo* pMateInfo = TGM_BamInStreamLiteGetMateInfo(pBamInStreamLite, index);

                    const bam1_t* pUpAlgn = NULL;
                    const bam1_t* pDownAlgn = NULL;

                    if (retNum == 1)
                    {
                        if (pMateInfo != NULL)        // sorted by coordinate without ZA
                            pDownAlgn = pAlgns[0];
                        else                          // sorted by coordinate with ZA
                            pUpAlgn = pAlgns[0];
                    }
                    else if (retNum == 2)             // sorted by query name or this is a split bam file
                    {
                        pUpAlgn = pAlgns[0];
                        pDownAlgn = pAlgns[1];
                    }
                    else
                        continue;

                    // get the read pair information
                    TGM_Status readStatus = TGM_OK;
                    if (pMateInfo != NULL)
                        readStatus = TGM_LoadPairStats(&pairStats, pDownAlgn, pLibTable);
                    else
                        readStatus = TGM_LoadPairStats(&pairStats, pUpAlgn, pLibTable);

                    // update the read pair table
                    if (readStatus == TGM_OK)
                    {
                        if (pMateInfo != NULL) // no za tag
                        {
                            TGM_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, NULL, &pairStats, pMateInfo, pLibTable, pHistArray, pBuildPars);
                        }
                        else 
                        {
                            TGM_Status ZAstatus = TGM_LoadZAtag(&zaTag, pUpAlgn);

                            if (ZAstatus == TGM_OK)
                                TGM_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, &zaTag, &pairStats, NULL, pLibTable, pHistArray, pBuildPars);
                            else if (pDownAlgn != NULL)
                                TGM_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, NULL, &pairStats, NULL, pLibTable, pHistArray, pBuildPars);
                        }
                    }
                }

            }while(bamStatus == TGM_OK);

            // write out the read pair table for the current bam file
            TGM_ReadPairTableWrite(pReadPairTable, pFileHash);

            // clear the read pair table for the next bam file
            TGM_ReadPairTableClear(pReadPairTable);

            // close the bam file
            TGM_BamInStreamLiteClose(pBamInStreamLite);
            TGM_BamHeaderFree(pBamHeader);
        }

        // write the number of pairs into each read pair file
        TGM_ReadPairFilesWriteNum(pFileHash);

        if (pLibTable->size > 0)
        {
            // this function will check the empty condition or null pointer condition
            pSpecialID = TGM_SpecialPairTableGetID(pReadPairTable->pSpecialPairTable);
            writeLibInfo = TRUE;
        }

        // clean up
        fclose(histOutput);
        TGM_FragLenHistArrayFree(pHistArray);
        TGM_ReadPairTableFree(pReadPairTable);
    }

    // if we are going to use the split alignments
    if (pBuildPars->splitFileInput != NULL) 
    {
        TGM_Bool hasSplitTable = FALSE;
        TGM_SplitPairTable* pSplitPairTable = NULL;

        while (TGM_GetNextLine(bamFileName, TGM_MAX_LINE, pBuildPars->splitFileInput) == TGM_OK)
        {
            // open the bam file
            TGM_BamInStreamLiteOpen(pBamInStreamLite, bamFileName);

            // load the bam header
            TGM_BamHeader* pBamHeader = TGM_BamInStreamLiteLoadHeader(pBamInStreamLite);

            // load the bam header before read any alignments
            if (TGM_LibInfoTableSetRGSplit(pLibTable, pBamHeader, pBuildPars->specialPrefix, pBuildPars->prefixLen) != TGM_OK)
                TGM_ErrQuit("ERROR: Found an error when loading the split bam file.\n");

            // we have to change the fitler function here for split pairs
            TGM_BamInStreamLiteSetFilter(pBamInStreamLite, TGM_SplitFilter);
            
            // split filter does not need any parameters
            TGM_BamInStreamLiteSetFilterData(pBamInStreamLite, NULL);

            // set the sorting order
            TGM_SortMode sortMode = TGM_SORTED_SPLIT;
            TGM_BamInStreamLiteSetSortMode(pBamInStreamLite, sortMode);

            int retNum = 0;
            const bam1_t* pAlgns[3] = {NULL, NULL, NULL};
            TGM_Status bamStatus = TGM_OK;

            // create a split table object if we haven't yet
            if (!hasSplitTable)
            {
                uint32_t numChr = TGM_LibInfoTableCountNormalChr(pLibTable);
                pSplitPairTable = TGM_SplitPairTableAlloc(numChr, pBuildPars->detectSet);

                // open the read pair files 
                if (pFileHash == NULL)
                {
                    pFileHash = TGM_ReadPairFilesOpen(pLibTable, pBuildPars->detectSet, pBuildPars->workingDir);
                    TGM_ReadPairFilesWriteNum(pFileHash);
                }

                // initialize the number of split pairs
                TGM_ReadPairFilesWriteSplitNum(pFileHash);

                hasSplitTable = TRUE;
            }

            do
            {
                int64_t index = -1;
                bamStatus = TGM_BamInStreamLiteRead(pAlgns, &retNum, &index, pBamInStreamLite);

                if (retNum == 2)
                {
                    TGM_ZAtag zaTag;

                    const bam1_t* pFirstPartial = pAlgns[0];
                    const bam1_t* pSecondPartial = pAlgns[1];

                    // FIXME: ZA tag is not available in current split bam
                    // change this block back when the ZA tag is available
                    /*  
                    TGM_Status ZAstatus = TGM_LoadZAtag(&zaTag, pFirstPartial);

                    if (ZAstatus == TGM_OK)
                        TGM_SplitPairTableUpdate(pSplitPairTable, pLibTable, pFirstPartial, pSecondPartial, &zaTag, pBuildPars->minMQ);
                    */
                    TGM_SplitPairTableUpdate(pSplitPairTable, pLibTable, pFirstPartial, pSecondPartial, &zaTag, pBuildPars->minMQ);
                }

            }while(bamStatus == TGM_OK);

            // write the split pair inforamtion into the read pair files
            TGM_SplitPairTableWrite(pSplitPairTable, pFileHash);
            TGM_SplitPairTableClear(pSplitPairTable);

            // close the current split bam file and free the bam header
            TGM_BamInStreamLiteClose(pBamInStreamLite);
            TGM_BamHeaderFree(pBamHeader);
        }

        // write the number of split pairs of each type into the read pair file
        TGM_ReadPairFilesWriteSplitNum(pFileHash);

        // free the split pair table
        TGM_SplitPairTableFree(pSplitPairTable);
    }

    // write the library information into file
    TGM_LibInfoTableWrite(pLibTable, writeLibInfo, libTableOutput);
    // write the detect set into the library information file
    TGM_ReadPairTableWriteDetectSet(pBuildPars, libTableOutput);

    // write the special reference ID into the library information file
    if (pSpecialID != NULL)
    {
        TGM_SpecialPairTableWriteID(pSpecialID, libTableOutput);
        free(pSpecialID);
    }
    else
    {
        uint32_t zero = 0;
        fwrite(&zero, sizeof(uint32_t), 1, libTableOutput);
    }

    // clean up
    fclose(libTableOutput);
    TGM_ReadPairFilesClose(pFileHash);
    TGM_LibInfoTableFree(pLibTable);
    TGM_BamInStreamLiteFree(pBamInStreamLite);
}


/*  
// read the bam alignments from the primary bam files, 
// extract the SV-related information out and wirte it into 
// read pair files according to the read pair type and 
// location
void TGM_ReadPairBuild(const TGM_ReadPairBuildPars* pBuildPars)
{
    // some default capacity of the containers
    unsigned int buffCapacity = 100;
    unsigned int capAnchor = 150;
    unsigned int capSample = 20;
    unsigned int capReadGrp = 20;
    unsigned int capHist = 10;

    // phony parameters for the bam input stream
    unsigned int reportSize = 0;
    unsigned int numThread = 0;

    // initialize the library information table
    TGM_LibInfoTable* pLibTable = TGM_LibInfoTableAlloc(capAnchor, capSample, capReadGrp);
    TGM_LibInfoTableSetCutoff(pLibTable, pBuildPars->cutoff);
    TGM_LibInfoTableSetTrimRate(pLibTable, pBuildPars->trimRate);

    // this is the data used to filter the read pairs in the bam
    TGM_FilterDataRP* pFilterData = TGM_FilterDataRPAlloc(pLibTable->pAnchorInfo, pBuildPars->binLen);

    // set the stream mode to read pair filter
    TGM_StreamMode streamMode;
    TGM_SetStreamMode(&streamMode, TGM_ReadPairFilter, pFilterData, TGM_READ_PAIR_MODE);

    // structure initialization
    TGM_BamInStream* pBamInStream = TGM_BamInStreamAlloc(pBuildPars->binLen, numThread, buffCapacity, reportSize, &streamMode);
    TGM_FragLenHistArray* pHistArray = TGM_FragLenHistArrayAlloc(capHist);
    TGM_BamHeader* pBamHeader = NULL;

    // a flag to indicate if we have created the read pair table
    TGM_Bool hasReadPairTable = FALSE;
    TGM_ReadPairTable* pReadPairTable = NULL;

    // file hash
    khash_t(file)* pFileHash = NULL;

    // required arguments for bam in stream structure
    char* histOutputFile = TGM_CreateFileName(pBuildPars->workingDir, TGM_HistFileName);

    // open the fragment length histogram output file
    FILE* histOutput = fopen(histOutputFile, "w");
    if (histOutput == NULL)
        TGM_ErrQuit("ERROR: Cannot open fragment length histogram file: %s\n", histOutputFile);

    char bamFileName[TGM_MAX_LINE];
    while (TGM_GetNextLine(bamFileName, TGM_MAX_LINE, pBuildPars->fileListInput) == TGM_OK)
    {
        // update the fitler data with new bam file
        // this will affect the bam in stream
        TGM_FilterDataRPInit(pFilterData, bamFileName);
        TGM_FilterDataRPTurnOffCross(pFilterData);

        // open the bam file
        TGM_BamInStreamOpen(pBamInStream, bamFileName);

        // load the bam header before read any alignments
        pBamHeader = TGM_BamInStreamLoadHeader(pBamInStream);

        // get the file position of the bam aligments
        int64_t bamPos = TGM_BamInStreamTell(pBamInStream);

        // process the header information
        unsigned int oldSize = 0;
        if (TGM_LibInfoTableSetRG(pLibTable, &oldSize, pBamHeader) != TGM_OK)
            TGM_ErrQuit("ERROR: Found an error when loading the bam file.\n");

        // initialize the fragment length histogram array with the number of newly added libraries in the bam file
        TGM_FragLenHistArrayInit(pHistArray, pLibTable->size - oldSize);

        TGM_BamNode* pUpNode = NULL;
        TGM_BamNode* pDownNode = NULL;

        const bam1_t* pUpAlgn = NULL;
        const bam1_t* pDownAlgn = NULL;

        TGM_Status bamStatus = TGM_OK;
        while ((bamStatus = TGM_BamInStreamLoadPair(&pUpNode, &pDownNode, pBamInStream)) != TGM_EOF && bamStatus != TGM_ERR)
        {
            // we hit another chromosome
            if (bamStatus == TGM_OUT_OF_RANGE)
                continue;

            // load the read pairs
            if (!pFilterData->isFilled)
            {
                // usually we should catch the read pairs here
                // since their fragment length should be smaller than the bin length
                pUpAlgn = &(pUpNode->alignment);
                pDownAlgn = &(pDownNode->alignment);
            }
            else
            {
                // if the fragment length of the read pair is beyond our bin length
                // we should catch them here
                pUpAlgn = pFilterData->pUpAlgn;
                pDownAlgn = pFilterData->pDownAlgn;
            }

            // check if the incoming read pair is normal (unique-unique pair)
            // if yes, then update the corresponding fragment length histogram
            TGM_PairStats pairStats;
            unsigned int backHistIndex = 0;
            if (TGM_IsNormalPair(&pairStats, &backHistIndex, pUpAlgn, pDownAlgn, pLibTable, pBuildPars->minMQ))
                TGM_FragLenHistArrayUpdate(pHistArray, backHistIndex, pairStats.fragLen);

            // recycle those bam nodes that are allocated from the memory pool
            if (!pFilterData->isFilled)
            {
                TGM_BamInStreamRecycle(pBamInStream, pUpNode);
                TGM_BamInStreamRecycle(pBamInStream, pDownNode);
            }
        }

        // finish the process of the histogram and update the library information table
        TGM_FragLenHistArrayFinalize(pHistArray);
        TGM_LibInfoTableUpdate(pLibTable, pHistArray, oldSize);

        // write the fragment length histogram into the file
        TGM_FragLenHistArrayWrite(pHistArray, histOutput);


        // clear the bam in stream
        TGM_BamInStreamClear(pBamInStream);

        // rewind the bam file to the beginning of the alignments (right after the header)
        TGM_BamInStreamSeek(pBamInStream, bamPos, SEEK_SET);

        // we only have to create the read pair table and open the read pair files once
        if (!hasReadPairTable)
        {
            pReadPairTable = TGM_ReadPairTableAlloc(pLibTable->pAnchorInfo->size, pBuildPars->detectSet); 
            pFileHash = TGM_ReadPairFilesOpen(pLibTable, pBuildPars->detectSet, pBuildPars->workingDir);
            hasReadPairTable =TRUE;
        }

        // we need to load the cross pair if we want to detect inter-chromosome translocation
        if ((pBuildPars->detectSet & SV_INTER_CHR_TRNSLCTN) != 0)
            TGM_FilterDataRPTurnOnCross(pFilterData);

        while ((bamStatus = TGM_BamInStreamLoadPair(&pUpNode, &pDownNode, pBamInStream)) != TGM_EOF && bamStatus != TGM_ERR)
        {
            // we hit another chromosome
            if (bamStatus == TGM_OUT_OF_RANGE)
                continue;

            // load the read pairs
            if (!pFilterData->isFilled)
            {
                // usually we should catch the read pairs here
                // since their fragment length should be smaller than the bin length
                pUpAlgn = &(pUpNode->alignment);
                pDownAlgn = &(pDownNode->alignment);
            }
            else
            {
                // if the fragment length of the read pair is beyond our bin length
                // we should catch them here
                pUpAlgn = pFilterData->pUpAlgn;
                pDownAlgn = pFilterData->pDownAlgn;
            }

            TGM_ZAtag zaTag;
            TGM_PairStats pairStats;

            TGM_Status readStatus = TGM_LoadPairStats(&pairStats, pUpAlgn, pLibTable);

            if (readStatus == TGM_OK)
            {
                TGM_Status ZAstatus = TGM_LoadZAtag(&zaTag, pUpAlgn);

                if (ZAstatus == TGM_OK)
                    TGM_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, &zaTag, &pairStats, pLibTable, pHistArray, pBuildPars->minMQ);
                else
                    TGM_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, NULL, &pairStats, pLibTable, pHistArray, pBuildPars->minMQ);
            }


            // recycle those bam nodes that are allocated from the memory pool
            if (!pFilterData->isFilled)
            {
                TGM_BamInStreamRecycle(pBamInStream, pUpNode);
                TGM_BamInStreamRecycle(pBamInStream, pDownNode);
            }
        }

        TGM_ReadPairTableWrite(pReadPairTable, pFileHash);
        TGM_ReadPairTableClear(pReadPairTable);

        // close the bam file
        TGM_BamInStreamClose(pBamInStream);
        TGM_BamHeaderFree(pBamHeader);
    }

    // close all the read pair files
    TGM_ReadPairFilesClose(pFileHash);

    // get the library table output file name
    char* libTableOutputFile = TGM_CreateFileName(pBuildPars->workingDir, TGM_LibTableFileName);

    // write the library information table into the file
    FILE* libTableOutput = fopen(libTableOutputFile, "wb");
    if (libTableOutput == NULL)
        TGM_ErrQuit("ERROR: Cannot open the library output file: %s\n", libTableOutputFile);

    TGM_LibInfoTableWrite(pLibTable, libTableOutput);
    TGM_ReadPairTableWriteDetectSet(pReadPairTable, libTableOutput);
    TGM_SpecialPairTableWirteID(pReadPairTable->pSpecialPairTable, libTableOutput);

    // clean up
    fclose(histOutput);
    fclose(libTableOutput);

    free(libTableOutputFile);
    free(histOutputFile);

    TGM_FragLenHistArrayFree(pHistArray);
    TGM_LibInfoTableFree(pLibTable);
    TGM_FilterDataRPFree(pFilterData);
    TGM_ReadPairTableFree(pReadPairTable);
    TGM_BamInStreamFree(pBamInStream);
}

*/

// clear the read pair table
void TGM_ReadPairTableClear(TGM_ReadPairTable* pReadPairTable)
{
    if (pReadPairTable->pLongPairArray != NULL)
    {
        memset(pReadPairTable->pLongPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pLongPairArray->size = 0;
    }

    if (pReadPairTable->pShortPairArray != NULL)
    {
        memset(pReadPairTable->pShortPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pShortPairArray->size = 0;
    }

    if (pReadPairTable->pReversedPairArray != NULL)
    {
        memset(pReadPairTable->pReversedPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pReversedPairArray->size = 0;
    }

    if (pReadPairTable->pInvertedPairArray != NULL)
    {
        memset(pReadPairTable->pInvertedPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pInvertedPairArray->size = 0;
    }

    if (pReadPairTable->pCrossPairArray != NULL)
    {
        memset(pReadPairTable->pCrossPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pCrossPairArray->size = 0;
    }

    if (pReadPairTable->pSpecialPairTable != NULL)
        TGM_SpecialPairTableClear(pReadPairTable->pSpecialPairTable, pReadPairTable->numChr);
}

// update the read pair table with the incoming read pairs
void TGM_ReadPairTableUpdate(TGM_ReadPairTable* pReadPairTable, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const TGM_ZAtag* pZAtag, const TGM_PairStats* pPairStats, 
                            const TGM_MateInfo* pMateInfo, const TGM_LibInfoTable* pLibTable, const TGM_FragLenHistArray* pHistArray, const TGM_ReadPairBuildPars* pBuildPars)
{
    // check the read pair type
    unsigned char minMQ = pBuildPars->minMQ;
    unsigned char minSpMQ = pBuildPars->minSpMQ;

    SV_ReadPairType readPairType = PT_NORMAL;
    readPairType = TGM_CheckReadPairType(pZAtag, pPairStats, pDownAlgn, pMateInfo, minMQ, minSpMQ, pLibTable);

    // check the mapping quality for all read pairs except the special read pairs
    if (readPairType != PT_SPECIAL5 && readPairType != PT_SPECIAL3)
    {
        if (pDownAlgn != NULL)
        {
            if (pMateInfo != NULL)
            {
                if (pMateInfo->mapQ < minMQ || pDownAlgn->core.qual < minMQ)
                    readPairType = PT_UNKNOWN;
            }
            else
            {
                if (pUpAlgn->core.qual < minMQ || pDownAlgn->core.qual < minMQ)
                    readPairType = PT_UNKNOWN;
            }
        }
        else
        {
            if (pZAtag->bestMQ[0] < minMQ || pZAtag->bestMQ[1] < minMQ)
                readPairType = PT_UNKNOWN;
        }
    }

    switch (readPairType)
    {
        case PT_NORMAL:
        case PT_UNKNOWN:
            break;
        case PT_LONG:
            if (pReadPairTable->pLongPairArray != NULL)
                TGM_LocalPairArrayUpdate(pReadPairTable->pLongPairArray, pUpAlgn, pDownAlgn, pZAtag, pPairStats, pMateInfo, pLibTable, pHistArray, readPairType);
            break;
        case PT_SHORT:
            if (pReadPairTable->pShortPairArray != NULL)
                TGM_LocalPairArrayUpdate(pReadPairTable->pShortPairArray, pUpAlgn, pDownAlgn, pZAtag, pPairStats, pMateInfo, pLibTable, pHistArray, readPairType);
            break;
        case PT_REVERSED:
            if (pReadPairTable->pReversedPairArray != NULL)
                TGM_LocalPairArrayUpdate(pReadPairTable->pReversedPairArray, pUpAlgn, pDownAlgn, pZAtag, pPairStats, pMateInfo, pLibTable, pHistArray, readPairType);
            break;
        case PT_INVERTED3:
        case PT_INVERTED5:
            if (pReadPairTable->pInvertedPairArray != NULL)
                TGM_InvertedPairArrayUpdate(pReadPairTable->pInvertedPairArray, pUpAlgn, pDownAlgn, pZAtag, pPairStats, pMateInfo, readPairType);
            break;
        case PT_SPECIAL3:
        case PT_SPECIAL5:
            if (pReadPairTable->pSpecialPairTable != NULL) 
                TGM_SpecialPairTableUpdate(pReadPairTable->pSpecialPairTable, pUpAlgn, pDownAlgn, pLibTable, pHistArray, pPairStats, pZAtag, readPairType);
            break;
        case PT_CROSS:
            if (pReadPairTable->pCrossPairArray != NULL)
                TGM_CrossPairArrayUpdate(pReadPairTable->pCrossPairArray, pUpAlgn, pDownAlgn, pZAtag, pPairStats);
            break;
        default:
            break;
    }
}

// open a seriers read pair files for output. A read pair
// file will be open for each read pair type on each chromosome
void* TGM_ReadPairFilesOpen(const TGM_LibInfoTable* pLibTable, uint32_t detectSet, const char* workingDir)
{
    khash_t(file)* pFileHash = kh_init(file);

    // buffer to hold the reference ID
    char refBuff[4];

    // length of the working directory
    int dirLen = strlen(workingDir);

    // buffer to hold the full read pair file name
    char nameBuff[dirLen + 50];

    // copy the working directory name into the name buffer
    strncpy(nameBuff, workingDir, dirLen);
    char* namePos = nameBuff + dirLen;
    char* replacePos = namePos + 3;

    // the reference ID should not be larger than 1000 (arbitrary, but reasonable)
    const TGM_AnchorInfo* pAnchorInfo = pLibTable->pAnchorInfo;
    if (pAnchorInfo->size > 1000)
        TGM_ErrQuit("ERROR: reference ID overflow (0-999).\n");

    // this is the place holder put at the beginning of each read pair file
    // different SV types may have different file headers
    const uint64_t placeHolder = 0;
    TGM_ReadPairOutStream stream = {NULL, 0, 0};

    for (unsigned int i = 0; i != pAnchorInfo->size; ++i)
    {
        if (pAnchorInfo->pLength[i] <= 0
            || strncmp(pAnchorInfo->pAnchors[i], pAnchorInfo->pSpecialPrefix, pAnchorInfo->specialPrefixLen) == 0)
        {
            continue;
        }

        // copy the reference ID into the name buffer
        unsigned int offset = 2 - (unsigned int) log10(i);
        sprintf(refBuff, "000");
        sprintf(refBuff + offset, "%d", i);

        int fileID = 0;
        int ret = 0;
        khiter_t khIter = 0;
        unsigned int fileIndex = TGM_LONG_PAIR_FILE;

        for (unsigned int j = SV_DELETION; j <= SV_INTER_CHR_TRNSLCTN; ++j)
        {
            if ((detectSet & (1 << j)) != 0)
            {
                strcpy(namePos, TGM_READ_PAIR_FILE_NAME_TEMPLATE[fileIndex]);
                strncpy(replacePos, refBuff, 3);

                stream.output = fopen(nameBuff, "wb");
                if (stream.output == NULL)
                    TGM_ErrQuit("ERROR: Cannot open file: %d\n", nameBuff);

                fileID = i * TGM_NUM_RP_FILE + fileIndex; 
                khIter = kh_put(file, pFileHash, fileID, &ret);
                kh_value(pFileHash, khIter) = stream;
                // leave a place for the total number of read pairs
                fwrite(&placeHolder, sizeof(uint64_t), 1, stream.output);

                // for tandem duplication we have two read types to take care
                if (j == SV_TANDEM_DUP)
                {
                    ++fileIndex;

                    strcpy(namePos, TGM_READ_PAIR_FILE_NAME_TEMPLATE[fileIndex]);
                    strncpy(replacePos, refBuff, 3);

                    stream.output = fopen(nameBuff, "wb");
                    if (stream.output == NULL)
                        TGM_ErrQuit("ERROR: Cannot open file: %d\n", nameBuff);

                    fileID = i * TGM_NUM_RP_FILE + fileIndex; 
                    khIter = kh_put(file, pFileHash, fileID, &ret);
                    kh_value(pFileHash, khIter) = stream;

                    // leave a place for the total number of read pairs
                    fwrite(&placeHolder, sizeof(uint64_t), 1, stream.output);
                }
            }

            ++fileIndex;
        }
    }

    return ((void*) pFileHash);
}

// write the number of read pairs for each chromosome and SV types into the beginning of each read pair file
void TGM_ReadPairFilesWriteNum(void* pFileHash)
{
    khash_t(file)* pHash = (khash_t(file)*) pFileHash;

    for (khiter_t khIter = kh_begin(pHash); khIter != kh_end(pHash); ++khIter)
    {
        if (kh_exist(pHash, khIter))
        {
            kh_value(pHash, khIter).splitInfoPos = ftello(kh_value(pHash, khIter).output);

            fseeko(kh_value(pHash, khIter).output, 0, SEEK_SET);
            fwrite(&(kh_value(pHash, khIter).numPairs), sizeof(uint64_t), 1, kh_value(pHash, khIter).output);

            kh_value(pHash, khIter).numPairs = 0;
        }
    }
}

// close all the read pair files
void TGM_ReadPairFilesWriteSplitNum(void* pFileHash)
{
    khash_t(file)* pHash = (khash_t(file)*) pFileHash;

    for (khiter_t khIter = kh_begin(pHash); khIter != kh_end(pHash); ++khIter)
    {
        if (kh_exist(pHash, khIter))
        {
            fseeko(kh_value(pHash, khIter).output, kh_value(pHash, khIter).splitInfoPos, SEEK_SET);
            fwrite(&(kh_value(pHash, khIter).numPairs), sizeof(uint64_t), 1, kh_value(pHash, khIter).output);
        }
    }
}

// close all the read pair files
void TGM_ReadPairFilesClose(void* pFileHash)
{
    khash_t(file)* pHash = (khash_t(file)*) pFileHash;

    for (khiter_t khIter = kh_begin(pHash); khIter != kh_end(pHash); ++khIter)
    {
        if (kh_exist(pHash, khIter))
            fclose(kh_value(pHash, khIter).output);
    }

    kh_destroy(file, pHash);
}

// write the read pair table into the read pair files
void TGM_ReadPairTableWrite(const TGM_ReadPairTable* pReadPairTable, void* pHash)
{
    khash_t(file)* pFileHash = (khash_t(file)*) pHash;
    unsigned int numChr = pReadPairTable->numChr;

    int fileID = 0;
    khiter_t khIter = 0;

    int currPos[] = {0, 0, 0, 0, 0, 0};

    for (unsigned int i = 0; i != numChr; ++i)
    {
        int numPairs = pReadPairTable->pLongPairArray == NULL ? 0 : pReadPairTable->pLongPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_LONG_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pLongPairArray->data + currPos[TGM_LONG_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_LONG_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pShortPairArray == NULL ? 0 : pReadPairTable->pShortPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_SHORT_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pShortPairArray->data + currPos[TGM_SHORT_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_SHORT_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pReversedPairArray == NULL ? 0 : pReadPairTable->pReversedPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_REVERSED_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pReversedPairArray->data + currPos[TGM_REVERSED_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_REVERSED_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pInvertedPairArray == NULL ? 0 : pReadPairTable->pInvertedPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_INVERTED_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pInvertedPairArray->data + currPos[TGM_INVERTED_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_INVERTED_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pSpecialPairTable == NULL ? 0 : pReadPairTable->pSpecialPairTable->array.chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_SPEICAL_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pSpecialPairTable->array.data + currPos[TGM_SPEICAL_PAIR_FILE], sizeof(TGM_SpecialPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_SPEICAL_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pCrossPairArray == NULL ? 0 : pReadPairTable->pCrossPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_CROSS_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pCrossPairArray->data + currPos[TGM_CROSS_PAIR_FILE], sizeof(TGM_CrossPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_CROSS_PAIR_FILE] += numPairs;
        }
    }

    // take care of those special pairs whose anchor refID is not the same as its mate's refID.
    if (pReadPairTable->pSpecialPairTable != NULL)
    {
        const TGM_SpecialPairArray* pSpecialPairArray = &(pReadPairTable->pSpecialPairTable->crossArray);
        for (unsigned int i = 0; i != pSpecialPairArray->size; ++i)
        {
            fileID = pSpecialPairArray->data[i].refID[0] * TGM_NUM_RP_FILE + TGM_SPEICAL_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += 1;
                fwrite(pSpecialPairArray->data + i, sizeof(TGM_SpecialPair), 1, kh_value(pFileHash, khIter).output);
            }
        }
    }
}

void TGM_SplitPairTableClear(TGM_SplitPairTable* pSplitPairTable)
{
    if (pSplitPairTable->pSplitLongArray != NULL)
    {
        memset(pSplitPairTable->pSplitLongArray->chrCount, 0, sizeof(uint64_t) * pSplitPairTable->numChr);
        pSplitPairTable->pSplitLongArray->size = 0;
    }

    if (pSplitPairTable->pSplitShortArray != NULL)
    {
        memset(pSplitPairTable->pSplitShortArray->chrCount, 0, sizeof(uint64_t) * pSplitPairTable->numChr);
        pSplitPairTable->pSplitShortArray->size = 0;
    }

    if (pSplitPairTable->pSplitReversedArray != NULL)
    {
        memset(pSplitPairTable->pSplitReversedArray->chrCount, 0, sizeof(uint64_t) * pSplitPairTable->numChr);
        pSplitPairTable->pSplitReversedArray->size = 0;
    }

    if (pSplitPairTable->pSplitInvertedArray != NULL)
    {
        memset(pSplitPairTable->pSplitInvertedArray->chrCount, 0, sizeof(uint64_t) * pSplitPairTable->numChr);
        pSplitPairTable->pSplitInvertedArray->size = 0;
    }

    if (pSplitPairTable->pSplitSpecialArray != NULL)
    {
        memset(pSplitPairTable->pSplitSpecialArray->chrCount, 0, sizeof(uint64_t) * pSplitPairTable->numChr);
        pSplitPairTable->pSplitSpecialArray->size = 0;
    }
}

// write the read pair table into the read pair files
void TGM_SplitPairTableWrite(const TGM_SplitPairTable* pSplitPairTable, void* pHash)
{
    khash_t(file)* pFileHash = (khash_t(file)*) pHash;
    unsigned int numChr = pSplitPairTable->numChr;

    int fileID = 0;
    khiter_t khIter = 0;

    int currPos[] = {0, 0, 0, 0, 0};

    for (unsigned int i = 0; i != numChr; ++i)
    {
        int numPairs = pSplitPairTable->pSplitLongArray == NULL ? 0 : pSplitPairTable->pSplitLongArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_LONG_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pSplitPairTable->pSplitLongArray->data + currPos[TGM_LONG_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_LONG_PAIR_FILE] += numPairs;
        }

        numPairs = pSplitPairTable->pSplitShortArray == NULL ? 0 : pSplitPairTable->pSplitShortArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_SHORT_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pSplitPairTable->pSplitShortArray->data + currPos[TGM_SHORT_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_SHORT_PAIR_FILE] += numPairs;
        }

        numPairs = pSplitPairTable->pSplitReversedArray == NULL ? 0 : pSplitPairTable->pSplitReversedArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_REVERSED_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pSplitPairTable->pSplitReversedArray->data + currPos[TGM_REVERSED_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_REVERSED_PAIR_FILE] += numPairs;
        }

        numPairs = pSplitPairTable->pSplitInvertedArray == NULL ? 0 : pSplitPairTable->pSplitInvertedArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_INVERTED_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pSplitPairTable->pSplitInvertedArray->data + currPos[TGM_INVERTED_PAIR_FILE], sizeof(TGM_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_INVERTED_PAIR_FILE] += numPairs;
        }

        numPairs = pSplitPairTable->pSplitSpecialArray == NULL ? 0 : pSplitPairTable->pSplitSpecialArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * TGM_NUM_RP_FILE + TGM_SPEICAL_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pSplitPairTable->pSplitSpecialArray->data + currPos[TGM_SPEICAL_PAIR_FILE], sizeof(TGM_SpecialPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[TGM_SPEICAL_PAIR_FILE] += numPairs;
        }
    }
}

void TGM_ReadPairTableWriteDetectSet(const TGM_ReadPairBuildPars* pBuildPars, FILE* libOutput)
{
    fwrite(&(pBuildPars->detectSet), sizeof(uint32_t), 1, libOutput);
}

// write the special reference name into the end of the library information file
void TGM_SpecialPairTableWriteID(const char* pSpecialID, FILE* libOutput)
{
    uint32_t size = strlen(pSpecialID) / 2;

    // write the name of special reference at the end of the library file
    fwrite(&(size), sizeof(uint32_t), 1, libOutput);
    fwrite(pSpecialID, sizeof(char), size * 2, libOutput);
}

void TGM_SplitPairTableUpdate(TGM_SplitPairTable* pSplitPairTable, const TGM_LibInfoTable* pLibTable, const bam1_t* pFirstPartial, 
                             const bam1_t* pSecondPartial, const TGM_ZAtag* pZAtag, uint8_t minMQ)
{
    if (pFirstPartial->core.qual < minMQ || pSecondPartial->core.qual < minMQ)
        return;

    int32_t readGrpID = 0;

    static const char tagRG[2] = {'R', 'G'};
    uint8_t* rgPos = bam_aux_get(pFirstPartial, tagRG);
    if (rgPos != NULL)
    {
        const char* RG = bam_aux2Z(rgPos);
        TGM_Status status = TGM_LibInfoTableGetRGIndex(&readGrpID, pLibTable, RG);
        if (status != TGM_OK)
            return;
    }
    else
        return;

    TGM_SplitPairArray* pSplitPairArray = NULL;
    switch (pZAtag->readPairType)
    {
        case PT_LONG:
            pSplitPairArray = pSplitPairTable->pSplitLongArray;
            break;
        case PT_SHORT:
            pSplitPairArray = pSplitPairTable->pSplitShortArray;
            break;
        case PT_REVERSED:
            pSplitPairArray = pSplitPairTable->pSplitReversedArray;
            break;
        case PT_INVERTED3:
        case PT_INVERTED5:
            pSplitPairArray = pSplitPairTable->pSplitInvertedArray;
            break;
        case PT_SPECIAL3:
        case PT_SPECIAL5:
            pSplitPairArray = pSplitPairTable->pSplitSpecialArray;
            break;
        default:
            return;
            break;
    }

    TGM_SplitPairArrayUpdate(pSplitPairArray, pFirstPartial, pSecondPartial, readGrpID, pZAtag);
}
