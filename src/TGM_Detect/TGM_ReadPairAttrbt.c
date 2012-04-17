/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairData.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/28/2012 10:40:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "TGM_Error.h"
#include "TGM_Utilities.h"
#include "TGM_ReadPairAttrbt.h"

#define DEFAULT_RP_ATTRB_CAPACITY 50

static int CompareAttrbt(const void* a, const void* b)
{
    const TGM_ReadPairAttrbt* first = a;
    const TGM_ReadPairAttrbt* second = b;

    if (first->firstAttribute < second->firstAttribute)
    {
        return -1;
    }
    else if (first->firstAttribute > second->firstAttribute)
    {
        return 1;
    }
    else
    {
        if (first->secondAttribute < second->secondAttribute)
            return -1;
        else if (first->secondAttribute > second->secondAttribute)
            return 1;
        else
        {
            if (first->readGrpID <= second->readGrpID)
                return -1;
            else
                return 1;
        }
    }
}


TGM_ReadPairAttrbtArray* TGM_ReadPairAttrbtArrayAlloc(uint32_t numReadGrp)
{
    TGM_ReadPairAttrbtArray* pAttrbtArray = (TGM_ReadPairAttrbtArray*) malloc(sizeof(TGM_ReadPairAttrbtArray));
    if (pAttrbtArray == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a read pair attribute array object.\n");

    pAttrbtArray->data = (TGM_ReadPairAttrbt*) malloc(DEFAULT_RP_ATTRB_CAPACITY * sizeof(TGM_ReadPairAttrbt));
    if (pAttrbtArray->data == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the read pair attributes in the read pair attribute array object.\n");

    pAttrbtArray->pBoundaries = (double (*)[2]) malloc(sizeof(double) * 2 * numReadGrp);
    if (pAttrbtArray->pBoundaries == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of the boundaries in the read pair attribute array object.\n");

    pAttrbtArray->numReadGrp = numReadGrp;
    pAttrbtArray->size = 0;
    pAttrbtArray->capacity = DEFAULT_RP_ATTRB_CAPACITY;

    return pAttrbtArray;
}

void TGM_ReadPairAttrbtArrayFree(TGM_ReadPairAttrbtArray* pAttrbtArray)
{
    if (pAttrbtArray != NULL)
    {
        free(pAttrbtArray->data);
        free(pAttrbtArray->pBoundaries);
        free(pAttrbtArray);
    }
}


void TGM_ReadPairAttrbtArrayReInit(TGM_ReadPairAttrbtArray* pAttrbtArray, uint64_t newCapacity)
{
    pAttrbtArray->size = 0;

    if (newCapacity > pAttrbtArray->capacity)
    {
        free(pAttrbtArray->data);
        pAttrbtArray->data = (TGM_ReadPairAttrbt*) malloc(newCapacity * sizeof(TGM_ReadPairAttrbt));
        if (pAttrbtArray->data == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the read pair attributes in the read pair attribute array object.\n");

        pAttrbtArray->capacity = newCapacity;
    }
}

void TGM_ReadPairMakeLocal(TGM_ReadPairAttrbtArray* pAttrbtArray, const TGM_LocalPairArray* pLocalPairArray, const TGM_LibInfoTable* pLibTable, SV_ReadPairType readPairType)
{
    TGM_ReadPairAttrbtArrayReInit(pAttrbtArray, pLocalPairArray->size);
    pAttrbtArray->readPairType = readPairType;
    pAttrbtArray->size = pLocalPairArray->size;
    pAttrbtArray->numReadGrp = pLibTable->size;

    for (unsigned int i = 0; i != pLocalPairArray->size; ++i)
    {
        pAttrbtArray->data[i].origIndex = i;
        pAttrbtArray->data[i].readGrpID = pLocalPairArray->data[i].readGrpID;

        double median = pLibTable->pLibInfo[pAttrbtArray->data[i].readGrpID].fragLenMedian;

        pAttrbtArray->data[i].firstAttribute =  pLocalPairArray->data[i].upPos + (double) pLocalPairArray->data[i].fragLen / 2;
        pAttrbtArray->data[i].secondAttribute = pLocalPairArray->data[i].fragLen - median;
    }

    qsort(pAttrbtArray->data, pAttrbtArray->size, sizeof(pAttrbtArray->data[0]), CompareAttrbt);

    // neighbourhood scale factor (borrowed from Spanner)
    double boundScale[2] = {1.25, 0.5};
    for (unsigned int i = 0; i != pAttrbtArray->numReadGrp; ++i)
    {
        pAttrbtArray->pBoundaries[i][0] = (double) pLibTable->pLibInfo[i].fragLenMedian * boundScale[0];
        pAttrbtArray->pBoundaries[i][1] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale[1];
    }

    qsort(pAttrbtArray->data, pAttrbtArray->size, sizeof(pAttrbtArray->data[0]), CompareAttrbt);
}

void TGM_ReadPairMakeCross(TGM_ReadPairAttrbtArray* pAttrbtArray, const TGM_LibInfoTable* pLibTable, const TGM_CrossPairArray* pCrossPairArray)
{
    TGM_ReadPairAttrbtArrayReInit(pAttrbtArray, pCrossPairArray->size);
    pAttrbtArray->readPairType = PT_CROSS;
    pAttrbtArray->size = pCrossPairArray->size;
    pAttrbtArray->numReadGrp = pLibTable->size;

    for (unsigned int i = 0; i != pCrossPairArray->size; ++i)
    {
        pAttrbtArray->data[i].origIndex = i;
        pAttrbtArray->data[i].readGrpID = pCrossPairArray->data[i].readGrpID;

        pAttrbtArray->data[i].firstAttribute = pCrossPairArray->data[i].upRefID * 1e10 + pCrossPairArray->data[i].upPos;
        pAttrbtArray->data[i].secondAttribute = pCrossPairArray->data[i].downRefID * 1e10 + pCrossPairArray->data[i].downPos;
    }

    qsort(pAttrbtArray->data, pAttrbtArray->size, sizeof(pAttrbtArray->data[0]), CompareAttrbt);

    double boundScale = 1.25;
    for (unsigned int i = 0; i != pAttrbtArray->numReadGrp; ++i)
    {
        pAttrbtArray->pBoundaries[i][0] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale;
        pAttrbtArray->pBoundaries[i][1] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale;
    }
}

void TGM_ReadPairMakeInverted(TGM_ReadPairAttrbtArray* pAttrbtArrays[2], const TGM_LibInfoTable* pLibTable, const TGM_LocalPairArray* pInvertedPairArray)
{
    TGM_ReadPairAttrbtArrayReInit(pAttrbtArrays[0], pInvertedPairArray->size / 2);
    TGM_ReadPairAttrbtArrayReInit(pAttrbtArrays[1], pInvertedPairArray->size / 2);

    pAttrbtArrays[0]->size = 0;
    pAttrbtArrays[1]->size = 0;

    pAttrbtArrays[0]->readPairType = PT_INVERTED3;
    pAttrbtArrays[1]->readPairType = PT_INVERTED5;

    pAttrbtArrays[0]->numReadGrp = pLibTable->size;
    pAttrbtArrays[1]->numReadGrp = pLibTable->size;

    for (unsigned int i = 0; i != pInvertedPairArray->size; ++i)
    {
        const TGM_LocalPair* pInvertedPair = pInvertedPairArray->data + i;
        int arrayIndex = pInvertedPair->readPairType - PT_INVERTED3;

        if (pAttrbtArrays[arrayIndex]->size == pAttrbtArrays[arrayIndex]->capacity)
            TGM_ARRAY_RESIZE(pAttrbtArrays[arrayIndex], pAttrbtArrays[arrayIndex]->capacity * 2, TGM_ReadPairAttrbt);

        TGM_ReadPairAttrbt* pAttrbt = pAttrbtArrays[arrayIndex]->data + pAttrbtArrays[arrayIndex]->size;
        ++(pAttrbtArrays[arrayIndex]->size);

        pAttrbt->origIndex = i;
        pAttrbt->readGrpID = pInvertedPair->readGrpID;

        double median = pLibTable->pLibInfo[pInvertedPair->readGrpID].fragLenMedian;

        pAttrbt->firstAttribute =  pInvertedPair->upPos + (double) pInvertedPair->fragLen / 2;
        pAttrbt->secondAttribute = pInvertedPair->fragLen - median;

        ++(pAttrbtArrays[arrayIndex]->size);
    }

    qsort(pAttrbtArrays[0]->data, pAttrbtArrays[0]->size, sizeof(pAttrbtArrays[0]->data[0]), CompareAttrbt);
    qsort(pAttrbtArrays[1]->data, pAttrbtArrays[1]->size, sizeof(pAttrbtArrays[1]->data[0]), CompareAttrbt);

    // neighbourhood scale factor (borrowed from Spanner)
    double boundScale[2] = {1.25, 0.5};
    for (unsigned int i = 0; i != pAttrbtArrays[0]->numReadGrp; ++i)
    {
        pAttrbtArrays[0]->pBoundaries[i][0] = (double) pLibTable->pLibInfo[i].fragLenMedian * boundScale[0] * 2.0;
        pAttrbtArrays[0]->pBoundaries[i][1] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale[1];
    }

    memcpy(pAttrbtArrays[1]->pBoundaries, pAttrbtArrays[0]->pBoundaries, sizeof(double) * 2 * pAttrbtArrays[0]->numReadGrp);
}

void TGM_ReadPairMakeSpecial(TGM_ReadPairAttrbtArray* pAttrbtArrays[], int numArray, const TGM_SpecialPairArray* pSpecialPairArray, const TGM_LibInfoTable* pLibTable)
{
    // initialize all the attribute arrays
    for (unsigned int i = 0; i != numArray; i += 2)
    {
        TGM_ReadPairAttrbtArrayReInit(pAttrbtArrays[i], DEFAULT_RP_ATTRB_CAPACITY);
        TGM_ReadPairAttrbtArrayReInit(pAttrbtArrays[i + 1], DEFAULT_RP_ATTRB_CAPACITY);

        pAttrbtArrays[i]->readPairType = PT_SPECIAL3;
        pAttrbtArrays[i + 1]->readPairType = PT_SPECIAL5;
    }

    // fill the arrays with special pair information
    for (unsigned int i = 0; i != pSpecialPairArray->size; ++i)
    {
        const TGM_SpecialPair* pSpecialPair = pSpecialPairArray->data + i;
        int arrayIndex = pSpecialPair->specialID * 2 + pSpecialPair->readPairType - PT_SPECIAL3;

        if (pAttrbtArrays[arrayIndex]->size == pAttrbtArrays[arrayIndex]->capacity)
            TGM_ARRAY_RESIZE(pAttrbtArrays[arrayIndex], pAttrbtArrays[arrayIndex]->capacity * 2, TGM_ReadPairAttrbt);

        TGM_ReadPairAttrbt* pAttrbt = pAttrbtArrays[arrayIndex]->data + pAttrbtArrays[arrayIndex]->size;
        ++(pAttrbtArrays[arrayIndex]->size);

        pAttrbt->origIndex = i;
        pAttrbt->readGrpID = pSpecialPair->readGrpID;

        double halfMedian = pLibTable->pLibInfo[pSpecialPair->readGrpID].fragLenMedian / 2.0;
        double halfMedians[2] = {halfMedian, -halfMedian};
        uint32_t pos[2] = {pSpecialPair->pos[0], pSpecialPair->end[0]};

        int posIndex = pSpecialPair->readPairType - PT_SPECIAL3;

        pAttrbt->firstAttribute = pos[posIndex] + halfMedians[posIndex];
        pAttrbt->secondAttribute = 0;
    }

    // sort the attribute arrays according to the first attribute
    for (unsigned int i = 0; i != numArray; i += 2)
    {
        qsort(pAttrbtArrays[i]->data, pAttrbtArrays[i]->size, sizeof(pAttrbtArrays[i]->data[0]), CompareAttrbt);
        qsort(pAttrbtArrays[i + 1]->data, pAttrbtArrays[i + 1]->size, sizeof(pAttrbtArrays[i + 1]->data[0]), CompareAttrbt);
    }

    // set the boundary
    double boundScale = 1.25;
    for (unsigned int i = 0; i != pAttrbtArrays[0]->numReadGrp; ++i)
    {
        pAttrbtArrays[0]->pBoundaries[i][0] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale;
        pAttrbtArrays[0]->pBoundaries[i][1] = 1e-3;
    }

    for (unsigned int i = 1; i != numArray; ++i)
        memcpy(pAttrbtArrays[i]->pBoundaries, pAttrbtArrays[0]->pBoundaries, sizeof(double) * 2 * pAttrbtArrays[0]->numReadGrp);
}

