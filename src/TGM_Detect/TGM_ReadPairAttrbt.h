/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairData.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/28/2012 10:37:08 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_READPAIRDATA_H
#define  TGM_READPAIRDATA_H

#include <stdint.h>

#include "TGM_LibInfo.h"
#include "TGM_ReadPairBuild.h"

typedef struct TGM_ReadPairAttrbt
{
    double firstAttribute;

    double secondAttribute;

    uint64_t origIndex;

    int32_t readGrpID;

}TGM_ReadPairAttrbt;

typedef struct TGM_ReadPairAttrbtArray
{
    TGM_ReadPairAttrbt* data;

    double (*pBoundaries)[2];

    uint64_t size;

    uint64_t capacity;

    uint32_t numReadGrp;

    SV_ReadPairType readPairType;

}TGM_ReadPairAttrbtArray;

TGM_ReadPairAttrbtArray* TGM_ReadPairAttrbtArrayAlloc(uint32_t numReadGrp);

void TGM_ReadPairAttrbtArrayFree(TGM_ReadPairAttrbtArray* pAttrbtArray);

#define TGM_ReadPairAttrbtArrayGetFirstBound(pAttrbtArray, i) ((pAttrbtArray)->pBoundaries[(pAttrbtArray)->data[(i)].readGrpID][0])

#define TGM_ReadPairAttrbtArrayGetSecondBound(pAttrbtArray, i) ((pAttrbtArray)->pBoundaries[(pAttrbtArray)->data[(i)].readGrpID][1])

void TGM_ReadPairAttrbtArrayReInit(TGM_ReadPairAttrbtArray* pAttrbtArray, uint64_t newCapacity);

void TGM_ReadPairMakeLocal(TGM_ReadPairAttrbtArray* pAttrbtArray, const TGM_LocalPairArray* pLocalPairArray, const TGM_LibInfoTable* pLibTable, SV_ReadPairType readpairType);

void TGM_ReadPairMakeCross(TGM_ReadPairAttrbtArray* pAttrbtArray, const TGM_LibInfoTable* pLibTable, const TGM_CrossPairArray* pCrossPairArray);

void TGM_ReadPairMakeInverted(TGM_ReadPairAttrbtArray* pAttrbtArrays[2], const TGM_LibInfoTable* pLibTable, const TGM_LocalPairArray* pInvertedPairArray);

void TGM_ReadPairMakeSpecial(TGM_ReadPairAttrbtArray* pAttrbtArrays[], int numArray, const TGM_SpecialPairArray* pSpecialArray, const TGM_LibInfoTable* pLibTable);


#endif  /*TGM_READPAIRDATA_H*/
