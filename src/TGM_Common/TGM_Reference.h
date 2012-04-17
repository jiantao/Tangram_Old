/*
 * =====================================================================================
 *
 *       Filename:  TGM_Reference.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/03/2011 08:05:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_REFERENCE_H
#define  TGM_REFERENCE_H

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>

#include "TGM_Types.h"


//===============================
// Type and constant definition
//===============================

// the default number of chromosomes
#define DEFAULT_NUM_CHR 100

#define DEFAULT_NUM_SPECIAL_REF 30

// the default length of padding between two special references
#define DEFAULT_PADDING_LEN 300

// reset the reference object for next reading
#define TGM_ReferenceReset(pRef)               \
    do                                        \
    {                                         \
        (pRef)->seqLen = 0;                   \
                                              \
    }while(0) 

typedef struct TGM_SpecialRefInfo
{
    uint32_t* endPos;

    uint32_t numRefs;

    uint32_t capacity;

}TGM_SpecialRefInfo;

// reference header strcture
typedef struct TGM_RefHeader
{
    char** names;             // an array contains the name of each chromosome

    void* dict;               // a hash table of reference name and reference ID

    char* md5s;               // an array contains the md5 string of each chromosome

    int64_t* refFilePos;      // an array contains the file offset poisition of each chromosomes in the reference file

    int64_t* htFilePos;       // an array contains the file offset poisition of each chromosomes in the hash table file

    uint32_t numRefs;         // total number of chromosomes 

    uint32_t numSeqs;         // total number of reference sequences (special reference sequence may contain one or more chromosomes)

    uint32_t capacity;

    TGM_SpecialRefInfo* pSpecialRefInfo;

}TGM_RefHeader;

// an object holds the reference sequence of a chromosome
typedef struct TGM_Reference
{
    char* sequence;               // reference sequence

    int32_t  id;                  // id of the sequence

    uint32_t seqLen;              // length of the chromosome

    uint32_t seqCap;              // capacity of reference sequence

}TGM_Reference;

// an object holds the pointer to an existed reference
typedef struct TGM_RefView
{
    const char* sequence;       // a pointer to an existed reference

    int32_t id;                 // the id of that sequence

    uint32_t seqLen;            // the length of that sequence

}TGM_RefView;

// get the reference name given the reference ID
#define TGM_RefHeaderGetName(pRefHeader, refID) ((pRefHeader)->names[(refID)])

// get the md5 string (not null terminated) given the reference ID
#define TGM_RefHeaderGetMD5(pRefHeader, refID) ((pRefHeader)->md5s + (refID) * MD5_STR_LEN)

//===============================
// Constructors and Destructors
//===============================

// create a new reference object
TGM_Reference* TGM_ReferenceAlloc(void);

// free an existing reference object
void TGM_ReferenceFree(TGM_Reference* pRef);

TGM_RefHeader* TGM_RefHeaderAlloc(uint32_t refCapacity, uint32_t seqCapacity);

void TGM_RefHeaderFree(TGM_RefHeader* pRefHeader);

TGM_SpecialRefInfo* TGM_SpecialRefInfoAlloc(uint32_t capcity);

void TGM_SpecialRefInfoFree(TGM_SpecialRefInfo* pSpecialRefInfo);

TGM_RefView* TGM_RefViewAlloc(void);

void TGM_RefViewFree(TGM_RefView* pRefView);

//==========================================
// Interface functions related with input
//==========================================

//====================================================================
// function:
//      read the reference header from the reference file
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. refInput: a file pointer to the input reference file
// 
// return:
//      the file offset of the reference header structure in the 
//      reference file. It is used to check the compatibility between
//      the reference file and the hash table file
//====================================================================
TGM_RefHeader* TGM_RefHeaderRead(int64_t* pRefStart, FILE* refInput);

//===================================================================
// function:
//      get the reference ID given the reference name
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. refName: a reference name
// 
// return:
//      if the reference name is found, return the corresponding
//      reference ID. Otherwise, return -1
//===================================================================
int32_t TGM_RefHeaderGetRefID(const TGM_RefHeader* pRefHeader, const char* refName);

static inline int32_t TGM_RefHeaderGetSeqID(const TGM_RefHeader* pRefHeader, int32_t refID)
{
    if (refID < 0)
        return refID;
    else
        return (refID <= (int32_t) pRefHeader->numSeqs - 1 ? refID : pRefHeader->numSeqs - 1);
}

//====================================================================
// function:
//      read the spcail reference sequence from the input 
//      reference file 
//
// args:
//      1. pSpecialRef: a pointer to the reference sequence structure
//      2. pRefHeader:  a pointer to the reference header structure
//      3. refInput: a file pointer to the input reference file
//
//return:
//      if there is no special reference sequence, return TGM_ERR,
//      otherwise return TGM_OK.
//====================================================================
TGM_Status TGM_SpecialRefRead(TGM_Reference* pSpecialRef, const TGM_RefHeader* pRefHeader, FILE* refInput);

//====================================================================
// function:
//      jump to a certain chromosome given the reference ID
//
// args:
//      1. refInput: a file pointer to the input reference file
//      2. pRefHeader: a pointer to the reference header structure
//      3. refID: the ID of the reference we want to jump to
// 
// return:
//      TGM_OK: successfully jumped
//      TGM_ERR: jump failed
//====================================================================
TGM_Status TGM_ReferenceJump(FILE* refInput, const TGM_RefHeader* pRefHeader, int32_t refID);

//====================================================================
// function:
//      read the reference sequence from the input reference file 
//
// args:
//      1. pRef: a pointer to the reference sequence structure
//      2. refInput: a file pointer to the input reference file
// 
//====================================================================
void TGM_ReferenceRead(TGM_Reference* pRef, FILE* refInput);

//=====================================================================
// function:
//      get the reference ID and the real position from the position 
//      of the special sequence
//
// args:
//      1. pRefID: a pointer to the reference ID
//      2. pPos: a pointer to the real reference position
//      3. pRefHeader: a pointer to the reference header structure
//      4. specialPos: position on the special sequence
// 
// return:
//      return a pointer to the reference view on success, 
//      otherwise return NULL
//=====================================================================
TGM_Status TGM_GetRefFromSpecialPos(TGM_RefView* pRefView, int32_t* pRefID, uint32_t* pPos, const TGM_RefHeader* pRefHeader, const TGM_Reference* pSpecialRef, uint32_t specialPos);

//==========================================
// Interface functions related with output
//==========================================

//===================================================================
// function:
//      read the reference sequence in the fasta file line by line, 
//      one chromosome at each time
//
// args:
//      1. pRef: a pointer to the reference structure
//      2. pRefHeader: a pointer to the reference header structure
//      3. faInput: a file pointer to the input fasta file
// 
// return:
//      TGM_OK: successfully load the chromosome sequence
//      TGM_EOF: reach the end of the fasta file
//      TGM_ERR: find an error during loading
//===================================================================
TGM_Status TGM_ReferenceLoad(TGM_Reference* pRef, TGM_RefHeader* pRefHeader, FILE* faInput);

//=====================================================================
// function:
//      get the special reference ID from the reference ID
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. refID: reference ID
// 
// return:
//      return the spcail reference ID on success, otherwise return -1
//=====================================================================
static inline int32_t TGM_GetSpecialRefIDFromRefID(const TGM_RefHeader* pRefHeader, int32_t refID)
{
    if (pRefHeader->pSpecialRefInfo != NULL && refID >= (int32_t) pRefHeader->numSeqs - 1)
        return (refID - (pRefHeader->numSeqs - 1));

    return -1;
}

//===================================================================
// function:
//      Get the begin position of a certain special chromosome in
//      the special reference sequence
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. specialRefID: the special reference ID (start from zero)
// 
// return:
//      the start position of a certain special chromosome
//===================================================================
static inline uint32_t TGM_SpecialRefGetBeginPos(const TGM_RefHeader* pRefHeader, int32_t specialRefID)
{
    return (specialRefID == 0 ? 0 : pRefHeader->pSpecialRefInfo->endPos[specialRefID - 1] + DEFAULT_PADDING_LEN + 1);
}

//===================================================================
// function:
//      load all the special chromosomes at once and then combine
//      them into a long special reference sequence separated by
//      "XXXX..." padding
//
// args:
//      1. pRef: a pointer to the reference structure
//      2. pRefHeader: a pointer to the reference header structure
//      3. faInput: a file pointer to the input fasta file
// 
// return:
//      TGM_OK: successfully load the special chromosome sequence
//      TGM_ERR: find an error during loading
//===================================================================
TGM_Status TGM_SpecialRefLoad(TGM_Reference* pRef, TGM_RefHeader* pRefHeader, FILE* faInput);

//===================================================================
// function:
//      skip the reference sequence with unknown chromosome
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. faInput: a file pointer to the input fasta file
// 
// return:
//      TGM_OK: successfully skipped a chromosome sequence
//      TGM_EOF: reach the end of the fasta file
//      TGM_ERR: find an error during skipping
//===================================================================
TGM_Status TGM_ReferenceSkip(TGM_RefHeader* pRefHeader, FILE* faInput);

//===================================================================
// function:
//      leave enough space at the beginning of the output reference
//      output file to store the reference header position
//
// args:
//      1. refOutput: a file pointer to the output reference file
//===================================================================
void TGM_ReferenceLeaveStart(FILE* refOutput);

//===================================================================
// function:
//      set the reference header position
//
// args:
//      1. refHeaderPos: the reference header position
//      2. refOutput: a file pointer to the output reference file
//===================================================================
void TGM_ReferenceSetStart(int64_t refHeaderPos, FILE* refOutput);

//===================================================================
// function:
//      write a reference sequence into the reference output file
//
// args:
//      1. pRef: a pointer to the reference sequence structure
//      2. refOutput: a file pointer to the output reference file
//
// return:
//      the file offset of the current reference sequence
//===================================================================
int64_t TGM_ReferenceWrite(const TGM_Reference* pRef, FILE* refOutput);

//====================================================================
// function:
//      write the reference header into the output reference file
//
// args:
//      1. pRefHeader: a pointer to the reference sequence structure
//      2. refOutput: a file pointer to the output reference file
//
// return:
//      the file offset of the reference header
//====================================================================
int64_t TGM_RefHeaderWrite(const TGM_RefHeader* pRefHeader, FILE* refOutput);


#endif  /*TGM_REFERENCE_H*/
