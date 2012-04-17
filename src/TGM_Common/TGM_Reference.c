/*
 * =====================================================================================
 *
 *       Filename:  TGM_TGM_Reference.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/03/2011 08:06:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "md5.h"
#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Reference.h"


//===============================
// Type and constant definition
//===============================

// maximum number of character will be load in a line from the fasta file
#define MAX_REF_LINE 1024

// number of characters can be held in a reference object
// this value assure that the largest chromosome in the human reference, chromsome 1, can be
// load into the object without any reallocation.
#define DEFAULT_REF_CAP 300000000

// initialize the hash table for reference name
KHASH_MAP_INIT_STR(refName, int32_t);


//===================
// Static methods
//===================

// process a line of reference sequence
static const char* ProcessRefLine(unsigned short* len, char* buff)
{
    char* iter = buff;
    const char* refLine = buff;
    *len = 0;

    while (*refLine == ' ' || *refLine == '\t')
    {
        if (*refLine == '\n' || *refLine == '\0')
        {
            refLine = NULL;
            return refLine;
        }

        ++iter;
        ++refLine;
    }

    while (*iter != '\n' && *iter != ' ' && *iter != '\t' && *iter != '\0')
    {
        *iter = toupper(*iter);
        ++iter;
        ++(*len);
    }

    return refLine;
}

// process the header line in the fasta file to get the ID for the next chromosome
static void TGM_RefHeaderSetName(TGM_RefHeader* pRefHeader, const char* buff)
{
    // the end of chromosome ID is either detected with a space character, a tab character, a new line character or a null character.

    // skip the '>' character
    const char* header = buff + 1;
    while (*header != ' ' && *header != '\t' 
           && *header != '\n' && *header != '\0')
    {
        ++header;
    }

    if (pRefHeader->numRefs == pRefHeader->capacity)
    {
        pRefHeader->capacity *= 2;

        pRefHeader->names = (char**) realloc(pRefHeader->names, pRefHeader->capacity * sizeof(char*));
        if (pRefHeader->names == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the reference name.\n");

        // we need to allocate one more character for the trailing null character from sprintf
        pRefHeader->md5s = (char*) realloc(pRefHeader->md5s, (pRefHeader->capacity * MD5_STR_LEN + 1) * sizeof(char));
        if (pRefHeader->md5s == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for MD5 strings in a reference header object.\n");

        if (pRefHeader->pSpecialRefInfo == NULL)
        {
            pRefHeader->refFilePos = (int64_t*) realloc(pRefHeader->refFilePos, sizeof(int64_t) * pRefHeader->capacity);
            if (pRefHeader->refFilePos == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the storage of reference file positions in the reference header object.\n");

            pRefHeader->htFilePos = (int64_t*) realloc(pRefHeader->htFilePos, sizeof(int64_t) * pRefHeader->capacity);
            if (pRefHeader->htFilePos == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the storage of hash table file positions in the reference header object.\n");
        }
        else if (pRefHeader->numRefs == pRefHeader->numSeqs)
        {
            pRefHeader->refFilePos = (int64_t*) realloc(pRefHeader->refFilePos, sizeof(int64_t) * (pRefHeader->numSeqs + 1));
            if (pRefHeader->refFilePos == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the storage of reference file positions in the reference header object.\n");

            pRefHeader->htFilePos = (int64_t*) realloc(pRefHeader->htFilePos, sizeof(int64_t) * (pRefHeader->numSeqs + 1));
            if (pRefHeader->htFilePos == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the storage of hash table file positions in the reference header object.\n");
        }
    }

    unsigned int headerLen = header - buff - 1;
    if (headerLen != 0)
    {
        pRefHeader->names[pRefHeader->numRefs] = (char*) malloc((headerLen + 1) * sizeof(char));
        if (pRefHeader->names[pRefHeader->numRefs] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the reference name.\n");

        pRefHeader->names[pRefHeader->numRefs][headerLen] = '\0';
        strncpy(pRefHeader->names[pRefHeader->numRefs], buff + 1, headerLen);
    }
    else
        pRefHeader->names[pRefHeader->numRefs] = NULL;
}

static void TGM_RefHeaderSetMd5(TGM_RefHeader* pRefHeader, char* sequence, uint32_t seqLen)
{
    unsigned char MD5[MD5_CHECKSUM_LEN];
    memset(MD5, 0, MD5_CHECKSUM_LEN);

    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char*) sequence, seqLen);
    MD5Final(MD5, &context);

    char* md5String = pRefHeader->md5s + MD5_STR_LEN * pRefHeader->numRefs;
    for (unsigned int i = 0; i != MD5_CHECKSUM_LEN; ++i)
    {
        sprintf(md5String, "%02X", MD5[i]);
        md5String += 2;
    }
}


//===============================
// Constructors and Destructors
//===============================

// create a new reference object
TGM_Reference* TGM_ReferenceAlloc(void)
{
    TGM_Reference* newRef = (TGM_Reference*) malloc(sizeof(TGM_Reference));
    if (newRef == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a reference object.\n");

    newRef->sequence = (char*) malloc(sizeof(char) * DEFAULT_REF_CAP);
    if (newRef->sequence == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of sequence in a reference object.\n");

    newRef->id = 0;
    newRef->seqLen = 0;
    newRef->seqCap = DEFAULT_REF_CAP;

    return newRef;
}

// free an existing reference object
void TGM_ReferenceFree(TGM_Reference* pRef)
{
    if (pRef != NULL)
    {
        free(pRef->sequence);
        free(pRef);
    }
}

TGM_RefHeader* TGM_RefHeaderAlloc(uint32_t refCapacity, uint32_t seqCapacity)
{
    if (refCapacity < seqCapacity)
        TGM_ErrQuit("ERROR: the number of references should not less than that of sequences.\n");

    TGM_RefHeader* pRefHeader = (TGM_RefHeader*) calloc(1, sizeof(TGM_RefHeader));
    if (pRefHeader == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a reference header object.\n");

    pRefHeader->names = (char**) calloc(refCapacity, sizeof(char*));
    if (pRefHeader->names == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for reference names in a reference header object.\n");

    pRefHeader->dict = (void*) kh_init(refName);

    // we need to allocate one more character for the trailing null character from sprintf
    pRefHeader->md5s = (char*) calloc((refCapacity * MD5_STR_LEN + 1), sizeof(char));
    if (pRefHeader->md5s == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for MD5 strings in a reference header object.\n");

    pRefHeader->refFilePos = (int64_t*) malloc(sizeof(int64_t) * seqCapacity);
    if (pRefHeader->refFilePos == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of reference file positions in the reference header object.\n");

    pRefHeader->htFilePos = (int64_t*) malloc(sizeof(int64_t) * seqCapacity);
    if (pRefHeader->htFilePos == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of hash table file positions in the reference header object.\n");

    pRefHeader->pSpecialRefInfo = NULL;

    pRefHeader->numRefs = 0;
    pRefHeader->numSeqs = 0;
    pRefHeader->capacity = refCapacity;

    return pRefHeader;
}

void TGM_RefHeaderFree(TGM_RefHeader* pRefHeader)
{
    if (pRefHeader != NULL)
    {
        kh_destroy(refName, pRefHeader->dict);

        if (pRefHeader->names != NULL)
        {
            for (unsigned int i = 0; i != pRefHeader->numRefs; ++i)
                free(pRefHeader->names[i]);

            free(pRefHeader->names);
        }

        free(pRefHeader->md5s);
        free(pRefHeader->refFilePos);
        free(pRefHeader->htFilePos);
        TGM_SpecialRefInfoFree(pRefHeader->pSpecialRefInfo);

        free(pRefHeader);
    }
}

TGM_SpecialRefInfo* TGM_SpecialRefInfoAlloc(uint32_t capacity)
{
    TGM_SpecialRefInfo* pSpecialRefInfo = (TGM_SpecialRefInfo*) malloc(sizeof(TGM_SpecialRefInfo));
    if (pSpecialRefInfo == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the special reference information object.\n");

    pSpecialRefInfo->numRefs = 0;
    pSpecialRefInfo->capacity = capacity;

    pSpecialRefInfo->endPos = (uint32_t*) malloc(sizeof(uint32_t) * capacity);
    if (pSpecialRefInfo->endPos == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of end indices of special references.\n");

    return pSpecialRefInfo;
}

void TGM_SpecialRefInfoFree(TGM_SpecialRefInfo* pSpecialRefInfo)
{
    if (pSpecialRefInfo != NULL)
    {
        free(pSpecialRefInfo->endPos);
        free(pSpecialRefInfo);
    }
}

TGM_RefView* TGM_RefViewAlloc(void)
{
    TGM_RefView* pRefView = (TGM_RefView*) calloc(1, sizeof(TGM_RefView));
    if (pRefView == NULL)
        TGM_ErrQuit("Not enough memory for a reference view object.\n");

    return pRefView;
}

void TGM_RefViewFree(TGM_RefView* pRefView)
{
    free(pRefView);
}

//==========================================
// Interface functions related with input
//==========================================

// read the reference header from the reference file
TGM_RefHeader* TGM_RefHeaderRead(int64_t* pRefHeaderPos, FILE* refInput)
{
    size_t readSize = 0;
    int64_t refPos = 0;

    readSize = fread(pRefHeaderPos, sizeof(int64_t), 1, refInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the reference header position from the reference file.\n");

    if ((refPos = ftello(refInput)) < 0)
        TGM_ErrQuit("ERROR: Cannot get the offset of current file.\n");

    if (fseeko(refInput, *pRefHeaderPos, SEEK_SET) != 0)
        TGM_ErrQuit("ERROR: Cannot seek in the reference file.\n");

    uint32_t numRefs = 0;
    readSize = fread(&(numRefs), sizeof(uint32_t), 1, refInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of chromosomes from the reference file.\n");

    if (numRefs == 0)
        TGM_ErrQuit("ERROR: There are no reference sequences stored in the file.\n");


    uint32_t numSpecialRefs = 0;
    readSize = fread(&(numSpecialRefs), sizeof(uint32_t), 1, refInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of special references from the reference file.\n");

    TGM_RefHeader* pRefHeader = NULL;
    if (numSpecialRefs != 0)
    {
        pRefHeader = TGM_RefHeaderAlloc(numRefs, numRefs - numSpecialRefs + 1);
        pRefHeader->pSpecialRefInfo = TGM_SpecialRefInfoAlloc(numSpecialRefs);

        pRefHeader->pSpecialRefInfo->numRefs = numSpecialRefs;
        pRefHeader->numSeqs = numRefs - numSpecialRefs + 1;

        readSize = fread(pRefHeader->pSpecialRefInfo->endPos, sizeof(uint32_t), numSpecialRefs, refInput);
        if (readSize != pRefHeader->pSpecialRefInfo->numRefs)
            TGM_ErrQuit("ERROR: Cannot read the end positions of special references from the reference file.\n");
    }
    else
    {
        pRefHeader = TGM_RefHeaderAlloc(numRefs, numRefs - numSpecialRefs + 1);
        pRefHeader->pSpecialRefInfo = NULL;

        pRefHeader->numSeqs = numRefs;
    }

    pRefHeader->numRefs = numRefs;

    int khRet = 0;
    khiter_t iter;
    khash_t(refName)* hash = pRefHeader->dict;
    for (unsigned int i = 0; i != pRefHeader->numRefs; ++i)
    {
        uint32_t nameLen = 0;
        readSize = fread(&nameLen, sizeof(uint32_t), 1, refInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the reference name length from the reference file.\n");

        pRefHeader->names[i] = (char*) malloc(sizeof(char) * (nameLen + 1));
        if (pRefHeader->names[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the reference name.\n");

        pRefHeader->names[i][nameLen] = '\0';
        readSize = fread(pRefHeader->names[i], sizeof(char), nameLen, refInput);
        if (readSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot read the reference name from the reference file.\n");

        iter = kh_put(refName, hash, pRefHeader->names[i], &khRet);
        kh_value(hash, iter) = i;
    }

    readSize = fread(pRefHeader->md5s, sizeof(char), MD5_STR_LEN * pRefHeader->numRefs, refInput);
    if (readSize != MD5_STR_LEN * pRefHeader->numRefs)
        TGM_ErrQuit("ERROR: Cannot read the md5 strings from the reference file.\n");

    readSize = fread(pRefHeader->refFilePos, sizeof(int64_t), pRefHeader->numSeqs, refInput);
    if (readSize != pRefHeader->numSeqs)
        TGM_ErrQuit("ERROR: Cannot read the offset of chromosomes from the reference file.\n");

    readSize = fread(pRefHeader->htFilePos, sizeof(int64_t), pRefHeader->numSeqs, refInput);
    if (readSize != pRefHeader->numSeqs)
        TGM_ErrQuit("ERROR: Cannot read the offset of hash table from the reference file.\n");

    if (fseeko(refInput, refPos, SEEK_SET) != 0)
        TGM_ErrQuit("ERROR: Cannot seek in the reference file.\n");

    return pRefHeader;
}

// get the reference ID given the reference name
int32_t TGM_RefHeaderGetRefID(const TGM_RefHeader* pRefHeader, const char* refName)
{
    khiter_t iter;
    khash_t(refName)* hash = (khash_t(refName)*) pRefHeader->dict;
    iter = kh_get(refName, hash, refName);

    return iter == kh_end(hash)? -1 : kh_value(hash, iter);
}

TGM_Status TGM_SpecialRefRead(TGM_Reference* pSpecialRef, const TGM_RefHeader* pRefHeader, FILE* refInput)
{
    if (pRefHeader->pSpecialRefInfo != NULL)
    {
        int64_t refPos = ftello(refInput);
        if (refPos < 0)
            TGM_ErrQuit("ERROR: Cannot get the offset of current file.\n");

        int64_t specialRefPos = pRefHeader->refFilePos[pRefHeader->numSeqs - 1];
        if (fseeko(refInput, specialRefPos, SEEK_SET) != 0)
            TGM_ErrQuit("ERROR: Cannot seek in the reference file.\n");

        TGM_ReferenceRead(pSpecialRef, refInput);

        if (fseeko(refInput, refPos, SEEK_SET) != 0)
            TGM_ErrQuit("ERROR: Cannot seek in the reference file.\n");

        return TGM_OK;
    }

    return TGM_ERR;
}

// jump to a certain chromosome given the reference ID
TGM_Status TGM_ReferenceJump(FILE* refInput, const TGM_RefHeader* pRefHeader, int32_t refID)
{
    if (refID < 0)
        return TGM_ERR;

    int32_t seqID = TGM_RefHeaderGetSeqID(pRefHeader, refID);
    int64_t jumpPos = pRefHeader->refFilePos[seqID];

    if (fseeko(refInput, jumpPos, SEEK_SET) != 0)
        return TGM_ERR;

    return TGM_OK;
}

// read the reference sequence from the input reference file 
void TGM_ReferenceRead(TGM_Reference* pRef, FILE* refInput)
{
    size_t readSize = 0;

    readSize = fread(&(pRef->id), sizeof(int32_t), 1, refInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the pRef input file due to an error.\n");

    readSize = fread(&(pRef->seqLen), sizeof(uint32_t), 1, refInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read chromosome length from the reference file.\n");

    if (pRef->seqLen > pRef->seqCap)
    {
        pRef->seqCap = pRef->seqLen;

        free(pRef->sequence);
        pRef->sequence = (char*) malloc(sizeof(char) * pRef->seqCap);
        if (pRef->sequence == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for reference sequence.\n");
    }

    readSize = fread(pRef->sequence, sizeof(char), pRef->seqLen, refInput);
    if (readSize != pRef->seqLen)
        TGM_ErrQuit("ERROR: Cannot read chromosome sequence from the reference file.\n");
}


TGM_Status TGM_GetRefFromSpecialPos(TGM_RefView* pRefView, int32_t* pRefID, uint32_t* pPos, const TGM_RefHeader* pRefHeader, const TGM_Reference* pSpecialRef, uint32_t specialPos)
{
    if (pRefHeader->pSpecialRefInfo != NULL)
    {
        if (specialPos <= pRefHeader->pSpecialRefInfo->endPos[0])
        {
            *pRefID = pRefHeader->numSeqs - 1;
            *pPos    = specialPos;

            pRefView->id = *pRefID;
            pRefView->sequence = pSpecialRef->sequence;
            pRefView->seqLen = pRefHeader->pSpecialRefInfo->endPos[0] + 1;

            return TGM_OK;
        }

        for (unsigned int i = 0; i != pRefHeader->pSpecialRefInfo->numRefs - 1; ++i)
        {
            unsigned int beginPos = pRefHeader->pSpecialRefInfo->endPos[i] + DEFAULT_PADDING_LEN + 1;
            if (specialPos >= beginPos && specialPos <= pRefHeader->pSpecialRefInfo->endPos[i + 1])
            {
                *pRefID = pRefHeader->numSeqs + i;
                *pPos    = specialPos - beginPos;

                pRefView->id = *pRefID;
                pRefView->sequence = pSpecialRef->sequence + beginPos;
                pRefView->seqLen = pRefHeader->pSpecialRefInfo->endPos[i + 1] - beginPos + 1;

                return TGM_OK;
            }
        }
    }

    return TGM_ERR;
}


//==========================================
// Interface functions related with output
//==========================================

// read the reference sequence in the fasta file line by line, one chromosome at each time
TGM_Status TGM_ReferenceLoad(TGM_Reference* pRef, TGM_RefHeader* pRefHeader, FILE* faInput)
{
    char buff[MAX_REF_LINE];

    while (fgets(buff, MAX_REF_LINE, faInput) != NULL && buff[0] != '>')
    {
        // actual length of the pRef line
        // excluding the starting and trailing space and the new line character
        unsigned short length = 0;
        const char* refLine = ProcessRefLine(&length, buff);

        // we skip the following processes if we get an empty line
        if (refLine == NULL)
            continue;

        // we shouldn't get this step if we are handling human pRef
        // the default capacity(300Mbp) should be enough even for the largest chromosome in human genome
        if (length + pRef->seqLen > pRef->seqCap)
        {
            pRef->seqCap *= 2;
            pRef->sequence = (char*) realloc(pRef->sequence, sizeof(char) * pRef->seqCap);
            if (pRef->sequence == NULL) 
                TGM_ErrQuit("ERROR: Not enough memory for the storage of sequence in the pRef object.\n");
        }

        // copy the current line into the pRef object
        char* cpAddr= pRef->sequence + pRef->seqLen;
        strncpy(cpAddr, refLine, length);
        pRef->seqLen += length;
    }

    if (pRef->seqLen > 0)
    {
        pRef->id = pRefHeader->numRefs;
        TGM_RefHeaderSetMd5(pRefHeader, pRef->sequence, pRef->seqLen);
        ++(pRefHeader->numRefs);
        ++(pRefHeader->numSeqs);
    }
    else
    {
        free(pRefHeader->names[pRefHeader->numRefs]);
        pRefHeader->names[pRefHeader->numRefs] = NULL;
    }

    if (buff[0] == '>')
        TGM_RefHeaderSetName(pRefHeader, buff);
    else if (feof(faInput))
        return TGM_EOF;
    else
        return TGM_ERR;

    return TGM_OK;
}

TGM_Status TGM_SpecialRefLoad(TGM_Reference* pRef, TGM_RefHeader* pRefHeader, FILE* faInput)
{
    char buff[MAX_REF_LINE];
    char padding[DEFAULT_PADDING_LEN];

    for (unsigned int i = 0; i != DEFAULT_PADDING_LEN; ++i)
        padding[i] = 'X';

    TGM_SpecialRefInfo* pSpecialRefInfo = pRefHeader->pSpecialRefInfo;
    pRef->id = pRefHeader->numRefs;

    unsigned int currRefLen = 0;
    while (fgets(buff, MAX_REF_LINE, faInput) != NULL)
    {
        if (buff[0] != '>')
        {
            // actual length of the pRef line
            // excluding the starting and trailing space and the new line character
            unsigned short length = 0;
            const char* refLine = ProcessRefLine(&length, buff);

            if (refLine == NULL)
                continue;

            if (length + pRef->seqLen > pRef->seqCap)
            {
                pRef->seqCap = (length + pRef->seqLen) * 2;
                pRef->sequence = (char*) realloc(pRef->sequence, sizeof(char) * pRef->seqCap);
                if (pRef->sequence == NULL) 
                    TGM_ErrQuit("ERROR: Not enough memory for the storage of sequence in the pRef object.\n");
            }
            
            // copy the current line into the pRef object
            char* cpAddr= pRef->sequence + pRef->seqLen;
            strncpy(cpAddr, refLine, length);
            pRef->seqLen += length;
            currRefLen += length;
        }
        else
        {
            if (currRefLen > 0)
            {
                if (DEFAULT_PADDING_LEN + pRef->seqLen > pRef->seqCap)
                {
                    pRef->seqCap = (pRef->seqLen + DEFAULT_PADDING_LEN) * 2;
                    pRef->sequence = (char*) realloc(pRef->sequence, sizeof(char) * pRef->seqCap);
                    if (pRef->sequence == NULL) 
                        TGM_ErrQuit("ERROR: Not enough memory for the storage of sequence in the reference object.\n");
                }

                if (pSpecialRefInfo->numRefs == pSpecialRefInfo->capacity)
                {
                    pSpecialRefInfo->capacity *= 2;
                    pSpecialRefInfo->endPos = (uint32_t*) realloc(pSpecialRefInfo->endPos, pSpecialRefInfo->capacity * sizeof(uint32_t));
                    if (pSpecialRefInfo->endPos == NULL)
                        TGM_ErrQuit("ERROR: Not enough memory for the storage of the end positions in the special reference object.\n");
                }

                pSpecialRefInfo->endPos[pSpecialRefInfo->numRefs] = pRef->seqLen - 1;

                uint32_t beginPos = TGM_SpecialRefGetBeginPos(pRefHeader, pSpecialRefInfo->numRefs);
                uint32_t specialRefLen = pRef->seqLen - beginPos;
                TGM_RefHeaderSetMd5(pRefHeader, pRef->sequence + beginPos, specialRefLen);

                char* cpAddr= pRef->sequence + pRef->seqLen;
                strncpy(cpAddr, padding, DEFAULT_PADDING_LEN);
                pRef->seqLen += DEFAULT_PADDING_LEN;

                ++(pSpecialRefInfo->numRefs);
                ++(pRefHeader->numRefs);
            }
            else if (pSpecialRefInfo->numRefs > 0)
            {
                free(pRefHeader->names[pRefHeader->numRefs]);
                pRefHeader->names[pRefHeader->numRefs] = NULL;
            }

            TGM_RefHeaderSetName(pRefHeader, buff);
            currRefLen = 0;

            if (pRefHeader->names[pRefHeader->numRefs] == NULL)
                TGM_ReferenceSkip(pRefHeader, faInput);
        }
    }

    if (currRefLen > 0)
    {
        pSpecialRefInfo->endPos[pSpecialRefInfo->numRefs] = pRef->seqLen - 1;

        uint32_t beginPos = TGM_SpecialRefGetBeginPos(pRefHeader, pSpecialRefInfo->numRefs);
        uint32_t specialRefLen = pRef->seqLen - beginPos;
        TGM_RefHeaderSetMd5(pRefHeader, pRef->sequence + beginPos, specialRefLen);

        ++(pSpecialRefInfo->numRefs);
        ++(pRefHeader->numRefs);
    }
    else if (pSpecialRefInfo->numRefs > 0)
    {
        free(pRefHeader->names[pRefHeader->numRefs]);
        pRefHeader->names[pRefHeader->numRefs] = NULL;
    }

    if (pSpecialRefInfo->numRefs == 0)
    {
        TGM_SpecialRefInfoFree(pRefHeader->pSpecialRefInfo);
        pRefHeader->pSpecialRefInfo = NULL;
    }
    else
        ++(pRefHeader->numSeqs);

    if (!feof(faInput))
        return TGM_ERR;
    else
        return TGM_OK;
}

// skip the reference sequence with unknown chromosome ID
TGM_Status TGM_ReferenceSkip(TGM_RefHeader* pRefHeader, FILE* faInput)
{
    char buff[MAX_REF_LINE];

    while (fgets(buff, MAX_REF_LINE, faInput) != NULL && buff[0] != '>')
        continue;

    if (buff[0] == '>')
        TGM_RefHeaderSetName(pRefHeader, buff);
    else if (feof(faInput))
        return TGM_EOF;
    else
        return TGM_ERR;

    return TGM_OK;
}

// leave enough space at the beginning of the output reference
// output file to store the reference header position
void TGM_ReferenceLeaveStart(FILE* refOutput)
{
    int64_t emptyOffset = 0;
    size_t writeSize = 0;

    writeSize = fwrite(&emptyOffset, sizeof(int64_t), 1, refOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the offset of the reference header into the reference file.\n");
}

// set the reference header position
void TGM_ReferenceSetStart(int64_t refHeaderPos, FILE* refOutput)
{
    if (fseeko(refOutput, 0, SEEK_SET) != 0)
        TGM_ErrSys("ERROR: Cannot seek in the reference output file\n");

    size_t writeSize = 0;

    writeSize = fwrite(&refHeaderPos, sizeof(int64_t), 1, refOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the offset of the reference header into the reference file.\n");
}

// write a reference sequence into the reference output file
int64_t TGM_ReferenceWrite(const TGM_Reference* pRef, FILE* refOutput)
{
    int64_t offset = ftello(refOutput);
    if (offset < 0)
        TGM_ErrQuit("ERROR: Cannot get the offset of current file.\n");

    size_t writeSize = 0;

    writeSize = fwrite(&(pRef->id), sizeof(int32_t), 1, refOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write chromosome ID into the reference file.\n");

    writeSize = fwrite(&(pRef->seqLen), sizeof(uint32_t), 1, refOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write chromosome length into the reference file.\n");

    writeSize = fwrite(pRef->sequence, sizeof(char), pRef->seqLen, refOutput);
    if (writeSize != pRef->seqLen)
        TGM_ErrQuit("ERROR: Cannot write chromosome sequence into the reference file.\n");

    fflush(refOutput);

    return offset;
}

// write the reference header into the output reference file
int64_t TGM_RefHeaderWrite(const TGM_RefHeader* pRefHeader, FILE* refOutput)
{
    int64_t offset = ftello(refOutput);
    if (offset < 0)
        TGM_ErrQuit("ERROR: Cannot get the offset of current file.\n");

    size_t writeSize = 0;

    writeSize = fwrite(&(pRefHeader->numRefs), sizeof(uint32_t), 1, refOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the total number of chromosomes into the reference file.\n");

    if (pRefHeader->pSpecialRefInfo != NULL)
    {
        writeSize = fwrite(&(pRefHeader->pSpecialRefInfo->numRefs), sizeof(uint32_t), 1, refOutput);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the total number of special references into the reference file.\n");

        writeSize = fwrite(pRefHeader->pSpecialRefInfo->endPos, sizeof(uint32_t), pRefHeader->pSpecialRefInfo->numRefs, refOutput);
        if (writeSize != pRefHeader->pSpecialRefInfo->numRefs)
            TGM_ErrQuit("ERROR: Cannot write the end positions of special references into the reference file.\n");
    }
    else
    {
        uint32_t zeroSpecialRef = 0;
        writeSize = fwrite(&zeroSpecialRef, sizeof(uint32_t), 1, refOutput);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the total number of special references into the reference file.\n");
    }

    for (unsigned int i = 0; i != pRefHeader->numRefs; ++i)
    {
        uint32_t nameLen = strlen(pRefHeader->names[i]);
        writeSize = fwrite(&nameLen, sizeof(uint32_t), 1, refOutput);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the length of reference name into the reference file.\n");

        writeSize = fwrite(pRefHeader->names[i], sizeof(char),  nameLen, refOutput);
        if (writeSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot write the reference name into the reference file.\n");
    }

    writeSize = fwrite(pRefHeader->md5s, sizeof(char), MD5_STR_LEN * pRefHeader->numRefs, refOutput);
    if (writeSize != MD5_STR_LEN * pRefHeader->numRefs)
        TGM_ErrQuit("ERROR: Cannot write the MD5 strings into the reference file.\n");


    writeSize = fwrite(pRefHeader->refFilePos, sizeof(int64_t), pRefHeader->numSeqs, refOutput);
    if (writeSize != pRefHeader->numSeqs)
        TGM_ErrQuit("ERROR: Cannot write the file offsets of the reference into the reference file.\n");

    writeSize = fwrite(pRefHeader->htFilePos, sizeof(int64_t), pRefHeader->numSeqs, refOutput);
    if (writeSize != pRefHeader->numSeqs)
        TGM_ErrQuit("ERROR: Cannot write the file offsets of the hash table into the reference file.\n");

    fflush(refOutput);

    return offset;
}

