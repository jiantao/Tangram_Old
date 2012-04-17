/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairBuild.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/08/2011 11:44:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_READPAIRBUILD_H
#define  TGM_READPAIRBUILD_H

#include "TGM_LibInfo.h"


//===============================
// Type and constant definition
//===============================

// parameters used for read pair build
typedef struct
{
    double cutoff;                 // fragment length cutoff p-value

    double trimRate;               // trim rate from fragment length distribution

    unsigned int binLen;           // bin length for searching mate

    uint32_t detectSet;            // bit set used to store the SV types that should be detected 

    unsigned char minMQ;           // minimum mapping quality for a read pair

    unsigned char minSpMQ;         // minimum mapping quality for the anchor of the special pair

    FILE* fileListInput;           // input stream of a file list containing all the bam file names

    FILE* splitFileInput;          // input stream of a file list containing all the split bam file names

    const char* workingDir;        // working directory for the detector

    const char* specialPrefix;     // prefix of the special reference

    uint32_t prefixLen;            // length of the prefix of the special reference

}TGM_ReadPairBuildPars;

// local pair structure(for deletion, tademn duplication and inversion)
typedef struct TGM_LocalPairInfo
{
    int32_t readGrpID;                                        // read group ID

    int32_t refID:16, upMapQ:8, downMapQ:8;                   // reference ID, up mate mapping quality, down mate mapping quality

    int32_t upPos;                                            // alignment position of the up mate

    int32_t upEnd;                                            // alignment end of the up mate

    int32_t downPos;                                          // alignment position of the down mate

    int32_t fragLen;                                          // fragment length of the read pair

    int32_t upNumMM:16, downNumMM:16;                         // number of mismatches of the up mate, number of mismatches of the down mate

    int32_t fragLenQual:16, readPairType:8, pairMode:8;       // fragment length quality(-10log(p-value)), read pairt type, pair mode

}TGM_LocalPair;

// cross pair structure(two mate hit different reference uniquely)
typedef struct TGM_CrossPair
{
    int32_t readGrpID;                                            // read group ID
                                                                                                                                          
    int32_t upRefID:16, downRefID:16;                             // up reference ID, down reference ID
                                                                                                                                          
    int32_t upMapQ:16, downMapQ:16;                               // up mapping quality, down mapping quality
                                                                  
    int32_t upPos;                                                // alignment position of the up mate
                                                                                                                                          
    int32_t upEnd;                                                // alignment end of the up mate
                                                                                                                                          
    int32_t downPos;                                              // alignment position of the down mate
                                                                                                                                          
    int32_t downEnd;                                              // alignment end of the down mate
                                                                                                                                          
    int32_t upNumMM:16, downNumMM:16;                             // number of mismatches of the up mate, number of mismatches of the down
                                                                                                                                          
    int32_t fragLenQual:16, readPairType:8, pairMode:8;           // fragment length quality(-10log(p-value)), read pairt type, pair mode

}TGM_CrossPair;

// special pair structure (one mate is unique the other mate hits the special reference)
typedef struct TGM_SpecialPair
{
    int32_t readGrpID;                                           // read group ID
                                                                                                                                         
    int16_t refID[2];                                            // up reference ID, down reference ID
                                                                                                                                         
    int32_t pos[2];                                              // alignment position of the up and down mate
                                                                 
    int32_t end[2];                                              // alignment end of the up and down mate

    int16_t numMM[2];                                            // number of mismatches of the up mate, number of mismatches of the down
                                                                                                                                         
    uint8_t bestMQ[2];                                           // best mapping quality of the up and down mate

    uint16_t numSpeicalHits;                                     // number of hits on the special reference

    int32_t specialID:16, readPairType:8, pairMode:8;            // fragment length quality(-10log(p-value)), read pairt type, pair mode

}TGM_SpecialPair;

// local pair array
typedef struct TGM_LocalPairArray
{
    TGM_LocalPair* data;     // local pairs

    uint64_t* chrCount;     // number of pairs in each chromosome

    uint64_t size;          // numer of local pairs

    uint64_t capacity;      // capacity of the array

}TGM_LocalPairArray;

// cross pair array
typedef struct TGM_CrossPairArray
{
    TGM_CrossPair* data;     // cross pairs

    uint64_t* chrCount;     // number of pairs in each chromosomes(up mate)

    uint64_t size;          // number of cross pairs

    uint64_t capacity;      // capacity of the array

}TGM_CrossPairArray;

// special pair array
typedef struct TGM_SpecialPairArray
{
    TGM_SpecialPair* data;    // special pairs

    uint64_t* chrCount;      // number of pairs in each chromosomes(anchor mate)

    uint64_t size;           // number of special pairs

    uint64_t capacity;       // capacity of the array

}TGM_SpecialPairArray;

// special pair table
typedef struct TGM_SpecialPairTable
{
    char (*names)[3];                    // array of special reference names(from ZA tags)

    void* nameHash;                      // special reference name hash

    uint32_t size;                       // number of the special reference names

    uint32_t capacity;                   // capacity of the name array

    TGM_SpecialPairArray array;           // array for local special pairs

    TGM_SpecialPairArray crossArray;      // array for cross special pairs

}TGM_SpecialPairTable;

// read pair table
typedef struct TGM_ReadPairTable
{
    TGM_LocalPairArray* pLongPairArray;          // long pair array(deletion)

    TGM_LocalPairArray* pShortPairArray;         // short pair array(tadem dup)

    TGM_LocalPairArray* pReversedPairArray;      // reveresed pair array (tadem dup)

    TGM_LocalPairArray* pInvertedPairArray;      // inverted pair array (inversion)

    TGM_SpecialPairTable* pSpecialPairTable;     // special pair table (retro insertion)

    TGM_CrossPairArray* pCrossPairArray;         // cross pair array (inter-chromosome translocation)

    uint32_t detectSet;                         // SV detect set

    uint32_t numChr;                            // number of references

    uint64_t numPairs;                          // total number of pairs in the read pair table

}TGM_ReadPairTable;

typedef struct TGM_SplitPair
{
    int32_t readGrpID;             // read group ID
                                                                                                           
    int16_t refID[2];              // up reference ID, down reference ID
                                                                                                           
    int32_t pos[2];                // alignment position of the up and down mate
                                   
    int32_t end[2];                // alignment end of the up and down mate

    int16_t numMM[2];              // number of mismatches of the up mate, number of mismatches of the down
                                                                                                           
    uint8_t bestMQ[2];             // best mapping quality of the up and down mate

    int16_t readPairType;          // fragment length quality(-10log(p-value)), read pairt type, pair mode

}TGM_SplitPair;

typedef struct TGM_SplitPairArray
{
    TGM_SplitPair* data;      // data array

    uint64_t size;           // number of split pairs stored in the array

    uint64_t capacity;       // capacity of the split pair array

    uint64_t* chrCount;      // number of pairs in each chromosomes(anchor mate)

}TGM_SplitPairArray;

typedef struct TGM_SplitPairTable
{
    TGM_SplitPairArray* pSplitLongArray;            // long split pair array(deletion)

    TGM_SplitPairArray* pSplitShortArray;           // short split pair array(medium-size insertion)

    TGM_SplitPairArray* pSplitReversedArray;        // reversed split pair array(tandem duplication)

    TGM_SplitPairArray* pSplitInvertedArray;        // inverted split pair array(inversion)

    TGM_SplitPairArray* pSplitSpecialArray;         // special split pair array(MEI)

    uint32_t detectSet;                            // detect set(indicates which SV we should detect)

    uint32_t numChr;                               // number of chromosomes in the bam

}TGM_SplitPairTable;


//===============================
// Constructors and Destructors
//===============================

TGM_SpecialPairTable* TGM_SpecialPairTableAlloc(unsigned int numChr);

void TGM_SpecialPairTableFree(TGM_SpecialPairTable* pSpecialPairTable);

TGM_ReadPairTable* TGM_ReadPairTableAlloc(uint32_t numChr, uint32_t detectSet);

void TGM_ReadPairTableFree(TGM_ReadPairTable* pInfoTable);

TGM_SplitPairTable* TGM_SplitPairTableAlloc(uint32_t numChr, uint32_t detectSet);

void TGM_SplitPairTableFree(TGM_SplitPairTable* pSplitPairTable);

//======================
// Interface functions
//======================

//===============================================================
// function:
//      read the bam alignments from the primary bam files, 
//      extract the SV-related information out and wirte it into 
//      read pair files according to the read pair type and 
//      location
//
// args:
//      1. pBuildPars: a pointer to the struct of parameters 
//      for read pair build
// 
//=============================================================== 
void TGM_ReadPairBuild(const TGM_ReadPairBuildPars* pBuildPars);

//===============================================================
// function:
//      clear the read pair table
//
// args:
//      1. pReadPairTable: a pointer to the read pair table
//                         structure
// 
//=============================================================== 
void TGM_ReadPairTableClear(TGM_ReadPairTable* pReadPairTable);

//======================================================================
// function:
//      update the read pair table with the incoming read pairs
//
// args:
//      1. pReadPairTable: a pointer to the read pair table structure
//      2. pUpAlgn: a pointor to the bam alignment of the up mate
//      3. pDownAlgn: a pointor to the bam alignment of the down mate
//      4. pZAtag: a pointer to the ZA tag structure
//      5. pPairStats: a pointer to the pair stats structure
//      6. pMateInfo: a pointer to the mate information structure
//      7. pLibTable: a pointer to the library information table
//      8. pHistArray: a pointer to the fragment length histogram array
//      9. minMQ: minimum mapping quality for a read pair
//======================================================================
void TGM_ReadPairTableUpdate(TGM_ReadPairTable* pReadPairTable, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const TGM_ZAtag* pZAtag, const TGM_PairStats* pPairStats, 
                            const TGM_MateInfo* pMateInfo, const TGM_LibInfoTable* pLibTable, const TGM_FragLenHistArray* pHistArray, const TGM_ReadPairBuildPars* pBuildPars);

//================================================================
// function:
//      open a seriers read pair files for output. A read pair
//      file will be open for each read pair type on each 
//      chromosome
//
// args:
//      1. pLibTable: a pointer to the library information table
//      2. detectSet: a bit set indicating which SV event shoud
//                    be detected
//      3. workingDir: the working directory for the detector
//
// return:
//      a pointer ot a hash table of read pair files
//================================================================
void* TGM_ReadPairFilesOpen(const TGM_LibInfoTable* pLibTable, uint32_t detectSet, const char* workingDir);

//==================================================================
// function:
//      close all the read pair files
//
// args:
//      1. pFileHash: a pointer to a hash table of read pair files
//==================================================================
void TGM_ReadPairFilesClose(void* pFileHash);

//=================================================================
// function:
//      write the read pair table into the read pair files
//
// args:
//      1. pReadPairTable: a pointer to a read pair table
//      2. pFileHash: a pointer to a hash table of read pair files
//=================================================================
void TGM_ReadPairTableWrite(const TGM_ReadPairTable* pReadPairTable, void* pFileHash);

//=================================================================
// function:
//      write the number of read pairs for each chromosome and 
//      SV types into the beginning of each read pair file
//
// args:
//      1. pFileHash: a pointer to a hash table of read pair files
//=================================================================
void TGM_ReadPairFilesWriteNum(void* pFileHash);

//=================================================================
// function:
//      write the detect set to the end of the library information
//      file
//
// args:
//      1. pReadPairTable: a pointer to a read pair table
//      2. libOutput: a file pointer to a library information file
//=================================================================
void TGM_ReadPairTableWriteDetectSet(const TGM_ReadPairBuildPars* pBuildPars, FILE* libOutput);

//=================================================================
// function:
//      write the special reference name into the end of the
//      library information file
//
// args:
//      1. pSpecialTable: a pointer to a special pair table
//      2. libOutput: a file pointer to a library information file
//=================================================================
void TGM_SpecialPairTableWriteID(const char* pSpecialID, FILE* libOutput);

//=================================================================================
// function:
//      update the split pair table with the incoming read pair
//
// args:
//      1. pReadPairTable: a pointer to the read pair table structure
//      2. pLibTable: a pointer to the library information table
//      3. pFirstPartial: a pointor to the bam alignment of the first partial
//      4. pSecondPartial: a pointor to the bam alignment of the second partial
//      5. pZAtag: a pointer to the ZA tag structure
//      6. minMQ: minimum mapping quality for a read pair
//=================================================================================
void TGM_SplitPairTableUpdate(TGM_SplitPairTable* pSplitPairTable, const TGM_LibInfoTable* pLibTable, const bam1_t* pFirstPartial, 
                             const bam1_t* pSecondPartial, const TGM_ZAtag* pZAtag, uint8_t minMQ);

//=================================================================
// function:
//      write the split pair table into the read pair files
//
// args:
//      1. pSplitPairTable: a pointer to a read pair table
//      2. pFileHash: a pointer to a hash table of read pair files
//=================================================================
void TGM_SplitPairTableWrite(const TGM_SplitPairTable* pSplitPairTable, void* pFileHash);

void TGM_SplitPairTableClear(TGM_SplitPairTable* pSplitPairTable);

void TGM_ReadPairFilesWriteSplitNum(void* pFileHash);


#endif  /*TGM_READPAIRBUILD_H*/
