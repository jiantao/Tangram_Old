/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairScan.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/11/2012 02:44:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_READPAIRSCAN_H
#define  TGM_READPAIRSCAN_H

#include "TGM_LibInfo.h"

// parameters used for read pair build
typedef struct
{
    double cutoff;                 // fragment length cutoff p-value

    double trimRate;               // trim rate from fragment length distribution

    unsigned char minMQ;           // minimum mapping quality for a read pair

    FILE* fileListInput;           // input stream of a file list containing all the bam file names

    const char* workingDir;        // working directory for the detector

    const char* specialPrefix;     // prefix of the special reference

    uint32_t prefixLen;            // length of the prefix of the special reference

}TGM_ReadPairScanPars;

void TGM_ReadPairScan(const TGM_ReadPairScanPars* pScanPars);

void TGM_SpecialIDUpdate(TGM_SpecialID* pSpecialID, const TGM_ZAtag* pZAtag);

#endif  /*TGM_READPAIRSCAN_H*/
