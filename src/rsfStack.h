////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.4
////
////  Copyright 2013, Cleveland Clinic Foundation
////
////  This program is free software; you can redistribute it and/or
////  modify it under the terms of the GNU General Public License
////  as published by the Free Software Foundation; either version 2
////  of the License, or (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
////  You should have received a copy of the GNU General Public
////  License along with this program; if not, write to the Free
////  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
////  Boston, MA  02110-1301, USA.
////
////  Written by:
////    Hemant Ishwaran, Ph.D.
////    Director of Statistical Methodology
////    Professor, Division of Biostatistics
////    Clinical Research Building, Room 1058
////    1120 NW 14th Street
////    University of Miami, Miami FL 33136
////
////    email:  hemant.ishwaran@gmail.com
////    URL:    http://web.ccs.miami.edu/~hishwaran
////    --------------------------------------------------------------
////    Udaya B. Kogalur, Ph.D.
////    Adjunct Staff
////    Dept of Quantitative Health Sciences
////    Cleveland Clinic Foundation
////    
////    Kogalur & Company, Inc.
////    5425 Nestleway Drive, Suite L1
////    Clemmons, NC 27012
////
////    email:  commerce@kogalur.com
////    URL:    http://www.kogalur.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************

#ifndef RSFSTACK_H
#define RSFSTACK_H
#include "node.h"
void stackPreDefinedCommonArrays();
void unstackPreDefinedCommonArrays();
void stackPreDefinedGrowthArrays();
void unstackPreDefinedGrowthArrays();
void stackPreDefinedPredictArrays();
void unstackPreDefinedPredictArrays();
void stackPreDefinedInteractionArrays();
void unstackPreDefinedInteractionArrays();
void initializeArrays(char mode);
void initializeFactorArrays(char mode);
char stackCompetingArrays(char mode);
void getEventTypeSize(uint     obsSize, 
                      double  *status, 
                      uint    *mRecordMap, 
                      int    **mvSign,  
                      char     overWriteFlag,
                      uint    *eventTypeSize,
                      uint    *msize,
                      uint    *eventType);
void unstackCompetingArrays(char mode);
void stackFactorArrays(char mode);
void unstackFactorArrays(char mode);
char stackMissingSignatures (uint     obsSize, 
                             double  *statusPtr, 
                             double  *timePtr, 
                             double **predictorPtr,
                             uint    *recordMap,
                             uint     recordSize, 
                             uint   **p_recordIndex, 
                             uint    *p_vSize,
                             int   ***p_vSign, 
                             int    **p_vIndex,
                             int   ***p_vForestSign,
                             uint    *p_mFactorSize,
                             uint   **p_mFactorIndex);
void unstackMissingSignatures(uint      obsSize, 
                              double   *statusPtr, 
                              double   *timePtr, 
                              double  **predictorPtr,
                              uint     *recordMap,
                              uint      recordSize, 
                              uint     *recordIndex, 
                              uint      vSize,
                              int     **vSign, 
                              int      *vIndex,
                              int     **vForestSign,
                              uint      mFactorSize,
                              uint     *mFactorIndex);
char stackMissingArrays(char       mode,
                        char    ***p_dmRecordBootFlag,
                        double ****p_dmvImputation);
void unstackMissingArrays(char      mode,
                          char    **dmRecordBootFlag,
                          double ***dmvImputation);
uint stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***p_root,
                               double  **p_oobEnsemble,
                               double  **p_fullEnsemble,
                               double  **p_performance,
                               uint    **p_leafCount,
                               uint    **p_proximity,
                               double  **p_varImportance,
                               int     **p_seed,
                               double  **p_oobImputation,
                               double  **p_imputation,
                               double  **p_sumImputeStatusPtr,
                               double  **p_sumImputeTimePtr,
                               double ***p_sumImputePredictorPtr,
                               double  **p_sumOOBimputeStatusPtr,
                               double  **p_sumOOBimputeTimePtr,
                               double ***p_sumOOBimputePredictorPtr,
                               uint    **p_varUsed,
                               uint   ***p_varUsedPtr,
                               double  **p_splitDepth,
                               double ***p_splitDepthPtr,
                               double  **localSplitDepthPtr,
                               double  **p_oobPOEEnsemble,
                               double  **p_fullPOEEnsemble,
                               uint     *stackCount,
                               SEXP     *sexpVector);
void unstackDefinedOutputObjects(char      mode,
                                 Node    **root,
                                 double   *localSplitDepthPtr);
uint stackVariableOutputObjects(uint     totalNodeCount,
                                uint     totalMWCPCount,
                                uint   **p_treeID,
                                uint   **p_nodeID,
                                uint   **p_parmID,
                                double **p_contPT,
                                uint   **p_mwcpSZ,
                                uint   **p_mwcpPT,
                                uint     sexpLength,
                                char   **sexpString,
                                SEXP    *sexpVector);
#endif
