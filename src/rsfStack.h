//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.2.2
//
//  Copyright 2008, Cleveland Clinic Foundation
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public
//  License along with this program; if not, write to the Free
//  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
//  Boston, MA  02110-1301, USA.
//
//  Project funded by:
//    National Institutes of Health, HL072771-01
//
//    Michael S Lauer, MD, FACC, FAHA
//    Cleveland Clinic Lerner College of Medicine of CWRU
//    9500 Euclid Avenue
//    Cleveland, OH 44195
//
//    email:  lauerm@ccf.org
//    phone:   216-444-6798
//
//  Written by:
//    Hemant Ishwaran, Ph.D.
//    Dept of Quantitative Health Sciences/Wb4
//    Cleveland Clinic Foundation
//    9500 Euclid Avenue
//    Cleveland, OH 44195
//
//    email:  hemant.ishwaran@gmail.com
//    phone:  216-444-9932
//    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
//    --------------------------------------------------------------
//    Udaya B. Kogalur, Ph.D.
//    Kogalur Shear Corporation
//    5425 Nestleway Drive, Suite L1
//    Clemmons, NC 27012
//
//    email:  ubk2101@columbia.edu
//    phone:  919-824-9825
//    URL:    www.kogalur-shear.com
//
//**********************************************************************
//**********************************************************************

#ifndef RSFSTACK_H
#define RSFSTACK_H
#include "extern.h"
void stackPreDefinedCommonArrays();
void unstackPreDefinedCommonArrays();
void stackPreDefinedGrowthArrays();
void unstackPreDefinedGrowthArrays();
void stackPreDefinedPredictArrays();
void unstackPreDefinedPredictArrays();
void stackPreDefinedInteractionArrays();
void unstackPreDefinedInteractionArrays();
void initializeArrays(uint mode,
                      uint *sortedTimeInterestsize);
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
                          int   ***p_vForestSign);
char stackMissingArrays(uint       mode,
                        char    ***p_dmRecordBootFlag,
                        double ****p_dmvImputation);
void unstackMissingArrays(uint      mode,
                          char    **dmRecordBootFlag,
                          double ***dmvImputation);
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
                              int     **vForestSign);
uint stackDefinedOutputObjects(uint      mode,
                               uint      sortedTimeInterestSize,
                               char    **sexpString,
                               Node   ***p_root,
                               double  **p_oobEnsemble,
                               double  **p_fullEnsemble,
                               double  **p_performance,
                               uint    **p_leafCount,
                               uint    **p_proximity,
                               double  **p_varImportance,
                               int     **p_seed,
                               double  **p_imputation,
                               double  **p_oobImputation,
                               double  **p_sumImputeStatusPtr,
                               double  **p_sumImputeTimePtr,
                               double ***p_sumImputePredictorPtr,
                               double  **p_sumOOBimputeStatusPtr,
                               double  **p_sumOOBimputeTimePtr,
                               double ***p_sumOOBimputePredictorPtr,
                               uint    **p_varUsed,
                               uint   ***p_varUsedPtr,
                               uint     *stackCount,
                               SEXP     *sexpVector);
void unstackDefinedOutputObjects(uint      mode,
                                 uint      sortedTimeInterestSize,
                                 Node    **root);
uint stackVariableOutputObjects(uint     totalNodeCount,
                                uint   **p_treeID,
                                uint   **p_nodeID,
                                uint   **p_parmID,
                                double **p_spltPT,
                                uint     sexpLength,
                                char   **sexpString,
                                SEXP    *sexpVector);
#endif
