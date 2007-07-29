//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.0.0
//
//  Copyright 2007, Cleveland Clinic
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
void stackPreDefinedCommonArrays(uint **p_oobSampleSize);
void unstackPreDefinedCommonArrays(uint *oobSampleSize);
void stackPreDefinedGrowthArrays(double ***p_masterSplit,
                                 uint    **p_masterSplitSize,
                                 uint   ***p_masterSplitOrder);
void unstackPreDefinedGrowthArrays(double **masterSplit,
                                   uint    *masterSplitSize,
                                   uint   **masterSplitOrder);
void stackPreDefinedPredictArrays();
void unstackPreDefinedPredictArrays();
char stackAndInitializeMissingArrays(uint       mode,
                                     char      *p_mTimeIndexFlag,
                                     char    ***p_mRecordBootFlag,
                                     double ****p_mvImputation);
void unstackMissingArrays(uint      mode,
                          char    **mRecordBootFlag,
                          double ***mvImputation);
uint stackDefinedOutputObjects(uint      mode,
                               uint      sortedTimeInterestSize,
                               char    **sexpString,
                               Node   ***p_root,
                               double ***p_oobEnsemblePtr,
                               double ***p_fullEnsemblePtr,
                               double  **p_ensembleRun,
                               uint    **p_ensembleDen,
                               double  **p_oobEnsemble,
                               double  **p_fullEnsemble,
                               double  **p_performance,
                               uint    **p_leafCount,
                               uint    **p_proximity,
                               double  **p_varImportance,
                               double ***p_vimpEnsembleRun,
                               int     **p_seed,
                               double  **p_imputation,
                               double  **p_sImputeStatusPtr,
                               double  **p_sImputeTimePtr,
                               double ***p_sImputePredictorPtr,
                               uint     *stackCount,
                               SEXP     *sexpVector);
uint stackVariableOutputObjects(uint     totalNodeCount,
                                uint   **p_treeID,
                                uint   **p_nodeID,
                                uint   **p_parmID,
                                double **p_spltPT,
                                uint     sexpLength,
                                char   **sexpString,
                                SEXP    *sexpVector);
void unstackDefinedOutputObjects(uint      mode,
                                 Node    **root,
                                 uint      sortedTimeInterestSize,
                                 double  **oobEnsemblePtr,
                                 double  **fullEnsemblePtr,
                                 double   *ensembleRun,
                                 uint     *ensembleDen,
                                 double  **vimpEnsembleRun,
                                 double  **sImputePredictorPtr);
void initializeArrays(uint mode,
                      uint *sortedTimeInterestsize);
#endif
