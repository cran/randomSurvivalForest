//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 2.0.0
//
//  Copyright 2006, Cleveland Clinic
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

#ifndef RSFUTIL_H
#define RSFUTIL_H
#include "node.h"
extern uint     stackCount;
extern uint      _traceFlagDiagLevel;
extern uint      _traceFlagIterValue;
extern uint      _traceFlagToggler;
extern uint      _mup;
extern uint      _forestSize;
extern uint      _splitRule;
extern double   *_randomCovariateWeight;
extern uint      _observationSize;
extern double   *_time;
extern uint     *_status;
extern uint      _timeInterestSize;
extern double   *_timeInterest;
extern uint      _xSize;
extern double   *_xData;
extern double  **_observation;
extern Node    **_nodeMembership;
extern uint     *_bootMembershipIndex;
extern char     *_bootMembershipFlag;
extern uint      _fobservationSize;
extern double   *_ftime;
extern uint     *_fstatus;
extern double  **_fobservation;
extern Node    **_fnodeMembership;
extern double   *_fxData;
extern uint   *_masterTimeIndex;
extern double *_masterTime;
extern int    *_seedPtr;
uint updateTimeStamp(uint before);
uint getTraceFlag();
void updateTraceFlag(char reset);
void stackPreDefinedCommonArrays(uint **p_oobSampleSize);
void unstackPreDefinedCommonArrays(uint *oobSampleSize);
void stackPreDefinedGrowthArrays(double ***p_masterSplit,
                                 uint **p_masterSplitSize,
                                 uint ***p_masterSplitOrder);
void unstackPreDefinedGrowthArrays(double **masterSplit,
                                   uint *masterSplitSize,
                                   uint **masterSplitOrder);
void stackPreDefinedPredictArrays();
void unstackPreDefinedPredictArrays();
uint stackDefinedOutputObjects(uint     mode,
                               uint     sortedTimeInterestSize,
                               char   **sexpString,
                               Node  ***p_root,
                               double ***p_oobEnsemblePtr,
                               double ***p_fullEnsemblePtr,
                               double **p_ensembleRun,
                               uint   **p_ensembleDen,
                               double **p_oobEnsemble,
                               double **p_fullEnsemble,
                               double **p_performance,
                               uint   **p_leafCount,
                               uint   **p_proximity,
                               double **p_varImportance,
                               double ***p_vimpEnsembleRun,
                               int    **p_seed,
                               SEXP    *sexpVector);
uint stackVariableOutputObjects(uint     sexpLength,
                                uint     forestNodeCount,
                                char   **sexpString,
                                uint   **p_treeID,
                                uint   **p_nodeID,
                                uint   **p_parmID,
                                double **p_spltPT,
                                SEXP    *sexpVector);
void unstackDefinedOutputObjects(uint      mode,
                          Node    **root,
                          uint      sortedTimeInterestSize,
                          double  **oobEnsemblePtr,
                          double  **fullEnsemblePtr,
                          double   *ensembleRun,
                          uint     *ensembleDen,
                          double  **vimpEnsembleRun);
void initializeArrays(uint mode,
                      uint *masterTimeSize,
                      uint *sortedTimeInterestsize,
                      uint *masterDeathTimeSize);
void freeTree(Node *parent);
char restoreTree(uint    b,
                 Node   *parent,
                 uint   *leafCount,
                 uint   *offset,
                 uint   *treeID,
                 uint   *nodeID,
                 uint   *parmID,
                 double *spltPT);
void saveTree(uint    b,
              Node   *parent,
              uint   *offset,
              uint   *treeID,
              uint   *nodeID,
              uint   *parmID,
              double *spltPT);
Node* getMembership(Node *parent,
                    double **predictor,
                    uint  index);
void getCumulativeHazardEstimate(double **cumulativeHazard,
  Node   *parent,
  uint    sortedTimeInterestSize,
  uint    masterTimesize
);
double getConcordanceIndex(uint size, 
                           double *timePtr, 
                           uint *statusPtr, 
                           double *mortality,
                           uint *oobCount
			   );
double getPerformance(
  uint mode,
  uint sortedTimeInterestSize,
  uint leafCount,
  uint    masterTimeSize,
  double **oobEnsemblePtr,
  double **fullEnsemblePtr,
  double *ensembleRun,
  uint   *ensembleDen,
  uint    oobSampleSize,
  Node   *rootPtr,
  double **vimpEnsembleRun
);
#endif
