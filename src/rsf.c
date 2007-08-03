//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.0.1
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

#include    "global.h"
#include    "nrutil.h"
#include  "node_ops.h"
#include "rsfImpute.h"
#include  "rsfSplit.h"
#include   "rsfUtil.h"
#include  "rsfStack.h"
#include       "rsf.h"
uint stackCount;
char *sexpString[RSF_SEXP_CNT] = {
  "",              
  "",              
  "fullEnsemble",  
  "oobEnsemble",   
  "performance",   
  "leafCount",     
  "proximity",     
  "importance",    
  "treeID",        
  "nodeID",        
  "parmID",        
  "spltPT",        
  "seed",          
  "imputation"     
};
SEXP sexpVector[RSF_SEXP_CNT];
uint     *_treeID_;
uint     *_nodeID_;
uint     *_parmID_;
double   *_spltPT_;
int      *_seed_;
double   *_fullEnsemble_;
double   *_oobEnsemble_;
double   *_performance_;
uint     *_leafCount_;
uint     *_proximity_;
double   *_varImportance_;
double   *_imputation_;
uint      _mup;
uint      _splitRule;
uint      _forestSize;
uint      _minimumDeathCount;
uint      _randomCovariateCount;
double   *_randomCovariateWeight;
uint      _observationSize;
uint      _xSize;
double   *_time;
double   *_status;
double   *_xData;
uint      _fobservationSize;
double   *_ftime;
double   *_fstatus;
double   *_fxData;
uint      _timeInterestSize;
double   *_timeInterest;
SEXP      _sexp_xType;
char     **_xType;
double  **_observation;
double  **_fobservation;
double   *_masterTime;
uint     *_masterTimeIndex;
uint      _masterTimeSize;
uint     *_mRecordMap;
uint     *_fmRecordMap;
uint      _mRecordSize;
uint      _fmRecordSize;
uint     *_mRecordIndex;
uint     *_fmRecordIndex;
uint      _mvSize;
uint      _fmvSize;
int     **_mvSign;
int     **_fmvSign;
int      *_mvIndex;
int      *_fmvIndex;
int     **_mvForestSign;
int     **_fmvForestSign;
double   *_mStatus;
double   *_fmStatus;
double   *_mTime;
double   *_fmTime;
int      *_seed1Ptr;
int      *_seed2Ptr;
Node    **_nodeMembership;
uint     *_bootMembershipIndex;
char     *_bootMembershipFlag;
Node    **_fnodeMembership;
uint      _traceFlagDiagLevel;
uint      _traceFlagIterValue;
uint      _traceFlagToggler;
uint getTraceFlag() {
  uint result;
  if (_traceFlagToggler == 1) {
    result = _traceFlagDiagLevel;
  }
  else {
    result = 0;
  }
  return result;
}
void updateTraceFlag(char reset) {
  if (reset == TRUE) {
    _traceFlagToggler = 1;
  }
  else {
    if (_traceFlagToggler == 1) {
      _traceFlagToggler = _traceFlagIterValue;
    }
    else {
      _traceFlagToggler --;
    }
  }
}
char getBestSplit(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  double **masterSplit,
  uint   **masterSplitOrder) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetBestSplit() ENTRY ...\n");
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nAttempting to split node:  %10d", parent -> leafCount);
  }
  if (_splitRule == LOG_RANK) {
    return logRank(parent,
                   splitParameterMax,
                   splitValueMax,
                   masterSplit);
  }
  else if (_splitRule == CONSERVE_EVENTS) {
    return conserveEvents(parent,
                          splitParameterMax,
                          splitValueMax,
                          masterSplit);
  }
  else if (_splitRule == LOG_RANK_SCORE) {
    return logRankScore(parent,
                        splitParameterMax,
                        splitValueMax,
                        masterSplit,
                        masterSplitOrder);
  }
  else if (_splitRule == LOG_RANK_APPROX) {
    return logRankApprox(parent,
                          splitParameterMax,
                          splitValueMax,
                          masterSplit);
  }
  else {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid split rule:  %10d", _splitRule);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetBestSplit() EXIT ...\n");
  }
}
char forkAndUpdate(
  uint  *leafCount,
  Node  *parent,
  uint   splitParameter,
  uint   splitValueIndex,
  double splitValue) {
  uint i;
  uint leftDeathCount, rightDeathCount;
  char result;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nforkAndUpdate() ENTRY ...\n");
  }
  result = forkNode(parent, splitParameter, splitValueIndex, splitValue);
  if (result == TRUE) {
    leftDeathCount = rightDeathCount = 0;
    (*leafCount)++;
    for (i = 1; i <= _observationSize; i++) {
      if (_nodeMembership[i] == parent) {
        if (_observation[splitParameter][i] <= splitValue) {
          _nodeMembership[i] = parent -> left;
          ((parent -> left) -> leafCount) = (parent -> leafCount);
        }
        else {
          _nodeMembership[i] = parent -> right;
          ((parent -> right) -> leafCount) = *leafCount;
        }
      }
    }
    for (i = 1; i <= _observationSize; i++) {
      if (_nodeMembership[_bootMembershipIndex[i]] == parent -> left) {
        if (_status[_bootMembershipIndex[i]] == 1) leftDeathCount ++;
      }
      else if (_nodeMembership[_bootMembershipIndex[i]] == parent -> right) {
        if (_status[_bootMembershipIndex[i]] == 1) rightDeathCount ++;
      }
    }
    if (leftDeathCount  < (2 * (_minimumDeathCount))) {
      (parent -> left) -> splitFlag = FALSE;
    }
    if (rightDeathCount <= (2 * (_minimumDeathCount))) {
      (parent -> right)-> splitFlag = FALSE;
    }
  }
  else {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  forkNode() failed.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nforkAndUpdate() EXIT ...\n");
  }
  return result;
}
char makeTree (uint     b,
               Node    *parent,
               double **masterSplit,
               uint   **masterSplitOrder,
               char     mTimeIndexFlag) {
  char splitResult, forkResult;
  uint splitParameterMax, splitValueMax;
  uint i;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nmakeTree() ENTRY ...\n");
  }
  if (_mRecordSize > 0) {
    imputeNode(RSF_GROW, b, parent);
    if (mTimeIndexFlag == TRUE) {
      updateTimeIndexArray(parent);
    }
  }
  splitResult = getBestSplit(parent,
                             & splitParameterMax,
                             & splitValueMax,
                             masterSplit,
                             masterSplitOrder);
  if (splitResult == TRUE) {
    forkResult = forkAndUpdate(_leafCount_ + b,
                               parent,
                               splitParameterMax,
                               splitValueMax,
                               masterSplit[splitParameterMax][splitValueMax]);
    if (forkResult == TRUE) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nNode Membership:  \n");
        for (i=1; i <=  _observationSize; i++) {
          Rprintf("%10d %10d \n", i, _nodeMembership[i] -> leafCount);
        }
      }
      if ((parent -> left) -> splitFlag == TRUE) {
        if (getTraceFlag() & DL1_TRACE) {
          Rprintf("\nmakeTree() LEFT:  \n");
        }
        makeTree (b,
                  parent -> left,
                  masterSplit,
                  masterSplitOrder,
                  mTimeIndexFlag);
      }
      else {
        if (_mRecordSize > 0) {
          imputeNode(RSF_GROW, b, parent -> left);
          if (mTimeIndexFlag == TRUE) {
            updateTimeIndexArray(parent -> left);
          }
        }
      }
      if ((parent -> right) -> splitFlag == TRUE) {
        if (getTraceFlag() & DL1_TRACE) {
          Rprintf("\nmakeTree() RIGHT:  \n");
        }
        makeTree (b,
                  parent -> right,
                  masterSplit,
                  masterSplitOrder,
                  mTimeIndexFlag);
      }
      else {
        if (_mRecordSize > 0) {
          imputeNode(RSF_GROW, b, parent -> right);
          if (mTimeIndexFlag == TRUE) {
            updateTimeIndexArray(parent -> right);
          }
        }
      }
    }
    else {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  forkAndUpdate() failed.");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
    }
  }
  else {
    parent -> splitFlag = FALSE;
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\ngetBestSplit() FAILED ...\n");
    }
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nmakeTree() EXIT ...\n");
  }
  return splitResult;
}
char bootstrap (uint    mode,
                uint    b,
                Node    *rootPtr,
                uint    *oobSampleSize,
                double **masterSplit,
                uint    *masterSplitSize,
                uint   **masterSplitOrder,
                char   **mRecordBootFlag) {
  char mOutcomeFlag;
  uint **masterSplitBounds;
  uint *leadIndex;
  uint ibSampleSize;
  char result;
  uint i,j,k;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nbootstrap() ENTRY ...\n");
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Bootstrap sample:  %10d ", b);  
  }
  result = TRUE;
  mOutcomeFlag = FALSE;
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nStart of seed chain for bootstrap:  ");
    Rprintf("\n%10d %20d ", b, *_seed1Ptr);
  }
  for (i=1; i <= _observationSize; i++) {
    _nodeMembership[i] = rootPtr;
    _bootMembershipFlag[i] = FALSE;
  }
  for (i=1; i <= _observationSize; i++) {
    k = ceil(ran1(_seed1Ptr)*((_observationSize)*1.0));
    _bootMembershipFlag[k] = TRUE;
    _bootMembershipIndex[i] = k;
  }
  if (mode == RSF_PRED) {
    for (i=1; i <= _fobservationSize; i++) {
      _fnodeMembership[i] = rootPtr;
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nIn-Bag Membership:  ");
    Rprintf("\n     index    IBindex\n");
    for (i=1; i <=  _observationSize; i++) {
      Rprintf("%10d %10d \n", i, _bootMembershipIndex[i]);
    } 
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Bootstrapping and initialization of root node complete.");  
  }
  if (mode == RSF_GROW) {
    if (_mRecordSize > 0) {
      for (i = 1; i <= _mRecordSize; i++) {
        if (_bootMembershipFlag[_mRecordIndex[i]] == TRUE) {
          mRecordBootFlag[b][i] = TRUE;
        }
        else {
          mRecordBootFlag[b][i] = FALSE;
        }
      }
    }  
  }  
  for (i=1; i <= _observationSize; i++) {
    if (_bootMembershipFlag[i] == FALSE) {
      oobSampleSize[b] ++;
    }
  }
  ibSampleSize = _observationSize - oobSampleSize[b];
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nOOB Size:  %10d ", oobSampleSize[b]);
  }
  if (_mRecordSize > 0) {
    result = getForestSign(mode, b);
    if (result == FALSE) {
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Status or Time values are all missing in the sample.  Bootstrap sample has been discarded.");  
      }
    }
  }
  if ((result == TRUE) && (mode == RSF_GROW)) {
    leadIndex = uivector(1, _xSize);
    masterSplitBounds = uimatrix(1, _xSize, 1, 2);
    for (j=1; j <= _xSize; j++) {
      masterSplitBounds[j][1] = 0;
      masterSplitBounds[j][2] = 0;
      masterSplitSize[j] = 0;
      for (i = 1; i <= _observationSize; i++) {
        masterSplit[j][i] = 0;
      }
    }
    for (j=1; j <= _xSize; j++) {
        leadIndex[j] = 0;
        for (i=1; i <= _observationSize; i++) {
          if (_bootMembershipFlag[i] == TRUE) {
            if (_mRecordMap[i] == 0) {
              leadIndex[j]++;
              masterSplit[j][leadIndex[j]] = _observation[j][i];
            }
            else {
              if(_mvSign[j+2][_mRecordMap[i]] == 0) {
                leadIndex[j]++;
                masterSplit[j][leadIndex[j]] = _observation[j][i];
              }
            }
          }
        }
    }  
    if (getTraceFlag() & DL3_TRACE) {
      Rprintf("\n\nIn-Bag Raw MasterSplit Predcitor Data:  ");
      Rprintf("\n          ");
      for (j=1; j <= _xSize; j++) {
        Rprintf(" %12d", j);
      }
      for (i=1; i <= ibSampleSize; i++) {
        Rprintf("\n%10d", i);
        for (j=1; j <= _xSize; j++) {
          Rprintf(" %12.4f", masterSplit[j][i]);
        }
      }
    }
    for (j=1; j <= _xSize; j++) {
      if (leadIndex[j] > 0) {
        masterSplitBounds[j][1] = 1;
        masterSplitSize[j] = 1;
        hpsort(masterSplit[j], leadIndex[j]);
        for (i = 2; i <= leadIndex[j]; i++) {
          if (masterSplit[j][i] > masterSplit[j][masterSplitSize[j]]) {
            masterSplitSize[j]++;
            masterSplit[j][masterSplitSize[j]] = masterSplit[j][i];
          }
        }
        masterSplitBounds[j][2] = masterSplitSize[j];
        for (i = masterSplitSize[j]+1; i <= leadIndex[j]; i++) {
          masterSplit[j][i] = 0;
        }
      }
    }  
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\n\nIn-Bag Sorted MasterSplit Data:  ");
      Rprintf("\n          ");
      for (j=1; j <= _xSize; j++) {
        Rprintf(" %12d", j);
      }
      for (i=1; i <= ibSampleSize; i++) {
        Rprintf("\n%10d", i);
        for (j=1; j <= _xSize; j++) {
          Rprintf(" %12.4f", masterSplit[j][i]);
        }
      }
      Rprintf("\n\nSize of Permissible Splits:  ");
      Rprintf("\n          ");
      for (j=1; j <= _xSize; j++) {
        Rprintf(" %12d", masterSplitSize[j]);
      }
      Rprintf("\n");
      Rprintf("\nMaster Permissible Split Boundaries:  \n");
      for (i=1; i <= _xSize; i++) {
        Rprintf("%10d %10d %10d \n", i, masterSplitBounds[i][1], masterSplitBounds[i][2]);
      }
    }
    if (_splitRule == LOG_RANK_SCORE) {
      double *predictorValue = dvector(1, _observationSize);
      for (j=1; j <= _xSize; j++) {
        for (i=1; i <= _observationSize; i++) {
          masterSplitOrder[j][i] = 0;
          predictorValue[i] = _observation[j][_bootMembershipIndex[i]];
        }
        indexx(_observationSize, predictorValue, masterSplitOrder[j]);
      }
      free_dvector(predictorValue, 1, _observationSize);
      if (getTraceFlag() & DL3_TRACE) {
        Rprintf("\n\nBootstraped Ordered Predictor Matrix:  ");
        Rprintf("\n     index   observations -> \n");
        for (i = 1; i <= _observationSize; i++) {
          Rprintf("%10d", i);
          for (j=1; j <= _xSize; j++) {
            Rprintf(" %10d", masterSplitOrder[j][i]);
          }
          Rprintf("\n");
        }
      }
    }
    nrCopyMatrix(rootPtr -> permissibleSplit, masterSplitBounds, _xSize, 2);
    free_uivector(leadIndex, 1, _xSize);
    free_uimatrix(masterSplitBounds, 1, _xSize, 1, 2);
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Initialization of master split data complete.");  
    }
    result = TRUE;
  }  
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nbootstrap() EXIT ...\n");
  }
  return result;
}
SEXP rsfGrow(SEXP traceFlag,
             SEXP mup,  
             SEXP seedPtr,  
             SEXP splitRule,  
             SEXP randomCovariateCount,  
             SEXP forestSize,  
             SEXP minimumDeathCount,
             SEXP observationSize,
             SEXP time,
             SEXP status,
             SEXP xSize,
             SEXP xData,
             SEXP timeInterestSize,
             SEXP timeInterest,
             SEXP randomCovariateWeight,
             SEXP xType) {
  int seed1Value         = INTEGER(seedPtr)[0];
  int seed2Value         = INTEGER(seedPtr)[0];
  _seed1Ptr              = &seed1Value;
  _seed2Ptr              = &seed2Value;
  _splitRule            = INTEGER(splitRule)[0];
  _randomCovariateCount = INTEGER(randomCovariateCount)[0];
  _forestSize           = INTEGER(forestSize)[0];
  _minimumDeathCount    = INTEGER(minimumDeathCount)[0];
  _observationSize      = INTEGER(observationSize)[0];
  _time                 =    REAL(time);  _time--;
  _status               =    REAL(status);  _status--;
  _xSize                = INTEGER(xSize)[0];
  _xData                =    REAL(xData);
  _timeInterestSize     = INTEGER(timeInterestSize)[0];
  _timeInterest         =    REAL(timeInterest);  _timeInterest--;
  _randomCovariateWeight =   REAL(randomCovariateWeight);  _randomCovariateWeight--;
  _sexp_xType = xType;
  _mup                  = INTEGER(mup)[0] | MUP_PERF;  
  if (*_seed1Ptr >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Random seed must be less than zero.  \n");
    return R_NilValue;
  }
  if ( (_splitRule != LOG_RANK) && (_splitRule != CONSERVE_EVENTS) && (_splitRule != LOG_RANK_SCORE) && (_splitRule != LOG_RANK_APPROX) ) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Invalid split rule:  %10d \n", _splitRule);
    return R_NilValue;
  }
  if ( ((_randomCovariateCount < 1) || (_randomCovariateCount > _xSize)) ) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of random covariate parameters must be greater");
    Rprintf("\nRSF:  than zero and less than the total number of covariates:  %10d \n", _randomCovariateCount);
    return R_NilValue;
  }
  if (_minimumDeathCount < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Minimum number of deaths must be greater than zero:  %10d \n", _minimumDeathCount);
    return R_NilValue;
  }
  _treeID_ = NULL;
  _nodeID_ = NULL;
  _parmID_ = NULL;
  _spltPT_ = NULL;
  return rsf(RSF_GROW, INTEGER(traceFlag)[0]);
}
SEXP rsfPredict(SEXP traceFlag,
                SEXP mup,
                SEXP seedPtr,  
                SEXP forestSize, 
                SEXP observationSize,
                SEXP time,
                SEXP status,
                SEXP xSize,
                SEXP xData,
                SEXP fobservationSize,
                SEXP ftime,
                SEXP fstatus,
                SEXP fxData,
                SEXP timeInterestSize,
                SEXP timeInterest,
                SEXP treeID,
                SEXP nodeID,
                SEXP parmID,
                SEXP spltPT,
                SEXP seed,
                SEXP xType) {
  int seed1Value         = INTEGER(seedPtr)[0];
  int seed2Value         = INTEGER(seedPtr)[0];
  _seed1Ptr              = &seed1Value;
  _seed2Ptr              = &seed2Value;
  _forestSize           = INTEGER(forestSize)[0];
  _observationSize      = INTEGER(observationSize)[0];
  _time                 =    REAL(time);  _time --;
  _status               =    REAL(status);  _status --;
  _xSize                = INTEGER(xSize)[0];
  _xData                =    REAL(xData);
  _fobservationSize     = INTEGER(fobservationSize)[0];
  _ftime                =    REAL(ftime);  _ftime --;
  _fstatus              =    REAL(fstatus);  _fstatus --;
  _fxData               =    REAL(fxData);
  _timeInterestSize     = INTEGER(timeInterestSize)[0];
  _timeInterest         =    REAL(timeInterest);  _timeInterest --;
  _treeID_              = (uint*) INTEGER(treeID);  _treeID_ --;
  _nodeID_              = (uint*) INTEGER(nodeID);  _nodeID_ --;
  _parmID_              = (uint*) INTEGER(parmID);  _parmID_ --;
  _spltPT_              =    REAL(spltPT);  _spltPT_ --;
  _seed_                = INTEGER(seed);
  _sexp_xType = xType;
  _mup                  = INTEGER(mup)[0] & (~MUP_TREE);  
  if (_fobservationSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of observations in prediction must be at least one:  %10d \n", _fobservationSize);
    return R_NilValue;
  }
  if (*_seed1Ptr >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Random seed must be less than zero.  \n");
    return R_NilValue;
  }
  return rsf(RSF_PRED, INTEGER(traceFlag)[0]);
}
SEXP rsf(uint mode, uint traceFlag) {
  uint sexpIndex;
  uint sortedTimeInterestSize;
  char      mTimeIndexFlag;   
  char    **mRecordBootFlag;  
  double ***mvImputation;     
  double **masterSplit;        
  uint    *masterSplitSize;    
  uint   **masterSplitOrder;   
  Node   **root;
  double  *sImputeStatusPtr;
  double  *sImputeTimePtr;
  double **sImputePredictorPtr;
  double **oobEnsemblePtr;
  double **fullEnsemblePtr;
  double  *ensembleRun;
  uint    *ensembleDen;
  double **vimpEnsembleRun;
  uint    *oobSampleSize;
  uint     forestNodeCounter;
  uint     rejectedTreeCount;  
  uint i, j, k, b;
  Node *rootPtr;
  uint obsSize;
  double *statusPtr;
  double *timePtr;
  double concordanceIndex;
  char result;
  clock_t    start;
  clock_t    now;
  _traceFlagDiagLevel = (traceFlag | mode) & 0xFF;
  _traceFlagIterValue = traceFlag >> 8;
  updateTraceFlag(TRUE);
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Native code rsf() entry. \n");
  }
  start = clock();
  if (_forestSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of bootstrap iterations must be greater than zero:  %10d \n", _forestSize);
    return R_NilValue;
  }
  if (_observationSize < 2) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of observations must be greater than one:  %10d \n", _observationSize);
    return R_NilValue;
  }
  if (_timeInterestSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of time points of interest must be greater than zero:  %10d \n", _timeInterestSize);
    return R_NilValue;
  }
  if (_xSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of parameters must be greater than zero:  %10d \n", _xSize);
    return R_NilValue;
  }
  if (mode == RSF_GROW) {
    for (i = 1; i <= _xSize; i++) {
      if(_randomCovariateWeight[i] < 0) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Parameter verification failed.");
        Rprintf("\nRSF:  Random covariate weight elements must be greater than or equal to zero:  %12.4f \n", _randomCovariateWeight[i]);
        return R_NilValue;
      }
    }
  }
  if (mode == RSF_PRED) {
      if(*_seed_ >= 0) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Parameter verification failed.");
        Rprintf("\nRSF:  Random seed element must be less than zero:  %10d \n", *_seed_);
        return R_NilValue;
      }
  }
  if (getTraceFlag() & DL0_TRACE) {
    if (mode == RSF_GROW) {
      Rprintf("\nRSF:  Mode is GROW.");
      Rprintf("\nRSF:  Split rule is:                %10d", _splitRule);
      Rprintf("\nRSF:  Number of GROW observations:  %10d", _observationSize);
    }
    if (mode == RSF_PRED) {
      Rprintf("\nRSF:  Mode is PRED.");
      Rprintf("\nRSF:  Number of PRED observations:  %10d", _fobservationSize);
    }
    Rprintf("\nRSF:  Number of predictors:         %10d", _xSize);
  }
  stackPreDefinedCommonArrays(& oobSampleSize);
  if (mode == RSF_GROW) {
    stackPreDefinedGrowthArrays(
                                & masterSplit,
                                & masterSplitSize,
                                & masterSplitOrder);
  }
  else {
    stackPreDefinedPredictArrays();
  }
  initializeArrays(mode, &sortedTimeInterestSize);
  result =  stackAndInitializeMissingArrays(mode,
                                            & mTimeIndexFlag,
                                            & mRecordBootFlag,
                                            & mvImputation);
  if (result == FALSE) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    return R_NilValue;
  }
  sexpIndex = stackDefinedOutputObjects(mode,
                                        sortedTimeInterestSize,
                                        sexpString,
                                        & root,
                                        & oobEnsemblePtr,
                                        & fullEnsemblePtr,
                                        & ensembleRun,
                                        & ensembleDen,
                                        & _oobEnsemble_,
                                        & _fullEnsemble_,
                                        & _performance_,
                                        & _leafCount_,
                                        & _proximity_,
                                        & _varImportance_,
                                        & vimpEnsembleRun,
                                        & _seed_,
                                        & _imputation_,
                                        & sImputeStatusPtr,
                                        & sImputeTimePtr,
                                        & sImputePredictorPtr,
                                        & stackCount,
                                        sexpVector
                                        );
  if (mode == RSF_GROW) {
    ran1(_seed1Ptr);  
    *_seed1Ptr = -abs(*_seed1Ptr);
    if (_mup & MUP_TREE) {
      *_seed_ = *_seed1Ptr;
    }
  }
  else {
   *_seed1Ptr = *_seed_;
  }
  ran2(_seed2Ptr);
  ran2(_seed2Ptr);
  *_seed2Ptr = -abs(*_seed2Ptr);
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nRandom seed for ran1():  %20d", *_seed1Ptr);
    Rprintf("\nRandom seed for ran2():  %20d", *_seed2Ptr);
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Number of trees in forest:  %10d \n", _forestSize);
  }
  if (mode == RSF_GROW) {
    forestNodeCounter = 0;
  }
  else {
    forestNodeCounter = 1;
  }
  for (b = 1; b <= _forestSize; b++) {
    if ( (b == 1) || (b == _forestSize)) {
      updateTraceFlag(TRUE);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nStart of Tree:  %10d", b);
    }
    if (getTraceFlag()) {
      if (_mRecordSize > 0) {
        unImpute (RSF_GROW);
      }
      if (mode == RSF_PRED) {
        if (_fmRecordSize > 0) {
          unImpute (RSF_PRED);
        }
      }
    }
    rootPtr = makeNode();  
    rootPtr -> parent = rootPtr;
    rootPtr -> left = NULL;
    rootPtr -> right = NULL;
    rootPtr -> splitValueIndex = 0;
    rootPtr -> splitValue      = 0.0;
    rootPtr -> splitParameter = 0;
    rootPtr -> splitFlag = TRUE;
    rootPtr -> leafCount = 1;
    if (_mup & MUP_TREE) {
      root[b] = rootPtr;
    }
    if (mode == RSF_GROW) {
      result = bootstrap (RSF_GROW,
                          b,
                          rootPtr,
                          oobSampleSize,
                          masterSplit,
                          masterSplitSize,
                          masterSplitOrder,
                          mRecordBootFlag); 
      if (result) {
        _leafCount_[b] = 1;
        makeTree (b, 
                  rootPtr,
                  masterSplit,
                  masterSplitOrder,
                  mTimeIndexFlag);  
      }
      else {
        if (getTraceFlag() & DL0_TRACE) {
          Rprintf("\nRSF:  Tree rejected:  %10d", b);  
        }
      }
    }  
    else {
      result = bootstrap (RSF_PRED,
                          b,
                          rootPtr,
                          oobSampleSize,
                          masterSplit,      
                          masterSplitSize,  
                          masterSplitOrder, 
                          mRecordBootFlag); 
      if (result) {
        _leafCount_[b] = 1;
        if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
          Rprintf("\nIncoming Tree:  ");
          Rprintf("\n      tree       parm       node         splt \n");
        }
        restoreTree(b,
                    rootPtr,
                    _leafCount_ + b,
                    & forestNodeCounter,
                    _treeID_,
                    _nodeID_,
                    _parmID_,
                    _spltPT_);
        if (getTraceFlag() & DL0_TRACE) {
          Rprintf("\nRSF:  Tree construction complete:  %10d", b);  
          Rprintf("\nRSF:  Final leaf count:  %10d\n", _leafCount_[b]);
        }
        imputeTree(mode, b, rootPtr);
        if (mTimeIndexFlag == TRUE) {
          updateTimeIndexArray(NULL);
        }
      }  
      else {
        forestNodeCounter ++;
        if (getTraceFlag() & DL0_TRACE) {
          Rprintf("\nRSF:  Tree rejected:  %10d", b);  
        }
      }
    }  
    if (result) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nFinal GROW Membership (all data):  %10d", b);
        Rprintf("\n     index       leaf    IBindex\n");
        for (i=1; i <=  _observationSize; i++) {
          Rprintf("%10d %10d %10d \n", i, _nodeMembership[i] -> leafCount, _bootMembershipIndex[i]);
        }
      }
      if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
        Rprintf("\nFinal PRED Membership (all data):  %10d", b);
        Rprintf("\n     index       leaf\n");
        for (i=1; i <=  _fobservationSize; i++) {
          Rprintf("%10d %10d \n", i, _fnodeMembership[i] -> leafCount);
        }
      }
      if (_mup & MUP_PROX) {
        if (mode == RSF_GROW) {
          k = 0;
          for (i = 1; i <= _observationSize; i++) {
            k += i - 1;
            for (j = 1; j <= i; j++) {
              if ( (_nodeMembership[i] -> leafCount) == (_nodeMembership[j] -> leafCount) ) {
                _proximity_[k + j] ++;
              }
            }
          }
        }
        else {
          k = 0;
          for (i = 1; i <= _fobservationSize; i++) {
            k += i - 1;
            for (j = 1; j <= i; j++) {
              if ( (_fnodeMembership[i] -> leafCount) == (_fnodeMembership[j] -> leafCount) ) {
                _proximity_[k + j] ++;
              }
            }
          }
        }
        if (getTraceFlag() & DL0_TRACE) {
          Rprintf("\nRSF:  Proximity matrix calculation complete.");  
        }
      }  
      if (mode == RSF_GROW) {
        if (_mRecordSize > 0) {
          imputeUpdate(RSF_GROW, mvImputation[b]);
        }
      }
      else {
        if (_fmRecordSize > 0) {
          imputeUpdate(RSF_PRED, mvImputation[b]);
        }
      }
      _performance_[b] = getPerformance(mode,
                                        sortedTimeInterestSize,
                                        _leafCount_[b],
                                        oobEnsemblePtr,
                                        fullEnsemblePtr,
                                        ensembleRun,
                                        ensembleDen,
                                        oobSampleSize[b],
                                        rootPtr,
                                        vimpEnsembleRun,
                                        b,
                                        mRecordBootFlag,
                                        mvImputation);
    }  
    if (!(_mup & MUP_TREE)) {
      freeTree(rootPtr);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nEnd   of Tree:  %10d", b);
    }
    updateTraceFlag(FALSE);
  }  
  updateTraceFlag(TRUE);
  if (getTraceFlag()) {
    if (_mRecordSize > 0) {
      unImpute (RSF_GROW);
    }
    if (mode == RSF_PRED) {
      if (_fmRecordSize > 0) {
        unImpute (RSF_PRED);
      }
    }
  }
  rejectedTreeCount = k = 0;
  for (b = 1; b <= _forestSize; b++) {
    if (_leafCount_[b] == 0) {
      rejectedTreeCount ++;
    }
    if (_leafCount_[b] == 1) {
      k ++;
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nFinal Tree Leaf Counts:  ");
    Rprintf("\n         tree    leafCount");
    for (b = 1; b <= _forestSize; b++) {
      Rprintf("\n %12d %12d", b, _leafCount_[b]);
    }
    Rprintf("\n");
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Trees rejected:  %10d ", rejectedTreeCount);
    Rprintf("\nRSF:  Trees stumped:   %10d ", k);
    Rprintf("\nRSF:  Trees (total):   %10d \n", _forestSize);
  }
  if (getTraceFlag() & DL1_TRACE) {
    if (_mup & MUP_PROX) {
      Rprintf("\nProximity Matrix:  \n");
      k = 0;
      for (i = 1; i <= _observationSize; i++) {
        k += i - 1;
        for (j = 1; j <= i; j++) {
          Rprintf("%10d ", _proximity_[k + j]);
        }
        Rprintf("\n");
      }
    }
  }
  if (mode == RSF_GROW) {
    if (_mRecordSize > 0) {
      imputeSummary(mode,
                    mRecordBootFlag,
                    mvImputation,
                    sImputeStatusPtr,
                    sImputeTimePtr,
                    sImputePredictorPtr);
    }
  }
  else {
    if (_fmRecordSize > 0) {
      imputeSummary(mode,
                    mRecordBootFlag,
                    mvImputation,
                    sImputeStatusPtr,
                    sImputeTimePtr,
                    sImputePredictorPtr);
    }
  }
  if (rejectedTreeCount < _forestSize) {
    if (_mup & MUP_VIMP) {
      if (mode == RSF_GROW) {
        obsSize = _observationSize;
        statusPtr = _status;
        timePtr = _time;
      }
      else {
        obsSize = _fobservationSize;
        statusPtr = _fstatus;
        timePtr = _ftime;
      }
      for (k=1; k <= _xSize; k++) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDen[i] != 0) {
            vimpEnsembleRun[k][i] = vimpEnsembleRun[k][i] / ensembleDen[i];
          }
        }
        if (getTraceFlag() & DL1_TRACE) {
          Rprintf("\nConcordance Calculation for VIMP covariate:  %10d", k);  
        }
        concordanceIndex = getConcordanceIndex(obsSize, 
                                               statusPtr,
                                               timePtr,
                                               vimpEnsembleRun[k], 
                                               ensembleDen);
        if (ISNA(concordanceIndex)) {
          _varImportance_[k] = NA_REAL;
        }
        else {
          _varImportance_[k] = 1 - concordanceIndex;
        }
      }
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Variable Importance Measure: \n");
        Rprintf("          ");  
        for (k=1; k <= _xSize; k++) {
          Rprintf("%10d", k);
        }
        Rprintf("\n          ");
        for (k=1; k <= _xSize; k++) {
          Rprintf("%10.4f", _varImportance_[k]);
        }
        Rprintf("\n");
      }
    }
    if (mode == RSF_GROW) {
      for (i = 1; i <= _observationSize; i++) {
        for (j=1; j <= sortedTimeInterestSize; j++) {
          if (ensembleDen[i] != 0) {
            oobEnsemblePtr[j][i] = oobEnsemblePtr[j][i] / ensembleDen[i];
          }
          fullEnsemblePtr[j][i] = fullEnsemblePtr[j][i] / (_forestSize - rejectedTreeCount);
        }
      }
    }
    else {
      for (i = 1; i <= _fobservationSize; i++) {
        for (j=1; j <= sortedTimeInterestSize; j++) {
          fullEnsemblePtr[j][i] = fullEnsemblePtr[j][i] / (_forestSize - rejectedTreeCount);
        }
      }
    }
  }  
  else {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Insufficient trees for analysis.  \n");
  }
  forestNodeCounter = 0;
  if (mode == RSF_GROW) {
    if (_mup & MUP_TREE) {
      for (b = 1; b <= _forestSize; b++) {
        if (_leafCount_[b] > 0) {
          forestNodeCounter += (2 * _leafCount_[b]) - 1;
        }
        else {
          forestNodeCounter ++;
        }
      }
    }
  }
  sexpIndex = 
    stackVariableOutputObjects(forestNodeCounter,  
                               & _treeID_,         
                               & _nodeID_,         
                               & _parmID_,         
                               & _spltPT_,         
                               sexpIndex, 
                               sexpString,
                               sexpVector);
  if (mode == RSF_GROW) {
    if (_mup & MUP_TREE) {
      forestNodeCounter = 1;
      for (b = 1; b <= _forestSize; b++) {
        saveTree(b, root[b], & forestNodeCounter, _treeID_, _nodeID_, _parmID_, _spltPT_);
      }
      forestNodeCounter --;
      if (getTraceFlag() & DL3_TRACE) {
        Rprintf("\nGrown Forest Node Count:  %10d", forestNodeCounter);
        Rprintf("\nGrown Forest Output:  ");
        Rprintf("\n   TREE_ID    NODE_ID    PARM_ID      SPLT_PT");
        for (i = 1; i <= forestNodeCounter; i++) {
          Rprintf("\n%10d %10d %10d %12.4f",
                  _treeID_[i],
                  _nodeID_[i],
                  _parmID_[i],
                  _spltPT_[i]);
        }
        Rprintf("\n");
        Rprintf("\nGrown Forest Membership by Observation:  ");
        Rprintf("\n    OBSERV    TREE_ID    NODE_ID");
        for (b = 1; b <= _forestSize; b++) {
          for (k = 1; k <= _observationSize; k++) {
            Rprintf("\n%10d %10d %10d ", k, b, getMembership(root[b], _observation, k) -> leafCount);
          }
        }
        Rprintf("\n");
      }
    }  
  }  
  unstackPreDefinedCommonArrays(oobSampleSize);
  if (mode == RSF_GROW) {
    unstackPreDefinedGrowthArrays(masterSplit,
                                  masterSplitSize,
                                  masterSplitOrder
                                  );
  }
  else {
    unstackPreDefinedPredictArrays();
  }
  unstackMissingArrays(mode,
                       mRecordBootFlag,
                       mvImputation);
  unstackDefinedOutputObjects(mode,
                              root,
                              sortedTimeInterestSize,
                              oobEnsemblePtr,
                              fullEnsemblePtr,
                              ensembleRun,
                              ensembleDen,
                              vimpEnsembleRun,
                              sImputePredictorPtr);
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Native code rsf() exit. \n");
    now = updateTimeStamp(start);
  }
  UNPROTECT(stackCount + 2);
  return sexpVector[RSF_OUTP_ID];
}
