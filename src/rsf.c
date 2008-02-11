//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.2.0
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

#include    "global.h"
#include    "nrutil.h"
#include  "node_ops.h"
#include "rsfImpute.h"
#include  "rsfSplit.h"
#include   "rsfUtil.h"
#include   "rsfBootstrap.h"
#include  "rsfStack.h"
#include       "rsf.h"
uint stackCount;
char *sexpString[RSF_SEXP_CNT] = {
  "",              
  "",              
  "fullEnsemble",  
  "oobEnsemble",   
  "performance",   
  "proximity",     
  "leafCount",     
  "treeID",        
  "nodeID",        
  "parmID",        
  "spltPT",        
  "seed",          
  "importance",    
  "imputation",    
  "oobImputation", 
  "varUsed"        
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
double   *_importance_;
double   *_imputation_;
double   *_oobImputation_;
uint     *_varUsed_;
uint      _opt;
uint      _splitRule;
uint      _imputeSize;
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
uint      _intrPredictorSize;
uint     *_intrPredictor;
uint     *_intrObservation;
char     **_xType;
double  **_observation;
double  **_fobservation;
double   *_masterTime;
uint     *_masterTimeIndex;
uint      _masterTimeSize;
char     *_importanceFlag;
char      _mTimeIndexFlag; 
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
double **_oobEnsemblePtr;
double **_fullEnsemblePtr;
double  *_ensembleRun;
double **_vimpEnsembleRun;  
uint    *_oobEnsembleDen;
uint    *_fullEnsembleDen;
double  *_sImputeStatusPtr;
double  *_sImputeTimePtr;
double **_sImputePredictorPtr;
double  *_sOOBImputeStatusPtr;
double  *_sOOBImputeTimePtr;
double **_sOOBImputePredictorPtr;
uint   **_varUsedPtr;
int      *_seed1Ptr;
int      *_seed2Ptr;
Node    **_nodeMembership;
uint     *_bootMembershipIndex;
char     *_bootMembershipFlag;
uint     *_oobSampleSize;
Node    **_fnodeMembership;
char     *_fbootMembershipFlag;
uint     *_foobSampleSize;
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
  double **masterSplit) {
  char result;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetBestSplit() ENTRY ...\n");
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nAttempting to split node:  %10d", parent -> leafCount);
  }
  switch(_splitRule) {
  case LOG_RANK:
    result = logRank(parent,
                     splitParameterMax,
                     splitValueMax,
                     masterSplit);
    break;
  case CONSERVE_EVENTS:
    result = conserveEvents(parent,
                            splitParameterMax,
                            splitValueMax,
                            masterSplit);
    break;
  case LOG_RANK_SCORE:
    result = logRankScore(parent,
                          splitParameterMax,
                          splitValueMax,
                          masterSplit);
    break;
  case LOG_RANK_APPROX:
    result = logRankApprox(parent,
                           splitParameterMax,
                           splitValueMax,
                           masterSplit);
    break;
  case RANDOM_SPLIT:
    result = randomSplit(parent,
                         splitParameterMax,
                         splitValueMax,
                         masterSplit);
    break;
  case LOG_RANK_RANDOM:
    result = logRankRandom(parent,
                           splitParameterMax,
                           splitValueMax,
                           masterSplit);
    break;
default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid split rule:  %10d", _splitRule);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetBestSplit() EXIT ...\n");
  }
  return result;
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
char makeTree (char     multipleImputeFlag,
               uint     b,
               Node    *parent,
               double **masterSplit) {
  char splitResult, forkResult;
  uint splitParameterMax, splitValueMax;
  uint i;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nmakeTree() ENTRY ...\n");
  }
  if (multipleImputeFlag == FALSE) {
    if (_mRecordSize > 0) {
      imputeNode(RSF_GROW,
                 TRUE,
                 b, 
                 parent);
      if (_mTimeIndexFlag == TRUE) {
        updateTimeIndexArray(parent);
      }
    }
  }
  splitResult = getBestSplit(parent,
                             & splitParameterMax,
                             & splitValueMax,
                             masterSplit);
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
        makeTree (multipleImputeFlag,
                  b,
                  parent -> left,
                  masterSplit);
      }
      else {
        if (multipleImputeFlag == FALSE) {
          if (_mRecordSize > 0) {
            imputeNode(RSF_GROW,
                       TRUE,
                       b, 
                       parent -> left);
            if (_mTimeIndexFlag == TRUE) {
              updateTimeIndexArray(parent -> left);
            }
          }
        }
      }
      if ((parent -> right) -> splitFlag == TRUE) {
        if (getTraceFlag() & DL1_TRACE) {
          Rprintf("\nmakeTree() RIGHT:  \n");
        }
        makeTree (multipleImputeFlag,
                  b,
                  parent -> right,
                  masterSplit);
      }
      else {
        if (multipleImputeFlag == FALSE) {
          if (_mRecordSize > 0) {
            imputeNode(RSF_GROW,
                       TRUE,
                       b, 
                       parent -> right);
            if (_mTimeIndexFlag == TRUE) {
              updateTimeIndexArray(parent -> right);
            }
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
SEXP rsfGrow(SEXP traceFlag,
             SEXP opt,  
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
             SEXP xType,
             SEXP imputeSize) {
  uint i;
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
  _imputeSize           = INTEGER(imputeSize)[0];
  _sexp_xType = xType;
  _opt                  = INTEGER(opt)[0];
  _opt                  = _opt | OPT_FENS;  
  _opt                  = _opt | OPT_OENS;  
  _opt                  = _opt | OPT_PERF;  
  _opt                  = _opt | OPT_LEAF;  
  _opt                  = _opt | OPT_MISS;  
  _opt                  = _opt | OPT_OMIS;  
  if (_opt & OPT_TREE) {
    _opt = _opt | OPT_SEED;
  }
  else {
    _opt = _opt & (~OPT_SEED);
  }
  if (*_seed1Ptr >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Random seed must be less than zero.  \n");
    return R_NilValue;
  }
  if ( (_splitRule != LOG_RANK) && 
       (_splitRule != CONSERVE_EVENTS) && 
       (_splitRule != LOG_RANK_SCORE) && 
       (_splitRule != LOG_RANK_APPROX) && 
       (_splitRule != RANDOM_SPLIT) &&
       (_splitRule != LOG_RANK_RANDOM) ) {
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
  for (i = 1; i <= _xSize; i++) {
    if(_randomCovariateWeight[i] < 0) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Parameter verification failed.");
      Rprintf("\nRSF:  Random covariate weight elements must be greater than or equal to zero:  %12.4f \n", _randomCovariateWeight[i]);
      return R_NilValue;
    }
  }
  _treeID_ = NULL;
  _nodeID_ = NULL;
  _parmID_ = NULL;
  _spltPT_ = NULL;
  return rsf(RSF_GROW, INTEGER(traceFlag)[0]);
}
SEXP rsfPredict(SEXP traceFlag,
                SEXP opt,
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
  _opt                  = INTEGER(opt)[0];
  _opt                  = _opt | OPT_FENS;  
  _opt                  = _opt & (~OPT_OENS);  
  _opt                  = _opt | OPT_LEAF;  
  _opt = _opt & (~OPT_TREE);  
  _opt = _opt & (~OPT_SEED);  
  _opt                  = _opt | OPT_MISS;  
  _opt                  = _opt & (~OPT_OMIS);  
  _imputeSize = 1;
  if (*_seed1Ptr >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  User random seed must be less than zero.  \n");
    return R_NilValue;
  }
  if(*_seed_ >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Forest random seed element must be less than zero:  %10d \n", *_seed_);
    return R_NilValue;
  }
  if (_fobservationSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of individuals in prediction must be at least one:  %10d \n", _fobservationSize);
    return R_NilValue;
  }
  return rsf(RSF_PRED, INTEGER(traceFlag)[0]);
}
SEXP rsfInteraction(SEXP traceFlag,
                    SEXP opt,
                    SEXP seedPtr,  
                    SEXP forestSize, 
                    SEXP observationSize,
                    SEXP time,
                    SEXP status,
                    SEXP xSize,
                    SEXP xData,
                    SEXP timeInterestSize,
                    SEXP timeInterest,
                    SEXP treeID,
                    SEXP nodeID,
                    SEXP parmID,
                    SEXP spltPT,
                    SEXP seed,
                    SEXP xType,
                    SEXP intrPredictorSize,
                    SEXP intrPredictor,
                    SEXP fobservationSize,
                    SEXP intrObservation) {
  uint i;
  uint leadingIndex;
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
  _timeInterestSize     = INTEGER(timeInterestSize)[0];
  _timeInterest         =    REAL(timeInterest);  _timeInterest --;
  _treeID_              = (uint*) INTEGER(treeID);  _treeID_ --;
  _nodeID_              = (uint*) INTEGER(nodeID);  _nodeID_ --;
  _parmID_              = (uint*) INTEGER(parmID);  _parmID_ --;
  _spltPT_              =    REAL(spltPT);  _spltPT_ --;
  _seed_                = INTEGER(seed);
  _sexp_xType = xType;
  _intrPredictorSize         = INTEGER(intrPredictorSize)[0];
  _intrPredictor             = (uint*) INTEGER(intrPredictor);  _intrPredictor --;
  _fobservationSize     = INTEGER(fobservationSize)[0];
  _intrObservation      = (uint*) INTEGER(intrObservation); _intrObservation --;
  _opt                  = INTEGER(opt)[0];
  _opt                  = _opt & (~OPT_FENS);  
  _opt                  = _opt | OPT_OENS;  
  _opt                  = _opt | OPT_PERF;  
  _opt                  = _opt | OPT_LEAF;  
  _opt                  = _opt & (~OPT_MISS);  
  _opt                  = _opt | OPT_OMIS;  
  _opt = _opt & (~OPT_TREE);  
  _opt = _opt & (~OPT_SEED);  
  _imputeSize = 1;
  if (*_seed1Ptr >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  User random seed must be less than zero.  \n");
    return R_NilValue;
  }
  if(*_seed_ >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Forest random seed element must be less than zero:  %10d \n", *_seed_);
    return R_NilValue;
  }
  if((_intrPredictorSize <= 0) || (_intrPredictorSize > _xSize)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of predictors to be perturbed must be greater than zero and less than %10d:  %10d \n", _xSize, _intrPredictorSize);
    return R_NilValue;
  }
  uint *intrPredictorCopy = uivector(1, _intrPredictorSize);
  for (i=1; i <= _intrPredictorSize; i++) {
    intrPredictorCopy[i] = _intrPredictor[i];
  }
  hpsortui(intrPredictorCopy, _intrPredictorSize);
  leadingIndex = 1;
  for (i=2; i <= _intrPredictorSize; i++) {
    if (intrPredictorCopy[i] > intrPredictorCopy[leadingIndex]) {
      leadingIndex++;
    }
  }
  if (_intrPredictorSize != leadingIndex) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Interaction terms are not unique.");
    Rprintf("\nRSF:  Only %10d of %10d are unique.", leadingIndex, _intrPredictorSize);
  }
  free_uivector(intrPredictorCopy, 1, _intrPredictorSize);
  for (i=1; i <= _intrPredictorSize; i++) {
    if (_intrPredictor[i] > _xSize) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Parameter verification failed.");
      Rprintf("\nRSF:  Interaction terms are not coherent.");
      Rprintf("\nRSF:  Predictor encountered is %10d, maximum allowable is %10d.", _intrPredictor[i], _xSize);
    }
  }
  if((_fobservationSize <= 0) || (_fobservationSize > _observationSize)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of individuals in INTR data set must be greater than zero and less than %10d:  %10d \n", _observationSize, _fobservationSize);
    return R_NilValue;
  }
  hpsortui(_intrObservation, _fobservationSize);
  leadingIndex = 1;
  for (i=2; i <= _fobservationSize; i++) {
    if (_intrObservation[i] > _intrObservation[leadingIndex]) {
      leadingIndex++;
    }
  }
  if (_fobservationSize != leadingIndex) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Individuals in INTR data subset are not unique.");
    Rprintf("\nRSF:  Only %10d of %10d are unique.", leadingIndex, _fobservationSize);
  }
  for (i=1; i <= _fobservationSize; i++) {
    if (_intrObservation[i] > _observationSize) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Parameter verification failed.");
      Rprintf("\nRSF:  Individuals in INTR data subset are not coherent.");
      Rprintf("\nRSF:  Individual encountered is %10d, maximum allowable is %10d.", _intrObservation[i], _observationSize);
    }
  }
  return rsf(RSF_INTR, INTEGER(traceFlag)[0]);
}
SEXP rsf(char mode, uint traceFlag) {
  uint sexpIndex;
  uint sortedTimeInterestSize;
  char      multipleImputeFlag;     
  char      concordanceImputeFlag;  
  char      updateFlag;       
  char    **dmRecordBootFlag;  
  double ***dmvImputation;     
  double **masterSplit;        
  uint    *masterSplitSize;    
  Node   **root;
  uint     forestNodeCounter;
  uint     rejectedTreeCount;  
  Node *rootPtr;
  uint obsSize;
  uint varSize;
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  uint    *ensembleDenPtr;
  double concordanceIndex;
  char result;
  clock_t    start;
  clock_t    now;
  uint i, j, k, b, r, s;
  statusPtr = NULL;     
  timePtr   = NULL;     
  predictorPtr = NULL;  
  _traceFlagDiagLevel = (traceFlag | mode) & 0xFF;
  _traceFlagIterValue = traceFlag >> 8;
  updateTraceFlag(TRUE);
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Native code rsf() entry. \n");
    Rprintf("\nRSF:  OPT %10x \n", _opt); 
  }
  start = clock();
  if (_imputeSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number imputations must be greater than zero:  %10d \n", _forestSize);
    return R_NilValue;
  }
  if (_forestSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of bootstrap iterations must be greater than zero:  %10d \n", _forestSize);
    return R_NilValue;
  }
  if (_observationSize < 2) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of individuals must be greater than one:  %10d \n", _observationSize);
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
  if (getTraceFlag() & DL0_TRACE) {
    switch (mode) {
    case RSF_GROW:
      Rprintf("\nRSF:  Mode is GROW.");
      Rprintf("\nRSF:  Split rule is:               %10d", _splitRule);
      Rprintf("\nRSF:  Number of GROW individuals:  %10d", _observationSize);
      Rprintf("\nRSF:  Number of impute iterations: %10d", _imputeSize);
      break;
    case RSF_PRED:
      Rprintf("\nRSF:  Mode is PRED.");
      Rprintf("\nRSF:  Number of PRED individuals:  %10d", _fobservationSize);
      break;
    case RSF_INTR:
      Rprintf("\nRSF:  Mode is INTR.");
      Rprintf("\nRSF:  Number of INTR individuals:  %10d", _fobservationSize);
      break;
    default:
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Unknown case in switch encountered. ");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
      break;
    }
    Rprintf("\nRSF:  Number of predictors:        %10d", _xSize);
    Rprintf("\nRSF:  Number of trees in forest:   %10d \n", _forestSize);
  }
  stackPreDefinedCommonArrays();
  switch (mode) {
  case RSF_GROW:
    stackPreDefinedGrowthArrays(& masterSplit,
                                & masterSplitSize);
    break;
  case RSF_PRED:
    stackPreDefinedPredictArrays();
    break;
  case RSF_INTR:
    stackPreDefinedInteractionArrays();
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  initializeArrays(mode, &sortedTimeInterestSize);
  result =  stackMissingArrays(mode,
                               & dmRecordBootFlag,
                               & dmvImputation);
  if (result == FALSE) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Missingness Parameter verification failed.");
    return R_NilValue;
  }
  sexpIndex = stackDefinedOutputObjects(mode,
                                        sortedTimeInterestSize,
                                        sexpString,
                                        & root,
                                        & _oobEnsemble_,
                                        & _fullEnsemble_,
                                        & _performance_,
                                        & _leafCount_,
                                        & _proximity_,
                                        & _importance_,
                                        & _seed_,
                                        & _imputation_,
                                        & _oobImputation_,
                                        & _sImputeStatusPtr,
                                        & _sImputeTimePtr,
                                        & _sImputePredictorPtr,
                                        & _sOOBImputeStatusPtr,
                                        & _sOOBImputeTimePtr,
                                        & _sOOBImputePredictorPtr,
                                        & _varUsed_,
                                        & _varUsedPtr,
                                        & stackCount,
                                        sexpVector
                                        );
  if (mode == RSF_GROW) {
    ran1(_seed1Ptr);  
    *_seed1Ptr = -abs(*_seed1Ptr);
  }
  else {
    *_seed1Ptr = *_seed_;
  }
  ran2(_seed2Ptr);
  ran2(_seed2Ptr);
  *_seed2Ptr = -abs(*_seed2Ptr);
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nStart of random seed chain for ran1():  %20d", *_seed1Ptr);
    Rprintf("\nStart of random seed chain for ran2():  %20d", *_seed2Ptr);
  }
  if (mode == RSF_GROW) {
    forestNodeCounter = 0;
  }
  else {
    forestNodeCounter = 1;
  }
  multipleImputeFlag = FALSE;
  for (r = 1; r <= _imputeSize; r++) {
    if (getTraceFlag() & DL1_TRACE) {
      if (mode == RSF_GROW) {
        Rprintf("\nStart of impute iteration:  %10d", r);
      }
    }
      if (r == _imputeSize) {
        if (_opt & OPT_TREE) { 
        *_seed1Ptr = -abs(*_seed1Ptr);
        *_seed_ = *_seed1Ptr;
      }
    }
    if (mode == RSF_GROW) {
      if (r > 1) {
        multipleImputeFlag = TRUE;
      } 
    }
    for (b = 1; b <= _forestSize; b++) {
      if (b == _forestSize) {
        updateTraceFlag(TRUE);
      }
      if (getTraceFlag() & DL1_TRACE) {
        Rprintf("\nStart of iteration:  (%10d, %10d)", r, b);
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nTree random seed for ran1():  %20d", *_seed1Ptr);
        Rprintf("\nTree random seed for ran2():  %20d", *_seed2Ptr);
      }
      if (TRUE) {
        unImpute (mode);
      }
        if (mode == RSF_GROW) {
          if (_imputeSize > 1) {
          if (r == 2) {
            if (_mRecordSize > 0) {
              imputeUpdateShadow(RSF_GROW, 
                                 FALSE, 
                                 dmvImputation, 
                                 _status, 
                                 _time, 
                                 _observation);
            }
          }  
          else {
            if (_mRecordSize > 0) {
              imputeUpdateShadow(RSF_GROW, 
                                 ACTIVE, 
                                 dmvImputation, 
                                 _status, 
                                 _time, 
                                 _observation);
            }
          }  
        }
      }
      rootPtr = makeNode();  
      rootPtr -> parent = rootPtr;
      rootPtr -> left  = NULL;
      rootPtr -> right = NULL;
      rootPtr -> splitValueIndex = 0;
      rootPtr -> splitValue      = 0.0;
      rootPtr -> splitParameter  = 0;
      rootPtr -> splitFlag       = TRUE;
      rootPtr -> leafCount = 1;
      if ((_opt & OPT_TREE) && (r == _imputeSize)) {
        root[b] = rootPtr;
      }
      if (mode == RSF_GROW) {
        result = bootstrap (multipleImputeFlag,
                            mode,
                            b,
                            rootPtr,
                            masterSplit,
                            masterSplitSize,
                            dmRecordBootFlag); 
        if (result) {
          _leafCount_[b] = 1;
          if (getTraceFlag() & DL1_TRACE) {
            Rprintf("\nStart of makeTree():  (%10d, %10d)", r, b);
          }
          makeTree (multipleImputeFlag,
                    b, 
                    rootPtr,
                    masterSplit);
          if (getTraceFlag() & DL1_TRACE) {
            Rprintf("\nEnd of makeTree():  (%10d, %10d)", r, b);
          }
          if (getTraceFlag() & DL0_TRACE) {
            Rprintf("\nRSF:  Tree construction complete:  (%10d, %10d)", r, b);  
            Rprintf("\nRSF:  Final leaf count:  %10d\n", _leafCount_[b]);
          }
          if ((_mRecordSize > 0) && (r > 1)) {
            for (j = 1; j <= _leafCount_[b]; j++) {
              k = 0;
              for (i = 1; i <= _observationSize; i++) {
                if ((_nodeMembership[i] -> leafCount) == j) {
                  k = i;
                  i = _observationSize;
                }
              }
              if (k == 0) {
                Rprintf("\nRSF:  *** ERROR *** ");
                Rprintf("\nRSF:  Multiply imputed node membership not found.");
                Rprintf("\nRSF:  Please Contact Technical Support.");
                Rprintf("\nRSF:  The application will now exit.\n");
                Rprintf("\n   imputation         tree         leaf \n");
                Rprintf(" %12d %12d %12d \n", r, b, j);
                Rprintf("\nDiagnostic Trace of (individual, boot, node, leaf) vectors in data set:  ");
                Rprintf("\n        index         boot         node         leaf \n");
                for (s = 1; s <= _observationSize; s++) {
                  Rprintf(" %12d %12d %12x %12d \n", s, _bootMembershipFlag[s], _nodeMembership[s], _nodeMembership[s] -> leafCount);
                }
                exit(TRUE);
              }
              else {
                imputeNode(mode, 
                           FALSE, 
                           b, 
                           _nodeMembership[k]);
              }  
            }  
          }  
        }  
        else {
          if (getTraceFlag() & DL0_TRACE) {
            Rprintf("\nRSF:  Tree rejected:  %10d", b);  
          }
        }
      }  
      else {
        result = bootstrap (multipleImputeFlag,
                            mode,
                            b,
                            rootPtr,
                            masterSplit,      
                            masterSplitSize,  
                            dmRecordBootFlag); 
        if (result) {
          _leafCount_[b] = 1;
          if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
            Rprintf("\nIncoming Tree:  ");
            Rprintf("\n      tree       parm       node         splt \n");
          }
          if (getTraceFlag() & DL1_TRACE) {
            Rprintf("\nStart of restore process:  (%10d, %10d)", r, b);
          }
          restoreTree(b,
                      rootPtr,
                      _leafCount_ + b,
                      & forestNodeCounter,
                      _treeID_,
                      _nodeID_,
                      _parmID_,
                      _spltPT_);
          if (getTraceFlag() & DL1_TRACE) {
            Rprintf("\nEnd of restore process:  (%10d, %10d)", r, b);
          }
          if (getTraceFlag() & DL0_TRACE) {
            Rprintf("\nRSF:  Tree construction complete:  (%10d, %10d)", r, b);  
            Rprintf("\nRSF:  Final leaf count:  %10d\n", _leafCount_[b]);
          }
          imputeTree(mode, b, rootPtr);
          if (_mTimeIndexFlag == TRUE) {
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
          Rprintf("\n     index       leaf  bootIndex\n");
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
        if ((getTraceFlag() & RSF_INTR) && (getTraceFlag() & DL2_TRACE)) {
          Rprintf("\nFinal INTR Membership (all data):  %10d", b);
          Rprintf("\n     index       leaf\n");
          for (i=1; i <=  _fobservationSize; i++) {
            Rprintf("%10d %10d \n", i, _fnodeMembership[i] -> leafCount);
          }
        }
        if (r == _imputeSize) {
          if (_opt & OPT_PROX) {
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
              Rprintf("\nRSF:  Proximity calculation complete.");  
            }
          }  
        }  
        updateFlag = FALSE;
        switch (mode) {
        case RSF_GROW:
          if ((r == 1) || (r < _imputeSize)) {
            if (_mRecordSize > 0) {
              updateFlag = TRUE;
              statusPtr = _status;
              timePtr = _time;
              predictorPtr = _observation;
            }
          }
          break;
        case RSF_PRED:
          if (_fmRecordSize > 0) {
            updateFlag = TRUE;
            statusPtr = _fstatus;
            timePtr = _ftime;
            predictorPtr = _fobservation;
          }
          break;
        case RSF_INTR:
          if (_fmRecordSize > 0) {
            updateFlag = TRUE;
            statusPtr = _fstatus;
            timePtr = _ftime;
            predictorPtr = _fobservation;
          }
          break;
        default:
          Rprintf("\nRSF:  *** ERROR *** ");
          Rprintf("\nRSF:  Unknown case in switch encountered. ");
          Rprintf("\nRSF:  Please Contact Technical Support.");
          Rprintf("\nRSF:  The application will now exit.\n");
          exit(TRUE);
          break;
        }
        if (updateFlag == TRUE) {
          imputeUpdateSummary(mode, 
                              statusPtr, 
                              timePtr, 
                              predictorPtr, 
                              dmvImputation[b]);
        }
        if (r == _imputeSize) {
          updateEnsembleEvents(multipleImputeFlag,
                               mode,
                               sortedTimeInterestSize,
                               rootPtr,
                               b,
                               dmRecordBootFlag,
                               dmvImputation);
        }
      }  
      if (!(_opt & OPT_TREE) || !(r == _imputeSize)) {
        freeTree(rootPtr);
      }
      if (getTraceFlag() & DL1_TRACE) {
        Rprintf("\nEnd of iteration:  (%10d, %10d)", r, b);
      }
      updateTraceFlag(FALSE);
    }  
    updateTraceFlag(TRUE);
    if (mode == RSF_GROW) {
      if (r < _imputeSize) {
        if (r == 1) {
          if (_mRecordSize > 0) {
            imputeSummary(RSF_GROW,
                          FALSE,
                          dmRecordBootFlag,
                          dmvImputation);
          }  
        }  
        else {
          if (_mRecordSize > 0) {
            imputeSummary(RSF_GROW,
                          ACTIVE,
                          dmRecordBootFlag,
                          dmvImputation);
          }
        }  
      }  
    } 
    if (getTraceFlag() & DL1_TRACE) {
      if (getTraceFlag() & DL1_TRACE) {
        Rprintf("\nEnd of impute iteration:  %10d", r);
      }
    }
  }  
  if (TRUE) {
    unImpute (mode);
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
    if (_opt & OPT_PROX) {
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
  switch (mode) {
  case RSF_GROW:
    if (_imputeSize == 1) {
      if (_mRecordSize > 0) {
        imputeSummary(RSF_GROW,
                      TRUE,
                      dmRecordBootFlag,
                      dmvImputation);
        imputeSummary(RSF_GROW,
                      FALSE,
                      dmRecordBootFlag,
                      dmvImputation);
      }
    }  
    else {
      if (_mRecordSize > 0) {
        imputeSummary(RSF_GROW,
                      ACTIVE,
                      dmRecordBootFlag,
                      dmvImputation);
      }
    }  
    break;
  case RSF_PRED:
    if (_fmRecordSize > 0) {
      imputeSummary(RSF_PRED,
                    ACTIVE,
                    dmRecordBootFlag,
                    dmvImputation);
    }
    break;
  case RSF_INTR:
    if (_fmRecordSize > 0) {
      imputeSummary(RSF_INTR,
                    FALSE,
                    dmRecordBootFlag,
                    dmvImputation);
    }
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  if (rejectedTreeCount < _forestSize) {
    if (_opt & OPT_VIMP) {
      concordanceImputeFlag = FALSE;
      switch (mode) {
      case RSF_GROW:
        obsSize = _observationSize;
        varSize = _xSize;
        statusPtr = _status;
        timePtr = _time;
        ensembleDenPtr = _oobEnsembleDen;
        if (_mRecordSize > 0) {
          concordanceImputeFlag = TRUE;
        }
        break;
      case RSF_PRED:
        obsSize = _fobservationSize;
        varSize = _xSize;
        statusPtr = _fstatus;
        timePtr = _ftime;
        ensembleDenPtr = _fullEnsembleDen;
        if (_fmRecordSize > 0) {
          concordanceImputeFlag = TRUE;
        }
        break;
      case RSF_INTR:
        obsSize = _fobservationSize;
        if ((_opt & OPT_VIMP) && (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN)) {
          varSize = 1;
        }
        else {
          varSize = _intrPredictorSize;
        }
        statusPtr = _fstatus;
        timePtr = _ftime;
        ensembleDenPtr = _oobEnsembleDen;
        if (_fmRecordSize > 0) {
          concordanceImputeFlag = TRUE;
        }
        break;
      default:
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Unknown case in switch encountered. ");
        Rprintf("\nRSF:  Please Contact Technical Support.");
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
        break;
      }
      for (k=1; k <= varSize; k++) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDenPtr[i] != 0) {
            _vimpEnsembleRun[k][i] = _vimpEnsembleRun[k][i] / ensembleDenPtr[i];
          }
        }
        if (getTraceFlag() & DL1_TRACE) {
          Rprintf("\nConcordance Calculation for VIMP (covariate):  %10d", k);  
        }
        if (concordanceImputeFlag == TRUE) {
          imputeConcordance(mode,
                            _forestSize,
                            dmRecordBootFlag,
                            dmvImputation,
                            statusPtr,
                            timePtr);
        }
        concordanceIndex = getConcordanceIndex(obsSize, 
                                               statusPtr,
                                               timePtr,
                                               _vimpEnsembleRun[k], 
                                               ensembleDenPtr);
        if (ISNA(concordanceIndex)) {
          _importance_[k] = NA_REAL;
        }
        else {
          _importance_[k] = 1 - concordanceIndex;
        }
      }
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Variable Importance Measure: \n");
        Rprintf("          ");  
        for (k=1; k <= varSize; k++) {
          Rprintf("%10d", k);
        }
        Rprintf("\n          ");
        for (k=1; k <= varSize; k++) {
          Rprintf("%10.4f", _importance_[k]);
        }
        Rprintf("\n");
      }
    }
    if (mode == RSF_GROW) {
      for (i = 1; i <= _observationSize; i++) {
        for (j=1; j <= sortedTimeInterestSize; j++) {
          if (_oobEnsembleDen[i] != 0) {
            _oobEnsemblePtr[j][i] = _oobEnsemblePtr[j][i] / _oobEnsembleDen[i];
          }
          _fullEnsemblePtr[j][i] = _fullEnsemblePtr[j][i] / _fullEnsembleDen[i];
        }
      }
    }
    if (mode == RSF_PRED) {
      for (i = 1; i <= _fobservationSize; i++) {
        for (j=1; j <= sortedTimeInterestSize; j++) {
          _fullEnsemblePtr[j][i] = _fullEnsemblePtr[j][i] / _fullEnsembleDen[i];
        }
      }
    }
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Ensemble outputs finalized. \n");
    }
  }  
  else {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Insufficient trees for analysis.  \n");
  }
  forestNodeCounter = 0;
  if (mode == RSF_GROW) {
    if (_opt & OPT_TREE) {
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
    if (_opt & OPT_TREE) {
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
  unstackPreDefinedCommonArrays();
  if (mode == RSF_GROW) {
    unstackPreDefinedGrowthArrays(masterSplit,
                                  masterSplitSize);
  }
  else {
    unstackPreDefinedPredictArrays();
  }
  unstackMissingArrays(mode,
                       dmRecordBootFlag,
                       dmvImputation);
  unstackDefinedOutputObjects(mode,
                              sortedTimeInterestSize,
                              root);
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Native code rsf() exit. \n");
    now = updateTimeStamp(start);
  }
  UNPROTECT(stackCount + 2);
  return sexpVector[RSF_OUTP_ID];
}
