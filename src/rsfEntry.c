//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.5.1
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

#include   "global.h"
#include   "nrutil.h"
#include   "rsfEntry.h"
extern uint getTraceFlag();
extern SEXP rsf(uint mode, uint traceFlag);
SEXP rsfGrow(SEXP traceFlag,
             SEXP opt,  
             SEXP seedPtr,  
             SEXP splitRule,  
             SEXP splitRandomRule,  
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
  _splitRandomRule      = INTEGER(splitRandomRule)[0];
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
       (_splitRule != RANDOM_SPLIT) ) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Invalid split rule:  %10d \n", _splitRule);
    return R_NilValue;
  }
  if (_splitRule == RANDOM_SPLIT) {
    _splitRandomRule = 1;
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
  _contPT_ = NULL;
  _mwcpSZ_ = NULL;
  _mwcpPT_ = NULL;
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
                SEXP contPT,
                SEXP mwcpSZ,
                SEXP mwcpPT,
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
  _contPT_              =    REAL(contPT);  _contPT_ --;
  _mwcpSZ_              = (uint*) INTEGER(mwcpSZ);  _mwcpSZ_ --;
  _mwcpPT_              = (uint*) INTEGER(mwcpPT);  _mwcpPT_ --;
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
                    SEXP contPT,
                    SEXP mwcpSZ,
                    SEXP mwcpPT,
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
  _contPT_              =    REAL(contPT);  _contPT_ --;
  _mwcpSZ_              = (uint*) INTEGER(mwcpSZ);  _mwcpSZ_ --;
  _mwcpPT_              = (uint*) INTEGER(mwcpPT);  _mwcpPT_ --;
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
