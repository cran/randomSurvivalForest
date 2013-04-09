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

#include        "global.h"
#include         "trace.h"
#include        "nrutil.h"
#include        "factor.h"
#include      "rsfStack.h"
#include       "nodeOps.h"
#include       "rsfTree.h"
#include  "rsfBootstrap.h"
#include     "rsfImpute.h"
#include "rsfImportance.h"
#include       "rsfUtil.h"
#include           "rsf.h"
uint stackCount;
char  *sexpString[RSF_SEXP_CNT] = {
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
  "contPT",        
  "mwcpSZ",        
  "mwcpPT",        
  "seed",          
  "importance",    
  "imputation",    
  "oobImputation", 
  "varUsed",       
  "splitDepth",    
  "fullPOE",       
  "oobPOE"         
};
SEXP sexpVector[RSF_SEXP_CNT];
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
double   *_splitDepth_;
double   *_fullEnsemblePOE_;
double   *_oobEnsemblePOE_;
uint     *_treeID_;
uint     *_nodeID_;
uint     *_parmID_;
uint     *_mwcpSZ_;
double   *_contPT_;
uint     *_mwcpPT_;
uint      _opt;
uint      _splitRule;
uint      _splitRandomRule;
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
double   *rsf_ftime;
double   *_fstatus;
double   *_fxData;
uint      _timeInterestSize;
double   *_timeInterest;
SEXP      _sexp_xType;
uint      _intrPredictorSize;
uint     *_intrPredictor;
uint     *_intrIndividual;
uint     *_individualIndex;
uint     *_predictorIndex;
char    **_xType;
double  **_observation;
double  **_fobservation;
uint     *_eventType;
uint      _eventTypeSize;
uint     *_eventTypeIndex;
uint     *_eventTypeIndexSize;
double   *_masterTime;
uint     *_masterTimeIndex;
uint      _masterTimeSize;
char     *_importanceFlag;
uint      _sortedTimeInterestSize;
uint      _factorCount;
uint     *_factorMap;
uint     *_factorIndex;
uint     *_factorSize;
uint      _maxFactorLevel;
Factor  **_factorList;
uint      _mFactorSize;
uint      _fmFactorSize;
uint     *_mFactorIndex;
uint     *_fmFactorIndex;
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
double **_performancePtr;
double ***_oobEnsemblePtr;
double ***_fullEnsemblePtr;
uint     *_oobEnsembleDen;
uint     *_fullEnsembleDen;
uint    **_oobVimpInvalidDen;
double   **_importancePtr;
double   **_vimpMortality;
double ****_crVimpEnsemble;
double  ***_crVimpPOE;
uint      _mStatusSize; 
uint     *_eIndividualSize;
uint     *_meIndividualSize;
uint    **_eIndividual;
double ***_oobSubSurvivalPtr;
double ***_fullSubSurvivalPtr;
double ***_oobSubDistributionPtr;
double ***_fullSubDistributionPtr;
double  **_oobPOEPtr;
double  **_fullPOEPtr;
double  *_sImputeStatusPtr;
double  *_sImputeTimePtr;
double **_sImputePredictorPtr;
double  *_sOOBImputeStatusPtr;
double  *_sOOBImputeTimePtr;
double **_sOOBImputePredictorPtr;
uint   **_varUsedPtr;
double **_splitDepthPtr;
uint     _totalMWCPCount;
int      *_seed1Ptr;
int      *_seed2Ptr;
Node    **_nodeMembership;
uint     *_bootMembershipIndex;
char     *_bootMembershipFlag;
uint     *_oobSampleSize;
char     *_genericMembershipFlag;
Node    **_fnodeMembership;
uint     *_foobSampleSize;
double   _splitValueMaxCont;
uint     _splitValueMaxFactSize;
uint    *_splitValueMaxFactPtr;
clock_t _benchTime;
clock_t _splitTime;
clock_t _hazrdTime;
clock_t _chazrdTime;
clock_t _ensblTime;
clock_t _censblTime;
clock_t _censblTimeSub;
clock_t _vimprTime;
clock_t _cindxTime;
SEXP rsf(char mode, uint traceFlag) {
  uint sexpIndex;
  char       multipleImputeFlag;     
  char       updateFlag;             
  char     **dmRecordBootFlag;  
  double  ***dmvImputation;     
  Node   **root;
  uint     forestNodeCounter;
  double  *localSplitDepth;
  uint *maximumDepthPtr;
  uint maximumDepth;
  Node  *rootPtr;
  Node  *terminalNode;
  uint   rejectedTreeCount;  
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  uint **mwcpPtrPtr;
  uint  *mwcpPtr;
  char  result;
  uint i, j, k, b, r;
  statusPtr    = NULL;  
  timePtr      = NULL;  
  predictorPtr = NULL;  
  mwcpPtrPtr   = NULL;  
  mwcpPtr      = NULL;  
  _benchTime  = 0;
  _splitTime   = 0;
  _hazrdTime   = 0;
  _chazrdTime  = 0;
  _ensblTime   = 0;
  _censblTime  = 0;
  _censblTimeSub  = 0;
  _vimprTime   = 0;
  _cindxTime   = 0;
  setTraceFlag(traceFlag);
  updateTraceFlag(TRUE);
  if (_imputeSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number imputations must be greater than zero:  %10d \n", _forestSize);
    Rprintf("\nRSF:  The application will now exit.\n");
    return R_NilValue;
  }
  if (_forestSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of bootstrap iterations must be greater than zero:  %10d \n", _forestSize);
    Rprintf("\nRSF:  The application will now exit.\n");
    return R_NilValue;
  }
  if (_observationSize < 2) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of individuals must be greater than one:  %10d \n", _observationSize);
    Rprintf("\nRSF:  The application will now exit.\n");
    return R_NilValue;
  }
  if (_xSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of parameters must be greater than zero:  %10d \n", _xSize);
    Rprintf("\nRSF:  The application will now exit.\n");
    return R_NilValue;
  }
  stackPreDefinedCommonArrays();
  switch (mode) {
  case RSF_GROW:
    stackPreDefinedGrowthArrays();
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
    error("\nRSF:  The application will now exit.\n");
    break;
  }
  initializeArrays(mode);
  stackFactorArrays(mode);
  result =  stackMissingArrays(mode,
                               & dmRecordBootFlag,
                               & dmvImputation);
  if (result == FALSE) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Missingness Parameter verification failed.");
    Rprintf("\nRSF:  The application will now exit.\n");
    return R_NilValue;
  }
  result = stackCompetingArrays(mode);
  if (result == FALSE) {
    Rprintf("\nRSF:  The application will now exit.\n");
    return R_NilValue;
  }
  if (_eventTypeSize > 1) {
    _opt = _opt & (~OPT_POUT_TYPE);
  }
  sexpIndex = stackDefinedOutputObjects(mode,
                                        sexpString,
                                        & root,
                                        & _oobEnsemble_,
                                        & _fullEnsemble_,
                                        & _performance_,
                                        & _leafCount_,
                                        & _proximity_,
                                        & _importance_,
                                        & _seed_,
                                        & _oobImputation_,
                                        & _imputation_,
                                        & _sImputeStatusPtr,
                                        & _sImputeTimePtr,
                                        & _sImputePredictorPtr,
                                        & _sOOBImputeStatusPtr,
                                        & _sOOBImputeTimePtr,
                                        & _sOOBImputePredictorPtr,
                                        & _varUsed_,
                                        & _varUsedPtr,
                                        & _splitDepth_,
                                        & _splitDepthPtr,
                                        & localSplitDepth,
                                        & _oobEnsemblePOE_,
                                        & _fullEnsemblePOE_,
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
  *_seed2Ptr = -abs(*_seed2Ptr) * 251;
  if (mode == RSF_GROW) {
    if (_opt & OPT_TREE) {
      _totalMWCPCount = 0;
    }
  }
  else {
    forestNodeCounter = 1;
    mwcpPtr = _mwcpPT_;
    mwcpPtrPtr = & mwcpPtr;
  }
  multipleImputeFlag = FALSE;
  for (r = 1; r <= _imputeSize; r++) {
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
      if (TRUE) {
        unImpute (mode);
      }
      if (mode == RSF_GROW) {
        if (_imputeSize > 1) {
          if (r > 1) {
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
            if (_mTimeIndexFlag == TRUE) {
              updateTimeIndexArray(NULL);
            }
          }
        }
      }
      rootPtr = makeNode(_xSize);  
      rootPtr -> parent = NULL;
      rootPtr -> leafCount          = 1;
      for (i = 1; i <= _xSize; i++) {
        (rootPtr -> permissibleSplit)[i] = TRUE;
      }
      if ((_opt & OPT_TREE) && (r == _imputeSize)) {
        root[b] = rootPtr;
      }
      maximumDepth = 0;
      maximumDepthPtr = &maximumDepth;
      if (mode == RSF_GROW) {
        result = bootstrap (mode,
                            b,
                            rootPtr,
                            dmRecordBootFlag); 
        if (result) {
          _leafCount_[b] = 1;
          makeTree (multipleImputeFlag, b, rootPtr, 0, maximumDepthPtr);
          if ((_imputeSize > 1) && (r > 1) && (r < _imputeSize) ) {
            if (_mRecordSize > 0) {
              for (j = 1; j <= _leafCount_[b]; j++) {
                imputeNode(mode, FALSE, b, getTerminalNode(mode, j));
              }  
            }  
          }  
        }  
        else {
        }
      }  
      else {
        result = bootstrap (mode,
                            b,
                            rootPtr,
                            dmRecordBootFlag); 
        if (result) {
          _leafCount_[b] = 1;
          restoreTree(b,
                      rootPtr,
                      _leafCount_ + b,
                      & forestNodeCounter,
                      _treeID_,
                      _nodeID_,
                      _parmID_,
                      _contPT_,
                      _mwcpSZ_,
                      mwcpPtrPtr,
                      0,
                      maximumDepthPtr);
          imputeTree(mode, b, rootPtr, TRUE);
          if (_mTimeIndexFlag == TRUE) {
          }
        }  
        else {
          forestNodeCounter ++;
        }
      }  
      if (result) {
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
          }  
          if (_opt & OPT_SPLT_DPTH) {
            for (i = 1; i <= _observationSize; i++) {
              if (TRUE) {
                for (j = 1; j <= _xSize; j++) {
                  localSplitDepth[j] = NA_REAL;
                }
                terminalNode = getProxyMember(rootPtr, _observation, i);
                for (k = 1; k <= terminalNode -> depth; k++) {
                  if (ISNA(localSplitDepth[(terminalNode -> splitDepth)[k]])) {
                     localSplitDepth[(terminalNode -> splitDepth)[k]] = (double) k;
                  }
                }
                for (j = 1; j <= _xSize; j++) {
                  if (ISNA(localSplitDepth[j])) {
                    localSplitDepth[j] = (double) *maximumDepthPtr + 1;
                  }
                }
                for (j = 1; j <= _xSize; j++) {
                  _splitDepthPtr[j][i] += localSplitDepth[j];
                }
              }
            }
          }
        }  
        updateFlag = FALSE;
        switch (mode) {
        case RSF_PRED:
          if (_fmRecordSize > 0) {
            updateFlag = TRUE;
            statusPtr = _fstatus;
            timePtr = rsf_ftime;
            predictorPtr = _fobservation;
          }
          break;
        case RSF_INTR:
          if (_fmRecordSize > 0) {
            updateFlag = TRUE;
            statusPtr = _fstatus;
            timePtr = rsf_ftime;
            predictorPtr = _fobservation;
          }
          break;
        default:
          if ((r == 1) || (r < _imputeSize)) {
            if (_mRecordSize > 0) {
              updateFlag = TRUE;
              statusPtr = _status;
              timePtr = _time;
              predictorPtr = _observation;
            }
          }
          break;
        }
        if (updateFlag == TRUE) {
          imputeUpdateSummary(mode, 
                              statusPtr, 
                              timePtr, 
                              predictorPtr,
                              b, 
                              dmvImputation);
        }
        if (!(_opt & OPT_IMPU_ONLY)) {
          if (r == _imputeSize) {
            updateEnsembleCalculations(multipleImputeFlag,
                                       mode,
                                       rootPtr,
                                       b,
                                       dmRecordBootFlag,
                                       dmvImputation);
          }
        }
      }  
      if (!(_opt & OPT_TREE) || !(r == _imputeSize)) {
        freeTree(rootPtr);
      }
      updateTraceFlag(FALSE);
    }  
    updateTraceFlag(TRUE);
    if (mode == RSF_GROW) {
      if (r == 1) {
        if (_mRecordSize > 0) {
          imputeSummary(mode,
                        FALSE,
                        dmRecordBootFlag,
                        dmvImputation);
          imputeMultipleTime(FALSE);
        }
      }  
      else {
        if (r < _imputeSize) {
          if (_mRecordSize > 0) {
            imputeSummary(mode,
                          ACTIVE,
                          dmRecordBootFlag,
                          dmvImputation);
            imputeMultipleTime(ACTIVE);
          }
        }
      }  
    } 
  }  
  if (FALSE) {
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
  switch (mode) {
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
    if (_imputeSize == 1) {
      if (_mRecordSize > 0) {
        imputeSummary(RSF_GROW,
                      TRUE,
                      dmRecordBootFlag,
                      dmvImputation);
      }
    }  
    else {
    }  
    break;
  }
  if (rejectedTreeCount < _forestSize) {
    if (_opt & OPT_VIMP) {
      finalizeVariableImportance(mode,
                                 rejectedTreeCount,
                                 dmRecordBootFlag,
                                 dmvImputation);
    }
    if (mode == RSF_GROW) {
      if (_opt & OPT_OENS) {
        for (i = 1; i <= _observationSize; i++) {
          if (_opt & (OPT_POUT_TYPE)) {
            if (_oobEnsembleDen[i] != 0) {
              _oobEnsemblePtr[1][1][i] = _oobEnsemblePtr[1][1][i] / _oobEnsembleDen[i];
            }
          }
          else {
            for (k=1; k <= _sortedTimeInterestSize; k++) {
              if (_oobEnsembleDen[i] != 0) {
                _oobEnsemblePtr[1][k][i] = _oobEnsemblePtr[1][k][i] / _oobEnsembleDen[i];
              }
            }
            if (_eventTypeSize > 1) {
              if (_oobEnsembleDen[i] != 0) {
                for (j = 1; j <= _eventTypeSize; j++) {
                  _oobPOEPtr[j][i] = _oobPOEPtr[j][i] / _oobEnsembleDen[i];
                }
              }
            }
          }
        }
      }  
      if (_opt & OPT_FENS) {
        for (i = 1; i <= _observationSize; i++) {
          if (_opt & (OPT_POUT_TYPE)) {
            _fullEnsemblePtr[1][1][i] = _fullEnsemblePtr[1][1][i] / _fullEnsembleDen[i];
          }
          else {
            for (k=1; k <= _sortedTimeInterestSize; k++) {
              _fullEnsemblePtr[1][k][i] = _fullEnsemblePtr[1][k][i] / _fullEnsembleDen[i];
            }
            if (_eventTypeSize > 1) {
              for (j = 1; j <= _eventTypeSize; j++) {
                _fullPOEPtr[j][i] = _fullPOEPtr[j][i] / _fullEnsembleDen[i];
              }
            }
          }
        }
      }
    }  
    if (mode == RSF_PRED) {
      if (_opt & OPT_FENS) {
        for (i = 1; i <= _fobservationSize; i++) {
          if (_opt & (OPT_POUT_TYPE)) {
            _fullEnsemblePtr[1][1][i] = _fullEnsemblePtr[1][1][i] / _fullEnsembleDen[i];
          }
          else {
            for (k=1; k <= _sortedTimeInterestSize; k++) {
              _fullEnsemblePtr[1][k][i] = _fullEnsemblePtr[1][k][i] / _fullEnsembleDen[i];
            }
            if (_eventTypeSize > 1) {
              for (j = 1; j <= _eventTypeSize; j++) {
                _fullPOEPtr[j][i] = _fullPOEPtr[j][i] / _fullEnsembleDen[i];
              }
            }
          }
        }
      }
    }
    if (_opt & OPT_SPLT_DPTH) {
      for (j = 1; j <= _xSize; j++) {
        for (i = 1; i <= _observationSize; i++) {
          _splitDepthPtr[j][i] = _splitDepthPtr[j][i] / (_forestSize - rejectedTreeCount);
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
                               _totalMWCPCount,    
                               & _treeID_,         
                               & _nodeID_,         
                               & _parmID_,         
                               & _contPT_,         
                               & _mwcpSZ_,         
                               & _mwcpPT_,         
                               sexpIndex, 
                               sexpString,
                               sexpVector);
  if (mode == RSF_GROW) {
    if (_opt & OPT_TREE) {
      mwcpPtr = _mwcpPT_;
      mwcpPtrPtr = & mwcpPtr;
      forestNodeCounter = 1;
      for (b = 1; b <= _forestSize; b++) {
        saveTree(b, 
                 root[b], 
                 & forestNodeCounter, 
                 _treeID_, 
                 _nodeID_, 
                 _parmID_, 
                 _contPT_,
                 _mwcpSZ_,
                 mwcpPtrPtr);
      }
      forestNodeCounter --;
    }  
  }  
  unstackDefinedOutputObjects(mode,
                              root,
                              localSplitDepth);
  unstackCompetingArrays(mode);
  unstackMissingArrays(mode,
                       dmRecordBootFlag,
                       dmvImputation);
  unstackFactorArrays(mode);
  switch (mode) {
  case RSF_GROW:
    unstackPreDefinedGrowthArrays();
    break;
  case RSF_PRED:
    unstackPreDefinedPredictArrays();
    break;
  case RSF_INTR:
    unstackPreDefinedInteractionArrays();
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    error("\nRSF:  The application will now exit.\n");
    break;
  }
  unstackPreDefinedCommonArrays();
  UNPROTECT(stackCount + 2);
  return sexpVector[RSF_OUTP_ID];
}
