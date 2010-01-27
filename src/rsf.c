////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.1
////
////  Copyright 2009, Cleveland Clinic Foundation
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
////  ----------------------------------------------------------------
////  Project Partially Funded By:
////    --------------------------------------------------------------
////    National Institutes of Health,  Grant HHSN268200800026C/0001
////
////    Michael S. Lauer, M.D., FACC, FAHA 
////    National Heart, Lung, and Blood Institute
////    6701 Rockledge Dr, Room 10122
////    Bethesda, MD 20892
////
////    email:  lauerm@nhlbi.nih.gov
////
////    --------------------------------------------------------------
////    Case Western Reserve University/Cleveland Clinic  
////    CTSA Grant:  UL1 RR024989, National Center for
////    Research Resources (NCRR), NIH
////
////    --------------------------------------------------------------
////    Dept of Defense Era of Hope Scholar Award, Grant W81XWH0910339
////    Andy Minn, M.D., Ph.D.
////    Department of Radiation and Cellular Oncology, and
////    Ludwig Center for Metastasis Research
////    The University of Chicago, Jules F. Knapp Center, 
////    924 East 57th Street, Room R318
////    Chicago, IL 60637
//// 
////    email:  aminn@radonc.uchicago.edu
////
////    --------------------------------------------------------------
////    Bryan Lau, Ph.D.
////    Department of Medicine, Johns Hopkins School of Medicine,
////    Baltimore, Maryland 21287
////
////    email:  blau1@jhmi.edu
////
////  ----------------------------------------------------------------
////  Written by:
////    --------------------------------------------------------------
////    Hemant Ishwaran, Ph.D.
////    Dept of Quantitative Health Sciences/Wb4
////    Cleveland Clinic Foundation
////    9500 Euclid Avenue
////    Cleveland, OH 44195
////
////    email:  hemant.ishwaran@gmail.com
////    phone:  216-444-9932
////    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
////
////    --------------------------------------------------------------
////    Udaya B. Kogalur, Ph.D.
////    Dept of Quantitative Health Sciences/Wb4
////    Cleveland Clinic Foundation
////    
////    Kogalur Shear Corporation
////    5425 Nestleway Drive, Suite L1
////    Clemmons, NC 27012
////
////    email:  ubk2101@columbia.edu
////    phone:  919-824-9825
////    URL:    www.kogalur-shear.com
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
  uint   obsSize;
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  uint **mwcpPtrPtr;
  uint  *mwcpPtr;
  char  result;
  uint i, j, k, b, r;
  clock_t totalTime, startTime;
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
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\n\nRSF:  Native code rsf() entry:  bld20100127 \n");
  }
  startTime = clock();
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
  if (getTraceFlag() & SUMM_USR_TRACE) {
    switch (mode) {
    case RSF_GROW:
      Rprintf("\nRSF:  Mode is GROW.");
      Rprintf("\nRSF:  Split rule is:               %10d", _splitRule);
      Rprintf("\nRSF:  Split random rule is:        %10d", _splitRandomRule);
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
  if (getTraceFlag() & SUMM_USR_TRACE) {
    int endian = 1;
    if (*(char *) &endian == 1) {
      Rprintf("\nRSF:  System is little-endian.  ");
    }
    else {
      Rprintf("\nRSF:  System is big-endian.  ");
    }
    Rprintf("\nRSF:  Size of (int)    is:  %10d   ", sizeof(int));
    Rprintf("\nRSF:  Size of (double) is:  %10d \n", sizeof(double));
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
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
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
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nStart of random seed chain ran1():  %20d", *_seed1Ptr);
    Rprintf("\nStart of random seed chain ran2():  %20d", *_seed2Ptr);
  }
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
    if (getTraceFlag() & SUMM_LOW_TRACE) {
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
      if (getTraceFlag() & SUMM_LOW_TRACE) {
        Rprintf("\nStart of iteration:  (%10d, %10d)", r, b);
      }
      if (getTraceFlag() & SUMM_MED_TRACE) {
        Rprintf("\nTree random seed ran1():  %20d", *_seed1Ptr);
        Rprintf("\nTree random seed ran2():  %20d", *_seed2Ptr);
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
          if (getTraceFlag() & SUMM_LOW_TRACE) {
            Rprintf("\nStart of makeTree():  (MI iter %10d, tree %10d)", r, b);
          }
          if (getTraceFlag() & TIME_DEF_TRACE) {
            _benchTime = clock();
          }
          makeTree (multipleImputeFlag, b, rootPtr, 0, maximumDepthPtr);
          if (getTraceFlag() & TIME_DEF_TRACE) {
            _splitTime = _splitTime + (clock() - _benchTime);
          }
          if (getTraceFlag() & SUMM_LOW_TRACE) {
            Rprintf("\nEnd of makeTree():  (MI iter %10d, tree %10d)", r, b);
          }
          if (getTraceFlag() & SUMM_USR_TRACE) {
            Rprintf("\nRSF:  Tree construction complete:  (MI iter %10d, tree %10d)", r, b);  
            Rprintf("\nRSF:  Final leaf count:  %10d\n", _leafCount_[b]);
          }
          if ((_imputeSize > 1) && (r > 1) && (r < _imputeSize) ) {
            if (_mRecordSize > 0) {
              for (j = 1; j <= _leafCount_[b]; j++) {
                imputeNode(mode, FALSE, b, getTerminalNode(mode, j));
              }  
            }  
          }  
        }  
        else {
          if (getTraceFlag() & SUMM_USR_TRACE) {
            Rprintf("\nRSF:  Tree rejected:  %10d", b);  
          }
        }
      }  
      else {
        result = bootstrap (mode,
                            b,
                            rootPtr,
                            dmRecordBootFlag); 
        if (result) {
          _leafCount_[b] = 1;
          if ((mode == RSF_PRED) && (getTraceFlag() & SUMM_MED_TRACE)) {
            Rprintf("\nIncoming Tree:  ");
            Rprintf("\n      tree       parm       node         splt \n");
          }
          if (getTraceFlag() & SUMM_LOW_TRACE) {
            Rprintf("\nStart of restore process:  (%10d, %10d)", r, b);
          }
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
          if (getTraceFlag() & SUMM_LOW_TRACE) {
            Rprintf("\nEnd of restore process:  (%10d, %10d)", r, b);
          }
          if (getTraceFlag() & SUMM_USR_TRACE) {
            Rprintf("\nRSF:  Tree construction complete:  (%10d, %10d)", r, b);  
            Rprintf("\nRSF:  Final leaf count:  %10d\n", _leafCount_[b]);
          }
          imputeTree(mode, b, rootPtr, TRUE);
          if (_mTimeIndexFlag == TRUE) {
          }
        }  
        else {
          forestNodeCounter ++;
          if (getTraceFlag() & SUMM_USR_TRACE) {
            Rprintf("\nRSF:  Tree rejected:  %10d", b);  
          }
        }
      }  
      if (result) {
        if (getTraceFlag() & SUMM_MED_TRACE) {
          Rprintf("\nFinal TRAINING Membership (all data, in-bag and OOB):  %10d", b);
          Rprintf("\n     index       leaf  bootIndex\n");
          for (i=1; i <=  _observationSize; i++) {
            Rprintf("%10d %10d %10d \n", i, _nodeMembership[i] -> leafCount, _bootMembershipIndex[i]);
          }
          uint  *tnmCount   = uivector(1, _leafCount_[b]);
          uint **eventCount = uimatrix(1, _eventTypeSize + 1, 1, _leafCount_[b]);
          Rprintf("\nFinal TRAINING Diagnostics (in-bag data):  %10d", b);
          Rprintf("\n      leaf   tnmCount totEvCount");
          for (i=1; i <= _leafCount_[b]; i++) {
            tnmCount[i] = 0;
            eventCount[1][i] = 0;
          }
          Rprintf("\n");
          for (i=1; i <= _observationSize; i++) {
            tnmCount[_nodeMembership[_bootMembershipIndex[i]] -> leafCount] ++;
            if (_status[_bootMembershipIndex[i]] > 0) {
              eventCount[1][_nodeMembership[_bootMembershipIndex[i]] -> leafCount] ++;
            }
          }
          for (i=1; i <= _leafCount_[b]; i++) {
            Rprintf("%10d %10d %10d", i, tnmCount[i], eventCount[1][i]);
            Rprintf("\n");
          }
          free_uivector(tnmCount, 1, _leafCount_[b]);
          free_uimatrix(eventCount, 1, _eventTypeSize + 1, 1, _leafCount_[b]);
        }
        if ((mode == RSF_PRED) && (getTraceFlag() & SUMM_MED_TRACE)) {
          Rprintf("\nFinal PRED Membership (all data):  %10d", b);
          Rprintf("\n     index       leaf\n");
          for (i=1; i <=  _fobservationSize; i++) {
            Rprintf("%10d %10d \n", i, _fnodeMembership[i] -> leafCount);
          }
        }
        if ((mode == RSF_INTR) && (getTraceFlag() & SUMM_MED_TRACE)) {
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
            if (getTraceFlag() & SUMM_USR_TRACE) {
              Rprintf("\nRSF:  Proximity calculation complete.");  
            }
          }  
          if (_opt & OPT_SPLT_DPTH) {
            if (getTraceFlag() & SUMM_MED_TRACE) {
              Rprintf("\nTerminal Node Predictor Split Depths for Tree:  %10d", b);
              Rprintf("\n      leaf      count     splits --> \n");
              for (k=1; k <= _leafCount_[b]; k++) {
                terminalNode = getTerminalNode(mode, k);
                if (terminalNode != NULL) {
                  Rprintf("%10d %10d ", k, terminalNode -> depth);
                  for (j = 1; j <= terminalNode -> depth; j++) {
                    Rprintf("%10d ", (terminalNode -> splitDepth)[j]);
                  }
                }
                else {
                  Rprintf("%10d ", k);
                }
                Rprintf("\n");
              }
            }
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
        if (getTraceFlag() & FORK_DEF_TRACE) {
          Rprintf("\nDe-allocating treeID:  %10d", b);          
        }
        freeTree(rootPtr);
      }
      if (getTraceFlag() & TIME_DEF_TRACE) {
        totalTime = clock() - startTime;
      }
      if (getTraceFlag() & SUMM_LOW_TRACE) {
        Rprintf("\nEnd of iteration:    (%10d, %10d)", r, b);
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
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nEnd of impute iteration:  %10d", r);
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
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nFinal Tree Leaf Counts:  ");
    Rprintf("\n         tree    leafCount");
    for (b = 1; b <= _forestSize; b++) {
      Rprintf("\n %12d %12d", b, _leafCount_[b]);
    }
    Rprintf("\n");
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Trees rejected:  %10d ", rejectedTreeCount);
    Rprintf("\nRSF:  Trees stumped:   %10d ", k);
    Rprintf("\nRSF:  Trees (total):   %10d \n", _forestSize);
  }
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    if (_opt & OPT_PROX) {
      switch (mode) {
      case RSF_GROW:
        obsSize = _observationSize;
      case RSF_PRED:
        obsSize = _fobservationSize;
        break;
      default:
        obsSize = 0;
      }
      Rprintf("\nProximity Matrix:  \n");
      k = 0;
      for (i = 1; i <= obsSize; i++) {
        k += i - 1;
        for (j = 1; j <= i; j++) {
          Rprintf("%10d ", _proximity_[k + j]);
        }
        Rprintf("\n");
      }
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
    if (getTraceFlag() & SUMM_USR_TRACE) {
      if ((_opt & OPT_OENS) || (_opt & OPT_FENS)) {
        Rprintf("\nRSF:  Ensemble outputs finalized. \n");
      }
    }
    if (getTraceFlag() & ENSB_LOW_TRACE) {
      switch (mode) {
      case RSF_GROW:
        obsSize = _observationSize;
        break;
      case RSF_PRED:
        obsSize = _fobservationSize;
        break;
      case RSF_INTR:
        obsSize = _fobservationSize;
        break;
      default:
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Unknown case in switch encountered. ");
        Rprintf("\nRSF:  Please Contact Technical Support.");
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
        break;
      }
      if (_opt & OPT_OENS) {
        for (j = 1; j <= 1; j++) {
          Rprintf("\nOOB Ensemble CHF:  \n");
          Rprintf("          ");
          for (i=1; i <= obsSize; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n");
          for (k=1; k <= _sortedTimeInterestSize; k++) {
            Rprintf("%10d", k);
            for (i=1; i <= obsSize; i++) {
              Rprintf("%10.4f", _oobEnsemblePtr[1][k][i]);
            }
            Rprintf("\n");
          }
        }
        if (_eventTypeSize > 1) {
          for (j = 1; j <= _eventTypeSize; j++) {
            Rprintf("\nOOB Ensemble Conditional CHF:  Event %10d \n", j);
            Rprintf("          ");
            for (i=1; i <= obsSize; i++) {
              Rprintf("%10d", i);
            }
            Rprintf("\n");
            for (k=1; k <= _sortedTimeInterestSize; k++) {
              Rprintf("%10d", k);
              for (i=1; i <= obsSize; i++) {
                Rprintf("%10.4f", _oobEnsemblePtr[j+1][k][i]);
              }
              Rprintf("\n");
            }
          }
          Rprintf("\nOOB Ensemble POE:  \n");
          Rprintf("          ");
          for (i=1; i <= obsSize; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n");
          for (j = 1; j <= _eventTypeSize; j++) {
            Rprintf("%10d", j);
            for (i=1; i <= obsSize; i++) {
              Rprintf("%10.4f", _oobPOEPtr[j][i]);
            }
            Rprintf("\n");
          }
        }
      }
      if (_opt & OPT_FENS) {
        for (j = 1; j <= 1; j++) {
          Rprintf("\nFULL Ensemble CHF:  \n");
          Rprintf("          ");
          for (i=1; i <= obsSize; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n");
          for (k=1; k <= _sortedTimeInterestSize; k++) {
            Rprintf("%10d", k);
            for (i=1; i <= obsSize; i++) {
              Rprintf("%10.4f", _fullEnsemblePtr[1][k][i]);
            }
            Rprintf("\n");
          }
        }
        if (_eventTypeSize > 1) {
          for (j = 1; j <= _eventTypeSize; j++) {
            Rprintf("\nFULL Ensemble Conditional CHF:  Event %10d \n", j);
            Rprintf("          ");
            for (i=1; i <= obsSize; i++) {
              Rprintf("%10d", i);
            }
            Rprintf("\n");
            for (k=1; k <= _sortedTimeInterestSize; k++) {
              Rprintf("%10d", k);
            for (i=1; i <= obsSize; i++) {
                Rprintf("%10.4f", _fullEnsemblePtr[j+1][k][i]);
              }
              Rprintf("\n");
            }
          }
          Rprintf("\nFULL Ensemble POE:  \n");
          Rprintf("          ");
          for (i=1; i <= obsSize; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n");
          for (j = 1; j <= _eventTypeSize; j++) {
            Rprintf("%10d", j);
            for (i=1; i <= obsSize; i++) {
              Rprintf("%10.4f", _fullPOEPtr[j][i]);
            }
            Rprintf("\n");
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
      if (getTraceFlag() & SUMM_MED_TRACE) {
        Rprintf("\nRSF:  Predictor Split Depths by Individual: \n");
        Rprintf("\n   predictor   individual -->");
        Rprintf("\n             ");
        for (i = 1; i <= _observationSize; i++) {
          Rprintf("%12d ", i);
        }
        Rprintf("\n");
        for (j=1; j <= _xSize; j++) {
          Rprintf("%12d ", j);
          for (i = 1; i <= _observationSize; i++) {
            Rprintf("%12.4f ", _splitDepthPtr[j][i]);
          }
          Rprintf("\n");
        }
      }
      if (getTraceFlag() & SUMM_USR_TRACE) {
        Rprintf("\nRSF:  Predictor split depth output finalized. \n");
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
      if (getTraceFlag() & SUMM_HGH_TRACE) {
        Rprintf("\nGrown Forest Node Count:  %10d", forestNodeCounter);
        Rprintf("\nGrown Forest Output:  ");
        Rprintf("\n    treeID     nodeID     parmID       contPT     mwcpSZ \n");
        for (i = 1; i <= forestNodeCounter; i++) {
          Rprintf("%10d %10d %10d %12.4f %10d \n", _treeID_[i], _nodeID_[i], _parmID_[i], _contPT_[i], _mwcpSZ_[i]);
        }
        Rprintf("\n");
        Rprintf("\n     index       mwcpPT \n");
        for (i = 1; i <= _totalMWCPCount; i++) {
          Rprintf("%10d %12d \n",i, _mwcpPT_[i]);
        }
        Rprintf("\n");
        Rprintf("\nGrown Forest Membership by Observation:  ");
        Rprintf("\n    OBSERV    TREE_ID    NODE_ID");
        for (b = 1; b <= _forestSize; b++) {
          for (k = 1; k <= _observationSize; k++) {
            Rprintf("\n%10d %10d %10d ", k, b, getProxyMember(root[b], _observation, k) -> leafCount);
          }
        }
        Rprintf("\n");
      }
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
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  unstackPreDefinedCommonArrays();
  if (getTraceFlag() & TIME_DEF_TRACE) {
    totalTime = clock() - startTime;
    Rprintf("\n  _splitTime  = %20d ",  _splitTime);
    Rprintf("\n  _hazrdTime  = %20d ",  _hazrdTime);
    Rprintf("\n _chazrdTime  = %20d ",  _chazrdTime);
    Rprintf("\n  _ensblTime  = %20d ",  _ensblTime);
    Rprintf("\n _censblTime  = %20d ",  _censblTime);
    Rprintf("\n _censblTimeSub  = %20d ",  _censblTimeSub);
    Rprintf("\n  _vimprTime  = %20d ",  _vimprTime);
    Rprintf("\n  _cindxTime  = %20d ",  _cindxTime);
    Rprintf("\n");
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Native code rsf() exit. \n");
  }
  UNPROTECT(stackCount + 2);
  return sexpVector[RSF_OUTP_ID];
}
