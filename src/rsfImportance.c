////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.2
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
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include  "rsfFactorOps.h"
#include     "rsfImpute.h"
#include       "rsfUtil.h"
#include "rsfImportance.h"
Node *getProxyMember(Node    *parent, 
                     double **predictor, 
                     uint     index) {
  char daughterFlag;
  uint i;
  Node *result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    daughterFlag = RIGHT;
    if (strcmp(_xType[parent -> splitParameter], "C") == 0) {
      daughterFlag = splitOnFactor((uint) predictor[parent -> splitParameter][index], parent -> splitValueFactPtr);
    }
    else {
      if (predictor[parent -> splitParameter][index] <= (parent -> splitValueCont)) {
        daughterFlag = LEFT;
      }
    }
    if (daughterFlag == LEFT) {
      result = getProxyMember(parent ->  left, predictor, index);
    }
    else {
      result = getProxyMember(parent -> right, predictor, index);
    }
  }
  return result;
}
Node *randomizeMembership(Node    *parent, 
                          double **predictor, 
                          uint     individual, 
                          uint     splitParameter) {
  char daughterFlag;
  char randomSplitFlag;
  Node *result;
  result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    randomSplitFlag = FALSE;
    if (splitParameter > 0) {
      if ((parent -> splitParameter) == splitParameter) {
        randomSplitFlag = TRUE;
      }
    }
    else {
      if(_importanceFlag[parent -> splitParameter] == TRUE) {
        randomSplitFlag = TRUE;
      }
    }
    if(randomSplitFlag == TRUE) {
      if (ran2(_seed2Ptr) <= 0.5) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter);
      }
    }
    else {
      daughterFlag = RIGHT;
      if (strcmp(_xType[parent -> splitParameter], "C") == 0) {
        daughterFlag = splitOnFactor((uint) predictor[parent -> splitParameter][individual], parent -> splitValueFactPtr);
      }
      else {
        if (predictor[parent -> splitParameter][individual] <= (parent -> splitValueCont)) {
          daughterFlag = LEFT;
        }
      }
      if (daughterFlag == LEFT) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter);
      }
    }
  }
  return result;
}
void permute(uint n, uint *indx) {
  uint i,j,k;
  for (i=1; i<= n; i++) {
    indx[i] = 0;
  }
  for (i=n; i > 0; i--) {
    k = (uint) ceil(ran2(_seed2Ptr)*(i*1.0));
    for (j = 1; k > 0; j++) {
      if (indx[j] == 0) {
        k--;
      }
    }
    indx[j-1] = i;
  }
}
void getVariableImportance (uint      mode,
                            uint      leafCount,
                            Node     *rootPtr,
                            uint      b) {
  uint j, k, n, p;
  uint obsSize;
  uint varCount;
  double **predictorPtr;
  char selectionFlag;
  char result;
  if (!(_opt & OPT_VIMP)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to compute variable importance though not requested.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  selectionFlag = ACTIVE;  
  result = TRUE;
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    varCount = _xSize;
    predictorPtr = _observation;
    selectionFlag = FALSE;
    if (_oobSampleSize[b] > 0) {
      result = TRUE;
    }
    break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    varCount = _xSize;
    predictorPtr = _fobservation;
    result = TRUE;
    selectionFlag = ACTIVE;
    break;
  case RSF_INTR:
    obsSize = _fobservationSize;
    varCount = _intrPredictorSize;
    predictorPtr = _fobservation;
    selectionFlag = FALSE;
    if (_foobSampleSize[b] > 0) {
      result = TRUE;
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
  if (result == TRUE) {    
    if (_opt & (~OPT_VIMP) & (~OPT_VIMP_JOIN) & OPT_VIMP_TYPE) {
      getVimpRandom(mode, rootPtr, predictorPtr, obsSize, varCount, selectionFlag);
    }
    else if (!(_opt & (~OPT_VIMP) & (~OPT_VIMP_JOIN) & OPT_VIMP_TYPE)) {
      getVimpPermute(mode, rootPtr, predictorPtr, b, obsSize, varCount, selectionFlag);
    }
    else {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Unknown VIMP perturbation type encountered. ");
      Rprintf("\nRSF:  Option flag is:  %10x", _opt);
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
    } 
  }  
  else {
  }
}
void getVimpRandom (uint      mode,
                    Node     *rootPtr,
                    double  **predictorPtr,
                    uint      obsSize,
                    uint      varCount,
                    char      selectionFlag) {
  Node *terminalNode;
  uint i, j, k, p;
  if (!(_opt & (~OPT_VIMP) & OPT_VIMP_JOIN)) {
    for (p=1; p <= varCount; p++) {
      for (i=1; i <= obsSize; i++) {
        if ((_genericMembershipFlag[_individualIndex[i]] == selectionFlag) || (selectionFlag == ACTIVE)) {
          terminalNode = randomizeMembership(rootPtr, predictorPtr, i, _predictorIndex[p]);
          if (!ISNA(terminalNode -> mortality)) {
            _vimpMortality[p][i] += terminalNode -> mortality;
            if (_eventTypeSize > 1) {
              for (j=1; j<= _eventTypeSize; j++) {
                _crVimpPOE[p][j][i] += (double) (terminalNode -> poe)[j] / (terminalNode -> eventCount);
                for (k=1; k <= _sortedTimeInterestSize; k++) {
                  _crVimpEnsemble[p][j][k][i] += terminalNode -> subSurvival[j][k];
                }
              }
            }
          }
          else {
            if (_opt & OPT_VOUT_TYPE) {
              _oobVimpInvalidDen[p][i] ++;
            }
            else {
              Rprintf("\nRSF:  *** ERROR *** ");
              Rprintf("\nRSF:  NA encountered for mortality in VIMP.");
              Rprintf("\nRSF:  Please Contact Technical Support.");
              Rprintf("\nRSF:  The application will now exit.\n");
              exit(TRUE);
            }
          }
        }
      }
    }  
  }
  else {
    for (i=1; i <= obsSize; i++) {
      if ((_genericMembershipFlag[_individualIndex[i]] == selectionFlag) || (selectionFlag == ACTIVE)) {
        terminalNode = randomizeMembership(rootPtr, predictorPtr, i, 0);
        if (!ISNA(terminalNode -> mortality)) {
          _vimpMortality[1][i] += terminalNode -> mortality;
          if (_eventTypeSize > 1) {
            for (j=1; j<= _eventTypeSize; j++) {
              _crVimpPOE[1][j][i] += (double) (terminalNode -> poe)[j] / (terminalNode -> eventCount);
              for (k=1; k <= _sortedTimeInterestSize; k++) {
                _crVimpEnsemble[1][j][k][i] += terminalNode -> subSurvival[j][k];
              }
            }
          }
        }
        else {
          if (_opt & OPT_VOUT_TYPE) {
            _oobVimpInvalidDen[1][i] ++;
          }
          else {
            Rprintf("\nRSF:  *** ERROR *** ");
            Rprintf("\nRSF:  NA encountered for mortality in VIMP.");
            Rprintf("\nRSF:  Please Contact Technical Support.");
            Rprintf("\nRSF:  The application will now exit.\n");
            exit(TRUE);
          }
        }
      }
    }
  }
}
void getVimpPermute(uint      mode,
                    Node     *rootPtr,
                    double  **predictorPtr,
                    uint      b,
                    uint      obsSize,
                    uint      varCount,
                    char      selectionFlag) {
  Node *terminalNode;
  uint permuteObsSize;
  uint    *indexVIMP;
  uint    *permuteVIMP;
  uint i, j, k, p;
  switch (mode) {
  case RSF_GROW:
    permuteObsSize = _oobSampleSize[b];
    break;
  case RSF_PRED:
    permuteObsSize = _fobservationSize;
    break;
  case RSF_INTR:
    permuteObsSize = _foobSampleSize[b];
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  indexVIMP = uivector(1, permuteObsSize);
  permuteVIMP = uivector(1, permuteObsSize);
  k = 0;
  for (i=1; i <= obsSize; i++) {
    if ((_genericMembershipFlag[_individualIndex[i]] == selectionFlag) || (selectionFlag == ACTIVE)) {
      k++;
      indexVIMP[k] = i;
    }
  }
  if (k != permuteObsSize) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  VIMP candidate selection failed.");
    Rprintf("\nRSF:  %10d available, %10d selected.", permuteObsSize, k);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (!(_opt & (~OPT_VIMP) & OPT_VIMP_JOIN)) {
    double  *originalVIMP = dvector(1, permuteObsSize);
    for (p=1; p <= varCount; p++) {
      for (k=1; k<= permuteObsSize; k++) {
        originalVIMP[k] = predictorPtr[_predictorIndex[p]][indexVIMP[k]];
      }
      permute(permuteObsSize, permuteVIMP);
      for (k=1; k <= permuteObsSize; k++) {
        predictorPtr[_predictorIndex[p]][indexVIMP[k]] = originalVIMP[permuteVIMP[k]];
      }
      for (i=1; i <= obsSize; i++) {
        if ((_genericMembershipFlag[_individualIndex[i]] == selectionFlag) || (selectionFlag == ACTIVE)) {
          terminalNode = getProxyMember(rootPtr, predictorPtr, i);
          if (!ISNA(terminalNode -> mortality)) {
            _vimpMortality[p][i] += terminalNode -> mortality;
            if (_eventTypeSize > 1) {
              for (j=1; j<= _eventTypeSize; j++) {
                _crVimpPOE[p][j][i] += (double) (terminalNode -> poe)[j] / (terminalNode -> eventCount);
                for (k=1; k <= _sortedTimeInterestSize; k++) {
                  _crVimpEnsemble[p][j][k][i] += terminalNode -> subSurvival[j][k];
                }
              }
            }
          }
          else {
            if (_opt & OPT_VOUT_TYPE) {
              _oobVimpInvalidDen[p][i] ++;
            }
            else {
              Rprintf("\nRSF:  *** ERROR *** ");
              Rprintf("\nRSF:  NA encountered for mortality in VIMP.");
              Rprintf("\nRSF:  Please Contact Technical Support.");
              Rprintf("\nRSF:  The application will now exit.\n");
              exit(TRUE);
            }
          }
        }
      }
      for (k=1; k <= permuteObsSize; k++) {
        predictorPtr[_predictorIndex[p]][indexVIMP[k]] = originalVIMP[k];
      }
    }  
    free_dvector(originalVIMP, 1, permuteObsSize);
  }
  else {
    double **intrOriginalVIMP = dmatrix(1, _intrPredictorSize, 1, permuteObsSize);
    for (p=1; p <= _intrPredictorSize; p++) {
      for (k=1; k<= permuteObsSize; k++) {
        intrOriginalVIMP[p][k] = predictorPtr[_intrPredictor[p]][indexVIMP[k]];
      }
      permute(permuteObsSize, permuteVIMP);
      for (k=1; k <= permuteObsSize; k++) {
        predictorPtr[_intrPredictor[p]][indexVIMP[k]] = intrOriginalVIMP[p][permuteVIMP[k]];
      }
    }
    for (i=1; i <= _fobservationSize; i++) {
      if ( _bootMembershipFlag[_intrIndividual[i]] == FALSE ) {
        terminalNode = getProxyMember(rootPtr, predictorPtr, i);
        if (!ISNA(terminalNode -> mortality)) {
          _vimpMortality[1][i] += terminalNode -> mortality;
          if (_eventTypeSize > 1) {
            for (j=1; j<= _eventTypeSize; j++) {
              _crVimpPOE[1][j][i] += (double) (terminalNode -> poe)[j] / (terminalNode -> eventCount);
              for (k=1; k <= _sortedTimeInterestSize; k++) {
                _crVimpEnsemble[1][j][k][i] += terminalNode -> subSurvival[j][k];
              }
            }
          }
        }
        else {
          if (_opt & OPT_VOUT_TYPE) {
            _oobVimpInvalidDen[1][i] ++;
          }
          else {
            Rprintf("\nRSF:  *** ERROR *** ");
            Rprintf("\nRSF:  NA encountered for mortality in VIMP.");
            Rprintf("\nRSF:  Please Contact Technical Support.");
            Rprintf("\nRSF:  The application will now exit.\n");
            exit(TRUE);
          }
        }
      }
    }
    for (p=1; p <= _intrPredictorSize; p++) {
      for (k=1; k <= permuteObsSize; k++) {
        predictorPtr[_intrPredictor[p]][indexVIMP[k]] = intrOriginalVIMP[p][k];
      }
    }
    free_dmatrix(intrOriginalVIMP, 1, _intrPredictorSize, 1, permuteObsSize);
  }  
  free_uivector(indexVIMP, 1, permuteObsSize);
  free_uivector(permuteVIMP, 1, permuteObsSize);
} 
void finalizeVariableImportance(uint       mode,
                                uint       rejectedTreeCount, 
                                char     **dmRecordBootFlag,
                                double  ***dmvImputation) {
  uint      obsSize;
  uint      varCount;
  double   *statusPtr;
  double   *timePtr;
  uint     *ensembleDenPtr;
  double  **poePtr;
  double    concordanceIndex;
  int       concordancePolarity;
  char      concordanceImputeFlag;  
  double   *crPerformanceVector;
  double ***crVimpMortality;
  double    value;
  uint     *denominatorCount;
  uint i, j, k, n, p;
  if (!(rejectedTreeCount < _forestSize)) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Insufficient trees for VIMP analysis.  \n");
    return;
  }
  if (!(_opt & OPT_VIMP)) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  VIMP analysis requested while OPT bit not set.  \n");
    return;
  }
  crPerformanceVector = NULL;  
  crVimpMortality     = NULL;  
  if (_opt & (OPT_POUT_TYPE)) {
    concordancePolarity = -1;
  }
  else {
    concordancePolarity = 1;
  }
  concordanceImputeFlag = FALSE;
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    varCount = _xSize;
    statusPtr = _status;
    timePtr = _time;
    ensembleDenPtr = _oobEnsembleDen;
    poePtr = _oobPOEPtr;
    if (_mRecordSize > 0) {
      concordanceImputeFlag = TRUE;
    }
    break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    varCount = _xSize;
    statusPtr = _fstatus;
    timePtr = rsf_ftime;
    ensembleDenPtr = _fullEnsembleDen;
    poePtr = _fullPOEPtr;
    if (_fmRecordSize > 0) {
      concordanceImputeFlag = TRUE;
    }
    break;
  case RSF_INTR:
    obsSize = _fobservationSize;
    if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
      varCount = 1;
    }
    else {
      varCount = _intrPredictorSize;
    }
    statusPtr = _fstatus;
    timePtr = rsf_ftime;
    ensembleDenPtr = _oobEnsembleDen;
    poePtr = _oobPOEPtr;
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
  if (_opt & OPT_VOUT_TYPE) {
    denominatorCount = uivector(1, obsSize);
  }
  else {
    denominatorCount = ensembleDenPtr;
  }
  if (_eventTypeSize > 1) {
    crVimpMortality = dmatrix3(1, varCount, 1, _eventTypeSize, 1, obsSize);
    crPerformanceVector = dvector(1, _eventTypeSize);
    for (p=1; p <= varCount; p++) {
      for (i = 1; i <= obsSize; i++) {
        for (j = 1; j <= _eventTypeSize; j++) {
          for (k = 1; k <= _sortedTimeInterestSize; k++) {
            if(_crVimpEnsemble[p][j][k][i] > 0) {
              if (_crVimpPOE[p][j][i] > 0) {
                value = _crVimpEnsemble[p][j][k][i] / _crVimpPOE[p][j][i];
                value = (value <= 1.0) ? value : 1.0;
                _crVimpEnsemble[p][j][k][i] = - log (value);
              }
              else {
                value = _crVimpEnsemble[p][j][k][i] / 1.0;
                value = (value <= 1.0) ? value : 1.0;
                _crVimpEnsemble[p][j][k][i] = - log (value);
              }
            }
            else {
              if (_crVimpPOE[p][j][i] > 0) {
                if (k > 1) {
                  _crVimpEnsemble[p][j][k][i] = _crVimpEnsemble[p][j][k-1][i];
                }
                else {
                  _crVimpEnsemble[p][j][k][i] = 0.0;
                }
              }
              else {
                _crVimpEnsemble[p][j][k][i] = 1.0;
              }
            }
          }
        }
      }  
    }  
    for (p = 1; p <= varCount; p++) {
      for (j = 1; j <= _eventTypeSize; j++) {
        for (i = 1; i <= obsSize; i++) {
          crVimpMortality[p][j][i] = 0.0;
            for (k = 1; k <= _sortedTimeInterestSize; k++) {            
              crVimpMortality[p][j][i] += _crVimpEnsemble[p][j][k][i];
            }
        }
      }
    }
  }  
  if (concordanceImputeFlag == TRUE) {
    imputeConcordance(mode,
                      _forestSize,
                      dmRecordBootFlag,
                      dmvImputation,
                      statusPtr,
                      timePtr);
  }
  for (p=1; p <= varCount; p++) {
    for (i = 1; i <= obsSize; i++) {
      if (_opt & OPT_VOUT_TYPE) {
        denominatorCount[i] = ensembleDenPtr[i] - _oobVimpInvalidDen[p][i];
      }
      if (denominatorCount[i] != 0) {
        _vimpMortality[p][i] = _vimpMortality[p][i] / denominatorCount[i];
      }
    }
    concordanceIndex = getConcordanceIndex(concordancePolarity,
                                           obsSize, 
                                           statusPtr,
                                           timePtr,
                                           _vimpMortality[p], 
                                           denominatorCount);
    if (ISNA(concordanceIndex)) {
      _importancePtr[1][p] = NA_REAL;
    }
    else {
      _importancePtr[1][p] = 1 - concordanceIndex;
    }
    if (_eventTypeSize > 1) {
      getConditionalPerformance(mode,
                                concordancePolarity, 
                                obsSize, 
                                statusPtr,
                                timePtr,
                                crVimpMortality[p], 
                                denominatorCount,
                                crPerformanceVector);
      for (j=1; j <=_eventTypeSize; j++) {
        _importancePtr[1+j][p] = crPerformanceVector[j];
      }
    }  
  }  
  if (_eventTypeSize > 1) {
    free_dvector(crPerformanceVector, 1, _eventTypeSize);
    free_dmatrix3(crVimpMortality, 1, varCount, 1, _eventTypeSize, 1, obsSize);
  }
  if (_opt & OPT_VOUT_TYPE) {
    free_uivector(denominatorCount, 1, obsSize);
  }
} 
