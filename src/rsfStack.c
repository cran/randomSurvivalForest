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
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include       "nodeOps.h"
#include  "rsfFactorOps.h"
#include     "rsfImpute.h"
#include       "rsfTree.h"
#include      "rsfStack.h"
void stackPreDefinedCommonArrays() {
  _nodeMembership = nodePtrVector(1, _observationSize);
  _bootMembershipIndex = uivector(1, _observationSize);
  _bootMembershipFlag = cvector(1, _observationSize);
  _masterTime  = dvector(1, _observationSize);
  _masterTimeIndex  = uivector(1, _observationSize);
  _oobSampleSize = uivector(1, _forestSize);
}
void unstackPreDefinedCommonArrays() {
  free_nodePtrVector(_nodeMembership, 1, _observationSize);
  free_uivector(_bootMembershipIndex, 1, _observationSize);
  free_cvector(_bootMembershipFlag, 1, _observationSize);
  free_dvector(_masterTime, 1, _observationSize);
  free_uivector(_masterTimeIndex, 1, _observationSize);
  free_uivector(_oobSampleSize, 1, _forestSize);
}
void stackPreDefinedGrowthArrays() {
  uint i;
  _individualIndex = uivector(1, _observationSize);
  for (i = 1; i <= _observationSize; i++) {
    _individualIndex[i] = i;
  }
  _predictorIndex  = uivector(1, _xSize);
  for (i = 1; i <= _xSize; i++) {
    _predictorIndex[i] = i;
  }
  _genericMembershipFlag = _bootMembershipFlag;
}
void unstackPreDefinedGrowthArrays() {
  free_uivector(_individualIndex, 1, _observationSize);
  free_uivector(_predictorIndex, 1, _xSize);
}
void stackPreDefinedPredictArrays() {
  uint i;
  _fnodeMembership = nodePtrVector(1, _fobservationSize);
  _individualIndex = uivector(1, _fobservationSize);
  for (i = 1; i <= _fobservationSize; i++) {
    _individualIndex[i] = i;
  }
  _predictorIndex = uivector(1, _xSize);
  for (i = 1; i <= _xSize; i++) {
    _predictorIndex[i] = i;
  }
  _genericMembershipFlag = cvector(1, _fobservationSize);
  for (i = 1; i <= _fobservationSize; i++) {
    _genericMembershipFlag[i] = ACTIVE;
  }
}
void unstackPreDefinedPredictArrays() {
  free_nodePtrVector(_fnodeMembership, 1, _fobservationSize);
  free_uivector(_individualIndex, 1, _fobservationSize);
  free_uivector(_predictorIndex, 1, _xSize);
  free_cvector(_genericMembershipFlag, 1, _fobservationSize);
}
void stackPreDefinedInteractionArrays() {
  uint i;
  _fnodeMembership = nodePtrVector(1, _fobservationSize);
  _foobSampleSize = uivector(1, _forestSize);
  _importanceFlag = cvector(1, _xSize);
  for (i = 1; i <= _xSize; i++) {
    _importanceFlag[i] = FALSE;
  }
  for (i = 1; i <= _intrPredictorSize; i++) {
    _importanceFlag[_intrPredictor[i]] = TRUE;
  }
  _individualIndex = _intrIndividual;
  _predictorIndex = _intrPredictor;
  _genericMembershipFlag = _bootMembershipFlag;
}
void unstackPreDefinedInteractionArrays() {
  free_nodePtrVector(_fnodeMembership, 1, _fobservationSize);
  free_uivector(_foobSampleSize, 1, _forestSize);
  free_cvector(_importanceFlag, 1, _xSize);
}
void initializeArrays(char mode) {
  uint i,j;
  uint leadingIndex;
  _masterTimeSize = 0;
  for (j = 1; j <= _observationSize; j++) {
    if (!ISNA(_time[j])) {
      _masterTimeSize ++;
      _masterTime[_masterTimeSize] = _time[j];
    }
  }
  hpsort(_masterTime, _masterTimeSize);
  leadingIndex = 1;
  for (i=2; i <= _masterTimeSize; i++) {
    if (_masterTime[i] > _masterTime[leadingIndex]) {
      leadingIndex++;
      _masterTime[leadingIndex] = _masterTime[i];
    }
  }
  _masterTimeSize = leadingIndex;
  for (i= _masterTimeSize + 1; i <= _observationSize; i++) {
    _masterTime[i] = 0;
  }
  if (!(_opt & OPT_IMPU_ONLY)) {
    hpsort(_timeInterest, _timeInterestSize);
    _sortedTimeInterestSize = 1;
    for (i=2; i <= _timeInterestSize; i++) {
      if (_timeInterest[i] > _timeInterest[_sortedTimeInterestSize]) {
        (_sortedTimeInterestSize) ++;
        _timeInterest[_sortedTimeInterestSize] = _timeInterest[i];
      }
    }
    if (_sortedTimeInterestSize != _timeInterestSize) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Time points of interest are not unique.");
      Rprintf("\nRSF:  The ensemble estimate output matrix is being");
      Rprintf("\nRSF:  resized as [N'] x [n], where N' is the");
      Rprintf("\nRSF:  unique time points of interest and n is");
      Rprintf("\nRSF:  number of observations in the data.");
    }
    for (i= (_sortedTimeInterestSize) + 1; i <= _timeInterestSize; i++) {
      _timeInterest[i] = 0;
    }
  }
}
void initializeFactorArrays(char mode) {
  uint i, j;
  uint factorLevel;
  if (_factorCount < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to initialize factorness in its absence.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    error("\nRSF:  The application will now exit.\n");
  }
  _maxFactorLevel = 0;
  for (j = 1; j <= _factorCount; j++) {
    _factorSize[j] = 0;
    for (i = 1; i <= _observationSize; i++) {
      if (!ISNA(_observation[_factorIndex[j]][i])) {
        if (_observation[_factorIndex[j]][i] >= 1) {
          factorLevel = (uint) _observation[_factorIndex[j]][i];
          if (_factorSize[j] < factorLevel) {
            _factorSize[j] = factorLevel;
          }
        }
        else {
          Rprintf("\nRSF:  *** ERROR *** ");
          Rprintf("\nRSF:  Factor level less than one (1):  %10.4f", _observation[_factorIndex[j]][i]);
          error("\nRSF:  The application will now exit.\n");
        }
      }
      if (_maxFactorLevel < _factorSize[j]) {
        _maxFactorLevel = _factorSize[j];
      }
    }
  }
  if ( mode != RSF_GROW) {
    for (j = 1; j <= _factorCount; j++) {
      factorLevel = 0;
      for (i = 1; i <= _fobservationSize; i++) {
        if (!ISNA(_fobservation[_factorIndex[j]][i])) {
          if (_fobservation[_factorIndex[j]][i] >= 1) {
            factorLevel = (factorLevel > (uint) _fobservation[_factorIndex[j]][i]) ? factorLevel : ((uint) _fobservation[_factorIndex[j]][i]);
          }
          else {
            Rprintf("\nRSF:  *** ERROR *** ");
            Rprintf("\nRSF:  Factor level less than one (1):  %10.4f", _fobservation[_factorIndex[j]][i]);
            error("\nRSF:  The application will now exit.\n");
          }
        }
      }
      if (factorLevel > _factorSize[j]) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  !GROW factor level greater than maximum GROW factor level:  %10d vs. %10d", factorLevel, _factorSize[j]);
        error("\nRSF:  The application will now exit.\n");
      }
    }
  }
  _factorList = factorPtrVector(1, _maxFactorLevel);
  for (j = 1; j <= _maxFactorLevel; j++) {
    _factorList[j] = NULL;
  }
  for (j = 1; j <= _factorCount; j++) {
    _factorList[_factorSize[j]] = makeFactor(_factorSize[j], FALSE);
  }
}
char stackCompetingArrays(char mode) {
  uint obsSize;
  double  *statusPtr;
  uint    *mRecordMap;
  int    **mvSign;
  char eventAnalysisFlag, consistencyFlag, overWriteFlag;
  char statusFlag;
  uint *eventCounter;
  uint i, j, jgrow;
  _eventType = uivector(1, _observationSize);
  overWriteFlag = FALSE;
  if (mode == RSF_GROW) {
    if (_splitRule < SPLIT_RULE_CR) {
      overWriteFlag = TRUE;
    }
  }
  getEventTypeSize(_observationSize, 
                   _status, 
                   _mRecordMap, 
                   _mvSign,  
                   overWriteFlag, 
                   & _eventTypeSize,
                   & _mStatusSize,
                   _eventType);
  if (mode == RSF_GROW) {
    if (_splitRule > SPLIT_RULE_CR) {
      if (_eventTypeSize == 1) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Split rule specified is for Competing Risk scenarios only.");
        Rprintf("\nRSF:  The data set does not contain multiple events.");
        return FALSE;
      }
    }
  }
  if (_eventTypeSize > 1) {
    _eventTypeIndex  = uivector(1, _eventType[_eventTypeSize]);
    for (j = 1; j <= _eventType[_eventTypeSize]; j++) {
      _eventTypeIndex[j] = 0;
    }
    for (j = 1; j <= _eventTypeSize; j++) {
      _eventTypeIndex[_eventType[j]] = j;
    }
  }
  switch (mode) {
  case RSF_PRED:
    if ((_opt & OPT_PERF) || (_opt & OPT_VIMP)) {
      eventAnalysisFlag = TRUE;        
    }
    else {
      eventAnalysisFlag = FALSE;
    }
    break;
  case RSF_INTR:
    eventAnalysisFlag = TRUE;
    break;
  default:
    eventAnalysisFlag = FALSE;
    break;
  }  
  if (eventAnalysisFlag == TRUE) {
      uint feventTypeSize;
      uint *feventType;
      feventType = uivector(1, _fobservationSize);
      feventTypeSize = 0;
      consistencyFlag = TRUE;
      overWriteFlag = FALSE;
      if (_eventTypeSize == 1) {
        overWriteFlag = TRUE;
      }
      getEventTypeSize(_fobservationSize, 
                       _fstatus, 
                       _fmRecordMap, 
                       _fmvSign,  
                       overWriteFlag, 
                       & feventTypeSize,
                       & _mStatusSize,
                       feventType);
      if (_eventTypeSize > 1) {
        for (j = 1; j <= feventTypeSize; j++) {
          for (jgrow = 1; jgrow <= _eventTypeSize; jgrow++) {
            if (feventType[j] != _eventType[jgrow]) {
              if (jgrow == _eventTypeSize) {
                consistencyFlag = FALSE;
              }
            }
            else {
              jgrow = _eventTypeSize;
            }
          }
        }
      }
      free_uivector(feventType, 1, _fobservationSize);
      if (consistencyFlag == FALSE) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Unknown event type encountered in !GROW mode. ");
        Rprintf("\nRSF:  Please Contact Technical Support.");
        return FALSE;
      }
  }  
  if (_eventTypeSize > 1) {
    if (mode == RSF_PRED) {
      if ((_opt & OPT_PERF) || (_opt & OPT_VIMP)) {
        eventAnalysisFlag = TRUE;        
      }
      else {
        eventAnalysisFlag = FALSE;
      }
    }
    else {
      eventAnalysisFlag = TRUE;
    }
  }
  else {
    eventAnalysisFlag = FALSE;
  }
  if (eventAnalysisFlag == TRUE) {
    if (mode == RSF_GROW) {
      obsSize = _observationSize;
      statusPtr = _status;
      mvSign = _mvSign;
      mRecordMap = _mRecordMap;
    }
    else {
      obsSize = _fobservationSize;
      statusPtr = _fstatus;
      mvSign = _fmvSign;
      mRecordMap = _fmRecordMap;
    }
    _eIndividualSize  = uivector(1, _eventTypeSize);
    for (j = 1; j <= _eventTypeSize; j++) {
      _eIndividualSize[j] = 0;
    }
    for (i = 1; i <= obsSize; i++) {
      statusFlag = FALSE;
      if (mRecordMap[i] == 0) { 
        statusFlag = TRUE;
      }
      else {
        if (mvSign[(uint) abs(CENS_IDX)][mRecordMap[i]] == 0) {
          statusFlag = TRUE;
        }
      }
      if (statusFlag == TRUE) {
        if ((uint) statusPtr[i] > 0) {
          _eIndividualSize[_eventTypeIndex[(uint) statusPtr[i]]] ++;
        } 
        else {
          for (j=1; j <= _eventTypeSize; j++) {
            _eIndividualSize[j] ++;
          }
        }
      }
    }  
    _eIndividual = puivector(1, _eventTypeSize);
    for (j = 1; j <= _eventTypeSize; j++) {
      _eIndividual[j] = uivector(1, _eIndividualSize[j] + _mStatusSize + 1);
    }
    eventCounter = uivector(1, _eventTypeSize);
    for (j = 1; j <= _eventTypeSize; j++) {
      eventCounter[j] = 0;
    }
    for (i = 1; i <= obsSize; i++) {
      statusFlag = FALSE;
      if (mRecordMap[i] == 0) { 
        statusFlag = TRUE;
      }
      else {
        if (mvSign[(uint) abs(CENS_IDX)][mRecordMap[i]] == 0) {
          statusFlag = TRUE;
        }
      }
      if (statusFlag == TRUE) {
        if ((uint) statusPtr[i] > 0) {
          j = _eventTypeIndex[(uint) statusPtr[i]];
          eventCounter[j] ++;
          _eIndividual[j][eventCounter[j]] = i;
        } 
        else {
          for (j=1; j <= _eventTypeSize; j++) {
            eventCounter[j] ++;
            _eIndividual[j][eventCounter[j]] = i;
          }
        }
      }
    }
    free_uivector(eventCounter, 1, _eventTypeSize);
    if (_mStatusSize > 0) {
      _meIndividualSize  = uivector(1, _eventTypeSize);
    }
    else {
      _meIndividualSize  = _eIndividualSize;
    }
  }  
  return TRUE;
}
void getEventTypeSize(uint     obsSize, 
                      double  *status, 
                      uint    *mRecordMap, 
                      int    **mvSign,  
                      char     overWriteFlag,
                      uint    *eventTypeSize,
                      uint    *msize,
                      uint    *eventType) {
  uint statusFlag;
  uint leadingIndex;
  uint i;
  *eventTypeSize = *msize = 0;
  for (i = 1; i <= obsSize; i++) {
    statusFlag = FALSE;
    if (mRecordMap[i] == 0) { 
      statusFlag = TRUE;
    }
    else {
      if (mvSign[(uint) abs(CENS_IDX)][mRecordMap[i]] == 0) {
        statusFlag = TRUE;
      }
    }
    if (statusFlag == TRUE) {
      if ((uint) status[i] > 0) {
        if (overWriteFlag == TRUE) { 
          status[i] = 1;
        }
        else {
          (*eventTypeSize) ++;
          eventType[*eventTypeSize] = (uint) status[i];
        }
      }  
      else {
      }
    }
    else {
      (*msize) ++;
    }
  }  
  if (overWriteFlag == TRUE) {
    *eventTypeSize = 1;
  }
  else {
    hpsortui(eventType, *eventTypeSize);
    leadingIndex = 1;
    for (i=2; i <= *eventTypeSize; i++) {
      if (eventType[i] > eventType[leadingIndex]) {
        leadingIndex++;
        eventType[leadingIndex] = eventType[i];
      }
    }
    *eventTypeSize = leadingIndex;
    for (i= *eventTypeSize + 1; i <= obsSize; i++) {
      eventType[i] = 0;
    }
  }
}
void unstackCompetingArrays(char mode) {
  char eventAnalysisFlag;
  uint j;
  free_uivector(_eventType, 1, _observationSize);
  if (_eventTypeSize > 1) {
    if (mode == RSF_PRED) {
      if ((_opt & OPT_PERF) || (_opt & OPT_VIMP)) {
        eventAnalysisFlag = TRUE;        
      }
      else {
        eventAnalysisFlag = FALSE;
      }
    }
    else {
      eventAnalysisFlag = TRUE;
    }
  }
  else {
    eventAnalysisFlag = FALSE;
  }
  if (eventAnalysisFlag == TRUE) {
    if (_mStatusSize > 0) {
      free_uivector(_meIndividualSize, 1, _eventTypeSize);
    }
    for (j = 1; j <= _eventTypeSize; j++) {
      free_uivector(_eIndividual[j], 1, _eIndividualSize[j] + _mStatusSize + 1);
    }
    free_puivector(_eIndividual, 1, _eventTypeSize);
    free_uivector(_eIndividualSize, 1, _eventTypeSize);
    free_uivector(_eventTypeIndex, 1, _eventType[_eventTypeSize]);
  }  
}
void stackFactorArrays(char mode) {
  uint j, p;
  _xType = pcvector(1, _xSize);
  for (p = 1; p <= _xSize; p++) {
    _xType[p] = (char*) CHAR(STRING_ELT(AS_CHARACTER(_sexp_xType), p-1));
    if ((strcmp(_xType[p], "C") != 0) && (strcmp(_xType[p], "I") != 0) && (strcmp(_xType[p], "R") != 0)) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Invalid predictor type:  [%10d] = %2s", p, _xType[p]);
      Rprintf("\nRSF:  Type must be 'C', 'I', or 'R'.");
      error("\nRSF:  The application will now exit.\n");
    }
  }
  _factorMap = uivector(1, _xSize);
  _factorCount = 0;
  for (p = 1; p <= _xSize; p++) {
    _factorMap[p] = 0;
    if (strcmp(_xType[p], "C") == 0) {
      _factorCount ++;
      _factorMap[p] = _factorCount;
    }
  }
  if (_factorCount > 0) {
    _factorIndex = uivector(1, _factorCount);
    j = 0;
    for (p = 1; p <= _xSize; p++) {
      if (_factorMap[p] > 0) {
        _factorIndex[++j] = p;
      }
    }
    _factorSize = uivector(1, _factorCount);
  }
}
void unstackFactorArrays(char mode) {
  uint j;
  free_pcvector(_xType, 1, _xSize);
  free_uivector(_factorMap, 1, _xSize);
  if (_factorCount > 0) {
    free_uivector(_factorIndex, 1, _factorCount);
    free_uivector(_factorSize, 1, _factorCount);
    for (j = 1; j <= _maxFactorLevel; j++) {
      if (_factorList[j] != NULL) {
        free_Factor(_factorList[j]);
      }
    }
    free_factorPtrVector(_factorList, 1, _maxFactorLevel);
  }
}
char stackMissingSignatures(uint     obsSize, 
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
                            uint   **p_mFactorIndex) {
  char result;
  uint i, j, p;
  if (recordSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to allocate for missingness in its absence.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    error("\nRSF:  The application will now exit.\n");
  }
  *p_recordIndex = uivector(1, recordSize);
  i = 0;
  for (j = 1; j <= obsSize; j++) {
    if (recordMap[j] > 0) {
      i++;
      (*p_recordIndex)[i] = j;
    }
  }
  *p_vSign = imatrix(1, _xSize+2, 1, recordSize);
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= _xSize+2; i++) {
      (*p_vSign)[i][j] = 0;
    }
  }
  for (j = 1; j <= recordSize; j++) {
    if (ISNA(statusPtr[(*p_recordIndex)[j]])) {
      (*p_vSign)[(uint) abs(CENS_IDX)][j] = 1;
    }
    if (ISNA(timePtr[(*p_recordIndex)[j]])) {
      (*p_vSign)[(uint) abs(TIME_IDX)][j] = 1;
    }
    for (i = 1; i <= _xSize; i++) {
      if (ISNA(predictorPtr[i][(*p_recordIndex)[j]])) {
        (*p_vSign)[i+2][j] = 1;
      }
    }
  }
  result = FALSE;
  *p_vIndex = ivector(1, _xSize+2);
  *p_vSize = 0;
  for (i = 1; i <= _xSize+2; i++) {
    (*p_vIndex)[i] = 0;
    switch (i) {
    case (uint) (-CENS_IDX):
      for (j = 1; j <= recordSize; j++) {
        if ((*p_vSign)[i][j] == 1) {
          (*p_vSize) ++;
          (*p_vIndex)[*p_vSize] = CENS_IDX;
          j = recordSize;
        }
      }
      break;
    case (uint) (-TIME_IDX):
      for (j = 1; j <= recordSize; j++) {
        if ((*p_vSign)[i][j] == 1) {
          (*p_vSize) ++;
          (*p_vIndex)[*p_vSize] = TIME_IDX;
          result = TRUE;
          j = recordSize;
        }
      }
      break;
    default:
      for (j = 1; j <= recordSize; j++) {
        if ((*p_vSign)[i][j] == 1) {
          (*p_vSize) ++;
          (*p_vIndex)[*p_vSize] = i-2;
          j = recordSize;
        }
      }
      break;
    }
  }  
  *p_vForestSign = imatrix(1, _forestSize,  1, *p_vSize);
  for (j = 1; j <= _forestSize; j++) {
    for (p = 1; p <= (*p_vSize); p++) {
      (*p_vForestSign)[j][p] = -1;
    }
  }
  *p_mFactorIndex = uivector(1, _xSize);
  for (p = 1; p <= _xSize; p++) {
    (*p_mFactorIndex)[p] = 0;
  }
  *p_mFactorSize = 0;
  for (p = 1; p <= *p_vSize; p++) {
    switch ((*p_vIndex)[p]) {
    case CENS_IDX:
      break;
    case TIME_IDX:
      break;
    default:
      if (strcmp(_xType[(*p_vIndex)[p]], "C") == 0) {
        (*p_mFactorSize) ++;
        (*p_mFactorIndex)[*p_mFactorSize] = (*p_vIndex)[p];
      }
      break;
    }
  }
  return result;
}
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
                              uint     *mFactorIndex) {
  free_dvector(statusPtr, 1, obsSize);
  free_dvector(timePtr, 1, obsSize);
  free_dmatrix(predictorPtr, 1, _xSize, 1, obsSize);
  free_uivector(recordIndex, 1, recordSize);
  free_imatrix(vSign, 1, _xSize+2, 1, recordSize);
  free_ivector(vIndex, 1, _xSize+2);
  free_imatrix(vForestSign, 1, _forestSize, 1, vSize);
  free_uivector(mFactorIndex, 1, _xSize);
}
char stackMissingArrays(char       mode,
                        char    ***p_dmRecordBootFlag,
                        double ****p_dmvImputation) {
  char result;
  char mFlag;
  char dualUseFlag;
  uint recordSize;
  uint vSize;
  uint i,j,k,p;
  result = TRUE;
  for (i = 1 ; i <= _observationSize; i++) {
    if (!ISNA(_time[i])) {
      if (_time[i] <= 0) {
        result = FALSE;
        Rprintf("\nRSF:  TRAINING time elements must be greater than zero or NA:  [%10d] = %12.4f \n", i, _time[i]);
      }
    }
  }
  if (result == FALSE) {
    return result;
  }
  mFlag = FALSE;
  for (i = 1 ; i <= _observationSize; i++) {
    if (!ISNA(_time[i])) {
      mFlag = TRUE;
      i = _observationSize;
    }
  }
  if (mFlag == FALSE) {
    Rprintf("\nRSF:  All TRAINING time elements are missing. \n");
    result = FALSE;
  }
  mFlag = FALSE;
  for (i = 1 ; i <= _observationSize; i++) {
    if (!ISNA(_status[i]) && (_status[i] > 0)) {
      mFlag = TRUE;
      i = _observationSize;
    }
  }
  if (mFlag == FALSE) {
    Rprintf("\nRSF:  All TRAINING status elements are are censored or missing. \n");
    result = FALSE;
  }
  if (result == FALSE) {
    return result;
  }
  if (mode == RSF_PRED) {
    if ((_opt & OPT_PERF) || (_opt & OPT_VIMP)) {
      for (i = 1 ; i <= _fobservationSize; i++) {
        if (!ISNA(rsf_ftime[i])) {
          if (rsf_ftime[i] <= 0) {
            result = FALSE;
            Rprintf("\nRSF:  PRED time elements must be greater than zero or NA when PERF or VIMP is requested:  [%10d] = %12.4f \n", i, rsf_ftime[i]);
          }
        }
      }
    }
    else {
      for (i = 1 ; i <= _fobservationSize; i++) {
        _fstatus[i] = rsf_ftime[i] = -1;
      }
    }
    if (result == FALSE) {
      return result;
    }
  }  
  _mRecordMap = uivector(1, _observationSize);
  _mRecordSize = getRecordMap(_mRecordMap, 
                              _observationSize, 
                              _status, 
                              _time, 
                              _xData);
  if (_mRecordSize == 0) {
    _mTimeIndexFlag = FALSE;
    _observation = pdvector(1, _xSize);
    for (i=1; i <= _xSize; i++) {
      _observation[i] = (_xData + ((i-1)*(_observationSize)) - 1);
    }
    if (mode == RSF_GROW) {
      _imputeSize = 1;
    }
  }
  else {
    _mStatus = dvector(1, _observationSize);
    _mTime   = dvector(1, _observationSize);
    for (i = 1; i <= _observationSize; i++) {
      _mStatus[i] = _status[i];
      _mTime[i] = _time[i];
    }
    _status = _mStatus;
    _time   = _mTime;
    _observation = dmatrix(1, _xSize, 1, _observationSize);
    for (p=1; p <= _xSize; p++) {
      for (i = 1; i <= _observationSize; i++) {
        _observation[p][i] = (_xData+((p-1)*(_observationSize)))[i-1];
      }
    }
    _mTimeIndexFlag = stackMissingSignatures(_observationSize,
                                             _status,
                                             _time,
                                             _observation,
                                             _mRecordMap,
                                             _mRecordSize,
                                             & _mRecordIndex,
                                             & _mvSize,
                                             & _mvSign,
                                             & _mvIndex,
                                             & _mvForestSign,
                                             & _mFactorSize,
                                             & _mFactorIndex);
  }  
  if (mode == RSF_PRED) {
    _fmRecordMap = uivector(1, _fobservationSize);
    _fmRecordSize = getRecordMap(_fmRecordMap, 
                                 _fobservationSize, 
                                 _fstatus, 
                                 rsf_ftime, 
                                 _fxData);
    if (_fmRecordSize == 0) {
      _fobservation = pdvector(1, _xSize);
      for (j=1; j <= _xSize; j++) {
        _fobservation[j] = (_fxData + ((j-1)*(_fobservationSize)) - 1);
      }
    }  
    else {
      _fmStatus = dvector(1, _fobservationSize);
      _fmTime   = dvector(1, _fobservationSize);
      for (i = 1; i <= _fobservationSize; i++) {
        _fmStatus[i] = _fstatus[i];
        _fmTime[i] = rsf_ftime[i];
      }
      _fstatus = _fmStatus;
      rsf_ftime   = _fmTime;
      _fobservation = dmatrix(1, _xSize, 1, _fobservationSize);
      for (p=1; p <= _xSize; p++) {
        for (i = 1; i <= _fobservationSize; i++) {
          _fobservation[p][i] = (_fxData+((p-1)*(_fobservationSize)))[i-1];
        }
      }
      stackMissingSignatures(_fobservationSize,
                             _fstatus,
                             rsf_ftime,
                             _fobservation,
                             _fmRecordMap,
                             _fmRecordSize,
                             & _fmRecordIndex,
                             & _fmvSize,
                             & _fmvSign,
                             & _fmvIndex,
                             & _fmvForestSign,
                             & _fmFactorSize,
                             & _fmFactorIndex);
    }  
  }  
  if (mode == RSF_INTR) {
    _fstatus = dvector(1, _fobservationSize);
    rsf_ftime   = dvector(1, _fobservationSize);
    _fobservation = dmatrix(1, _xSize, 1, _fobservationSize);
    for (j = 1; j <= _fobservationSize; j++) {
      rsf_ftime[j] = _time[_intrIndividual[j]];
      _fstatus[j] = _status[_intrIndividual[j]];
      for (i=1; i <= _xSize; i++) {
        _fobservation[i][j] = (_xData+((i-1)*(_observationSize)))[_intrIndividual[j]-1];
      }
    }
    _fmRecordMap = uivector(1, _fobservationSize);
    _fmRecordSize = getRecordMap(_fmRecordMap, 
                                 _fobservationSize, 
                                 _fstatus, 
                                 rsf_ftime, 
                                 &_fobservation[1][1]);
    if (_fmRecordSize == 0) {
    }
    else {
      stackMissingSignatures(_fobservationSize,
                             _fstatus,
                             rsf_ftime,
                             _fobservation,
                             _fmRecordMap,
                             _fmRecordSize,
                             & _fmRecordIndex,
                             & _fmvSize,
                             & _fmvSign,
                             & _fmvIndex,
                             & _fmvForestSign,
                             & _fmFactorSize,
                             & _fmFactorIndex);
    }  
  }  
  dualUseFlag = FALSE;
  switch (mode) {
  case RSF_PRED:
    if (_fmRecordSize > 0) {
      recordSize = _fmRecordSize;
      vSize = _fmvSize;
      dualUseFlag = TRUE;
      mFlag = ACTIVE;
    }
    else {
      _opt = _opt & (~OPT_MISS);
    }
    break;
  case RSF_INTR:
    if (_mRecordSize > 0) {
      recordSize = _fmRecordSize;
      vSize = _fmvSize;
      dualUseFlag = TRUE;
      mFlag = FALSE;
    }
    else {
      _opt = _opt & (~OPT_MISS);
      _opt = _opt & (~OPT_OMIS);
    }
    break;
  default:
    if (_mRecordSize > 0) {
      recordSize = _mRecordSize;
      vSize = _mvSize;
      dualUseFlag = TRUE;
      mFlag = FALSE;
    }
    else {
      _opt = _opt & (~OPT_MISS);
      _opt = _opt & (~OPT_OMIS);
    }
    break;
  }  
  if (dualUseFlag == TRUE) {
    *p_dmRecordBootFlag = cmatrix(1, _forestSize, 1, recordSize);
    for (j = 1; j <= _forestSize; j++) {
      for (i = 1; i <= recordSize; i++) {
        (*p_dmRecordBootFlag)[j][i] = mFlag;
      }
    }
    *p_dmvImputation = dmatrix3(1, _forestSize, 1, recordSize, 1, vSize);
    for (i = 1; i <= _forestSize; i++) {
      for (j = 1; j <= recordSize; j++) {
        for (k = 1; k <= vSize; k++) {
          (*p_dmvImputation)[i][j][k] = NA_REAL;
        }
      }
    }
  }
  for (i=1; i <= _observationSize; i++) {
    if (!ISNA(_time[i])) {
      k = 1;
      while (k <= _masterTimeSize) {
        if (_time[i] == _masterTime[k]) {
          _masterTimeIndex[i] = k;
          k = _masterTimeSize;
        }
        k++;
      }
    }
    else {
      _masterTimeIndex[i] = 0;
    }
  }
  if (_factorCount > 0) {
    initializeFactorArrays(mode);
  }
  return result;
}
void unstackMissingArrays(char   mode,
                          char  **dmRecordBootFlag,
                          double ***dmvImputation) {
  char dualUseFlag;
  uint recordSize;
  uint vSize;
  free_uivector(_mRecordMap, 1, _observationSize);
  if (_mRecordSize == 0) {
    free_pdvector(_observation, 1, _xSize);
  }
  else {
    unstackMissingSignatures(_observationSize,
                             _status,
                             _time,
                             _observation,
                             _mRecordMap,
                             _mRecordSize,
                             _mRecordIndex,
                             _mvSize,
                             _mvSign,
                             _mvIndex,
                             _mvForestSign,
                             _mFactorSize,
                             _mFactorIndex);
  }  
  if (mode == RSF_PRED) {
    free_uivector(_fmRecordMap, 1, _fobservationSize);
    if (_fmRecordSize == 0) {
      free_pdvector(_fobservation, 1, _xSize);
    }
    else {
      unstackMissingSignatures(_fobservationSize,
                               _fstatus,
                               rsf_ftime,
                               _fobservation,
                               _fmRecordMap,
                               _fmRecordSize,
                               _fmRecordIndex,
                               _fmvSize,
                               _fmvSign,
                               _fmvIndex,
                               _fmvForestSign,
                               _fmFactorSize,
                               _fmFactorIndex);
    }
  }  
  if (mode == RSF_INTR) {
    free_uivector(_fmRecordMap, 1, _fobservationSize);
    if (_fmRecordSize == 0) {
      free_dvector(_fstatus, 1, _fobservationSize);
      free_dvector(rsf_ftime, 1, _fobservationSize);
      free_dmatrix(_fobservation, 1, _xSize, 1, _fobservationSize);
    }
    else {
      unstackMissingSignatures(_fobservationSize,
                               _fstatus,
                               rsf_ftime,
                               _fobservation,
                               _fmRecordMap,
                               _fmRecordSize,
                               _fmRecordIndex,
                               _fmvSize,
                               _fmvSign,
                               _fmvIndex,
                               _fmvForestSign,
                               _fmFactorSize,
                               _fmFactorIndex);
   }
  }  
  dualUseFlag = FALSE;
  switch (mode) {
  case RSF_PRED:
    if (_fmRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = _fmRecordSize;
      vSize = _fmvSize;
    }
    break;
  case RSF_INTR:
    if (_fmRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = _fmRecordSize;
      vSize = _fmvSize;
    }
    break;
  default:
    if (_mRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = _mRecordSize;
      vSize = _mvSize;
    }
    break;
  }  
  if (dualUseFlag == TRUE) {
    free_cmatrix(dmRecordBootFlag, 1, _forestSize, 1, recordSize);
    free_dmatrix3(dmvImputation, 1, _forestSize, 1, recordSize, 1, vSize);
  }
}
uint stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***p_root,
                               double  **p_oobEnsemble,
                               double  **p_fullEnsemble,
                               double  **p_performance,
                               uint    **p_leafCount,
                               uint    **p_proximity,
                               double  **p_importance,
                               int     **p_seed,
                               double  **p_oobImputation,
                               double  **p_imputation,
                               double  **p_sImputeStatusPtr,
                               double  **p_sImputeTimePtr,
                               double ***p_sImputePredictorPtr,
                               double  **p_sOOBImputeStatusPtr,
                               double  **p_sOOBImputeTimePtr,
                               double ***p_sOOBImputePredictorPtr,
                               uint    **p_varUsed,
                               uint   ***p_varUsedPtr,
                               double  **p_splitDepth,
                               double ***p_splitDepthPtr,
                               double  **localSplitDepthPtr,
                               double  **p_oobEnsemblePOE,
                               double  **p_fullEnsemblePOE,
                               uint     *stackCount,
                               SEXP     *sexpVector) {
  uint sexpIndex;
  uint ensembleSize;
  uint performanceSize;
  uint proximitySize;
  uint imputationSize;
  uint importanceSize;
  uint importanceMult;
  uint varUsedSize;
  uint  obsSize;
  uint  mRecordSize;
  uint *mRecordIndex;
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  uint adjEventSize;
  uint i,j,k,p;
  sexpIndex      = 0;  
  ensembleSize   = 0;  
  performanceSize= 0;  
  proximitySize  = 0;  
  imputationSize = 0;  
  importanceSize = 0;  
  importanceMult = 0;  
  varUsedSize    = 0;  
  obsSize        = 0;  
  mRecordSize    = 0;  
  statusPtr      = NULL;  
  timePtr        = NULL;  
  predictorPtr   = NULL;  
  mRecordIndex   = NULL;  
  if (_eventTypeSize > 1) {
    adjEventSize = _eventTypeSize + 1;
  }
  else {
    adjEventSize = 1;
  }
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    mRecordSize = _mRecordSize;
    statusPtr = _status;
    timePtr = _time;
    predictorPtr = _observation;
    mRecordIndex = _mRecordIndex;
    if (_opt & OPT_IMPU_ONLY) {
      *stackCount = 1;
    }
    else {
      *stackCount = 4;
    }
    if (_opt & OPT_FENS) {
      if (_eventTypeSize > 1) {
        (*stackCount) += 1;
      }
    }
    if (_opt & OPT_OENS) {
      if (_eventTypeSize > 1) {
        (*stackCount) += 1;
      }
    }
    if (_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2; 
      (*stackCount) += 1;
    }
    if (_opt & OPT_SEED) {
      if (_opt & OPT_TREE) {
        (*stackCount) += 7;
      }
      else {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  SEXP TREE output request inconsistent.");
        Rprintf("\nRSF:  Please Contact Technical Support.");
        error("\nRSF:  The application will now exit.\n");
      }
    }
    if ((_opt & OPT_MISS) && (_opt & OPT_OMIS)) {
      imputationSize = (_xSize + 3) * mRecordSize;
      (*stackCount) += 2;
    }
    if (_opt & OPT_VUSE) {
      if (_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
        varUsedSize = _forestSize;
      }
      else {
        varUsedSize = 1;
      }
      (*stackCount) += 1;
    }
    if (_opt & OPT_SPLT_DPTH) {
      (*stackCount) += 1;
    }
    if (_opt & OPT_VIMP) {
      (*stackCount) += 1;
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceMult = 1;
        importanceSize = adjEventSize;
      }
      else {
        importanceMult = _xSize;
        importanceSize = adjEventSize * _xSize;
      }
    }
    break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    mRecordSize = _fmRecordSize;
    statusPtr = _fstatus;
    timePtr = rsf_ftime;
    predictorPtr = _fobservation;
    mRecordIndex = _fmRecordIndex;
    *stackCount = 2;
    if (_opt & OPT_FENS) {
      if (_eventTypeSize > 1) {
        (*stackCount) += 1;
      }
    }
    if (_opt & OPT_PERF) {
      (*stackCount) += 1;
    }
    if (_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2; 
      (*stackCount) += 1;
    }
    if (_opt & OPT_MISS) {
      imputationSize = (_xSize + 3) * mRecordSize;
      (*stackCount) += 1;
    }
    if (_opt & OPT_SPLT_DPTH) {
      (*stackCount) += 1;
    }
    if (_opt & OPT_VIMP) {
      (*stackCount) += 1;
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceMult = 1;
        importanceSize = adjEventSize;
      }
      else {
        importanceMult = _xSize;
        importanceSize = adjEventSize * _xSize;
      }
    }
    break;
  case RSF_INTR:
    obsSize = _fobservationSize;
    mRecordSize = _fmRecordSize;
    statusPtr = _fstatus;
    timePtr = rsf_ftime;
    predictorPtr = _fobservation;
    mRecordIndex = _fmRecordIndex;
    *stackCount = 3;
    if (_opt & OPT_OENS) {
      if (_eventTypeSize > 1) {
        (*stackCount) += 1;
      }
    }
    if (_opt & OPT_OMIS) {
      imputationSize = (_xSize + 3) * mRecordSize;
      (*stackCount) += 1;
    }
    if (_opt & OPT_SPLT_DPTH) {
      (*stackCount) += 1;
    }
    if (_opt & OPT_VIMP) {
      (*stackCount) += 1;
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceMult = 1;
        importanceSize = adjEventSize;
      }
      else {
        importanceMult = _intrPredictorSize;
        importanceSize = adjEventSize * _intrPredictorSize;
      }
    }
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    error("\nRSF:  The application will now exit.\n");
    break;
  }
  performanceSize = adjEventSize * _forestSize;
  ensembleSize = adjEventSize * _sortedTimeInterestSize * obsSize;
  PROTECT(sexpVector[RSF_OUTP_ID] = allocVector(VECSXP, *stackCount));
  PROTECT(sexpVector[RSF_STRG_ID] = allocVector(STRSXP, *stackCount));
  setAttrib(sexpVector[RSF_OUTP_ID], R_NamesSymbol, sexpVector[RSF_STRG_ID]);
  sexpIndex = 0;
  if (_opt & OPT_FENS) {
    PROTECT(sexpVector[RSF_FENS_ID] = NEW_NUMERIC(ensembleSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_FENS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_FENS_ID]));
    *p_fullEnsemble = NUMERIC_POINTER(sexpVector[RSF_FENS_ID]);
    sexpIndex ++;
    _fullEnsemblePtr = ppdvector(1, adjEventSize);
    _fullEnsembleDen = uivector(1, obsSize);
    for (j = 1; j <= adjEventSize; j++) {
      _fullEnsemblePtr[j] = pdvector(1, _sortedTimeInterestSize);
      for (k = 1; k <= _sortedTimeInterestSize; k++) {
        _fullEnsemblePtr[j][k]  = (*p_fullEnsemble) + ((j-1) * _sortedTimeInterestSize * obsSize) + ((k-1) * obsSize) - 1;
      }
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= adjEventSize; j++) {
        for (k = 1; k <= _sortedTimeInterestSize; k++) {
          _fullEnsemblePtr[j][k][i] = 0.0;
        }
      }
      _fullEnsembleDen[i] = 0;
    }
    if (_eventTypeSize > 1) {
      PROTECT(sexpVector[RSF_FPOE_ID] = NEW_NUMERIC(_eventTypeSize * obsSize));
      SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_FPOE_ID]);
      SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_FPOE_ID]));
      *p_fullEnsemblePOE = NUMERIC_POINTER(sexpVector[RSF_FPOE_ID]);
      sexpIndex ++;
      _fullPOEPtr = pdvector(1, _eventTypeSize);
      for (j = 1; j <= _eventTypeSize; j++) {
          _fullPOEPtr[j]  = (*p_fullEnsemblePOE) + ((j-1) * obsSize) - 1;
      }
      for (i = 1; i <= obsSize; i++) {
        for (j = 1; j <= _eventTypeSize; j++) {
            _fullPOEPtr[j][i] = 0.0;
        }
      }
      _fullSubSurvivalPtr     = dmatrix3(1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      _fullSubDistributionPtr = dmatrix3(1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        for (j = 1; j <= _eventTypeSize; j++) {
          for (k = 1; k <= _sortedTimeInterestSize; k++) {
            _fullSubSurvivalPtr[j][k][i] = 0.0;
            _fullSubDistributionPtr[j][k][i] = 0.0;
          }
        }
      }
    }
  }
  if (_opt & OPT_OENS) {
    PROTECT(sexpVector[RSF_OENS_ID] = NEW_NUMERIC(ensembleSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_OENS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_OENS_ID]));
    *p_oobEnsemble = NUMERIC_POINTER(sexpVector[RSF_OENS_ID]);
    sexpIndex ++;
    _oobEnsemblePtr  = ppdvector(1, adjEventSize);
    _oobEnsembleDen  = uivector(1, obsSize);
    for (j = 1; j <= adjEventSize; j++) {
      _oobEnsemblePtr[j] = pdvector(1, _sortedTimeInterestSize);
      for (k = 1; k <= _sortedTimeInterestSize; k++) {
        _oobEnsemblePtr[j][k]  = (*p_oobEnsemble) + ((j-1) * _sortedTimeInterestSize * obsSize) + ((k-1) * obsSize) - 1;
      }
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= adjEventSize; j++) {
        for (k = 1; k <= _sortedTimeInterestSize; k++) {
          _oobEnsemblePtr[j][k][i] = 0.0;
        }
      }
      _oobEnsembleDen[i] = 0;
    }
    if (_eventTypeSize > 1) {
      PROTECT(sexpVector[RSF_OPOE_ID] = NEW_NUMERIC(_eventTypeSize * obsSize));
      SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_OPOE_ID]);
      SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_OPOE_ID]));
      *p_oobEnsemblePOE = NUMERIC_POINTER(sexpVector[RSF_OPOE_ID]);
      sexpIndex ++;
      _oobPOEPtr  = pdvector(1, _eventTypeSize);
      for (j = 1; j <= _eventTypeSize; j++) {
          _oobPOEPtr[j] = (*p_oobEnsemblePOE) + ((j-1) * obsSize) - 1;
      }
      for (i = 1; i <= obsSize; i++) {
        for (j = 1; j <= _eventTypeSize; j++) {
            _oobPOEPtr[j][i] = 0.0;
        }
      }
      _oobSubSurvivalPtr     = dmatrix3(1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      _oobSubDistributionPtr = dmatrix3(1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        for (j = 1; j <= _eventTypeSize; j++) {
          for (k = 1; k <= _sortedTimeInterestSize; k++) {
            _oobSubSurvivalPtr[j][k][i] = 0.0;
            _oobSubDistributionPtr[j][k][i] = 0.0;
          }
        }
      }
    }  
  }
  if (_opt & OPT_PERF) {
    PROTECT(sexpVector[RSF_PERF_ID] = NEW_NUMERIC(performanceSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PERF_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_PERF_ID]));
    *p_performance = NUMERIC_POINTER(sexpVector[RSF_PERF_ID]);
    _performancePtr = pdvector(1, adjEventSize);
    for (i = 1; i <= adjEventSize; i++) {
      _performancePtr[i]  = (*p_performance)  + ((i-1) * _forestSize) - 1;
    }
    for (j = 1; j <= adjEventSize; j++) {
      for (k = 1; k <= _forestSize; k++) {
        _performancePtr[j][k] = NA_REAL;
      }
    }
    sexpIndex ++;
  }
  if (_opt & OPT_PROX) {
    PROTECT(sexpVector[RSF_PROX_ID] = NEW_INTEGER(proximitySize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PROX_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_PROX_ID]));
    *p_proximity = (uint*) INTEGER_POINTER(sexpVector[RSF_PROX_ID]);
    (*p_proximity) --;
    for (i = 1; i <= proximitySize; i++) {
      (*p_proximity)[i] = 0;
    }
    sexpIndex++;
  }
  if (_opt & OPT_LEAF) {
    PROTECT(sexpVector[RSF_LEAF_ID] = NEW_INTEGER(_forestSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_LEAF_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_LEAF_ID]));
    *p_leafCount    = (uint*) INTEGER_POINTER(sexpVector[RSF_LEAF_ID]);
    (*p_leafCount) --;
    for (i = 1; i <= _forestSize; i++) {
      (*p_leafCount)[i] = 0;
    }
    sexpIndex ++;
  }
  if (_opt & OPT_SEED) {
    PROTECT(sexpVector[RSF_SEED_ID] = NEW_INTEGER(1));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_SEED_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_SEED_ID]));
    *p_seed = INTEGER_POINTER(sexpVector[RSF_SEED_ID]);
    sexpIndex++;
  }
  if (_opt & OPT_TREE) {
    *p_root = nodePtrVector(1, _forestSize);
    for (i = 1; i <= _forestSize; i++) {
      (*p_root)[i] = NULL;
    }
  }
  if (_opt & OPT_MISS) {
    PROTECT(sexpVector[RSF_MISS_ID] = NEW_NUMERIC(imputationSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_MISS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_MISS_ID]));
    *p_imputation = NUMERIC_POINTER(sexpVector[RSF_MISS_ID]);
    *p_sImputePredictorPtr = pdvector(1, _xSize);
    *p_sImputeStatusPtr = (*p_imputation)  + (1 * mRecordSize) - 1;
    *p_sImputeTimePtr   = (*p_imputation)  + (2 * mRecordSize) - 1;
    for (i = 1; i <= _xSize; i++) {
      (*p_sImputePredictorPtr)[i]  = (*p_imputation)  + ((i+2) * mRecordSize) - 1;
    }
    for (i = 1; i <= mRecordSize; i++) {
      (*p_imputation)[i-1] = (double) mRecordIndex[i];
      (*p_sImputeStatusPtr)[i] = statusPtr[mRecordIndex[i]];
      (*p_sImputeTimePtr)[i]   = timePtr[mRecordIndex[i]];
      for (j = 1; j <= _xSize; j++) {
        (*p_sImputePredictorPtr)[j][i] = predictorPtr[j][mRecordIndex[i]];
      }
    }
    sexpIndex++;
  }
  if (_opt & OPT_OMIS) {
    PROTECT(sexpVector[RSF_OMIS_ID] = NEW_NUMERIC(imputationSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_OMIS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_OMIS_ID]));
    *p_oobImputation = NUMERIC_POINTER(sexpVector[RSF_OMIS_ID]);
    *p_sOOBImputePredictorPtr = pdvector(1, _xSize);
    *p_sOOBImputeStatusPtr = (*p_oobImputation)  + (1 * mRecordSize) - 1;
    *p_sOOBImputeTimePtr   = (*p_oobImputation)  + (2 * mRecordSize) - 1;
    for (i = 1; i <= _xSize; i++) {
      (*p_sOOBImputePredictorPtr)[i]  = (*p_oobImputation)  + ((i+2) * mRecordSize) - 1;
    }
    for (i = 1; i <= mRecordSize; i++) {
      (*p_oobImputation)[i-1] = (double) mRecordIndex[i];
      (*p_sOOBImputeStatusPtr)[i] = statusPtr[mRecordIndex[i]];
      (*p_sOOBImputeTimePtr)[i]   = timePtr[mRecordIndex[i]];
      for (j = 1; j <= _xSize; j++) {
        (*p_sOOBImputePredictorPtr)[j][i] = predictorPtr[j][mRecordIndex[i]];
      }
    }
    sexpIndex++;
  }
  if (_opt & OPT_VIMP) {
    PROTECT(sexpVector[RSF_VIMP_ID] = NEW_NUMERIC(importanceSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_VIMP_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_VIMP_ID]));
    *p_importance = NUMERIC_POINTER(sexpVector[RSF_VIMP_ID]);
    _importancePtr = pdvector(1, adjEventSize);
    for (i = 1; i <= adjEventSize; i++) {
      _importancePtr[i]  = (*p_importance)  + ((i-1) * importanceMult) - 1;
    }
    for (j = 1; j <= adjEventSize; j++) {
      for (k = 1; k <= importanceMult; k++) {
        _importancePtr[j][k] = NA_REAL;
      }
    }
    _vimpMortality = dmatrix(1, importanceMult, 1, obsSize);
    for (k=1; k <= importanceMult; k++) {
      for (i = 1; i <= obsSize; i++) {
        _vimpMortality[k][i] = 0.0;
      }
    }
    if (_opt & OPT_VOUT_TYPE) {
      _oobVimpInvalidDen  = uimatrix(1, importanceMult, 1, obsSize);
      for (k=1; k <= importanceMult; k++) {
      for (i = 1; i <= obsSize; i++) {
        _oobVimpInvalidDen[k][i] = 0;
      }
    }
    }
    if (_eventTypeSize > 1) {
      _crVimpEnsemble = dmatrix4(1, importanceMult, 1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      _crVimpPOE      = dmatrix3(1, importanceMult, 1, _eventTypeSize, 1, obsSize);
      for (p = 1; p <= importanceMult; p++) {
        for (j = 1; j <= _eventTypeSize; j++) {
          for (i = 1; i <= obsSize; i++) {
            _crVimpPOE[p][j][i] = 0.0;
            for (k = 1; k <= _sortedTimeInterestSize; k++) {
              _crVimpEnsemble[p][j][k][i] = 0.0;
            }
          }
        }
      }
    }
    sexpIndex++;
  }  
  if (_opt & OPT_VUSE) {
    PROTECT(sexpVector[RSF_VUSE_ID] = NEW_INTEGER(varUsedSize * _xSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_VUSE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_VUSE_ID]));
    *p_varUsed = (uint*) INTEGER_POINTER(sexpVector[RSF_VUSE_ID]);
    *p_varUsedPtr = puivector(1, varUsedSize);
    for (i = 1; i <= varUsedSize; i++) {
      (*p_varUsedPtr)[i] = (*p_varUsed) + ((i-1)*(_xSize)) - 1;
    }
    for (i = 1; i <= varUsedSize; i++) {
      for (j = 1; j <= _xSize; j++) {
        (*p_varUsedPtr)[i][j] = 0;
      }
    }
    sexpIndex++;
  }
  if (_opt & OPT_SPLT_DPTH) {
    PROTECT(sexpVector[RSF_DPTH_ID] = NEW_NUMERIC(_xSize * _observationSize));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_DPTH_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_DPTH_ID]));
    *p_splitDepth = NUMERIC_POINTER(sexpVector[RSF_DPTH_ID]);
    *p_splitDepthPtr  = pdvector(1, _xSize);
    for (j = 1; j <= _xSize; j++) {
      (*p_splitDepthPtr)[j] = (*p_splitDepth) + ((j-1)*(_observationSize)) - 1;
    }
    for (j = 1; j <= _xSize; j++) {
      for (i = 1; i <= _observationSize; i++) {
        (*p_splitDepthPtr)[j][i] = 0;
      }
    }
    *localSplitDepthPtr = dvector(1, _xSize);
    sexpIndex ++;
  }
  return (sexpIndex);
}
void unstackDefinedOutputObjects(char      mode,
                                 Node    **root,
                                 double   *localSplitDepthPtr) {
  uint obsSize;
  uint importanceMult;
  uint varUsedSize;
  uint adjEventSize;
  uint i, j;
  obsSize        = 0;  
  importanceMult = 0;  
  varUsedSize    = 0;  
  if (_eventTypeSize > 1) {
    adjEventSize = _eventTypeSize + 1;
  }
  else {
    adjEventSize = 1;
  }
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    if (_opt & OPT_VUSE) {
      if (_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
        varUsedSize = _forestSize;
      }
      else {
        varUsedSize = 1;
      }
    }
    if (_opt & OPT_VIMP) {
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceMult = 1;
      }
      else {
        importanceMult = _xSize;
      }
    }
    break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    if (_opt & OPT_VIMP) {
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceMult = 1;
      }
      else {
        importanceMult = _xSize;
      }
    }
    break;
  case RSF_INTR:
    obsSize = _fobservationSize;
    if (_opt & OPT_VIMP) {
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceMult = 1;
      }
      else {
        importanceMult = _intrPredictorSize;
      }
    }
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    error("\nRSF:  The application will now exit.\n");
    break;
  }
  if (_opt & OPT_FENS) {
    for (j = 1; j <= adjEventSize; j++) {
      free_pdvector(_fullEnsemblePtr[j], 1, _sortedTimeInterestSize);
    }
    free_ppdvector(_fullEnsemblePtr, 1, adjEventSize);
    free_uivector(_fullEnsembleDen, 1, obsSize);
    if (_eventTypeSize > 1) {
      free_pdvector(_fullPOEPtr, 1, _eventTypeSize);
      free_dmatrix3(_fullSubSurvivalPtr, 1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      free_dmatrix3(_fullSubDistributionPtr, 1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
    }
  }
  if (_opt & OPT_OENS) {
    for (j = 1; j <= adjEventSize; j++) {
      free_pdvector(_oobEnsemblePtr[j], 1, _sortedTimeInterestSize);
    }
    free_ppdvector(_oobEnsemblePtr, 1, adjEventSize);
    free_uivector(_oobEnsembleDen, 1, obsSize);
    if (_eventTypeSize > 1) {
      free_pdvector(_oobPOEPtr, 1, _eventTypeSize);
      free_dmatrix3(_oobSubSurvivalPtr, 1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      free_dmatrix3(_oobSubDistributionPtr, 1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
    }
  }
  if (_opt & OPT_PERF) {
    free_pdvector(_performancePtr, 1, adjEventSize);
  }
  if (_opt & OPT_TREE) {
    for (i = 1; i <= _forestSize; i++) {
      freeTree(root[i]);
    }
    free_nodePtrVector(root, 1, _forestSize);
  }
  if (_opt & OPT_MISS) {
    free_pdvector(_sImputePredictorPtr, 1, _xSize);
  }
  if (_opt & OPT_OMIS) {
    free_pdvector(_sOOBImputePredictorPtr, 1, _xSize);
  }
  if (_opt & OPT_VIMP) {
    free_pdvector(_importancePtr, 1, adjEventSize);
    free_dmatrix(_vimpMortality, 1, importanceMult, 1, obsSize);
    if (_opt & OPT_VOUT_TYPE) {
      free_uimatrix(_oobVimpInvalidDen, 1, importanceMult, 1, obsSize);
    }
    if (_eventTypeSize > 1) {
      free_dmatrix4(_crVimpEnsemble, 1, importanceMult, 1, _eventTypeSize, 1, _sortedTimeInterestSize, 1, obsSize);
      free_dmatrix3(_crVimpPOE, 1, importanceMult, 1, _eventTypeSize, 1, obsSize);
    }
  }
  if (_opt & OPT_VUSE) {
    free_puivector(_varUsedPtr, 1, varUsedSize);
  }  
  if (_opt & OPT_SPLT_DPTH) {
    free_pdvector(_splitDepthPtr, 1, _xSize);
    free_dvector(localSplitDepthPtr, 1, _xSize);
  }  
}
uint stackVariableOutputObjects(uint     totalNodeCount,
                                uint     totalMWCPCount,
                                uint   **p_treeID,
                                uint   **p_nodeID,
                                uint   **p_parmID,                                   
                                double **p_contPT,
                                uint   **p_mwcpSZ,
                                uint   **p_mwcpPT,
                                uint     sexpIndex,
                                char   **sexpString,
                                SEXP    *sexpVector) {
  if (_opt & OPT_TREE) {
    PROTECT(sexpVector[RSF_TREE_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RSF_NODE_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RSF_PARM_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RSF_CONT_PT] = NEW_NUMERIC(totalNodeCount));    
    PROTECT(sexpVector[RSF_MWCP_SZ] = NEW_INTEGER(totalNodeCount));  
    PROTECT(sexpVector[RSF_MWCP_PT] = NEW_INTEGER(totalMWCPCount));  
    *p_treeID = (uint*) INTEGER_POINTER(sexpVector[RSF_TREE_ID]);
    *p_nodeID = (uint*) INTEGER_POINTER(sexpVector[RSF_NODE_ID]);
    *p_parmID = (uint*) INTEGER_POINTER(sexpVector[RSF_PARM_ID]);
    *p_contPT = NUMERIC_POINTER(sexpVector[RSF_CONT_PT]);
    *p_mwcpSZ = (uint*) INTEGER_POINTER(sexpVector[RSF_MWCP_SZ]);
    *p_mwcpPT = (uint*) INTEGER_POINTER(sexpVector[RSF_MWCP_PT]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_TREE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_TREE_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_NODE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_NODE_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PARM_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_PARM_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_CONT_PT]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_CONT_PT]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_MWCP_SZ]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_MWCP_SZ]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_MWCP_PT]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_MWCP_PT]));
    (*p_treeID) --;
    (*p_nodeID) --;
    (*p_parmID) --;
    (*p_contPT) --;
    (*p_mwcpSZ) --;
    (*p_mwcpPT) --;
  }
  return (sexpIndex);
}
