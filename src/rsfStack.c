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
#include   "node_ops.h"
#include   "rsfFactorOps.h"
#include   "rsfImpute.h"
#include   "rsfUtil.h"
#include   "rsfStack.h"
extern uint getTraceFlag();
void stackPreDefinedCommonArrays() {
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedCommonArrays() ENTRY ...\n");
  }
  _nodeMembership = nodePtrVector(1, _observationSize);
  _bootMembershipIndex = uivector(1, _observationSize);
  _bootMembershipFlag = cvector(1, _observationSize);
  _masterTime  = dvector(1, _observationSize);
  _masterTimeIndex  = uivector(1, _observationSize);
  _oobSampleSize = uivector(1, _forestSize);
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedCommonArrays() EXIT ...\n");
  }
}
void unstackPreDefinedCommonArrays() {
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedCommonArrays() ENTRY ...\n");
  }
  free_nodePtrVector(_nodeMembership, 1, _observationSize);
  free_uivector(_bootMembershipIndex, 1, _observationSize);
  free_cvector(_bootMembershipFlag, 1, _observationSize);
  free_dvector(_masterTime, 1, _observationSize);
  free_uivector(_masterTimeIndex, 1, _observationSize);
  free_uivector(_oobSampleSize, 1, _forestSize);
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedCommonArrays() EXIT ...\n");
  }
}
void stackPreDefinedGrowthArrays() {
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedGrowthArrays() ENTRY ...\n");
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedGrowthArrays() EXIT ...\n");
  }
}
void unstackPreDefinedGrowthArrays() {
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedGrowthArrays() ENTRY ...\n");
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedGrowthArrays() EXIT ...\n");
  }
}
void stackPreDefinedPredictArrays() {
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedPredictArrays() ENTRY ...\n");
  }
  _fnodeMembership = nodePtrVector(1, _fobservationSize);
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedPredictArrays() EXIT ...\n");
  }
}
void unstackPreDefinedPredictArrays() {
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedPredictArrays() ENTRY ...\n");
  }
  free_nodePtrVector(_fnodeMembership, 1, _fobservationSize);
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedPredictArrays() EXIT ...\n");
  }
}
void stackPreDefinedInteractionArrays() {
  uint i;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedInteractionArrays() ENTRY ...\n");
  }
  _fnodeMembership = nodePtrVector(1, _fobservationSize);
  _foobSampleSize = uivector(1, _forestSize);
  _importanceFlag = cvector(1, _xSize);
  for (i = 1; i <= _xSize; i++) {
    _importanceFlag[i] = FALSE;
  }
  for (i = 1; i <= _intrPredictorSize; i++) {
    _importanceFlag[_intrPredictor[i]] = TRUE;
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackPreDefinedInteractionArrays() EXIT ...\n");
  }
}
void unstackPreDefinedInteractionArrays() {
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedInteractionArrays() ENTRY ...\n");
  }
  free_nodePtrVector(_fnodeMembership, 1, _fobservationSize);
  free_cvector(_importanceFlag, 1, _xSize);
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackPreDefinedInteractionArrays() EXIT ...\n");
  }
}
void initializeArrays(char mode, uint *sortedTimeInterestSize) {
  uint i,j;
  uint leadingIndex;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nCommon Incoming Parameters:  ");
    Rprintf("\n            traceFlag:  %10d", getTraceFlag());
    Rprintf("\n                  opt:  %10x", _opt);
    Rprintf("\n           forestSize:  %10d", _forestSize);
    Rprintf("\n      observationSize:  %10d", _observationSize);
    Rprintf("\n     timeInterestSize:  %10d", _timeInterestSize);
    Rprintf("\n                xSize:  %10d", _xSize);
    Rprintf("\n randomCovariateCount:  %10d", _randomCovariateCount);
    Rprintf("\n");
  }
  if ((mode == RSF_GROW) && (getTraceFlag() & SUMM_LOW_TRACE)) {
    Rprintf("\nIncoming Random Covariate Weights:  ");
    Rprintf("\n     index       weight");
    for (j=1; j <= _xSize; j++) {
      Rprintf("\n%10d  %12.4f", j, _randomCovariateWeight[j]);
    }
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nIncoming GROW Data:  ");
    Rprintf("\n       index       status         time   observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= _xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (j = 1; j <= _observationSize; j++) {
      Rprintf("%12d %12.4f %12.4f", j, _status[j], _time[j]);
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12.4f", (_xData+((i-1)*(_observationSize)))[j-1]);
      }
      Rprintf("\n");
    }
  }
  if ((mode == RSF_PRED) && (getTraceFlag() & SUMM_MED_TRACE)) {
    Rprintf("\nIncoming PRED Data:  ");
    Rprintf("\n       index       status         time   observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= _xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (j = 1; j <= _fobservationSize; j++) {
      Rprintf("%12d %12.4f %12.4f", j, _fstatus[j], _ftime[j]);
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12.4f", (_fxData+((i-1)*(_fobservationSize)))[j-1]);
      }
      Rprintf("\n");
    }
  }
  if ((mode == RSF_INTR) && (getTraceFlag() & SUMM_MED_TRACE)) {
    Rprintf("\nIncoming INTR Data (implied):  ");
    Rprintf("\n       index       status         time   observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= _xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (j = 1; j <= _fobservationSize; j++) {
      Rprintf("%12d %12.4f %12.4f", _intrObservation[j], _status[_intrObservation[j]], _time[_intrObservation[j]]);
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12.4f", (_xData+((i-1)*(_observationSize)))[_intrObservation[j]-1]);
      }
      Rprintf("\n");
    }
    Rprintf("\nIncoming INTR Predictors:  ");
    Rprintf("\n       index    predictor");
    for (i=1; i <= _intrPredictorSize; i++) {
      Rprintf("\n%12d %12d ", i, _intrPredictor[i]);
    }
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Initial read complete.");  
  }
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
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\n\nSorted Unique Event Times:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10.4f \n", i, _masterTime[i]);
    }
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Initialization of master time data complete.");  
  }
  hpsort(_timeInterest, _timeInterestSize);
  *sortedTimeInterestSize = 1;
  for (i=2; i <= _timeInterestSize; i++) {
    if (_timeInterest[i] > _timeInterest[*sortedTimeInterestSize]) {
      (*sortedTimeInterestSize) ++;
      _timeInterest[*sortedTimeInterestSize] = _timeInterest[i];
    }
  }
  if (*sortedTimeInterestSize != _timeInterestSize) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Time points of interest are not unique.");
    Rprintf("\nRSF:  The ensemble estimate output matrix is being");
    Rprintf("\nRSF:  resized as [N'] x [n], where N' is the");
    Rprintf("\nRSF:  unique time points of interest and n is");
    Rprintf("\nRSF:  number of observations in the data.");
  }
  for (i= (*sortedTimeInterestSize) + 1; i <= _timeInterestSize; i++) {
    _timeInterest[i] = 0;
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\n\nSorted Distinct Times of Interest:  \n");
    for (i=1; i <= *sortedTimeInterestSize; i++) {
      Rprintf("%10d %10.4f \n", i, _timeInterest[i]);
    }
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Initialization of time interest data complete.");  
  }
}
void initializeFactorArrays(char mode) {
  uint i, j;
  uint factorLevel;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\ninitializeFactorArrays() ENTRY ...\n");
  }
  if (_factorCount < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to initialize factorness in its absence.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
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
          Rprintf("\nRSF:  The application will now exit.\n");
          exit(TRUE);
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
            Rprintf("\nRSF:  The application will now exit.\n");
            exit(TRUE);
          }
        }
      }
      if (factorLevel > _factorSize[j]) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  !GROW factor level greater than maximum GROW factor level:  %10d vs. %10d", factorLevel, _factorSize[j]);
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
      }
    }
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nGROW Factor data:  ");
    Rprintf("\n     index    factors -> ");
    Rprintf("\n          ");
    for (j = 1; j <= _factorCount; j++) {
      Rprintf(" %10d", _factorIndex[j]);
    }
    Rprintf("\n\n");
    for (i = 1; i <= _observationSize; i++) {
      Rprintf("%10d", i);
      for (j = 1; j <= _factorCount; j++) {
        Rprintf(" %10.0f", _observation[_factorIndex[j]][i]);
      }
      Rprintf("\n");
    }
    Rprintf("\nMaximum factor levels:  ");
    Rprintf("\n          ");
    for (j = 1; j <= _factorCount; j++) {
      Rprintf(" %10d", _factorSize[j]);
    }
    Rprintf("\n");
    if (mode != RSF_GROW) {
      Rprintf("\n!GROW Factor data:  ");
      Rprintf("\n     index    factors -> ");
      Rprintf("\n          ");
      for (j = 1; j <= _factorCount; j++) {
        Rprintf(" %10d", _factorIndex[j]);
      }
      Rprintf("\n\n");
      for (i = 1; i <= _fobservationSize; i++) {
        Rprintf("%10d", i);
        for (j = 1; j <= _factorCount; j++) {
          Rprintf(" %10.0f", _fobservation[_factorIndex[j]][i]);
        }
        Rprintf("\n");
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\ninitializeFactorArrays() EXIT ...\n");
  }
}
void stackFactorArrays(char mode) {
  uint j, p;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackFactorArrays() ENTRY ...\n");
  }
  _xType = pcvector(1, _xSize);
  for (p = 1; p <= _xSize; p++) {
    _xType[p] = (char*) CHAR(STRING_ELT(AS_CHARACTER(_sexp_xType), p-1));
    if ((strcmp(_xType[p], "C") != 0) && (strcmp(_xType[p], "I") != 0) && (strcmp(_xType[p], "R") != 0)) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Invalid predictor type:  [%10d] = %2s", p, _xType[p]);
      Rprintf("\nRSF:  Type must be 'C', 'I', or 'R'.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
    }
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nIncoming Predictor Types:  ");
    Rprintf("\n     index  predictor \n");
    for (p = 1; p <= _xSize; p++) {
      Rprintf("%10d %10s \n", p, _xType[p]);
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
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nFactor Index Mapping:  ");
      Rprintf("\n     index  predictor \n");
      for (j = 1; j <= _factorCount; j++) {
        Rprintf("%10d %10d \n", j, _factorIndex[j]);
      }
    }
    _factorSize = uivector(1, _factorCount);
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackFactorArrays() EXIT ...\n");
  }
}
void unstackFactorArrays(char mode) {
  uint j;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackFactorArrays() ENTRY ...\n");
  }
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackFactorArrays() EXIT ...\n");
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackMissingSignatures() ENTRY ...\n");
  }
  if (recordSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to allocate for missingness in its absence.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
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
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nIndex of Individuals with any Missing Outcomes or Predictors:  ");
    Rprintf("\n    mIndex   orgIndex \n");
    for (i = 1; i <= recordSize; i++) {
      Rprintf("%10d %10d \n", i, (*p_recordIndex)[i]);
    }
    Rprintf("\nIncoming Indices of Missing Outcomes and Predictors:  ");
    Rprintf("\n   element      index \n");
    for (i = 1; i <= (*p_vSize); i++) {
      Rprintf("%10d %10d \n", i, (*p_vIndex)[i]);
    }
    Rprintf("\nIncoming Signatures of Missing Outcomes and Predictors:  ");
    Rprintf("\n     index   signatures -> \n");
    Rprintf(  "                  ");
    for (i=1; i <= _xSize; i++) {
      Rprintf("%3d", i);
    }
    Rprintf("\n");
    Rprintf(  "              C  T");
    for (i=1; i <= _xSize; i++) {
      Rprintf("%3s", _xType[i]);
    }
    Rprintf("\n");
    for (i = 1; i <= recordSize; i++) {
      Rprintf("%10d  ", (*p_recordIndex)[i]);
      for (j=1; j <= _xSize+2; j++) {
        Rprintf("%3d", (*p_vSign)[j][i]);
      }
      Rprintf("\n");
    }
    Rprintf("\nLinking of factor data structures to missing data structures:  ");
    Rprintf("\n       index   mFactorIdx \n");
    for (i = 1; i <= (*p_mFactorSize); i++) {
      Rprintf("%12d %12d \n", i, (*p_mFactorIndex)[i]);
    }
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackMissingSignatures() EXIT ...\n");
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackMissingSignatures() ENTRY ...\n");
  }
  free_dvector(statusPtr, 1, obsSize);
  free_dvector(timePtr, 1, obsSize);
  free_dmatrix(predictorPtr, 1, _xSize, 1, obsSize);
  free_uivector(recordIndex, 1, recordSize);
  free_imatrix(vSign, 1, _xSize+2, 1, recordSize);
  free_ivector(vIndex, 1, _xSize+2);
  free_imatrix(vForestSign, 1, _forestSize, 1, vSize);
  free_uivector(mFactorIndex, 1, _xSize);
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackMissingSignatures() EXIT ...\n");
  }
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackMissingArrays() ENTRY ...\n");
  }
  result = TRUE;
  for (i = 1 ; i <= _observationSize; i++) {
    if (!ISNA(_time[i])) {
      if (_time[i] <= 0) {
        result = FALSE;
        Rprintf("\nRSF:  GROW time elements must be greater than zero or NA:  [%10d] = %12.4f \n", i, _time[i]);
      }
    }
    if ( (_status[i] != 0) && (_status[i] != 1) && (!ISNA(_status[i])) ) {
      result = FALSE;
      Rprintf("\nRSF:  GROW status elements must be equal to one, zero, or NA:  [%10d] = %12.4f \n", i, _status[i]);
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
    Rprintf("\nRSF:  All GROW time elements are missing. \n");
    result = FALSE;
  }
  mFlag = FALSE;
  for (i = 1 ; i <= _observationSize; i++) {
    if (!ISNA(_status[i]) || (_status[i] == 1)) {
      mFlag = TRUE;
      i = _observationSize;
    }
  }
  if (mFlag == FALSE) {
    Rprintf("\nRSF:  All GROW status elements are are censored or missing. \n");
    result = FALSE;
  }
  if (result == FALSE) {
    return result;
  }
  if (mode == RSF_PRED) {
    if ((_opt & OPT_PERF) || (_opt & OPT_VIMP)) {
      for (i = 1 ; i <= _fobservationSize; i++) {
        if (!ISNA(_ftime[i])) {
          if (_ftime[i] <= 0) {
            result = FALSE;
            Rprintf("\nRSF:  PRED time elements must be greater than zero or NA when PERF or VIMP is requested:  [%10d] = %12.4f \n", i, _ftime[i]);
          }
        }
        if ( (_fstatus[i] != 0) && (_fstatus[i] != 1) && (!ISNA(_fstatus[i])) ) {
          result = FALSE;
          Rprintf("\nRSF:  PRED status elements must be equal to one, zero, or NA when PERF or VIMP is requested:  [%10d] = %12.4f \n", i, _fstatus[i]);
        }
      }
    }
    else {
      for (i = 1 ; i <= _fobservationSize; i++) {
        _fstatus[i] = _ftime[i] = -1;
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
    if (getTraceFlag() & SUMM_USR_TRACE) {
      Rprintf("\nRSF:  Missing data analysis of GROW complete -- none found.");  
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
    if (getTraceFlag() & SUMM_USR_TRACE) {
      Rprintf("\nRSF:  Missing data analysis of GROW complete -- some found.");  
    }
  }  
  if (mode == RSF_PRED) {
    _fmRecordMap = uivector(1, _fobservationSize);
    _fmRecordSize = getRecordMap(_fmRecordMap, 
                                 _fobservationSize, 
                                 _fstatus, 
                                 _ftime, 
                                 _fxData);
    if (_fmRecordSize == 0) {
      _fobservation = pdvector(1, _xSize);
      for (j=1; j <= _xSize; j++) {
        _fobservation[j] = (_fxData + ((j-1)*(_fobservationSize)) - 1);
      }
      if (getTraceFlag() & SUMM_USR_TRACE) {
        Rprintf("\nRSF:  Missing data analysis of PRED complete -- none found.");  
      }
    }  
    else {
      _fmStatus = dvector(1, _fobservationSize);
      _fmTime   = dvector(1, _fobservationSize);
      for (i = 1; i <= _fobservationSize; i++) {
        _fmStatus[i] = _fstatus[i];
        _fmTime[i] = _ftime[i];
      }
      _fstatus = _fmStatus;
      _ftime   = _fmTime;
      _fobservation = dmatrix(1, _xSize, 1, _fobservationSize);
      for (p=1; p <= _xSize; p++) {
        for (i = 1; i <= _fobservationSize; i++) {
          _fobservation[p][i] = (_fxData+((p-1)*(_fobservationSize)))[i-1];
        }
      }
      stackMissingSignatures(_fobservationSize,
                             _fstatus,
                             _ftime,
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
      if (getTraceFlag() & MISS_LOW_TRACE) {
        Rprintf("\nPRED Missing Signatures For Forest (cleared). \n");
        Rprintf(  " Outcome or Predictor: ");
        for (p=1; p <= _fmvSize; p++) {
          Rprintf("%3d", _fmvIndex[p]);
        }
        Rprintf("\n");
        for (i = 1; i <= _forestSize; i++) {
        Rprintf(  "            %10d ", i);
          for (p=1; p <= _fmvSize; p++) {
            Rprintf("%3d", _fmvForestSign[i][p]);
          }
          Rprintf("\n");
        }
      }
      if (getTraceFlag() & SUMM_USR_TRACE) {
        Rprintf("\nRSF:  Missing data analysis of PRED complete -- some found.");  
      }
    }  
  }  
  if (mode == RSF_INTR) {
    _fstatus = dvector(1, _fobservationSize);
    _ftime   = dvector(1, _fobservationSize);
    _fobservation = dmatrix(1, _xSize, 1, _fobservationSize);
    for (j = 1; j <= _fobservationSize; j++) {
      _ftime[j] = _time[_intrObservation[j]];
      _fstatus[j] = _status[_intrObservation[j]];
      for (i=1; i <= _xSize; i++) {
        _fobservation[i][j] = (_xData+((i-1)*(_observationSize)))[_intrObservation[j]-1];
      }
    }
    if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nIncoming INTR Data (actual):  ");
      Rprintf("\n       index       status         time   observations -> \n");
      Rprintf("\n                                      ");
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12d", i);
      }
      Rprintf("\n");
      for (j = 1; j <= _fobservationSize; j++) {
        Rprintf("%12d %12.4f %12.4f", j, _status[j], _time[j]);
        for (i=1; i <= _xSize; i++) {
          Rprintf(" %12.4f", (_fobservation[i][j]));
        }
        Rprintf("\n");
      }
    }
    if (getTraceFlag() & SUMM_USR_TRACE) {
      Rprintf("\nRSF:  Initialization of INTR data structures complete.");  
    }
    _fmRecordMap = uivector(1, _fobservationSize);
    _fmRecordSize = getRecordMap(_fmRecordMap, 
                                 _fobservationSize, 
                                 _fstatus, 
                                 _ftime, 
                                 &_fobservation[1][1]);
    if (_fmRecordSize == 0) {
    }
    else {
      stackMissingSignatures(_fobservationSize,
                             _fstatus,
                             _ftime,
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
      if (getTraceFlag() & MISS_LOW_TRACE) {
        Rprintf("\nINTR Missing Signatures For Forest (cleared). \n");
        Rprintf(  " Outcome or Predictor: ");
        for (p=1; p <= _fmvSize; p++) {
          Rprintf("%3d", _fmvIndex[p]);
        }
        Rprintf("\n");
        for (i = 1; i <= _forestSize; i++) {
        Rprintf(  "            %10d ", i);
          for (p=1; p <= _fmvSize; p++) {
            Rprintf("%3d", _fmvForestSign[i][p]);
          }
          Rprintf("\n");
        }
      }
      if (getTraceFlag() & SUMM_USR_TRACE) {
        Rprintf("\nRSF:  Missing data analysis of INTR complete -- some found.");  
      }
    }  
  }  
  dualUseFlag = FALSE;
  switch (mode) {
  case RSF_GROW:
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
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackMissingArrays() EXIT ...\n");
  }
  return result;
}
void unstackMissingArrays(char   mode,
                          char  **dmRecordBootFlag,
                          double ***dmvImputation) {
  char dualUseFlag;
  uint recordSize;
  uint vSize;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackMissingArrays() ENTRY ...\n");
  }
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
                               _ftime,
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
      free_pdvector(_fobservation, 1, _xSize);
    }
    else {
      unstackMissingSignatures(_fobservationSize,
                               _fstatus,
                               _ftime,
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
  case RSF_GROW:
    if (_mRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = _mRecordSize;
      vSize = _mvSize;
    }
    break;
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
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }  
  if (dualUseFlag == TRUE) {
    free_cmatrix(dmRecordBootFlag, 1, _forestSize, 1, recordSize);
    free_dmatrix3(dmvImputation, 1, _forestSize, 1, recordSize, 1, vSize);
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackMissingArrays() EXIT ...\n");
  }
}
uint stackDefinedOutputObjects(char      mode,
                               uint      sortedTimeInterestSize,
                               char    **sexpString,
                               Node   ***p_root,
                               double  **p_oobEnsemble,
                               double  **p_fullEnsemble,
                               double  **p_performance,
                               uint    **p_leafCount,
                               uint    **p_proximity,
                               double  **p_importance,
                               int     **p_seed,
                               double  **p_imputation,
                               double  **p_oobImputation,
                               double  **p_sImputeStatusPtr,
                               double  **p_sImputeTimePtr,
                               double ***p_sImputePredictorPtr,
                               double  **p_sOOBImputeStatusPtr,
                               double  **p_sOOBImputeTimePtr,
                               double ***p_sOOBImputePredictorPtr,
                               uint    **p_varUsed,
                               uint   ***p_varUsedPtr,
                               uint     *stackCount,
                               SEXP     *sexpVector) {
  uint sexpIndex;
  uint ensembleSize;
  uint performanceSize;
  uint proximitySize;
  uint imputationSize;
  uint importanceSize;
  uint varUsedSize;
  uint  obsSize;
  uint  mRecordSize;
  uint *mRecordIndex;
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  uint i,j,p;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackDefinedOutputObjects() ENTRY ...\n");
  }
  sexpIndex      = 0;  
  ensembleSize   = 0;  
  performanceSize= 0;  
  proximitySize  = 0;  
  imputationSize = 0;  
  importanceSize = 0;  
  varUsedSize    = 0;  
  obsSize        = 0;  
  mRecordSize    = 0;  
  statusPtr      = NULL;  
  timePtr        = NULL;  
  predictorPtr   = NULL;  
  mRecordIndex   = NULL;  
  _ensembleRun = NULL;
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    mRecordSize = _mRecordSize;
    performanceSize = _forestSize;
    statusPtr = _status;
    timePtr = _time;
    predictorPtr = _observation;
    mRecordIndex = _mRecordIndex;
    *stackCount = 4;
    if (_opt & (OPT_POUT_TYPE)) {
      ensembleSize = 1;
    }
    else {
      ensembleSize = sortedTimeInterestSize;
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
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
      }
    }
    if ((_opt & OPT_MISS) && (_opt & OPT_OMIS)) {
      imputationSize = (_xSize + 3) * mRecordSize;
      (*stackCount) += 2;
    }
    if (_opt & OPT_VIMP) {
      importanceSize = _xSize;
      (*stackCount) += 1;
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
    break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    mRecordSize = _fmRecordSize;
    performanceSize = _forestSize;
    statusPtr = _fstatus;
    timePtr = _ftime;
    predictorPtr = _fobservation;
    mRecordIndex = _fmRecordIndex;
    *stackCount = 2;
    if (_opt & (OPT_POUT_TYPE)) {
      ensembleSize = 1;
    }
    else {
      ensembleSize = sortedTimeInterestSize;
    }
    if (_opt & OPT_PERF) {
      (*stackCount) += 1;
    }
    if (_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2; 
      (*stackCount) += 1;
    }
    if (_opt & OPT_VIMP) {
      importanceSize = _xSize;
      (*stackCount) += 1;
    }
    if (_opt & OPT_MISS) {
      imputationSize = (_xSize + 3) * mRecordSize;
      (*stackCount) += 1;
    }
    break;
  case RSF_INTR:
    obsSize = _fobservationSize;
    mRecordSize = _fmRecordSize;
    performanceSize = _forestSize;
    statusPtr = _fstatus;
    timePtr = _ftime;
    predictorPtr = _fobservation;
    mRecordIndex = _fmRecordIndex;
    *stackCount = 3;
    if (_opt & OPT_VIMP) {
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceSize = 1;
      }
      else {
        importanceSize = _intrPredictorSize;
      }
      (*stackCount) += 1;
    }
    if (_opt & (OPT_POUT_TYPE)) {
      ensembleSize = 1;
    }
    else {
      ensembleSize = sortedTimeInterestSize;
    }
    if (_opt & OPT_OMIS) {
      imputationSize = (_xSize + 3) * mRecordSize;
      (*stackCount) += 1;
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
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nTotal Stack Count:  %12d", *stackCount);  
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nSEXP VECTOR mapped at: %20x", sexpVector);  
  }    
  PROTECT(sexpVector[RSF_OUTP_ID] = allocVector(VECSXP, *stackCount));
  PROTECT(sexpVector[RSF_STRG_ID] = allocVector(STRSXP, *stackCount));
  setAttrib(sexpVector[RSF_OUTP_ID], R_NamesSymbol, sexpVector[RSF_STRG_ID]);
  sexpIndex = 0;
  if (_opt & OPT_FENS) {
    PROTECT(sexpVector[RSF_FENS_ID] = NEW_NUMERIC(ensembleSize * obsSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nFENS memory mapped at: %20x", sexpVector + RSF_FENS_ID);  
      Rprintf("\nFENS memory:           %20x", sexpVector[RSF_FENS_ID]);  
      Rprintf("\nFENS memory sized:     %20d", ensembleSize * obsSize);
      Rprintf("\nFENS SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_FENS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_FENS_ID]));
    *p_fullEnsemble = NUMERIC_POINTER(sexpVector[RSF_FENS_ID]);
    _fullEnsemblePtr = pdvector(1, ensembleSize);
    _fullEnsembleDen = uivector(1, obsSize);
    if (_ensembleRun == NULL) {
      _ensembleRun     = dvector(1, obsSize);
    }
    for (i = 1; i <= ensembleSize; i++) {
      _fullEnsemblePtr[i] = (*p_fullEnsemble) + ((i-1)*(obsSize)) - 1;
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= ensembleSize; j++) {
        _fullEnsemblePtr[j][i] = 0.0;
      }
      _fullEnsembleDen[i] = 0;
    }
    sexpIndex ++;
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT FENS complete.");  
    }
  }
  if (_opt & OPT_OENS) {
    PROTECT(sexpVector[RSF_OENS_ID] = NEW_NUMERIC(ensembleSize * obsSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nOENS memory mapped at: %20x", sexpVector + RSF_OENS_ID);  
      Rprintf("\nOENS memory:           %20x", sexpVector[RSF_OENS_ID]);  
      Rprintf("\nOENS memory sized:     %20d", ensembleSize * obsSize);
      Rprintf("\nOENS SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_OENS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_OENS_ID]));
    *p_oobEnsemble = NUMERIC_POINTER(sexpVector[RSF_OENS_ID]);
    _oobEnsemblePtr  = pdvector(1, ensembleSize);
    _oobEnsembleDen  = uivector(1, obsSize);
    if (_ensembleRun == NULL) {
      _ensembleRun     = dvector(1, obsSize);
    }
    for (i = 1; i <= ensembleSize; i++) {
      _oobEnsemblePtr[i]  = (*p_oobEnsemble)  + ((i-1)*(obsSize)) - 1;
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= ensembleSize; j++) {
        _oobEnsemblePtr[j][i]  = 0.0;
      }
      _oobEnsembleDen[i]  = 0;
    }
    sexpIndex ++;
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT OENS complete.");  
    }
  }
  if (_opt & OPT_PERF) {
    PROTECT(sexpVector[RSF_PERF_ID] = NEW_NUMERIC(performanceSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nPERF memory mapped at: %20x", sexpVector + RSF_PERF_ID);  
      Rprintf("\nPERF memory:           %20x", sexpVector[RSF_PERF_ID]);  
      Rprintf("\nPERF memory sized:     %20d", performanceSize);  
      Rprintf("\nPERF SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PERF_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_PERF_ID]));
    *p_performance  = NUMERIC_POINTER(sexpVector[RSF_PERF_ID]);
    (*p_performance) --;
    for (i = 1; i <= performanceSize; i++) {
      (*p_performance)[i] = NA_REAL;
    }
    sexpIndex ++;
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT PERF complete.");  
    }
  }
  if (_opt & OPT_PROX) {
    PROTECT(sexpVector[RSF_PROX_ID] = NEW_INTEGER(proximitySize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nPROX memory mapped at: %20x", sexpVector + RSF_PROX_ID);  
      Rprintf("\nPROX memory:           %20x", sexpVector[RSF_PROX_ID]);  
      Rprintf("\nPROX memory sized:     %20d", proximitySize);  
      Rprintf("\nPROX SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PROX_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_PROX_ID]));
    *p_proximity = (uint*) INTEGER_POINTER(sexpVector[RSF_PROX_ID]);
    (*p_proximity) --;
    for (i = 1; i <= proximitySize; i++) {
      (*p_proximity)[i] = 0;
    }
    sexpIndex++;
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT PROX complete.");  
    }
  }
  if (_opt & OPT_LEAF) {
    PROTECT(sexpVector[RSF_LEAF_ID] = NEW_INTEGER(_forestSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nLEAF memory mapped at: %20x", sexpVector + RSF_LEAF_ID);  
      Rprintf("\nLEAF memory:           %20x", sexpVector[RSF_LEAF_ID]);  
      Rprintf("\nLEAF memory sized:     %20d", _forestSize);  
      Rprintf("\nLEAF SEXP index:       %20d", sexpIndex);  
    }
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_LEAF_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_LEAF_ID]);
    *p_leafCount    = (uint*) INTEGER_POINTER(sexpVector[RSF_LEAF_ID]);
    (*p_leafCount) --;
    for (i = 1; i <= _forestSize; i++) {
      (*p_leafCount)[i] = 0;
    }
    sexpIndex ++;
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT LEAF complete.");  
    }
  }
  if (_opt & OPT_SEED) {
    PROTECT(sexpVector[RSF_SEED_ID] = NEW_INTEGER(1));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nSEED memory mapped at: %20x", sexpVector + RSF_SEED_ID);  
      Rprintf("\nSEED memory:           %20x", sexpVector[RSF_SEED_ID]);  
      Rprintf("\nSEED memory sized:     %20d", 1);  
      Rprintf("\nSEED SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_SEED_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_SEED_ID]));
    *p_seed = INTEGER_POINTER(sexpVector[RSF_SEED_ID]);
    sexpIndex++;
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT SEED complete.");  
    }
  }
  if (_opt & OPT_TREE) {
    *p_root = nodePtrVector(1, _forestSize);
    for (i = 1; i <= _forestSize; i++) {
      (*p_root)[i] = NULL;
    }
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT TREE complete.");  
    }
  }
  if (_opt & OPT_MISS) {
    PROTECT(sexpVector[RSF_MISS_ID] = NEW_NUMERIC(imputationSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nMISS memory mapped at: %20x", sexpVector + RSF_MISS_ID);  
      Rprintf("\nMISS memory:           %20x", sexpVector[RSF_MISS_ID]);  
      Rprintf("\nMISS memory sized:     %20d", imputationSize);  
      Rprintf("\nMISS SEXP index:       %20d", sexpIndex);  
    }
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
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT MISS complete.");  
    }
  }
  if (_opt & OPT_OMIS) {
    PROTECT(sexpVector[RSF_OMIS_ID] = NEW_NUMERIC(imputationSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nOMIS memory mapped at: %20x", sexpVector + RSF_OMIS_ID);  
      Rprintf("\nOMIS memory:           %20x", sexpVector[RSF_OMIS_ID]);  
      Rprintf("\nOMIS memory sized:     %20d", imputationSize);  
      Rprintf("\nOMIS SEXP index:       %20d", sexpIndex);  
    }
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
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT OMIS complete.");  
    }
  }
  if (_opt & OPT_VIMP) {
    PROTECT(sexpVector[RSF_VIMP_ID] = NEW_NUMERIC(importanceSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nVIMP memory mapped at: %20x", sexpVector + RSF_VIMP_ID);  
      Rprintf("\nVIMP memory:           %20x", sexpVector[RSF_VIMP_ID]);  
      Rprintf("\nVIMP memory sized:     %20d", importanceSize);  
      Rprintf("\nVIMP SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_VIMP_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_VIMP_ID]));
    *p_importance = NUMERIC_POINTER(sexpVector[RSF_VIMP_ID]);
    (*p_importance) --;
    for (i = 1; i <= importanceSize; i++) {
      (*p_importance)[i] = NA_REAL;
    }
    _vimpEnsembleRun = dmatrix(1, importanceSize, 1, obsSize);
    for (i = 1; i <= importanceSize; i++) {
      (*p_importance)[i] = NA_REAL;
      for (j = 1; j <= obsSize; j++) {
        _vimpEnsembleRun[i][j] = 0;
      }
    }
    sexpIndex++;
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT VIMP complete.");  
    }
  }
  if (_opt & OPT_VUSE) {
    PROTECT(sexpVector[RSF_VUSE_ID] = NEW_INTEGER(varUsedSize * _xSize));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nVUSE memory mapped at: %20x", sexpVector + RSF_VUSE_ID);  
      Rprintf("\nVUSE memory:           %20x", sexpVector[RSF_VUSE_ID]);  
      Rprintf("\nVUSE memory sized:     %20d", varUsedSize * _xSize);  
      Rprintf("\nVUSE SEXP index:       %20d", sexpIndex);  
    }
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
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nAllocating for OPT VUSE complete.");  
    }
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Allocation of defined output objects complete:  %10d", sexpIndex);  
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    if (_opt & OPT_MISS) {
      Rprintf("\nImputed Data Output Object:  (at initialization)");
      Rprintf("\n       index   imputation -> \n");
      Rprintf(  "             %12d %12d", CENS_IDX, TIME_IDX);
      for (p=1; p <= _xSize; p++) {
        Rprintf(" %12d", p);
      }
      Rprintf("\n");
      for (i = 1; i <= mRecordSize; i++) {
        Rprintf("%12d", mRecordIndex[i]);
        Rprintf(" %12.4f", (*p_sImputeStatusPtr)[i]);
        Rprintf(" %12.4f", (*p_sImputeTimePtr)[i]);
        for (p = 1; p <= _xSize; p++) {
          Rprintf(" %12.4f", (*p_sImputePredictorPtr)[p][i]);
        }
        Rprintf("\n");
      }
    }
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackDefinedOutputObjects() EXIT ...\n");
  }
  return (sexpIndex);
}
void unstackDefinedOutputObjects(char      mode,
                                 uint      sortedTimeInterestSize,
                                 Node    **root) {
  uint  obsSize;
  uint mRecordSize;
  uint importanceSize;
  uint varUsedSize;
  uint i;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackDefinedOutputObjects() ENTRY ...\n");
  }
  obsSize        = 0;  
  mRecordSize    = 0;  
  importanceSize = 0;  
  varUsedSize    = 0;  
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    mRecordSize = _mRecordSize;
    if (_opt & OPT_VIMP) {
      importanceSize = _xSize;
    }
    if (_opt & OPT_VUSE) {
      if (_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
        varUsedSize = _forestSize;
      }
      else {
        varUsedSize = 1;
      }
    }
    break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    mRecordSize = _fmRecordSize;
    if (_opt & OPT_VIMP) {
      importanceSize = _xSize;
    }
    break;
  case RSF_INTR:
    obsSize = _observationSize;
    mRecordSize = _fmRecordSize;
    if (_opt & OPT_VIMP) {
      if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
        importanceSize = 1;
      }
      else {
        importanceSize = _intrPredictorSize;
      }
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
  if (_opt & OPT_FENS) {
    free_pdvector(_fullEnsemblePtr, 1, sortedTimeInterestSize);
    free_uivector(_fullEnsembleDen, 1, obsSize);
    if (_ensembleRun != NULL) {
      free_dvector(_ensembleRun, 1, obsSize);
      _ensembleRun = NULL;
    }
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nDe-allocating for OPT FENS complete.");  
    }
  }
  if (_opt & OPT_OENS) {
    free_uivector(_oobEnsembleDen, 1, obsSize);
    free_pdvector(_oobEnsemblePtr, 1, sortedTimeInterestSize);
    if (_ensembleRun != NULL) {
      free_dvector(_ensembleRun, 1, obsSize);
      _ensembleRun = NULL;
    }
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nDe-allocating for OPT OENS complete.");  
    }
  }
  if (_opt & OPT_TREE) {
    for (i = 1; i <= _forestSize; i++) {
      freeTree(root[i]);
    }
    free_nodePtrVector(root, 1, _forestSize);
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nDe-allocating for OPT TREE complete.");  
    }
  }
  if (_opt & OPT_MISS) {
    free_pdvector(_sImputePredictorPtr, 1, _xSize);
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nDe-allocating for OPT MISS complete.");  
    }
  }
  if (_opt & OPT_OMIS) {
    free_pdvector(_sOOBImputePredictorPtr, 1, _xSize);
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nDe-allocating for OPT OMIS complete.");  
    }
  }
  if (_opt & OPT_VIMP) {
    free_dmatrix(_vimpEnsembleRun, 1, importanceSize, 1, obsSize);
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nDe-allocating for OPT VIMP complete.");  
    }
  }
  if (_opt & OPT_VUSE) {
    free_puivector(_varUsedPtr, 1, varUsedSize);
    if (getTraceFlag() & SUMM_LOW_TRACE) {
      Rprintf("\nDe-allocating for OPT VUSE complete.");  
    }
  }  
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nunstackDefinedOutputObjects() EXIT ...\n");
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackVariableOutputObjects() ENTRY ...\n");
  }
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
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nTREE memory mapped at: %20x", sexpVector + RSF_TREE_ID);  
      Rprintf("\nTREE memory:           %20x", sexpVector[RSF_TREE_ID]);  
      Rprintf("\nTREE memory sized:     %20d", totalNodeCount);  
      Rprintf("\nTREE SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_TREE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_TREE_ID]));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nNODE memory mapped at: %20x", sexpVector + RSF_NODE_ID);  
      Rprintf("\nNODE memory:           %20x", sexpVector[RSF_NODE_ID]);  
      Rprintf("\nNODE memory sized:     %20d", totalNodeCount);  
      Rprintf("\nNODE SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_NODE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_NODE_ID]));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nPARM memory mapped at: %20x", sexpVector + RSF_PARM_ID);  
      Rprintf("\nPARM memory:           %20x", sexpVector[RSF_PARM_ID]);  
      Rprintf("\nPARM memory sized:     %20d", totalNodeCount);  
      Rprintf("\nPARM SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PARM_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_PARM_ID]));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nCONT memory mapped at: %20x", sexpVector + RSF_CONT_PT);  
      Rprintf("\nCONT memory:           %20x", sexpVector[RSF_CONT_PT]);  
      Rprintf("\nCONT memory sized:     %20d", totalNodeCount);  
      Rprintf("\nCONT SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_CONT_PT]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_CONT_PT]));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nMWSZ memory mapped at: %20x", sexpVector + RSF_MWCP_SZ);  
      Rprintf("\nMWSZ memory:           %20x", sexpVector[RSF_MWCP_SZ]);  
      Rprintf("\nMWSZ memory sized:     %20d", totalNodeCount);  
      Rprintf("\nMWSZ SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_MWCP_SZ]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_MWCP_SZ]));
    if (getTraceFlag() & SUMM_MED_TRACE) {
      Rprintf("\nMWPT memory mapped at: %20x", sexpVector + RSF_MWCP_PT);  
      Rprintf("\nMWPT memory:           %20x", sexpVector[RSF_MWCP_PT]);  
      Rprintf("\nMWPT memory sized:     %20d", totalMWCPCount);  
      Rprintf("\nMWPT SEXP index:       %20d", sexpIndex);  
    }
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_MWCP_PT]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_MWCP_PT]));
    (*p_treeID) --;
    (*p_nodeID) --;
    (*p_parmID) --;
    (*p_contPT) --;
    (*p_mwcpSZ) --;
    (*p_mwcpPT) --;
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Allocation of variable output objects complete:  %10d", sexpIndex);  
  }
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nstackVariableOutputObjects() EXIT ...\n");
  }
  return (sexpIndex);
}
