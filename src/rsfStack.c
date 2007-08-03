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

#include   "global.h"
#include   "nrutil.h"
#include   "node_ops.h"
#include   "rsfImpute.h"
#include   "rsfUtil.h"
#include   "rsfStack.h"
extern uint getTraceFlag();
void stackPreDefinedCommonArrays(uint **p_oobSampleSize) {
  uint i;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackPreDefinedCommonArrays() ENTRY ...\n");
  }
  _nodeMembership = nodePtrVector(1, _observationSize);
  _bootMembershipIndex = uivector(1, _observationSize);
  _bootMembershipFlag = cvector(1, _observationSize);
  _masterTime  = dvector(1, _observationSize);
  _masterTimeIndex  = uivector(1, _observationSize);
  *p_oobSampleSize = uivector(1, _forestSize);
  for (i = 1; i <= _forestSize; i++) {
    (*p_oobSampleSize)[i] = 0;
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackPreDefinedCommonArrays() EXIT ...\n");
  }
}
void unstackPreDefinedCommonArrays(uint *oobSampleSize) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackPreDefinedCommonArrays() ENTRY ...\n");
  }
  free_nodePtrVector(_nodeMembership, 1, _observationSize);
  free_uivector(_bootMembershipIndex, 1, _observationSize);
  free_cvector(_bootMembershipFlag, 1, _observationSize);
  free_dvector(_masterTime, 1, _observationSize);
  free_uivector(_masterTimeIndex, 1, _observationSize);
  free_uivector(oobSampleSize, 1, _forestSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackPreDefinedCommonArrays() EXIT ...\n");
  }
}
void stackPreDefinedPredictArrays() {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackPreDefinedPredictArrays() ENTRY ...\n");
  }
  _fnodeMembership = nodePtrVector(1, _fobservationSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackPreDefinedPredictArrays() EXIT ...\n");
  }
}
void unstackPreDefinedPredictArrays() {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackPreDefinedPredictArrays() ENTRY ...\n");
  }
  free_nodePtrVector(_fnodeMembership, 1, _fobservationSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackPreDefinedPredictArrays() EXIT ...\n");
  }
}
void stackPreDefinedGrowthArrays(double ***p_masterSplit,
                                 uint **p_masterSplitSize,
                                 uint ***p_masterSplitOrder) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackPreDefinedGrowthArrays() ENTRY ...\n");
  }
  *p_masterSplit = dmatrix(1, _xSize, 1, _observationSize);
  *p_masterSplitSize = uivector(1, _xSize);
  if (_splitRule == LOG_RANK_SCORE) {
    *p_masterSplitOrder = uimatrix(1, _xSize, 1, _observationSize);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackPreDefinedGrowthArrays() EXIT ...\n");
  }
}
void unstackPreDefinedGrowthArrays(double **masterSplit,
                                   uint *masterSplitSize,
                                   uint **masterSplitOrder) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackPreDefinedGrowthArrays() ENTRY ...\n");
  }
  free_dmatrix(masterSplit, 1, _xSize, 1, _observationSize);
  free_uivector(masterSplitSize, 1, _xSize);
  if (_splitRule == LOG_RANK_SCORE) {
    free_uimatrix(masterSplitOrder, 1, _xSize, 1, _observationSize);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackPreDefinedGrowthArrays() EXIT ...\n");
  }
}
char stackAndInitializeMissingArrays(uint       mode,
                                     char      *p_mTimeIndexFlag,
                                     char    ***p_mRecordBootFlag,
                                     double ****p_mvImputation) {
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  char result;
  char mFlag;
  uint i,j,k,p;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackAndInitializeMissingArrays() ENTRY ...\n");
  }
  result = TRUE;
  *p_mTimeIndexFlag = FALSE;
  _xType = pcvector(1, _xSize);
  for (p = 1; p <= _xSize; p++) {
    _xType[p] = (char*) CHAR(STRING_ELT(AS_CHARACTER(_sexp_xType), p-1));
    if ((strcmp(_xType[p], "I") != 0) && (strcmp(_xType[p], "R") != 0)) {
      Rprintf("\nRSF:  Invalid predictor type:  [%10d] = %2s", p, _xType[p]);
      Rprintf("\nRSF:  Type must be 'I', or 'R'.");
      result = FALSE;
    }
  }
  if (result == FALSE) {
    return result;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nIncoming Predictor Types:  ");
    Rprintf("\n     index  predictor \n");
    for (p = 1; p <= _xSize; p++) {
      Rprintf("%10d %10s \n", p, _xType[p]);
    }
  }
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
    if (_mup & MUP_PERF) {
      for (i = 1 ; i <= _fobservationSize; i++) {
        if (!ISNA(_ftime[i])) {
          if (_ftime[i] <= 0) {
            result = FALSE;
            Rprintf("\nRSF:  PRED time elements must be greater than zero or NA:  [%10d] = %12.4f \n", i, _ftime[i]);
          }
        }
        if ( (_fstatus[i] != 0) && (_fstatus[i] != 1) && (!ISNA(_fstatus[i])) ) {
          result = FALSE;
          Rprintf("\nRSF:  PRED status elements must be equal to one, zero, or NA:  [%10d] = %12.4f \n", i, _fstatus[i]);
        }
      }
    }
    else {
      for (i = 1 ; i <= _fobservationSize; i++) {
        _fstatus[i] = _ftime[i] = NA_REAL;
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
    _observation = pdvector(1, _xSize);
    for (i=1; i <= _xSize; i++) {
      _observation[i] = (_xData + ((i-1)*(_observationSize)) - 1);
    }
    if (getTraceFlag() & DL0_TRACE) {
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
    _mRecordIndex = uivector(1, _mRecordSize);
    i = 0;
    for (j = 1; j <= _observationSize; j++) {
      if (_mRecordMap[j] > 0) {
        i++;
        _mRecordIndex[i] = j;
      }
    }
    _mvSign = imatrix(1, _xSize+2, 1, _mRecordSize);
    for (j = 1; j <= _mRecordSize; j++) {
      for (i = 1; i <= _xSize+2; i++) {
        _mvSign[i][j] = 0;
      }
    }
    statusPtr = _status;
    timePtr = _time;
    predictorPtr = _observation;
    for (j = 1; j <= _mRecordSize; j++) {
      if (ISNA(statusPtr[_mRecordIndex[j]])) {
        _mvSign[(uint) abs(CENS_IDX)][j] = 1;
      }
      if (ISNA(timePtr[_mRecordIndex[j]])) {
        _mvSign[(uint) abs(TIME_IDX)][j] = 1;
      }
      for (i = 1; i <= _xSize; i++) {
        if (ISNA(predictorPtr[i][_mRecordIndex[j]])) {
          _mvSign[i+2][j] = 1;
        }
      }
    }
    _mvIndex = ivector(1, _xSize+2);
    _mvSize = 0;
    for (i = 1; i <= _xSize+2; i++) {
      _mvIndex[i] = 0;
      switch (i) {
      case (uint) (-CENS_IDX):
        for (j = 1; j <= _mRecordSize; j++) {
          if (_mvSign[i][j] == 1) {
            _mvSize ++;
            _mvIndex[_mvSize] = CENS_IDX;
            j = _mRecordSize;
          }
        }
        break;
      case (uint) (-TIME_IDX):
        for (j = 1; j <= _mRecordSize; j++) {
          if (_mvSign[i][j] == 1) {
            _mvSize ++;
            _mvIndex[_mvSize] = TIME_IDX;
            *p_mTimeIndexFlag = TRUE;
            j = _mRecordSize;
          }
        }
        break;
      default:
        for (j = 1; j <= _mRecordSize; j++) {
          if (_mvSign[i][j] == 1) {
            _mvSize ++;
            _mvIndex[_mvSize] = i-2;
            j = _mRecordSize;
          }
        }
        break;
      }
    }  
    _mvForestSign = imatrix(1, _forestSize,  1, _mvSize);
    for (j = 1; j <= _forestSize; j++) {
      for (p = 1; p <= _mvSize; p++) {
        _mvForestSign[j][p] = 0;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nMapping of GROW Individuals to Missing Record Index:  ");
      Rprintf("\n    record      index \n");
      for (i = 1; i <= _observationSize; i++) {
        if (_mRecordMap[i] != 0) {
          Rprintf("%10d %10d \n", i, _mRecordMap[i]);
        }
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nIndex of GROW Individuals with any Missing Outcomes or Predictors:  ");
      Rprintf("\n   element      index \n");
      for (i = 1; i <= _mRecordSize; i++) {
          Rprintf("%10d %10d \n", i, _mRecordIndex[i]);
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nIncoming Indices of Missing GROW Outcomes and Predictors:  ");
      Rprintf("\n   element      index \n");
      for (i = 1; i <= _mvSize; i++) {
        Rprintf("%10d %10d \n", i, _mvIndex[i]);
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nIncoming Signatures of Missing GROW Outcomes and Predictors:  ");
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
    }
    if (getTraceFlag() & DL2_TRACE) {
      for (i = 1; i <= _mRecordSize; i++) {
        Rprintf("%10d  ", _mRecordIndex[i]);
        for (j=1; j <= _xSize+2; j++) {
          Rprintf("%3d", _mvSign[j][i]);
        }
        Rprintf("\n");
      }
    }
    if (getTraceFlag() & DL0_TRACE) {
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
      if (getTraceFlag() & DL0_TRACE) {
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
      _fmRecordIndex = uivector(1, _fmRecordSize);
      i = 0;
      for (j = 1; j <= _fobservationSize; j++) {
        if (_fmRecordMap[j] > 0) {
          i++;
          _fmRecordIndex[i] = j;
        }
      }
      _fmvSign = imatrix(1, _xSize+2, 1, _fmRecordSize);
      for (j = 1; j <= _fmRecordSize; j++) {
        for (i = 1; i <= _xSize+2; i++) {
          _fmvSign[i][j] = 0;
        }
      }
      statusPtr = _fstatus;
      timePtr = _ftime;
      predictorPtr = _fobservation;
      for (j = 1; j <= _fmRecordSize; j++) {
        if (ISNA(statusPtr[_fmRecordIndex[j]])) {
          _fmvSign[(uint) abs(CENS_IDX)][j] = 1;
        }
        if (ISNA(timePtr[_fmRecordIndex[j]])) {
          _fmvSign[(uint) abs(TIME_IDX)][j] = 1;
        }
        for (i = 1; i <= _xSize; i++) {
          if (ISNA(predictorPtr[i][_fmRecordIndex[j]])) {
            _fmvSign[i+2][j] = 1;
          }
        }
      }
      _fmvIndex = ivector(1, _xSize+2);
      _fmvSize = 0;
      for (i = 1; i <= _xSize+2; i++) {
        _fmvIndex[i] = 0;
        switch (i) {
        case (uint) (-CENS_IDX):
          for (j = 1; j <= _fmRecordSize; j++) {
            if (_fmvSign[i][j] == 1) {
              _fmvSize ++;
              _fmvIndex[_fmvSize] = CENS_IDX;
              j = _fmRecordSize;
            }
          }
          break;
        case (uint) (-TIME_IDX):
          for (j = 1; j <= _fmRecordSize; j++) {
            if (_fmvSign[i][j] == 1) {
              _fmvSize ++;
              _fmvIndex[_fmvSize] = TIME_IDX;
              j = _fmRecordSize;
            }
          }
          break;
        default:
          for (j = 1; j <= _fmRecordSize; j++) {
            if (_fmvSign[i][j] == 1) {
              _fmvSize ++;
              _fmvIndex[_fmvSize] = i-2;
              j = _fmRecordSize;
            }
          }
          break;
        }
      }  
      _fmvForestSign = imatrix(1, _forestSize,  1, _fmvSize);
      for (j = 1; j <= _forestSize; j++) {
        for (p = 1; p <= _fmvSize; p++) {
          _fmvForestSign[j][p] = 1;
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nMapping of PRED Individuals to Missing Record Index:  ");
        Rprintf("\n    record      index \n");
        for (i = 1; i <= _fobservationSize; i++) {
          if (_fmRecordMap[i] != 0) {
            Rprintf("%10d %10d \n", i, _fmRecordMap[i]);
          }
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nIndex of PRED Individuals with any Missing Outcomes or Predictors:  ");
        Rprintf("\n   element      index \n");
        for (i = 1; i <= _fmRecordSize; i++) {
          Rprintf("%10d %10d \n", i, _fmRecordIndex[i]);
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nIncoming Indices of Missing PRED Outcomes and Predictors:  ");
        Rprintf("\n   element      index \n");
        for (i = 1; i <= _fmvSize; i++) {
          Rprintf("%10d %10d \n", i, _fmvIndex[i]);
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nIncoming Signatures of Missing PRED Outcomes and Predictors:  ");
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
      }
      if (getTraceFlag() & DL2_TRACE) {
        for (i = 1; i <= _fmRecordSize; i++) {
          Rprintf("%10d  ", _fmRecordIndex[i]);
          for (j=1; j <= _xSize+2; j++) {
            Rprintf("%3d", _fmvSign[j][i]);
          }
          Rprintf("\n");
        }
      }
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Missing data analysis of PRED complete -- some found.");  
      }
    }  
  }  
  if (mode == RSF_GROW) {
    if (_mRecordSize > 0) {
      *p_mRecordBootFlag = cmatrix(1, _forestSize, 1, _mRecordSize);
      for (j = 1; j <= _forestSize; j++) {
        for (i = 1; i <= _mRecordSize; i++) {
          (*p_mRecordBootFlag)[j][i] = FALSE;
        }
      }
      *p_mvImputation = dmatrix3(1, _forestSize, 1, _mRecordSize, 1, _mvSize);
      for (i = 1; i <= _forestSize; i++) {
        for (j = 1; j <= _mRecordSize; j++) {
          for (k = 1; k <= _mvSize; k++) {
            (*p_mvImputation)[i][j][k] = NA_REAL;
          }
        }
      }
    }
  }
  else {
    if (_fmRecordSize > 0) {
      *p_mRecordBootFlag = cmatrix(1, _forestSize, 1, _fmRecordSize);
      for (j = 1; j <= _forestSize; j++) {
        for (i = 1; i <= _fmRecordSize; i++) {
          (*p_mRecordBootFlag)[j][i] = ACTIVE;
        }
      }
      *p_mvImputation = dmatrix3(1, _forestSize, 1, _fmRecordSize, 1, _fmvSize);
      for (i = 1; i <= _forestSize; i++) {
        for (j = 1; j <= _fmRecordSize; j++) {
          for (k = 1; k <= _fmvSize; k++) {
            (*p_mvImputation)[i][j][k] = NA_REAL;
          }
        }
      }
    }
  }
  if (*p_mTimeIndexFlag == FALSE) {
      updateTimeIndexArray(NULL);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackAndInitializeMissingArrays() EXIT ...\n");
  }
  return result;
}
void unstackMissingArrays(uint   mode,
                          char  **mRecordBootFlag,
                          double ***mvImputation) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackMissingArrays() ENTRY ...\n");
  }
  free_pcvector(_xType, 1, _xSize);
  free_uivector(_mRecordMap, 1, _observationSize);
  if (_mRecordSize == 0) {
    free_pdvector(_observation, 1, _xSize);
  }
  else {
    free_dvector(_mStatus, 1, _observationSize);
    free_dvector(_mTime, 1, _observationSize);
    free_dmatrix(_observation, 1, _xSize, 1, _observationSize);
    free_uivector(_mRecordIndex, 1, _mRecordSize);
    free_imatrix(_mvSign, 1, _xSize+2, 1, _mRecordSize);
    free_ivector(_mvIndex, 1, _xSize+2);
    free_imatrix(_mvForestSign, 1, _forestSize, 1, _mvSize);
  }  
  if (mode == RSF_PRED) {
    free_uivector(_fmRecordMap, 1, _fobservationSize);
    if (_fmRecordSize == 0) {
      free_pdvector(_fobservation, 1, _xSize);
    }
    else {
      free_dvector(_fmStatus, 1, _fobservationSize);
      free_dvector(_fmTime, 1, _fobservationSize);
      free_dmatrix(_fobservation, 1, _xSize, 1, _fobservationSize);
      free_uivector(_fmRecordIndex, 1, _fmRecordSize);
      free_imatrix(_fmvSign, 1, _xSize+2, 1, _fmRecordSize);
      free_ivector(_fmvIndex, 1, _xSize+2);
      free_imatrix(_fmvForestSign, 1, _forestSize, 1, _fmvSize);
    }  
  }  
  if (mode == RSF_GROW) {
    if (_mRecordSize > 0) {
      free_cmatrix(mRecordBootFlag, 1, _forestSize, 1, _mRecordSize);
      free_dmatrix3(mvImputation, 1, _forestSize, 1, _mRecordSize, 1, _mvSize);
    }
  }
  else {
    if (_fmRecordSize > 0) {
      free_cmatrix(mRecordBootFlag, 1, _forestSize, 1, _fmRecordSize);
      free_dmatrix3(mvImputation, 1, _forestSize, 1, _fmRecordSize, 1, _fmvSize);
    }
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackMissingArrays() EXIT ...\n");
  }
}
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
                               SEXP     *sexpVector) {
  uint sexpIndex;
  uint ensembleSize;
  uint proximitySize;
  uint imputationSize;
  uint  obsSize;
  uint  mRecordSize;
  uint *mRecordIndex;
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  uint i,j,p;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackDefinedOutputObjects() ENTRY ...\n");
  }
  proximitySize = 0;  
  imputationSize = 0;  
  *stackCount = 3;
  if (mode == RSF_GROW) {
    obsSize = _observationSize;
    mRecordSize = _mRecordSize;
    (*stackCount) += 1;
  }
  else {
    obsSize = _fobservationSize;
    mRecordSize = _fmRecordSize;
  }
  ensembleSize = sortedTimeInterestSize * obsSize;
  if (_mup & MUP_PROX) {
    proximitySize = ((obsSize + 1)  * obsSize) / 2; 
    (*stackCount) += 1;
  }
  if (_mup & MUP_VIMP) {
    (*stackCount) += 1;
  }
  if (mRecordSize > 0) {
    imputationSize = (_xSize + 3) * mRecordSize;
    (*stackCount) += 1;
  }
  if (_mup & MUP_TREE) (*stackCount) += 5;
  PROTECT(sexpVector[RSF_OUTP_ID] = allocVector(VECSXP, *stackCount));
  PROTECT(sexpVector[RSF_STRG_ID] = allocVector(STRSXP, *stackCount));
  setAttrib(sexpVector[RSF_OUTP_ID], R_NamesSymbol, sexpVector[RSF_STRG_ID]);
  PROTECT(sexpVector[RSF_FENS_ID] = NEW_NUMERIC(ensembleSize));
  PROTECT(sexpVector[RSF_PERF_ID] = NEW_NUMERIC(_forestSize));
  PROTECT(sexpVector[RSF_LEAF_ID] = NEW_INTEGER(_forestSize));
  *p_fullEnsemble = NUMERIC_POINTER(sexpVector[RSF_FENS_ID]);
  *p_performance  = NUMERIC_POINTER(sexpVector[RSF_PERF_ID]);
  *p_leafCount    = (uint*) INTEGER_POINTER(sexpVector[RSF_LEAF_ID]);
  SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], 0, sexpVector[RSF_FENS_ID]);
  SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], 1, sexpVector[RSF_PERF_ID]);
  SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], 2, sexpVector[RSF_LEAF_ID]);
  SET_STRING_ELT(sexpVector[RSF_STRG_ID], 0, mkChar(sexpString[RSF_FENS_ID]));
  SET_STRING_ELT(sexpVector[RSF_STRG_ID], 1, mkChar(sexpString[RSF_PERF_ID]));
  SET_STRING_ELT(sexpVector[RSF_STRG_ID], 2, mkChar(sexpString[RSF_LEAF_ID]));
  sexpIndex = 3;
  if (mode == RSF_GROW) {
    PROTECT(sexpVector[RSF_OENS_ID] = NEW_NUMERIC(ensembleSize));
    *p_oobEnsemble = NUMERIC_POINTER(sexpVector[RSF_OENS_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_OENS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_OENS_ID]));
    sexpIndex++;
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for GROW mode information complete.");  
    }
  }
  if (_mup & MUP_PROX) {
    PROTECT(sexpVector[RSF_PROX_ID] = NEW_INTEGER(proximitySize));
    *p_proximity = (uint*) INTEGER_POINTER(sexpVector[RSF_PROX_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PROX_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_PROX_ID]));
    sexpIndex++;
    (*p_proximity) --;
    for (i = 1; i <= proximitySize; i++) {
      (*p_proximity)[i] = 0;
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for MUP proximity information complete.");  
    }
  }
  if (_mup & MUP_VIMP) {
    PROTECT(sexpVector[RSF_VIMP_ID] = NEW_NUMERIC(_xSize));
    *p_varImportance = NUMERIC_POINTER(sexpVector[RSF_VIMP_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_VIMP_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_VIMP_ID]));
    sexpIndex++;
    (*p_varImportance) --;
    (*p_vimpEnsembleRun) = dmatrix(1, _xSize, 1, obsSize);
    for (i = 1; i <= _xSize; i++) {
      (*p_varImportance)[i] = NA_REAL;
      for (j = 1; j <= obsSize; j++) {
        (*p_vimpEnsembleRun)[i][j] = 0;
      }
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for MUP variable importance complete.");  
    }
  }
  if (mRecordSize > 0) {
    PROTECT(sexpVector[RSF_MISS_ID] = NEW_NUMERIC(imputationSize));
    *p_imputation = NUMERIC_POINTER(sexpVector[RSF_MISS_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_MISS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_MISS_ID]));
    sexpIndex++;
    *p_sImputePredictorPtr = pdvector(1, _xSize);
    *p_sImputeStatusPtr = (*p_imputation)  + (1 * mRecordSize) - 1;
    *p_sImputeTimePtr   = (*p_imputation)  + (2 * mRecordSize) - 1;
    for (i = 1; i <= _xSize; i++) {
      (*p_sImputePredictorPtr)[i]  = (*p_imputation)  + ((i+2) * mRecordSize) - 1;
    }
    for (i = 1; i <= mRecordSize; i++) {
      if (mode == RSF_GROW) {
        statusPtr = _status;
        timePtr = _time;
        predictorPtr = _observation;
        mRecordIndex = _mRecordIndex;
      }
      else {
        statusPtr = _fstatus;
        timePtr = _ftime;
        predictorPtr = _fobservation;
        mRecordIndex = _fmRecordIndex;
      }
      (*p_imputation)[i-1] = (double) mRecordIndex[i];
      (*p_sImputeStatusPtr)[i] = statusPtr[mRecordIndex[i]];
      (*p_sImputeTimePtr)[i]   = timePtr[mRecordIndex[i]];
      for (j = 1; j <= _xSize; j++) {
        (*p_sImputePredictorPtr)[j][i] = predictorPtr[j][mRecordIndex[i]];
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nImputed Data Output Object:  (at initialization)");
      Rprintf("\n     index   imputation -> \n");
      Rprintf(  "           %12d %12d", CENS_IDX, TIME_IDX);
      for (p=1; p <= _xSize; p++) {
        Rprintf(" %12d", p);
      }
      Rprintf("\n");
      for (i = 1; i <= mRecordSize; i++) {
        Rprintf("%10d", mRecordIndex[i]);
        Rprintf(" %12.4f", (*p_sImputeStatusPtr)[i]);
        Rprintf(" %12.4f", (*p_sImputeTimePtr)[i]);
        for (p = 1; p <= _xSize; p++) {
          Rprintf(" %12.4f", (*p_sImputePredictorPtr)[p][i]);
        }
        Rprintf("\n");
      }
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for imputed data output complete.");  
    }
  }
  if (mode == RSF_GROW) {
    if (_mup & MUP_TREE) {
      *p_root = nodePtrVector(1, _forestSize);
      PROTECT(sexpVector[RSF_SEED_ID] = NEW_INTEGER(1));
      *p_seed = INTEGER_POINTER(sexpVector[RSF_SEED_ID]);
      SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_SEED_ID]);
      SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_SEED_ID]));
      sexpIndex++;
      if (getTraceFlag() & DL1_TRACE) {
        Rprintf("\nAllocating for MUP forest information complete.");  
      }
    }
  }
  *p_ensembleRun = dvector(1, obsSize);
  *p_ensembleDen = uivector(1, obsSize);
  *p_fullEnsemblePtr = pdvector(1, sortedTimeInterestSize);
  if (mode == RSF_GROW) {
    *p_oobEnsemblePtr = pdvector(1, sortedTimeInterestSize);
    for (i = 1; i <= sortedTimeInterestSize; i++) {
      (*p_oobEnsemblePtr)[i]  = (*p_oobEnsemble)  + ((i-1)*(obsSize)) - 1;
      (*p_fullEnsemblePtr)[i] = (*p_fullEnsemble) + ((i-1)*(obsSize)) - 1;
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= sortedTimeInterestSize; j++) {
        (*p_oobEnsemblePtr)[j][i] = 0.0;
        (*p_fullEnsemblePtr)[j][i] = 0.0;
      }
      (*p_ensembleDen)[i] = 0;
    }
  }
  if (mode == RSF_PRED) {
    for (i = 1; i <= sortedTimeInterestSize; i++) {
      (*p_fullEnsemblePtr)[i] = (*p_fullEnsemble) + ((i-1)*(obsSize)) - 1;
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= sortedTimeInterestSize; j++) {
        (*p_fullEnsemblePtr)[j][i] = 0.0;
      }
      (*p_ensembleDen)[i] = 0;
    }
  }
  (*p_performance) --;
  (*p_leafCount)   --;
  for (i = 1; i <= _forestSize; i++) {
    (*p_leafCount)[i] = 0;
    (*p_performance)[i] = NA_REAL;
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Allocation of defined output objects complete:  %10d", sexpIndex);  
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackDefinedOutputObjects() EXIT ...\n");
  }
  return (sexpIndex);
}
uint stackVariableOutputObjects(uint     totalNodeCount,
                                uint   **p_treeID,
                                uint   **p_nodeID,
                                uint   **p_parmID,                                   
                                double **p_spltPT,
                                uint     sexpIndex,
                                char   **sexpString,
                                SEXP    *sexpVector) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackVariableOutputObjects() ENTRY ...\n");
  }
  if (_mup & MUP_TREE) {
    PROTECT(sexpVector[RSF_TREE_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RSF_NODE_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RSF_PARM_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RSF_SPLT_PT] = NEW_NUMERIC(totalNodeCount));    
    *p_treeID = (uint*) INTEGER_POINTER(sexpVector[RSF_TREE_ID]);
    *p_nodeID = (uint*) INTEGER_POINTER(sexpVector[RSF_NODE_ID]);
    *p_parmID = (uint*) INTEGER_POINTER(sexpVector[RSF_PARM_ID]);
    *p_spltPT = NUMERIC_POINTER(sexpVector[RSF_SPLT_PT]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_TREE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_TREE_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_NODE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_NODE_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PARM_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_PARM_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_SPLT_PT]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_SPLT_PT]));
    (*p_treeID) --;
    (*p_nodeID) --;
    (*p_parmID) --;
    (*p_spltPT) --;
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Allocation of variable output objects complete:  %10d", sexpIndex);  
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nstackVariableOutputObjects() EXIT ...\n");
  }
  return (sexpIndex);
}
void unstackDefinedOutputObjects(uint      mode,
                                 Node    **root,
                                 uint      sortedTimeInterestSize,
                                 double  **oobEnsemblePtr,
                                 double  **fullEnsemblePtr,
                                 double   *ensembleRun,
                                 uint     *ensembleDen,
                                 double  **vimpEnsembleRun,
                                 double  **sImputePredictorPtr) {
  uint mRecordSize;
  uint i;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackDefinedOutputObjects() ENTRY ...\n");
  }
  if (_mup & MUP_TREE) {
    for (i = 1; i <= _forestSize; i++) {
      freeTree(root[i]);
    }
    free_nodePtrVector(root, 1, _forestSize);
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nDe-allocating for forest information.");  
    }
  }
  if (mode == RSF_GROW) {
    mRecordSize = _mRecordSize;
  }
  else {
    mRecordSize = _fmRecordSize;
  }
  if (mRecordSize > 0) {
    free_pdvector(sImputePredictorPtr, 1, _xSize);
  }
  free_pdvector(fullEnsemblePtr, 1, sortedTimeInterestSize);
  if (mode == RSF_GROW) {
    free_dvector(ensembleRun, 1, _observationSize);
    free_uivector(ensembleDen, 1, _observationSize);
    free_pdvector(oobEnsemblePtr, 1, sortedTimeInterestSize);
  }
  else {
    free_dvector(ensembleRun, 1, _fobservationSize);
    free_uivector(ensembleDen, 1, _fobservationSize);
  }
  if (_mup & MUP_VIMP) {
    if (mode == RSF_GROW) {
      free_dmatrix(vimpEnsembleRun, 1, _xSize, 1, _observationSize);
    }
    else {
      free_dmatrix(vimpEnsembleRun, 1, _xSize, 1, _fobservationSize);
    }
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nunstackDefinedOutputObjects() EXIT ...\n");
  }
}
void initializeArrays(uint mode, uint *sortedTimeInterestSize) {
  uint i,j;
  uint leadingIndex;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nCommon Incoming Parameters:  ");
    Rprintf("\n         traceFlag:  %10d", getTraceFlag());
    Rprintf("\n               mup:  %10d", _mup);
    Rprintf("\n        forestSize:  %10d", _forestSize);
    Rprintf("\n   observationSize:  %10d", _observationSize);
    Rprintf("\n  timeInterestSize:  %10d", _timeInterestSize);
    Rprintf("\n             xSize:  %10d", _xSize);
    Rprintf("\n");
  }
  if ((getTraceFlag() & RSF_GROW) && (getTraceFlag() & DL1_TRACE)) {
    Rprintf("\nIncoming Random Covariate Weights:  ");
    Rprintf("\n     index       weight");
    for (j=1; j <= _xSize; j++) {
      Rprintf("\n%10d  %12.4f", j, _randomCovariateWeight[j]);
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nIncoming GROW Data:  ");
    Rprintf("\n     index         time     status   predictors -> \n");
  }
  _masterTimeSize = 0;
  for (j = 1; j <= _observationSize; j++) {
    if (!ISNA(_time[j])) {
      _masterTimeSize ++;
      _masterTime[_masterTimeSize] = _time[j];
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    for (j = 1; j <= _observationSize; j++) {
      Rprintf("%10d %12.4f %12.4f", j, _time[j], _status[j]);
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12.4f", (_xData+((i-1)*(_observationSize)))[j-1]);
      }
      Rprintf("\n");
    }
  }
  if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
    Rprintf("\nIncoming PRED Data:  ");
    Rprintf("\n     index         time     status   observations -> \n");
    for (j = 1; j <= _fobservationSize; j++) {
      Rprintf("%10d %12.4f %12.4f", j, _ftime[j], _fstatus[j]);
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12.4f", (_fxData+((i-1)*(_fobservationSize)))[j-1]);
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Initial read complete.");  
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
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nSorted Unique Event Times:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10.4f \n", i, _masterTime[i]);
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
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
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nSorted Distinct Times of Interest:  \n");
    for (i=1; i <= *sortedTimeInterestSize; i++) {
      Rprintf("%10d %10.4f \n", i, _timeInterest[i]);
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of time interest data complete.");  
  }
}
