////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.0
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
#include       "rsfTree.h"
#include  "rsfFactorOps.h"
#include     "rsfImpute.h"
void imputeInteraction (uint treeID, Node *parent) {
  double *valuePtr, *imputePtr;
  uint unsignedIndex;
  uint i,p;
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nimputeInteraction() ENTRY ...\n");
  }
  if (!(_fmRecordSize > 0)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to impute interaction with no missingness.  ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  for (p = 1; p <= _fmvSize; p++) {
    if (_fmvForestSign[treeID][p] != -1) {
      switch (_fmvIndex[p]) {
      case CENS_IDX:
        unsignedIndex = abs(_fmvIndex[p]);
        valuePtr = _status;
        imputePtr = _fstatus;
        break;
      case TIME_IDX:
        unsignedIndex = abs(_fmvIndex[p]);
        valuePtr = _time;
        imputePtr = _ftime;
        break;
      default:
        unsignedIndex = (uint) _fmvIndex[p] + 2;
        valuePtr = _observation[(uint) _fmvIndex[p]];
        imputePtr = _fobservation[(uint) _fmvIndex[p]];
        break;
      }
      for (i = 1; i <= _fobservationSize; i++) {
        if (_fnodeMembership[i] == parent) {
          if (_fmRecordMap[i] > 0) {
            if(_fmvSign[unsignedIndex][_fmRecordMap[i]] == 1) {
              imputePtr[i] = valuePtr[_intrIndividual[i]];
              if (getTraceFlag() & MISS_HGH_TRACE) {
                Rprintf("\nImputed Interaction Value for:  ");
                Rprintf("\n[indv, coord] = [%10d, %10d] \n", i, _fmvIndex[p]);
                Rprintf("%12.4f \n", imputePtr[i]);
              }
            }  
          }  
        }  
      }  
    }  
  }  
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nimputeInteraction() EXIT ...\n");
  }
}
char imputeNode (uint     type,
                 char     seedChainFlag,
                 uint     treeID, 
                 Node    *parent) {
  uint     obsSize;
  double  *status;
  double  *time;
  double **predictor;
  Node   **nodeMembership;
  uint    *mRecordMap;
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  int    **mvForestSign;
  double *valuePtr, *imputePtr;
  char mPredictorFlag;
  uint unsignedIndex;
  char result;
  uint i,p;
  uint localDistributionSize;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeNode() ENTRY ...\n");
  }
  result = FALSE;
  switch (type) {
  case RSF_PRED:
    if (_fmRecordSize > 0) {
      obsSize = _fobservationSize;
      status = _fstatus;
      time = _ftime;
      predictor = _fobservation;
      nodeMembership = _fnodeMembership;
      mRecordMap = _fmRecordMap;
      mRecordSize = _fmRecordSize;
      mRecordIndex = _fmRecordIndex;
      mvSize = _fmvSize;
      mvSign = _fmvSign;
      mvIndex = _fmvIndex;
      mvForestSign = _fmvForestSign;
      result = TRUE;
    }
    break;
  default:
    if (_mRecordSize > 0) {
      obsSize = _observationSize;
      status = _status;
      time = _time;
      predictor = _observation;
      nodeMembership = _nodeMembership;
      mRecordMap = _mRecordMap;
      mRecordSize = _mRecordSize;
      mRecordIndex = _mRecordIndex;
      mvSize = _mvSize;
      mvSign = _mvSign;
      mvIndex = _mvIndex;
      mvForestSign = _mvForestSign;
      result = TRUE;
    }
    break;
  }
  if (result == FALSE) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to impute node with no missingness in type:  %10d", type);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
   if (getTraceFlag() & MISS_LOW_TRACE) {
    if (type == RSF_GROW) {
      Rprintf("\nData for type GROW (before imputation) in leaf:  ");
    }
    else {
      Rprintf("\nData for type NON-GROW (before imputation) in leaf:  ");
    }
    Rprintf("%10d", parent -> leafCount);
    Rprintf("\n     index   imputation -> \n");
    Rprintf(  "          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      if (nodeMembership[mRecordIndex[i]] == parent) {
        Rprintf("%10d", mRecordIndex[i]);
        for (p = 1; p <= mvSize; p++) {
          switch (mvIndex[p]) {
          case CENS_IDX:
            valuePtr = status;
            break;
          case TIME_IDX:
            valuePtr = time;
            break;
          default:
            valuePtr = predictor[(uint) mvIndex[p]];
            break;
          }
          Rprintf(" %12.4f", valuePtr[mRecordIndex[i]]);
        }
        Rprintf("\n");
      }
    }
  }
  double *localDistribution = dvector(1, _observationSize);
  for (p = 1; p <= mvSize; p++) {
    if (mvForestSign[treeID][p] != -1) {
      switch (mvIndex[p]) {
      case CENS_IDX:
        unsignedIndex = abs(mvIndex[p]);
        valuePtr = _status;
        imputePtr = status;
        break;
      case TIME_IDX:
        unsignedIndex = abs(mvIndex[p]);
        valuePtr = _time;
        imputePtr = time;
        break;
      default:
        unsignedIndex = (uint) mvIndex[p] + 2;
        valuePtr = _observation[(uint) mvIndex[p]];
        imputePtr = predictor[(uint) mvIndex[p]];
        break;
      }
      localDistributionSize = 0;
      for (i = 1; i <= _observationSize; i++) {
        if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
          mPredictorFlag = TRUE;
          if (_mRecordMap[_bootMembershipIndex[i]] == 0) {
            mPredictorFlag = FALSE;
          }
          else if (_mvSign[unsignedIndex][_mRecordMap[_bootMembershipIndex[i]]] == 0) {
            mPredictorFlag = FALSE;
          }
          if (mPredictorFlag == FALSE) {
            localDistributionSize ++;
            localDistribution[localDistributionSize] = valuePtr[_bootMembershipIndex[i]];
          }
        }
      }  
      if (localDistributionSize > 0) {
        for (i = 1; i <= obsSize; i++) {
          if (nodeMembership[i] == parent) {
            if (mRecordMap[i] > 0) {
              if(mvSign[unsignedIndex][mRecordMap[i]] == 1) {
                imputePtr[i] = getSampleValue(localDistribution, localDistributionSize, seedChainFlag);
                if (getTraceFlag() & MISS_HGH_TRACE) {
                  Rprintf("\nNode Imputed Value for:  ");
                  Rprintf("\n[indv, coord] = [%10d, %10d] \n", i, mvIndex[p]);
                  Rprintf("%12.4f \n", imputePtr[i]);
                }
              }  
            }  
          }  
        }  
      }  
      else {
        if (getTraceFlag() & MISS_MED_TRACE) {
          Rprintf("\nNode Value Not Imputable for:  ");
          Rprintf("\n[coord] = [%10d] \n", mvIndex[p]);
          Rprintf("\nPrevious value retained. \n");
        }
      }
    }  
  }  
  free_dvector(localDistribution, 1, _observationSize);
  if (getTraceFlag() & MISS_LOW_TRACE) {
    if (type == RSF_GROW) {
      Rprintf("\nData For GROW (after imputation) in leaf:  ");
    }
    else {
      Rprintf("\nData For PRED (after imputation) in leaf:  ");
    }
    Rprintf("%10d", parent -> leafCount);
    Rprintf("\n     index   imputation -> \n");
    Rprintf(  "          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      if (nodeMembership[mRecordIndex[i]] == parent) {
        Rprintf("%10d", mRecordIndex[i]);
        for (p = 1; p <= mvSize; p++) {
          switch (mvIndex[p]) {
          case CENS_IDX:
            valuePtr = status;
            break;
          case TIME_IDX:
            valuePtr = time;
            break;
          default:
            valuePtr = predictor[(uint) mvIndex[p]];
            break;
          }
          Rprintf(" %12.4f", valuePtr[mRecordIndex[i]]);
        }
        Rprintf("\n");
      }
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeNode() EXIT ...\n");
  }
  return TRUE;
}
void imputeTree(uint mode, uint b, Node *parent, char rootFlag) {
  char daughterFlag;
  char result;
  uint i;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeTree() ENTRY ...\n");
  }
  switch (mode) {
  case RSF_GROW:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unused GROW case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  case RSF_PRED:
    if (_mRecordSize > 0) {
      if (rootFlag == TRUE) {
        result = TRUE;
      }
      else {
        result = testNodeSize(parent);
      }
      if (result) {
        imputeNode(RSF_GROW, 
                   TRUE, 
                   b, 
                   parent);
        if (_mTimeIndexFlag == TRUE) {
          updateTimeIndexArray(parent);
        }
      }
    }
    if (_fmRecordSize > 0) {
      imputeNode(RSF_PRED, 
                 FALSE, 
                 b, 
                 parent);
    }
    break;
  case RSF_INTR:
    if (_mRecordSize > 0) {
      if (rootFlag == TRUE) {
        result = TRUE;
      }
      else {
        result = testNodeSize(parent);
      }
      if (result) {
        imputeNode(RSF_GROW, 
                   TRUE, 
                   b, 
                   parent);
        if (_mTimeIndexFlag == TRUE) {
          updateTimeIndexArray(parent);
        }
      }
    }
    if (_fmRecordSize > 0) {
      imputeInteraction(b, parent);
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
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nForking On:  ");
      Rprintf("\n   LeafCnt   SpltParm  ");
      Rprintf("\n%10d %10d ", parent -> leafCount, parent -> splitParameter);
      if (parent -> splitValueFactSize > 0) {
        for (i = parent -> splitValueFactSize; i >= 1; i--) {
          Rprintf("%8x ", (parent -> splitValueFactPtr)[i]);
        }
        Rprintf("\n");
      }
      else {
        Rprintf(" %20.4f \n", parent -> splitValueCont);
      }
    }
    for (i=1; i <= _observationSize; i++) {
      if (_nodeMembership[i] == parent) {
        daughterFlag = RIGHT;
        if (strcmp(_xType[parent -> splitParameter], "C") == 0) {
          if (getTraceFlag() & MISS_LOW_TRACE) {
            Rprintf("\nNode Membership:  %10d %16d", i, (uint) _observation[parent -> splitParameter][i]);
          }
          daughterFlag = splitOnFactor((uint) _observation[parent -> splitParameter][i], parent -> splitValueFactPtr);
        }
        else {
          if (getTraceFlag() & MISS_LOW_TRACE) {
            Rprintf("\nNode Membership:  %10d %16.4f, %16.4f", i, _observation[parent -> splitParameter][i], parent -> splitValueCont);
          }
          if ( _observation[parent -> splitParameter][i] <= (parent -> splitValueCont) ) {
            daughterFlag = LEFT;
          }
        }
        if (daughterFlag == LEFT) {
          _nodeMembership[i] = parent -> left;
          if (getTraceFlag() & MISS_LOW_TRACE) {
            Rprintf(" --> LEFT ");
          }
        }
        else {
          _nodeMembership[i] = parent -> right;
          if (getTraceFlag() & MISS_LOW_TRACE) {
            Rprintf(" --> RGHT ");
          }
        }
      }
    }
    if (mode != RSF_GROW) {
      for (i=1; i <= _fobservationSize; i++) {
        if (_fnodeMembership[i] == parent) {
          daughterFlag = RIGHT;
          if (strcmp(_xType[parent -> splitParameter], "C") == 0) {
            daughterFlag = splitOnFactor((uint) _fobservation[parent -> splitParameter][i], parent -> splitValueFactPtr);
          }
          else {
            if ( _fobservation[parent -> splitParameter][i] <= (parent -> splitValueCont) ) {
              daughterFlag = LEFT;
            }
          }
          if (daughterFlag == LEFT) {
            _fnodeMembership[i] = parent -> left;
          }
          else {
            _fnodeMembership[i] = parent -> right;
          }
        }
      }
    }
    imputeTree(mode, b, parent -> left, FALSE);
    imputeTree(mode, b, parent -> right, FALSE);
  }  
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeTree() EXIT ...\n");
  }
}
void imputeUpdateSummary (uint      mode, 
                          double  *statusPtr, 
                          double  *timePtr, 
                          double **predictorPtr, 
                          uint treeID, 
                          double ***dmvImputationPtr) {
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  int    **mvForestSign;
  double  *valuePtr;
  uint     unsignedIndex;
  uint     obsSize;
  char result;
  uint i, p;
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeUpdateSummary(%2d) ENTRY ...\n", mode);
  }
  result = FALSE;
  if (mode == RSF_GROW) {
    if (_mRecordSize > 0) {
      mRecordSize = _mRecordSize;
      mRecordIndex = _mRecordIndex;
      mvSize = _mvSize;
      mvSign = _mvSign;
      mvIndex = _mvIndex;
      mvForestSign = _mvForestSign;
      obsSize = _observationSize;
      result = TRUE;
    }
  }
  else {
    if (_fmRecordSize > 0) {
      mRecordSize = _fmRecordSize;
      mRecordIndex = _fmRecordIndex;
      mvSize = _fmvSize;
      mvSign = _fmvSign;
      mvIndex = _fmvIndex;
      mvForestSign = _fmvForestSign;
      obsSize = _fobservationSize;
      result = TRUE;
    }
  }
  if (result == FALSE) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to update forest impute data with no missingness in mode:  %10d", mode);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nShadowed data for tree:  ");
    Rprintf("\n       index       status         time observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= _xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (i = 1; i <= obsSize; i++) {
      Rprintf("%12d %12.4f %12.4f", i, statusPtr[i], timePtr[i]);
      for (p = 1; p <= _xSize; p++) {
        Rprintf(" %12.4f", predictorPtr[p][i]);
      }
      Rprintf("\n");
    }
    Rprintf("\nForest Impute Data Structure For Tree (before update):  ");
    Rprintf("\n       index   imputation -> \n");
    Rprintf(  "          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      Rprintf("%12d", mRecordIndex[i]);
      for (p = 1; p <= mvSize; p++) {
        Rprintf(" %12.4f", dmvImputationPtr[treeID][i][p]);
      }
      Rprintf("\n");
    }
  }
  for (p = 1; p <= mvSize; p++) {
    if (mvForestSign[treeID][p] != -1) {
      for (i = 1; i <= mRecordSize; i++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          unsignedIndex = abs(mvIndex[p]);
          valuePtr = statusPtr;
          break;
        case TIME_IDX:
          unsignedIndex = abs(mvIndex[p]);
          valuePtr = timePtr;
          break;
        default:
          unsignedIndex = (uint) mvIndex[p] + 2;
          valuePtr = predictorPtr[(uint) mvIndex[p]];
          break;
        }
        if (mvSign[unsignedIndex][i] == 1) {
          if (ISNA(valuePtr[mRecordIndex[i]])) {
            Rprintf("\nRSF:  *** ERROR *** ");
            Rprintf("\nRSF:  Attempt to update forest impute data with invalid shadowed value, NA. ");
            Rprintf("\nRSF:  Invalid value for:  [indv][outcome/predictor] = [%10d][%10d] ", mRecordIndex[i], mvIndex[p]);
            Rprintf("\nRSF:  Please Contact Technical Support.");
            Rprintf("\nRSF:  The application will now exit.\n");
            Rprintf("\nDiagnostic Trace of Shadowed Data:  ");
            Rprintf("\n       index   imputation -> \n");
            Rprintf(  "            ");
            for (p=1; p <= mvSize; p++) {
              Rprintf(" %12d", mvIndex[p]);
            }
            Rprintf("\n");
            for (i = 1; i <= mRecordSize; i++) {
              Rprintf("%12d", mRecordIndex[i]);
              for (p = 1; p <= mvSize; p++) {
                switch (mvIndex[p]) {
                case CENS_IDX:
                  valuePtr = statusPtr;
                  break;
                case TIME_IDX:
                  valuePtr = timePtr;
                  break;
                default:
                  valuePtr = predictorPtr[(uint) mvIndex[p]];
                  break;
                }
                Rprintf(" %12.4f", valuePtr[mRecordIndex[i]]);
              }
              Rprintf("\n");
            }
            exit(TRUE);
          }  
          else {
            dmvImputationPtr[treeID][i][p] = valuePtr[mRecordIndex[i]];
          }
        }  
      }  
    }  
  }  
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nForest Impute Data Structure For Tree (after update):  ");
    Rprintf("\n       index   imputation -> \n");
    Rprintf(  "          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      Rprintf("%12d", mRecordIndex[i]);
      for (p = 1; p <= mvSize; p++) {
        Rprintf(" %12.4f", dmvImputationPtr[treeID][i][p]);
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeUpdateSummary(%2d) EXIT ...\n", mode);
  }
}
void imputeUpdateShadow (uint      mode, 
                         char      selectionFlag,
                         double ***dmvImputation,
                         double   *shadowStatus, 
                         double   *shadowTime, 
                         double  **shadowPredictor) {
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  double  *outTime;
  double  *outStatus;
  double **outPredictor;
  double  *valuePtr;
  double  *outputPtr;
  uint unsignedIndex;
  char outcomeFlag;
  uint i, p;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeUpdateShadow(%2d) ENTRY ...\n", mode);
  }
  mRecordSize  = 0;     
  mRecordIndex = NULL;  
  mvSize       = 0;     
  mvSign       = NULL;  
  mvIndex      = NULL;  
  outStatus    = NULL;  
  outTime      = NULL;  
  outPredictor = NULL;  
  valuePtr     = NULL;  
  outputPtr    = NULL;  
  unsignedIndex = 0;    
  switch (mode) {
  case RSF_PRED:
    mRecordSize = _fmRecordSize;
    mRecordIndex = _fmRecordIndex;
    mvSize = _fmvSize;
    mvSign = _fmvSign;
    mvIndex = _fmvIndex;
    outStatus    = _sImputeStatusPtr;
    outTime      = _sImputeTimePtr;
    outPredictor = _sImputePredictorPtr;
    break;
  case RSF_INTR:
    mRecordSize = _fmRecordSize;
    mRecordIndex = _fmRecordIndex;
    mvSize = _fmvSize;
    mvSign = _fmvSign;
    mvIndex = _fmvIndex;
    outStatus    = _sOOBImputeStatusPtr;
    outTime      = _sOOBImputeTimePtr;
    outPredictor = _sOOBImputePredictorPtr;
    break;
  default:
    mRecordSize = _mRecordSize;
    mRecordIndex = _mRecordIndex;
    mvSize = _mvSize;
    mvSign = _mvSign;
    mvIndex = _mvIndex;
    if ((selectionFlag == TRUE) || (selectionFlag == ACTIVE)) {
      outStatus    = _sImputeStatusPtr;
      outTime      = _sImputeTimePtr;
      outPredictor = _sImputePredictorPtr;
    }
    else {
      outStatus    = _sOOBImputeStatusPtr;
      outTime      = _sOOBImputeTimePtr;
      outPredictor = _sOOBImputePredictorPtr;
    }
    break;
  }
  if (mRecordSize == 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to update shadow data with no missingness in mode:  %10d", mode);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nImputed Shadow Data:  (before transfer)");
    Rprintf("\n       index   imputation -> \n");
    Rprintf(  "          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      Rprintf("%12d", mRecordIndex[i]);
      for (p = 1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          Rprintf(" %12.4f", shadowStatus[mRecordIndex[i]]);
          break;
        case TIME_IDX:
          Rprintf(" %12.4f", shadowTime[mRecordIndex[i]]);
          break;
        default:
          if (shadowPredictor != NULL) {
            Rprintf(" %12.4f", shadowPredictor[(uint) mvIndex[p]][mRecordIndex[i]]);
          }
          break;
        }
      }
      Rprintf("\n");
    }
  }
  outcomeFlag = TRUE;
  for (p = 1; p <= mvSize; p++) {
    for (i = 1; i <= mRecordSize; i++) {
      switch (mvIndex[p]) {
      case CENS_IDX:
        unsignedIndex = abs(mvIndex[p]);
        valuePtr = shadowStatus;
        outputPtr = outStatus;
        break;
      case TIME_IDX:
        unsignedIndex = abs(mvIndex[p]);
        valuePtr = shadowTime;
        outputPtr = outTime;
        break;
      default:
        if (shadowPredictor != NULL) {
          unsignedIndex = (uint) mvIndex[p] + 2;
          valuePtr = shadowPredictor[(uint) mvIndex[p]];
          outputPtr = outPredictor[(uint) mvIndex[p]];
        }
        outcomeFlag = FALSE;
        break;
      }
      if (outcomeFlag || (shadowPredictor != NULL)) {
        if (mvSign[unsignedIndex][i] == 1) {
          if (ISNA(outputPtr[i])) {
            if (getTraceFlag() & MISS_MED_TRACE) {
              Rprintf("\nShadow Data Receiving NA from SEXP Output:  ");
              Rprintf("\n[indv, OUT/PRED] = [%10d, %10d] \n", mRecordIndex[i], mvIndex[p]);
            }
          }
          valuePtr[mRecordIndex[i]] = outputPtr[i];
        }  
      }
    }  
  }  
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nImputed Shadow Data:  (after transfer)");
    Rprintf("\n       index   imputation -> \n");
    Rprintf(  "          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      Rprintf("%12d", mRecordIndex[i]);
      for (p = 1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          Rprintf(" %12.4f", shadowStatus[mRecordIndex[i]]);
          break;
        case TIME_IDX:
          Rprintf(" %12.4f", shadowTime[mRecordIndex[i]]);
          break;
        default:
          if (shadowPredictor != NULL) {
            Rprintf(" %12.4f", shadowPredictor[(uint) mvIndex[p]][mRecordIndex[i]]);
          }
          break;
        }
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeUpdateShadow(%2d) EXIT ...\n", mode);
  }
}
void imputeSummary(uint      mode,
                   char      selectionFlag,
                   char    **dmRecordBootFlag,
                   double ***dmvImputation) {
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeSummary(%2d) ENTRY ...\n", selectionFlag);
  }
  imputeCommon(mode,
               _forestSize,
               selectionFlag, 
               dmRecordBootFlag,
               dmvImputation,
               TRUE);
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeSummary(%2d) EXIT ...\n", selectionFlag);
  }
}
void imputeConcordance(uint      mode,
                       uint      b,
                       char    **dmRecordBootFlag,
                       double ***dmvImputation,
                       double   *tempStatus,
                       double   *tempTime) {
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeConcordance() ENTRY ...\n");
  }
  switch(mode) {
  case RSF_GROW:
    imputeCommon(mode,
                 b,
                 FALSE, 
                 dmRecordBootFlag,
                 dmvImputation,
                 FALSE);
    imputeUpdateShadow(mode,
                       FALSE, 
                       dmvImputation, 
                       tempStatus, 
                       tempTime, 
                       NULL);
    break;
  case RSF_PRED:
    imputeCommon(mode,
                 b,
                 ACTIVE, 
                 dmRecordBootFlag,
                 dmvImputation,
                 FALSE);
    imputeUpdateShadow(mode,
                       ACTIVE, 
                       dmvImputation, 
                       tempStatus, 
                       tempTime, 
                       NULL);
    break;
  case RSF_INTR:
    imputeCommon(mode,
                 b,
                 FALSE, 
                 dmRecordBootFlag,
                 dmvImputation,
                 FALSE);
    imputeUpdateShadow(mode,
                       FALSE, 
                       dmvImputation, 
                       tempStatus, 
                       tempTime, 
                       NULL);
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeConcordance() EXIT ...\n");
  }
}
void imputeCommon(uint      mode,
                  uint      b,
                  char      selectionFlag,
                  char    **dmRecordBootFlag,
                  double ***dmvImputation,
                  char      predictorFlag) {
  uint i,p,tree;
  uint q,s;
  char mPredictorFlag;
  char outcomeFlag;
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  int    **mvForestSign;
  double  *outTime;
  double  *outStatus;
  double **outPredictor;
  double *valuePtr;
  double *naivePtr;
  uint    unsignedIndex;
  double imputedValue;
  uint localDistributionSize;
  uint maxDistributionSize;
  char result;
  if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nimputeCommon(%2d) ENTRY ...\n", selectionFlag);
  }
  valuePtr      = NULL;  
  naivePtr      = NULL;  
  unsignedIndex = 0;     
  if ((b < 1) || (b > _forestSize)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid tree count in imputeCommon():  %10d", b);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if ((selectionFlag != TRUE) && (selectionFlag != FALSE) && (selectionFlag != ACTIVE)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid selectionFlag in imputeCommon():  %10d", selectionFlag);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  result = FALSE;
  switch (mode) {
  case RSF_PRED:
    if (_fmRecordSize > 0) {
      mRecordSize = _fmRecordSize;
      mRecordIndex = _fmRecordIndex;
      mvSize = _fmvSize;
      mvSign = _fmvSign;
      mvIndex = _fmvIndex;
      mvForestSign = _fmvForestSign;
      maxDistributionSize = ((_observationSize) > (_forestSize)) ? (_observationSize) : (_forestSize);
      outStatus    = _sImputeStatusPtr;
      outTime      = _sImputeTimePtr;
      outPredictor = _sImputePredictorPtr;
      result = TRUE;
    }
    break;
  case RSF_INTR:
    if (_fmRecordSize > 0) {
      mRecordSize = _fmRecordSize;
      mRecordIndex = _fmRecordIndex;
      mvSize = _fmvSize;
      mvSign = _fmvSign;
      mvIndex = _fmvIndex;
      mvForestSign = _fmvForestSign;
      maxDistributionSize = ((_observationSize) > (_forestSize)) ? (_observationSize) : (_forestSize);
      outStatus    = _sOOBImputeStatusPtr;
      outTime      = _sOOBImputeTimePtr;
      outPredictor = _sOOBImputePredictorPtr;
      result = TRUE;
    }
    break;
  default:
    if (_mRecordSize > 0) {
      mRecordSize = _mRecordSize;
      mRecordIndex = _mRecordIndex;
      mvSize = _mvSize;
      mvSign = _mvSign;
      mvIndex = _mvIndex;
      mvForestSign = _mvForestSign;
      maxDistributionSize = ((_observationSize) > (_forestSize)) ? (_observationSize) : (_forestSize);
      if ((selectionFlag == TRUE) || (selectionFlag == ACTIVE)) {
        outStatus    = _sImputeStatusPtr;
        outTime      = _sImputeTimePtr;
        outPredictor = _sImputePredictorPtr;
      }
      else {
        outStatus    = _sOOBImputeStatusPtr;
        outTime      = _sOOBImputeTimePtr;
        outPredictor = _sOOBImputePredictorPtr;
      }
      result = TRUE;
    }
    break;
  }
  if (result == FALSE) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to impute in imputeCommon() with no missingness in mode:  %10d", mode);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nOutputs (SEXP) (before update):  ");
    Rprintf("\n        index   imputation -> \n");
    Rprintf(  "             ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      Rprintf("%12d", mRecordIndex[i]);
      for (p = 1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          Rprintf(" %12.4f", outStatus[i]);
          break;
        case TIME_IDX:
          Rprintf(" %12.4f", outTime[i]);
          break;
        default:
          if (predictorFlag == TRUE) {
            Rprintf(" %12.4f", outPredictor[(uint) mvIndex[p]][i]);
          }
          break;
        }
      }
      Rprintf("\n");
    }
  }
  imputedValue = 0.0;  
  double *localDistribution = dvector(1, maxDistributionSize);
  char  *naiveMvFlag = cvector(1, mvSize);
  char **naiveSign = cmatrix(1, mRecordSize, 1, mvSize);
  for (p = 1; p <= mvSize; p++) {
    naiveMvFlag[p] = FALSE;
  }
  for (i = 1; i <= mRecordSize; i++) {
    outcomeFlag = TRUE;
    for (p = 1; p <= mvSize; p++) {
      naiveSign[i][p] = FALSE;  
      switch (mvIndex[p]) {
      case CENS_IDX:
        unsignedIndex = abs(mvIndex[p]);
        break;
      case TIME_IDX:
        unsignedIndex = abs(mvIndex[p]);
        break;
      default:
        if (predictorFlag == TRUE) {
          unsignedIndex = (uint) mvIndex[p] + 2;
        }
        outcomeFlag = FALSE;
        break;
      }
      if (outcomeFlag || (predictorFlag == TRUE)) {
        if (mvSign[unsignedIndex][i] == 1) {
          localDistributionSize = 0;
          for (tree = 1; tree <= b; tree++) {
            if (_leafCount_[tree] > 0) {
              if (mvForestSign[tree][p] == 1) {
                if ((dmRecordBootFlag[tree][i] == selectionFlag) || (selectionFlag == ACTIVE)) {
                  if (!ISNA(dmvImputation[tree][i][p])) {
                    localDistribution[++localDistributionSize] = dmvImputation[tree][i][p];
                  }
                  else {
                    Rprintf("\nRSF:  *** ERROR *** ");
                    Rprintf("\nRSF:  Invalid imputed value, NA:  [tree][i][p] = [%10d][%10d][%10d] ", tree, i, p);
                    Rprintf("\nRSF:  Please Contact Technical Support.");
                    Rprintf("\nRSF:  The application will now exit.\n");
                    Rprintf("\nDiagnostic Trace of (tree, outcome/predictor) pairs for [indv, outcome/predictor] = [%10d][%10d] ", mRecordIndex[i], mvIndex[p]);
                    Rprintf("\n      tree   imputation -> \n");
                    Rprintf(  "          ");
                    for (q=1; q <= mvSize; q++) {
                      Rprintf(" %12d", mvIndex[q]);
                    }
                    Rprintf("\n");
                    for (s = 1; s <= _forestSize; s++) {
                      Rprintf("%10d", s);
                      for (q = 1; q <= mvSize; q++) {
                        Rprintf(" %12.4f", dmvImputation[s][i][q]);
                      }
                      Rprintf("\n");
                    }
                    exit(TRUE);
                  }  
                }  
              }  
            }  
          }  
          if (localDistributionSize > 0) {
            switch (mvIndex[p]) {
            case CENS_IDX:
              imputedValue = getMaximalValue(localDistribution, localDistributionSize);
              outStatus[i] = imputedValue;
              if (getTraceFlag() & MISS_HGH_TRACE) {
                Rprintf("\nSummary Status Imputed Value for:  ");
                Rprintf("\n[indv, CENS] = [%10d, %10d] \n", mRecordIndex[i], mvIndex[p]);
                Rprintf("%10d \n", (int) outStatus[i]);
              }
              break;
            case TIME_IDX:
              imputedValue = getMeanValue(localDistribution, localDistributionSize);
              outTime[i] = imputedValue;
              if (getTraceFlag() & MISS_HGH_TRACE) {
                Rprintf("\nSummary Time Imputed Value for:  ");
                Rprintf("\n[indv, TIME] = [%10d, %10d] \n", mRecordIndex[i], mvIndex[p]);
                Rprintf("%12.4f \n", outTime[i]);
              }
              break;
            default:
              if (strcmp(_xType[(uint) mvIndex[p]], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize);
              }
              outPredictor[(uint) mvIndex[p]][i] = imputedValue;
              if (getTraceFlag() & MISS_HGH_TRACE) {
                Rprintf("\nSummary Predictor Imputed Value for:  ");
                Rprintf("\n[indv, pred] = [%10d, %10d] \n", mRecordIndex[i], mvIndex[p]);
                Rprintf("%12.4f \n", outPredictor[(uint) mvIndex[p]][i]);
              }
              break;
            }
          }  
          else {
            naiveMvFlag[p] = TRUE;
            naiveSign[i][p] = TRUE;
            if (getTraceFlag() & MISS_HGH_TRACE) {
              Rprintf("\n[Indv, Pred] requiring naive imputation:  (%10d, %10d) \n", mRecordIndex[i], p);
              Rprintf("\nDiagnostic Trace of (tree, outcome/predictor) pairs for [indv, outcome/predictor] = [%10d][%10d] ", mRecordIndex[i], mvIndex[p]);
              Rprintf("\n      tree   imputation -> \n");
              Rprintf(  "          ");
              for (q=1; q <= mvSize; q++) {
                Rprintf(" %12d", mvIndex[q]);
              }
              Rprintf("\n");
              for (s = 1; s <= _forestSize; s++) {
                Rprintf("%10d", s);
                for (q = 1; q <= mvSize; q++) {
                  Rprintf(" %12.4f", dmvImputation[s][i][q]);
                }
                Rprintf("\n");
              }
            }          
          }
        }  
      }  
      else {
        p = mvSize;
      }
    }  
  }  
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nOutputs (SEXP) (before naive imputation):  ");
    Rprintf("\n       index   imputation -> \n");
    Rprintf(  "            ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      Rprintf("%12d", mRecordIndex[i]);
      for (p = 1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          Rprintf(" %12.4f", outStatus[i]);
          break;
        case TIME_IDX:
          Rprintf(" %12.4f", outTime[i]);
          break;
        default:
          if (predictorFlag == TRUE) {
            Rprintf(" %12.4f", outPredictor[(uint) mvIndex[p]][i]);
          }
          break;
        }
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nNaive Imputation Flag of Outcomes and Predictors:  ");
    Rprintf("\n     index flags -> \n");
    Rprintf(  "          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf("%4d", mvIndex[p]);
    }
    Rprintf("\n          ");
    for (p=1; p <= mvSize; p++) {
      Rprintf("%4d", naiveMvFlag[p]);
    }
    for (i=1; i <= mRecordSize; i++) {
      Rprintf("\n");
      Rprintf("%10d", mRecordIndex[i]);
      for (p=1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          Rprintf("%4d", naiveSign[i][p]);
          break;
        case TIME_IDX:
          Rprintf("%4d", naiveSign[i][p]);
          break;
        default:
          if (predictorFlag == TRUE) {
            Rprintf("%4d", naiveSign[i][p]);
          }
          break;
        }
      }
    }
    Rprintf("\n");
  }
  outcomeFlag = TRUE;
  for (p = 1; p <= mvSize; p++) {
    switch (mvIndex[p]) {
    case CENS_IDX:
      unsignedIndex = abs(mvIndex[p]);
      valuePtr = _status;
      naivePtr = outStatus;
      break;
    case TIME_IDX:
      unsignedIndex = abs(mvIndex[p]);
      valuePtr = _time;
      naivePtr = outTime;
      break;
    default:
      if (predictorFlag == TRUE) {
        unsignedIndex = (uint) mvIndex[p] + 2;
        valuePtr = _observation[(uint) mvIndex[p]];
        naivePtr = outPredictor[(uint) mvIndex[p]];
      }
      outcomeFlag = FALSE;
      break;
    }
    if (outcomeFlag || (predictorFlag == TRUE)) {
      if (naiveMvFlag[p] == TRUE) {
        localDistributionSize = 0;
        for (i=1; i <= _observationSize; i++) {
          mPredictorFlag = TRUE;
          if (_mRecordMap[i] == 0) {
            mPredictorFlag = FALSE;
          }
          else if (_mvSign[unsignedIndex][_mRecordMap[i]] == 0) {
            mPredictorFlag = FALSE;
          }
          if (mPredictorFlag == FALSE) {
            localDistribution[++localDistributionSize] = valuePtr[i];
          }
        }  
        if (localDistributionSize > 0) {
          for (i=1; i <= mRecordSize; i++) {
            if (naiveSign[i][p] == TRUE) {
              naivePtr[i] = getSampleValue(localDistribution, localDistributionSize, FALSE);
              if (getTraceFlag() & MISS_HGH_TRACE) {
                Rprintf("\nSummary Imputed Value for:  ");
                Rprintf("\n[indv, outcome/pred] = [%10d, %10d] \n", mRecordIndex[i], mvIndex[p]);
                Rprintf("%12.4f \n", naivePtr[i]);
              }
            }
          }
        }  
        else {
          if (getTraceFlag() & MISS_HGH_TRACE) {
            Rprintf("\nNaive imputation failed for predictor:  %10d", mvIndex[p]);
          }
        }
      }  
    }  
    else {
      p = mvSize;
    }
  }  
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nOutputs (SEXP) (after naive imputation):  ");
    Rprintf("\n        index   imputation -> \n");
    Rprintf(  "             ");
    for (p=1; p <= mvSize; p++) {
      Rprintf(" %12d", mvIndex[p]);
    }
    Rprintf("\n");
    for (i = 1; i <= mRecordSize; i++) {
      Rprintf(" %12d", mRecordIndex[i]);
      for (p = 1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          Rprintf(" %12.4f", outStatus[i]);
          break;
        case TIME_IDX:
          Rprintf(" %12.4f", outTime[i]);
          break;
        default:
          if (predictorFlag == TRUE) {
            Rprintf(" %12.4f", outPredictor[(uint) mvIndex[p]][i]);
          }
          break;
        }
      }
      Rprintf("\n");
    }
  }
  free_dvector(localDistribution, 1, maxDistributionSize);
  free_cvector(naiveMvFlag, 1, mvSize);
  free_cmatrix(naiveSign, 1, mRecordSize, 1, mvSize);
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeCommon() EXIT ...\n");
  }
}
void unImpute (uint mode) {
  double  *statusPtr;
  double  *timePtr;
  double **predictorPtr;
  uint    mRecordSize; 
  uint   *mRecordIndex;
  int   **mvSign;
  uint    mvSize;
  int    *mvIndex;
  double *valuePtr;
  uint unsignedIndex;
  char result;
  uint p, i;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nunImpute() ENTRY ...\n");
  }
  result = FALSE;
  if (_mRecordSize > 0) {
    statusPtr = _status;
    timePtr = _time;
    predictorPtr = _observation;
    mRecordSize = _mRecordSize;
    mRecordIndex = _mRecordIndex;
    mvSize = _mvSize;
    mvSign = _mvSign;
    mvIndex = _mvIndex;
    for (p = 1; p <= mvSize; p++) {
      switch (mvIndex[p]) {
      case CENS_IDX:
        unsignedIndex = abs(mvIndex[p]);
        valuePtr = statusPtr;
        break;
      case TIME_IDX:
        unsignedIndex = abs(mvIndex[p]);
        valuePtr = timePtr;
        break;
      default:
        unsignedIndex = (uint) mvIndex[p] + 2;
        valuePtr = predictorPtr[(uint) mvIndex[p]];
        break;
      }
      for (i = 1; i <= mRecordSize; i++) {
        if (mvSign[unsignedIndex][i] == 1) {
          valuePtr[mRecordIndex[i]] = NA_REAL;
        }
      }
    }
  }
  if (mode != RSF_GROW) {
    if (_fmRecordSize > 0) {
      statusPtr = _fstatus;
      timePtr = _ftime;
      predictorPtr = _fobservation;
      mRecordSize = _fmRecordSize;
      mRecordIndex = _fmRecordIndex;
      mvSize = _fmvSize;
      mvSign = _fmvSign;
      mvIndex = _fmvIndex;
      for (p = 1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          unsignedIndex = abs(mvIndex[p]);
          valuePtr = statusPtr;
          break;
        case TIME_IDX:
          unsignedIndex = abs(mvIndex[p]);
          valuePtr = timePtr;
          break;
        default:
          unsignedIndex = (uint) mvIndex[p] + 2;
          valuePtr = predictorPtr[(uint) mvIndex[p]];
          break;
        }
        for (i = 1; i <= mRecordSize; i++) {
          if (mvSign[unsignedIndex][i] == 1) {
            valuePtr[mRecordIndex[i]] = NA_REAL;
          }
        }
      }
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nunImpute() EXIT ...\n");
  }
}
void imputeMultipleTime (char selectionFlag) {
  double  *outTime;
  double  *outStatus;
  double meanValue;
  double leftDistance, rightDistance;
  uint minimumIndex;
  char result;
  uint i,j;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeMultipleTime() ENTRY ...\n");
  }
  result = FALSE;
  if (_mRecordSize > 0) {
    result = TRUE;
  }
  if (result == FALSE) {
  }
  if (selectionFlag == FALSE) {
    outTime   = _sOOBImputeTimePtr;
    outStatus = _sOOBImputeStatusPtr;
  }
  else {
    outTime  = _sImputeTimePtr;
    outStatus  = _sImputeStatusPtr;
  }    
  for (i=1; i <= _mRecordSize; i++) {
    if(_mvSign[abs(TIME_IDX)][i] == 1) { 
      meanValue = outTime[i];
      if ((meanValue < _masterTime[1]) || (meanValue > _masterTime[_masterTimeSize])) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  The summary mean value for time is out of range:  indv %10d, value %12.4f", _mRecordIndex[i], meanValue);
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
      }
      leftDistance = rightDistance = minimumIndex = 0;
      for (j = 1; j <= _masterTimeSize; j++) {
        if (meanValue <= _masterTime[j]) {
          minimumIndex = j;
          j = _masterTimeSize;
        }
      }
      if (minimumIndex == 1) {
      }
      else {
        leftDistance = meanValue - _masterTime[minimumIndex-1];
        rightDistance = _masterTime[minimumIndex] - meanValue;
        if (leftDistance < rightDistance) {
          minimumIndex = j-1;            
        }
        else {
          if (abs(leftDistance -rightDistance) < EPSILON) {
            if (ran2(_seed2Ptr) <= 0.5) {
              minimumIndex = j-1;
            }
          }
        }
      }
      if (getTraceFlag() & MISS_LOW_TRACE) {
        if ((leftDistance > 0) && (rightDistance > 0)) {
          Rprintf("\nSummary time overridden with nearest neighbour for indv idx:  %10d at master time idx %10d ", _mRecordIndex[i], minimumIndex);
          Rprintf("\nSummary -> NearestNhbr:  %12.4f -> %12.4f  \n", outTime[i], _masterTime[minimumIndex]);
        }
      }
      outTime[i] = _masterTime[minimumIndex];
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nimputeMultipleTime EXIT ...\n");
  }
}
double getMaximalValue(double *value, uint size) {
  double result;
  uint classCount, maximalClassSize, maximalClassCount;
  uint randomIndex;
  uint j;
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetMaximalValue() ENTRY ...\n");
  }
  uint   *classSize  = uivector(1, size);
  for (j = 1; j <= size; j++) {
    classSize[j] = 0;
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMaximal Class Values (raw):  ");
    Rprintf("\n     index       value \n");
    for (j=1; j <= size; j++) {
      Rprintf("%10d %10d \n", j, (int) value[j]);
    }
  }
  hpsort(value, size);
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMaximal Class Values (sorted):  ");
    Rprintf("\n     index       value \n");
    for (j=1; j <= size; j++) {
      Rprintf("%10d %10d \n", j, (int) value[j]);
    }
  }
  classCount = 1;
  classSize[1] = 1;
  for (j = 2; j <= size; j++) {
    if (value[j] > value[classCount]) {
      classCount ++;
      value[classCount] = value[j];
    }
    classSize[classCount] ++;
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMaximal Class Values (summary):  ");
    Rprintf("\n      size      value \n");
    for (j=1; j <= classCount; j++) {
      Rprintf("%10d %10d \n", classSize[j], (int) value[j]);
    }
  }
  maximalClassSize = maximalClassCount = 0;
  for (j=1; j <= classCount; j++) {
    if (classSize[j] > maximalClassSize) {
      maximalClassSize = classSize[j];
    }
  }
  for (j=1; j <= classCount; j++) {
    if (classSize[j] == maximalClassSize) {
      maximalClassCount ++;
    }
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMaximal Class (size, number):  ");
    Rprintf("\n(%10d, %10d) \n", maximalClassSize, maximalClassCount);
  }
  if (maximalClassCount > 1) {
    if (getTraceFlag() & MISS_HGH_TRACE) {
      Rprintf("\nUsing unchained generator to break class size tie. ");
    }
    if (getTraceFlag() & MISS_HGH_TRACE) {
      Rprintf("\nMaximal value random seed for ran2():  %20d", *_seed2Ptr);
    }
    randomIndex = (uint) ceil(ran2(_seed2Ptr)*((maximalClassCount)*1.0));
  }
  else {
    randomIndex = 1;
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMaximal Class Random Index:  ");
    Rprintf("\n%10d \n", randomIndex);
  }
  j = 0;
  while (randomIndex > 0) {
    j++;
    if (classSize[j] == maximalClassSize) {
      randomIndex --;
    }
  }
  result = value[j];
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMaximal Class Value:  ");
    Rprintf("\n%10d \n", (int) result);
  }
  free_uivector(classSize, 1, size);
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetMaximalValue() EXIT ...\n");
  }
  return result;
}
double getMedianValue(double *value, uint size) {
  double result;
  uint medianIndex;
  uint j;
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetMedianValue() ENTRY ...\n");
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMedian Values (raw):  \n");
    for (j=1; j <= size; j++) {
      Rprintf("%10d %12.4f \n", j, value[j]);
    }
  }
  hpsort(value, size);
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMedian Values (sorted):  \n");
    for (j=1; j <= size; j++) {
      Rprintf("%10d %12.4f \n", j, value[j]);
    }
  }
  if (size > 1) {
    medianIndex = (uint) ceil(size/2);
  }
  else {
    medianIndex = 1;
  }
  result = value[medianIndex];
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMedian Index Index:  ");
    Rprintf("\n%10d \n", medianIndex);
    Rprintf("\nMedian Value:  ");
    Rprintf("\n%12.4f \n", result);
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetMedianValue() EXIT ...\n");
  }
  return result;
}
double getMeanValue(double *value, uint size) {
  double result;
  uint j;
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetMeanValue() ENTRY ...\n");
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMean Values (raw):  \n");
    for (j=1; j <= size; j++) {
      Rprintf("%10d %12.4f \n", j, value[j]);
    }
  }
  result = 0.0;
  for (j = 1; j <= size; j++) {
    result = result + value[j];
  }
  result = result / size;
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nMean Value:  \n");
    Rprintf("%12.4f \n", result);
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetMeanValue() EXIT ...\n");
  }
  return result;
}
double getSampleValue(double *value, uint size, char chainFlag) {
  uint randomIndex;
  uint j;
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetSampleValue() ENTRY ...\n");
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\nSample Values (raw):  \n");
    for (j=1; j <= size; j++) {
      Rprintf("%10d %12.4f \n", j, value[j]);
    }
  }
  if (chainFlag) {
    if (getTraceFlag() & MISS_HGH_TRACE) {
      Rprintf("\nSample value random seed for ran1():  %20d", *_seed1Ptr);
    }
    randomIndex = (uint) ceil(ran1(_seed1Ptr)*((size)*1.0));
  }
  else {
    if (getTraceFlag() & MISS_HGH_TRACE) {
      Rprintf("\nSample value random seed for ran2():  %20d", *_seed2Ptr);
    }
    randomIndex = (uint) ceil(ran2(_seed2Ptr)*((size)*1.0));
  }
  if (getTraceFlag() & MISS_HGH_TRACE) {
    Rprintf("\ngetSampleValue() EXIT ...\n");
  }
  return value[randomIndex];
}
uint getRecordMap(uint *map, 
                  uint size, 
                  double *status, 
                  double *time, 
                  double *data) {
  uint i, p;
  uint mSize;
  char mFlag;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\ngetRecordMap() ENTRY ...\n");
  }
  mSize  = 0;
  for (i = 1; i <= size; i++) {
    mFlag = FALSE;
    if (ISNA(status[i])) {
      mFlag = TRUE;
    }
    if (ISNA(time[i])) {
      mFlag = TRUE;
    }
    if (mFlag == FALSE) {
      for (p = 1; p <= _xSize; p++) {
        if (ISNA((data+((p-1)*(size)))[i-1])) {
          mFlag = TRUE;
          p = _xSize;
        }
      }
    }
    if (mFlag == TRUE) {
     mSize ++;
     if (map != NULL) {
       map[i] = mSize;
     }
    }
    else {
      if (map != NULL) {
        map[i] = 0;
      }
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nMapping of Individuals to Missing Record Index:  ");
    Rprintf("\n  orgIndex     mIndex \n");
    for (i = 1; i <= size; i++) {
      if (map[i] != 0) {
        Rprintf("%10d %10d \n", i, map[i]);
      }
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\ngetRecordMap() EXIT ...\n");
  }
  return mSize;
}
char getForestSign (uint mode, uint b) {
  char result;
  uint i,p,q,m;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\ngetForestSign() ENTRY ...\n");
  }
  result = TRUE;
  if (_mRecordSize > 0) {
    int **mvBootstrapSign = imatrix(1, _mvSize, 1, _observationSize);
    for (p = 1; p <= _mvSize; p++) {
      for (i = 1; i <= _observationSize; i++) {
        mvBootstrapSign[p][i] = 0;
      }
    }
    for (p = 1; p <= _mvSize; p++) {
      _mvForestSign[b][p] = 0;
    }
    for (i=1; i <= _observationSize; i++) {
      m = _bootMembershipIndex[i];
      if (_mRecordMap[m] != 0) {
        for (p = 1; p <= _mvSize; p++) {
          switch (_mvIndex[p]) {
          case CENS_IDX:
            mvBootstrapSign[p][i] = _mvSign[(uint) abs(CENS_IDX)][_mRecordMap[m]];
            break;
          case TIME_IDX:
            mvBootstrapSign[p][i] = _mvSign[(uint) abs(TIME_IDX)][_mRecordMap[m]];
            break;
          default:
            mvBootstrapSign[p][i] = _mvSign[(uint) _mvIndex[p]+2][_mRecordMap[m]];
            break;
          }
        }
      }
      else {
        for (p=1; p <= _mvSize; p++) {
          mvBootstrapSign[p][i] = 0;
        }
      }
      for (p = 1; p <= _mvSize; p++) {
        _mvForestSign[b][p] = _mvForestSign[b][p] + mvBootstrapSign[p][i];
      }
    }
    m = 0;
    for (p = 1; p <= _mvSize; p++) {
      if (_mvForestSign[b][p] > 0) {
        if (_mvForestSign[b][p] == _observationSize) {
          _mvForestSign[b][p] = -1;
        }
        else {
          _mvForestSign[b][p] = 1;
        }
      }
      switch (_mvIndex[p]) {
      case CENS_IDX:
        if (_mvForestSign[b][p] == -1) result = FALSE;
        break;
      case TIME_IDX:
        if (_mvForestSign[b][p] == -1) result = FALSE;
        break;
      default:
        if (_mvForestSign[b][p] == -1) m ++;
        break;
      }
    }  
    if (m == _mvSize) {
      result = FALSE;
    }
    if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\n\nGROW Individual Sample Signature For Tree:  %10d", b);
      Rprintf("\n bootIndex   signatures ->\n");
      Rprintf(  "            ");
      for (i=1; i <= _mvSize; i++) {
        Rprintf("%3d", _mvIndex[i]);
      }
      Rprintf("\n");
      for (i=1; i <=  _observationSize; i++) {
        Rprintf("%10d  ", i);
        for (p=1; p <= _mvSize; p++) {
          Rprintf("%3d", mvBootstrapSign[p][i]);
        }
        Rprintf("\n");
      }
    }
    if (getTraceFlag() & MISS_LOW_TRACE) {
      Rprintf("\nGROW Missing Sample Signature For Tree:  %10d \n", b);
      Rprintf(  "            ");
      for (p=1; p <= _mvSize; p++) {
        Rprintf("%3d", _mvIndex[p]);
      }
      Rprintf("\n");
      Rprintf(  "            ");
      for (p=1; p <= _mvSize; p++) {
        Rprintf("%3d", _mvForestSign[b][p]);
      }
      Rprintf("\n");
    }
    free_imatrix(mvBootstrapSign, 1, _mvSize, 1, _observationSize);
  }
  if (mode != RSF_GROW) {
    if (_fmRecordSize > 0) {
      for (p = 1; p <= _fmvSize; p++) {
        _fmvForestSign[b][p] = 1;
      }
      if (_mRecordSize > 0) {
        p = q = 1;
        while ((q <= _fmvSize) && (p <= _mvSize)) {
          if (_mvIndex[p] == _fmvIndex[q]) {
            if (_mvForestSign[b][p] == -1) {
              _fmvForestSign[b][q] = -1;
            }
            p++;
            q++;
          }
          else {
            switch (_fmvIndex[q]) {
            case CENS_IDX:
              if (_mvIndex[p] > 0) {
                q++;
              }
              else {
                if (abs(_mvIndex[q]) < abs(_mvIndex[p])) {
                  q++;
                }
                else {
                  p++;
                }
              }
              break;
            case TIME_IDX:
              if (_mvIndex[p] > 0) {
                q++;
              }
              else {
                if (abs(_mvIndex[q]) < abs(_mvIndex[p])) {
                  q++;
                }
                else {
                  p++;
                }
              }
              break;
            default:
              if (_fmvIndex[q] < _mvIndex[p]) {
                q++;
              }
              else {
                p++;
              }
              break;
            }
          }
        }
      }  
      if (getTraceFlag() & MISS_LOW_TRACE) {
        Rprintf("\nPRED Missing Signature For Tree:  %10d \n", b);
        Rprintf(  " Outcome or Predictor: ");
        for (p=1; p <= _fmvSize; p++) {
          Rprintf("%3d", _fmvIndex[p]);
        }
        Rprintf("\n");
        Rprintf(  "                       ");
        for (p=1; p <= _fmvSize; p++) {
          Rprintf("%3d", _fmvForestSign[b][p]);
        }
        Rprintf("\n");
      }
    }  
  }  
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\ngetForestSign() EXIT ...\n");
  }
  return result;
}
void updateTimeIndexArray(Node *parent) {
  char nodeFlag;
  char idxFoundFlag;
  uint i,k;
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nupdateTimeIndexArray() ENTRY ...\n");
  }
  for (i=1; i <= _observationSize; i++) {
    if (parent != NULL) {    
      nodeFlag = FALSE;
      if (_nodeMembership[i] == parent) {
        nodeFlag = TRUE;
      }
    }
    else {
      nodeFlag = TRUE;
    }
    if (nodeFlag) {
      idxFoundFlag = FALSE;
      if (!ISNA(_time[i])) {
        k = 1;
        while (k <= _masterTimeSize) {
          if (_time[i] == _masterTime[k]) {
            _masterTimeIndex[i] = k;
            idxFoundFlag = TRUE;
            k = _masterTimeSize;
          }
          k++;
        }
      }
      else {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Invalid event time encountered:  %12.4f", _time[i]);
        Rprintf("\nRSF:  Please Contact Technical Support.");
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
      }
      if (idxFoundFlag == FALSE) {
        _masterTimeIndex[i] = 0;
      }
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nMaster Time Index for node:  \n");
    Rprintf("\n      Indv      Index         Time:  \n");
    for (i=1; i <= _observationSize; i++) {
      if (parent != NULL) {    
        nodeFlag = FALSE;
        if (_nodeMembership[i] == parent) {
          nodeFlag = TRUE;
        }
      }
      else {
        nodeFlag = TRUE;
      }
      if (nodeFlag) {
        if(_masterTimeIndex[i] > 0) {
          Rprintf("%10d %10d %12.4f \n", i, _masterTimeIndex[i], _masterTime[_masterTimeIndex[i]]);
        }
        else {
          Rprintf("%10d %10d           XX \n", i, _masterTimeIndex[i]);
        }
      }
    }
  }
  if (getTraceFlag() & MISS_LOW_TRACE) {
    Rprintf("\nupdateTimeIndexArray() EXIT ...\n");
  }
}
void updateEventTypeSubsets(double *summaryStatus, 
                            uint    mRecordSize,
                            int   **mvSign,
                            uint   *mRecordIndex) {
  uint i, j, n;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEventTypeSubsets() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  if (_eventTypeSize == 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to update event type subsets in a non-CR analysis.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (_mStatusSize > 0) {
    uint *eventCounter = uivector(1, _eventTypeSize);
    for (j = 1; j <= _eventTypeSize; j++) {
      eventCounter[j] = _eIndividualSize[j]; 
    }
    for (i = 1; i <= mRecordSize; i++) {
      if (mvSign[(uint) abs(CENS_IDX)][i] == 1) {
        if ((uint) summaryStatus[mRecordIndex[i]] > 0) {
          j = _eventTypeIndex[(uint) summaryStatus[mRecordIndex[i]]];
          eventCounter[j] ++;
          _eIndividual[j][eventCounter[j]] = mRecordIndex[i];
        }
        else {
          for (j=1; j <= _eventTypeSize; j++) {
            eventCounter[j] ++;
            _eIndividual[j][eventCounter[j]] = mRecordIndex[i];
          }
        }
      }
    }
    for (j = 1; j <= _eventTypeSize; j++) {
      _meIndividualSize[j] = eventCounter[j];
    }
    free_uivector(eventCounter, 1, _eventTypeSize);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nEvent Type Subset Sizes:  \n");
    for (j=1; j <= _eventTypeSize; j++) {
      Rprintf("%10d %10d %10d \n", j, _eventType[j], _meIndividualSize[j]);
    }
    Rprintf("\nEvent Type Subsets:  \n");
    Rprintf("          ");
    for (n=1; n <= _observationSize; n++) {
      Rprintf("%10d", n);
    }
    Rprintf("\n");
    for (j=1; j <= _eventTypeSize; j++) {
      Rprintf("%10d", j);
      for (n=1; n <= _meIndividualSize[j]; n++) {
        Rprintf("%10d", _eIndividual[j][n]);
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _censblTime = _censblTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEventTypeSubsets() EXIT ...\n");  
  }
}
