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

#include   "global.h"
#include   "nrutil.h"
#include   "rsfImpute.h"
#include   "rsfBootstrap.h"
extern uint getTraceFlag();
char bootstrap (char     multipleImputeFlag,
                uint     mode,
                uint     b,
                Node    *rootPtr,
                double **masterSplit,
                uint    *masterSplitSize,
                char   **dmRecordBootFlag) {
  char mOutcomeFlag;
  uint **masterSplitBounds;
  uint *leadIndex;
  uint ibgSampleSize;
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
    Rprintf("\nStart of forest seed chain for bootstrap:  ");
    Rprintf("\n%10d %20d ", b, *_seed1Ptr);
  }
  for (i=1; i <= _observationSize; i++) {
    _bootMembershipFlag[i] = FALSE;
  }
  for (i=1; i <= _observationSize; i++) {
    k = ceil(ran1(_seed1Ptr)*((_observationSize)*1.0));
    _bootMembershipFlag[k] = TRUE;
    _bootMembershipIndex[i] = k;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nIn-Bag Membership:  ");
    Rprintf("\n bootIndex   orgIndex (in original data) \n");
    for (i=1; i <=  _observationSize; i++) {
      Rprintf("%10d %10d \n", i, _bootMembershipIndex[i]);
    } 
  }
  for (i=1; i <= _observationSize; i++) {
    _nodeMembership[i] = rootPtr;
  }
  _oobSampleSize[b] = 0;
  for (i=1; i <= _observationSize; i++) {
    _nodeMembership[i] = rootPtr;
    if (_bootMembershipFlag[i] == FALSE) {
      _oobSampleSize[b] ++;
    }
  }
  ibgSampleSize = _observationSize - _oobSampleSize[b];
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nIBG Size:  %10d ", ibgSampleSize);
    Rprintf(  "\nOOB Size:  %10d ", _oobSampleSize[b]);
  }
  result = getForestSign(mode, b);
  if (result == FALSE) {
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Status or Time values are all missing in the sample.  Bootstrap sample has been discarded.");  
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Bootstrapping and all mode initialization complete.");  
  }
  if (result == TRUE) {
    if (mode == RSF_GROW) {
      if (_mRecordSize > 0) {
        for (i = 1; i <= _mRecordSize; i++) {
          if (_bootMembershipFlag[_mRecordIndex[i]] == TRUE) {
            dmRecordBootFlag[b][i] = TRUE;
          }
          else {
            dmRecordBootFlag[b][i] = FALSE;
          }
        }
      }  
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
            if (multipleImputeFlag == TRUE) {
              leadIndex[j]++;
              masterSplit[j][leadIndex[j]] = _observation[j][i];
            }
            else {
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
      }  
      if (getTraceFlag() & DL3_TRACE) {
        Rprintf("\n\nIn-Bag Raw MasterSplit Predcitor Data:  ");
        Rprintf("\n          ");
        for (j=1; j <= _xSize; j++) {
          Rprintf(" %12d", j);
        }
        for (i=1; i <= ibgSampleSize; i++) {
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
        for (i=1; i <= ibgSampleSize; i++) {
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
      nrCopyMatrix(rootPtr -> permissibleSplit, masterSplitBounds, _xSize, 2);
      free_uivector(leadIndex, 1, _xSize);
      free_uimatrix(masterSplitBounds, 1, _xSize, 1, 2);
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Initialization of master split data complete.");  
      }
    }  
    if (mode == RSF_PRED) {
      for (i=1; i <= _fobservationSize; i++) {
        _fnodeMembership[i] = rootPtr;
      }
    }
    if (mode == RSF_INTR) {
      for (i=1; i <= _fobservationSize; i++) {
        _fnodeMembership[i] = rootPtr;
      }
      _foobSampleSize[b] = 0;
      for (i=1; i <= _fobservationSize; i++) {
        if (_bootMembershipFlag[_intrObservation[i]] == FALSE) {
          _foobSampleSize[b] ++;
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf(  "\nINTR OOB Size:  %10d ", _foobSampleSize[b]);
      }
      if (_fmRecordSize > 0) {
        for (i = 1; i <= _fmRecordSize; i++) {
          if (_bootMembershipFlag[_intrObservation[_fmRecordIndex[i]]] == TRUE) {
            dmRecordBootFlag[b][i] = TRUE;
          }
          else {
            dmRecordBootFlag[b][i] = FALSE;
          }
        }
      }
    }
  }  
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nBootstrap sample is invalid:  %10d \n", b);      
    }
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nbootstrap() EXIT ...\n");
  }
  return result;
}
