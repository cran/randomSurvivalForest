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
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include     "rsfImpute.h"
#include  "rsfBootstrap.h"
char bootstrap (uint     mode,
                uint     b,
                Node    *rootPtr,
                char   **dmRecordBootFlag) {
  char mOutcomeFlag;
  uint ibgSampleSize;
  char result;
  uint i,k;
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\nbootstrap() ENTRY ...\n");
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
    Rprintf("\nRSF:  Bootstrap sample:  %10d ", b);  
  }
  result = TRUE;
  mOutcomeFlag = FALSE;
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\nBootstrap random seed chain ran1():  ");
    Rprintf("\n%10d %20d ", b, *_seed1Ptr);
  }
  for (i=1; i <= _observationSize; i++) {
    _bootMembershipFlag[i] = FALSE;
  }
  for (i=1; i <= _observationSize; i++) {
    k = (uint) ceil(ran1(_seed1Ptr)*((_observationSize)*1.0));
    _bootMembershipFlag[k] = TRUE;
    _bootMembershipIndex[i] = k;
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
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
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\n\nIBG Size:  %10d ", ibgSampleSize);
    Rprintf(  "\nOOB Size:  %10d ", _oobSampleSize[b]);
  }
  result = getForestSign(mode, b);
  if (result == FALSE) {
    if (getTraceFlag() & SUMM_USR_TRACE) {
      Rprintf("\nRSF:  Status or Time values are all missing in the sample.  Bootstrap sample has been discarded.");  
    }
  }
  if (getTraceFlag() & SUMM_USR_TRACE) {
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
        if (_bootMembershipFlag[_intrIndividual[i]] == FALSE) {
          _foobSampleSize[b] ++;
        }
      }
      if (getTraceFlag() & SUMM_HGH_TRACE) {
        Rprintf(  "\nINTR OOB Size:  %10d ", _foobSampleSize[b]);
      }
      if (_fmRecordSize > 0) {
        for (i = 1; i <= _fmRecordSize; i++) {
          if (_bootMembershipFlag[_intrIndividual[_fmRecordIndex[i]]] == TRUE) {
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
    if (getTraceFlag() & SUMM_HGH_TRACE) {
      Rprintf("\nBootstrap sample is invalid:  %10d \n", b);      
    }
  }
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\nbootstrap() EXIT ...\n");
  }
  return result;
}
