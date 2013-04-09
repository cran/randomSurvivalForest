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
#include     "rsfImpute.h"
#include  "rsfBootstrap.h"
char bootstrap (uint     mode,
                uint     b,
                Node    *rootPtr,
                char   **dmRecordBootFlag) {
  char result;
  uint i,k;
  result = TRUE;
  for (i=1; i <= _observationSize; i++) {
    _bootMembershipFlag[i] = FALSE;
  }
  for (i=1; i <= _observationSize; i++) {
    k = (uint) ceil(ran1(_seed1Ptr)*((_observationSize)*1.0));
    _bootMembershipFlag[k] = TRUE;
    _bootMembershipIndex[i] = k;
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
  result = getForestSign(mode, b);
  if (result == FALSE) {
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
  }
  return result;
}
