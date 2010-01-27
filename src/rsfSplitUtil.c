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
#include  "rsfFactorOps.h"
#include  "rsfSplitUtil.h"
void updateMaximumSplit(double  delta, 
                        uint    randomCovariate,
                        uint    index,
                        char    factorFlag,
                        uint    mwcpSizeAbsolute,
                        double *deltaMax,
                        uint   *splitParameterMax,
                        void   *permissibleSplitPtr) {
  uint k;
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nupdateMaximumSplit() ENTRY ...\n");
  }
  if (delta > *deltaMax) {
    *deltaMax = delta;
    *splitParameterMax = randomCovariate;
    if (getTraceFlag() & SPLT_MED_TRACE) {
      Rprintf("\n\nUpdated Running Split Statistics: \n");
      Rprintf("  SplitParm  SplitValIdx        Delta \n");
      Rprintf(" %10d %12d %12.4f \n", randomCovariate, index, *deltaMax);
    }
    if (factorFlag == TRUE) {
      if (_splitValueMaxFactSize > 0) {
        if (_splitValueMaxFactSize != mwcpSizeAbsolute) {
          free_uivector(_splitValueMaxFactPtr, 1, _splitValueMaxFactSize);
          _splitValueMaxFactSize = mwcpSizeAbsolute;
          _splitValueMaxFactPtr = uivector(1, _splitValueMaxFactSize);
        }
      }
      else {
        _splitValueMaxFactSize = mwcpSizeAbsolute;
        _splitValueMaxFactPtr = uivector(1, _splitValueMaxFactSize);
      }
      _splitValueMaxCont = NA_REAL;
      for (k=1; k <= _splitValueMaxFactSize; k++) {
        _splitValueMaxFactPtr[k] = 
          ((uint*) permissibleSplitPtr + ((index - 1) * _splitValueMaxFactSize))[k];
      }
      if (getTraceFlag() & SPLT_MED_TRACE) {
        Rprintf(" at MWCPsize= %2d, mwcp= ", _splitValueMaxFactSize);
        for (k = _splitValueMaxFactSize; k >= 1; k--) {
          Rprintf("%8x ", _splitValueMaxFactPtr[k]);
        }
        Rprintf("\n");
      }
    }
    else {
      if (_splitValueMaxFactSize > 0) {
        free_uivector(_splitValueMaxFactPtr, 1, _splitValueMaxFactSize);
        _splitValueMaxFactSize = 0;
        _splitValueMaxFactPtr = NULL;
      }
      else {
      }
      _splitValueMaxCont = ((double*) permissibleSplitPtr)[index];
      if (getTraceFlag() & SPLT_MED_TRACE) {
        Rprintf(" at %12.4f \n", _splitValueMaxCont);
      }
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nupdateMaximumSplit() EXIT ...\n");
  }
}
uint stackAndSelectRandomCovariates(Node     *parent,
                                    uint      nodeSize,
                                    uint     *nodeIndex,
                                    uint    **covariateIndex,
                                    double ***permissibleSplit,
                                    uint    **permissibleSplitSize) {
  uint i;
  uint actualCovariateCount;
  uint candidateCovariate;
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nstackAndSelectRandomCovariates() ENTRY ...\n");
  }
  *covariateIndex = uivector(1, _xSize);
  *permissibleSplit = dmatrix(1, _randomCovariateCount, 1, nodeSize);
  *permissibleSplitSize = uivector(1, _randomCovariateCount);
  char *randomSplitVector = cvector(1, _xSize);
  if (nodeSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid nodeSize encountered in stackAndSelectRandomCovariates():  %10d", nodeSize);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  nrCopyVector(randomSplitVector, parent -> permissibleSplit, _xSize);
  for(i=1; i <= _randomCovariateCount; i++) {
    (*covariateIndex)[i] = 0;
  }
  actualCovariateCount =  1;
  candidateCovariate   = -1;
  while ((actualCovariateCount  <= _randomCovariateCount) && (candidateCovariate != 0)) {
    candidateCovariate = getSelectableElement(_xSize, randomSplitVector, _randomCovariateWeight);
    if (candidateCovariate != 0) {
      if (getTraceFlag() & SPLT_HGH_TRACE) {
        Rprintf("\nCandidate covariate (index, candidate):  %10d %10d \n", actualCovariateCount, candidateCovariate);
      }
      for (i=1; i <= nodeSize; i++) {
        (*permissibleSplit)[actualCovariateCount][i] = _observation[candidateCovariate][nodeIndex[i]];
      }
      (*permissibleSplitSize)[actualCovariateCount] = 1;
      hpsort((*permissibleSplit)[actualCovariateCount], nodeSize);
      for (i = 2; i <= nodeSize; i++) {
        if ((*permissibleSplit)[actualCovariateCount][i] > (*permissibleSplit)[actualCovariateCount][(*permissibleSplitSize)[actualCovariateCount]]) {
          (*permissibleSplitSize)[actualCovariateCount] ++;
          (*permissibleSplit)[actualCovariateCount][(*permissibleSplitSize)[actualCovariateCount]] = (*permissibleSplit)[actualCovariateCount][i];
        }
      }
      if (getTraceFlag() & SPLT_HGH_TRACE) {
        Rprintf("\nCandidate covariate size:  %10d \n", (*permissibleSplitSize)[actualCovariateCount]);
      }
      if((*permissibleSplitSize)[actualCovariateCount] >= 2) {
        randomSplitVector[candidateCovariate] = ACTIVE;
        (*covariateIndex)[actualCovariateCount] = candidateCovariate;
        if (getTraceFlag() & SPLT_HGH_TRACE) {
          Rprintf("\nCovariate selected:  %10d \n", candidateCovariate);
          Rprintf("\nSplit points:       index        value");
          for (i = 1; i <= (*permissibleSplitSize)[actualCovariateCount]; i++) {
            Rprintf("\n               %10d %12.4f", i, (*permissibleSplit)[actualCovariateCount][i]);
          }
        }
        actualCovariateCount ++;
      }
      else {
        (parent -> permissibleSplit)[candidateCovariate] = FALSE;
        randomSplitVector[candidateCovariate] = FALSE;
        if (getTraceFlag() & SPLT_HGH_TRACE) {
          Rprintf("\nCovariate rejected:  %10d \n", candidateCovariate);
        }
      }
    }
  }
  actualCovariateCount --;
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nCovariate Random Selection:  \n");
    for (i=1; i <= actualCovariateCount; i++) {
      Rprintf("%10d %10d \n", i, (*covariateIndex)[i]);
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nPermissible Split Sizes: \n");
    Rprintf(" covariate       size \n");
    for (i=1; i <= actualCovariateCount; i++) {
      Rprintf("%10d %10d \n", (*covariateIndex)[i], (*permissibleSplitSize)[i]);
    }
  }
  free_cvector(randomSplitVector, 1, _xSize);
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nstackAndSelectRandomCovariates() EXIT ...\n");
  }
  return actualCovariateCount;
}
void unstackRandomCovariates(uint     nodeSize, 
                             uint    *covariateIndex,
                             double **permissibleSplit,
                             uint    *permissibleSplitSize) {
  free_uivector(covariateIndex, 1, _xSize);
  free_dmatrix(permissibleSplit, 1, nodeSize, 1, _randomCovariateCount);
  free_uivector(permissibleSplitSize, 1, _randomCovariateCount);
}
uint getSelectableElement (uint    length,
                           char   *permissible,
                           double *weight) {
  char   *localPermissible = NULL;  
  double *cdf = NULL;  
  uint selectableCount;
  uint covariateIndex;
  double randomValue;
  uint i, j, k, p, index;
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\ngetSelectableElement() ENTRY ...\n");
  }
  if (length > 0) {
    localPermissible = cvector(1, length);
    cdf = dvector(1, length);
  }
  selectableCount = 0;
  for (i=1; i <= length; i++) {
    if (permissible[i] == TRUE) {
      if (weight != NULL) {
        if (weight[i] > 0) {
          localPermissible[i] = TRUE;
          selectableCount ++;
        }
        else {
          localPermissible[i] = FALSE;
        }
      }
      else {
        localPermissible[i] = TRUE;
        selectableCount ++;
      }
    }
    else {
      localPermissible[i] = FALSE;
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nSelectable vs Total Count: \n");  
    Rprintf("%10d %10d \n", selectableCount, length);
  }
  if (selectableCount > 0) {
    if (weight != NULL) { 
      covariateIndex = 0;
      for (k=1; k <= _xSize; k++) {
        if (localPermissible[k] == TRUE) {
          cdf[++covariateIndex] = weight[k];
        }
      }
      for (k=2; k <= covariateIndex; k++) {
        cdf[k] += cdf[k-1];
      }
      if (getTraceFlag() & SPLT_HGH_TRACE) {
        Rprintf("\nUpdated CDF based on weights:  ");
        Rprintf("\n     index          cdf");
        for (k=1; k <= covariateIndex; k++) {
          Rprintf("\n%10d  %12.4f", k, cdf[k]);
        }
        Rprintf("\n");
      }
      randomValue = ran2(_seed2Ptr)*cdf[covariateIndex];
      j=1;
      while (randomValue > cdf[j]) {
        j++;
      }
      for (index = 1; j > 0; index++) {
        if (localPermissible[index] == TRUE) {
          j--;
        }
      }
      index --;
    }
    else {
      p = (uint) ceil(ran2(_seed2Ptr)*(selectableCount*1.0));
      index = 1;
      while (p > 0) {
        if (permissible[index] == TRUE) {
          p --;
        }
        index ++;
      }
      index --;
    }
  }  
  else {
    index = 0;
  }
  if (length > 0) {
    free_cvector(localPermissible, 1, length);
    free_dvector(cdf, 1, length);
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nSelected Index:  %10d", index);
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\ngetSelectableElement() EXIT ...\n");
  }
  return index;
}
char getDeathCount(Node *parent, 
                   uint *localMembershipIndex, 
                   uint *localDeathTimeCount, 
                   uint *localDeathTimeIndex,
                   uint *localMembershipSize,
                   uint *localDeathTimeSize) {
  uint parentDeathCount = 0;
  uint i;
  char result;
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\ngetDeathCount() ENTRY ...\n");
  }
  *localMembershipSize = 0;
  *localDeathTimeSize  = 0;
  result = FALSE;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++(*localMembershipSize)] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] > 0) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount ++;
      }
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for Parent Node: \n");
    for (i=1; i <= (*localMembershipSize); i++) {
      Rprintf("%10d %10d %10d %10d %12.4f \n", i, localMembershipIndex[i], _masterTimeIndex[localMembershipIndex[i]], (uint) _status[localMembershipIndex[i]], _time[localMembershipIndex[i]]);
    }
    Rprintf("\nLocal Death Time Counts:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10d %12.4f \n", i, localDeathTimeCount[i], _masterTime[i]);
    }
  }
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeIndex[++(*localDeathTimeSize)] = i;
      }
    }
    if (getTraceFlag() & SPLT_HGH_TRACE) {
      Rprintf("\nLocal Death Times (i, _masterTimeIndex): \n");
      for (i=1; i <= (*localDeathTimeSize); i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if ((*localDeathTimeSize) >= _minimumDeathCount) {
      result = TRUE;
    }
    else {
      if (getTraceFlag() & SPLT_HGH_TRACE) {
        Rprintf("\nMinimum unique deaths not acheived:  %10d versus %10d ", *localDeathTimeSize, _minimumDeathCount);
        Rprintf("\nNode will not be split.  \n");
      }
    }
  }
  else {
    if (getTraceFlag() & SPLT_HGH_TRACE) {
      Rprintf("\nLess than twice the minimum number of deaths encountered.  ");
      Rprintf("\nNode will not be split.  \n");
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\ngetDeathCount() EXIT ...\n");
  }
  return result;
}
void stackSplit(uint **localMembershipIndex, 
                uint **localDeathTimeCount, 
                uint **localDeathTimeIndex) {
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nstackSplit() ENTRY ...\n");
  }
  *localMembershipIndex = uivector(1, _observationSize);
  *localDeathTimeCount = uivector(1, _masterTimeSize);
  *localDeathTimeIndex = uivector(1, _masterTimeSize);
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nstackSplit() EXIT ...\n");
  }
}
void unstackSplit(uint *localMembershipIndex, 
                  uint *localDeathTimeCount, 
                  uint *localDeathTimeIndex) {
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nunstackSplit() ENTRY ...\n");
  }
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, _masterTimeSize);
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nunstackSplit() EXIT ...\n");
  }
}
void stackSplitCompact(uint   deathTimeSize,
                       uint **nodeParentDeath,
                       uint **nodeParentAtRisk,
                       uint **nodeLeftDeath,
                       uint **nodeLeftAtRisk,
                       uint **nodeRightDeath,
                       uint **nodeRightAtRisk,
                       uint   nodeSize,
                       char **localSplitIndicator) {
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nstackSplitCompact() ENTRY ...\n");
  }
  *nodeParentDeath  = uivector(1, deathTimeSize);
  *nodeParentAtRisk = uivector(1, deathTimeSize);
  *nodeLeftDeath  = uivector(1, deathTimeSize);
  *nodeLeftAtRisk = uivector(1, deathTimeSize);
  *nodeRightDeath  = uivector(1, deathTimeSize);
  *nodeRightAtRisk = uivector(1, deathTimeSize);
  if ((nodeSize > 0) && (localSplitIndicator != NULL)) {
    *localSplitIndicator = cvector(1, nodeSize);
  } 
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nstackSplitCompact() EXIT ...\n");
  }
}
void unstackSplitCompact(uint  deathTimeSize,
                         uint *nodeParentDeath,
                         uint *nodeParentAtRisk,
                         uint *nodeLeftDeath,
                         uint *nodeLeftAtRisk,
                         uint *nodeRightDeath,
                         uint *nodeRightAtRisk,
                         uint  nodeSize,
                         char *localSplitIndicator) {
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nunstackSplitCompact() ENTRY ...\n");
  }
  free_uivector(nodeParentDeath, 1, deathTimeSize);
  free_uivector(nodeParentAtRisk, 1, deathTimeSize);
  free_uivector(nodeLeftDeath, 1, deathTimeSize);
  free_uivector(nodeLeftAtRisk, 1, deathTimeSize);
  free_uivector(nodeRightDeath, 1, deathTimeSize);
  free_uivector(nodeRightAtRisk, 1, deathTimeSize);
  if ((nodeSize > 0) && (localSplitIndicator != NULL)) {
    free_cvector(localSplitIndicator, 1, nodeSize);
  } 
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nunstackSplitCompact() EXIT ...\n");
  }
}
void getAtRisk(uint *localMembershipIndex,
               uint *localDeathTimeCount,
               uint *localDeathTimeIndex,
               uint  localMembershipSize,
               uint  localDeathTimeSize,
               uint *nodeParentDeath,
               uint *nodeParentAtRisk) {
  uint i, j;
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\ngetAtRisk() ENTRY ...\n");
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nLocal Death Counts (timIdx, deaths): \n");
  }
  for (i=1; i <= localDeathTimeSize; i++) {
    nodeParentAtRisk[i] = 0;
    nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
    if (getTraceFlag() & SPLT_HGH_TRACE) {
      Rprintf("%10d %10d \n", i, nodeParentDeath[i]);
    }
    for (j=1; j <= localMembershipSize; j++) {
      if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
        nodeParentAtRisk[i] ++;
      }
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nLocal At Risk Counts (timIdx, at risk): \n");
    for (i=1; i <= localDeathTimeSize; i++) {
      Rprintf("%10d %10d \n", i, nodeParentAtRisk[i]);
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\ngetAtRisk() EXIT ...\n");
  }
}
uint stackAndConstructSplitVector (uint     localMembershipSize,
                                   uint     randomCovariateIndex,
                                   double  *permissibleSplit,
                                   uint     permissibleSplitSize,
                                   char    *factorFlag,
                                   char    *deterministicSplitFlag,
                                   uint    *mwcpSizeAbsolute,
                                   void   **permissibleSplitPtr) {
  uint j, k, j2, k2;
  uint factorSizeAbsolute;
  uint offset;
  uint splitLength;
  uint relativePair;
  if (getTraceFlag() & SPLT_MED_TRACE) {
    Rprintf("\nstackAndConstructSplitVector() ENTRY ...\n");
  }
  splitLength = 0;  
  (*permissibleSplitPtr) = NULL;  
  if (getTraceFlag() & SPLT_MED_TRACE) {
    Rprintf("\nSplitting on (parameter, of size):  %10d %10d \n", randomCovariateIndex, permissibleSplitSize);
  }
  if (strcmp(_xType[randomCovariateIndex], "C") == 0) {
    *factorFlag = TRUE;
    if(_factorList[permissibleSplitSize] == NULL) {
      _factorList[permissibleSplitSize] = makeFactor(permissibleSplitSize, FALSE);
    }
    factorSizeAbsolute = _factorSize[_factorMap[randomCovariateIndex]];
    *mwcpSizeAbsolute = _factorList[factorSizeAbsolute] -> mwcpSize;
    if (getTraceFlag() & SPLT_MED_TRACE) {
      Rprintf("\n(Absolute Factor Size, Absolute MWCP Size):  (%10d, %10d)", factorSizeAbsolute, *mwcpSizeAbsolute);
    }
    if (_splitRule == RANDOM_SPLIT) {
      splitLength = 1 + ((_splitRandomRule <= localMembershipSize) ? _splitRandomRule : localMembershipSize);
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(_splitRandomRule == 0) {
        *deterministicSplitFlag = TRUE;
        if ((_factorList[permissibleSplitSize] -> r) > MAX_EXACT_LEVEL) {
          *deterministicSplitFlag = FALSE;
        }
        else {
          if ( *((uint *) _factorList[permissibleSplitSize] -> complementaryPairCount) >= localMembershipSize ) {
            *deterministicSplitFlag = FALSE;
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = localMembershipSize + 1;
          if (getTraceFlag() & SPLT_MED_TRACE) {
            if (permissibleSplitSize <= MAX_EXACT_LEVEL) {
              Rprintf("\nFactor override to random (pSplit, nSplit, ndSize):  (%10d, %10d, %10d) \n", 
                      *((uint*) _factorList[permissibleSplitSize] -> complementaryPairCount), 
                      _splitRandomRule,
                      localMembershipSize);
            }
            else {
              Rprintf("\nFactor override to random (pSplit, nSplit, ndSize):  (%24.0f, %10d, %10d) \n", 
                      *((double*) _factorList[permissibleSplitSize] -> complementaryPairCount), 
                      _splitRandomRule,
                      localMembershipSize);
            }
          }
        }
        else {
          splitLength = *((uint*) _factorList[permissibleSplitSize] -> complementaryPairCount) + 1;
        }
      }
      else {
        *deterministicSplitFlag = FALSE;
        if ((_factorList[permissibleSplitSize] -> r) <= MAX_EXACT_LEVEL) {
          if (*((uint*) _factorList[permissibleSplitSize] -> complementaryPairCount) <= ((_splitRandomRule <= localMembershipSize) ? _splitRandomRule : localMembershipSize)) {
            splitLength = *((uint*) _factorList[permissibleSplitSize] -> complementaryPairCount) + 1;
            *deterministicSplitFlag = TRUE;
            if (getTraceFlag() & SPLT_MED_TRACE) {
              Rprintf("\nFactor override to determ (pSplit, nSplit, ndSize):  (%10d, %10d, %10d) \n", 
                      *((uint*) _factorList[permissibleSplitSize] -> complementaryPairCount), 
                      _splitRandomRule,
                      localMembershipSize);
            }
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = 1 + ((_splitRandomRule <= localMembershipSize) ? _splitRandomRule : localMembershipSize);
        }
      }  
    }  
    (*permissibleSplitPtr) = uivector(1, splitLength * (*mwcpSizeAbsolute));
    for (offset = 1; offset <= *mwcpSizeAbsolute; offset++) {
      ((uint*) (*permissibleSplitPtr) + ((splitLength - 1) * (*mwcpSizeAbsolute)))[offset] = 0;
    }
    if (*deterministicSplitFlag) {
      bookFactor(_factorList[permissibleSplitSize]);
      j2 = 0;
      for (j = 1; j <= _factorList[permissibleSplitSize] -> cardinalGroupCount; j++) {
        for (k2 = 1; k2 <= ((uint*) _factorList[permissibleSplitSize] -> cardinalGroupSize)[j]; k2++) {
          ++j2;
          relativePair = (_factorList[permissibleSplitSize] -> cardinalGroupBinary)[j][k2];
          convertRelToAbsBinaryPair(permissibleSplitSize, 
                                    factorSizeAbsolute, 
                                    relativePair,
                                    permissibleSplit,
                                    (uint*) (*permissibleSplitPtr) + ((j2 - 1) * (*mwcpSizeAbsolute)));
        }
      }
    }  
    else {
      for (j = 1; j < splitLength; j++) {
        getRandomPair(permissibleSplitSize, factorSizeAbsolute, permissibleSplit, (uint*) (*permissibleSplitPtr) + ((j - 1) * (*mwcpSizeAbsolute)));
      }
    }
  }  
  else {
    *factorFlag = FALSE;
    if (_splitRule == RANDOM_SPLIT) {
      splitLength = 1 + ((_splitRandomRule <= localMembershipSize) ? _splitRandomRule : localMembershipSize);
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(_splitRandomRule == 0) {
        splitLength = permissibleSplitSize;
        (*permissibleSplitPtr) = permissibleSplit;
        *deterministicSplitFlag = TRUE;
      }
      else {
        if (permissibleSplitSize <= _splitRandomRule) {
          splitLength = permissibleSplitSize;
          (*permissibleSplitPtr) = permissibleSplit;
          *deterministicSplitFlag = TRUE;
          if (getTraceFlag() & SPLT_MED_TRACE) {
            Rprintf("\nContinuous override to determ (pSplit, nSplit, ndSize):  (%10d, %10d, %10d) \n", 
                    permissibleSplitSize,
                    _splitRandomRule,
                    localMembershipSize);
          }
        }
        else {
          splitLength = _splitRandomRule + 1;
          *deterministicSplitFlag = FALSE;
        }
      }  
    }  
    if (*deterministicSplitFlag == FALSE) {
      (*permissibleSplitPtr) = dvector(1, splitLength);
      ((double*) (*permissibleSplitPtr))[splitLength] = 0;
      for (j = 1; j < splitLength; j++) {
        k = permissibleSplitSize - 1;
        ((double*) (*permissibleSplitPtr))[j]  = permissibleSplit[(uint) ceil(ran2(_seed2Ptr)*(k*1.0))];
      }
    }  
  }  
  if (getTraceFlag() & SPLT_MED_TRACE) {
    Rprintf("\nstackAndConstructSplitVector(%10d) EXIT ...\n", splitLength);
  }
  return splitLength;
}
void unstackSplitVector(uint   permissibleSplitSize,
                        uint   splitLength,
                        char   factorFlag,
                        char   deterministicSplitFlag,
                        void  *permissibleSplitPtr) {
  if (factorFlag == TRUE) {
    free_uivector(permissibleSplitPtr, 1, splitLength * (_factorList[permissibleSplitSize] -> mwcpSize));
    if (deterministicSplitFlag == FALSE) {
      if (permissibleSplitSize > SAFE_FACTOR_SIZE) {
        unBookFactor(_factorList[permissibleSplitSize]);
      }
    }
  }
  else {
    if (deterministicSplitFlag == FALSE) {
      free_dvector(permissibleSplitPtr, 1, splitLength);
    }
  }
}
void virtuallySplitNode (uint  localMembershipSize,
                         char  factorFlag,
                         uint  mwcpSizeAbsolute,
                         uint  randomCovariate,
                         uint *localMembershipIndex,
                         void *permissibleSplitPtr,
                         uint  offset,
                         uint  localDeathTimeSize,
                         uint *localDeathTimeIndex,
                         uint *nodeParentAtRisk,
                         uint *nodeParentDeath,
                         uint *nodeLeftAtRisk,
                         uint *nodeLeftDeath,
                         uint *leftDeathTimeSize,
                         uint *nodeRightAtRisk,
                         uint *nodeRightDeath,
                         uint *rightDeathTimeSize,
                         char *localSplitIndicator) {
  char daughterFlag;
  uint k, m, index;
  if (getTraceFlag() & SPLT_MED_TRACE) {
    Rprintf("\nvirtuallySplitNode() ENTRY ...\n");
  }
  if (getTraceFlag() & SPLT_MED_TRACE) {
    if (factorFlag == TRUE) {
      Rprintf("\nSplitting on (index, factor (mwcp)):  ");
      for (k = mwcpSizeAbsolute; k >= 1; k--) {
        Rprintf("( %10d, %8x)", offset, ((uint*) permissibleSplitPtr + ((offset - 1) * mwcpSizeAbsolute))[k]);
      }
    }
    else {
      Rprintf("\nSplitting on (index, value):  ( %10d, %12.4f) \n", offset, ((double*) permissibleSplitPtr)[offset]);
    }
  }
  *leftDeathTimeSize = *rightDeathTimeSize = 0;
  for (k=1; k <= localDeathTimeSize; k++) {
    nodeLeftDeath[k] = nodeLeftAtRisk[k] = 0;
  }
  for (k=1; k <= localMembershipSize; k++) {
    daughterFlag = RIGHT;
    if (factorFlag == TRUE) {
      daughterFlag = splitOnFactor((uint) _observation[randomCovariate][localMembershipIndex[k]], (uint*) permissibleSplitPtr + ((offset - 1) * mwcpSizeAbsolute));
    }
    else {
      if (_observation[randomCovariate][localMembershipIndex[k]] <= ((double*) permissibleSplitPtr)[offset]) {
        daughterFlag = LEFT;
      }
    }
    if (localSplitIndicator != NULL) {
      localSplitIndicator[k] = daughterFlag;
    }
    if (daughterFlag == LEFT) {
      if (getTraceFlag() & SPLT_HGH_TRACE) {
        Rprintf("\nMember of LEFT Daughter (index):  %10d %10d \n", k, localMembershipIndex[k]);
      }
      index = 0;  
      for (m = 1; m <= localDeathTimeSize; m++) {
        if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
          if (getTraceFlag() & TURN_OFF_TRACE) {
            Rprintf("NLAR:  %10d %10d %10d \n", m, localDeathTimeIndex[m], _masterTimeIndex[localMembershipIndex[k]]);
          }
          nodeLeftAtRisk[m] ++;
          index = m;
        }
        else {
          m = localDeathTimeSize;
        }
      }
      if (_status[localMembershipIndex[k]] > 0) {
        nodeLeftDeath[index] ++;
      }
    }  
    else {
      if (getTraceFlag() & SPLT_HGH_TRACE) {
        Rprintf("\nMember of RGHT Daughter (index):  %10d %10d \n", k, localMembershipIndex[k]);
      }
    }
  }  
  for (k=1; k <= localDeathTimeSize; k++) {
    nodeRightDeath[k] = nodeParentDeath[k] - nodeLeftDeath[k];
    nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k];
    if (nodeLeftDeath[k] > 0) {
      (*leftDeathTimeSize) ++;
    }
    if (nodeRightDeath[k] > 0) {
      (*rightDeathTimeSize) ++;
    }
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nRunning Split Risk Counts: \n");
    Rprintf("     timIdx    PARrisk    LFTrisk    RGTrisk    PARdeath   LFTdeath   RGTdeath\n");
    for (k=1; k <=  localDeathTimeSize; k++) {
      Rprintf(" %10d %10d %10d %10d %10d %10d %10d\n", k,
              nodeParentAtRisk[k], nodeLeftAtRisk[k], nodeRightAtRisk[k],
              nodeParentDeath[k], nodeLeftDeath[k], nodeRightDeath[k]);
    }
    uint totalDeathCount = 0;
    uint leftDeathCount  = 0;
    uint rightDeathCount = 0;
    for (k=1; k <=  localDeathTimeSize; k++) {
      totalDeathCount += nodeParentDeath[k];
      leftDeathCount  += nodeLeftDeath[k];
      rightDeathCount += nodeRightDeath[k];
    }
    Rprintf("\nRunning Split Total LFT & RGT Deaths: \n");
    Rprintf("                                             %10d %10d %10d \n", totalDeathCount, leftDeathCount, rightDeathCount);
    Rprintf("\nRunning Split Total LFT & RGT Unique Death Times: \n");
    Rprintf(" %10d %10d \n", *leftDeathTimeSize, *rightDeathTimeSize);
  }          
  if (getTraceFlag() & SPLT_MED_TRACE) {
    Rprintf("\nvirtuallySplitNode() EXIT ...\n");
  }
}
void getReweightedRandomPair (uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\ngetReweightedRandomPair() ENTRY ...\n");
  }
  if(_factorList[relativeFactorSize] == NULL) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Factor not allocated for size:  %10d", relativeFactorSize);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  randomGroupIndex = (uint) ceil(ran2(_seed2Ptr)*((_factorList[relativeFactorSize] -> cardinalGroupCount) *1.0));
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nRandomly Selected Group Index:  %10d", randomGroupIndex);
  }
  createRandomBinaryPair(relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\ngetReweightedRandomPair() EXIT ...\n");
  }
}
void getRandomPair (uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  double randomValue;
  uint k;
  if (getTraceFlag() & FACT_LOW_TRACE) {
    Rprintf("\ngetRandomPair() ENTRY ...\n");
  }
  if(_factorList[relativeFactorSize] == NULL) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Factor not allocated for size:  %10d", relativeFactorSize);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  double *cdf = dvector(1, _factorList[relativeFactorSize] -> cardinalGroupCount);
  if (relativeFactorSize <= MAX_EXACT_LEVEL) {
    for (k=1; k <= _factorList[relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = (double) ((uint*) _factorList[relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  else {
    for (k=1; k <= _factorList[relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = ((double*) _factorList[relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  if (getTraceFlag() & FACT_LOW_TRACE) {
    Rprintf("\nFactor (relativeFactorSize, cardinalGroupCount):  (%10d, %10d) \n", relativeFactorSize, _factorList[relativeFactorSize] -> cardinalGroupCount);
  }
  for (k=2; k <= _factorList[relativeFactorSize] -> cardinalGroupCount; k++) {
    cdf[k] += cdf[k-1];
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nUpdated CDF based on cardinal group size:  ");
    Rprintf("\n     index          cdf");
    for (k=1; k <= _factorList[relativeFactorSize] -> cardinalGroupCount; k++) {
      Rprintf("\n%10d  %12f", k, cdf[k]);
    }
    Rprintf("\n");
  }
  randomValue = ceil((ran2(_seed2Ptr) * cdf[_factorList[relativeFactorSize] -> cardinalGroupCount]));
  randomGroupIndex = 1;
  while (randomValue > cdf[randomGroupIndex]) {
    randomGroupIndex ++;
  }
  free_dvector(cdf, 1, _factorList[relativeFactorSize] -> cardinalGroupCount);
  if (getTraceFlag() & FACT_LOW_TRACE) {
    Rprintf("\nRandomly Selected Group Index:  %10d", randomGroupIndex);
  }
  createRandomBinaryPair(relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
  if (getTraceFlag() & FACT_LOW_TRACE) {
    Rprintf("\ngetRandomPair() EXIT ...\n");
  }
}
void createRandomBinaryPair(uint    relativeFactorSize, 
                            uint    absoluteFactorSize,
                            uint    groupIndex, 
                            double *absoluteLevel, 
                            uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint k, offset;
  if (getTraceFlag() & FACT_LOW_TRACE) {
    Rprintf("\ncreateRandomBinaryPair() ENTRY ...\n");
  }
  if (getTraceFlag() & FACT_LOW_TRACE) {
    Rprintf("\nrelativeSize absoluteSize   groupIndex ");
    Rprintf("\n%12d %12d %12d", relativeFactorSize, absoluteFactorSize, groupIndex);
    Rprintf("\n");
  }
  mwcpSizeAbsolute = _factorList[absoluteFactorSize] -> mwcpSize;
  char *localPermissible = cvector(1, relativeFactorSize);
  uint *randomLevel = uivector(1, groupIndex);
  for (k = 1; k <= relativeFactorSize; k++) {
    localPermissible[k] = TRUE;
  }
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = getSelectableElement(relativeFactorSize, localPermissible, NULL);
    localPermissible[randomLevel[k]] = FALSE;
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nAbsolute Levels:  ");
    Rprintf("\n     index      level");
    for (k=1; k <= relativeFactorSize; k++) {
      Rprintf("\n%10d  %10d", k, (uint) absoluteLevel[k]);
    }
    Rprintf("\n");
    Rprintf("\nRandomly Selected Levels Prior to Remapping:  ");
    Rprintf("\n     index      level");
    for (k=1; k <= groupIndex; k++) {
      Rprintf("\n%10d  %10d", k, randomLevel[k]);
    }
    Rprintf("\n");
  }
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = (uint) absoluteLevel[randomLevel[k]];
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nRandomly Selected Levels After Remapping:  ");
    Rprintf("\n     index      level");
    for (k=1; k <= groupIndex; k++) {
      Rprintf("\n%10d  %10d", k, randomLevel[k]);
    }
    Rprintf("\n");
  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
      pair[offset] = 0;
  }
  for (k = 1; k <= groupIndex; k++) {
    mwcpLevelIdentifier = (randomLevel[k] >> (3 + ulog2(SIZE_OF_INTEGER))) + ((randomLevel[k] & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
    if (getTraceFlag() & FACT_HGH_TRACE) {
      Rprintf("\n MWCP Level Identifier:   %10d ", mwcpLevelIdentifier);
      Rprintf("\n upower() bit:  %10d ", randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
    }
    pair[mwcpLevelIdentifier] += upower(2, randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  }
  free_cvector(localPermissible, 1, relativeFactorSize);
  free_uivector(randomLevel, 1, groupIndex);
  if (getTraceFlag() & FACT_LOW_TRACE) {
    Rprintf("\ncreateRandomBinaryPair() EXIT ...\n");
  }
}
void convertRelToAbsBinaryPair(uint    relativeFactorSize, 
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel, 
                               uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint coercedAbsoluteLevel;
  uint k, offset;
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nconvertRelToAbsBinaryPair() ENTRY ...\n");
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nrelativeSize absoluteSize      relPair");
    Rprintf("\n%12d %12d %12x", relativeFactorSize, absoluteFactorSize, relativePair);
    Rprintf("\n");
  }
  mwcpSizeAbsolute = _factorList[absoluteFactorSize] -> mwcpSize;
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nAbsolute Levels:  ");
    Rprintf("\n     index      level");
    for (k=1; k <= relativeFactorSize; k++) {
      Rprintf("\n%10d  %10d", k, (uint) absoluteLevel[k]);
    }
    Rprintf("\n");
    Rprintf("\nRelative Pair:  %8x \n", relativePair);
  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
      pair[offset] = 0;
  }
  for (k = 1; k <= relativeFactorSize; k++) {
    if (relativePair & ((uint) 0x01)) {
      coercedAbsoluteLevel = (uint) absoluteLevel[k];
      mwcpLevelIdentifier = (coercedAbsoluteLevel >> (3 + ulog2(SIZE_OF_INTEGER))) + ((coercedAbsoluteLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      pair[mwcpLevelIdentifier] += upower(2, coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
      if (getTraceFlag() & FACT_HGH_TRACE) {
        Rprintf("\n MWCP Level Identifier:   %10d ", mwcpLevelIdentifier);
        Rprintf("\n upower() bit:  %10d ", coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1);
      }
    }
    relativePair = relativePair >> 1;
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nconvertRelToAbsBinaryPair() EXIT ...\n");
  }
}
char summarizeSplitResult(uint splitParameterMax, double deltaMax) {
  char result;
  uint k;
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nsummarizeSplitResult() ENTRY ...\n");
  }
  if (splitParameterMax > 0) {
    result = TRUE;
    if (getTraceFlag() & SPLT_LOW_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm        Delta \n");
      Rprintf(" %10d %12.4f \n", splitParameterMax, deltaMax);
      if (strcmp(_xType[splitParameterMax], "C") == 0) {
        Rprintf(" at MWCPsize= %2d, mwcp= ", _splitValueMaxFactSize);
        for (k = _splitValueMaxFactSize; k >= 1; k--) {
          Rprintf("%8x ", _splitValueMaxFactPtr[k]);
        }
        Rprintf("\n");
      }
      else {
        Rprintf(" at %12.4f \n", _splitValueMaxCont);
      }
    }
  }
  else {
    result = FALSE;
  }
  if (getTraceFlag() & SPLT_HGH_TRACE) {
    Rprintf("\nsummarizeSplitResult(%1d) EXIT ...\n", result);
  }
  return result;
}
