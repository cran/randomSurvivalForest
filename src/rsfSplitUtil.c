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
  if (delta > *deltaMax) {
    *deltaMax = delta;
    *splitParameterMax = randomCovariate;
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
    }
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
      if((*permissibleSplitSize)[actualCovariateCount] >= 2) {
        randomSplitVector[candidateCovariate] = ACTIVE;
        (*covariateIndex)[actualCovariateCount] = candidateCovariate;
        actualCovariateCount ++;
      }
      else {
        (parent -> permissibleSplit)[candidateCovariate] = FALSE;
        randomSplitVector[candidateCovariate] = FALSE;
      }
    }
  }
  actualCovariateCount --;
  free_cvector(randomSplitVector, 1, _xSize);
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
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeIndex[++(*localDeathTimeSize)] = i;
      }
    }
    if ((*localDeathTimeSize) >= _minimumDeathCount) {
      result = TRUE;
    }
    else {
    }
  }
  else {
  }
  return result;
}
void stackSplit(uint **localMembershipIndex, 
                uint **localDeathTimeCount, 
                uint **localDeathTimeIndex) {
  *localMembershipIndex = uivector(1, _observationSize);
  *localDeathTimeCount = uivector(1, _masterTimeSize);
  *localDeathTimeIndex = uivector(1, _masterTimeSize);
}
void unstackSplit(uint *localMembershipIndex, 
                  uint *localDeathTimeCount, 
                  uint *localDeathTimeIndex) {
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, _masterTimeSize);
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
  *nodeParentDeath  = uivector(1, deathTimeSize);
  *nodeParentAtRisk = uivector(1, deathTimeSize);
  *nodeLeftDeath  = uivector(1, deathTimeSize);
  *nodeLeftAtRisk = uivector(1, deathTimeSize);
  *nodeRightDeath  = uivector(1, deathTimeSize);
  *nodeRightAtRisk = uivector(1, deathTimeSize);
  if ((nodeSize > 0) && (localSplitIndicator != NULL)) {
    *localSplitIndicator = cvector(1, nodeSize);
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
  free_uivector(nodeParentDeath, 1, deathTimeSize);
  free_uivector(nodeParentAtRisk, 1, deathTimeSize);
  free_uivector(nodeLeftDeath, 1, deathTimeSize);
  free_uivector(nodeLeftAtRisk, 1, deathTimeSize);
  free_uivector(nodeRightDeath, 1, deathTimeSize);
  free_uivector(nodeRightAtRisk, 1, deathTimeSize);
  if ((nodeSize > 0) && (localSplitIndicator != NULL)) {
    free_cvector(localSplitIndicator, 1, nodeSize);
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
  for (i=1; i <= localDeathTimeSize; i++) {
    nodeParentAtRisk[i] = 0;
    nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
    for (j=1; j <= localMembershipSize; j++) {
      if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
        nodeParentAtRisk[i] ++;
      }
    }
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
  splitLength = 0;  
  (*permissibleSplitPtr) = NULL;  
  if (strcmp(_xType[randomCovariateIndex], "C") == 0) {
    *factorFlag = TRUE;
    if(_factorList[permissibleSplitSize] == NULL) {
      _factorList[permissibleSplitSize] = makeFactor(permissibleSplitSize, FALSE);
    }
    factorSizeAbsolute = _factorSize[_factorMap[randomCovariateIndex]];
    *mwcpSizeAbsolute = _factorList[factorSizeAbsolute] -> mwcpSize;
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
      index = 0;  
      for (m = 1; m <= localDeathTimeSize; m++) {
        if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
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
}
void getReweightedRandomPair (uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  if(_factorList[relativeFactorSize] == NULL) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Factor not allocated for size:  %10d", relativeFactorSize);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  randomGroupIndex = (uint) ceil(ran2(_seed2Ptr)*((_factorList[relativeFactorSize] -> cardinalGroupCount) *1.0));
  createRandomBinaryPair(relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
}
void getRandomPair (uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  double randomValue;
  uint k;
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
  for (k=2; k <= _factorList[relativeFactorSize] -> cardinalGroupCount; k++) {
    cdf[k] += cdf[k-1];
  }
  randomValue = ceil((ran2(_seed2Ptr) * cdf[_factorList[relativeFactorSize] -> cardinalGroupCount]));
  randomGroupIndex = 1;
  while (randomValue > cdf[randomGroupIndex]) {
    randomGroupIndex ++;
  }
  free_dvector(cdf, 1, _factorList[relativeFactorSize] -> cardinalGroupCount);
  createRandomBinaryPair(relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
}
void createRandomBinaryPair(uint    relativeFactorSize, 
                            uint    absoluteFactorSize,
                            uint    groupIndex, 
                            double *absoluteLevel, 
                            uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint k, offset;
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
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = (uint) absoluteLevel[randomLevel[k]];
  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
      pair[offset] = 0;
  }
  for (k = 1; k <= groupIndex; k++) {
    mwcpLevelIdentifier = (randomLevel[k] >> (3 + ulog2(SIZE_OF_INTEGER))) + ((randomLevel[k] & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
    pair[mwcpLevelIdentifier] += upower(2, randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  }
  free_cvector(localPermissible, 1, relativeFactorSize);
  free_uivector(randomLevel, 1, groupIndex);
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
  mwcpSizeAbsolute = _factorList[absoluteFactorSize] -> mwcpSize;
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
      pair[offset] = 0;
  }
  for (k = 1; k <= relativeFactorSize; k++) {
    if (relativePair & ((uint) 0x01)) {
      coercedAbsoluteLevel = (uint) absoluteLevel[k];
      mwcpLevelIdentifier = (coercedAbsoluteLevel >> (3 + ulog2(SIZE_OF_INTEGER))) + ((coercedAbsoluteLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      pair[mwcpLevelIdentifier] += upower(2, coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
    }
    relativePair = relativePair >> 1;
  }
}
char summarizeSplitResult(uint splitParameterMax, double deltaMax) {
  char result;
  uint k;
  if (splitParameterMax > 0) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  return result;
}
