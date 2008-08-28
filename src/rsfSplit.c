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

#include       "global.h"
#include       "nrutil.h"
#include "rsfFactorOps.h"
#include     "rsfSplit.h"
#include "rsfSplitUtil.h"
extern uint getTraceFlag();
char randomSplit(Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  double delta, deltaNum, deltaDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j;
  if (getTraceFlag() & SPLT_MED_TRACE) {
    Rprintf("\nrandomSplit() ENTRY ...\n");
  }
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  _splitValueMaxFactSize = 0;
  _splitValueMaxFactPtr  = NULL;
  _splitValueMaxCont     = NA_REAL;
  deltaMax               = -EPSILON;
  stackSplit(& localMembershipIndex, 
             & localDeathTimeCount, 
             & localDeathTimeIndex);
  result = getDeathCount(parent, 
                         localMembershipIndex, 
                         localDeathTimeCount, 
                         localDeathTimeIndex,
                         & localMembershipSize,
                         & localDeathTimeSize);
  if (result) {
    char *covariateStatus = NULL;  
    uint *nodeParentDeath, *nodeLeftDeath, *nodeRightDeath;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    stackSplitCompact(localDeathTimeSize,
                      & nodeParentDeath,
                      & nodeParentAtRisk,
                      & nodeLeftDeath,
                      & nodeLeftAtRisk,
                      & nodeRightDeath,
                      & nodeRightAtRisk);
    getAtRisk(localMembershipIndex,
              localDeathTimeCount,
              localDeathTimeIndex,
              localMembershipSize,
              localDeathTimeSize,
              nodeParentDeath,
              nodeParentAtRisk);
    uint actualCovariateCount = stackAndSelectRandomCovariates(parent, 
                                                               localMembershipSize,
                                                               localMembershipIndex,
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize);
    if (actualCovariateCount > 0) {
      covariateStatus = cvector(1, actualCovariateCount);
    }
    for (i = 1; i <= actualCovariateCount; i++) {
      covariateStatus[i] = TRUE;
    }
    i = getSelectableElement(actualCovariateCount, covariateStatus, NULL);
    while ((i != 0) && ((*splitParameterMax) == 0)) {
      splitLength = stackAndConstructSplitVector(localMembershipSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      for (j = 1; j < splitLength; j++) {
        virtuallySplitNode(localMembershipSize,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           localMembershipIndex,
                           permissibleSplitPtr,
                           j,
                           localDeathTimeSize,
                           localDeathTimeIndex,
                           nodeParentAtRisk,
                           nodeParentDeath,
                           nodeLeftAtRisk,
                           nodeLeftDeath,
                           & leftDeathTimeSize,
                           nodeRightAtRisk,
                           nodeRightDeath,
                           & rightDeathTimeSize);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          delta = deltaNum = deltaDen =  0.0;
          updateMaximumSplit(delta,
                             randomCovariateIndex[i],
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             & deltaMax,
                             splitParameterMax,
                             permissibleSplitPtr);
        }  
      }  
      unstackSplitVector(permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         permissibleSplitPtr);
      if(*splitParameterMax == 0) {
        covariateStatus[i] = FALSE;
        i = getSelectableElement(actualCovariateCount, covariateStatus, NULL);
      }
    }  
    unstackRandomCovariates(localMembershipSize, 
                            randomCovariateIndex, 
                            permissibleSplit, 
                            permissibleSplitSize);
    if (actualCovariateCount > 0) {
      free_cvector(covariateStatus, 1, actualCovariateCount);
    }
    unstackSplitCompact(localDeathTimeSize,
                        nodeParentDeath,
                        nodeParentAtRisk,
                        nodeLeftDeath,
                        nodeLeftAtRisk,
                        nodeRightDeath,
                        nodeRightAtRisk);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  if (getTraceFlag() & SPLT_LOW_TRACE) {
    Rprintf("\nrandomSplit(%1d) EXIT ...\n", result);
  }
  return result;
}
char logRankScore(Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  double delta, deltaNum, deltaDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  if (getTraceFlag() & SPLT_LOW_TRACE) {
    Rprintf("\nlogrankScore() ENTRY ...\n");
  }
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  _splitValueMaxFactSize = 0;
  _splitValueMaxFactPtr  = NULL;
  _splitValueMaxCont     = NA_REAL;
  deltaMax               = -EPSILON;
  stackSplit(& localMembershipIndex, 
             & localDeathTimeCount, 
             & localDeathTimeIndex);
  result = getDeathCount(parent, 
                         localMembershipIndex, 
                         localDeathTimeCount, 
                         localDeathTimeIndex,
                         & localMembershipSize,
                         & localDeathTimeSize);
  if (result) {
    double meanSurvRank, varSurvRank;
    uint leftMembershipSize;
    char daughterFlag;
    double *predictorValue = dvector(1, localMembershipSize);
    uint   *localSplitRank  = uivector(1, localMembershipSize);
    uint *survivalTimeIndexRank = uivector(1, localMembershipSize);
    double *survivalRank = dvector(1, localMembershipSize);
    uint *nodeParentDeath, *nodeLeftDeath, *nodeRightDeath;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    stackSplitCompact(localDeathTimeSize,
                      & nodeParentDeath,
                      & nodeParentAtRisk,
                      & nodeLeftDeath,
                      & nodeLeftAtRisk,
                      & nodeRightDeath,
                      & nodeRightAtRisk);
    getAtRisk(localMembershipIndex,
              localDeathTimeCount,
              localDeathTimeIndex,
              localMembershipSize,
              localDeathTimeSize,
              nodeParentDeath,
              nodeParentAtRisk);
    uint actualCovariateCount = stackAndSelectRandomCovariates(parent, 
                                                               localMembershipSize,
                                                               localMembershipIndex,
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize);
    for (i = 1; i <= actualCovariateCount; i++) {
      if (getTraceFlag() & SPLT_MED_TRACE) {
        Rprintf("\nSplitting on (index, parameter, size):  %10d %10d %10d \n", i, randomCovariateIndex[i], permissibleSplitSize[i]);
      }
      splitLength = stackAndConstructSplitVector(localMembershipSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      for (k=1; k <= localMembershipSize; k++) {
        predictorValue[k] = _observation[randomCovariateIndex[i]][localMembershipIndex[k]];
      }
      indexx(localMembershipSize, predictorValue, localSplitRank);
      for (k = 1; k <= localMembershipSize; k++) {
        survivalTimeIndexRank[k] = 0;
        for (j = 1; j <= localMembershipSize; j++) {
          if ( _masterTimeIndex[localMembershipIndex[j]]  <= _masterTimeIndex[localMembershipIndex[localSplitRank[k]]] ) {
            survivalTimeIndexRank[k] ++;
          }
        }
      }
      if (getTraceFlag() & SPLT_MED_TRACE) {
        Rprintf("\nLocal Membership Information for Parent Node: \n");
        Rprintf("       RANK    INDVidx       SPLTval    TIMEidx -> SORTidx \n");
        for (k=1; k <=  localMembershipSize; k++) {
          Rprintf(" %10d %10d %10.4f %10d %10d\n", k,
                  localMembershipIndex[localSplitRank[k]], _observation[i][localMembershipIndex[localSplitRank[k]]],
                  _masterTimeIndex[localMembershipIndex[localSplitRank[k]]], survivalTimeIndexRank[k]);
        }
      }
      meanSurvRank = varSurvRank = 0;
      for (k = 1; k <= localMembershipSize; k++) {
        survivalRank[k] = 0;
        for (j = 1; j <= survivalTimeIndexRank[k]; j++) {
          survivalRank[k] = survivalRank[k] + (_status[localMembershipIndex[localSplitRank[j]]] / (localMembershipSize - survivalTimeIndexRank[j] + 1) );
        }
        survivalRank[k] = _status[localMembershipIndex[localSplitRank[k]]] - survivalRank[k];
        meanSurvRank = meanSurvRank + survivalRank[k];
        varSurvRank = varSurvRank +  pow(survivalRank[k], 2.0);
      }
      varSurvRank = ( varSurvRank - (pow(meanSurvRank, 2.0) / localMembershipSize) ) / (localMembershipSize - 1);
      meanSurvRank = meanSurvRank / localMembershipSize;
      for (j = 1; j < splitLength; j++) {
        virtuallySplitNode(localMembershipSize,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           localMembershipIndex,
                           permissibleSplitPtr,
                           j,
                           localDeathTimeSize,
                           localDeathTimeIndex,
                           nodeParentAtRisk,
                           nodeParentDeath,
                           nodeLeftAtRisk,
                           nodeLeftDeath,
                           & leftDeathTimeSize,
                           nodeRightAtRisk,
                           nodeRightDeath,
                           & rightDeathTimeSize);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          delta = deltaNum = deltaDen =  0.0;
          leftMembershipSize = 0;
          for (k=1; k <= localMembershipSize; k++) {
            daughterFlag = RIGHT;
            if (factorFlag == TRUE) {
              daughterFlag = splitOnFactor((uint) _observation[randomCovariateIndex[i]][localMembershipIndex[localSplitRank[k]]], (uint*) permissibleSplitPtr + ((j - 1) * mwcpSizeAbsolute));
            }
            else {
              if (_observation[randomCovariateIndex[i]][localMembershipIndex[localSplitRank[k]]] <= ((double*) permissibleSplitPtr)[j]) {
                daughterFlag = LEFT;
              }
            }
            if (daughterFlag == LEFT) {
              leftMembershipSize ++;
              deltaNum = deltaNum + survivalRank[k];
            }
          }  
          deltaNum  = deltaNum - (leftMembershipSize * meanSurvRank);
          deltaDen = leftMembershipSize * (1.0 - (leftMembershipSize / localMembershipSize)) * varSurvRank;
          deltaNum = fabs(deltaNum);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNum <= EPSILON) {
              delta = 0.0;
            }
          }
          else {
            delta = deltaNum / deltaDen;
          }
          updateMaximumSplit(delta,
                             randomCovariateIndex[i],
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             & deltaMax,
                             splitParameterMax,
                             permissibleSplitPtr);
        }  
      }  
      unstackSplitVector(permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         permissibleSplitPtr);
    }  
    unstackRandomCovariates(localMembershipSize, 
                            randomCovariateIndex, 
                            permissibleSplit, 
                            permissibleSplitSize);
    unstackSplitCompact(localDeathTimeSize,
                        nodeParentDeath,
                        nodeParentAtRisk,
                        nodeLeftDeath,
                        nodeLeftAtRisk,
                        nodeRightDeath,
                        nodeRightAtRisk);
    free_dvector(predictorValue, 1, localMembershipSize);
    free_uivector(localSplitRank, 1, localMembershipSize);
    free_uivector(survivalTimeIndexRank, 1, localMembershipSize);
    free_dvector(survivalRank, 1, localMembershipSize);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  if (getTraceFlag() & SPLT_LOW_TRACE) {
    Rprintf("\nlogRankScore(%1d) EXIT ...\n", result);
  }
  return result;
}
char conserveEvents(Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  double delta, deltaNum, deltaDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  if (getTraceFlag() & SPLT_LOW_TRACE) {
    Rprintf("\nconserveEvents() ENTRY ...\n");
  }
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  _splitValueMaxFactSize = 0;
  _splitValueMaxFactPtr  = NULL;
  _splitValueMaxCont     = NA_REAL;
  deltaMax               = -EPSILON;
  stackSplit(& localMembershipIndex, 
             & localDeathTimeCount, 
             & localDeathTimeIndex);
  result = getDeathCount(parent, 
                         localMembershipIndex, 
                         localDeathTimeCount, 
                         localDeathTimeIndex,
                         & localMembershipSize,
                         & localDeathTimeSize);
  if (result) {
    double nelsonAalenSumLeft, nelsonAalenSumRight;
    double *nelsonAalenLeft  = dvector(1, localDeathTimeSize);
    double *nelsonAalenRight = dvector(1, localDeathTimeSize);
    uint *nodeParentDeath, *nodeLeftDeath, *nodeRightDeath;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    stackSplitCompact(localDeathTimeSize,
                      & nodeParentDeath,
                      & nodeParentAtRisk,
                      & nodeLeftDeath,
                      & nodeLeftAtRisk,
                      & nodeRightDeath,
                      & nodeRightAtRisk);
    getAtRisk(localMembershipIndex,
              localDeathTimeCount,
              localDeathTimeIndex,
              localMembershipSize,
              localDeathTimeSize,
              nodeParentDeath,
              nodeParentAtRisk);
    uint actualCovariateCount = stackAndSelectRandomCovariates(parent, 
                                                               localMembershipSize,
                                                               localMembershipIndex,
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize);
    for (i = 1; i <= actualCovariateCount; i++) {
      splitLength = stackAndConstructSplitVector(localMembershipSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      for (j = 1; j < splitLength; j++) {
        virtuallySplitNode(localMembershipSize,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           localMembershipIndex,
                           permissibleSplitPtr,
                           j,
                           localDeathTimeSize,
                           localDeathTimeIndex,
                           nodeParentAtRisk,
                           nodeParentDeath,
                           nodeLeftAtRisk,
                           nodeLeftDeath,
                           & leftDeathTimeSize,
                           nodeRightAtRisk,
                           nodeRightDeath,
                           & rightDeathTimeSize);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          delta = deltaNum = deltaDen =  0.0;
          nelsonAalenSumLeft = nelsonAalenSumRight = 0.0;
          for (k=1; k <= localDeathTimeSize; k++) {
            if (nodeLeftAtRisk[k] != 0) {
              nelsonAalenLeft[k] = (double) nodeLeftDeath[k]/nodeLeftAtRisk[k];
            }
            else {
              nelsonAalenLeft[k] = 0;
            }
            if (nodeRightAtRisk[k] != 0) {
              nelsonAalenRight[k] = (double) nodeRightDeath[k]/nodeRightAtRisk[k];
            }
            else {
              nelsonAalenRight[k] = 0;
            }
          }
          for (k=2; k <= localDeathTimeSize; k++) {
            nelsonAalenLeft[k] += nelsonAalenLeft[k-1];
            nelsonAalenRight[k] += nelsonAalenRight[k-1];
          }
          for (k=1; k <= localDeathTimeSize-1; k++) {
            nelsonAalenSumLeft += (nodeLeftAtRisk[k] - nodeLeftAtRisk[k+1]) * nodeLeftAtRisk[k+1] * nelsonAalenLeft[k];
            nelsonAalenSumRight += (nodeRightAtRisk[k] - nodeRightAtRisk[k+1]) * nodeRightAtRisk[k+1] * nelsonAalenRight[k];
          }
          delta = ((nodeLeftAtRisk[1] * nelsonAalenSumLeft) + (nodeRightAtRisk[1] * nelsonAalenSumRight)) / (nodeLeftAtRisk[1] + nodeRightAtRisk[1]);
          delta = 1.0 / (1.0 + delta);
          updateMaximumSplit(delta,
                             randomCovariateIndex[i],
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             & deltaMax,
                             splitParameterMax,
                             permissibleSplitPtr);
        }  
      }  
      unstackSplitVector(permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         permissibleSplitPtr);
    }  
    unstackRandomCovariates(localMembershipSize, 
                            randomCovariateIndex, 
                            permissibleSplit, 
                            permissibleSplitSize);
    unstackSplitCompact(localDeathTimeSize,
                        nodeParentDeath,
                        nodeParentAtRisk,
                        nodeLeftDeath,
                        nodeLeftAtRisk,
                        nodeRightDeath,
                        nodeRightAtRisk);
    free_dvector(nelsonAalenLeft, 1, localDeathTimeSize);
    free_dvector(nelsonAalenRight, 1, localDeathTimeSize);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  if (getTraceFlag() & SPLT_LOW_TRACE) {
    Rprintf("\nconserveEvents(%1d) EXIT ...\n", result);
  }
  return result;
}
char logRank (Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  double delta, deltaNum, deltaDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  if (getTraceFlag() & SPLT_LOW_TRACE) {
    Rprintf("\nlogRank() ENTRY ...\n");
  }
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  _splitValueMaxFactSize = 0;
  _splitValueMaxFactPtr  = NULL;
  _splitValueMaxCont     = NA_REAL;
  deltaMax               = -EPSILON;
  stackSplit(& localMembershipIndex, 
             & localDeathTimeCount, 
             & localDeathTimeIndex);
  result = getDeathCount(parent, 
                         localMembershipIndex, 
                         localDeathTimeCount, 
                         localDeathTimeIndex,
                         & localMembershipSize,
                         & localDeathTimeSize);
  if (result) {
    uint *nodeParentDeath, *nodeLeftDeath, *nodeRightDeath;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    stackSplitCompact(localDeathTimeSize,
                      & nodeParentDeath,
                      & nodeParentAtRisk,
                      & nodeLeftDeath,
                      & nodeLeftAtRisk,
                      & nodeRightDeath,
                      & nodeRightAtRisk);
    getAtRisk(localMembershipIndex,
              localDeathTimeCount,
              localDeathTimeIndex,
              localMembershipSize,
              localDeathTimeSize,
              nodeParentDeath,
              nodeParentAtRisk);
    uint actualCovariateCount = stackAndSelectRandomCovariates(parent, 
                                                               localMembershipSize,
                                                               localMembershipIndex,
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize);
    for (i = 1; i <= actualCovariateCount; i++) {
      splitLength = stackAndConstructSplitVector(localMembershipSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      for (j = 1; j < splitLength; j++) {
        virtuallySplitNode(localMembershipSize,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           localMembershipIndex,
                           permissibleSplitPtr,
                           j,
                           localDeathTimeSize,
                           localDeathTimeIndex,
                           nodeParentAtRisk,
                           nodeParentDeath,
                           nodeLeftAtRisk,
                           nodeLeftDeath,
                           & leftDeathTimeSize,
                           nodeRightAtRisk,
                           nodeRightDeath,
                           & rightDeathTimeSize);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          delta = deltaNum = deltaDen =  0.0;
          for (k=1; k <= localDeathTimeSize; k++) {
            deltaNum = deltaNum + ((double) nodeLeftDeath[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentDeath[k]) / nodeParentAtRisk[k]));
            if (getTraceFlag() & TURN_OFF_TRACE) {
              Rprintf("\nPartial Sum deltaDen:  %10d %10.4f", k, deltaDen);
            }
            if (nodeParentAtRisk[k] >= 2) {
              deltaDen = deltaDen + (
                                     ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                                     (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                                     ((double) (nodeParentAtRisk[k] - nodeParentDeath[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentDeath[k]
                                     );
              if (getTraceFlag() & TURN_OFF_TRACE) {
                Rprintf("\nPartial Sum deltaDen:  %10d %10.4f", k, deltaDen);
              }
            }
          }
          deltaNum = fabs(deltaNum);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNum <= EPSILON) {
              delta = 0.0;
            }
          }
          else {
            delta = deltaNum / deltaDen;
          }
          updateMaximumSplit(delta,
                             randomCovariateIndex[i],
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             & deltaMax,
                             splitParameterMax,
                             permissibleSplitPtr);
        }  
      }  
      unstackSplitVector(permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         permissibleSplitPtr);
    }  
    unstackRandomCovariates(localMembershipSize, 
                            randomCovariateIndex, 
                            permissibleSplit, 
                            permissibleSplitSize);
    unstackSplitCompact(localDeathTimeSize,
                        nodeParentDeath,
                        nodeParentAtRisk,
                        nodeLeftDeath,
                        nodeLeftAtRisk,
                        nodeRightDeath,
                        nodeRightAtRisk);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  if (getTraceFlag() & SPLT_LOW_TRACE) {
    Rprintf("\nlogRank(%1d) EXIT ...\n", result);
  }
  return result;
}
