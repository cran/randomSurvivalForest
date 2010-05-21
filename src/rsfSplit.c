////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.3
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
#include      "rsfSplit.h"
char getBestSplit(Node *parent, 
                  uint *splitParameterMax) {
  char  result;
  switch(_splitRule) {
  case LOG_RANK:
    result = logRank(parent, splitParameterMax);
    break;
  case CONSERVE_EVENTS:
    result = conserveEvents(parent, splitParameterMax);
    break;
  case LOG_RANK_SCORE:
    result = logRankScore(parent, splitParameterMax);
    break;
  case RANDOM_SPLIT:
    result = randomSplit(parent, splitParameterMax);
    break;
  case LOG_RANK_LAU_CR:
    result = logRankLauCR(parent, splitParameterMax);
    break;
  case LOG_RANK_CR:
    result = logRankCR(parent, splitParameterMax);
    break;
  case SUB_CUM_HAZ_CR:
    result = subCumHazCR(parent, splitParameterMax);
    break;
default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid split rule:  %10d", _splitRule);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  return result;
}
char randomSplit(Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  uint saveMinimumDeathCount;
  double delta, deltaNum, deltaDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j;
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  _splitValueMaxFactSize = 0;
  _splitValueMaxFactPtr  = NULL;
  _splitValueMaxCont     = NA_REAL;
  deltaMax               = -EPSILON;
  stackSplit(& localMembershipIndex, 
             & localDeathTimeCount, 
             & localDeathTimeIndex);
  saveMinimumDeathCount = _minimumDeathCount;
  _minimumDeathCount = 1;
  result = getDeathCount(parent, 
                         localMembershipIndex, 
                         localDeathTimeCount, 
                         localDeathTimeIndex,
                         & localMembershipSize,
                         & localDeathTimeSize);
  _minimumDeathCount = saveMinimumDeathCount;
  if(result) {
    if (localDeathTimeSize >= _minimumDeathCount) {
      result = TRUE;
    }
    else {
      result = FALSE;
    }
  }
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
                      & nodeRightAtRisk,
                      0,
                      NULL);
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
                           & rightDeathTimeSize,
                           NULL);
        if ((leftDeathTimeSize  >= 1) && (rightDeathTimeSize  >= 1)) {
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
                        nodeRightAtRisk,
                        0,
                        NULL);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
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
                      & nodeRightAtRisk,
                      0,
                      NULL);
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
                           & rightDeathTimeSize,
                           NULL);
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
                        nodeRightAtRisk,
                        0,
                        NULL);
    free_dvector(predictorValue, 1, localMembershipSize);
    free_uivector(localSplitRank, 1, localMembershipSize);
    free_uivector(survivalTimeIndexRank, 1, localMembershipSize);
    free_dvector(survivalRank, 1, localMembershipSize);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
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
                      & nodeRightAtRisk,
                      0,
                      NULL);
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
                           & rightDeathTimeSize,
                           NULL);
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
                        nodeRightAtRisk,
                        0,
                        NULL);
    free_dvector(nelsonAalenLeft, 1, localDeathTimeSize);
    free_dvector(nelsonAalenRight, 1, localDeathTimeSize);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
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
                      & nodeRightAtRisk,
                      0,
                      NULL);
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
                           & rightDeathTimeSize,
                           NULL);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          delta = deltaNum = deltaDen =  0.0;
          for (k=1; k <= localDeathTimeSize; k++) {
            deltaNum = deltaNum + ((double) nodeLeftDeath[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentDeath[k]) / nodeParentAtRisk[k]));
            if (nodeParentAtRisk[k] >= 2) {
              deltaDen = deltaDen + (
                                     ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                                     (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                                     ((double) (nodeParentAtRisk[k] - nodeParentDeath[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentDeath[k]
                                     );
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
                        nodeRightAtRisk,
                        0,
                        NULL);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  return result;
}
char subCumHazCR (Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  double delta, deltaNum, deltaMax;
  double sumValueParent, sumValueLeft, sumValueRight;
  double denValueParent, denValueLeft, denValueRight;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint index;
  uint i, j, k, m, q;
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
    char *localSplitIndicator;
    stackSplitCompact(localDeathTimeSize,
                      & nodeParentDeath,
                      & nodeParentAtRisk,
                      & nodeLeftDeath,
                      & nodeLeftAtRisk,
                      & nodeRightDeath,
                      & nodeRightAtRisk,
                      localMembershipSize,
                      & localSplitIndicator);
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
                           & rightDeathTimeSize,
                           localSplitIndicator);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          uint **nodeParentEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeLeftEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeRightEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint *weight  = uivector(1, _eventTypeSize);
          for (m=1; m <= localDeathTimeSize; m++) {
            for (q = 1; q <= _eventTypeSize; q++) {
              nodeParentEvent[q][m] = nodeLeftEvent[q][m] = nodeRightEvent[q][m] = 0;
            }
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_status[localMembershipIndex[k]] > 0) {
              index = 0;  
              for (m = 1; m <= localDeathTimeSize; m++) {
                if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                  index = m;
                }
                else {
                  m = localDeathTimeSize;
                }
              }
              nodeParentEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              if (localSplitIndicator[k] == LEFT) {
                nodeLeftEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              }
              else {
                nodeRightEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              }
            }
          }
          for (q = 1; q <= _eventTypeSize; q++) {
            weight[q] = 1.0;
          }
          double **subDensityParent = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double **subDensityLeft   = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double **subDensityRight  = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double *survivalParent = dvector(1, localDeathTimeSize);
          double *survivalLeft   = dvector(1, localDeathTimeSize);
          double *survivalRight  = dvector(1, localDeathTimeSize);
          double **subDistributionParent = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double **subDistributionLeft   = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double **subDistributionRight  = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double **subDistributionCHFParent = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double **subDistributionCHFLeft   = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double **subDistributionCHFRight  = dmatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          double estimate;
          for (q=1; q <= localDeathTimeSize; q++) {
            estimate = 1.0;
            if (nodeParentDeath[q] > 0) {
              if (nodeParentAtRisk[q] >= 1) {
                estimate = 1.0 - ((double) nodeParentDeath[q] / nodeParentAtRisk[q]);
              }
            }
            survivalParent[q] = estimate;
            estimate = 1.0;
            if (nodeLeftDeath[q] > 0) {
              if (nodeLeftAtRisk[q] >= 1) {
                estimate = 1.0 - ((double) nodeLeftDeath[q] / nodeLeftAtRisk[q]);
              }
            }
            survivalLeft[q] = estimate;
            estimate = 1.0;
            if (nodeRightDeath[q] > 0) {
              if (nodeRightAtRisk[q] >= 1) {
                estimate = 1.0 - ((double) nodeRightDeath[q] / nodeRightAtRisk[q]);
              }
            }
            survivalRight[q] = estimate;
          }  
          for (q=2; q <= localDeathTimeSize; q++) {
            survivalParent[q] *= survivalParent[q-1];
            survivalLeft[q] *= survivalLeft[q-1];
            survivalRight[q] *= survivalRight[q-1];
          }
          for (q=1; q <= _eventTypeSize; q++) {      
            subDensityParent[q][1] = 0.0;
            if (nodeParentEvent[q][1] > 0) {
              if (nodeParentAtRisk[1] >= 1) {
                subDensityParent[q][1] = 1.0 * ((double) nodeParentEvent[q][1] / nodeParentAtRisk[1]);
              }
            }
            subDensityLeft[q][1] = 0.0;
            if (nodeLeftEvent[q][1] > 0) {
              if (nodeLeftAtRisk[1] >= 1) {
                subDensityLeft[q][1] = 1.0 * ((double) nodeLeftEvent[q][1] / nodeLeftAtRisk[1]);
              }
            }
            subDensityRight[q][1] = 0.0;
            if (nodeRightEvent[q][1] > 0) {
              if (nodeRightAtRisk[1] >= 1) {
                subDensityRight[q][1] = 1.0 * ((double) nodeRightEvent[q][1] / nodeRightAtRisk[1]);
              }
            }
            for (m=2; m <= localDeathTimeSize; m++) {
              subDensityParent[q][m] = 0.0;
              if (nodeParentEvent[q][m] > 0) {
                if (nodeParentAtRisk[m] >= 1) {
                  subDensityParent[q][m] = survivalParent[m-1] * ((double) nodeParentEvent[q][m] / nodeParentAtRisk[m]);
                }
              }
            }
            for (m=2; m <= localDeathTimeSize; m++) {
              subDensityLeft[q][m] = 0.0;
              if (nodeLeftEvent[q][m] > 0) {
                if (nodeLeftAtRisk[m] >= 1) {
                  subDensityLeft[q][m] = survivalLeft[m-1] * ((double) nodeLeftEvent[q][m] / nodeLeftAtRisk[m]);
                }
              }
            }
            for (m=2; m <= localDeathTimeSize; m++) {
              subDensityRight[q][m] = 0.0;
              if (nodeRightEvent[q][m] > 0) {
                if (nodeRightAtRisk[m] >= 1) {
                  subDensityRight[q][m] = survivalRight[m-1] * ((double) nodeRightEvent[q][m] / nodeRightAtRisk[m]);
                }
              }
            }
          }
          for (q=1; q <= _eventTypeSize; q++) {
            subDistributionParent[q][1] = subDensityParent[q][1];
            subDistributionLeft[q][1] = subDensityLeft[q][1];
            subDistributionRight[q][1] = subDensityRight[q][1];
            sumValueParent = sumValueLeft = sumValueRight = 0.0;
            for (m=2; m <= localDeathTimeSize; m++) {
              sumValueParent += subDensityParent[q][m-1];
              denValueParent = 1.0 - sumValueParent;
              sumValueLeft += subDensityLeft[q][m-1];
              denValueLeft = 1.0 - sumValueLeft;
              sumValueRight += subDensityRight[q][m-1];
              denValueRight = 1.0 - sumValueRight;
              if (fabs(denValueParent) < EPSILON) {
                subDistributionParent[q][m] = 0.0;
              }
              else {
                subDistributionParent[q][m] = subDensityParent[q][m] / denValueParent;
              }
              if (fabs(denValueLeft) < EPSILON) {
                subDistributionLeft[q][m] = 0.0;
              }
              else {
                subDistributionLeft[q][m] = subDensityLeft[q][m] / denValueLeft;
              }
              if (fabs(denValueRight) < EPSILON) {
                subDistributionRight[q][m] = 0.0;
              }
              else {
                subDistributionRight[q][m] = subDensityRight[q][m] / denValueRight;
              }
            }
          }
          for (q=1; q <= _eventTypeSize; q++) {      
            subDistributionCHFParent[q][1] = subDistributionParent[q][1];
            subDistributionCHFLeft[q][1] = subDistributionLeft[q][1];
            subDistributionCHFRight[q][1] = subDistributionRight[q][1];
            for (m=2; m <= localDeathTimeSize; m++) {
              subDistributionCHFParent[q][m] += subDistributionCHFParent[q][m-1] + subDistributionParent[q][m];
              subDistributionCHFLeft[q][m] += subDistributionCHFLeft[q][m-1] + subDistributionLeft[q][m];
              subDistributionCHFRight[q][m] += subDistributionCHFRight[q][m-1] + subDistributionRight[q][m];
            }
          }
          delta = deltaNum = 0.0;
          for (q = 1; q <= _eventTypeSize; q++) {
            deltaNum = 0;
            for (m=1; m <= localDeathTimeSize; m++) {
              deltaNum += pow(subDistributionCHFLeft[q][m] - subDistributionCHFParent[q][m], 2.0);
            }
            delta += weight[q] * deltaNum;
          }
          updateMaximumSplit(delta,
                             randomCovariateIndex[i],
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             & deltaMax,
                             splitParameterMax,
                             permissibleSplitPtr);
          free_dmatrix(subDensityParent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dmatrix(subDensityLeft, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dmatrix(subDensityRight, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dvector(survivalParent, 1, localDeathTimeSize);
          free_dvector(survivalLeft, 1, localDeathTimeSize);
          free_dvector(survivalRight, 1, localDeathTimeSize);
          free_dmatrix(subDistributionParent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dmatrix(subDistributionLeft, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dmatrix(subDistributionRight, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dmatrix(subDistributionCHFParent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dmatrix(subDistributionCHFLeft, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_dmatrix(subDistributionCHFRight, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeParentEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeLeftEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeRightEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uivector(weight, 1, _eventTypeSize);
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
                        nodeRightAtRisk,
                        localMembershipSize,
                        localSplitIndicator);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  return result;
}
char logRankCR (Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  double delta, deltaNum, deltaSubNum, deltaDen, deltaSubDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint index;
  uint i, j, k, m, q;
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
    char *localSplitIndicator;
    stackSplitCompact(localDeathTimeSize,
                      & nodeParentDeath,
                      & nodeParentAtRisk,
                      & nodeLeftDeath,
                      & nodeLeftAtRisk,
                      & nodeRightDeath,
                      & nodeRightAtRisk,
                      localMembershipSize,
                      & localSplitIndicator);
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
                           & rightDeathTimeSize,
                           localSplitIndicator);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          uint **nodeParentEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeLeftEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeRightEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint *weight  = uivector(1, _eventTypeSize);
          for (m=1; m <= localDeathTimeSize; m++) {
            for (q = 1; q <= _eventTypeSize; q++) {
              nodeParentEvent[q][m] = nodeLeftEvent[q][m] = nodeRightEvent[q][m] = 0;
            }
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_status[localMembershipIndex[k]] > 0) {
              index = 0;  
              for (m = 1; m <= localDeathTimeSize; m++) {
                if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                  index = m;
                }
                else {
                  m = localDeathTimeSize;
                }
              }
              nodeParentEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              if (localSplitIndicator[k] == LEFT) {
                nodeLeftEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              }
              else {
                nodeRightEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              }
            }
          }
          for (q = 1; q <= _eventTypeSize; q++) {
            weight[q] = 1.0;
          }
          delta = deltaNum = deltaDen =  0.0;
          for (q = 1; q <= _eventTypeSize; q++) {
            deltaSubNum = 0;
            for (m=1; m <= localDeathTimeSize; m++) {
              deltaSubNum = deltaSubNum + (nodeLeftEvent[q][m] - (nodeParentEvent[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])));
            }
            deltaNum = deltaNum + (weight[q] * deltaSubNum);
          }
          for (q = 1; q <= _eventTypeSize; q++) {
            deltaSubDen = 0.0;
            for (m=1; m <= localDeathTimeSize; m++) {
              if (nodeParentAtRisk[m] >= 2) {
                deltaSubDen = deltaSubDen  + (
                                              (nodeParentEvent[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                              (1.0 - ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                              ((double) (nodeParentAtRisk[m] - nodeParentEvent[q][m]) / (nodeParentAtRisk[m] - 1))
                                              );
              }
            }
            deltaDen += (weight[q] * weight[q] * deltaSubDen);
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
          free_uimatrix(nodeParentEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeLeftEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeRightEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uivector(weight, 1, _eventTypeSize);
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
                        nodeRightAtRisk,
                        localMembershipSize,
                        localSplitIndicator);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  return result;
}
char logRankLauCR (Node *parent, uint *splitParameterMax) {
  uint *localMembershipIndex;
  uint *localDeathTimeCount;
  uint *localDeathTimeIndex;
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  double delta, deltaNum, deltaSubNum, deltaDen, deltaSubDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint index;
  uint i, j, k, m, q, r, s;
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
    char *localSplitIndicator;
    stackSplitCompact(localDeathTimeSize,
                      & nodeParentDeath,
                      & nodeParentAtRisk,
                      & nodeLeftDeath,
                      & nodeLeftAtRisk,
                      & nodeRightDeath,
                      & nodeRightAtRisk,
                      localMembershipSize,
                      & localSplitIndicator);
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
                           & rightDeathTimeSize,
                           localSplitIndicator);
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          uint **nodeParentEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeLeftEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeRightEvent = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint *weight  = uivector(1, _eventTypeSize);
          uint **nodeParentInclusiveAtRisk = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeLeftInclusiveAtRisk = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          uint **nodeRightInclusiveAtRisk = uimatrix(1, _eventTypeSize, 1, localDeathTimeSize);
          for (m=1; m <= localDeathTimeSize; m++) {
            for (q = 1; q <= _eventTypeSize; q++) {
              nodeParentEvent[q][m] = nodeLeftEvent[q][m] = nodeRightEvent[q][m] = 0;
            }
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_status[localMembershipIndex[k]] > 0) {
              index = 0;  
              for (m = 1; m <= localDeathTimeSize; m++) {
                if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                  index = m;
                }
                else {
                  m = localDeathTimeSize;
                }
              }
              nodeParentEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              if (localSplitIndicator[k] == LEFT) {
                nodeLeftEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              }
              else {
                nodeRightEvent[_eventTypeIndex[(uint) _status[localMembershipIndex[k]]]][index] ++;
              }
            }
          }
          for (m=1; m <= localDeathTimeSize; m++) {
            for (q = 1; q <= _eventTypeSize; q++) {
              nodeParentInclusiveAtRisk[q][m] = nodeParentAtRisk[m];
              nodeLeftInclusiveAtRisk[q][m]   = nodeLeftAtRisk[m];
              nodeRightInclusiveAtRisk[q][m]  = nodeRightAtRisk[m];
              for (s = 1; s < m; s++) {
                for (r = 1; r <= _eventTypeSize; r++) {
                  if (q != r) {
                    nodeParentInclusiveAtRisk[q][m]  += nodeParentEvent[r][s];
                    nodeLeftInclusiveAtRisk[q][m]    += nodeLeftEvent[r][s];
                    nodeRightInclusiveAtRisk[q][m]   += nodeRightEvent[r][s];
                  }
                }
              }
            }
          }
          for (q = 1; q <= _eventTypeSize; q++) {
            weight[q] = 1.0;
          }
          delta = deltaNum = deltaDen =  0.0;
          for (q = 1; q <= _eventTypeSize; q++) {
            deltaSubNum = 0;
            for (m=1; m <= localDeathTimeSize; m++) {
              deltaSubNum = deltaSubNum + (nodeLeftEvent[q][m] - (nodeParentEvent[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])));
            }
            deltaNum = deltaNum + (weight[q] * deltaSubNum);
          }
          for (q = 1; q <= _eventTypeSize; q++) {
            deltaSubDen = 0;
            for (m=1; m <= localDeathTimeSize; m++) {
              if (nodeParentAtRisk[m] >= 2) {
                deltaSubDen = deltaSubDen  + (
                                              (nodeParentEvent[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                              (1.0 - ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                              ((double) (nodeParentInclusiveAtRisk[q][m] - nodeParentEvent[q][m]) / (nodeParentInclusiveAtRisk[q][m] - 1))
                                              );
              }
            }
            deltaDen = deltaDen + (weight[q] * weight[q] * deltaSubDen);
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
          free_uimatrix(nodeParentInclusiveAtRisk, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeLeftInclusiveAtRisk, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeRightInclusiveAtRisk, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeParentEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeLeftEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uimatrix(nodeRightEvent, 1, _eventTypeSize, 1, localDeathTimeSize);
          free_uivector(weight, 1, _eventTypeSize);
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
                        nodeRightAtRisk,
                        localMembershipSize,
                        localSplitIndicator);
  }  
  unstackSplit(localMembershipIndex, 
               localDeathTimeCount, 
               localDeathTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, deltaMax);
  return result;
}
