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
#include "rsfSplit.h"
extern uint getTraceFlag();
char logRankRandom (Node    *parent,
                    uint    *splitParameterMax,
                    uint    *splitValueMax,
                    double **masterSplit) {
  double delta, deltaNum, deltaDen;
  double deltaMax;
  uint i,j,k,m, index;
  uint actualCovariateCount;
  uint parentDeathCount;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nlogRankRandom() ENTRY ...\n");
  }
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeLeftDeath    = uivector(1, _masterTimeSize);
  uint *nodeRightDeath   = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  uint *nodeLeftAtRisk   = uivector(1, _masterTimeSize);
  uint *nodeRightAtRisk  = uivector(1, _masterTimeSize);
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, _masterTimeSize);
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = deltaMax = -EPSILON;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for Parent Node: \n");
    for (i=1; i <= localMembershipSize; i++) {
      Rprintf("%10d %10d \n", i, localMembershipIndex[i]);
    }
    Rprintf("\nLocal Death Time Counts:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10d \n", i, localDeathTimeCount[i]);
    }
  }
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeIndex[++localDeathTimeSize] = i;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nLocal Death Times (i, _masterTimeIndex): \n");
      for (i=1; i <= localDeathTimeSize; i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if (localDeathTimeSize >= _minimumDeathCount) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal Death Counts (i, deaths): \n");
      }
      for (i=1; i <= localDeathTimeSize; i++) {
        nodeParentAtRisk[i] = 0;
        nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
        if (getTraceFlag() & DL2_TRACE) {
          Rprintf("%10d %10d \n", i, nodeParentDeath[i]);
        }
        for (j=1; j <= localMembershipSize; j++) {
          if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
            nodeParentAtRisk[i] ++;
          }
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal At Risk Counts (i, at risk): \n");
        for (j=1; j <= localDeathTimeSize; j++) {
          Rprintf("%10d %10d \n", j, nodeParentAtRisk[j]);
        }
      }
      actualCovariateCount = selectRandomCovariates(parent, randomCovariateIndex);
      for (i=1; i <= actualCovariateCount; i++) {
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
                  );
        }
        k = (parent -> permissibleSplit)[randomCovariateIndex[i]][2] - (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
        j = floor(ran2(_seed2Ptr)*(k*1.0));
        j = j + (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
        }
        deltaNum = deltaDen =  0.0;
        leftDeathTimeSize = rightDeathTimeSize = 0;
        for (k=1; k <= localDeathTimeSize; k++) {
          nodeLeftDeath[k] = nodeLeftAtRisk[k] = 0;
        }
        for (k=1; k <= localMembershipSize; k++) {
          if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] <= masterSplit[randomCovariateIndex[i]][j]) {
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[k]);
            }
            index = 0;  
            for (m = 1; m <= localDeathTimeSize; m++) {
              if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                if (getTraceFlag() & DL3_TRACE) {
                  Rprintf("NLAR:  %10d %10d %10d \n", m, localDeathTimeIndex[m], _masterTimeIndex[localMembershipIndex[k]]);
                }
                nodeLeftAtRisk[m] ++;
                index = m;
              }
              else {
                m = localDeathTimeSize;
              }
            }
            if (_status[localMembershipIndex[k]] == 1) {
              nodeLeftDeath[index] ++;
            }
          }
          else {
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\nMember of Right Daughter (index):  %10d %10d \n", k, localMembershipIndex[k]);
            }
          }
        }  
        for (k=1; k <= localDeathTimeSize; k++) {
          nodeRightDeath[k] = nodeParentDeath[k] - nodeLeftDeath[k];
          nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k];
          if (nodeLeftDeath[k] > 0) {
            leftDeathTimeSize ++;
          }
          if (nodeRightDeath[k] > 0) {
            rightDeathTimeSize ++;
          }
        }
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nRunning Split Risk Counts: \n");
          Rprintf("     index    PARrisk    LFTrisk    RGTrisk    PARdeath   LFTdeath   RGTdeath\n");
          for (k=1; k <=  localDeathTimeSize; k++) {
            Rprintf("%10d %10d %10d %10d %10d %10d %10d\n", k,
                    nodeParentAtRisk[k], nodeLeftAtRisk[k], nodeRightAtRisk[k],
                    nodeParentDeath[k], nodeLeftDeath[k], nodeRightDeath[k]);
          }
          Rprintf("\nRunning Split Total LFT & RGT Unique Death Times: \n");
          Rprintf("%10d %10d \n", leftDeathTimeSize, rightDeathTimeSize);
        }          
        if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
          for (k=1; k <= localDeathTimeSize; k++) {
            deltaNum = deltaNum + (
                                   (double) nodeLeftDeath[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentDeath[k]) / nodeParentAtRisk[k])
                                   );
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\nPartial Sum deltaNum:  %10d %10.4f", k, deltaNum);
            }
            if (nodeParentAtRisk[k] >= 2) {
              deltaDen = deltaDen + (
                                     ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                                     (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                                     ((double) (nodeParentAtRisk[k] - nodeParentDeath[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentDeath[k]
                                     );
              if (getTraceFlag() & DL3_TRACE) {
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
          if (delta > deltaMax) {
            deltaMax = delta;
            *splitParameterMax = randomCovariateIndex[i];
            *splitValueMax = j;
          }
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\n\nRunning Split Statistics: \n");
            Rprintf("  SplitParm SplitIndex        Numer        Denom        Delta \n");
            Rprintf(" %10d %10d %12.4f %12.4f %12.4f \n", i, j, deltaNum, deltaDen, delta);
          }
        }  
      }  
    }  
    else {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nMinimum unique deaths not acheived:  %10d versus %10d \n", localDeathTimeSize, _minimumDeathCount);
      }
    }
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum number of deaths");
      Rprintf("\nRSF:  encountered in logRank().  Node not split.  \n");
    }
  }
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeLeftDeath,   1, _masterTimeSize);
  free_uivector(nodeRightDeath,  1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
  free_uivector(nodeLeftAtRisk,   1, _masterTimeSize);
  free_uivector(nodeRightAtRisk,  1, _masterTimeSize);
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, _masterTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
      Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRankRandom() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRankRandom() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
char randomSplit(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  double **masterSplit) {
  double delta, deltaMax;
  uint i,j,k;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nrandomSplit() ENTRY ...\n");
  }
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *leftDeathTimeCount = uivector(1, _masterTimeSize);
  uint *rightDeathTimeCount = uivector(1, _masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, _masterTimeSize);
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = deltaMax = -EPSILON;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for Parent Node: \n");
    for (i=1; i <= localMembershipSize; i++) {
      Rprintf("%10d %10d \n", i, localMembershipIndex[i]);
    }
  }
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeIndex[++localDeathTimeSize] = i;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nLocal Death Times (i, _masterTimeIndex): \n");
      for (i=1; i <= localDeathTimeSize; i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if (localDeathTimeSize >= _minimumDeathCount) {
      actualCovariateCount = selectRandomCovariates(parent, randomCovariateIndex);
      char *covariateStatus = cvector(1, actualCovariateCount);
      for (i = 1; i <= actualCovariateCount; i++) {
        covariateStatus[i] = TRUE;
      }
      i = selectPermissibleElements(actualCovariateCount, covariateStatus);
      while ((i != 0) && ((*splitParameterMax) == 0) && ((*splitValueMax) == 0)) {
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
                  );
        }
        k = (parent -> permissibleSplit)[randomCovariateIndex[i]][2] - (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
        j = floor(ran2(_seed2Ptr)*(k*1.0));
        j = j + (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
        }
        leftDeathCount = 0;
        leftDeathTimeSize = rightDeathTimeSize = 0;
        for (k=1; k <= _masterTimeSize; k++) {
          leftDeathTimeCount[k] = 0;
        }
        for (k=1; k <= localMembershipSize; k++) {
          if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] <= masterSplit[randomCovariateIndex[i]][j]) {
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[k]);
            }
            if (_status[localMembershipIndex[k]] == 1) {
              leftDeathCount ++;
              leftDeathTimeCount[_masterTimeIndex[localMembershipIndex[k]]] ++;
            }
          }
        }  
        for (k=1; k <= _masterTimeSize; k++) {
          rightDeathTimeCount[k] = localDeathTimeCount[k] - leftDeathTimeCount[k];
          if (leftDeathTimeCount[k] > 0) {
            leftDeathTimeSize ++;
          }
          if (rightDeathTimeCount[k] > 0) {
            rightDeathTimeSize ++;
          }
        }
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nRunning Split Totals: \n");
          Rprintf("   PARUdeath   LFTUdeath  RGTUdeath\n");
          Rprintf("%10d %10d %10d \n", localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize);
        }
        if ((leftDeathTimeSize >= (_minimumDeathCount)) && (rightDeathTimeSize >= (_minimumDeathCount))) {
          *splitParameterMax = randomCovariateIndex[i];
          *splitValueMax = j;
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\n\nRunning Split Statistics: \n");
            Rprintf(" SplitParm SplitIndex   SplitValue        Delta \n");
            Rprintf("%10d %10d %12.4f %12.4f \n", i, j, masterSplit[randomCovariateIndex[i]][j], delta);
          }
        }  
        else {
          covariateStatus[i] = FALSE;
          i = selectPermissibleElements(actualCovariateCount, covariateStatus);
        }
      }  
      free_cvector(covariateStatus, 1, actualCovariateCount);
    }  
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum number of deaths");
      Rprintf("\nRSF:  encountered in logRankApprox().  Node not split.  \n");
    }
  }
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(leftDeathTimeCount, 1, _masterTimeSize);
  free_uivector(rightDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, _masterTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
      Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nrandomSplit() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nrandomSplit() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
char logRankApprox(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  double **masterSplit) {
  double delta, deltaNum, deltaDen;
  double eOneMono;
  double deltaMax;
  uint i,j,k, leadingIndex;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nlogRankApprox() ENTRY ...\n");
  }
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *leftDeathTimeCount = uivector(1, _masterTimeSize);
  uint *rightDeathTimeCount = uivector(1, _masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, _masterTimeSize);
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = deltaMax = -EPSILON;
  leadingIndex = 0;  
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for Parent Node: \n");
    for (i=1; i <= localMembershipSize; i++) {
      Rprintf("%10d %10d \n", i, localMembershipIndex[i]);
    }
  }
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeIndex[++localDeathTimeSize] = i;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nLocal Death Times (i, _masterTimeIndex): \n");
      for (i=1; i <= localDeathTimeSize; i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if (localDeathTimeSize >= _minimumDeathCount) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal Death Counts (i, deaths): \n");
      }
      for (i=1; i <= localDeathTimeSize; i++) {
        nodeParentAtRisk[i] = 0;
        nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
        if (getTraceFlag() & DL2_TRACE) {
          Rprintf("%10d %10d \n", i, nodeParentDeath[i]);
        }
        for (j=1; j <= localMembershipSize; j++) {
          if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
            nodeParentAtRisk[i] ++;
          }
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal At Risk Counts (i, at risk): \n");
        for (j=1; j <= localDeathTimeSize; j++) {
          Rprintf("%10d %10d \n", j, nodeParentAtRisk[j]);
        }
      }
      actualCovariateCount = selectRandomCovariates(parent, randomCovariateIndex);
      double *nelsonAalen = dvector(1, localDeathTimeSize);
      uint *survivalTimeIndex = uivector(1, localMembershipSize);
      for (k=1; k <= localDeathTimeSize; k++) {
        if (nodeParentAtRisk[k] != 0) {
          nelsonAalen[k] = (double) nodeParentDeath[k]/nodeParentAtRisk[k];
        }
        else {
          nelsonAalen[k] = 0;
        }
      }
      for (k=2; k <= localDeathTimeSize; k++) {
        nelsonAalen[k] += nelsonAalen[k-1];
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\n     index     NA_par\n");
        for (k=1; k <= localDeathTimeSize-1; k++) {
          Rprintf("%10d %10.4f \n", k, nelsonAalen[k]);
        }
      }
      for (k = 1; k <= localMembershipSize; k++) {
        for (j = 1; j <= localDeathTimeSize; j++) {
          if ( _masterTimeIndex[localMembershipIndex[k]] <= localDeathTimeIndex[j] ) {
            leadingIndex = j;
          }
          else {
            j = localDeathTimeSize;
          }
        }
        survivalTimeIndex[k] = leadingIndex;
      }
      for (i=1; i <= actualCovariateCount; i++) {
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
                  );
        }
        eOneMono = 0.0;
        for (j = (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
             j < (parent -> permissibleSplit)[randomCovariateIndex[i]][2];
             j++) {
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
          }
          leftDeathCount = 0;
          leftDeathTimeSize = rightDeathTimeSize = 0;
          for (k=1; k <= _masterTimeSize; k++) {
            leftDeathTimeCount[k] = 0;
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] <= masterSplit[randomCovariateIndex[i]][j]) {
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[k]);
              }
              if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] == masterSplit[randomCovariateIndex[i]][j]) {
                eOneMono = eOneMono + nelsonAalen[survivalTimeIndex[k]];
              }
              if (_status[localMembershipIndex[k]] == 1) {
                leftDeathCount ++;
                leftDeathTimeCount[_masterTimeIndex[localMembershipIndex[k]]] ++;
              }
            }
          }  
          for (k=1; k <= _masterTimeSize; k++) {
            rightDeathTimeCount[k] = localDeathTimeCount[k] - leftDeathTimeCount[k];
            if (leftDeathTimeCount[k] > 0) {
              leftDeathTimeSize ++;
            }
            if (rightDeathTimeCount[k] > 0) {
              rightDeathTimeSize ++;
            }
          }
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nRunning Split Totals: \n");
            Rprintf("   PARUdeath   LFTUdeath  RGTUdeath\n");
            Rprintf("%10d %10d %10d \n", localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize);
          }
          if ((leftDeathTimeSize >= (_minimumDeathCount)) && (rightDeathTimeSize >= (_minimumDeathCount))) {
            deltaNum = fabs(leftDeathCount - eOneMono);
            if (eOneMono > (parentDeathCount - eOneMono)) {
              deltaDen = eOneMono;
            }
            else {
              deltaDen = (parentDeathCount - eOneMono);
            }
            if (deltaDen <= EPSILON) {
              if (deltaNum <= EPSILON) {
                delta = 0.0;
              }
            }
            else {
              delta = deltaNum / deltaDen;
            }
            if (delta > deltaMax) {
              deltaMax = delta;
              *splitParameterMax = randomCovariateIndex[i];
              *splitValueMax = j;
            }
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\n\nRunning Split Statistics: \n");
              Rprintf(" SplitParm SplitIndex   SplitValue        Delta \n");
              Rprintf("%10d %10d %12.4f %12.4f \n", i, j, masterSplit[randomCovariateIndex[i]][j], delta);
            }
          }  
        }  
      }  
      free_dvector(nelsonAalen , 1, localDeathTimeSize);
      free_uivector (survivalTimeIndex, 1, localMembershipSize);
    }  
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum number of deaths");
      Rprintf("\nRSF:  encountered in logRankApprox().  Node not split.  \n");
    }
  }
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(leftDeathTimeCount, 1, _masterTimeSize);
  free_uivector(rightDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, _masterTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
      Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRankApprox() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRankApprox() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
char logRankScore(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  double **masterSplit) {
  double delta, deltaNum, deltaDen;
  double meanSurvRank, varSurvRank;
  double deltaMax;
  uint i,j,k;
  uint actualCovariateCount;
  uint parentDeathCount;
  uint localMembershipSize, leftMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nlogrankScore() ENTRY ...\n");
  }
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *leftDeathTimeCount = uivector(1, _masterTimeSize);
  uint *rightDeathTimeCount = uivector(1, _masterTimeSize);
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = deltaMax = -EPSILON;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nLocal Membership Size:  %10d \n", localMembershipSize);
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
  }
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        ++localDeathTimeSize;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nLocal Death Time Size:  %10d \n", localDeathTimeSize);
    }
    if (localDeathTimeSize >= _minimumDeathCount) {
      actualCovariateCount = selectRandomCovariates(parent, randomCovariateIndex);
      double *predictorValue = dvector(1, localMembershipSize);
      uint   *localSplitRank  = uivector(1, localMembershipSize);
      uint *survivalTimeIndexRank = uivector(1, localMembershipSize);
      double *survivalRank = dvector(1, localMembershipSize);
      for (i=1; i <= actualCovariateCount; i++) {
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
        if (getTraceFlag() & DL3_TRACE) {
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
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\n\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
                  (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
                  );
        }
        for (j = (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
             j < (parent -> permissibleSplit)[randomCovariateIndex[i]][2];
             j++) {
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
          }
          deltaNum = 0.0;
          leftMembershipSize = leftDeathTimeSize = rightDeathTimeSize = 0;
          for (k=1; k <= _masterTimeSize; k++) {
            leftDeathTimeCount[k] = 0;
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_observation[randomCovariateIndex[i]][localMembershipIndex[localSplitRank[k]]] <= masterSplit[randomCovariateIndex[i]][j]) {
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[localSplitRank[k]]);
              }
              leftMembershipSize ++;
              if (_status[localMembershipIndex[localSplitRank[k]]] == 1) {
                leftDeathTimeCount[_masterTimeIndex[localMembershipIndex[localSplitRank[k]]]] ++;
              }
              deltaNum = deltaNum + survivalRank[k];
            }
          }  
          for (k=1; k <= _masterTimeSize; k++) {
            rightDeathTimeCount[k] = localDeathTimeCount[k] - leftDeathTimeCount[k];
            if (leftDeathTimeCount[k] > 0) {
              leftDeathTimeSize ++;
            }
            if (rightDeathTimeCount[k] > 0) {
              rightDeathTimeSize ++;
            }
          }
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nRunning Split Totals: \n");
            Rprintf("    PARcnt      LFTcnt     RGTcnt  PARUdeath   LFTUdeath  RGTUdeath\n");
            Rprintf("%10d %10d %10d %10d %10d %10d \n", localMembershipSize, leftMembershipSize, localMembershipSize - leftMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize);
          }
          if ((leftDeathTimeSize >= (_minimumDeathCount)) && (rightDeathTimeSize >= (_minimumDeathCount))) {
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
            if (delta > deltaMax) {
              deltaMax = delta;
              *splitParameterMax = randomCovariateIndex[i];
              *splitValueMax = j;
            }
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\n\nRunning Split Statistics: \n");
              Rprintf(" SplitParm SplitIndex   SplitValue        Delta \n");
              Rprintf("%10d %10d %12.4f %12.4f \n", i, j, masterSplit[randomCovariateIndex[i]][j], delta);
            }
          }  
        }  
      }  
      free_dvector(predictorValue, 1, localMembershipSize);
      free_uivector(localSplitRank, 1, localMembershipSize);
      free_uivector(survivalTimeIndexRank, 1, localMembershipSize);
      free_dvector(survivalRank, 1, localMembershipSize);
    }  
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum number of deaths");
      Rprintf("\nRSF:  encountered in logRankScore().  Node not split.  \n");
    }
  }
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(leftDeathTimeCount, 1, _masterTimeSize);
  free_uivector(rightDeathTimeCount, 1, _masterTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
      Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRankScore() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRankScore() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
char conserveEvents(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  double **masterSplit) {
  double delta, nelsonAalenSumLeft, nelsonAalenSumRight;
  double deltaMax;
  uint i,j,k,m, index;
  uint actualCovariateCount;
  uint parentDeathCount;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nconserveEvents() ENTRY ...\n");
  }
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeLeftDeath    = uivector(1, _masterTimeSize);
  uint *nodeRightDeath    = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  uint *nodeLeftAtRisk   = uivector(1, _masterTimeSize);
  uint *nodeRightAtRisk   = uivector(1, _masterTimeSize);
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, _masterTimeSize);
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = deltaMax = -EPSILON;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for Parent Node: \n");
    for (i=1; i <= localMembershipSize; i++) {
      Rprintf("%10d %10d \n", i, localMembershipIndex[i]);
    }
    Rprintf("\nLocal Death Time Counts:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10d \n", i, localDeathTimeCount[i]);
    }
  }
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeIndex[++localDeathTimeSize] = i;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nLocal Death Times (i, _masterTimeIndex): \n");
      for (i=1; i <= localDeathTimeSize; i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if (localDeathTimeSize >= _minimumDeathCount) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal Death Counts (i, deaths): \n");
      }
      for (i=1; i <= localDeathTimeSize; i++) {
        nodeParentAtRisk[i] = 0;
        nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
        if (getTraceFlag() & DL2_TRACE) {
          Rprintf("%10d %10d \n", i, nodeParentDeath[i]);
        }
        for (j=1; j <= localMembershipSize; j++) {
          if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
            nodeParentAtRisk[i] ++;
          }
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal At Risk Counts (i, at risk): \n");
        for (j=1; j <= localDeathTimeSize; j++) {
          Rprintf("%10d %10d \n", j, nodeParentAtRisk[j]);
        }
      }
      actualCovariateCount = selectRandomCovariates(parent, randomCovariateIndex);
      double *nelsonAalenLeft  = dvector(1, localDeathTimeSize);
      double *nelsonAalenRight = dvector(1, localDeathTimeSize);
      for (i=1; i <= actualCovariateCount; i++) {
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\n\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
            (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
            (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
          );
        }
        for (j = (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
             j < (parent -> permissibleSplit)[randomCovariateIndex[i]][2];
             j++) {
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
          }
          nelsonAalenSumLeft = nelsonAalenSumRight = 0.0;
          leftDeathTimeSize = rightDeathTimeSize = 0;
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeLeftDeath[k] = nodeLeftAtRisk[k] = 0;
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] <= masterSplit[randomCovariateIndex[i]][j]) {
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[k]);
              }
              index = 0;  
              for (m = 1; m <= localDeathTimeSize; m++) {
                if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                  if (getTraceFlag() & DL3_TRACE) {
                    Rprintf("NLAR:  %10d %10d %10d \n", m, localDeathTimeIndex[m], _masterTimeIndex[localMembershipIndex[k]]);
                  }
                  nodeLeftAtRisk[m] ++;
                  index = m;
                }
                else {
                  m = localDeathTimeSize;
                }
              }
              if (_status[localMembershipIndex[k]] == 1) {
                nodeLeftDeath[index] ++;
              }
            }
            else {
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nMember of Right Daughter (index):  %10d %10d \n", k, localMembershipIndex[k]);
              }
            }
          }  
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeRightDeath[k] = nodeParentDeath[k] - nodeLeftDeath[k];
            nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k];
            if (nodeLeftDeath[k] > 0) {
              leftDeathTimeSize ++;
            }
            if (nodeRightDeath[k] > 0) {
              rightDeathTimeSize ++;
            }
          }
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nRunning Split Risk Counts: \n");
            Rprintf("     index     PARrisk    LFTrisk    RGTrisk    PARdeath   LFTdeath   RGTdeath\n");
            for (k=1; k <=  localDeathTimeSize; k++) {
              Rprintf("%10d %10d %10d %10d %10d %10d %10d\n", k,
                      nodeParentAtRisk[k], nodeLeftAtRisk[k], nodeRightAtRisk[k],
                      nodeParentDeath[k], nodeLeftDeath[k], nodeRightDeath[k]);
            }
            Rprintf("\nRunning Split Total LFT & RGT Unique Death Times: \n");
            Rprintf("%10d %10d \n", leftDeathTimeSize, rightDeathTimeSize);
          }
          if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
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
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\n     index     LFTrsk     RGTrsk     NA_lft     NA_rgt\n");
              for (k=1; k <= localDeathTimeSize-1; k++) {
                Rprintf("%10d %10d %10d %10.4f %10.4f \n", k, nodeLeftAtRisk[k], nodeRightAtRisk[k], nelsonAalenLeft[k], nelsonAalenRight[k]);
              }
            }
            for (k=1; k <= localDeathTimeSize-1; k++) {
              nelsonAalenSumLeft += (nodeLeftAtRisk[k] - nodeLeftAtRisk[k+1]) * nodeLeftAtRisk[k+1] * nelsonAalenLeft[k];
              nelsonAalenSumRight += (nodeRightAtRisk[k] - nodeRightAtRisk[k+1]) * nodeRightAtRisk[k+1] * nelsonAalenRight[k];
            }
            delta = ((nodeLeftAtRisk[1] * nelsonAalenSumLeft) + (nodeRightAtRisk[1] * nelsonAalenSumRight)) / (nodeLeftAtRisk[1] + nodeRightAtRisk[1]);
            delta = 1.0 / (1.0 + delta);
            if (delta > deltaMax) {
              deltaMax = delta;
              *splitParameterMax = randomCovariateIndex[i];
              *splitValueMax = j;
            }
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\n\nRunning Split Statistics: \n");
              Rprintf(" SplitParm SplitIndex   SplitValue        Delta \n");
              Rprintf("%10d %10d %12.4f %12.4f \n", i, j, masterSplit[randomCovariateIndex[i]][j], delta);
            }
          }  
        }  
      }  
      free_dvector(nelsonAalenLeft, 1, localDeathTimeSize);
      free_dvector(nelsonAalenRight, 1, localDeathTimeSize);
    }  
    else {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nMinimum unique deaths not acheived:  %10d versus %10d \n", localDeathTimeSize, _minimumDeathCount);
      }
    }
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum number of deaths");
      Rprintf("\nRSF:  encountered in conserveEvents().  Node not split.  \n");
    }
  }
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeLeftDeath,   1, _masterTimeSize);
  free_uivector(nodeRightDeath,  1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
  free_uivector(nodeLeftAtRisk,   1, _masterTimeSize);
  free_uivector(nodeRightAtRisk,  1, _masterTimeSize);
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, _masterTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
      Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nconserveEvents() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nconserveEvents() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
char logRank (
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  double **masterSplit) {
  double delta, deltaNum, deltaDen;
  double deltaMax;
  uint i,j,k,m, index;
  uint actualCovariateCount;
  uint parentDeathCount;
  uint localMembershipSize, localDeathTimeSize, leftDeathTimeSize, rightDeathTimeSize;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nlogRank() ENTRY ...\n");
  }
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeLeftDeath    = uivector(1, _masterTimeSize);
  uint *nodeRightDeath   = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  uint *nodeLeftAtRisk   = uivector(1, _masterTimeSize);
  uint *nodeRightAtRisk  = uivector(1, _masterTimeSize);
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, _masterTimeSize);
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = deltaMax = -EPSILON;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for Parent Node: \n");
    for (i=1; i <= localMembershipSize; i++) {
      Rprintf("%10d %10d \n", i, localMembershipIndex[i]);
    }
    Rprintf("\nLocal Death Time Counts:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10d \n", i, localDeathTimeCount[i]);
    }
  }
  if (parentDeathCount >= (2 * (_minimumDeathCount))) {
    for (i=1; i <= _masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeIndex[++localDeathTimeSize] = i;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nLocal Death Times (i, _masterTimeIndex): \n");
      for (i=1; i <= localDeathTimeSize; i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if (localDeathTimeSize >= _minimumDeathCount) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal Death Counts (i, deaths): \n");
      }
      for (i=1; i <= localDeathTimeSize; i++) {
        nodeParentAtRisk[i] = 0;
        nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
        if (getTraceFlag() & DL2_TRACE) {
          Rprintf("%10d %10d \n", i, nodeParentDeath[i]);
        }
        for (j=1; j <= localMembershipSize; j++) {
          if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
            nodeParentAtRisk[i] ++;
          }
        }
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nLocal At Risk Counts (i, at risk): \n");
        for (j=1; j <= localDeathTimeSize; j++) {
          Rprintf("%10d %10d \n", j, nodeParentAtRisk[j]);
        }
      }
      actualCovariateCount = selectRandomCovariates(parent, randomCovariateIndex);
      for (i=1; i <= actualCovariateCount; i++) {
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
            (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
            (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
          );
        }
        for (j = (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
             j < (parent -> permissibleSplit)[randomCovariateIndex[i]][2];
             j++) {
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
          }
          deltaNum = deltaDen =  0.0;
          leftDeathTimeSize = rightDeathTimeSize = 0;
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeLeftDeath[k] = nodeLeftAtRisk[k] = 0;
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] <= masterSplit[randomCovariateIndex[i]][j]) {
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[k]);
              }
              index = 0;  
              for (m = 1; m <= localDeathTimeSize; m++) {
                if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                  if (getTraceFlag() & DL3_TRACE) {
                    Rprintf("NLAR:  %10d %10d %10d \n", m, localDeathTimeIndex[m], _masterTimeIndex[localMembershipIndex[k]]);
                  }
                  nodeLeftAtRisk[m] ++;
                  index = m;
                }
                else {
                  m = localDeathTimeSize;
                }
              }
              if (_status[localMembershipIndex[k]] == 1) {
                nodeLeftDeath[index] ++;
              }
            }
            else {
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nMember of Right Daughter (index):  %10d %10d \n", k, localMembershipIndex[k]);
              }
            }
          }  
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeRightDeath[k] = nodeParentDeath[k] - nodeLeftDeath[k];
            nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k];
            if (nodeLeftDeath[k] > 0) {
              leftDeathTimeSize ++;
            }
            if (nodeRightDeath[k] > 0) {
              rightDeathTimeSize ++;
            }
          }
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nRunning Split Risk Counts: \n");
            Rprintf("     index    PARrisk    LFTrisk    RGTrisk    PARdeath   LFTdeath   RGTdeath\n");
            for (k=1; k <=  localDeathTimeSize; k++) {
              Rprintf("%10d %10d %10d %10d %10d %10d %10d\n", k,
                      nodeParentAtRisk[k], nodeLeftAtRisk[k], nodeRightAtRisk[k],
                      nodeParentDeath[k], nodeLeftDeath[k], nodeRightDeath[k]);
            }
            Rprintf("\nRunning Split Total LFT & RGT Unique Death Times: \n");
            Rprintf("%10d %10d \n", leftDeathTimeSize, rightDeathTimeSize);
          }          
          if ((leftDeathTimeSize  >= (_minimumDeathCount)) && (rightDeathTimeSize  >= (_minimumDeathCount))) {
            for (k=1; k <= localDeathTimeSize; k++) {
              deltaNum = deltaNum + (
                (double) nodeLeftDeath[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentDeath[k]) / nodeParentAtRisk[k])
              );
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nPartial Sum deltaNum:  %10d %10.4f", k, deltaNum);
              }
              if (nodeParentAtRisk[k] >= 2) {
                deltaDen = deltaDen + (
                  ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                  (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                  ((double) (nodeParentAtRisk[k] - nodeParentDeath[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentDeath[k]
                );
                if (getTraceFlag() & DL3_TRACE) {
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
            if (delta > deltaMax) {
              deltaMax = delta;
              *splitParameterMax = randomCovariateIndex[i];
              *splitValueMax = j;
            }
            if (getTraceFlag() & DL3_TRACE) {
              Rprintf("\n\nRunning Split Statistics: \n");
              Rprintf("  SplitParm SplitIndex        Numer        Denom        Delta \n");
              Rprintf(" %10d %10d %12.4f %12.4f %12.4f \n", i, j, deltaNum, deltaDen, delta);
            }
          }  
        }  
      }  
    }  
    else {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nMinimum unique deaths not acheived:  %10d versus %10d \n", localDeathTimeSize, _minimumDeathCount);
      }
    }
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum number of deaths");
      Rprintf("\nRSF:  encountered in logRank().  Node not split.  \n");
    }
  }
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeLeftDeath,   1, _masterTimeSize);
  free_uivector(nodeRightDeath,  1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
  free_uivector(nodeLeftAtRisk,   1, _masterTimeSize);
  free_uivector(nodeRightAtRisk,  1, _masterTimeSize);
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, _masterTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
      Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRank() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nlogRank() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
uint selectRandomCovariates (
  Node *parent,
  uint *covariateIndex) {
  uint i,j,k;
  double randomValue;
  uint maxCovariateCount;
  uint actualCovariateCount;
  uint unselectedCovariateCount;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nselectRandomCovariate() ENTRY ...\n");
  }
  char *randomSplitVector = cvector(1, _xSize);
  double *cdf = dvector(1, _xSize);
  maxCovariateCount = 0;
  for (i=1; i <= _xSize; i++) {
    if (_randomCovariateWeight[i] == 0) {
      randomSplitVector[i] = FALSE;
    }
    else if ((parent -> permissibleSplit)[i][1] == (parent -> permissibleSplit)[i][2]) {
      randomSplitVector[i] = FALSE;
    }
    else {
      randomSplitVector[i] = TRUE;
      maxCovariateCount++;
    }
  }
  if (_randomCovariateCount > maxCovariateCount) {
    actualCovariateCount = maxCovariateCount;
  }
  else {
    actualCovariateCount = _randomCovariateCount;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nCovariate Counts:  actual vs allowed \n");  
    Rprintf("%10d %10d \n", actualCovariateCount, maxCovariateCount);
  }
  for (i=1; i <= actualCovariateCount; i++) {
    unselectedCovariateCount = 0;
    for (k=1; k <= _xSize; k++) {
      if (randomSplitVector[k] == TRUE) {
        cdf[++unselectedCovariateCount] = _randomCovariateWeight[k];
      }
    }
    for (k=2; k <= unselectedCovariateCount; k++) {
      cdf[k] += cdf[k-1];
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nUpdated CDF of Covariate Weights:  ");
      Rprintf("\n     index          cdf");
      for (k=1; k <= unselectedCovariateCount; k++) {
        Rprintf("\n%10d  %12.4f", k, cdf[k]);
      }
      Rprintf("\n");
    }
    randomValue = ran2(_seed2Ptr)*cdf[unselectedCovariateCount];
    j=1;
    while (randomValue > cdf[j]) {
      j++;
    }
    for (k = 1; j > 0; k++) {
      if (randomSplitVector[k] == TRUE) {
        j--;
      }
    }
    randomSplitVector[k-1] = ACTIVE;
    covariateIndex[i] = k-1;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nCovariate Random Selection:  \n");
    for (i=1; i <= actualCovariateCount; i++) {
      Rprintf("%10d %10d \n", i, covariateIndex[i]);
    }
  }
  free_cvector(randomSplitVector, 1, _xSize);
  free_dvector(cdf, 1, _xSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nselectRandomCovariate() EXIT ...\n");
  }
  return actualCovariateCount;
}
uint selectPermissibleElements (uint length,
                                char *permissible) {
  uint permissibleCount;
  uint i, p, index;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nselectPermissibleElements() ENTRY ...\n");
  }
  permissibleCount = 0;
  for (i=1; i <= length; i++) {
    if (permissible[i] == TRUE) {
      permissibleCount ++;
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nPermissible Count:  actual vs maximum \n");  
    Rprintf("%10d %10d \n", permissibleCount, length);
  }
  if (permissibleCount > 0) {
    p = ceil(ran2(_seed2Ptr)*(permissibleCount*1.0));
    index = 1;
    while (p > 0) {
      if (permissible[index] == TRUE) {
        p --;
      }
      index ++;
    }
    index --;
  }
  else {
    index = 0;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nSelected Index:  %10d", index);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nselectPermissibleElements() EXIT ...\n");
  }
  return index;
}
