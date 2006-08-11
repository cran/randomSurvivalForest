//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 1.0.0
//
//  Copyright 2006, Cleveland Clinic
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
#include "node_ops.h"
#include      "rsf.h"
uint      _observationSize;
uint      _xSize;
int      *_seedPtr;
uint     *_status;
uint      _randomCovariateCount;
uint      _minimumDeathCount;
uint      _traceFlag;
uint      _masterTimeSize;
double  **_masterSplit;
double  **_observation;
uint     *_masterTimeIndex;
uint   *_bootMembershipIndex;
uchar  *_bootMembershipFlag;
clock_t start;
clock_t now;
uint updateTimeStamp(uint before) {
  uint stamp;
  double cpuTimeUsed;
  stamp = clock();
  cpuTimeUsed = ((double) (stamp - before)) / CLOCKS_PER_SEC;
  Rprintf("\nRSF:  CPU process time:  %20.3f \n", cpuTimeUsed);
  return stamp;
}
void freeTree(Node *parent) {
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    freeTree(parent -> left);
    freeTree(parent -> right);
  }
  free_Node(parent);
}
Node* getMembership(Node *parent,
                    uint index) {
  uint i;
  Node *result = parent;
  for (i=1; i <= _xSize; i++) {
    if ( (_masterSplit[i][(parent -> permissibleSplit)[i][1]] > _observation[i][index]) ||
         (_masterSplit[i][(parent -> permissibleSplit)[i][2]] < _observation[i][index]) ) {
      Rprintf("\nRSF:  *** WARNING *** \n");
      Rprintf("\nRSF:  Inconsistent call to getMembership().  ");
      Rprintf("\nRSF:  Observation is not in tree.\n");
      return NULL;
    }
  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    if (_observation[parent -> splitParameter][index] <=
        _masterSplit[parent -> splitParameter][parent -> splitValue]) {
      result = getMembership(parent -> left, index);
    }
    else {
      result = getMembership(parent -> right, index);
    }
  }
  return result;
}
char conserveEvents(Node *parent,
                    uint *splitParameterMax,
                    uint *splitValueMax,
                    Node **nodeMembership,
                    uint masterDeathTimeSize) {
  double delta, nelsonAalenSumLeft, nelsonAalenSumRight;
  double deltaMax = -EPSILON;
  uint i,j,k,m, index;
  uint maxCovariateCount = 0;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount, rightDeathCount;
  uint localMembershipSize, localDeathTimeSize;
  uint *nodeParentDeath  = uivector(1, masterDeathTimeSize);
  uint *nodeLeftDeath    = uivector(1, masterDeathTimeSize);
  uint *nodeRightDeath    = uivector(1, masterDeathTimeSize);
  uint *nodeParentAtRisk = uivector(1, masterDeathTimeSize);
  uint *nodeLeftAtRisk   = uivector(1, masterDeathTimeSize);
  uint *nodeRightAtRisk   = uivector(1, masterDeathTimeSize);
  uchar *randomSplitVector = ucvector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, masterDeathTimeSize);
  if (_traceFlag & DL1_TRACE) {
    Rprintf("\nconserveEvents() ENTRY ...\n");
  }
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for parent node: \n");
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
    if (_traceFlag & DL2_TRACE) {
      Rprintf("\nLocal Death Times (i, masterTimeIndex): \n");
      for (i=1; i <= localDeathTimeSize; i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if (localDeathTimeSize > 1) {
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nLocal Death Counts (i, deaths): \n");
      }
      for (i=1; i <= localDeathTimeSize; i++) {
        nodeParentAtRisk[i] = 0;
        nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
        if (_traceFlag & DL2_TRACE) {
          Rprintf("%10d %10d \n", i, nodeParentDeath[i]);
        }
        for (j=1; j <= localMembershipSize; j++) {
          if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
            nodeParentAtRisk[i] ++;
          }
        }
      }
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nLocal At Risk Counts (i, at risk): \n");
        for (j=1; j <= localDeathTimeSize; j++) {
          Rprintf("%10d %10d \n", j, nodeParentAtRisk[j]);
        }
      }
      for (i=1; i <= _xSize; i++) {
        if ((parent -> permissibleSplit)[i][1] == (parent -> permissibleSplit)[i][2]) {
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
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nCovariate Counts:  actual vs allowed \n");  
        Rprintf("%10d %10d \n", actualCovariateCount, maxCovariateCount);
      }
      uint *randomCovariateIndex = uivector(1, actualCovariateCount);
      double *nelsonAalenLeft  = dvector(1, localDeathTimeSize);
      double *nelsonAalenRight = dvector(1, localDeathTimeSize);
      m = actualCovariateCount;
      for (i=maxCovariateCount; m > 0; i--) {
        k = ceil(ran1(_seedPtr)*(i*1.0));
        for (j = 1; k > 0; j++) {
          if (randomSplitVector[j] == TRUE) {
            k--;
          }
        }
        randomSplitVector[j-1] = ACTIVE;
        randomCovariateIndex[m] = j-1;
        m--;
      }
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nCovariate Random Selection:  \n");
        for (i=1; i <= actualCovariateCount; i++) {
          Rprintf("%10d %10d \n", i, randomCovariateIndex[i]);
        }
      }
      for (i=1; i <= actualCovariateCount; i++) {
        if (_traceFlag & DL3_TRACE) {
          Rprintf("\n\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
            (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
            (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
          );
        }
        for (j = (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
             j < (parent -> permissibleSplit)[randomCovariateIndex[i]][2];
             j++) {
          if (_traceFlag & DL3_TRACE) {
            Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
          }
          nelsonAalenSumLeft = nelsonAalenSumRight = 0.0;
          leftDeathCount = rightDeathCount = 0;
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeLeftDeath[k] = nodeLeftAtRisk[k] = 0;
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] <= _masterSplit[randomCovariateIndex[i]][j]) {
              if (_traceFlag & DL3_TRACE) {
                Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[k]);
              }
              index = 0;  
              for (m = 1; m <= localDeathTimeSize; m++) {
                if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                  if (_traceFlag & DL2_TRACE) {
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
          }  
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeRightDeath[k] = nodeParentDeath[k] - nodeLeftDeath[k];
            nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k];
            leftDeathCount += nodeLeftDeath[k];
          }
          rightDeathCount = parentDeathCount - leftDeathCount;
          if (_traceFlag & DL3_TRACE) {
            Rprintf("\nRunning Split Risk Counts: \n");
            Rprintf("     index     PARrisk    LFTrisk    RGTrisk    PARdeath   LFTdeath   RGTdeath\n");
            for (k=1; k <=  localDeathTimeSize; k++) {
              Rprintf("%10d %10d %10d %10d %10d %10d %10d\n", k,
                      nodeParentAtRisk[k], nodeLeftAtRisk[k], nodeRightAtRisk[k],
                      nodeParentDeath[k], nodeLeftDeath[k], nodeRightDeath[k]);
            }
            Rprintf("\nRunning Split Total LFT & RGT Death Counts: \n");
            Rprintf("%10d %10d \n", leftDeathCount, rightDeathCount);
          }
          if ((leftDeathCount  >= (_minimumDeathCount)) && (rightDeathCount  >= (_minimumDeathCount))) {
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
            if (_traceFlag & DL3_TRACE) {
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
            if (_traceFlag & DL3_TRACE) {
              Rprintf("\n\nRunning Split Statistics: \n");
              Rprintf(" SplitParm SplitIndex SplitValue      Delta \n");
              Rprintf("%10d %10d %10.4f %10.4f \n", i, j, _masterSplit[randomCovariateIndex[i]][j], delta);
            }
          }  
        }  
      }  
      free_uivector(randomCovariateIndex, 1, actualCovariateCount);
      free_dvector(nelsonAalenLeft, 1, localDeathTimeSize);
      free_dvector(nelsonAalenRight, 1, localDeathTimeSize);
    }   
  }  
  else {
    if (_traceFlag & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum threshold of deaths");
      Rprintf("\nRSF:  encountered in conserveEvents().  This condition");
      Rprintf("\nRSF:  can occur on the root split and is potentially");
      Rprintf("\nRSF:  nominal behaviour.  \n");
    }
  }
    if (_traceFlag & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex        Delta \n");
      Rprintf(" %10d %10d %12.4f \n", *splitParameterMax, *splitValueMax, deltaMax);
    }
  free_uivector(nodeParentDeath, 1, masterDeathTimeSize);
  free_uivector(nodeLeftDeath,   1, masterDeathTimeSize);
  free_uivector(nodeRightDeath,  1, masterDeathTimeSize);
  free_uivector(nodeParentAtRisk, 1, masterDeathTimeSize);
  free_uivector(nodeLeftAtRisk,   1, masterDeathTimeSize);
  free_uivector(nodeRightAtRisk,  1, masterDeathTimeSize);
  free_ucvector(randomSplitVector, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, masterDeathTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (_traceFlag & DL1_TRACE) {
      Rprintf("\nconserveEvents() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (_traceFlag & DL1_TRACE) {
      Rprintf("\nconserveEvents() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
char logRank(Node *parent,
             uint *splitParameterMax,
             uint *splitValueMax,
             Node **nodeMembership,
             uint masterDeathTimeSize) {
  double delta, deltaNum, deltaDen;
  double deltaMax = -EPSILON;
  uint i,j,k,m, index;
  uint maxCovariateCount = 0;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount, rightDeathCount;
  uint localMembershipSize, localDeathTimeSize;
  uint *nodeParentDeath  = uivector(1, masterDeathTimeSize);
  uint *nodeLeftDeath    = uivector(1, masterDeathTimeSize);
  uint *nodeRightDeath   = uivector(1, masterDeathTimeSize);
  uint *nodeParentAtRisk = uivector(1, masterDeathTimeSize);
  uint *nodeLeftAtRisk   = uivector(1, masterDeathTimeSize);
  uint *nodeRightAtRisk  = uivector(1, masterDeathTimeSize);
  uchar *randomSplitVector = ucvector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, masterDeathTimeSize);
  if (_traceFlag & DL1_TRACE) {
    Rprintf("\nlogRank() ENTRY ...\n");
  }
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = -1.0;  
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (nodeMembership[_bootMembershipIndex[i]] == parent) {
      localMembershipIndex[++localMembershipSize] = _bootMembershipIndex[i];
      if (_status[_bootMembershipIndex[i]] == 1) {
        localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        parentDeathCount++;
      }
    }
  }
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\nParent Death Count:  %10d \n", parentDeathCount);
    Rprintf("\nLocal Membership Index for parent node: \n");
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
    if (_traceFlag & DL2_TRACE) {
      Rprintf("\nLocal Death Times (i, masterTimeIndex): \n");
      for (i=1; i <= localDeathTimeSize; i++) {
        Rprintf("%10d %10d \n", i, localDeathTimeIndex[i]);
      }
    }
    if (localDeathTimeSize > 1) {
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nLocal Death Counts (i, deaths): \n");
      }
      for (i=1; i <= localDeathTimeSize; i++) {
        nodeParentAtRisk[i] = 0;
        nodeParentDeath[i] = localDeathTimeCount[localDeathTimeIndex[i]];
        if (_traceFlag & DL2_TRACE) {
          Rprintf("%10d %10d \n", i, nodeParentDeath[i]);
        }
        for (j=1; j <= localMembershipSize; j++) {
          if (localDeathTimeIndex[i] <= _masterTimeIndex[localMembershipIndex[j]]) {
            nodeParentAtRisk[i] ++;
          }
        }
      }
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nLocal At Risk Counts (i, at risk): \n");
        for (j=1; j <= localDeathTimeSize; j++) {
          Rprintf("%10d %10d \n", j, nodeParentAtRisk[j]);
        }
      }
      for (i=1; i <= _xSize; i++) {
        if ((parent -> permissibleSplit)[i][1] == (parent -> permissibleSplit)[i][2]) {
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
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nCovariate Counts:  actual vs allowed \n");  
        Rprintf("%10d %10d \n", actualCovariateCount, maxCovariateCount);
      }
      uint *randomCovariateIndex = uivector(1, actualCovariateCount);
      m = actualCovariateCount;
      for (i=maxCovariateCount; m > 0; i--) {
        k = ceil(ran1(_seedPtr)*(i*1.0));
        for (j = 1; k > 0; j++) {
          if (randomSplitVector[j] == TRUE) {
            k--;
          }
        }
        randomSplitVector[j-1] = ACTIVE;
        randomCovariateIndex[m] = j-1;
        m--;
      }
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nCovariate Random Selection:  \n");
        for (i=1; i <= actualCovariateCount; i++) {
          Rprintf("%10d %10d \n", i, randomCovariateIndex[i]);
        }
      }
      for (i=1; i <= actualCovariateCount; i++) {
        if (_traceFlag & DL3_TRACE) {
          Rprintf("\nSplitting on parameter:  %10d %10d", i, randomCovariateIndex[i]);
          Rprintf("\n           with limits:  %10d %10d \n",
            (parent -> permissibleSplit)[randomCovariateIndex[i]][1], 
            (parent -> permissibleSplit)[randomCovariateIndex[i]][2]
          );
        }
        for (j = (parent -> permissibleSplit)[randomCovariateIndex[i]][1];
             j < (parent -> permissibleSplit)[randomCovariateIndex[i]][2];
             j++) {
          if (_traceFlag & DL3_TRACE) {
            Rprintf("\nSplitting on parameter, value:  %10d %10d \n", i, j);
          }
          deltaNum = deltaDen =  0.0;
          leftDeathCount = rightDeathCount = 0;
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeLeftDeath[k] = nodeLeftAtRisk[k] = 0;
          }
          for (k=1; k <= localMembershipSize; k++) {
            if (_observation[randomCovariateIndex[i]][localMembershipIndex[k]] <= _masterSplit[randomCovariateIndex[i]][j]) {
              if (_traceFlag & DL3_TRACE) {
                Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, localMembershipIndex[k]);
              }
              index = 0;  
              for (m = 1; m <= localDeathTimeSize; m++) {
                if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                  if (_traceFlag & DL3_TRACE) {
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
              if (_traceFlag & DL3_TRACE) {
                Rprintf("\nMember of Right Daughter (index):  %10d %10d \n", k, localMembershipIndex[k]);
                index = 0;  
                for (m = 1; m <= localDeathTimeSize; m++) {
                  if (localDeathTimeIndex[m] <= _masterTimeIndex[localMembershipIndex[k]]) {
                    Rprintf("NRAR:  %10d %10d %10d \n", m, localDeathTimeIndex[m], _masterTimeIndex[localMembershipIndex[k]]);
                    nodeRightAtRisk[m] ++;
                    index = m;
                  }
                  else {
                    m = localDeathTimeSize;
                  }
                }
                if (_status[localMembershipIndex[k]] == 1) {
                  nodeRightDeath[index] ++;
                }
              }
            }
          }  
          for (k=1; k <= localDeathTimeSize; k++) {
            leftDeathCount += nodeLeftDeath[k];
          }
          rightDeathCount = parentDeathCount - leftDeathCount;
          if (_traceFlag & DL3_TRACE) {
            Rprintf("\nRunning Split Risk Counts: \n");
            Rprintf("     index    PARrisk    LFTrisk    PARdeath   LFTdeath   \n");
            for (k=1; k <=  localDeathTimeSize; k++) {
              Rprintf("%10d %10d %10d %10d %10d %10d %10d\n", k,
                      nodeParentAtRisk[k], nodeLeftAtRisk[k], nodeRightAtRisk[k],
                      nodeParentDeath[k], nodeLeftDeath[k], nodeRightDeath[k]);
            }
            Rprintf("\nRunning Split Total LFT & RGT Death Counts: \n");
            Rprintf("%10d %10d \n", leftDeathCount, rightDeathCount);
          }          
          if ((leftDeathCount  >= (_minimumDeathCount)) && (rightDeathCount  >= (_minimumDeathCount))) {
            for (k=1; k <= localDeathTimeSize; k++) {
              deltaNum = deltaNum + (
                (double) nodeLeftDeath[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentDeath[k]) / nodeParentAtRisk[k])
              );
              if (_traceFlag & DL3_TRACE) {
                Rprintf("\nPartial Sum deltaNum:  %10d %10.4f", k, deltaNum);
              }
              if (nodeParentAtRisk[k] >= 2) {
                deltaDen = deltaDen + (
                  ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                  (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                  ((double) (nodeParentAtRisk[k] - nodeParentDeath[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentDeath[k]
                );
                if (_traceFlag & DL3_TRACE) {
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
            if (_traceFlag & DL3_TRACE) {
              Rprintf("\n\nRunning Split Statistics: \n");
              Rprintf("  SplitParm SplitIndex        Numer        Denom        Delta \n");
              Rprintf(" %10d %10d %12.4f %12.4f %12.4f \n", i, j, deltaNum, deltaDen, delta);
            }
          }  
        }  
      }  
      free_uivector(randomCovariateIndex, 1, actualCovariateCount);
    }  
  }  
  else {
    if (_traceFlag & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum threshold of deaths");
      Rprintf("\nRSF:  encountered in logRank().  This condition");
      Rprintf("\nRSF:  can occur on the root split and is potentially");
      Rprintf("\nRSF:  nominal behaviour.  \n");
    }
  }
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\nBest Split Statistics: \n");
    Rprintf("  SplitParm SplitIndex        Delta \n");
    Rprintf(" %10d %10d %12.4f \n", *splitParameterMax, *splitValueMax, deltaMax);
  }
  free_uivector(nodeParentDeath, 1, masterDeathTimeSize);
  free_uivector(nodeLeftDeath,   1, masterDeathTimeSize);
  free_uivector(nodeParentAtRisk, 1, masterDeathTimeSize);
  free_uivector(nodeLeftAtRisk,   1, masterDeathTimeSize);
  free_ucvector(randomSplitVector, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, masterDeathTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
    if (_traceFlag & DL1_TRACE) {
      Rprintf("\nlogRank() EXIT (TRUE) ...\n");
    }
    return TRUE;
  }
  else {
    if (_traceFlag & DL1_TRACE) {
      Rprintf("\nlogRank() EXIT (FALSE) ...\n");
    }
    return FALSE;
  }
}
char getBestSplit(Node *parent,
                 uint *splitParameterMax,
                 uint *splitValueMax,
                 Node **nodeMembership,
                 uint masterDeathTimeSize,
                 uint *splitRule) {
  if (*splitRule == LOG_RANK) {
    return logRank(parent,
                   splitParameterMax,
                   splitValueMax,
                   nodeMembership,
                   masterDeathTimeSize);
  }
  else if (*splitRule == CONSERVE_EVENTS) {
    return conserveEvents(parent,
                          splitParameterMax,
                          splitValueMax,
                          nodeMembership,
                          masterDeathTimeSize);
  }
  else {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid split rule:  %10d", *splitRule);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
}
char forkAndUpdate(Node **nodeMembership,
                   uint *leafCount,
                   Node *parent,
                   uint splitParameter,
                   uint splitValue) {
  uint i;
  uint leftDeathCount, rightDeathCount;
  char result;
  if (_traceFlag & DL1_TRACE) {
    Rprintf("\ngetforkAndUpdate() ENTRY ...\n");
  }
  result = forkNode(parent, splitParameter, splitValue);
  if (result == TRUE) {
    leftDeathCount = rightDeathCount = 0;
    (*leafCount)++;
    for (i = 1; i <= _observationSize; i++) {
      if (nodeMembership[i] == parent) {
        if (_observation[splitParameter][i] <= _masterSplit[splitParameter][splitValue]) {
          nodeMembership[i] = parent -> left;
          ((parent -> left) -> leafCount) = (parent -> leafCount);
        }
        else {
          nodeMembership[i] = parent -> right;
          ((parent -> right) -> leafCount) = *leafCount;
        }
      }
    }
    for (i = 1; i <= _observationSize; i++) {
      if (nodeMembership[_bootMembershipIndex[i]] == parent -> left) {
        if (_status[_bootMembershipIndex[i]] == 1) leftDeathCount ++;
      }
      else if (nodeMembership[_bootMembershipIndex[i]] == parent -> right) {
        if (_status[_bootMembershipIndex[i]] == 1) rightDeathCount ++;
      }
    }
    if (leftDeathCount  < (2 * (_minimumDeathCount))) {
      (parent -> left) -> splitFlag = FALSE;
    }
    if (rightDeathCount <= (2 * (_minimumDeathCount))) {
      (parent -> right)-> splitFlag = FALSE;
    }
  }
  else {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  forkNode() failed.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (_traceFlag & DL1_TRACE) {
    Rprintf("\ngetforkAndUpdate() EXIT ...\n");
  }
  return result;
}
char makeTree(Node *parent,
              Node **nodeMembership,
              uint *leafCount, 
              uint masterDeathTimeSize,
              uint *splitRule) {
  char splitResult, forkResult;
  uint splitParameterMax, splitValueMax;
  uint i;
  if (_traceFlag & DL1_TRACE) {
    Rprintf("\nmakeTree() ENTRY ...\n");
  }
  splitResult = getBestSplit(parent,
                             &splitParameterMax,
                             &splitValueMax,
                             nodeMembership,
                             masterDeathTimeSize,
                             splitRule);
  if (splitResult == TRUE) {
    forkResult = forkAndUpdate(nodeMembership,
                               leafCount,
                               parent,
                               splitParameterMax,
                               splitValueMax);
    if (forkResult == TRUE) {
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nNode Membership:  \n");
        for (i=1; i <=  _observationSize; i++) {
          Rprintf("%10d %10d \n", i, nodeMembership[i] -> leafCount);
        }
      }
      if ((parent -> left) -> splitFlag == TRUE) {
        if (_traceFlag & DL1_TRACE) {
          Rprintf("\nmakeTree() LEFT:  \n");
        }
        makeTree(parent -> left,
                 nodeMembership,
                 leafCount,
                 masterDeathTimeSize,
                 splitRule);
      }
      if ((parent -> right) -> splitFlag == TRUE) {
        if (_traceFlag & DL1_TRACE) {
          Rprintf("\nmakeTree() RIGHT:  \n");
        }
        makeTree(parent -> right,
                 nodeMembership,
                 leafCount,
                 masterDeathTimeSize,
                 splitRule);
      }
    }
    else {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  forkAndUpdate() failed.");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
    }
  }
  else {
    parent -> splitFlag = FALSE;
    if (_traceFlag & DL1_TRACE) {
      Rprintf("\ngetBestSplit() FAILED ...\n");
    }
  }
  if (_traceFlag & DL1_TRACE) {
    Rprintf("\nmakeTree() EXIT ...\n");
  }
  return splitResult;
}
void getCumulativeHazardEstimate(double **cumulativeHazard,
                                 Node *parent,
                                 double *timeInterest,
                                 uint sortedTimeInterestSize,
                                 Node **nodeMembership,
                                 double *masterTime) {
  uint i, j;
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  uint priorTimePointIndex, currentTimePointIndex;
  double estimate;
  for (i=1; i <= _masterTimeSize; i++) {
    nodeParentAtRisk[i] = nodeParentDeath[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (nodeMembership[_bootMembershipIndex[i]] == parent) {
      for (j=1; j <= _masterTimeIndex[_bootMembershipIndex[i]]; j++) {
        nodeParentAtRisk[j] ++;
      }
      if (_status[_bootMembershipIndex[i]] == 1) {
        nodeParentDeath[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
      }
    }
  }
  priorTimePointIndex = 0;
  currentTimePointIndex = 1;
  for (j=1; j <= sortedTimeInterestSize; j++) {
    for (i = priorTimePointIndex + 1; i <= _masterTimeSize; i++) {
      if (masterTime[i] <= timeInterest[j]) {
        currentTimePointIndex = i;
      }
      else {
        i = _masterTimeSize;
      }
    }
    if (_traceFlag & DL2_TRACE) {
      Rprintf("\nCumulative hazard at risk and death counts for node:  %10d \n", parent -> leafCount);
      Rprintf("  with time point index:  %10d \n", currentTimePointIndex);
      for (i=1; i <= _masterTimeSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (i=1; i <= _masterTimeSize; i++) {
        Rprintf("%10d", nodeParentAtRisk[i]);
      }
      Rprintf("\n");
      for (i=1; i <= _masterTimeSize; i++) {
        Rprintf("%10d", nodeParentDeath[i]);
      }
      Rprintf("\n");
    }
    estimate = 0.0;
    for (i = priorTimePointIndex + 1; i <= currentTimePointIndex; i++) {
      if (nodeParentDeath[i] > 0) {
        if (nodeParentAtRisk[i] >= 1) {
          estimate = estimate + ((double) nodeParentDeath[i] / nodeParentAtRisk[i]);
        }
      }
    }
    cumulativeHazard[j][parent -> leafCount] = estimate;
    priorTimePointIndex = currentTimePointIndex;
  }  
  for (j=2; j <= sortedTimeInterestSize; j++) {
    cumulativeHazard[j][parent -> leafCount] += cumulativeHazard[j-1][parent -> leafCount];
  }
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
}
void rsf(uint *traceFlag,
         int  *memoryUseProtocol,
         int  *seed,
         uint *splitRule,
         uint *randomCovariateCount,
         uint *bootstrapSampleCount,
         uint *minimumDeathCount,
         uint *observationSizePtr,
         double *time,
         uint *status,
         uint *timeInterestSize,
         double *timeInterest,
         uint *xSizePtr,
         double *xData,
         double *ensembleEstimator,
         double *performanceMeasure,
         uint *leafCount,
         uint *proximity) {
  uint i,j,k,b;
  uint leadingIndex;
  Node **root;
  double **ensembleEstimatorPtr;
  uint sortedTimeInterestSize;
  uint *masterSplitSize;
  double *masterTime;
  uint masterDeathTimeSize;
  uint **masterSplitBounds;
  uint proximitySize;
  Node **nodeMembership;
  double **cumulativeHazard;
  uint inBootstrapSampleSize;
  ulong concordancePairSize;
  ulong concordanceWorseCount;
  double *ensembleEstimatorRun;
  uint *ensembleEstimatorDen;
  uint *oobSampleSize;
  uint discardedSampleCount;
  _traceFlag = *traceFlag;
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Native code normal entry. \n");
  }
  start = clock();
  if (*seed >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Random seed must be less than zero.  \n");
    return;
  }
  else if ((*splitRule != LOG_RANK) && (*splitRule != CONSERVE_EVENTS)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Invalid split rule:  %10d  \n", *splitRule);
    return;
  }
  else if (*observationSizePtr < 2) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of observations must be greater than one.  \n");
    return;
  }
  else if (*xSizePtr < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of parameters must be greater than zero.  \n");
    return;
  }
  else if ((*randomCovariateCount < 1) || (*randomCovariateCount > *xSizePtr)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of random covariate parameters must be greater");
    Rprintf("\nRSF:  than zero and less than the total number of covariates.  \n");
    return;
  }
  else if (*minimumDeathCount < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Minimum number of deaths must be greater than zero.  \n");
    return;
  }
  else if (*bootstrapSampleCount < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of bootstrap iterations must be greater than zero.  \n");
    return;
  }
  else if (*timeInterestSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of time points of interest must be greater than zero.  \n");
    return;
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Number of observations:  %10d", *observationSizePtr);
    Rprintf("\nRSF:  Number of predictors:    %10d \n", *xSizePtr);
  }
  _observationSize = *observationSizePtr;
  _xSize = *xSizePtr;
  _seedPtr = seed;
  _status = status;
  _randomCovariateCount = *randomCovariateCount;
  _minimumDeathCount = *minimumDeathCount;
  root = NULL;  
  if (*memoryUseProtocol != FALSE) {
    root = nodePtrVector(1, *bootstrapSampleCount);
  }
  nodeMembership = nodePtrVector(1, _observationSize);
  _bootMembershipIndex = uivector(1, _observationSize);
  _bootMembershipFlag = ucvector(1, _observationSize);
  _observation = dmatrix(1, _xSize, 1, _observationSize);
  _masterSplit = dmatrix(1, _xSize, 1, _observationSize);
  masterSplitBounds = uimatrix(1, _xSize, 1, 2);
  masterTime  = dvector(1, _observationSize);
  _masterTimeIndex  = uivector(1, _observationSize);
  masterSplitSize = uivector(1, _xSize);
  oobSampleSize = uivector(1, *bootstrapSampleCount);
  time--;
  _status--;
  timeInterest--;
  performanceMeasure--;
  leafCount--;
  proximity--;
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\nIncoming Data:  ");
    Rprintf("\n     index       time     status   observations -> \n");
  }
  for (i=1; i <= _observationSize; i++) {
    if (_traceFlag & DL2_TRACE) {
      Rprintf("%10d %10.4f %10d", i, time[i], _status[i]);
      for (j=1; j <= _xSize; j++) {
        Rprintf(" %10.4f", (xData+((j-1)*(_observationSize)))[i-1]);
      }
      Rprintf("\n");
    }
    masterTime[i] = time[i];
    for (j=1; j <= _xSize; j++) {
      _observation[j][i] = (xData+((j-1)*(_observationSize)))[i-1];
    }
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Initial read complete.");  
  }
  hpsort(masterTime, _observationSize);
  _masterTimeSize = leadingIndex = 1;
  for (i=2; i <= _observationSize; i++) {
    if (masterTime[i] > masterTime[leadingIndex]) {
      _masterTimeSize++;
      leadingIndex++;
      masterTime[leadingIndex] = masterTime[i];
    }
  }
  for (i=_masterTimeSize+1; i <= _observationSize; i++) {
    masterTime[i] = 0;
  }
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\n\nSorted Distinct Event Times:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10.4f \n", i, masterTime[i]);
    }
  }
  for (j=1; j <= _observationSize; j++) {
    k = 1;
    while (k <= _masterTimeSize) {
      if (time[j] == masterTime[k]) {
        _masterTimeIndex[j] = k;
        k = _masterTimeSize;
      }
      k++;
    }
  }
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\nMaster Time Index:  \n");
    for (i=1; i <= _observationSize; i++) {
      Rprintf("%10d %10d \n", i, _masterTimeIndex[i]);
    }
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of master time data complete.");  
  }
  hpsort(timeInterest, *timeInterestSize);
  sortedTimeInterestSize = leadingIndex = 1;
  for (i=2; i <= *timeInterestSize; i++) {
    if (timeInterest[i] > timeInterest[leadingIndex]) {
      sortedTimeInterestSize++;
      leadingIndex++;
      timeInterest[leadingIndex] = timeInterest[i];
    }
  }
  if (sortedTimeInterestSize != *timeInterestSize) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Time points of interest are not unique.");
    Rprintf("\nRSF:  The ensemble estimate output matrix is being");
    Rprintf("\nRSF:  resized as [N'] x [n], where N' is the");
    Rprintf("\nRSF:  unique time points of interest and n is");
    Rprintf("\nRSF:  number of observations in the data.");
  }
  for (i=sortedTimeInterestSize+1; i <= *timeInterestSize; i++) {
    timeInterest[i] = 0;
  }
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\n\nSorted Distinct Times of Interest:  \n");
    for (i=1; i <= sortedTimeInterestSize; i++) {
      Rprintf("%10d %10.4f \n", i, timeInterest[i]);
    }
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of time interest data complete.");  
  }
  uchar *masterDeathTimeIndicator = ucvector(1, _masterTimeSize);
  for (i=1; i <= _masterTimeSize; i++) {
    masterDeathTimeIndicator[i] = FALSE;
  }
  k = 0;
  for (i=1; i <= _observationSize; i++) {
    if (_status[i] == 0) {
      k++;
    }
    else {
      masterDeathTimeIndicator[_masterTimeIndex[i]] = TRUE;
    }
  }
  if (k == _observationSize) {
    for (i=1; i <= _observationSize; i++) {
      _status[i] = 1;
    }
    for (i=1; i <= _masterTimeSize; i++) {
      masterDeathTimeIndicator[i] = TRUE;
    }
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  All observations were censored, and will be considered to be deaths.");
  }
  masterDeathTimeSize = 0;
  for (i=1; i<= _masterTimeSize; i++) {
    if (masterDeathTimeIndicator[i] == TRUE) {
      masterDeathTimeSize ++;
    }
  }
  if (_traceFlag & DL2_TRACE) {
    Rprintf("\n\nMaster Death Time Indicator:  \n");
    for (i=1; i <= _masterTimeSize; i++) {
      Rprintf("%10d %10.4f %10d \n", i, masterTime[i], masterDeathTimeIndicator[i]);
    }
  }
  free_ucvector(masterDeathTimeIndicator, 1, _masterTimeSize);
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of death time data complete");  
  }
  ensembleEstimatorRun = dvector(1, _observationSize);
  ensembleEstimatorDen = uivector(1, _observationSize);
  ensembleEstimatorPtr = pdvector(1, sortedTimeInterestSize);
  for (i=1; i <= sortedTimeInterestSize; i++) {
    ensembleEstimatorPtr[i] = ensembleEstimator+((i-1)*(_observationSize)) - 1;
  }
  for (i=1; i <= _observationSize; i++) {
    for (j=1; j <= sortedTimeInterestSize; j++) {
      ensembleEstimatorPtr[j][i] = 0.0;
    }
    ensembleEstimatorDen[i] = 0;
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of ensemble estimates complete.");  
  }
  if (*memoryUseProtocol != FALSE) {
    proximitySize = ((_observationSize + 1)  * _observationSize) / 2; 
    for (i=1; i <= proximitySize; i++) {
      proximity[i] = 0;
    }
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of proximity information complete.");  
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Number of bootstrap samples:  %10d \n", *bootstrapSampleCount);
  }
  discardedSampleCount = 0;
  for (b=1; b <= *bootstrapSampleCount; b++) {
    Node *rootPtr;
    char result;
    leafCount[b] = 1;
    oobSampleSize[b] = 0;
    performanceMeasure[b] = 0;
    if (_traceFlag & DL0_TRACE) {
      Rprintf("\nRSF:  Bootstrap iteration:  %10d ", b);  
    }
    rootPtr = makeNode();
    if (*memoryUseProtocol != FALSE) {
      root[b] = rootPtr;
    }
    rootPtr -> parent = rootPtr;
    rootPtr -> left = NULL; rootPtr -> right = NULL;
    rootPtr -> splitValue = 0;
    rootPtr -> splitParameter = 0;
    rootPtr -> splitFlag = TRUE;
    rootPtr -> leafCount = 1;
    for (i=1; i <= _observationSize; i++) {
      nodeMembership[i] = rootPtr;
      _bootMembershipFlag[i] = FALSE;
    }
    for (i=1; i <= _observationSize; i++) {
      k = ceil(ran1(_seedPtr)*((_observationSize)*1.0));
      _bootMembershipFlag[k] = TRUE;
      _bootMembershipIndex[i] = k;
    }
    for (i=1; i <= _observationSize; i++) {
      if (_bootMembershipFlag[i] == FALSE) {
        oobSampleSize[b] ++;
      }
    }
    if (_traceFlag & DL2_TRACE) {
      Rprintf("\n\nSample Size:  %10d ", oobSampleSize[b]);
      Rprintf("\n\nRoot Membership:  \n");
      for (i=1; i <=  _observationSize; i++) {
        Rprintf("%10d %10d %10d \n", i, _bootMembershipIndex[i], nodeMembership[i] -> leafCount);
      }
    }
    if (_traceFlag & DL0_TRACE) {
      Rprintf("\nRSF:  Initialization of bootstrap root node complete.");  
    }
    if (oobSampleSize[b] != 0) {
      inBootstrapSampleSize = _observationSize - oobSampleSize[b];
      if (_traceFlag & DL3_TRACE) {
        Rprintf("\nRaw Bootstrap MasterSplit Data:  ");  
      }
      k = 1;
      for (i=1; i <= _observationSize; i++) {
        if (_bootMembershipFlag[i] == TRUE) {
          for (j=1; j <= _xSize; j++) {
            _masterSplit[j][k] = (xData+((j-1)*(_observationSize)))[i-1];
          }
          if (_traceFlag & DL3_TRACE) {
            Rprintf("\n%10d ", i);  
            for (j=1; j <= _xSize; j++) {
              Rprintf("%10.4f ", _masterSplit[j][k]);
            }
          }
          k++;
        }
      }
      for (i=1; i <= _xSize; i++) {
        hpsort(_masterSplit[i], inBootstrapSampleSize);
      }
      if (_traceFlag & DL3_TRACE) {
        Rprintf("\n\nRaw Sorted Bootstrap MasterSplit Data:  ");
        for (i=1; i <= inBootstrapSampleSize; i++) {
          Rprintf("\n%10d", i);
          for (j=1; j <= _xSize; j++) {
            Rprintf(" %10.4f", _masterSplit[j][i]);
          }
        }
      }
      for (i=1; i <= _xSize; i++) {
        masterSplitBounds[i][1] = 1;
        masterSplitSize[i] = leadingIndex = 1;
        for (j=2; j <= inBootstrapSampleSize; j++) {
          if (_masterSplit[i][j] > _masterSplit[i][leadingIndex]) {
            masterSplitSize[i]++;
            leadingIndex++;
            _masterSplit[i][leadingIndex] = _masterSplit[i][j];
          }
        }
        masterSplitBounds[i][2] = leadingIndex;
        for (j=masterSplitSize[i]+1; j <= inBootstrapSampleSize; j++) {
          _masterSplit[i][j] = 0;
        }
      }
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\n\nBootstrap MasterSplit Data:  ");
        for (i=1; i <= inBootstrapSampleSize; i++) {
          Rprintf("\n%10d", i);
          for (j=1; j <= _xSize; j++) {
            Rprintf(" %10.4f", _masterSplit[j][i]);
          }
        }
        Rprintf("\n\nSize of Permissible Splits:  ");
        Rprintf("\n          ");
        for (j=1; j <= _xSize; j++) {
          Rprintf(" %10d", masterSplitSize[j]);
        }
        Rprintf("\n");
      }
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nMaster Permissible Split Boundaries:  \n");
        for (i=1; i <= _xSize; i++) {
          Rprintf("%10d %10d %10d \n", i, masterSplitBounds[i][1], masterSplitBounds[i][2]);
        }
      }
      copyMatrix(rootPtr -> permissibleSplit, masterSplitBounds, _xSize, 2);
      if (_traceFlag & DL0_TRACE) {
        Rprintf("\nRSF:  Initialization of master split data complete.");  
      }
      result = makeTree(rootPtr, nodeMembership, &leafCount[b],
                               masterDeathTimeSize, splitRule);
      if (_traceFlag & DL0_TRACE) {
        Rprintf("\nRSF:  Bootstrap iteration tree complete.");  
        Rprintf("\nRSF:  Final leaf count:  %10d ", leafCount[b]);
        now = updateTimeStamp(now);
      }
      if (*memoryUseProtocol != FALSE) {
        k = 0;
        for (i = 1; i <= _observationSize; i++) {
          k += i - 1;
          for (j = 1; j <= i; j++) {
            if ( (nodeMembership[i] -> leafCount) == (nodeMembership[j] -> leafCount) ) {
              proximity[k + j] ++;
            }
          }
        }
        if (_traceFlag & DL0_TRACE) {
          Rprintf("\nRSF:  Proximity matrix calculation complete.");  
        }
      }
      if (_traceFlag & DL2_TRACE) {
        Rprintf("\nFinal Membership for Bootstrap Sample:  \n");
        for (i=1; i <=  _observationSize; i++) {
          Rprintf("%10d %10d %10d \n", i, _bootMembershipIndex[i], nodeMembership[i] -> leafCount);
        }
      }
      if (result == TRUE) {
        cumulativeHazard = dmatrix(1, sortedTimeInterestSize, 1, leafCount[b]);
        for (i=1; i <= leafCount[b]; i++) {
          for (j = 1; j <= _observationSize; j++) {
            if (_bootMembershipFlag[j] == TRUE) {
              if ((nodeMembership[j] -> leafCount) == i) {
                k = j;
                j = _observationSize;
              }
            }
          }
          getCumulativeHazardEstimate(cumulativeHazard,
                                      nodeMembership[k],
                                      timeInterest,
                                      sortedTimeInterestSize,
                                      nodeMembership,
                                      masterTime
          );
        }
        if (_traceFlag & DL0_TRACE) {
          Rprintf("\nRSF:  CHE calculation complete.");  
        }
        if (_traceFlag & DL2_TRACE) {
          Rprintf("\nTree specific cumulative hazard calculation: \n");
          Rprintf("          ");
          for (i=1; i <= leafCount[b]; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n");
          for (j=1; j <= sortedTimeInterestSize; j++) {
            Rprintf("%10d", j);
            for (i=1; i <= leafCount[b]; i++) {
              Rprintf("%10.4f", cumulativeHazard[j][i]);
            }
            Rprintf("\n");
          }
        }
        for (i=1; i <= _observationSize; i++) {
          if (_bootMembershipFlag[i] == FALSE) {
            k = nodeMembership[i] -> leafCount;
            for (j=1; j <= sortedTimeInterestSize; j++) {
              ensembleEstimatorPtr[j][i] += cumulativeHazard[j][k];
            }
            ensembleEstimatorDen[i] ++;
          }
        }
        free_dmatrix(cumulativeHazard, 1, sortedTimeInterestSize, 1, leafCount[b]);
        if (_traceFlag & DL2_TRACE) {
          Rprintf("\nEnsemble Estimator Numerator calculation: \n");
          Rprintf("          ");
          for (i=1; i <= _observationSize; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n");
          for (j=1; j <= sortedTimeInterestSize; j++) {
            Rprintf("%10d", j);
            for (i=1; i <= _observationSize; i++) {
              Rprintf("%10.4f", ensembleEstimatorPtr[j][i]);
            }
            Rprintf("\n");
          }
        }
        for (i = 1; i <= _observationSize; i++) {
          ensembleEstimatorRun[i] = 0.0;
          for (j=1; j <= sortedTimeInterestSize; j++) {
            ensembleEstimatorRun[i] += ensembleEstimatorPtr[j][i];
          }
          if (ensembleEstimatorDen[i] != 0) {
            ensembleEstimatorRun[i] = ensembleEstimatorRun[i] / ensembleEstimatorDen[i];
          }
        }
        if (_traceFlag & DL0_TRACE) {
          Rprintf("\nRSF:  Bootstrap iteration ensemble estimate complete.");  
        }
        if (_traceFlag & DL2_TRACE) {
          Rprintf("\nEnsemble Estimator Denominator calculation: \n");
          Rprintf("          ");  
          for (i=1; i <= _observationSize; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n          ");
          for (i=1; i <= _observationSize; i++) {
            Rprintf("%10d", ensembleEstimatorDen[i]);
          }
          Rprintf("\n");
          Rprintf("\nEnsemble Estimator:  %10d \n", b);
          Rprintf("          ");
          for (i=1; i <= _observationSize; i++) {
            Rprintf("%10d", i);
          }
          Rprintf("\n          ");
          for (i=1; i <= _observationSize; i++) {
            Rprintf("%10.4f", ensembleEstimatorRun[i]);
          }
          Rprintf("\n");
        }
        concordancePairSize = concordanceWorseCount = 0;
        for (i=1; i < _observationSize; i++) {
          for (j=i+1; j <= _observationSize; j++) {
            if (ensembleEstimatorDen[i] != 0  && ensembleEstimatorDen[j] != 0) {
              if ( (_masterTimeIndex[i] >  _masterTimeIndex[j] && _status[j] == 1) ||
                   (_masterTimeIndex[i] == _masterTimeIndex[j] && _status[j] == 1 && _status[i] == 0) ) {
                concordancePairSize += 2;
                if (ensembleEstimatorRun[j] > ensembleEstimatorRun[i]) {              
                    concordanceWorseCount += 2;
                  }
                else if (fabs(ensembleEstimatorRun[j] - ensembleEstimatorRun[i]) < EPSILON) {              
                  concordanceWorseCount += 1;
                }  
              }
              else if ( (_masterTimeIndex[j] >  _masterTimeIndex[i] && _status[i] == 1) ||
                        (_masterTimeIndex[j] == _masterTimeIndex[i] && _status[i] == 1 && _status[j] == 0) ) {
                concordancePairSize += 2;
                if (ensembleEstimatorRun[i] > ensembleEstimatorRun[j]) {
                  concordanceWorseCount += 2;
                }
                else if (fabs(ensembleEstimatorRun[i] - ensembleEstimatorRun[j]) < EPSILON) {
                  concordanceWorseCount += 1;
                }  
              }
            }
          }
        }
        performanceMeasure[b] = 1.0 - ((double) concordanceWorseCount / concordancePairSize);
        if (_traceFlag & DL0_TRACE) {
          Rprintf("\nRSF:  Concordance pair and error update complete.");  
          Rprintf("\nRSF:  Number of concordance pairs:         %20d", concordancePairSize);
          Rprintf("\nRSF:  Number of pairs with worse outcome:  %20d", concordanceWorseCount);
          Rprintf("\nRSF:  Error Rate:                          %20.4f ", performanceMeasure[b]);
          Rprintf("\n");
        }
        if (_traceFlag & DL2_TRACE) {
          Rprintf("\nDiagnostic Classification of OOB elements:   \n");
          Rprintf("\n       Leaf     Observ         EE    Data ->  \n");
          for (i=1; i <= leafCount[b]; i++) {
            for (j = 1; j <= _observationSize; j++) {
              if (_bootMembershipFlag[j] == FALSE) {
                if ((nodeMembership[j] -> leafCount) == i) {
                  Rprintf("\n %10d %10d %10.4f", i, j, ensembleEstimatorRun[j]);
                  for (k = 1; k <= *xSizePtr; k++) {
                    Rprintf("%10.4f", _observation[k][j]);
                  }
                }
              }
            }
            Rprintf("\n");
          }
        }
      }
      else {
        discardedSampleCount++;
        if (b > 1) {
          performanceMeasure[b] = performanceMeasure[b-1];
        }
      }  
    }  
    else {
      discardedSampleCount++;
      if (b > 1) {
        performanceMeasure[b] = performanceMeasure[b-1];
      }
    }
    if (*memoryUseProtocol == FALSE) {
      freeTree(rootPtr);
    }
  }  
  if (_traceFlag & DL2_TRACE) {
    if (*memoryUseProtocol != FALSE) {
      Rprintf("\nProximity Matrix:  \n");
      k = 0;
      for (i = 1; i <= _observationSize; i++) {
        k += i - 1;
        for (j = 1; j <= i; j++) {
          Rprintf("%10d ", proximity[k + j]);
        }
        Rprintf("\n");
      }
    }
  }
  k = 0;
  for (b=1; b <= *bootstrapSampleCount; b++) {
    if (oobSampleSize[b] == 0) {
      k++;
    }
  }
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Bootstrap samples:                                                  %10d ", *bootstrapSampleCount);
    Rprintf("\nRSF:  Bootstrap samples rejected (total):                                 %10d " , discardedSampleCount);
    Rprintf("\nRSF:  Bootstrap samples rejected (OOB subset size was zero):              %10d ", k);
    Rprintf("\nRSF:  Bootstrap samples rejected (threshold of minimum deaths too high):  %10d \n", discardedSampleCount - k);
    Rprintf("\n");
  }
  if (discardedSampleCount == *bootstrapSampleCount) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Insufficient trees for analysis.  \n");
  }
  for (i = 1; i <= _observationSize; i++) {
    if (ensembleEstimatorDen[i] != 0) {
      for (j=1; j <= sortedTimeInterestSize; j++) {
        ensembleEstimatorPtr[j][i] = ensembleEstimatorPtr[j][i] / ensembleEstimatorDen[i];
      }
    }
  }
  if (*memoryUseProtocol != FALSE) {
    for (b=1; b <= *bootstrapSampleCount; b++) {
      freeTree(root[b]);
    }
    free_nodePtrVector(root, 1, *bootstrapSampleCount);
  }
  free_nodePtrVector(nodeMembership, 1, _observationSize);
  free_uivector(_bootMembershipIndex, 1, _observationSize);
  free_ucvector(_bootMembershipFlag, 1, _observationSize);
  free_dmatrix(_observation, 1, _xSize, 1, _observationSize);
  free_dmatrix(_masterSplit, 1, _xSize, 1, _observationSize);
  free_uimatrix(masterSplitBounds, 1, _xSize, 1, 2);
  free_dvector(masterTime, 1, _observationSize);
  free_uivector(_masterTimeIndex, 1, _observationSize);
  free_uivector(masterSplitSize, 1, _xSize);
  free_uivector(oobSampleSize, 1, *bootstrapSampleCount);
  free_dvector(ensembleEstimatorRun, 1, _observationSize);
  free_uivector(ensembleEstimatorDen, 1, _observationSize);
  free_pdvector(ensembleEstimatorPtr, 1, sortedTimeInterestSize);
  if (_traceFlag & DL0_TRACE) {
    Rprintf("\nRSF:  Native code nominal exit.  ");
    now = updateTimeStamp(start);
  }
}  
