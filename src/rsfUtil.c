//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.2.2
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
#include   "node_ops.h"
#include   "rsfImpute.h"
#include   "rsfUtil.h"
extern uint getTraceFlag();
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
char restoreTree(uint    b,
                 Node   *parent,
                 uint   *leafCount,
                 uint   *offset,
                 uint   *treeID,
                 uint   *nodeID,
                 uint   *parmID,
                 double *spltPT) {
  uint i;
  char splitResult;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nrestoreTree() ENTRY ...\n");
  }
  if (b != treeID[*offset]) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid forest input record at line:  %10d", b);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    Rprintf("\nDiagnostic Trace of Tree Record:  \n");
    Rprintf("\n    treeID     nodeID     parmID       spltPT \n");
    Rprintf("%10d %10d %10d %12.4f \n", treeID[*offset], nodeID[*offset], parmID[*offset], spltPT[*offset]);
    exit(TRUE);
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n    treeID     nodeID     parmID       spltPT \n");
    Rprintf("%10d %10d %10d %12.4f \n", treeID[*offset], nodeID[*offset], parmID[*offset], spltPT[*offset]);
  }
  parent -> parent = parent;
  parent -> left  = NULL;
  parent -> right = NULL;
  for (i = 1; i <= _xSize; i++) {
    parent -> permissibleSplit[i] = FALSE;
  }
  parent -> splitFlag = FALSE;
  parent -> leafCount = nodeID[*offset];
  parent -> splitParameter = parmID[*offset];
  parent -> splitValue = spltPT[*offset];
  (*offset) ++;
  if ((parent -> splitParameter) != 0) {
    splitResult = TRUE;
    (*leafCount)++;
    parent -> left  = makeNode();
    ((parent -> left) -> leafCount) = (parent -> leafCount);
    restoreTree(b, 
                parent -> left, 
                leafCount, 
                offset, 
                treeID, 
                nodeID, 
                parmID, 
                spltPT);
    parent -> right = makeNode();
    ((parent -> right) -> leafCount) = *leafCount;
    restoreTree(b, 
                parent -> right, 
                leafCount, 
                offset, 
                treeID, 
                nodeID, 
                parmID, 
                spltPT);
  }
  else {
    splitResult = FALSE;
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nrestoreTree() EXIT ...\n");
  }
  return splitResult;
}
void saveTree(uint    b,
              Node   *parent,
              uint   *offset,
              uint   *treeID,
              uint   *nodeID,
              uint   *parmID,
              double *spltPT) {
  treeID[*offset] = b;
  nodeID[*offset] = parent -> leafCount;
  parmID[*offset] = parent -> splitParameter;
  spltPT[*offset] = parent -> splitValue;
  (*offset) ++;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    saveTree(b, parent ->  left, offset, treeID, nodeID, parmID, spltPT);
    saveTree(b, parent -> right, offset, treeID, nodeID, parmID, spltPT);
  }
}
char testNodeSize(Node *parent) {
  uint localDeathTimeSize;
  char result;
  uint i;
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ntestNodeSize() ENTRY ...\n");
  }
  localDeathTimeSize = 0;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      if (_status[_bootMembershipIndex[i]] == 1) {
        if (_mRecordMap[_bootMembershipIndex[i]] == 0) {
          localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        }
        else if (_mvSign[abs(CENS_IDX)][_mRecordMap[_bootMembershipIndex[i]]] == 0) {
          localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        }
      }
    }
  }
  for (i=1; i <= _masterTimeSize; i++) {
    if (localDeathTimeCount[i] > 0) {
      ++localDeathTimeSize;
    }
  }
  free_uivector(localDeathTimeCount, 1, _masterTimeSize);
  if (localDeathTimeSize >= _minimumDeathCount) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\ntestNodeSize() is:  %2d\n", result);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ntestNodeSize() EXIT ...\n");
  }
    return result;
}
Node* getMembership(Node *parent, double **predictor, uint index) {
  Node *result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    if ( (predictor[parent -> splitParameter][index]) <= (parent -> splitValue) ) {
      result = getMembership(parent ->  left, predictor, index);
    }
    else {
      result = getMembership(parent -> right, predictor, index);
    }
  }
  return result;
}
Node *getTerminalNode(uint leaf) {
  uint i, j;
  Node *parent;
  parent = NULL;
  for (j = 1; j <= _observationSize; j++) {
    if ((_nodeMembership[j] -> leafCount) == leaf) {
      parent = _nodeMembership[j];
      j = _observationSize;
    }
  }
  if (parent == NULL) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Proxy member for node %12d not found.", leaf);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    Rprintf("\nDiagnostic Trace of (individual, boot, node, leaf) vectors in data set:  ");
    Rprintf("\n        index         boot         node         leaf \n");
    for (i = 1; i <= _observationSize; i++) {
      Rprintf(" %12d %12d %12x %12d \n", i, _bootMembershipFlag[i], _nodeMembership[i], _nodeMembership[i] -> leafCount);
    }
    exit(TRUE);
  }
  return parent;
}
Node* randomizeMembership(Node    *parent, 
                          double **predictor, 
                          uint     individual, 
                          uint     splitParameter) {
  char randomSplitFlag;
  Node *result;
  result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    randomSplitFlag = FALSE;
    if (splitParameter > 0) {
      if ((parent -> splitParameter) == splitParameter) {
        randomSplitFlag = TRUE;
      }
    }
    else {
      if(_importanceFlag[parent -> splitParameter] == TRUE) {
        randomSplitFlag = TRUE;
      }
    }
    if(randomSplitFlag == TRUE) {
      if (ran2(_seed2Ptr) <= 0.5) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter);
      }
    }
    else {
      if ( (predictor[parent -> splitParameter][individual]) <= (parent -> splitValue) ) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter);
      }
    }
  }
  return result;
}
double getConcordanceIndex(int     polarity,
                           uint    size, 
                           double *statusPtr, 
                           double *timePtr, 
                           double *predictedOutcome,
                           uint   *oobCount) {
  uint i,j;
  ulong concordancePairSize;
  ulong concordanceWorseCount;
  double result;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetConcordanceIndex() ENTRY ...\n");
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nStatus and Time used in Concordance Index Calculations:  ");
    Rprintf("\n       Status         Time");
    for (i=1; i <= size; i++) {
      Rprintf("\n %12d %12.4f", (uint) statusPtr[i], timePtr[i]);
    }
    Rprintf("\n");
  }
  concordancePairSize = concordanceWorseCount = 0;
  for (i=1; i < size; i++) {
    for (j=i+1; j <= size; j++) {
      if (oobCount[i] != 0  && oobCount[j] != 0) {
        if ( (timePtr[i] > timePtr[j] && statusPtr[j] == 1) ||
             (timePtr[i] == timePtr[j] && statusPtr[j] == 1 && statusPtr[i] == 0) ) {
          concordancePairSize += 2;
          if ((polarity * predictedOutcome[j]) > (polarity * predictedOutcome[i])) {              
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[j] - predictedOutcome[i]) < EPSILON) {
            concordanceWorseCount += 1;
          }  
        }
        else if ( (timePtr[j] > timePtr[i] && statusPtr[i] == 1) ||
                  (timePtr[j] == timePtr[i] && statusPtr[i] == 1 && statusPtr[j] == 0) ) {
          concordancePairSize += 2;
          if ((polarity * predictedOutcome[i]) > (polarity * predictedOutcome[j])) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 1;
          }  
        }
        else if ( timePtr[i] == timePtr[j] && statusPtr[i] == 1 && statusPtr[j] == 1) {
          concordancePairSize += 2;
          if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
           concordanceWorseCount += 1;
          }
          else {
            concordanceWorseCount += 2;
          }
        }
      }  
    }  
  }  
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nConcordance pair and error update complete:");  
    Rprintf("\nCount of concordance pairs:         %20d", concordancePairSize);
    Rprintf("\nCount of pairs with worse outcome:  %20d", concordanceWorseCount);
    Rprintf("\n");
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetConcordanceIndex() EXIT() ...\n");
  }
  if (concordancePairSize == 0) {
    result = NA_REAL;
  }
  else {
    result = ((double) concordanceWorseCount / (double) concordancePairSize);
  }
  return result;
}
void getNelsonAalenEstimate(double **nelsonAalen,
                            uint treeID,
                            uint sortedTimeInterestSize) {
  uint leaf, i, j;
  Node *parent;
  uint priorTimePointIndex, currentTimePointIndex;
  double estimate;
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  getNelsonAalenEstimate() ENTRY ...\n");  
  }
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  for (leaf=1; leaf <= _leafCount_[treeID]; leaf++) {
    for (i=1; i <= _masterTimeSize; i++) {
      nodeParentAtRisk[i] = nodeParentDeath[i] = 0;
    }
    parent = getTerminalNode(leaf);
    for (i=1; i <= _observationSize; i++) {
      if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
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
        if (_masterTime[i] <= _timeInterest[j]) {
          currentTimePointIndex = i;
        }
        else {
          i = _masterTimeSize;
        }
      }
      if (getTraceFlag() & DL3_TRACE) {
        Rprintf("\nNelson-Aalen at risk and death counts for node:  %10d \n", parent -> leafCount);
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
      nelsonAalen[j][parent -> leafCount] = estimate;
      priorTimePointIndex = currentTimePointIndex;
    }  
    for (j=2; j <= sortedTimeInterestSize; j++) {
      nelsonAalen[j][parent -> leafCount] += nelsonAalen[j-1][parent -> leafCount];
    }
  }  
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nTree specific Nelson-Aalen estimator matrix:  %10d \n", treeID);
    Rprintf("          ");
    for (i=1; i <= _leafCount_[treeID]; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n");
    for (j=1; j <= sortedTimeInterestSize; j++) {
      Rprintf("%10d", j);
      for (i=1; i <= _leafCount_[treeID]; i++) {
        Rprintf("%10.4f", nelsonAalen[j][i]);
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  getNelsonAalenEstimate() EXIT ...\n");  
  }
}
void updateEnsembleCHF(uint mode, 
                       uint sortedTimeInterestSize,
                       uint treeID,
                       double **cumulativeHazard) {
  uint i, j, k;
  uint obsSize;
  uint *genericEnsembleDenPtr;
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  updateEnsembleCHF() ENTRY ...\n");  
  }
  switch (mode) {
  case RSF_GROW:
    for (i=1; i <= _observationSize; i++) {
      k = _nodeMembership[i] -> leafCount;
      if (_oobSampleSize[treeID] > 0) {
        if ( _bootMembershipFlag[i] == FALSE ) {
          for (j=1; j <= sortedTimeInterestSize; j++) {
            _oobEnsemblePtr[j][i] += cumulativeHazard[j][k];
          }
          _oobEnsembleDen[i] ++;
        }
      }
      _ensembleRun[i] = 0.0;
      for (j=1; j <= sortedTimeInterestSize; j++) {
        _fullEnsemblePtr[j][i] += cumulativeHazard[j][k];
        _ensembleRun[i] += _oobEnsemblePtr[j][i];
      }
      _fullEnsembleDen[i] ++;
      if (_oobEnsembleDen[i] != 0) {
        _ensembleRun[i] = _ensembleRun[i] / _oobEnsembleDen[i];
      }
    }
    break;
  case RSF_PRED:
    for (i=1; i <= _fobservationSize; i++) {
      k = _fnodeMembership[i] -> leafCount;
      _ensembleRun[i] = 0.0;
      for (j=1; j <= sortedTimeInterestSize; j++) {
        _fullEnsemblePtr[j][i] += cumulativeHazard[j][k];
        _ensembleRun[i] += _fullEnsemblePtr[j][i];
      }
      _fullEnsembleDen[i] ++;
      if (_fullEnsembleDen[i] != 0) {
        _ensembleRun[i] = _ensembleRun[i] / _fullEnsembleDen[i];
      }
    }
    break;
  case RSF_INTR:
    for (i=1; i <= _fobservationSize; i++) {
      k = _nodeMembership[_intrObservation[i]] -> leafCount;
      if (_foobSampleSize[treeID] > 0) {
        if ( _bootMembershipFlag[_intrObservation[i]] == FALSE ) {
          for (j=1; j <= sortedTimeInterestSize; j++) {
            _oobEnsemblePtr[j][i] += cumulativeHazard[j][k];
          }
          _oobEnsembleDen[i] ++;
        }
      }
      _ensembleRun[i] = 0.0;
      for (j=1; j <= sortedTimeInterestSize; j++) {
        _ensembleRun[i] += _oobEnsemblePtr[j][i];
      }
      if (_oobEnsembleDen[i] != 0) {
        _ensembleRun[i] = _ensembleRun[i] / _oobEnsembleDen[i];
      }
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
  if (getTraceFlag() & DL2_TRACE) {
    switch (mode) {
    case RSF_GROW:
      obsSize = _observationSize;
      genericEnsembleDenPtr = _oobEnsembleDen;
      break;
    case RSF_PRED:
      obsSize = _fobservationSize;
      genericEnsembleDenPtr = _fullEnsembleDen;
      break;
    case RSF_INTR:
      obsSize = _fobservationSize;
      genericEnsembleDenPtr = _oobEnsembleDen;
      break;
    default:
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Unknown case in switch encountered. ");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
      break;
    }
    if (_opt & OPT_OENS) {
      Rprintf("\nOOB Ensemble Estimator Numerator calculation: \n");
      Rprintf("          ");
      for (i=1; i <= _observationSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (j=1; j <= sortedTimeInterestSize; j++) {
        Rprintf("%10d", j);
        for (i=1; i <= _observationSize; i++) {
          Rprintf("%10.4f", _oobEnsemblePtr[j][i]);
        }
        Rprintf("\n");
      }
    }
    if (_opt & OPT_FENS) {
      Rprintf("\nFull Ensemble Estimator Numerator calculation: \n");
      Rprintf("          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (j=1; j <= sortedTimeInterestSize; j++) {
        Rprintf("%10d", j);
        for (i=1; i <= obsSize; i++) {
          Rprintf("%10.4f", _fullEnsemblePtr[j][i]);
        }
        Rprintf("\n");
      }
    }
    if ((_opt & OPT_OENS) || (_opt & OPT_FENS)) {
      Rprintf("\nRunning Ensemble Estimator:  \n");
      Rprintf("          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10.4f", _ensembleRun[i]);
      }
      Rprintf("\n");
      Rprintf("\nEnsemble Estimator Denominator calculation: \n");
      Rprintf("          ");  
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", genericEnsembleDenPtr[i]);
      }
      Rprintf("\n");
    }
    if (getTraceFlag() & RSF_GROW) {
      Rprintf("\nClassification of OOB elements:  ");
      Rprintf("\n       Leaf     Indiv            EE      predictors ->  ");
      for (i=1; i <= _leafCount_[treeID]; i++) {
        for (j = 1; j <= _observationSize; j++) {
          if (_bootMembershipFlag[j] == FALSE) {
            if ((_nodeMembership[j] -> leafCount) == i) {
              Rprintf("\n %10d %10d %12.4f", i, j, _ensembleRun[j]);
              for (k = 1; k <= _xSize; k++) {
                Rprintf("%12.4f", _observation[k][j]);
              }
            }
          }
        }
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  updateEnsembleCHF() EXIT ...\n");  
  }
}
void getMeanSurvivalTime(double *meanSurvivalTime,
                         uint treeID) {
  uint leaf, i;
  double totalDeathTime;
  uint deathCount;
  Node *parent;
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  getMeanSurvivalTime() ENTRY ...\n");  
  }
  for (leaf=1; leaf <= _leafCount_[treeID]; leaf++) {
    totalDeathTime = 0;
    deathCount = 0;
    parent = getTerminalNode(leaf);
    for (i=1; i <= _observationSize; i++) {
      if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
        if (_status[_bootMembershipIndex[i]] == 1) {
          totalDeathTime += _masterTime[_masterTimeIndex[_bootMembershipIndex[i]]];
          deathCount ++;
        }
      }
    }
    if (deathCount > 0) {
      meanSurvivalTime[leaf] = totalDeathTime / deathCount;
    }
    else {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Zero death count encountered in node:  %10d", leaf);
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
    }
  }  
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nTree specific mean survival time vector:  %10d \n", treeID);
    Rprintf("          ");
    for (i=1; i <= _leafCount_[treeID]; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n");
    Rprintf("          ");
    for (i=1; i <= _leafCount_[treeID]; i++) {
      Rprintf("%10.4f", meanSurvivalTime[i]);
    }
    Rprintf("\n");
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  getMeanSurvivalTime() EXIT ...\n");  
  }
}
void updateEnsembleSurvivalTime(uint mode, 
                                uint treeID, 
                                double *meanSurvivalTime) {
  uint i, j, k;
  uint obsSize;
  uint *genericEnsembleDenPtr;
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  updateEnsembleSurvivalTime() ENTRY ...\n");  
  }
  switch (mode) {
  case RSF_GROW:
    for (i=1; i <= _observationSize; i++) {
      k = _nodeMembership[i] -> leafCount;
      if (_oobSampleSize[treeID] > 0) {
        if ( _bootMembershipFlag[i] == FALSE ) {
          _oobEnsemblePtr[1][i] += meanSurvivalTime[k];
          _oobEnsembleDen[i] ++;
        }
      }
      _fullEnsemblePtr[1][i] += meanSurvivalTime[k];
      _fullEnsembleDen[i] ++;
      if (_oobEnsembleDen[i] != 0) {
        _ensembleRun[i] = _oobEnsemblePtr[1][i];
        _ensembleRun[i] = _ensembleRun[i] / _oobEnsembleDen[i];
      }
      else {
        _ensembleRun[i] = 0;
      }
    }
    break;
  case RSF_PRED:
    for (i=1; i <= _fobservationSize; i++) {
      k = _fnodeMembership[i] -> leafCount;
      _fullEnsemblePtr[1][i] += meanSurvivalTime[k];
      _fullEnsembleDen[i] ++;
      if (_fullEnsembleDen[i] != 0) {
        _ensembleRun[i] = _fullEnsemblePtr[1][i];
        _ensembleRun[i] = _ensembleRun[i] / _fullEnsembleDen[i];
      }
      else {
        _ensembleRun[i] = 0;
      }
    }
    break;
  case RSF_INTR:
    for (i=1; i <= _fobservationSize; i++) {
      k = _nodeMembership[_intrObservation[i]] -> leafCount;
      if (_foobSampleSize[treeID] > 0) {
        if ( _bootMembershipFlag[_intrObservation[i]] == FALSE ) {
          _oobEnsemblePtr[1][i] += meanSurvivalTime[k];
          _oobEnsembleDen[i] ++;
        }
      }
      if (_oobEnsembleDen[i] != 0) {
        _ensembleRun[i] = _oobEnsemblePtr[1][i];
        _ensembleRun[i] = _ensembleRun[i] / _oobEnsembleDen[i];
      }
      else {
        _ensembleRun[i] = 0;
      }
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
  if (getTraceFlag() & DL2_TRACE) {
    switch (mode) {
    case RSF_GROW:
      obsSize = _observationSize;
      genericEnsembleDenPtr = _oobEnsembleDen;
      break;
    case RSF_PRED:
      obsSize = _fobservationSize;
      genericEnsembleDenPtr = _fullEnsembleDen;
      break;
    case RSF_INTR:
      obsSize = _fobservationSize;
      genericEnsembleDenPtr = _oobEnsembleDen;
      break;
    default:
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Unknown case in switch encountered. ");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
      break;
    }
    if (_opt & OPT_OENS) {
      Rprintf("\nOOB Ensemble Estimator Numerator calculation: \n");
      Rprintf("          ");
      for (i=1; i <= _observationSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (j=1; j <= 1; j++) {
        Rprintf("%10d", j);
        for (i=1; i <= _observationSize; i++) {
          Rprintf("%10.4f", _oobEnsemblePtr[j][i]);
        }
        Rprintf("\n");
      }
    }
    if (_opt & OPT_FENS) {
      Rprintf("\nFull Ensemble Estimator Numerator calculation: \n");
      Rprintf("          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (j=1; j <= 1; j++) {
        Rprintf("%10d", j);
        for (i=1; i <= obsSize; i++) {
          Rprintf("%10.4f", _fullEnsemblePtr[j][i]);
        }
        Rprintf("\n");
      }
    }
    if ((_opt & OPT_OENS) || (_opt & OPT_FENS)) {
      Rprintf("\nRunning Ensemble Estimator:  \n");
      Rprintf("          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10.4f", _ensembleRun[i]);
      }
      Rprintf("\n");
      Rprintf("\nEnsemble Estimator Denominator calculation: \n");
      Rprintf("          ");  
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", genericEnsembleDenPtr[i]);
      }
      Rprintf("\n");
    }
    if (getTraceFlag() & RSF_GROW) {
      Rprintf("\nClassification of OOB elements:  ");
      Rprintf("\n       Leaf     Indiv            EE      predictors ->  ");
      for (i=1; i <= _leafCount_[treeID]; i++) {
        for (j = 1; j <= _observationSize; j++) {
          if (_bootMembershipFlag[j] == FALSE) {
            if ((_nodeMembership[j] -> leafCount) == i) {
              Rprintf("\n %10d %10d %12.4f", i, j, _ensembleRun[j]);
              for (k = 1; k <= _xSize; k++) {
                Rprintf("%12.4f", _observation[k][j]);
              }
            }
          }
        }
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  updateEnsembleSurvivalTime() EXIT ...\n");  
  }
}
void getVariableImportance (char     mode,
                            uint     sortedTimeInterestSize,
                            uint     leafCount,
                            double **cumulativeHazard,
                            double  *meanSurvivalTime,
                            Node    *rootPtr,
                            uint     b) {
  uint i,j,k,p;
  uint obsSize;
  uint varSize;
  uint permuteObsSize;
  uint permuteVarSize;
  uint   *indexVIMP;
  uint   *permuteVIMP;
  double **originalVIMP;
  double **predictorPtr;
  char result;
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  getVariableImportance() ENTRY ...\n");  
  }
  if (!(_opt & OPT_VIMP)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to compute variable importance though not requested.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  obsSize = 0;          
  varSize = 0;          
  permuteObsSize = 0;   
  permuteVarSize = 0;   
  predictorPtr = NULL;  
  result = TRUE;
  switch (mode) {
  case RSF_GROW:
    predictorPtr = _observation;
    if (_oobSampleSize[b] == 0) {
      result = FALSE;
    }
    break;
  case RSF_PRED:
    predictorPtr = _fobservation;
    break;
  case RSF_INTR:
    predictorPtr = _fobservation;
    if (_foobSampleSize[b] == 0) {
      result = FALSE;
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
  if (result == TRUE) {    
    if (_opt & (~OPT_VIMP) & (~OPT_VIMP_JOIN) & OPT_VIMP_TYPE) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nType is VIMP_RAND.");  
        Rprintf("\nRandom seed for ran2():  %20d", *_seed2Ptr);
      }
      switch (mode) {
      case RSF_GROW:
        for (p=1; p <= _xSize; p++) {
          for (i=1; i <= _observationSize; i++) {
            if ( _bootMembershipFlag[i] == FALSE ) {
              k = randomizeMembership(rootPtr, predictorPtr, i, p) -> leafCount;
              if (_opt & (OPT_POUT_TYPE)) {
                _vimpEnsembleRun[p][i] += meanSurvivalTime[k];
              }
              else {
                for (j=1; j <= sortedTimeInterestSize; j++) {
                  _vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
                }
              }
            }
          }
        }
        break;
      case RSF_PRED:
        for (p=1; p <= _xSize; p++) {
          for (i=1; i <= _fobservationSize; i++) {
            k = randomizeMembership(rootPtr, predictorPtr, i, p) -> leafCount;
            if (_opt & (OPT_POUT_TYPE)) {
              _vimpEnsembleRun[p][i] += meanSurvivalTime[k];
            }
            else {
              for (j=1; j <= sortedTimeInterestSize; j++) {
                _vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
              }
            }
          }
        }
        break;
      case RSF_INTR:
        if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
          for (i=1; i <= _fobservationSize; i++) {
            if ( _bootMembershipFlag[_intrObservation[i]] == FALSE ) {
              k = randomizeMembership(rootPtr, predictorPtr, i, 0) -> leafCount;
              if (_opt & (OPT_POUT_TYPE)) {
                _vimpEnsembleRun[1][i] += meanSurvivalTime[k];
              }
              else {
                for (j=1; j <= sortedTimeInterestSize; j++) {
                  _vimpEnsembleRun[1][i] += cumulativeHazard[j][k];
                }
              }
            }
          }
        }  
        else {
          for (p=1; p <= _intrPredictorSize; p++) {
            for (i=1; i <= _fobservationSize; i++) {
              if ( _bootMembershipFlag[_intrObservation[i]] == FALSE ) {
                k = randomizeMembership(rootPtr, predictorPtr, i, _intrPredictor[p]) -> leafCount;
                if (_opt & (OPT_POUT_TYPE)) {
                  _vimpEnsembleRun[p][i] += meanSurvivalTime[k];
                }
                else {
                  for (j=1; j <= sortedTimeInterestSize; j++) {
                    _vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
                  }
                }
              }
            }
          }
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
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  VIMP random split calculation complete.");  
      }
    }  
    else if (!(_opt & (~OPT_VIMP) & (~OPT_VIMP_JOIN) & OPT_VIMP_TYPE)) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nType is VIMP_PERM.");  
        Rprintf("\nRandom seed for ran2():  %20d", *_seed2Ptr);
      }
      switch (mode) {
      case RSF_GROW:
        permuteObsSize = _oobSampleSize[b];
        permuteVarSize = 1;
        break;
      case RSF_PRED:
        permuteObsSize = _fobservationSize;
        permuteVarSize = 1;
        break;
      case RSF_INTR:
        permuteObsSize = _foobSampleSize[b];
        if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
          permuteVarSize = _intrPredictorSize;
        }
        else {
          permuteVarSize = 1;
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
      indexVIMP = uivector(1, permuteObsSize);
      permuteVIMP = uivector(1, permuteObsSize);
      originalVIMP = dmatrix(1, permuteVarSize, 1, permuteObsSize);
      k = 0;
      switch (mode) {
      case RSF_GROW:
        for (i=1; i <= _observationSize; i++) {
          if ( _bootMembershipFlag[i] == FALSE ) {
            k++;
            indexVIMP[k] = i;
          }
        }
        break;
      case RSF_PRED:
        for (i=1; i <= _fobservationSize; i++) {
          k++;
          indexVIMP[k] = i;
        }
        break;
      case RSF_INTR:
        for (i=1; i <= _fobservationSize; i++) {
          if ( _bootMembershipFlag[_intrObservation[i]] == FALSE ) {
            k++;
            indexVIMP[k] = i;
          }
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
      if (k != permuteObsSize) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  VIMP candidate selection failed.");
        Rprintf("\nRSF:  %10d available, %10d selected.", permuteObsSize, k);
        Rprintf("\nRSF:  Please Contact Technical Support.");
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
      }
      switch (mode) {
      case RSF_GROW:
        for (p=1; p <= _xSize; p++) {
          for (k=1; k<= permuteObsSize; k++) {
            originalVIMP[1][k] = predictorPtr[p][indexVIMP[k]];
          }
          permute(permuteObsSize, permuteVIMP);
          for (k=1; k <= permuteObsSize; k++) {
            predictorPtr[p][indexVIMP[k]] = originalVIMP[1][permuteVIMP[k]];
          }
          for (i=1; i <= _observationSize; i++) {
            if ( _bootMembershipFlag[i] == FALSE ) {
              k = getMembership(rootPtr, predictorPtr, i) -> leafCount;
              if (_opt & (OPT_POUT_TYPE)) {
                _vimpEnsembleRun[p][i] += meanSurvivalTime[k];
              }
              else {
                for (j=1; j <= sortedTimeInterestSize; j++) {
                  _vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
                }
              }
            }
          }
          for (k=1; k <= permuteObsSize; k++) {
            predictorPtr[p][indexVIMP[k]] = originalVIMP[1][k];
          }
        }
        break;
      case RSF_PRED:
        for (p=1; p <= _xSize; p++) {
          for (k=1; k<= permuteObsSize; k++) {
            originalVIMP[1][k] = predictorPtr[p][indexVIMP[k]];
          }
          permute(permuteObsSize, permuteVIMP);
          for (k=1; k <= permuteObsSize; k++) {
            predictorPtr[p][indexVIMP[k]] = originalVIMP[1][permuteVIMP[k]];
          }
          for (i=1; i <= _fobservationSize; i++) {
            k = getMembership(rootPtr, predictorPtr, i) -> leafCount;
            if (_opt & (OPT_POUT_TYPE)) {
              _vimpEnsembleRun[p][i] += meanSurvivalTime[k];
            }
            else {
              for (j=1; j <= sortedTimeInterestSize; j++) {
                _vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
              }
            }
          }
          for (k=1; k <= permuteObsSize; k++) {
            predictorPtr[p][indexVIMP[k]] = originalVIMP[1][k];
          }
        }
        break;
      case RSF_INTR:
        if (_opt & (~OPT_VIMP) & OPT_VIMP_JOIN) {
          for (p=1; p <= permuteVarSize; p++) {
            for (k=1; k<= permuteObsSize; k++) {
              originalVIMP[p][k] = predictorPtr[_intrPredictor[p]][indexVIMP[k]];
            }
            permute(permuteObsSize, permuteVIMP);
            for (k=1; k <= permuteObsSize; k++) {
              predictorPtr[_intrPredictor[p]][indexVIMP[k]] = originalVIMP[p][permuteVIMP[k]];
            }
          }
          for (i=1; i <= _fobservationSize; i++) {
            if ( _bootMembershipFlag[_intrObservation[i]] == FALSE ) {
              k = getMembership(rootPtr, predictorPtr, i) -> leafCount;
              if (_opt & (OPT_POUT_TYPE)) {
                _vimpEnsembleRun[1][i] += meanSurvivalTime[k];
              }
              else {
                for (j=1; j <= sortedTimeInterestSize; j++) {
                  _vimpEnsembleRun[1][i] += cumulativeHazard[j][k];
                }
              }
            }
          }
          for (p=1; p <= permuteVarSize; p++) {
            for (k=1; k <= permuteObsSize; k++) {
              predictorPtr[_intrPredictor[p]][indexVIMP[k]] = originalVIMP[p][k];
            }
          }
        }  
        else {
          for (p=1; p <= _intrPredictorSize; p++) {
            for (k=1; k<= permuteObsSize; k++) {
              originalVIMP[1][k] = predictorPtr[_intrPredictor[p]][indexVIMP[k]];
            }
            permute(permuteObsSize, permuteVIMP);
            for (k=1; k <= permuteObsSize; k++) {
              predictorPtr[_intrPredictor[p]][indexVIMP[k]] = originalVIMP[1][permuteVIMP[k]];
            }
            for (i=1; i <= _fobservationSize; i++) {
              if ( _bootMembershipFlag[_intrObservation[i]] == FALSE ) {
                k = getMembership(rootPtr, predictorPtr, i) -> leafCount;
                if (_opt & (OPT_POUT_TYPE)) {
                  _vimpEnsembleRun[p][i] += meanSurvivalTime[k];
                }
                else {
                  for (j=1; j <= sortedTimeInterestSize; j++) {
                    _vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
                  }
                }
              }
            }
            for (k=1; k <= permuteObsSize; k++) {
              predictorPtr[_intrPredictor[p]][indexVIMP[k]] = originalVIMP[1][k];
            }
          }
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
      free_uivector(indexVIMP, 1, permuteObsSize);
      free_uivector(permuteVIMP, 1, permuteObsSize);
      free_dmatrix(originalVIMP, 1, permuteVarSize, 1, permuteObsSize);
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  VIMP permutation calculation complete.");  
      }
    }  
    else {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Unknown VIMP perturbation type encountered. ");
      Rprintf("\nRSF:  Option flag is:  %10x", _opt);
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
    } 
    if (getTraceFlag() & DL2_TRACE) {
      switch (mode) {
      case RSF_GROW:
        obsSize = _observationSize;
        varSize = _xSize;
        break;
      case RSF_PRED:
        obsSize = _fobservationSize;
        varSize = _xSize;
        break;
      case RSF_INTR:
        obsSize = _fobservationSize;
        varSize = 1;
        break;
      default:
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Unknown case in switch encountered. ");
        Rprintf("\nRSF:  Please Contact Technical Support.");
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
        break;
      }
      Rprintf("\nTree specific variable importance calculation:  %10d \n", b);
      Rprintf("          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (p=1; p <= varSize; p++) {
        Rprintf("%10d", p);
        for (i=1; i <= obsSize; i++) {
          Rprintf("%10.4f", _vimpEnsembleRun[p][i]);
        }
        Rprintf("\n");
      }
    }
  }  
  else {
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nVIMP omitted since OOB sample size is zero. \n");      
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  getVariableImportance() EXIT ...\n");  
  }
}
void getVariablesUsed(Node *parent, uint *varUsedPtr) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetVariablesUsed() ENTRY ...\n");  
  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    varUsedPtr[parent -> splitParameter] ++;
     getVariablesUsed(parent ->  left, varUsedPtr);
     getVariablesUsed(parent -> right, varUsedPtr);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetVariablesUsed() EXIT ...\n");  
  }
  return;
}
void updateEnsembleEvents (char      multipleImputeFlag,
                           uint      mode,
                           uint      sortedTimeInterestSize,
                           Node     *rootPtr,
                           uint      b,
                           char    **dmRecordBootFlag,
                           double ***dmvImputation) {
  uint i;
  uint obsSize;
  double **cumulativeHazard;
  double  *meanSurvivalTime;
  uint    *genericEnsembleDenPtr;
  double *statusPtr, *orgStatusPtr, *mStatusPtr;
  double *timePtr, *orgTimePtr, *mTimePtr;
  uint *varUsedRow;
  char concordanceImputeFlag;
  int concordancePolarity;
  double concordanceIndex;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nupdateEnsembleEvents() ENTRY ...\n");
  }
  cumulativeHazard = NULL;  
  meanSurvivalTime = NULL;  
  statusPtr   = NULL;  
  timePtr     = NULL;  
  mStatusPtr  = NULL;  
  mTimePtr    = NULL;  
  genericEnsembleDenPtr = NULL;  
  obsSize = 0;  
  if (_leafCount_[b] == 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to compute performance on a rejected tree:  %10d", b);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (_opt & (OPT_POUT_TYPE)) {
    meanSurvivalTime = dvector(1, _leafCount_[b]);
    concordancePolarity = -1;
    getMeanSurvivalTime(meanSurvivalTime, b);
    updateEnsembleSurvivalTime(mode, b, meanSurvivalTime);
  }
  else {
    cumulativeHazard = dmatrix(1, sortedTimeInterestSize, 1, _leafCount_[b]);
    concordancePolarity = 1;
    getNelsonAalenEstimate(cumulativeHazard, b, sortedTimeInterestSize);
    updateEnsembleCHF(mode, sortedTimeInterestSize, b, cumulativeHazard);
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Predicted outcome calculation complete.");  
  }
  if (_opt & OPT_VIMP) {
    getVariableImportance(mode,
                          sortedTimeInterestSize,
                          _leafCount_[b],
                          cumulativeHazard,
                          meanSurvivalTime,
                          rootPtr,
                          b);
  }
  if (_opt & (OPT_POUT_TYPE)) {
    free_dvector(meanSurvivalTime, 1, _leafCount_[b]);
  }
  else {
    free_dmatrix(cumulativeHazard, 1, sortedTimeInterestSize, 1, _leafCount_[b]);
  }
  if (_opt & OPT_VUSE) {
    if (_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
      varUsedRow = _varUsedPtr[b];
    }
    else {
      varUsedRow = _varUsedPtr[1];
    }
    getVariablesUsed(rootPtr, varUsedRow);
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nVariables Used Counts by Parameter:  \n");
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12d", i);
      }
      Rprintf("\n");
      for (i=1; i <= _xSize; i++) {
        Rprintf(" %12d", varUsedRow[i]);
      }
      Rprintf("\n");
    }
  }
  if (_opt & OPT_PERF) {
    switch (mode) {
    case RSF_GROW:
      obsSize = _observationSize;
      genericEnsembleDenPtr = _oobEnsembleDen;
      break;
    case RSF_PRED:
      obsSize = _fobservationSize;
      genericEnsembleDenPtr = _fullEnsembleDen;
      break;
    case RSF_INTR:
      obsSize = _fobservationSize;
      genericEnsembleDenPtr = _oobEnsembleDen;
      break;
    default:
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Unknown case in switch encountered. ");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
      break;
    }
    concordanceImputeFlag = FALSE;
    if (mode == RSF_GROW) {
      statusPtr = orgStatusPtr = _status;
      timePtr = orgTimePtr = _time;
      if (multipleImputeFlag == FALSE) {
        if (_mRecordSize > 0) {
          concordanceImputeFlag = TRUE;
        }
      }
    } 
    else {
      statusPtr = orgStatusPtr = _fstatus;
      timePtr = orgTimePtr = _ftime;
      if (_fmRecordSize > 0) {
        concordanceImputeFlag = TRUE;
      }
    }  
    if (concordanceImputeFlag == TRUE) {
      mStatusPtr = dvector(1, obsSize);
      mTimePtr   = dvector(1, obsSize);
      for (i=1; i <= obsSize; i++) {
        mStatusPtr[i] = orgStatusPtr[i];
        mTimePtr[i]   = orgTimePtr[i];
      }
      imputeConcordance(mode,
                        b,
                        dmRecordBootFlag,
                        dmvImputation,
                        mStatusPtr,
                        mTimePtr);
      statusPtr = mStatusPtr;
      timePtr   = mTimePtr;
    }  
    concordanceIndex = getConcordanceIndex(concordancePolarity,
                                           obsSize, 
                                           statusPtr, 
                                           timePtr, 
                                           _ensembleRun, 
                                           genericEnsembleDenPtr);
    if (ISNA(concordanceIndex)) {
      _performance_[b] = NA_REAL;
    }
    else {
      _performance_[b] = 1.0 - concordanceIndex;
    }
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Error Rate:      %20.4f", _performance_[b]);
      Rprintf("\n");
    }
    if (concordanceImputeFlag > 0) {
      free_dvector(mStatusPtr, 1, obsSize);
      free_dvector(mTimePtr, 1, obsSize);
    }  
  }  
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nupdateEnsembleEvents() EXIT ...\n");
  }
}
