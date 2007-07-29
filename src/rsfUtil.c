//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.0.0
//
//  Copyright 2007, Cleveland Clinic
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
  char splitResult;
  uint i;
  if (getTraceFlag() & DL1_TRACE) {
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
    parent -> permissibleSplit[i][1] = 0;
    parent -> permissibleSplit[i][2] = 0;
  }
  parent -> splitValueIndex = 0;
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
double getConcordanceIndex(uint size, 
                           double *statusPtr, 
                           double *timePtr, 
                           double *mortality,
                           uint *oobCount
                          ) {
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
      Rprintf("\n %12d %12.4f", (uint) round(statusPtr[i]), timePtr[i]);
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
          if (mortality[j] > mortality[i]) {              
            concordanceWorseCount += 2;
          }
          else if (fabs(mortality[j] - mortality[i]) < EPSILON) {
            concordanceWorseCount += 1;
          }  
        }
        else if ( (timePtr[j] > timePtr[i] && statusPtr[i] == 1) ||
                  (timePtr[j] == timePtr[i] && statusPtr[i] == 1 && statusPtr[j] == 0) ) {
          concordancePairSize += 2;
          if (mortality[i] > mortality[j]) {
            concordanceWorseCount += 2;
          }
          else if (fabs(mortality[i] - mortality[j]) < EPSILON) {
            concordanceWorseCount += 1;
          }  
        }
        else if ( timePtr[i] == timePtr[j] && statusPtr[i] == 1 && statusPtr[j] == 1) {
          concordancePairSize += 2;
          if (fabs(mortality[i] - mortality[j]) < EPSILON) {
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
void getCumulativeHazardEstimate(double **cumulativeHazard,
                                 Node *parent,
                                 uint sortedTimeInterestSize) {
  uint i, j;
  uint priorTimePointIndex, currentTimePointIndex;
  double estimate;
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  for (i=1; i <= _masterTimeSize; i++) {
    nodeParentAtRisk[i] = nodeParentDeath[i] = 0;
  }
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
double getPerformance (uint      mode,
                       uint      sortedTimeInterestSize,
                       uint      leafCount,
                       double  **oobEnsemblePtr,
                       double  **fullEnsemblePtr,
                       double   *ensembleRun,
                       uint     *ensembleDen,
                       uint      oobSampleSize,
                       Node     *rootPtr,
                       double  **vimpEnsembleRun,
                       uint      b,
                       char    **mRecordBootFlag,
                       double ***mvImputation
) {
  uint i,j,k,p,r;
  uint obsSize;
  uint permuteSize;
  double **cumulativeHazard;
  double **genericEnsemblePtr;
  char   *randomVIMP;
  uint   *indexVIMP;
  uint   *permuteVIMP;
  double *originalVIMP;
  double *statusPtr, *orgStatusPtr, *mStatusPtr;
  double *timePtr, *orgTimePtr, *mTimePtr;
  uint    *mRecordMap;
  uint     mRecordSize;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  double **predictorPtr;
  double *orgOutcomePtr, *newOutcomePtr, *mOutcomePtr;
  char outcomeFlag, mPredictorFlag;
  char imputeFlag;
  uint unsignedIndex;
  double concordanceIndex;
  double performance;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetPerformance() ENTRY ...\n");
  }
  mStatusPtr = NULL;  
  mTimePtr   = NULL;  
  mRecordSize = NA_REAL;  
  mvSize      = 0;        
  mvIndex     = NULL;     
  if (leafCount == 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to compute performance on a rejected tree:  %10d", b);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  performance = NA_REAL;
  cumulativeHazard = dmatrix(1, sortedTimeInterestSize, 1, leafCount);
  if (mode == RSF_GROW) {
    obsSize = _observationSize;
  }
  else {
    obsSize = _fobservationSize;
  }
  for (i=1; i <= leafCount; i++) {
    k = 0;
    for (j = 1; j <= _observationSize; j++) {
      if (_bootMembershipFlag[j] == TRUE) {
        if ((_nodeMembership[j] -> leafCount) == i) {
          k = j;
          j = _observationSize;
        }
      }
    }
    if (k == 0) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  CHE Node membership not found for leaf:  %10d", i);
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      Rprintf("\nDiagnostic Trace of (individual, boot, node, leaf) vectors in data set:  ");
      Rprintf("\n        index         boot         node         leaf \n");
      for (r = 1; r <= _observationSize; r++) {
        Rprintf(" %12d %12d %12x %12d \n", r, _bootMembershipFlag[r], _nodeMembership[r], _nodeMembership[r] -> leafCount);
      }
      exit(TRUE);
    }
    getCumulativeHazardEstimate(cumulativeHazard,
                                _nodeMembership[k],
                                sortedTimeInterestSize);
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  CHE calculation complete.");  
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nTree specific cumulative hazard calculation: \n");
    Rprintf("          ");
    for (i=1; i <= leafCount; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n");
    for (j=1; j <= sortedTimeInterestSize; j++) {
      Rprintf("%10d", j);
      for (i=1; i <= leafCount; i++) {
        Rprintf("%10.4f", cumulativeHazard[j][i]);
      }
      Rprintf("\n");
    }
  }
  if (mode == RSF_GROW) {
    for (i=1; i <= _observationSize; i++) {
      k = _nodeMembership[i] -> leafCount;
      if ( _bootMembershipFlag[i] == FALSE ) {
        for (j=1; j <= sortedTimeInterestSize; j++) {
          oobEnsemblePtr[j][i] += cumulativeHazard[j][k];
        }
        ensembleDen[i] ++;
      }
      for (j=1; j <= sortedTimeInterestSize; j++) {
        fullEnsemblePtr[j][i] += cumulativeHazard[j][k];
      }
    }
  }
  else {
    for (i=1; i <= _fobservationSize; i++) {
      k = _fnodeMembership[i] -> leafCount;
      for (j=1; j <= sortedTimeInterestSize; j++) {
        fullEnsemblePtr[j][i] += cumulativeHazard[j][k];
      }
      ensembleDen[i] ++;
    }
  }
  if ((_mup & MUP_VIMP) && (oobSampleSize > 0)) {
    if (mode == RSF_GROW) {
      permuteSize = oobSampleSize;
      predictorPtr = _observation;
    }
    else {
      permuteSize = _fobservationSize;
      predictorPtr = _fobservation;
    }
    indexVIMP = uivector(1, permuteSize);
    randomVIMP = cvector(1, permuteSize);
    permuteVIMP = uivector(1, permuteSize);
    originalVIMP = dvector(1, permuteSize);
    k = 0;
    if (mode == RSF_GROW) {
      for (i=1; i <= _observationSize; i++) {
        if ( _bootMembershipFlag[i] == FALSE ) {
          k++;
          indexVIMP[k] = i;
        }
      }
    }
    else {
      for (i=1; i <= _fobservationSize; i++) {
        k++;
        indexVIMP[k] = i;
      }
    }
    if (k != permuteSize) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  VIMP candidate selection failed.");
      Rprintf("\nRSF:  %10d available, %10d selected.", permuteSize, k);
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
    }
    for (p=1; p <= _xSize; p++) {
      for (k=1; k<= permuteSize; k++) {
        originalVIMP[k] = predictorPtr[p][indexVIMP[k]];
        randomVIMP[k] = TRUE;
      }
      for (i=permuteSize; i > 0; i--) {
        k = ceil(ran2(_seed2Ptr)*(i*1.0));
        for (j = 1; k > 0; j++) {
          if (randomVIMP[j] == TRUE) {
            k--;
          }
        }
        randomVIMP[j-1] = ACTIVE;
        permuteVIMP[i] = j-1;
      }
      if (getTraceFlag() & DL3_TRACE) {
        Rprintf("\nCovariate Permutation:  %10d\n", p);
          Rprintf("     index   indexVIMP    permVIMP\n");
        for (i=1; i <= permuteSize; i++) {
          Rprintf("%10d %10d %10d \n", i, indexVIMP[i], permuteVIMP[i]);
        }
      }
      for (k=1; k <= permuteSize; k++) {
        predictorPtr[p][indexVIMP[k]] = originalVIMP[permuteVIMP[k]];
      }
      if (mode == RSF_GROW) {
        for (i=1; i <= _observationSize; i++) {
          if ( _bootMembershipFlag[i] == FALSE ) {
            k = getMembership(rootPtr, predictorPtr, i) -> leafCount;
            for (j=1; j <= sortedTimeInterestSize; j++) {
              vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
            }
          }
        }
      }
      else {
        for (i=1; i <= _fobservationSize; i++) {
          k = getMembership(rootPtr, predictorPtr, i) -> leafCount;
          for (j=1; j <= sortedTimeInterestSize; j++) {
            vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
          }
        }
      }
      for (k=1; k <= permuteSize; k++) {
        predictorPtr[p][indexVIMP[k]] = originalVIMP[k];
      }
    }  
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nTree specific variable importance calculation: \n");
      Rprintf("          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (p=1; p <= _xSize; p++) {
        Rprintf("%10d", p);
        for (i=1; i <= obsSize; i++) {
          Rprintf("%10.4f", vimpEnsembleRun[p][i]);
        }
        Rprintf("\n");
      }
    }
    free_uivector(indexVIMP, 1, permuteSize);
    free_cvector(randomVIMP, 1, permuteSize);
    free_uivector(permuteVIMP, 1, permuteSize);
    free_dvector(originalVIMP, 1, permuteSize);
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  VIMP calculation complete.");  
    }
  }  
  free_dmatrix(cumulativeHazard, 1, sortedTimeInterestSize, 1, leafCount);
  if ((getTraceFlag() & RSF_GROW) && (getTraceFlag() & DL2_TRACE)) {
    Rprintf("\nOOB Ensemble Estimator Numerator calculation: \n");
    Rprintf("          ");
    for (i=1; i <= _observationSize; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n");
    for (j=1; j <= sortedTimeInterestSize; j++) {
      Rprintf("%10d", j);
      for (i=1; i <= _observationSize; i++) {
        Rprintf("%10.4f", oobEnsemblePtr[j][i]);
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nFull Ensemble Estimator Numerator calculation: \n");
    Rprintf("          ");
    for (i=1; i <= obsSize; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n");
    for (j=1; j <= sortedTimeInterestSize; j++) {
      Rprintf("%10d", j);
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10.4f", fullEnsemblePtr[j][i]);
      }
      Rprintf("\n");
    }
  }
  if (mode == RSF_GROW) {
    genericEnsemblePtr = oobEnsemblePtr;
  }
  else {
    genericEnsemblePtr = fullEnsemblePtr;
  }
  for (i = 1; i <= obsSize; i++) {
    ensembleRun[i] = 0.0;
    for (j=1; j <= sortedTimeInterestSize; j++) {
      ensembleRun[i] += genericEnsemblePtr[j][i];
    }
    if (ensembleDen[i] != 0) {
      ensembleRun[i] = ensembleRun[i] / ensembleDen[i];
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Tree ensemble estimate complete.");  
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nEnsemble Estimator Denominator calculation: \n");
    Rprintf("          ");  
    for (i=1; i <= obsSize; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n          ");
    for (i=1; i <= obsSize; i++) {
      Rprintf("%10d", ensembleDen[i]);
    }
    Rprintf("\n");
    Rprintf("\nRunning Ensemble Estimator:  \n");
    Rprintf("          ");
    for (i=1; i <= obsSize; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n          ");
    for (i=1; i <= obsSize; i++) {
      Rprintf("%10.4f", ensembleRun[i]);
    }
    Rprintf("\n");
  }
  if ((getTraceFlag() & RSF_GROW) && (getTraceFlag() & DL2_TRACE)) {
    Rprintf("\nClassification of OOB elements:  ");
    Rprintf("\n       Leaf     Indiv            EE      predictors ->  ");
    for (i=1; i <= leafCount; i++) {
      for (j = 1; j <= _observationSize; j++) {
	if (_bootMembershipFlag[j] == FALSE) {
	  if ((_nodeMembership[j] -> leafCount) == i) {
	    Rprintf("\n %10d %10d %12.4f", i, j, ensembleRun[j]);
	    for (k = 1; k <= _xSize; k++) {
	      Rprintf("%12.4f", _observation[k][j]);
	    }
	  }
	}
      }
    }
    Rprintf("\n");
  }
  if (_mup & MUP_PERF) {
    imputeFlag = FALSE;
    if (mode == RSF_GROW) {
      statusPtr = orgStatusPtr = _status;
      timePtr = orgTimePtr = _time;
      if (_mRecordSize > 0) {
        imputeFlag = TRUE;
        obsSize = _observationSize;
        mRecordMap = _mRecordMap;
        mRecordSize = _mRecordSize;
        mvSize = _mvSize;
        mvSign = _mvSign;
        mvIndex = _mvIndex;
      }
    }
    else {
      statusPtr = orgStatusPtr = _fstatus;
      timePtr = orgTimePtr = _ftime;
      if (_fmRecordSize > 0) {
        imputeFlag = TRUE;
        obsSize = _fobservationSize;
        mRecordMap = _fmRecordMap;
        mRecordSize = _fmRecordSize;
        mvSize = _fmvSize;
        mvSign = _fmvSign;
        mvIndex = _fmvIndex;
      }
    }
    if (imputeFlag == TRUE) {
      outcomeFlag = TRUE;
      for (p=1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          statusPtr = dvector(1, obsSize);
          mStatusPtr = dvector(1, mRecordSize);
          unsignedIndex = abs(mvIndex[p]);
          orgOutcomePtr  = orgStatusPtr;
          newOutcomePtr = statusPtr;
          mOutcomePtr = mStatusPtr;
          break;
        case TIME_IDX:
          timePtr = dvector(1, obsSize);
          mTimePtr = dvector(1, mRecordSize);
          unsignedIndex = abs(mvIndex[p]);
          orgOutcomePtr = orgTimePtr;
          newOutcomePtr = timePtr;
          mOutcomePtr = mTimePtr;
          break;
        default:
          outcomeFlag = FALSE;
          break;
        }
        if (outcomeFlag == TRUE) {
          for (i=1; i <= obsSize; i++) {
            mPredictorFlag = TRUE;
            if (mRecordMap[i] == 0) {
              mPredictorFlag = FALSE;
              newOutcomePtr[i] = orgOutcomePtr[i];
            }
            else {
              if (mvSign[unsignedIndex][mRecordMap[i]] == 0) {
                mPredictorFlag = FALSE;
                newOutcomePtr[i] = orgOutcomePtr[i];
                mOutcomePtr[mRecordMap[i]] = orgOutcomePtr[i];
              }
            }
            if (mPredictorFlag == TRUE) {
              mOutcomePtr[mRecordMap[i]] = NA_REAL;
            }
          }  
        }  
        else {
          p = mvSize;
        }
      } 
      imputeConcordance(mode,
                        b,
                        mRecordBootFlag,
                        mvImputation,
                        mStatusPtr,
                        mTimePtr,
                        statusPtr,
                        timePtr);
    }  
    concordanceIndex = getConcordanceIndex(obsSize, 
                                           statusPtr, 
                                           timePtr, 
                                           ensembleRun, 
                                           ensembleDen);
    if (ISNA(concordanceIndex)) {
      performance = NA_REAL;
    }
    else {
      performance = 1.0 - concordanceIndex;
    }
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Error Rate:      %20.4f", performance);
      Rprintf("\n");
    }
    if (imputeFlag > 0) {
      outcomeFlag = TRUE;
      for (p=1; p <= mvSize; p++) {
        switch (mvIndex[p]) {
        case CENS_IDX:
          free_dvector(statusPtr, 1, obsSize);
          free_dvector(mStatusPtr, 1, mRecordSize);
          break;
        case TIME_IDX:
          free_dvector(timePtr, 1,  obsSize);
          free_dvector(mTimePtr, 1, mRecordSize);
          break;
        default:
          outcomeFlag = FALSE;
          break;
        }
        if (outcomeFlag == FALSE) {
          p = mvSize;
        }
      }  
    }  
  }  
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetPerformance() EXIT ...\n");
  }
  return performance;
}
