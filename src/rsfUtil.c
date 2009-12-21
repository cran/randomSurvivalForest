////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.0
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
#include       "rsfTree.h"
#include     "rsfImpute.h"
#include "rsfImportance.h"
#include       "rsfUtil.h"
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
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\ngetConcordanceIndex() ENTRY ...\n");
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  concordancePairSize = concordanceWorseCount = 0;
  for (i=1; i < size; i++) {
    for (j=i+1; j <= size; j++) {
      if (oobCount[i] != 0  && oobCount[j] != 0) {
        if ( (timePtr[i] > timePtr[j] && statusPtr[j] > 0) ||
             (timePtr[i] == timePtr[j] && statusPtr[j] > 0 && statusPtr[i] == 0) ) {
          concordancePairSize += 2;
          if ((polarity * predictedOutcome[j]) > (polarity * predictedOutcome[i])) { 
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[j] - predictedOutcome[i]) < EPSILON) {
            concordanceWorseCount += 1;
          }  
        }
        else if ( (timePtr[j] > timePtr[i] && statusPtr[i] > 0) ||
                  (timePtr[j] == timePtr[i] && statusPtr[i] > 0 && statusPtr[j] == 0) ) {
          concordancePairSize += 2;
          if ((polarity * predictedOutcome[i]) > (polarity * predictedOutcome[j])) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 1;
          }  
        }
        else if ( timePtr[i] == timePtr[j] && statusPtr[i] > 0 && statusPtr[j] > 0) {
          concordancePairSize += 2;
          if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 2;
          }
          else {
            concordanceWorseCount += 1;
          }
        }
      }  
    }  
  }  
  if (getTraceFlag() & ENSB_LOW_TRACE) {
    Rprintf("\nPredicted Outcome used in Concordance Index Calculations:  ");
    Rprintf("\n        count     OOBcount       Status         Time    Mortality");
    for (i=1; i <= size; i++) {
      Rprintf("\n %12d %12d %12d %12.4f, %12.4f", i, oobCount[i], (uint) statusPtr[i], timePtr[i], predictedOutcome[i]);
    }
    Rprintf("\n");
    Rprintf("\nConcordance pair and error update complete:");  
    Rprintf("\nCount of concordance pairs:         %20d", concordancePairSize);
    Rprintf("\nCount of pairs with worse outcome:  %20d", concordanceWorseCount);
    Rprintf("\n");
  }
  if (concordancePairSize == 0) {
    result = NA_REAL;
  }
  else {
    result = ((double) concordanceWorseCount / (double) concordancePairSize);
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _cindxTime = _cindxTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\ngetConcordanceIndex() EXIT() ...\n");
  }
  return result;
}
void getNelsonAalenEstimate(uint mode, double **nelsonAalen, uint treeID) {
  uint leaf, i, j;
  Node *parent;
  uint priorTimePointIndex, currentTimePointIndex;
  double estimate;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getNelsonAalenEstimate() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  for (leaf=1; leaf <= _leafCount_[treeID]; leaf++) {
    parent = getTerminalNode(mode, leaf);
    if (parent != NULL) { 
      for (i=1; i <= _masterTimeSize; i++) {
        nodeParentAtRisk[i] = nodeParentDeath[i] = 0;
      }
      parent -> mortality = 0.0;
      for (i=1; i <= _observationSize; i++) {
        if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
          for (j=1; j <= _masterTimeIndex[_bootMembershipIndex[i]]; j++) {
            nodeParentAtRisk[j] ++;
          }
          if (_status[_bootMembershipIndex[i]] > 0) {
            nodeParentDeath[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
          }
        }
      }
      priorTimePointIndex = 0;
      currentTimePointIndex = 1;
      for (j=1; j <= _sortedTimeInterestSize; j++) {
        for (i = priorTimePointIndex + 1; i <= _masterTimeSize; i++) {
          if (_masterTime[i] <= _timeInterest[j]) {
            currentTimePointIndex = i;
          }
          else {
            i = _masterTimeSize;
          }
        }
        if (getTraceFlag() & TURN_OFF_TRACE) {
          Rprintf("\nNelson-Aalen at risk and death counts for node:  %10d \n", leaf);
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
        nelsonAalen[j][leaf] = estimate;
        priorTimePointIndex = currentTimePointIndex;
      }  
      for (j=2; j <= _sortedTimeInterestSize; j++) {
        nelsonAalen[j][leaf] += nelsonAalen[j-1][leaf];
      }
      for (j = 1; j <= _sortedTimeInterestSize; j++) {
        parent -> mortality += nelsonAalen[j][leaf];
      }
    }
    else {
    }      
  }  
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
  if (getTraceFlag() & TURN_OFF_TRACE) {
    Rprintf("\nTree specific Nelson-Aalen estimator matrix:  %10d \n", treeID);
    Rprintf("          ");
    for (i=1; i <= _leafCount_[treeID]; i++) {
      Rprintf("%10d", i);
    }
    Rprintf("\n");
    for (j=1; j <= _sortedTimeInterestSize; j++) {
      Rprintf("%10d", j);
      for (i=1; i <= _leafCount_[treeID]; i++) {
        Rprintf("%10.4f", nelsonAalen[j][i]);
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _hazrdTime = _hazrdTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getNelsonAalenEstimate() EXIT ...\n");  
  }
}
void updateEnsembleCHF(uint      mode, 
                       uint      treeID,
                       double  **cumulativeHazard,
                       double   *mortality) {
  uint obsSize;
  unsigned char oobFlag, fullFlag, selectionFlag, mortalityFlag;
  double ***ensemblePtr;
  uint     *ensembleDenPtr;  
  Node    **nodeMembershipPtr;
  uint i, j, k, p, n;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEnsembleCHF() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    if (_oobSampleSize[treeID] > 0) {
      oobFlag = TRUE;
    }
    else {
      oobFlag = FALSE;
    }
    fullFlag = TRUE;
    mortalityFlag = TRUE;
    nodeMembershipPtr = _nodeMembership;
   break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    oobFlag = FALSE;
    fullFlag = TRUE;
    mortalityFlag = TRUE;
    nodeMembershipPtr = _fnodeMembership;
    break;
  case RSF_INTR:
    obsSize = _fobservationSize;
    if (_foobSampleSize[treeID] > 0) {
      oobFlag = TRUE;
    }
    else {
      oobFlag = FALSE;
    }
    fullFlag = FALSE;
    mortalityFlag = TRUE;
    nodeMembershipPtr = _nodeMembership;
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) { 
    if (oobFlag == TRUE) {
      ensemblePtr = _oobEnsemblePtr;
      ensembleDenPtr = _oobEnsembleDen;
      selectionFlag = FALSE;
    }
    else {
      if (fullFlag == TRUE) {
        ensemblePtr = _fullEnsemblePtr;
        ensembleDenPtr = _fullEnsembleDen;
        selectionFlag = ACTIVE;
      }
      else {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Unknown case in switch encountered. ");
        Rprintf("\nRSF:  Please Contact Technical Support.");
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
      }
    }
    if (mortalityFlag == TRUE) {
      for (i=1; i <= obsSize; i++) {
        mortality[i]  = 0.0;
      }
    }
    for (i=1; i <= obsSize; i++) {
      p = nodeMembershipPtr[_individualIndex[i]] -> leafCount;
      if ((_genericMembershipFlag[_individualIndex[i]] == selectionFlag) || (selectionFlag == ACTIVE)) {
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          ensemblePtr[1][k][i] += cumulativeHazard[k][p];
        }
        ensembleDenPtr[i] ++;
      }
      if (mortalityFlag == TRUE) {
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          mortality[i] += ensemblePtr[1][k][i];
        }
        if (ensembleDenPtr[i] != 0) {
          mortality[i] = mortality[i] / ensembleDenPtr[i];
        }
      }
    }  
    if (getTraceFlag() & ENSB_LOW_TRACE) {
      if (oobFlag == TRUE) {
        Rprintf("\nOOB Ensemble calculations follow: \n");        
      }
      else {
        if (fullFlag == TRUE) {
          Rprintf("\nFULL Ensemble calculations follow: \n");        
        }
      }
      Rprintf("\nEnsemble Estimator Numerator calculation: \n");
      Rprintf("          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (k=1; k <= _sortedTimeInterestSize; k++) {
        Rprintf("%10d", k);
        for (i=1; i <= obsSize; i++) {
          Rprintf("%10.4f", ensemblePtr[1][k][i]);
        }
        Rprintf("\n");
      }
      Rprintf("\nEnsemble Estimator Denominator calculation: \n");
      Rprintf("          ");  
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n          ");
      for (i=1; i <= obsSize; i++) {
        Rprintf("%10d", ensembleDenPtr[i]);
      }
      Rprintf("\n");
      if (mortalityFlag == TRUE) {
        Rprintf("\nMortality:  \n");
        Rprintf("          ");
        for (n=1; n <= obsSize; n++) {
          Rprintf("%10d", n);
        }
        Rprintf("\n          ");
        for (n=1; n <= obsSize; n++) {
          Rprintf("%10.4f", mortality[n]);
        }
        Rprintf("\n");
      }
      if (oobFlag == TRUE) {
        Rprintf("\nClassification of OOB elements:  ");
        Rprintf("\n       Leaf     Indiv            EE      predictors ->  ");
        for (p=1; p <= _leafCount_[treeID]; p++) {
          for (i = 1; i <= obsSize; i++) {
            if ( _bootMembershipFlag[_individualIndex[i]] == FALSE) {
              if ((nodeMembershipPtr[i] -> leafCount) == p) {
                Rprintf("\n %10d %10d %12.4f", p, i, mortality[i]);
                for (j = 1; j <= _xSize; j++) {
                  Rprintf("%12.4f", _observation[j][_individualIndex[i]]);
                }
              }
            }
          }
        }
        Rprintf("\n");
      }
    }
    if (mortalityFlag == TRUE) {
      mortalityFlag = FALSE;
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      if (fullFlag == TRUE) {
        fullFlag = FALSE;
      }
    }
  }  
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _ensblTime = _ensblTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEnsembleCHF() EXIT ...\n");  
  }
}
void getTreeSpecificSubSurvivalAndDistribution(uint mode, uint treeID) {
  uint leaf, i, j, k, q, m;
  uint nodeEventTimeSize;
  Node *parent;
  double estimate, headValue;
  uint priorTimePointIndex, currentTimePointIndex;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getTreeSpecificSubSurvivalAndDistribution() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  if (_eventTypeSize == 1) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Illegal sub-survival function call.");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
  }
  uint *nodeParentDeath  = uivector(1, _masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, _masterTimeSize);
  uint **nodeParentEvent  = uimatrix(1, _eventTypeSize, 1, _masterTimeSize);
  uint *nodeEventTimeIndex = uivector(1, _masterTimeSize);
  for (leaf=1; leaf <= _leafCount_[treeID]; leaf++) {
    parent = getTerminalNode(mode, leaf);
    if (parent != NULL) { 
      for (i=1; i <= _masterTimeSize; i++) {
        nodeParentAtRisk[i] = nodeParentDeath[i] = 0;
        for (j = 1; j <= _eventTypeSize; j++) {
          nodeParentEvent[j][i] = 0;
        }
        nodeEventTimeIndex[i] = 0;
      }
      for (j = 1; j <= _eventTypeSize; j++) {
        parent -> poe[j] = 0;
      }
      parent -> eventCount = 0;
      for (i=1; i <= _observationSize; i++) {
        if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
          for (k=1; k <= _masterTimeIndex[_bootMembershipIndex[i]]; k++) {
            nodeParentAtRisk[k] ++;
          }
          if (_status[_bootMembershipIndex[i]] > 0) {
            nodeParentDeath[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
            nodeParentEvent[_eventTypeIndex[(uint) _status[_bootMembershipIndex[i]]]][_masterTimeIndex[_bootMembershipIndex[i]]] ++;
            (parent -> poe)[_eventTypeIndex[(uint) _status[_bootMembershipIndex[i]]]] ++;
            parent -> eventCount += 1;
          }
        }
      }
      nodeEventTimeSize = 0;
      for (k=1; k <= _masterTimeSize; k++) {
        if (nodeParentDeath[k] > 0) {
          nodeEventTimeIndex[++nodeEventTimeSize] = k;
        }
      }
      if (getTraceFlag() & ENSB_HGH_TRACE) {
        Rprintf("\nSurvival at risk, death and event counts for node:  %10d  \n", leaf);
        Rprintf("mstTimIdx ");
        for (k=1; k <= _masterTimeSize; k++) {
          Rprintf("%10d", k);
        }
        Rprintf("\n");
        Rprintf("At Risk   ");
        for (k=1; k <= _masterTimeSize; k++) {
          Rprintf("%10d", nodeParentAtRisk[k]);
        }
        Rprintf("\n");
        Rprintf("Deaths    ");
        for (k=1; k <= _masterTimeSize; k++) {
          Rprintf("%10d", nodeParentDeath[k]);
        }
        Rprintf("\n");
        for (j=1; j <= _eventTypeSize; j++) {
          Rprintf("Ev %7d", j);
          for (k=1; k <= _masterTimeSize; k++) {
            Rprintf("%10d", nodeParentEvent[j][k]);
          }
          Rprintf("\n");
        }
        Rprintf("\n");
        Rprintf("\nNode specific POE counts of length [_eventTypeSize] for:  (tree, node) = (%10d, %10d)  \n", treeID, leaf);
        Rprintf("Event Types:  ");
        for (j=1; j <= _eventTypeSize; j++) {
          Rprintf("%10d ", _eventType[j]);
        }
        Rprintf("\n");
        Rprintf("              ");
        for (j=1; j <= _eventTypeSize; j++) {
          Rprintf("%10d ", (parent -> poe)[j]);
        }
        Rprintf("\n");
      }  
      double **subDensity = dmatrix(1, _eventTypeSize, 1, nodeEventTimeSize);
      double **survival = dmatrix(1, nodeEventTimeSize, 1, _leafCount_[treeID]);
      for (q=1; q <= nodeEventTimeSize; q++) {
        estimate = 1.0;
        if (nodeParentDeath[nodeEventTimeIndex[q]] > 0) {
          if (nodeParentAtRisk[nodeEventTimeIndex[q]] >= 1) {
            estimate = 1.0 - ((double) nodeParentDeath[nodeEventTimeIndex[q]] / nodeParentAtRisk[nodeEventTimeIndex[q]]);
          }
        }
        survival[q][leaf] = estimate;
      }  
      for (q=2; q <= nodeEventTimeSize; q++) {
        survival[q][leaf] *= survival[q-1][leaf];
      }
      for (j=1; j <= _eventTypeSize; j++) {      
        subDensity[j][1] = 0.0;
        if (nodeParentEvent[j][nodeEventTimeIndex[1]] > 0) {
          if (nodeParentAtRisk[nodeEventTimeIndex[1]] >= 1) {
            subDensity[j][1] = 1.0 * ((double) nodeParentEvent[j][nodeEventTimeIndex[1]] / nodeParentAtRisk[nodeEventTimeIndex[1]]);
          }
        }
        for (q=2; q <= nodeEventTimeSize; q++) {
          subDensity[j][q] = 0.0;
          if (nodeParentEvent[j][nodeEventTimeIndex[q]] > 0) {
            if (nodeParentAtRisk[nodeEventTimeIndex[q]] >= 1) {
              subDensity[j][q] = survival[q-1][leaf] * ((double) nodeParentEvent[j][nodeEventTimeIndex[q]] / nodeParentAtRisk[nodeEventTimeIndex[q]]);
            }
          }
        }
      }
      for (j = 1; j <= _eventTypeSize; j++) {
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          (parent -> subSurvival)[j][k] = 0.0;
        }
      }
      priorTimePointIndex = 0;
      currentTimePointIndex = 1;
      for (i = 1; i <= nodeEventTimeSize; i++) {
        for (k = priorTimePointIndex + 1; k <= _sortedTimeInterestSize; k++) {
          if (_timeInterest[k] <= _masterTime[nodeEventTimeIndex[i]] ) {
            currentTimePointIndex = k;
          }
          else {
            k = _sortedTimeInterestSize;
          }
        }
        for (j = 1; j <= _eventTypeSize; j++) {
          for (q = 1; q <= i; q++) {
            (parent -> subSurvival)[j][currentTimePointIndex] += subDensity[j][q];
          }
          if (i > 1) {
            for(k = priorTimePointIndex + 1; k < currentTimePointIndex; k++) {
              (parent -> subSurvival)[j][k] = (parent -> subSurvival)[j][priorTimePointIndex];
            }
            if (i == nodeEventTimeSize) {
              for(k = currentTimePointIndex + 1; k <= _sortedTimeInterestSize; k++) {
                (parent -> subSurvival)[j][k] = (parent -> subSurvival)[j][currentTimePointIndex];
              }
            }
          }
        }
        priorTimePointIndex = currentTimePointIndex;
      }
      if (getTraceFlag() & ENSB_HGH_TRACE) {
        Rprintf("\nNode specific survival function of length [nodeEventTimeSize] for:  (tree, node) = (%10d, %10d)  \n", treeID, leaf);
        Rprintf("              mTimIdx       time   survival \n");
        for (q=1; q <= nodeEventTimeSize; q++) {
          Rprintf("%10d %10d %10.4f %10.4f \n", q, nodeEventTimeIndex[q], _masterTime[nodeEventTimeIndex[q]], survival[q][leaf]);
        }
        Rprintf("\nNode specific sub-density function [_eventTypeSize] x [nodeEventTimeSize] for:  (tree, node) = (%10d, %10d)  \n", treeID, leaf);
        Rprintf("             mTimIdx       time ");
        for (j=1; j <= _eventTypeSize; j++) {
          Rprintf("%10d ", j);
        }
        Rprintf("\n");
        for (q=1; q <= nodeEventTimeSize; q++) {
          Rprintf("%10d %10d %10.4f ", q, nodeEventTimeIndex[q], _masterTime[nodeEventTimeIndex[q]]);
          for (j=1; j <= _eventTypeSize; j++) {
            Rprintf("%10.4f ", subDensity[j][q]);
          }
          Rprintf("\n");
        }
        Rprintf("\nNode specific pre sub-survival function (before head adjsutment) [_eventTypeSize] x [_sortedTimeInterestSize] x [_leafCount_[treeID]] for:  (tree, node) = (%10d, %10d)  \n", treeID, leaf);
        Rprintf("                 time ");
        for (j=1; j <= _eventTypeSize; j++) {
          Rprintf("%10d ", j);
        }
        Rprintf("\n");
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          Rprintf("%10d %10.4f ", k, _timeInterest[k]);
          for (j=1; j <= _eventTypeSize; j++) {
            Rprintf("%10.4f ", (parent -> subSurvival)[j][k]);
          }
          Rprintf("\n");
        }
      }
      for (j=1; j <= _eventTypeSize; j++) {
        headValue = (double) (parent -> poe)[j] / parent -> eventCount;
        headValue = (headValue > (parent -> subSurvival)[j][_sortedTimeInterestSize]) ? headValue : (parent -> subSurvival)[j][_sortedTimeInterestSize];
        if (headValue > 0) {
          for (k=1; k <= _sortedTimeInterestSize; k++) {
            (parent -> subSurvival)[j][k] = headValue - (parent -> subSurvival)[j][k];
            if ((parent -> subSurvival)[j][k] < - EPSILON) {
              Rprintf("\nRSF:  *** ERROR *** ");
              Rprintf("\nRSF:  Negative sub-survival value at (event, time, leaf, subSurvival) = (%4d, %12.4f, %4d, %12.4f).", j, _timeInterest[k], leaf, (parent -> subSurvival)[j][k]);
              Rprintf("\nRSF:  Please Contact Technical Support.");
              Rprintf("\nRSF:  The application will now exit.\n");
              exit(TRUE);
            }
          }
        }
        else {
          for (k=1; k <= _sortedTimeInterestSize; k++) {
            (parent -> subSurvival)[j][k] = 0.0;
          }
        }
      }
      if (getTraceFlag() & ENSB_HGH_TRACE) {
        Rprintf("\nNode specific sub-survival function (after head adjustment) [_eventTypeSize] x [_sortedTimeInterestSize] x [_leafCount_[treeID]] for:  (tree, node) = (%10d, %10d)  \n", treeID, leaf);
        Rprintf("                 time ");
        for (j=1; j <= _eventTypeSize; j++) {
          Rprintf("%10d ", j);
        }
        Rprintf("\n");
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          Rprintf("%10d %10.4f ", k, _timeInterest[k]);
          for (j=1; j <= _eventTypeSize; j++) {
            Rprintf("%10.4f ", (parent -> subSurvival)[j][k]);
          }
          Rprintf("\n");
        }
      }
      for (j = 1; j <= _eventTypeSize; j++) {
        if ((parent -> poe)[j] == 0) {
          for (k=1; k <= _sortedTimeInterestSize; k++) {
            for (m = 1; m <= _eventTypeSize; m++) {        
              (parent -> subSurvival)[j][k] = 0.0;
            }
          }
        }
      }
      free_dmatrix(subDensity, 1, _eventTypeSize, 1, nodeEventTimeSize);
      free_dmatrix(survival, 1, nodeEventTimeSize, 1, _leafCount_[treeID]);
    }
    else {
    }      
  }  
  free_uivector(nodeEventTimeIndex, 1, _masterTimeSize);
  free_uivector(nodeParentDeath, 1, _masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, _masterTimeSize);
  free_uimatrix(nodeParentEvent, 1, _eventTypeSize, 1, _masterTimeSize);
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _chazrdTime = _chazrdTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getTreeSpecificSubSurvivalAndDistribution() EXIT ...\n");  
  }
}
void updateEnsembleSubSurvivalAndDistribution(uint      mode, 
                                              uint      treeID,
                                              double  **conditionalMortality) {
  uint obsSize;
  unsigned char oobFlag, fullFlag, selectionFlag, mortalityFlag;
  double  **poePtr;
  double ***ensSubSurvivalPtr;
  double ***ensemblePtr;
  Node    **nodeMembershipPtr;
  double  value;
  uint i, j, k, n, p;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEnsembleSubSurvivalAndDistribution() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  if (_eventTypeSize == 1) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Illegal ensemble sub-survival call.");
      Rprintf("\nRSF:  Please Contact Technical Support.");
      Rprintf("\nRSF:  The application will now exit.\n");
      exit(TRUE);
  }
  poePtr            = NULL;  
  ensSubSurvivalPtr = NULL;  
  ensemblePtr       = NULL;  
  switch (mode) {
  case RSF_GROW:
    obsSize = _observationSize;
    if (_oobSampleSize[treeID] > 0) {
      oobFlag = TRUE;
    }
    else {
      oobFlag = FALSE;
    }
    fullFlag = TRUE;
    mortalityFlag = TRUE;
    nodeMembershipPtr = _nodeMembership;
    break;
  case RSF_PRED:
    obsSize = _fobservationSize;
    oobFlag = FALSE;
    fullFlag = TRUE;
    mortalityFlag = TRUE;
    nodeMembershipPtr = _fnodeMembership;
    break;
  case RSF_INTR:
    obsSize = _fobservationSize;
    if (_foobSampleSize[treeID] > 0) {
      oobFlag = TRUE;
    }
    else {
      oobFlag = FALSE;
    }
    fullFlag = FALSE;
    mortalityFlag = TRUE;
    nodeMembershipPtr = _nodeMembership;
    break;
  default:
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Unknown case in switch encountered. ");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
    break;
  }
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) { 
    if (oobFlag == TRUE) {
      poePtr = _oobPOEPtr;
      ensSubSurvivalPtr = _oobSubSurvivalPtr;
      ensemblePtr = _oobEnsemblePtr;
    }
    else {
      if (fullFlag == TRUE) {
        poePtr = _fullPOEPtr;
        ensSubSurvivalPtr = _fullSubSurvivalPtr;
        ensemblePtr = _fullEnsemblePtr;
      }
    }
    if (mortalityFlag == TRUE) {
      for (i=1; i <= obsSize; i++) {
        for (j = 1; j <= _eventTypeSize; j++) {
          conditionalMortality[j][i] = 0.0;
        }
      }
    }
    for (i=1; i <= obsSize; i++) {
      for (j = 1; j <= _eventTypeSize; j++) {
        for (k = 1; k <= _sortedTimeInterestSize; k++) {      
          ensemblePtr[j+1][k][i] = 0.0;
        }
      }
    }
    for (i=1; i <= obsSize; i++) {
      p = nodeMembershipPtr[_individualIndex[i]] -> leafCount;
      selectionFlag = TRUE;
      if (oobFlag == TRUE) {
        if (_bootMembershipFlag[_individualIndex[i]] == FALSE) {
          selectionFlag = TRUE;
        }
        else {
          selectionFlag = FALSE;
        }
      }
      if (selectionFlag) {
        for (j = 1; j <= _eventTypeSize; j++) {
          poePtr[j][i] += (double) (nodeMembershipPtr[_individualIndex[i]] -> poe)[j] / nodeMembershipPtr[_individualIndex[i]] -> eventCount;
        }
        for (j = 1; j <= _eventTypeSize; j++) {
          for (k=1; k <= _sortedTimeInterestSize; k++) {
            ensSubSurvivalPtr[j][k][i] += nodeMembershipPtr[_individualIndex[i]] -> subSurvival[j][k];
          }
        }
      }  
      for (j = 1; j <= _eventTypeSize; j++) {
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          if(ensSubSurvivalPtr[j][k][i] > 0) {
            if (getTraceFlag() & TIME_DEF_TRACE) {
            }
            if (poePtr[j][i] > 0) {
              value = ensSubSurvivalPtr[j][k][i] / poePtr[j][i];
              value = (value <= 1.0) ? value : 1.0;
              ensemblePtr[j+1][k][i] = - log (value);
            }
            else {
              value = ensSubSurvivalPtr[j][k][i] / 1.0;
              value = (value <= 1.0) ? value : 1.0;
              ensemblePtr[j+1][k][i] = - log (value);
            }
            if (getTraceFlag() & TIME_DEF_TRACE) {
            }
          }
          else {
            if (poePtr[j][i] > 0) {
              if (k > 1) {
                ensemblePtr[j+1][k][i] = ensemblePtr[j+1][k-1][i];
              }
              else {
                ensemblePtr[j+1][k][i] = 0.0;
              }
            }
            else {
              ensemblePtr[j+1][k][i] = 0.0;  
            }
          }
        }
      }
      if (mortalityFlag == TRUE) {
        for (j = 1; j <= _eventTypeSize; j++) {
          for (k=1; k <= _sortedTimeInterestSize; k++) {            
            conditionalMortality[j][i] += ensemblePtr[j + 1][k][i];
          }
        }
      }
    }  
    if (getTraceFlag() & ENSB_LOW_TRACE) {
      if (oobFlag == TRUE) {
        Rprintf("\nOOB Conditional Ensemble calculations follow: \n");        
      }
      else {
        if (fullFlag == TRUE) {
          Rprintf("\nFULL Conditional Ensemble calculations follow: \n");        
        }
      }
      for (j = 1; j <= _eventTypeSize; j++) {
        Rprintf("\nSub-SRV Numerator calculation:  Event %10d \n", j);
        Rprintf("          ");
        for (n=1; n <= obsSize; n++) {
          Rprintf("%10d", n);
        }
        Rprintf("\n");
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          Rprintf("%10d", k);
          for (n=1; n <= obsSize; n++) {
            Rprintf("%10.4f", ensSubSurvivalPtr[j][k][n]);
          }
          Rprintf("\n");
        }
      }
      Rprintf("\nProbability of Event:  \n");
      Rprintf("          ");
      for (n=1; n <= obsSize; n++) {
        Rprintf("%10d", n);
      }
      Rprintf("\n");
      for (j=1; j <= _eventTypeSize; j++) {
        Rprintf("%10d", j);
        for (n=1; n <= obsSize; n++) {
          Rprintf("%10.4f", poePtr[j][n]);
        }
        Rprintf("\n");
      }
      for (j = 1; j <= _eventTypeSize; j++) {
        Rprintf("\nEnsemble Conditional CHF:  Event %10d \n", j);
        Rprintf("          ");
        for (n=1; n <= obsSize; n++) {
          Rprintf("%10d", n);
        }
        Rprintf("\n");
        for (k=1; k <= _sortedTimeInterestSize; k++) {
          Rprintf("%10d", k);
          for (n=1; n <= obsSize; n++) {
            Rprintf("%10.4f", ensemblePtr[j+1][k][n]);
          }
          Rprintf("\n");
        }
      }
      if (mortalityFlag == TRUE) {
        Rprintf("\nConditional Mortality:  \n");
        Rprintf("          ");
        for (n=1; n <= obsSize; n++) {
          Rprintf("%10d", n);
        }
        Rprintf("\n");
        for (j = 1; j <= _eventTypeSize; j++) {       
          Rprintf("%10d", _eventType[j]);
          for (n=1; n <= obsSize; n++) {
            Rprintf("%10.4f", conditionalMortality[j][n]);
          }
          Rprintf("\n");
        }
      }
    }  
    if (mortalityFlag == TRUE) {
      mortalityFlag = FALSE;
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      if (fullFlag == TRUE) {
        fullFlag = FALSE;
      }
    }
  }  
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _censblTime = _censblTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEnsembleSubSurvivalAndDistribution() EXIT ...\n");  
  }
}
void getConditionalConcordanceArrays(uint     j, 
                                     double  *statusPtr, 
                                     double  *timePtr, 
                                     double  *mortalityPtr, 
                                     uint    *genericEnsembleDenPtr,
                                     double  *subsettedStatus,
                                     double  *subsettedTime,
                                     double  *subsettedMortality,
                                     uint    *subsettedEnsembleDen) {
  uint i;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getConditionalConcordanceArrays() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  if (_eventTypeSize == 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt to update event type subsets in a non-CR analysis.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  for (i = 1; i <= _meIndividualSize[j]; i++) {
    subsettedStatus[i]      = statusPtr[_eIndividual[j][i]];
    subsettedTime[i]        = timePtr[_eIndividual[j][i]];
    subsettedMortality[i]   = mortalityPtr[_eIndividual[j][i]];
    subsettedEnsembleDen[i] = genericEnsembleDenPtr[_eIndividual[j][i]];
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _censblTime = _censblTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getConditionalConcordanceArrays() EXIT ...\n");  
  }
}
void getMeanSurvivalTime(uint mode, double *meanSurvivalTime, uint treeID) {
  uint leaf, i;
  double totalDeathTime;
  uint deathCount;
  Node *parent;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getMeanSurvivalTime() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  for (leaf=1; leaf <= _leafCount_[treeID]; leaf++) {
    totalDeathTime = 0;
    deathCount = 0;
    parent = getTerminalNode(mode, leaf);
    if (parent != NULL) { 
      for (i=1; i <= _observationSize; i++) {
        if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
          if (_status[_bootMembershipIndex[i]] > 0) {
            totalDeathTime += _masterTime[_masterTimeIndex[_bootMembershipIndex[i]]];
            deathCount ++;
          }
        }
      }
      if (deathCount > 0) {
        meanSurvivalTime[leaf] = totalDeathTime / deathCount;
        parent -> mortality = meanSurvivalTime[leaf];
      }
      else {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Zero death count encountered in node:  %10d", leaf);
        Rprintf("\nRSF:  Please Contact Technical Support.");
        Rprintf("\nRSF:  The application will now exit.\n");
        exit(TRUE);
      }
    }
    else {
    }      
  }  
  if (getTraceFlag() & SUMM_HGH_TRACE) {
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
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _hazrdTime = _hazrdTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  getMeanSurvivalTime() EXIT ...\n");  
  }
}
void updateEnsembleSurvivalTime(uint     mode, 
                                uint     treeID, 
                                double  *meanSurvivalTime,
                                double *mortality) {
  uint i, j, k, p;
  uint obsSize;
  uint *genericEnsembleDenPtr;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEnsembleSurvivalTime() ENTRY ...\n");  
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _benchTime = clock();
  }
  switch (mode) {
  case RSF_GROW:
    for (i=1; i <= _observationSize; i++) {
      p = _nodeMembership[i] -> leafCount;
      if (_oobSampleSize[treeID] > 0) {
        if ( _bootMembershipFlag[i] == FALSE ) {
          _oobEnsemblePtr[1][1][i] += meanSurvivalTime[p];
          _oobEnsembleDen[i] ++;
        }
      }
      _fullEnsemblePtr[1][1][i] += meanSurvivalTime[p];
      _fullEnsembleDen[i] ++;
      if (_oobEnsembleDen[i] != 0) {
        mortality[i] = _oobEnsemblePtr[1][1][i];
        mortality[i] = mortality[i] / _oobEnsembleDen[i];
      }
      else {
        mortality[i] = 0;
      }
    }
    break;
  case RSF_PRED:
    for (i=1; i <= _fobservationSize; i++) {
      p = _fnodeMembership[i] -> leafCount;
      _fullEnsemblePtr[1][1][i] += meanSurvivalTime[p];
      _fullEnsembleDen[i] ++;
      if (_fullEnsembleDen[i] != 0) {
        mortality[i] = _fullEnsemblePtr[1][1][i];
        mortality[i] = mortality[i] / _fullEnsembleDen[i];
      }
      else {
        mortality[i] = 0;
      }
    }
    break;
  case RSF_INTR:
    for (i=1; i <= _fobservationSize; i++) {
      p = _nodeMembership[_intrIndividual[i]] -> leafCount;
      if (_foobSampleSize[treeID] > 0) {
        if ( _bootMembershipFlag[_intrIndividual[i]] == FALSE ) {
          _oobEnsemblePtr[1][1][i] += meanSurvivalTime[p];
          _oobEnsembleDen[i] ++;
        }
      }
      if (_oobEnsembleDen[i] != 0) {
        mortality[i] = _oobEnsemblePtr[1][1][i];
        mortality[i] = mortality[i] / _oobEnsembleDen[i];
      }
      else {
        mortality[i] = 0;
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
  if (getTraceFlag() & SUMM_HGH_TRACE) {
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
      for (k=1; k <= 1; k++) {
        Rprintf("%10d", k);
        for (i=1; i <= obsSize; i++) {
          Rprintf("%10.4f", _oobEnsemblePtr[1][k][i]);
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
      for (k=1; k <= 1; k++) {
        Rprintf("%10d", k);
        for (i=1; i <= obsSize; i++) {
          Rprintf("%10.4f", _fullEnsemblePtr[1][k][i]);
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
        Rprintf("%10.4f", mortality[i]);
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
    if (mode == RSF_GROW) {
      Rprintf("\nClassification of OOB elements:  ");
      Rprintf("\n       Leaf     Indiv            EE      predictors ->  ");
      for (j=1; j <= _leafCount_[treeID]; j++) {
        for (i = 1; i <= _observationSize; i++) {
          if (_bootMembershipFlag[i] == FALSE) {
            if ((_nodeMembership[i] -> leafCount) == j) {
              Rprintf("\n %10d %10d %12.4f", j, i, mortality[i]);
              for (k = 1; k <= _xSize; k++) {
                Rprintf("%12.4f", _observation[k][i]);
              }
            }
          }
        }
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & TIME_DEF_TRACE) {
    _ensblTime = _ensblTime + (clock() - _benchTime);
  }
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\nRSF:  updateEnsembleSurvivalTime() EXIT ...\n");  
  }
}
void getVariablesUsed(Node *parent, uint *varUsedPtr) {
  if (getTraceFlag() & OUTP_DEF_TRACE) {
    Rprintf("\ngetVariablesUsed() ENTRY ...\n");  
  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    varUsedPtr[parent -> splitParameter] ++;
     getVariablesUsed(parent ->  left, varUsedPtr);
     getVariablesUsed(parent -> right, varUsedPtr);
  }
  if (getTraceFlag() & OUTP_DEF_TRACE) {
    Rprintf("\ngetVariablesUsed() EXIT ...\n");  
  }
  return;
}
void updateEnsembleCalculations (char      multipleImputeFlag,
                                 uint      mode,
                                 Node     *rootPtr,
                                 uint      b,
                                 char    **dmRecordBootFlag,
                                 double ***dmvImputation) {
  uint i, j;
  uint obsSize;
  double  **cumulativeHazard;
  double   *meanSurvivalTime;
  double   *mortality;
  double  **conditionalMortality;
  uint     *genericEnsembleDenPtr;
  double *statusPtr, *orgStatusPtr, *mStatusPtr;
  double *timePtr, *orgTimePtr, *mTimePtr;
  uint *varUsedRow;
  char concordanceImputeFlag;
  int concordancePolarity;
  double concordanceIndex;
  double *crPerformanceVector;
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nupdateEnsembleCalculations() ENTRY ...\n");
  }
  cumulativeHazard     = NULL;  
  meanSurvivalTime     = NULL;  
  conditionalMortality = NULL;  
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
  mortality = dvector(1, obsSize);
  if (_opt & (OPT_POUT_TYPE)) {
    meanSurvivalTime = dvector(1, _leafCount_[b]);
    concordancePolarity = -1;
    getMeanSurvivalTime(mode, meanSurvivalTime, b);
    updateEnsembleSurvivalTime(mode, b, meanSurvivalTime, mortality);
  }
  else {
    cumulativeHazard = dmatrix(1, _sortedTimeInterestSize, 1, _leafCount_[b]);
    concordancePolarity = 1;
    getNelsonAalenEstimate(mode, cumulativeHazard, b);
    updateEnsembleCHF(mode, b, cumulativeHazard, mortality);
    if (_eventTypeSize > 1) {
      conditionalMortality = dmatrix(1, _eventTypeSize, 1, obsSize);
      getTreeSpecificSubSurvivalAndDistribution(mode, b);
      updateEnsembleSubSurvivalAndDistribution(mode, b, conditionalMortality);
    }
  }  
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nRSF:  Predicted outcome calculation complete.");  
  }
  if (_opt & OPT_VIMP) {
    getVariableImportance(mode,
                          _leafCount_[b],
                          rootPtr,
                          b);
  }
  if (_opt & (OPT_POUT_TYPE)) {
    free_dvector(meanSurvivalTime, 1, _leafCount_[b]);
  }
  else {
    free_dmatrix(cumulativeHazard, 1, _sortedTimeInterestSize, 1, _leafCount_[b]);
  }
  if (_opt & OPT_PERF) {
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
                                           mortality, 
                                           genericEnsembleDenPtr);
    if (ISNA(concordanceIndex)) {
      _performancePtr[1][b] = NA_REAL;
    }
    else {
      _performancePtr[1][b] = 1.0 - concordanceIndex;
    }
    if (_eventTypeSize > 1) {
      crPerformanceVector = dvector(1, _eventTypeSize);
      getConditionalPerformance(mode,
                                concordancePolarity, 
                                obsSize,
                                statusPtr, 
                                timePtr,
                                conditionalMortality,
                                genericEnsembleDenPtr,
                                crPerformanceVector);
      for (j=1; j <=_eventTypeSize; j++) {
        _performancePtr[1+j][b] = crPerformanceVector[j];
      }
      free_dvector(crPerformanceVector, 1, _eventTypeSize);
    }
    if (getTraceFlag() & SUMM_USR_TRACE) {
      Rprintf("\nRSF:  Error Rate for treeID:  %10d", b);
      if (_eventTypeSize > 1) {
        for (j = 1; j <= _eventTypeSize; j++) {
          Rprintf("                [%2d]", j);
        }
      }
      Rprintf("\n                   ");
      Rprintf("%20.4f", _performancePtr[1][b]);
      if (_eventTypeSize > 1) {
        for (j = 1; j <= _eventTypeSize; j++) {
          Rprintf("%20.4f", _performancePtr[1+j][b]);
        }
      }
      Rprintf("\n");
    }
    if (concordanceImputeFlag > 0) {
      free_dvector(mStatusPtr, 1, obsSize);
      free_dvector(mTimePtr, 1, obsSize);
    }  
  }  
  free_dvector(mortality, 1, obsSize);
  if (!(_opt & (OPT_POUT_TYPE))) {
    if (_eventTypeSize > 1) {
      free_dmatrix(conditionalMortality, 1, _eventTypeSize, 1, obsSize);
    }
  }
  if (_opt & OPT_VUSE) {
    if (_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
      varUsedRow = _varUsedPtr[b];
    }
    else {
      varUsedRow = _varUsedPtr[1];
    }
    getVariablesUsed(rootPtr, varUsedRow);
    if (getTraceFlag() & SUMM_LOW_TRACE) {
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
  if (getTraceFlag() & SUMM_LOW_TRACE) {
    Rprintf("\nupdateEnsembleCalculations() EXIT ...\n");
  }
}
void getConditionalPerformance (uint     mode,
                                uint     concordancePolarity, 
                                uint     obsSize,
                                double  *statusPtr, 
                                double  *timePtr,
                                double **conditionalMortality,
                                uint    *ensembleDenPtr,
                                double  *performanceVector) {
  uint   mRecordSize;
  int  **mvSign;
  uint  *mRecordIndex;
  double concordanceIndex;
  uint j;
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\ngetConditionalPerformance() ENTRY() ...\n");
  }
  if (_eventTypeSize == 1) { 
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Attempt at conditional performance updates in a non-CR analysis.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (_mStatusSize > 0) {
    if (mode == RSF_GROW) {
      mRecordSize = _mRecordSize;
      mvSign = _mvSign;
      mRecordIndex = _mRecordIndex;
    }
    else {
      mRecordSize = _fmRecordSize;
      mvSign = _fmvSign;
      mRecordIndex = _fmRecordIndex;
    }
    updateEventTypeSubsets(statusPtr, mRecordSize, mvSign, mRecordIndex);
  }
  double *subsettedStatus    = dvector(1, obsSize);
  double *subsettedTime      = dvector(1, obsSize);
  double *subsettedMortality = dvector(1, obsSize);
  uint *subsettedEnsembleDen = uivector(1, obsSize);
  for (j = 1; j <= _eventTypeSize; j++) {
    getConditionalConcordanceArrays(j, 
                                    statusPtr, 
                                    timePtr, 
                                    conditionalMortality[j],
                                    ensembleDenPtr,
                                    subsettedStatus,
                                    subsettedTime,
                                    subsettedMortality,
                                    subsettedEnsembleDen);
    concordanceIndex = getConcordanceIndex(concordancePolarity,
                                           _meIndividualSize[j], 
                                           subsettedStatus,
                                           subsettedTime,
                                           subsettedMortality,
                                           subsettedEnsembleDen);
    if (ISNA(concordanceIndex)) {
      performanceVector[j] = NA_REAL;
    }
    else {
      performanceVector[j] = 1.0 - concordanceIndex;
    }
  }
  free_dvector(subsettedStatus, 1, obsSize);
  free_dvector(subsettedTime, 1, obsSize);
  free_dvector(subsettedMortality, 1, obsSize);
  free_uivector(subsettedEnsembleDen, 1, obsSize);
  if (getTraceFlag() & SUMM_MED_TRACE) {
    Rprintf("\ngetConditionalPerformance() EXIT() ...\n");
  }
}
