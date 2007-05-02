//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 2.1.0
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
#include  "rsfUtil.h"
#include      "rsf.h"
SEXP sexpVector[RSF_SEXP_CNT];
uint stackCount;
char *sexpString[RSF_SEXP_CNT] = {
  "",              
  "",              
  "fullEnsemble",  
  "oobEnsemble",   
  "performance",   
  "leafCount",     
  "proximity",     
  "importance",    
  "treeID",        
  "nodeID",        
  "parmID",        
  "spltPT",        
  "seed"           
};
uint     *_treeID_;
uint     *_nodeID_;
uint     *_parmID_;
double   *_spltPT_;
int      *_seed_;
double   *_fullEnsemble_;
double   *_oobEnsemble_;
double   *_performance_;
uint     *_leafCount_;
uint     *_proximity_;
double   *_varImportance_;
uint      _mup;
int      *_seedPtr;
uint      _splitRule;
uint      _randomCovariateCount;
uint      _forestSize;
uint      _minimumDeathCount;
uint      _observationSize;
double   *_time;
uint     *_status;
uint      _xSize;
double   *_xData;
double   *_randomCovariateWeight;
uint      _fobservationSize;
double   *_ftime;
uint     *_fstatus;
double   *_fxData;
uint      _timeInterestSize;
double   *_timeInterest;
double  **_observation;
uint     *_masterTimeIndex;
double   *_masterTime;
double  **_fobservation;
Node    **_nodeMembership;
uint     *_bootMembershipIndex;
char     *_bootMembershipFlag;
Node    **_fnodeMembership;
uint      _traceFlagDiagLevel;
uint      _traceFlagIterValue;
uint      _traceFlagToggler;
clock_t    start;
clock_t    now;
char logRankApprox(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit) {
  double delta, deltaNum, deltaDen;
  double eOneMono;
  double deltaMax = -EPSILON;
  uint i,j,k, leadingIndex;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount, rightDeathCount;
  uint localMembershipSize, localDeathTimeSize;
  uint *nodeParentDeath  = uivector(1, masterDeathTimeSize);
  uint *nodeParentAtRisk = uivector(1, masterDeathTimeSize);
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, masterDeathTimeSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nlogRankApprox() ENTRY ...\n");
  }
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = 0;  
  leadingIndex = 0;  
  for (i=1; i <= masterTimeSize; i++) {
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
    for (i=1; i <= masterTimeSize; i++) {
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
    if (getTraceFlag() & DL3_TRACE) {
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
            }
          }
        }  
        rightDeathCount = parentDeathCount - leftDeathCount;
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nRunning Split Totals: \n");
          Rprintf("    PARdeath    LFTdeath   RGTdeath\n");
          Rprintf("%10d %10d %10d \n", parentDeathCount, leftDeathCount, rightDeathCount);
        }
        if ((leftDeathCount  >= (_minimumDeathCount)) && (rightDeathCount  >= (_minimumDeathCount))) {
          deltaNum = sqrt(parentDeathCount) * (leftDeathCount - eOneMono);
          deltaDen = eOneMono * (parentDeathCount - eOneMono);
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
    free_dvector(nelsonAalen , 1, localDeathTimeSize);
    free_uivector (survivalTimeIndex, 1, localMembershipSize);
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum threshold of deaths");
      Rprintf("\nRSF:  encountered in logRank().  This condition");
      Rprintf("\nRSF:  can occur on the root split and is potentially");
      Rprintf("\nRSF:  nominal behaviour.  \n");
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nBest Split Statistics: \n");
    Rprintf("  SplitParm SplitIndex        Delta \n");
    Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
  }
  free_uivector(nodeParentDeath, 1, masterDeathTimeSize);
  free_uivector(nodeParentAtRisk, 1, masterDeathTimeSize);
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, masterDeathTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
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
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit,
  uint   **masterSplitOrder) {
  double delta, deltaNum, deltaDen;
  double meanSurvRank, varSurvRank;
  double deltaMax = -EPSILON;
  uint i,j,k;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount, rightDeathCount;
  uint localMembershipSize, localDeathTimeSize, leftMembershipSize;
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, masterTimeSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nlogrankScore() ENTRY ...\n");
  }
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = 0;  
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
    for (i=1; i <= masterTimeSize; i++) {
      if (localDeathTimeCount[i] > 0) {
        localDeathTimeSize++;
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nLocal Death Time Size:  %10d \n", localDeathTimeSize);
    }
    if (localDeathTimeSize > 1) {
      actualCovariateCount = selectRandomCovariates(parent, randomCovariateIndex);
      uint   *localSplitRank  = uivector(1, localMembershipSize);
      uint *survivalTimeIndexRank = uivector(1, localMembershipSize);
      double *survivalRank = dvector(1, localMembershipSize);
      for (i=1; i <= actualCovariateCount; i++) {
        k = 0;
        for (j=1; j <= _observationSize; j++) {
          if (_nodeMembership[_bootMembershipIndex[masterSplitOrder[i][j]]] == parent) {
            localSplitRank[++k] = masterSplitOrder[i][j];
          }
        }
        for (k = 1; k <= localMembershipSize; k++) {
          survivalTimeIndexRank[k] = 0;
          for (j = 1; j <= localMembershipSize; j++) {
            if ( _masterTimeIndex[localMembershipIndex[j]]  <= _masterTimeIndex[_bootMembershipIndex[localSplitRank[k]]] ) {
              survivalTimeIndexRank[k] ++;
            }
          }
        }
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\nLocal Membership Information for Parent Node: \n");
          Rprintf("               MEMBidx    SPLTval    TIMEidx -> SORTidx \n");
          for (k=1; k <=  localMembershipSize; k++) {
            Rprintf("%10d %10d %10.4f %10d %10d\n", k,
                    _bootMembershipIndex[localSplitRank[k]], _observation[i][_bootMembershipIndex[localSplitRank[k]]],
                    _masterTimeIndex[_bootMembershipIndex[localSplitRank[k]]], survivalTimeIndexRank[k]);
          }
        }
        meanSurvRank = varSurvRank = 0;
        for (k = 1; k <= localMembershipSize; k++) {
          survivalRank[k] = 0;
          for (j = 1; j <= survivalTimeIndexRank[k]; j++) {
            survivalRank[k] = survivalRank[k] + ((double) _status[_bootMembershipIndex[localSplitRank[j]]] / (localMembershipSize - survivalTimeIndexRank[j] + 1) );
          }
          survivalRank[k] = _status[_bootMembershipIndex[localSplitRank[k]]] - survivalRank[k];
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
          leftDeathCount = 0;
          leftMembershipSize = 0;
          deltaNum = 0.0;
          for (k=1; k <= localMembershipSize; k++) {
            if (_observation[randomCovariateIndex[i]][_bootMembershipIndex[localSplitRank[k]]] <= masterSplit[randomCovariateIndex[i]][j]) {
              if (getTraceFlag() & DL3_TRACE) {
                Rprintf("\nMember of Left Daughter (index):   %10d %10d \n", k, _bootMembershipIndex[localSplitRank[k]]);
              }
              leftMembershipSize ++;
              if (_status[_bootMembershipIndex[localSplitRank[k]]] == 1) {
                leftDeathCount ++;
              }
              deltaNum = deltaNum + survivalRank[k];
            }
          }  
          rightDeathCount = parentDeathCount - leftDeathCount;
          if (getTraceFlag() & DL3_TRACE) {
            Rprintf("\nRunning Split Totals: \n");
            Rprintf("    PARcnt      LFTcnt     RGTcnt   PARdeath    LFTdeath   RGTdeath\n");
            Rprintf("%10d %10d %10d %10d %10d %10d \n", localMembershipSize, leftMembershipSize, localMembershipSize - leftMembershipSize, parentDeathCount, leftDeathCount, rightDeathCount);
          }
          if ((leftDeathCount  >= (_minimumDeathCount)) && (rightDeathCount  >= (_minimumDeathCount))) {
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
      free_uivector(localSplitRank, 1, localMembershipSize);
      free_uivector(survivalTimeIndexRank, 1, localMembershipSize);
      free_dvector(survivalRank, 1, localMembershipSize);
    }   
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum threshold of deaths");
      Rprintf("\nRSF:  encountered in conserveEvents().  This condition");
      Rprintf("\nRSF:  can occur on the root split and is potentially");
      Rprintf("\nRSF:  nominal behaviour.  \n");
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nBest Split Statistics: \n");
    Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
    Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
  }
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, masterTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
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
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit) {
  double delta, nelsonAalenSumLeft, nelsonAalenSumRight;
  double deltaMax = -EPSILON;
  uint i,j,k,m, index;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount, rightDeathCount;
  uint localMembershipSize, localDeathTimeSize;
  uint *nodeParentDeath  = uivector(1, masterDeathTimeSize);
  uint *nodeLeftDeath    = uivector(1, masterDeathTimeSize);
  uint *nodeRightDeath    = uivector(1, masterDeathTimeSize);
  uint *nodeParentAtRisk = uivector(1, masterDeathTimeSize);
  uint *nodeLeftAtRisk   = uivector(1, masterDeathTimeSize);
  uint *nodeRightAtRisk   = uivector(1, masterDeathTimeSize);
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, masterDeathTimeSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nconserveEvents() ENTRY ...\n");
  }
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  for (i=1; i <= masterTimeSize; i++) {
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
    for (i=1; i <= masterTimeSize; i++) {
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
    if (localDeathTimeSize > 1) {
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
          leftDeathCount = rightDeathCount = 0;
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
          }  
          for (k=1; k <= localDeathTimeSize; k++) {
            nodeRightDeath[k] = nodeParentDeath[k] - nodeLeftDeath[k];
            nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k];
            leftDeathCount += nodeLeftDeath[k];
          }
          rightDeathCount = parentDeathCount - leftDeathCount;
          if (getTraceFlag() & DL3_TRACE) {
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
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum threshold of deaths");
      Rprintf("\nRSF:  encountered in conserveEvents().  This condition");
      Rprintf("\nRSF:  can occur on the root split and is potentially");
      Rprintf("\nRSF:  nominal behaviour.  \n");
    }
  }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nBest Split Statistics: \n");
      Rprintf("  SplitParm SplitIndex   SplitValue        Delta \n");
      Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
    }
  free_uivector(nodeParentDeath, 1, masterDeathTimeSize);
  free_uivector(nodeLeftDeath,   1, masterDeathTimeSize);
  free_uivector(nodeRightDeath,  1, masterDeathTimeSize);
  free_uivector(nodeParentAtRisk, 1, masterDeathTimeSize);
  free_uivector(nodeLeftAtRisk,   1, masterDeathTimeSize);
  free_uivector(nodeRightAtRisk,  1, masterDeathTimeSize);
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, masterDeathTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
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
char logRank(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit) {
  double delta, deltaNum, deltaDen;
  double deltaMax = -EPSILON;
  uint i,j,k,m, index;
  uint actualCovariateCount;
  uint parentDeathCount, leftDeathCount, rightDeathCount;
  uint localMembershipSize, localDeathTimeSize;
  uint *nodeParentDeath  = uivector(1, masterDeathTimeSize);
  uint *nodeLeftDeath    = uivector(1, masterDeathTimeSize);
  uint *nodeRightDeath   = uivector(1, masterDeathTimeSize);
  uint *nodeParentAtRisk = uivector(1, masterDeathTimeSize);
  uint *nodeLeftAtRisk   = uivector(1, masterDeathTimeSize);
  uint *nodeRightAtRisk  = uivector(1, masterDeathTimeSize);
  uint *randomCovariateIndex = uivector(1, _xSize);
  uint *localMembershipIndex = uivector(1, _observationSize);
  uint *localDeathTimeCount = uivector(1, masterTimeSize);
  uint *localDeathTimeIndex = uivector(1, masterDeathTimeSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nlogRank() ENTRY ...\n");
  }
  *splitParameterMax = *splitValueMax = 0;
  localMembershipSize = localDeathTimeSize = 0;
  parentDeathCount = 0;
  delta = -1.0;  
  for (i=1; i <= masterTimeSize; i++) {
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
    for (i=1; i <= masterTimeSize; i++) {
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
    if (localDeathTimeSize > 1) {
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
          leftDeathCount = rightDeathCount = 0;
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
                  if (getTraceFlag() & DL2_TRACE) {
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
          if (getTraceFlag() & DL3_TRACE) {
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
  }  
  else {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRSF:  *** WARNING *** ");
      Rprintf("\nRSF:  Less than twice the minimum threshold of deaths");
      Rprintf("\nRSF:  encountered in logRank().  This condition");
      Rprintf("\nRSF:  can occur on the root split and is potentially");
      Rprintf("\nRSF:  nominal behaviour.  \n");
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nBest Split Statistics: \n");
    Rprintf("  SplitParm SplitIndex        Delta \n");
    Rprintf(" %10d %10d %12.4f %12.4f \n", *splitParameterMax, *splitValueMax, masterSplit[*splitParameterMax][*splitValueMax], deltaMax);
  }
  free_uivector(nodeParentDeath, 1, masterDeathTimeSize);
  free_uivector(nodeLeftDeath,   1, masterDeathTimeSize);
  free_uivector(nodeRightDeath,  1, masterDeathTimeSize);
  free_uivector(nodeParentAtRisk, 1, masterDeathTimeSize);
  free_uivector(nodeLeftAtRisk,   1, masterDeathTimeSize);
  free_uivector(nodeRightAtRisk,  1, masterDeathTimeSize);
  free_uivector(randomCovariateIndex, 1, _xSize);
  free_uivector(localMembershipIndex, 1, _observationSize);
  free_uivector(localDeathTimeCount, 1, masterTimeSize);
  free_uivector(localDeathTimeIndex, 1, masterDeathTimeSize);
  if (((*splitParameterMax) > 0) && ((*splitValueMax) > 0)) {
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
uint selectRandomCovariates(
  Node *parent,
  uint *covariateIndex) {
  uint i,j,k;
  double randomValue;
  uint maxCovariateCount;
  uint actualCovariateCount;
  uint unselectedCovariateCount;
  char *randomSplitVector = cvector(1, _xSize);
  double *cdf = dvector(1, _xSize);
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nselectRandomCovariate() ENTRY ...\n");
  }
  maxCovariateCount = 0;
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
    randomValue = ran1(_seedPtr)*cdf[unselectedCovariateCount];
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
char getBestSplit(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit,
  uint   **masterSplitOrder) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetBestSplit() ENTRY ...\n");
  }
  if (_splitRule == LOG_RANK) {
    return logRank(parent,
                   splitParameterMax,
                   splitValueMax,
                   masterDeathTimeSize,
                   masterTimeSize,
                   masterSplit);
  }
  else if (_splitRule == CONSERVE_EVENTS) {
    return conserveEvents(parent,
                          splitParameterMax,
                          splitValueMax,
                          masterDeathTimeSize,
                          masterTimeSize,
                          masterSplit);
  }
  else if (_splitRule == LOG_RANK_SCORE) {
    return logRankScore(parent,
                        splitParameterMax,
                        splitValueMax,
                        masterDeathTimeSize,
                        masterTimeSize,
                        masterSplit,
                        masterSplitOrder);
  }
  else if (_splitRule == LOG_RANK_APPROX) {
    return logRankApprox(parent,
                          splitParameterMax,
                          splitValueMax,
                          masterDeathTimeSize,
                          masterTimeSize,
                          masterSplit);
  }
  else {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid split rule:  %10d", _splitRule);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetBestSplit() EXIT ...\n");
  }
}
char forkAndUpdate(
  uint  *leafCount,
  Node  *parent,
  uint   splitParameter,
  uint   splitValueIndex,
  double splitValue) {
  uint i;
  uint leftDeathCount, rightDeathCount;
  char result;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetforkAndUpdate() ENTRY ...\n");
  }
  result = forkNode(parent, splitParameter, splitValueIndex, splitValue);
  if (result == TRUE) {
    leftDeathCount = rightDeathCount = 0;
    (*leafCount)++;
    for (i = 1; i <= _observationSize; i++) {
      if (_nodeMembership[i] == parent) {
        if (_observation[splitParameter][i] <= splitValue) {
          _nodeMembership[i] = parent -> left;
          ((parent -> left) -> leafCount) = (parent -> leafCount);
        }
        else {
          _nodeMembership[i] = parent -> right;
          ((parent -> right) -> leafCount) = *leafCount;
        }
      }
    }
    for (i = 1; i <= _observationSize; i++) {
      if (_nodeMembership[_bootMembershipIndex[i]] == parent -> left) {
        if (_status[_bootMembershipIndex[i]] == 1) leftDeathCount ++;
      }
      else if (_nodeMembership[_bootMembershipIndex[i]] == parent -> right) {
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
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\ngetforkAndUpdate() EXIT ...\n");
  }
  return result;
}
char makeTree(
  Node    *parent,
  uint    *leafCount, 
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit,
  uint   **masterSplitOrder) {
  char splitResult, forkResult;
  uint splitParameterMax, splitValueMax;
  uint i;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nmakeTree() ENTRY ...\n");
  }
  splitResult = getBestSplit(parent,
                             & splitParameterMax,
                             & splitValueMax,
                             masterDeathTimeSize,
                             masterTimeSize,
                             masterSplit,
                             masterSplitOrder);
  if (splitResult == TRUE) {
    forkResult = forkAndUpdate(leafCount,
                               parent,
                               splitParameterMax,
                               splitValueMax,
                               masterSplit[splitParameterMax][splitValueMax]);
    if (forkResult == TRUE) {
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nNode Membership:  \n");
        for (i=1; i <=  _observationSize; i++) {
          Rprintf("%10d %10d \n", i, _nodeMembership[i] -> leafCount);
        }
      }
      if ((parent -> left) -> splitFlag == TRUE) {
        if (getTraceFlag() & DL1_TRACE) {
          Rprintf("\nmakeTree() LEFT:  \n");
        }
        makeTree(parent -> left,
                 leafCount,
                 masterDeathTimeSize,
                 masterTimeSize,
                 masterSplit,
                 masterSplitOrder);
      }
      if ((parent -> right) -> splitFlag == TRUE) {
        if (getTraceFlag() & DL1_TRACE) {
          Rprintf("\nmakeTree() RIGHT:  \n");
        }
        makeTree(parent -> right,
                 leafCount,
                 masterDeathTimeSize,
                 masterTimeSize,
                 masterSplit,
                 masterSplitOrder);
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
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\ngetBestSplit() FAILED ...\n");
    }
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nmakeTree() EXIT ...\n");
  }
  return splitResult;
}
char bootstrap(
  uint     mode,
  uint     b,
  Node    *rootPtr,
  uint    *oobSampleSize,
  double **masterSplit,
  uint    *masterSplitSize,
  uint   **masterSplitOrder) {
  char result;
  uint **masterSplitBounds;
  uint ibSampleSize;
  uint leadingIndex;
  uint i,j,k;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nbootstrap() ENTRY ...\n");
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Bootstrap sample:  %10d ", b);  
  }
  result = FALSE;  
  if (_mup & MUP_TREE) {
    *_seedPtr = - abs(*_seedPtr);
    _seed_[b] = *_seedPtr;
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nSaved Seed for Bootstrap:  ");
      Rprintf("\n%10d %20d ", b, _seed_[b]);
      Rprintf("\n");
    }
  }
  else if (mode == RSF_PRED) {
    *_seedPtr = _seed_[b];
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nRestored Seed for Bootstrap:  ");
      Rprintf("\n%10d %20d ", b, *_seedPtr);
      Rprintf("\n");
    }
  }
  for (i=1; i <= _observationSize; i++) {
    _nodeMembership[i] = rootPtr;
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
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nOOB Size:  %10d ", oobSampleSize[b]);
    Rprintf("\n\nRoot Membership (all data):  ");
    Rprintf("\n     index       leaf    IBindex\n");
    for (i=1; i <=  _observationSize; i++) {
      Rprintf("%10d %10d %10d \n", i, _nodeMembership[i] -> leafCount, _bootMembershipIndex[i]);
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Bootstrapping and initialization of root node complete.");  
  }
  if (oobSampleSize[b] == 0) {
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  OOB sample size is zero.  Bootstrap sample has been discarded.");  
    }
    result = FALSE;
  }
  else if(mode == RSF_GROW) {
    if (getTraceFlag() & DL3_TRACE) {
      Rprintf("\nRaw Bootstrap MasterSplit Data:  ");  
    }
    k = 1;
    for (i=1; i <= _observationSize; i++) {
      if (_bootMembershipFlag[i] == TRUE) {
        for (j=1; j <= _xSize; j++) {
          masterSplit[j][k] = _observation[j][i];
        }
        if (getTraceFlag() & DL3_TRACE) {
          Rprintf("\n%10d ", i);  
          for (j=1; j <= _xSize; j++) {
            Rprintf("%10.4f ", masterSplit[j][k]);
          }
        }
        k++;
      }
    }
    ibSampleSize = _observationSize - oobSampleSize[b];
    for (i=1; i <= _xSize; i++) {
      hpsort(masterSplit[i], ibSampleSize);
   }
    if (getTraceFlag() & DL3_TRACE) {
      Rprintf("\n\nRaw Sorted Bootstrap MasterSplit Data:  ");
      for (i=1; i <= ibSampleSize; i++) {
        Rprintf("\n%10d", i);
        for (j=1; j <= _xSize; j++) {
          Rprintf(" %10.4f", masterSplit[j][i]);
        }
      }
    }
    masterSplitBounds = uimatrix(1, _xSize, 1, 2);
    for (i=1; i <= _xSize; i++) {
      masterSplitBounds[i][1] = 1;
      masterSplitSize[i] = leadingIndex = 1;
      for (j=2; j <= ibSampleSize; j++) {
        if (masterSplit[i][j] > masterSplit[i][leadingIndex]) {
          masterSplitSize[i]++;
          leadingIndex++;
          masterSplit[i][leadingIndex] = masterSplit[i][j];
        }
      }
      masterSplitBounds[i][2] = leadingIndex;
      for (j=masterSplitSize[i]+1; j <= ibSampleSize; j++) {
        masterSplit[i][j] = 0;
      }
    }
    if (_splitRule == LOG_RANK_SCORE) {
      double *predictorValue = dvector(1, _observationSize);
      for (j=1; j <= _xSize; j++) {
        for (i=1; i <= _observationSize; i++) {
          masterSplitOrder[j][i] = 0;
          predictorValue[i] = _observation[j][_bootMembershipIndex[i]];
        }
        indexx(_observationSize, predictorValue, masterSplitOrder[j]);
      }
      free_dvector(predictorValue, 1, _observationSize);
      if (getTraceFlag() & DL3_TRACE) {
        Rprintf("\n\nBootstraped Ordered Predictor Matrix:  ");
        Rprintf("\n     index   observations -> \n");
        for (i = 1; i <= _observationSize; i++) {
          Rprintf("%10d", i);
          for (j=1; j <= _xSize; j++) {
            Rprintf(" %10d", masterSplitOrder[j][i]);
          }
          Rprintf("\n");
        }
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\n\nBootstraped MasterSplit Data:  ");
      for (i=1; i <= ibSampleSize; i++) {
        Rprintf("\n%10d", i);
        for (j=1; j <= _xSize; j++) {
          Rprintf(" %10.4f", masterSplit[j][i]);
        }
      }
      Rprintf("\n\nSize of Permissible Splits:  ");
      Rprintf("\n          ");
      for (j=1; j <= _xSize; j++) {
        Rprintf(" %10d", masterSplitSize[j]);
      }
      Rprintf("\n");
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nMaster Permissible Split Boundaries:  \n");
      for (i=1; i <= _xSize; i++) {
        Rprintf("%10d %10d %10d \n", i, masterSplitBounds[i][1], masterSplitBounds[i][2]);
      }
    }
    nrCopyMatrix(rootPtr -> permissibleSplit, masterSplitBounds, _xSize, 2);
    free_uimatrix(masterSplitBounds, 1, _xSize, 1, 2);
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Initialization of master split data complete.");  
    }
    result = TRUE;
  }  
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nbootstrap() EXIT ...\n");
  }
  return result;
}
SEXP rsfGrow(SEXP traceFlag,
             SEXP mup,  
             SEXP seedPtr,  
             SEXP splitRule,  
             SEXP randomCovariateCount,  
             SEXP forestSize,  
             SEXP minimumDeathCount,
             SEXP observationSize,
             SEXP time,
             SEXP status,
             SEXP xSize,
             SEXP xData,
             SEXP timeInterestSize,
             SEXP timeInterest,
             SEXP randomCovariateWeight
            ) {
  int seedValue         = INTEGER(seedPtr)[0];
  _seedPtr              = &seedValue;
  _splitRule            = INTEGER(splitRule)[0];
  _randomCovariateCount = INTEGER(randomCovariateCount)[0];
  _forestSize           = INTEGER(forestSize)[0];
  _minimumDeathCount    = INTEGER(minimumDeathCount)[0];
  _observationSize      = INTEGER(observationSize)[0];
  _time                 =    REAL(time);  _time--;
  _status               = (uint*) INTEGER(status);  _status--;
  _xSize                = INTEGER(xSize)[0];
  _xData                =    REAL(xData);
  _timeInterestSize     = INTEGER(timeInterestSize)[0];
  _timeInterest         =    REAL(timeInterest);  _timeInterest--;
  _randomCovariateWeight =   REAL(randomCovariateWeight);  _randomCovariateWeight--;
  _mup                  = INTEGER(mup)[0] | MUP_PERF;  
  if (*_seedPtr >= 0) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Random seed must be less than zero.  \n");
    return R_NilValue;
  }
  if ( (_splitRule != LOG_RANK) && (_splitRule != CONSERVE_EVENTS) && (_splitRule != LOG_RANK_SCORE) && (_splitRule != LOG_RANK_APPROX) ) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Invalid split rule:  %10d \n", _splitRule);
    return R_NilValue;
  }
  if ( ((_randomCovariateCount < 1) || (_randomCovariateCount > _xSize)) ) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of random covariate parameters must be greater");
    Rprintf("\nRSF:  than zero and less than the total number of covariates:  %10d \n", _randomCovariateCount);
    return R_NilValue;
  }
  if (_minimumDeathCount < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Minimum number of deaths must be greater than zero:  %10d \n", _minimumDeathCount);
    return R_NilValue;
  }
  return rsf(RSF_GROW, INTEGER(traceFlag)[0]);
}
SEXP rsfPredict(SEXP traceFlag,
                SEXP mup,
                SEXP forestSize, 
                SEXP observationSize,
                SEXP time,
                SEXP status,
                SEXP xSize,
                SEXP xData,
                SEXP fobservationSize,
                SEXP ftime,
                SEXP fstatus,
                SEXP fxData,
                SEXP timeInterestSize,
                SEXP timeInterest,
                SEXP treeID,
                SEXP nodeID,
                SEXP parmID,
                SEXP spltPT,
                SEXP seed
               ) {
  int seedValue;
  _seedPtr              = &seedValue;
  _forestSize           = INTEGER(forestSize)[0];
  _observationSize      = INTEGER(observationSize)[0];
  _time                 =    REAL(time);  _time --;
  _status               = (uint*) INTEGER(status);  _status --;
  _xSize                = INTEGER(xSize)[0];
  _xData                =    REAL(xData);
  _fobservationSize     = INTEGER(fobservationSize)[0];
  _ftime                =    REAL(ftime);  _ftime --;
  _fstatus              = (uint*) INTEGER(fstatus);  _fstatus --;
  _fxData               =    REAL(fxData);
  _timeInterestSize     = INTEGER(timeInterestSize)[0];
  _timeInterest         =    REAL(timeInterest);  _timeInterest --;
  _treeID_              = (uint*) INTEGER(treeID);  _treeID_ --;
  _nodeID_              = (uint*) INTEGER(nodeID);  _nodeID_ --;
  _parmID_              = (uint*) INTEGER(parmID);  _parmID_ --;
  _spltPT_              =    REAL(spltPT);  _spltPT_ --;
  _seed_                = INTEGER(seed);  _seed_ --;
  _mup                  = INTEGER(mup)[0] & (~MUP_TREE);  
  _mup                  = _mup & (~MUP_VIMP);  
  if (_fobservationSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of prediction observations must be greater than one:  %10d \n", _fobservationSize);
    return R_NilValue;
  }
  return rsf(RSF_PRED, INTEGER(traceFlag)[0]);
}
SEXP rsf(uint mode, uint traceFlag) {
  uint sexpIndex;
  uint sortedTimeInterestSize;
  uint masterDeathTimeSize;
  uint masterTimeSize;
  double **masterSplit;        
  uint    *masterSplitSize;    
  uint   **masterSplitOrder;   
  Node   **root;
  double **oobEnsemblePtr;
  double **fullEnsemblePtr;
  double  *ensembleRun;
  uint    *ensembleDen;
  double **vimpEnsembleRun;
  uint    *oobSampleSize;
  uint     forestNodeCount;
  uint     rejectedTreeCount;  
  uint i, j, k, b;
  char splitResult;
  Node *rootPtr;
  char result;
  _traceFlagDiagLevel = (traceFlag | mode) & 0xFF;
  _traceFlagIterValue = traceFlag >> 8;
  updateTraceFlag(TRUE);
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Native code rsf() entry. \n");
  }
  start = clock();
  if (_forestSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of bootstrap iterations must be greater than zero:  %10d \n", _forestSize);
    return R_NilValue;
  }
  if (_observationSize < 2) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of observations must be greater than one:  %10d \n", _observationSize);
    return R_NilValue;
  }
  if (_timeInterestSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of time points of interest must be greater than zero:  %10d \n", _timeInterestSize);
    return R_NilValue;
  }
  if (_xSize < 1) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Parameter verification failed.");
    Rprintf("\nRSF:  Number of parameters must be greater than zero:  %10d \n", _xSize);
    return R_NilValue;
  }
  for (i =1 ; i <= _observationSize; i++) {
    if (_time[i] <= 0) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Parameter verification failed.");
      Rprintf("\nRSF:  Time vector must contain elements greater than zero:  %10d \n", _time[i]);
      return R_NilValue;
    }
    if ( (_status[i] != 0) && (_status[i] != 1) ) {
      Rprintf("\nRSF:  *** ERROR *** ");
      Rprintf("\nRSF:  Parameter verification failed.");
      Rprintf("\nRSF:  Status vector must contain elements equal to one or zero:  %10d \n", _status[i]);
      return R_NilValue;
    }
  }
  if (mode == RSF_GROW) {
    for (i = 1; i <= _xSize; i++) {
      if(_randomCovariateWeight[i] <= 0) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Parameter verification failed.");
        Rprintf("\nRSF:  Random covariate weight vector must contain elements greater than zero:  %12.4f \n", _randomCovariateWeight[i]);
        return R_NilValue;
      }
    }
  }
  if (mode == RSF_PRED) {
    for (i = 1; i <= _forestSize; i++) {
      if(_seed_[i] >= 0) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Parameter verification failed.");
        Rprintf("\nRSF:  Random seed vector must contain elements less than zero:  %10d \n", _seed_[i]);
        return R_NilValue;
      }
    }
    result = TRUE;
    for (i =1 ; i <= _fobservationSize; i++) {
      if ( (_ftime[i] <= 0) || ((_fstatus[i] != 0) && (_fstatus[i] != 1)) ) {
        result = FALSE;
      }
    }
    if (result == FALSE) {
      _mup = _mup & (~MUP_PERF);  
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Prediction data set time and status vector ignored.");  
        Rprintf("\nRSF:  No performance measure will be calculated.");  
      }
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    if (mode == RSF_GROW) {
      Rprintf("\nRSF:  Mode is GROW.");
      Rprintf("\nRSF:  Split rule is:                %10d", _splitRule);
      Rprintf("\nRSF:  Number of GROW observations:  %10d", _observationSize);
    }
    if (mode == RSF_PRED) {
      Rprintf("\nRSF:  Mode is PRED.");
      Rprintf("\nRSF:  Number of PRED observations:  %10d", _fobservationSize);
    }
    Rprintf("\nRSF:  Number of predictors:         %10d \n", _xSize);
  }
  stackPreDefinedCommonArrays(& oobSampleSize);
  if (mode == RSF_GROW) {
    stackPreDefinedGrowthArrays(
                                & masterSplit,
                                & masterSplitSize,
                                & masterSplitOrder);
  }
  else {
    stackPreDefinedPredictArrays();
  }
  initializeArrays(  mode,
                   & masterTimeSize,
                   & sortedTimeInterestSize,
                   & masterDeathTimeSize
                  );
  sexpIndex = stackDefinedOutputObjects(   mode,
                                           sortedTimeInterestSize,
                                           sexpString,
                                         & root,
                                         & oobEnsemblePtr,
                                         & fullEnsemblePtr,
                                         & ensembleRun,
                                         & ensembleDen,
                                         & _oobEnsemble_,
                                         & _fullEnsemble_,
                                         & _performance_,
                                         & _leafCount_,
                                         & _proximity_,
                                         & _varImportance_,
					 & vimpEnsembleRun,
                                         & _seed_,
                                         sexpVector
                                       );
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Number of trees in forest:  %10d \n", _forestSize);
  }
  forestNodeCount = 1;
  for (b = 1; b <= _forestSize; b++) {
    if ( (b == 1) || (b == _forestSize)) {
      updateTraceFlag(TRUE);
    }
    rootPtr = makeNode();  
    rootPtr -> parent = rootPtr;
    rootPtr -> left = NULL;
    rootPtr -> right = NULL;
    rootPtr -> splitValueIndex = 0;
    rootPtr -> splitParameter = 0;
    rootPtr -> splitFlag = TRUE;
    rootPtr -> leafCount = 1;
    if (_mup & MUP_TREE) {
      root[b] = rootPtr;
    }
    if (mode == RSF_GROW) {
      splitResult = bootstrap(mode,
                              b,
                              rootPtr,
                              oobSampleSize,
                              masterSplit,
                              masterSplitSize,
                              masterSplitOrder);
      if (splitResult) {
        splitResult = makeTree(rootPtr,
                               _leafCount_ + b,
                               masterDeathTimeSize,
                               masterTimeSize,
                               masterSplit,
                               masterSplitOrder);
      }
    }
    else {
      *_seedPtr = - _seed_[b];
      splitResult = bootstrap(mode,
                              b,
                              rootPtr,
                              oobSampleSize,
                              masterSplit,      
                              masterSplitSize,  
                              masterSplitOrder  
                              );
      if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
        Rprintf("\nIncoming Tree:  ");
        Rprintf("\n      tree       parm       node         splt \n");
      }
      splitResult = restoreTree(b,
                                rootPtr,
                                _leafCount_ + b,
                                & forestNodeCount,
                                _treeID_,
                                _nodeID_,
                                _parmID_,
                                _spltPT_);
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Tree construction complete:  %10d", b);  
        Rprintf("\nRSF:  Final leaf count:  %10d\n", _leafCount_[b]);
      }
      for (k = 1; k <= _observationSize; k++) {
        _nodeMembership[k] = getMembership(rootPtr, _observation, k);
      }
      for (k = 1; k <= _fobservationSize; k++) {
        _fnodeMembership[k] = getMembership(rootPtr,  _fobservation, k);
      }
    }
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nFinal GROW Membership (all data):  %10d", b);
      Rprintf("\n     index       leaf    IBindex\n");
      for (i=1; i <=  _observationSize; i++) {
        Rprintf("%10d %10d %10d \n", i, _nodeMembership[i] -> leafCount, _bootMembershipIndex[i]);
      }
    }
    if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
      Rprintf("\nFinal PRED Membership (all data):  %10d", b);
      Rprintf("\n     index       leaf\n");
      for (i=1; i <=  _fobservationSize; i++) {
        Rprintf("%10d %10d \n", i, _fnodeMembership[i] -> leafCount);
      }
    }
    if (_mup & MUP_PROX) {
      if (mode == RSF_GROW) {
        k = 0;
        for (i = 1; i <= _observationSize; i++) {
          k += i - 1;
          for (j = 1; j <= i; j++) {
            if ( (_nodeMembership[i] -> leafCount) == (_nodeMembership[j] -> leafCount) ) {
              _proximity_[k + j] ++;
            }
          }
        }
      }
      else {
        k = 0;
        for (i = 1; i <= _fobservationSize; i++) {
          k += i - 1;
          for (j = 1; j <= i; j++) {
            if ( (_fnodeMembership[i] -> leafCount) == (_fnodeMembership[j] -> leafCount) ) {
              _proximity_[k + j] ++;
            }
          }
        }
      }
      if (getTraceFlag() & DL0_TRACE) {
        Rprintf("\nRSF:  Proximity matrix calculation complete.");  
      }
    }
    if (splitResult) {
      _performance_[b] = getPerformance(mode,
                                        sortedTimeInterestSize,
                                        _leafCount_[b],
                                        masterTimeSize,
                                        oobEnsemblePtr,
                                        fullEnsemblePtr,
                                        ensembleRun,
                                        ensembleDen,
                                        oobSampleSize[b],
                                        rootPtr,
                                        vimpEnsembleRun);
    }
    else {
      if (b > 1) {
        _performance_[b] = _performance_[b-1];
      }
    }  
    if (!(_mup & MUP_TREE)) {
      freeTree(rootPtr);
    }
    updateTraceFlag(FALSE);
  }  
  updateTraceFlag(TRUE);
  if (_mup & MUP_VIMP) {
    for (k=1; k <= _xSize; k++) {
      for (i = 1; i <= _observationSize; i++) {
        if (ensembleDen[i] != 0) {
          vimpEnsembleRun[k][i] = vimpEnsembleRun[k][i] / ensembleDen[i];
        }
      }
      if (getTraceFlag() & DL1_TRACE) {
        Rprintf("\nConcordance Pairs for VIMP covariate:  %10d", k);  
      }
      _varImportance_[k] = 1.0 - getConcordanceIndex(_observationSize, 
                                                     _time, 
                                                     _status, 
                                                     vimpEnsembleRun[k], 
                                                     ensembleDen);
    }
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Variable Importance Measure: \n");
      Rprintf("          ");  
      for (k=1; k <= _xSize; k++) {
        Rprintf("%10d", k);
      }
      Rprintf("\n          ");
      for (k=1; k <= _xSize; k++) {
        Rprintf("%10.4f", _varImportance_[k]);
      }
      Rprintf("\n");
    }
  }
  if (getTraceFlag() & DL1_TRACE) {
    if (_mup & MUP_PROX) {
      Rprintf("\nProximity Matrix:  \n");
      k = 0;
      for (i = 1; i <= _observationSize; i++) {
        k += i - 1;
        for (j = 1; j <= i; j++) {
          Rprintf("%10d ", _proximity_[k + j]);
        }
        Rprintf("\n");
      }
    }
  }
  rejectedTreeCount = k = 0;
  for (b = 1; b <= _forestSize; b++) {
    if (mode == RSF_GROW) {
      if (oobSampleSize[b] == 0) {
        k ++;
      }
    }
    if (_leafCount_[b] == 1) {
      rejectedTreeCount ++;
    }
  }
  if ((getTraceFlag() & RSF_GROW) && (getTraceFlag() & DL0_TRACE)) {
    Rprintf("\nRSF:  Trees rejected (death threshold too high):  %10d ", rejectedTreeCount - k);
    Rprintf("\nRSF:  Trees rejected (OOB subset size was zero):  %10d ", k);
    Rprintf("\nRSF:  Trees rejected (total):                     %10d ", rejectedTreeCount);
    Rprintf("\nRSF:  Trees (total):                              %10d \n", _forestSize);
  }
  if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL0_TRACE)) {
    Rprintf("\nRSF:  Trees rejected (total):                     %10d ", rejectedTreeCount);
    Rprintf("\nRSF:  Trees (total):                              %10d \n", _forestSize);
  }
  if (rejectedTreeCount == _forestSize) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Insufficient trees for analysis.  \n");
  }
  if (mode == RSF_GROW) {
    for (i = 1; i <= _observationSize; i++) {
      for (j=1; j <= sortedTimeInterestSize; j++) {
        oobEnsemblePtr[j][i] = oobEnsemblePtr[j][i] / ensembleDen[i];
        fullEnsemblePtr[j][i] = fullEnsemblePtr[j][i] / _forestSize;
      }
    }
  }
  if (mode == RSF_PRED) {
    for (i = 1; i <= _fobservationSize; i++) {
      for (j=1; j <= sortedTimeInterestSize; j++) {
        fullEnsemblePtr[j][i] = fullEnsemblePtr[j][i] / _forestSize;
      }
    }
  }
  forestNodeCount = 0;
  for (b = 1; b <= _forestSize; b++) {
    forestNodeCount += (2 * _leafCount_[b]) - 1;
  }
  sexpIndex = stackVariableOutputObjects(
    sexpIndex, 
    forestNodeCount,
    sexpString,
    & _treeID_,
    & _nodeID_,
    & _parmID_,
    & _spltPT_,
    sexpVector);
  if (_mup & MUP_TREE) {
    forestNodeCount = 1;
    for (b = 1; b <= _forestSize; b++) {
      saveTree(b, root[b], & forestNodeCount, _treeID_, _nodeID_, _parmID_, _spltPT_);
    }
    forestNodeCount --;
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nGrown Forest Node Count:  %10d", forestNodeCount);
      Rprintf("\nGrown Forest Output:  ");
      Rprintf("\n   TREE_ID    NODE_ID    PARM_ID      SPLT_PT");
      for (i = 1; i <= forestNodeCount; i++) {
        Rprintf("\n%10d %10d %10d %12.4f",
          _treeID_[i],
          _nodeID_[i],
          _parmID_[i],
          _spltPT_[i]);
      }
      Rprintf("\n");
      Rprintf("\nGrown Forest Membership by Observation:  ");
      Rprintf("\n    OBSERV    TREE_ID    NODE_ID");
      for (b = 1; b <= _forestSize; b++) {
        for (k = 1; k <= _observationSize; k++) {
	  Rprintf("\n%10d %10d %10d ", k, b, getMembership(root[b], _observation, k) -> leafCount);
        }
      }
      Rprintf("\n");
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nGrown Forest Random Seed Vector:  ");
      Rprintf("\n   TREE_ID                 seed");
      for (b = 1; b <= _forestSize; b++) {
        Rprintf("\n%10d %20d ", b, _seed_[b]);
      }
      Rprintf("\n");
    }
  }
  unstackPreDefinedCommonArrays(oobSampleSize);
  if (mode == RSF_GROW) {
    unstackPreDefinedGrowthArrays(
                                  masterSplit,
                                  masterSplitSize,
                                  masterSplitOrder
                                 );
  }
  else {
    unstackPreDefinedPredictArrays();
  }
  unstackDefinedOutputObjects(mode,
                       root,
                       sortedTimeInterestSize,
                       oobEnsemblePtr,
                       fullEnsemblePtr,
                       ensembleRun,
                       ensembleDen,
                       vimpEnsembleRun);
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Native code rsf() exit. \n");
    now = updateTimeStamp(start);
  }
  UNPROTECT(stackCount + 2);
  return sexpVector[RSF_OUTP_ID];
}
