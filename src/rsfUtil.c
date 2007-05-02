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
#include   "node_ops.h"
#include   "rsfUtil.h"
uint updateTimeStamp(uint before) {
  uint stamp;
  double cpuTimeUsed;
  stamp = clock();
  cpuTimeUsed = ((double) (stamp - before)) / CLOCKS_PER_SEC;
  Rprintf("\nRSF:  CPU process time:  %20.3f \n", cpuTimeUsed);
  return stamp;
}
uint getTraceFlag() {
  uint result;
  if (_traceFlagToggler == 1) {
    result = _traceFlagDiagLevel;
  }
  else {
    result = 0;
  }
  return result;
}
void updateTraceFlag(char reset) {
  if (reset == TRUE) {
    _traceFlagToggler = 1;
  }
  else {
    if (_traceFlagToggler == 1) {
      _traceFlagToggler = _traceFlagIterValue;
    }
    else {
      _traceFlagToggler --;
    }
  }
}
void stackPreDefinedCommonArrays(uint **p_oobSampleSize) {
  uint i;
  _observation = pdvector(1, _xSize);
  _nodeMembership = nodePtrVector(1, _observationSize);
  _bootMembershipIndex = uivector(1, _observationSize);
  _bootMembershipFlag = cvector(1, _observationSize);
  _masterTime  = dvector(1, _observationSize);
  _masterTimeIndex  = uivector(1, _observationSize);
  *p_oobSampleSize = uivector(1, _forestSize);
  for (i = 1; i <= _forestSize; i++) {
    (*p_oobSampleSize)[i] = 0;
  }
}
void unstackPreDefinedCommonArrays(uint *oobSampleSize) {
  free_pdvector(_observation, 1, _xSize);
  free_nodePtrVector(_nodeMembership, 1, _observationSize);
  free_uivector(_bootMembershipIndex, 1, _observationSize);
  free_cvector(_bootMembershipFlag, 1, _observationSize);
  free_dvector(_masterTime, 1, _observationSize);
  free_uivector(_masterTimeIndex, 1, _observationSize);
  free_uivector(oobSampleSize, 1, _forestSize);
}
void stackPreDefinedPredictArrays() {
  _fobservation = pdvector(1, _xSize);
  _fnodeMembership = nodePtrVector(1, _fobservationSize);
}
void unstackPreDefinedPredictArrays() {
  free_pdvector(_fobservation, 1, _xSize);
  free_nodePtrVector(_fnodeMembership, 1, _fobservationSize);
}
void stackPreDefinedGrowthArrays(double ***p_masterSplit,
                                 uint **p_masterSplitSize,
                                 uint ***p_masterSplitOrder) {
  *p_masterSplit = dmatrix(1, _xSize, 1, _observationSize);
  *p_masterSplitSize = uivector(1, _xSize);
  if (_splitRule == LOG_RANK_SCORE) {
    *p_masterSplitOrder = uimatrix(1, _xSize, 1, _observationSize);
  }
}
void unstackPreDefinedGrowthArrays(double **masterSplit,
                                   uint *masterSplitSize,
                                   uint **masterSplitOrder) {
  free_dmatrix(masterSplit, 1, _xSize, 1, _observationSize);
  free_uivector(masterSplitSize, 1, _xSize);
  if (_splitRule == LOG_RANK_SCORE) {
    free_uimatrix(masterSplitOrder, 1, _xSize, 1, _observationSize);
  }
}
uint stackDefinedOutputObjects(uint      mode,
                               uint      sortedTimeInterestSize,
                               char    **sexpString,
                               Node   ***p_root,
                               double ***p_oobEnsemblePtr,
                               double ***p_fullEnsemblePtr,
                               double  **p_ensembleRun,
                               uint    **p_ensembleDen,
                               double  **p_oobEnsemble,
                               double  **p_fullEnsemble,
                               double  **p_performance,
                               uint    **p_leafCount,
                               uint    **p_proximity,
                               double  **p_varImportance,
                               double ***p_vimpEnsembleRun,
                               int     **p_seed,
                               SEXP     *sexpVector) {
  uint sexpIndex;
  uint obsSize;
  uint ensembleSize;
  uint proximitySize;
  uint i,j;
  proximitySize = 0;  
  stackCount = 3;
  if (mode == RSF_GROW) {
    obsSize = _observationSize;
    stackCount += 1;
  }
  else {
    obsSize = _fobservationSize;
  }
  ensembleSize = sortedTimeInterestSize * obsSize;
  if (_mup & MUP_PROX) {
    proximitySize = ((obsSize + 1)  * obsSize) / 2; 
    stackCount += 1;
  }
  if (_mup & MUP_VIMP) {
    stackCount += 1;
  }
  if (_mup & MUP_TREE) stackCount += 5;
  PROTECT(sexpVector[RSF_OUTP_ID] = allocVector(VECSXP, stackCount));
  PROTECT(sexpVector[RSF_STRG_ID] = allocVector(STRSXP, stackCount));
  setAttrib(sexpVector[RSF_OUTP_ID], R_NamesSymbol, sexpVector[RSF_STRG_ID]);
  PROTECT(sexpVector[RSF_FENS_ID] = NEW_NUMERIC(ensembleSize));
  PROTECT(sexpVector[RSF_PERF_ID] = NEW_NUMERIC(_forestSize));
  PROTECT(sexpVector[RSF_LEAF_ID] = NEW_INTEGER(_forestSize));
  *p_fullEnsemble = NUMERIC_POINTER(sexpVector[RSF_FENS_ID]);
  *p_performance  = NUMERIC_POINTER(sexpVector[RSF_PERF_ID]);
  *p_leafCount    = (uint*) INTEGER_POINTER(sexpVector[RSF_LEAF_ID]);
  SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], 0, sexpVector[RSF_FENS_ID]);
  SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], 1, sexpVector[RSF_PERF_ID]);
  SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], 2, sexpVector[RSF_LEAF_ID]);
  SET_STRING_ELT(sexpVector[RSF_STRG_ID], 0, mkChar(sexpString[RSF_FENS_ID]));
  SET_STRING_ELT(sexpVector[RSF_STRG_ID], 1, mkChar(sexpString[RSF_PERF_ID]));
  SET_STRING_ELT(sexpVector[RSF_STRG_ID], 2, mkChar(sexpString[RSF_LEAF_ID]));
  sexpIndex = 3;
  if (mode == RSF_GROW) {
    PROTECT(sexpVector[RSF_OENS_ID] = NEW_NUMERIC(ensembleSize));
    *p_oobEnsemble = NUMERIC_POINTER(sexpVector[RSF_OENS_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_OENS_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_OENS_ID]));
    sexpIndex++;
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for GROW mode information.");  
    }
  }
  if (_mup & MUP_PROX) {
    PROTECT(sexpVector[RSF_PROX_ID] = NEW_INTEGER(proximitySize));
    *p_proximity = (uint*) INTEGER_POINTER(sexpVector[RSF_PROX_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PROX_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_PROX_ID]));
    sexpIndex++;
    (*p_proximity) --;
    for (i = 1; i <= proximitySize; i++) {
      (*p_proximity)[i] = 0;
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for MUP proximity information.");  
    }
  }
  if (_mup & MUP_VIMP) {
    PROTECT(sexpVector[RSF_VIMP_ID] = NEW_NUMERIC(_xSize));
    *p_varImportance = NUMERIC_POINTER(sexpVector[RSF_VIMP_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_VIMP_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_VIMP_ID]));
    sexpIndex++;
    (*p_varImportance) --;
    (*p_vimpEnsembleRun) = dmatrix(1, _xSize, 1, _observationSize);
    for (i = 1; i <= _xSize; i++) {
      (*p_varImportance)[i] = 0;
      for (j = 1; j <= _observationSize; j++) {
        (*p_vimpEnsembleRun)[i][j] = 0;
      }
    }
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for MUP variable importance.");  
    }
  }
  if (_mup & MUP_TREE) {
    *p_root = nodePtrVector(1, _forestSize);
    PROTECT(sexpVector[RSF_SEED_ID] = NEW_INTEGER(_forestSize));
    *p_seed = INTEGER_POINTER(sexpVector[RSF_SEED_ID]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_SEED_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex, mkChar(sexpString[RSF_SEED_ID]));
    sexpIndex++;
    (*p_seed) --;
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nAllocating for MUP forest information.");  
    }
  }
  *p_ensembleRun = dvector(1, obsSize);
  *p_ensembleDen = uivector(1, obsSize);
  *p_fullEnsemblePtr = pdvector(1, sortedTimeInterestSize);
  if (mode == RSF_GROW) {
    *p_oobEnsemblePtr = pdvector(1, sortedTimeInterestSize);
    for (i = 1; i <= sortedTimeInterestSize; i++) {
      (*p_oobEnsemblePtr)[i]  = (*p_oobEnsemble)  + ((i-1)*(obsSize)) - 1;
      (*p_fullEnsemblePtr)[i] = (*p_fullEnsemble) + ((i-1)*(obsSize)) - 1;
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= sortedTimeInterestSize; j++) {
        (*p_oobEnsemblePtr)[j][i] = 0.0;
        (*p_fullEnsemblePtr)[j][i] = 0.0;
      }
      (*p_ensembleDen)[i] = 0;
    }
  }
  if (mode == RSF_PRED) {
    for (i = 1; i <= sortedTimeInterestSize; i++) {
      (*p_fullEnsemblePtr)[i] = (*p_fullEnsemble) + ((i-1)*(obsSize)) - 1;
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= sortedTimeInterestSize; j++) {
        (*p_fullEnsemblePtr)[j][i] = 0.0;
      }
      (*p_ensembleDen)[i] = _forestSize;
    }
  }
  (*p_performance) --;
  (*p_leafCount)   --;
  for (i = 1; i <= _forestSize; i++) {
    (*p_leafCount)[i] = 1;
    (*p_performance)[i] = 0;
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Allocation of defined output objects complete:  %10d", sexpIndex);  
  }
  return (sexpIndex);
}
uint stackVariableOutputObjects(uint     sexpIndex,
                                uint     forestNodeCount,
                                char   **sexpString,
                                uint   **p_treeID,
                                uint   **p_nodeID,
                                uint   **p_parmID,                                   
                                double **p_spltPT,
                                SEXP    *sexpVector) {
  if (_mup & MUP_TREE) {
    PROTECT(sexpVector[RSF_TREE_ID] = NEW_INTEGER(forestNodeCount));
    PROTECT(sexpVector[RSF_NODE_ID] = NEW_INTEGER(forestNodeCount));
    PROTECT(sexpVector[RSF_PARM_ID] = NEW_INTEGER(forestNodeCount));
    PROTECT(sexpVector[RSF_SPLT_PT] = NEW_NUMERIC(forestNodeCount));    
    *p_treeID = (uint*) INTEGER_POINTER(sexpVector[RSF_TREE_ID]);
    *p_nodeID = (uint*) INTEGER_POINTER(sexpVector[RSF_NODE_ID]);
    *p_parmID = (uint*) INTEGER_POINTER(sexpVector[RSF_PARM_ID]);
    *p_spltPT = NUMERIC_POINTER(sexpVector[RSF_SPLT_PT]);
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_TREE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_TREE_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_NODE_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_NODE_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_PARM_ID]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_PARM_ID]));
    SET_VECTOR_ELT(sexpVector[RSF_OUTP_ID], sexpIndex, sexpVector[RSF_SPLT_PT]);
    SET_STRING_ELT(sexpVector[RSF_STRG_ID], sexpIndex++, mkChar(sexpString[RSF_SPLT_PT]));
    (*p_treeID) --;
    (*p_nodeID) --;
    (*p_parmID) --;
    (*p_spltPT) --;
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Allocation of variable output objects complete:  %10d", sexpIndex);  
  }
  return (sexpIndex);
}
void unstackDefinedOutputObjects(uint      mode,
                          Node    **root,
                          uint      sortedTimeInterestSize,
                          double  **oobEnsemblePtr,
                          double  **fullEnsemblePtr,
                          double   *ensembleRun,
                          uint     *ensembleDen,
                          double  **vimpEnsembleRun) {
  uint i;
  if (_mup & MUP_TREE) {
    for (i = 1; i <= _forestSize; i++) {
      freeTree(root[i]);
    }
    free_nodePtrVector(root, 1, _forestSize);
    if (getTraceFlag() & DL1_TRACE) {
      Rprintf("\nDe-allocating for forest information.");  
    }
  }
  free_pdvector(fullEnsemblePtr, 1, sortedTimeInterestSize);
  if (mode == RSF_GROW) {
    free_dvector(ensembleRun, 1, _observationSize);
    free_uivector(ensembleDen, 1, _observationSize);
    free_pdvector(oobEnsemblePtr, 1, sortedTimeInterestSize);
  }
  else {
    free_dvector(ensembleRun, 1, _fobservationSize);
    free_uivector(ensembleDen, 1, _fobservationSize);
  }
  if (_mup & MUP_VIMP) {
    free_dmatrix(vimpEnsembleRun, 1, _xSize, 1, _observationSize);
  }
}
void initializeArrays(uint mode,
                      uint *masterTimeSize,
                      uint *sortedTimeInterestSize,
                      uint *masterDeathTimeSize
                     ) {
  uint i,j,k;
  uint leadingIndex;
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nCommon Incoming Parameters:  ");
    Rprintf("\n         traceFlag:  %10d", getTraceFlag());
    Rprintf("\n               mup:  %10d", _mup);
    Rprintf("\n        forestSize:  %10d", _forestSize);
    Rprintf("\n   observationSize:  %10d", _observationSize);
    Rprintf("\n  timeInterestSize:  %10d", _timeInterestSize);
    Rprintf("\n             xSize:  %10d", _xSize);
    Rprintf("\n");
  }
  if ((getTraceFlag() & RSF_GROW) && (getTraceFlag() & DL1_TRACE)) {
    Rprintf("\nIncoming Random Covariate Weights:  ");
    Rprintf("\n     index       weight");
    for (j=1; j <= _xSize; j++) {
      Rprintf("\n%10d  %12.4f", j, _randomCovariateWeight[j]);
    }
    Rprintf("\n");
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nIncoming GROW Data:  ");
    Rprintf("\n     index         time     status   observations -> \n");
  }
  for (i = 1; i <= _observationSize; i++) {
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("%10d %12.4f %10d", i, _time[i], _status[i]);
      for (j=1; j <= _xSize; j++) {
        Rprintf(" %12.4f", (_xData+((j-1)*(_observationSize)))[i-1]);
      }
      Rprintf("\n");
    }
    _masterTime[i] = _time[i];
  }
  for (j=1; j <= _xSize; j++) {
    _observation[j] = (_xData + ((j-1)*(_observationSize)) - 1);
  }
  if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
    Rprintf("\nIncoming PRED Data:  ");
    Rprintf("\n     index         time     status   observations -> \n");
    for (i = 1; i <= _fobservationSize; i++) {
      Rprintf("%10d %12.4f %10d", i, _ftime[i], _fstatus[i]);
      for (j=1; j <= _xSize; j++) {
        Rprintf(" %12.4f", (_fxData+((j-1)*(_fobservationSize)))[i-1]);
      }
      Rprintf("\n");
    }
  }
  if (mode == RSF_PRED) {
    for (j=1; j <= _xSize; j++) {
      _fobservation[j] = (_fxData + ((j-1)*(_fobservationSize)) - 1);
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Initial read complete.");  
  }
  hpsort(_masterTime, _observationSize);
  *masterTimeSize = leadingIndex = 1;
  for (i=2; i <= _observationSize; i++) {
    if (_masterTime[i] > _masterTime[leadingIndex]) {
      (*masterTimeSize) ++;
      leadingIndex++;
      _masterTime[leadingIndex] = _masterTime[i];
    }
  }
  for (i= (*masterTimeSize) + 1; i <= _observationSize; i++) {
    _masterTime[i] = 0;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nSorted Distinct Event Times:  \n");
    for (i=1; i <= *masterTimeSize; i++) {
      Rprintf("%10d %10.4f \n", i, _masterTime[i]);
    }
  }
  for (j=1; j <= _observationSize; j++) {
    k = 1;
    while (k <= *masterTimeSize) {
      if (_time[j] == _masterTime[k]) {
        _masterTimeIndex[j] = k;
        k = *masterTimeSize;
      }
      k++;
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nMaster Time Index:  \n");
    for (i=1; i <= _observationSize; i++) {
      Rprintf("%10d %10d \n", i, _masterTimeIndex[i]);
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of master time data complete.");  
  }
  hpsort(_timeInterest, _timeInterestSize);
  *sortedTimeInterestSize = leadingIndex = 1;
  for (i=2; i <= _timeInterestSize; i++) {
    if (_timeInterest[i] > _timeInterest[leadingIndex]) {
      (*sortedTimeInterestSize) ++;
      leadingIndex++;
      _timeInterest[leadingIndex] = _timeInterest[i];
    }
  }
  if (*sortedTimeInterestSize != _timeInterestSize) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Time points of interest are not unique.");
    Rprintf("\nRSF:  The ensemble estimate output matrix is being");
    Rprintf("\nRSF:  resized as [N'] x [n], where N' is the");
    Rprintf("\nRSF:  unique time points of interest and n is");
    Rprintf("\nRSF:  number of observations in the data.");
  }
  for (i= (*sortedTimeInterestSize) + 1; i <= _timeInterestSize; i++) {
    _timeInterest[i] = 0;
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nSorted Distinct Times of Interest:  \n");
    for (i=1; i <= *sortedTimeInterestSize; i++) {
      Rprintf("%10d %10.4f \n", i, _timeInterest[i]);
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of time interest data complete.");  
  }
  char *masterDeathTimeIndicator = cvector(1, *masterTimeSize);
  for (i=1; i <= *masterTimeSize; i++) {
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
    for (i=1; i <= *masterTimeSize; i++) {
      masterDeathTimeIndicator[i] = TRUE;
    }
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  All observations were censored, and will be considered to be deaths.");
  }
  *masterDeathTimeSize = 0;
  for (i=1; i<= *masterTimeSize; i++) {
    if (masterDeathTimeIndicator[i] == TRUE) {
      (*masterDeathTimeSize) ++;
    }
  }
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\n\nMaster Death Time Indicator:  \n");
    for (i=1; i <= *masterTimeSize; i++) {
      Rprintf("%10d %10.4f %10d \n", i, _masterTime[i], masterDeathTimeIndicator[i]);
    }
  }
  free_cvector(masterDeathTimeIndicator, 1, *masterTimeSize);
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Initialization of death time data complete.");  
  }
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
    exit(TRUE);
  }
  if ((getTraceFlag() & RSF_PRED) && (getTraceFlag() & DL2_TRACE)) {
    Rprintf("%10d %10d %10d %12.4f \n", treeID[*offset], parmID[*offset], nodeID[*offset], spltPT[*offset]);
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
    restoreTree(b, parent ->  left, leafCount, offset, treeID, nodeID, parmID, spltPT);
    parent -> right = makeNode();
    ((parent -> right) -> leafCount) = *leafCount;
    restoreTree(b, parent -> right, leafCount, offset, treeID, nodeID, parmID, spltPT);
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
void getCumulativeHazardEstimate(double **cumulativeHazard,
                                 Node *parent,
                                 uint sortedTimeInterestSize,
                                 uint masterTimeSize) {
  uint i, j;
  uint *nodeParentDeath  = uivector(1, masterTimeSize);
  uint *nodeParentAtRisk = uivector(1, masterTimeSize);
  uint priorTimePointIndex, currentTimePointIndex;
  double estimate;
  for (i=1; i <= masterTimeSize; i++) {
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
    for (i = priorTimePointIndex + 1; i <= masterTimeSize; i++) {
      if (_masterTime[i] <= _timeInterest[j]) {
        currentTimePointIndex = i;
      }
      else {
        i = masterTimeSize;
      }
    }
    if (getTraceFlag() & DL3_TRACE) {
      Rprintf("\nCumulative hazard at risk and death counts for node:  %10d \n", parent -> leafCount);
      Rprintf("  with time point index:  %10d \n", currentTimePointIndex);
      for (i=1; i <= masterTimeSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (i=1; i <= masterTimeSize; i++) {
        Rprintf("%10d", nodeParentAtRisk[i]);
      }
      Rprintf("\n");
      for (i=1; i <= masterTimeSize; i++) {
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
  free_uivector(nodeParentDeath, 1, masterTimeSize);
  free_uivector(nodeParentAtRisk, 1, masterTimeSize);
}
double getConcordanceIndex(uint size, 
                           double *timePtr, 
                           uint *statusPtr, 
                           double *mortality,
                           uint *oobCount
                          ) {
  uint i,j;
  ulong concordancePairSize;
  ulong concordanceWorseCount;
  concordancePairSize = concordanceWorseCount = 0;
  for (i=1; i < size; i++) {
    for (j=i+1; j <= size; j++) {
      if (oobCount[i] != 0  && oobCount[j] != 0) {
        if ( (timePtr[i] >  timePtr[j] && statusPtr[j] == 1) ||
             (timePtr[i] == timePtr[j] && statusPtr[j] == 1 && statusPtr[i] == 0) ) {
          concordancePairSize += 2;
          if (mortality[j] > mortality[i]) {              
            concordanceWorseCount += 2;
          }
          else if (fabs(mortality[j] - mortality[i]) < EPSILON) {              
            concordanceWorseCount += 1;
          }  
        }
        else if ( (timePtr[j] >  timePtr[i] && statusPtr[i] == 1) ||
                  (timePtr[j] == timePtr[i] && statusPtr[i] == 1 && statusPtr[j] == 0) ) {
          concordancePairSize += 2;
          if (mortality[i] > mortality[j]) {
            concordanceWorseCount += 2;
          }
          else if (fabs(mortality[i] - mortality[j]) < EPSILON) {
            concordanceWorseCount += 1;
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
  return ((double) concordanceWorseCount / (double) concordancePairSize);
}
double getPerformance(uint     mode,
                      uint     sortedTimeInterestSize,
                      uint     leafCount,
                      uint     masterTimeSize,
                      double **oobEnsemblePtr,
                      double **fullEnsemblePtr,
                      double  *ensembleRun,
                      uint    *ensembleDen,
                      uint     oobSampleSize,
                      Node    *rootPtr,
                      double **vimpEnsembleRun) {
  uint i,j,k,p;
  uint obsSize;
  double **cumulativeHazard;
  uint *indexOOB;
  char *randomOOB;
  uint *permuteOOB;
  double *originalOOB;
  double *timePtr;
  uint   *statusPtr;
  double performance;
  k = 0;  
  performance = 0;
  cumulativeHazard = dmatrix(1, sortedTimeInterestSize, 1, leafCount);
  if (mode == RSF_GROW) {
    obsSize = _observationSize;
  }
  else {
    obsSize = _fobservationSize;
  }
  for (i=1; i <= leafCount; i++) {
    for (j = 1; j <= _observationSize; j++) {
      if (_bootMembershipFlag[j] == TRUE) {
        if ((_nodeMembership[j] -> leafCount) == i) {
          k = j;
          j = _observationSize;
        }
      }
    }
    getCumulativeHazardEstimate(cumulativeHazard,
                                _nodeMembership[k],
                                sortedTimeInterestSize,
                                masterTimeSize
                               );
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
  if (mode == RSF_PRED) {
    for (i=1; i <= _fobservationSize; i++) {
      k = _fnodeMembership[i] -> leafCount;
      for (j=1; j <= sortedTimeInterestSize; j++) {
        fullEnsemblePtr[j][i] += cumulativeHazard[j][k];
      }
    }
  }
  if (_mup & MUP_VIMP) {
    indexOOB = uivector(1, oobSampleSize);
    randomOOB = cvector(1, oobSampleSize);
    permuteOOB = uivector(1, oobSampleSize);
    originalOOB = dvector(1, oobSampleSize);
    k = 1;
    for (i=1; i <= _observationSize; i++) {
      if ( _bootMembershipFlag[i] == FALSE ) {
        indexOOB[k] = i;
        k++;
      }
    }
    for (p=1; p <= _xSize; p++) {
      for (k=1; k<= oobSampleSize; k++) {
        originalOOB[k] = _observation[p][indexOOB[k]];
        randomOOB[k] = TRUE;
      }
      for (i=oobSampleSize; i > 0; i--) {
        k = ceil(ran1(_seedPtr)*(i*1.0));
        for (j = 1; k > 0; j++) {
          if (randomOOB[j] == TRUE) {
            k--;
          }
        }
        randomOOB[j-1] = ACTIVE;
        permuteOOB[i] = j-1;
      }
      if (getTraceFlag() & DL2_TRACE) {
        Rprintf("\nCovariate Permutation:  %10d\n", p);
          Rprintf("     index   indexOOB    permOOB\n");
        for (i=1; i <= oobSampleSize; i++) {
          Rprintf("%10d %10d %10d \n", i, indexOOB[i], permuteOOB[i]);
        }
      }
      for (k=1; k <= oobSampleSize; k++) {
        _observation[p][indexOOB[k]] = originalOOB[permuteOOB[k]];
      }
      for (i=1; i <= _observationSize; i++) {
          if ( _bootMembershipFlag[i] == FALSE ) {
          k = getMembership(rootPtr, _observation, i) -> leafCount;
          for (j=1; j <= sortedTimeInterestSize; j++) {
            vimpEnsembleRun[p][i] += cumulativeHazard[j][k];
          }
        }
      }
      for (k=1; k <= oobSampleSize; k++) {
        _observation[p][indexOOB[k]] = originalOOB[k];
      }
    }  
    if (getTraceFlag() & DL2_TRACE) {
      Rprintf("\nTree specific variable importance calculation: \n");
      Rprintf("          ");
      for (i=1; i <= _observationSize; i++) {
        Rprintf("%10d", i);
      }
      Rprintf("\n");
      for (p=1; p <= _xSize; p++) {
        Rprintf("%10d", p);
        for (i=1; i <= _observationSize; i++) {
          Rprintf("%10.4f", vimpEnsembleRun[p][i]);
        }
        Rprintf("\n");
      }
    }
    free_uivector(indexOOB, 1, oobSampleSize);
    free_cvector(randomOOB, 1, oobSampleSize);
    free_uivector(permuteOOB, 1, oobSampleSize);
    free_dvector(originalOOB, 1, oobSampleSize);
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
    for (i = 1; i <= _observationSize; i++) {
      ensembleRun[i] = 0.0;
      for (j=1; j <= sortedTimeInterestSize; j++) {
        ensembleRun[i] += oobEnsemblePtr[j][i];
      }
      if (ensembleDen[i] != 0) {
        ensembleRun[i] = ensembleRun[i] / ensembleDen[i];
      }
    }
  }
  if (mode == RSF_PRED) {
    for (i = 1; i <= _fobservationSize; i++) {
      ensembleRun[i] = 0.0;
      for (j=1; j <= sortedTimeInterestSize; j++) {
        ensembleRun[i] += fullEnsemblePtr[j][i];
      }
      if (ensembleDen[i] != 0) {
        ensembleRun[i] = ensembleRun[i] / _forestSize;
      }
    }
  }
  if (getTraceFlag() & DL0_TRACE) {
    Rprintf("\nRSF:  Bootstrap sample ensemble estimate complete.");  
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
    Rprintf("\nDiagnostic Classification of OOB elements:  ");
    Rprintf("\n       Leaf     Observ         EE    Data ->  ");
    for (i=1; i <= leafCount; i++) {
      for (j = 1; j <= _observationSize; j++) {
	if (_bootMembershipFlag[j] == FALSE) {
	  if ((_nodeMembership[j] -> leafCount) == i) {
	    Rprintf("\n %10d %10d %10.4f", i, j, ensembleRun[j]);
	    for (k = 1; k <= _xSize; k++) {
	      Rprintf("%10.4f", _observation[k][j]);
	    }
	  }
	}
      }
    }
    Rprintf("\n");
  }
  if (_mup & MUP_PERF) {
    if (mode == RSF_GROW) {
      timePtr = _time;
      statusPtr = _status;
    }
    else {
      timePtr = _ftime;
      statusPtr = _fstatus;
    }
    performance = 1.0 - getConcordanceIndex(obsSize, timePtr, statusPtr, ensembleRun, ensembleDen);
    if (getTraceFlag() & DL0_TRACE) {
      Rprintf("\nRSF:  Error Rate:                         %20.4f", performance);
      Rprintf("\n");
    }
  }
  return performance;
}
