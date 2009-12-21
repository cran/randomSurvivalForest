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
#include       "nodeOps.h"
#include  "rsfFactorOps.h"
#include     "rsfImpute.h"
#include      "rsfSplit.h"
#include       "rsfTree.h"
Node *getTerminalNode(uint mode, uint leaf) {
  uint i, j;
  Node *parent;
  parent = NULL;
  for (j = 1; j <= _observationSize; j++) {
    if ((_nodeMembership[j] -> leafCount) == leaf) {
      parent = _nodeMembership[j];
      j = _observationSize;
    }
  }
  if (!(_opt & OPT_VOUT_TYPE) && (parent == NULL)) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Proxy member for node %12d not found.", leaf);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    Rprintf("\nDiagnostic Trace of (individual, boot, node, leaf) vectors in data set:  ");
    Rprintf("\n        index         boot         node         leaf \n");
    for (i = 1; i <= _observationSize; i++) {
      Rprintf(" %12d %12d %12x %12d \n", i, 
              _bootMembershipFlag[i], _nodeMembership[i], 
              _nodeMembership[i] -> leafCount);
    }
    exit(TRUE);
  }
  return parent;
}
char testNodeSize(Node *parent) {
  uint localDeathTimeSize;
  char result;
  uint i;
  uint *localDeathTimeCount = uivector(1, _masterTimeSize);
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\ntestNodeSize() ENTRY ...\n");
  }
  localDeathTimeSize = 0;
  for (i=1; i <= _masterTimeSize; i++) {
    localDeathTimeCount[i] = 0;
  }
  for (i=1; i <= _observationSize; i++) {
    if (_nodeMembership[_bootMembershipIndex[i]] == parent) {
      if (_status[_bootMembershipIndex[i]] > 0) {
        if (_mRecordMap[_bootMembershipIndex[i]] == 0) {
          localDeathTimeCount[_masterTimeIndex[_bootMembershipIndex[i]]] ++;
        }
        else if ((_mvSign[abs(CENS_IDX)][_mRecordMap[_bootMembershipIndex[i]]] == 0) &&
                 (_mvSign[abs(TIME_IDX)][_mRecordMap[_bootMembershipIndex[i]]] == 0)) {
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
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\ntestNodeSize(%1d) EXIT ...\n", result);
  }
    return result;
}
char forkAndUpdate(uint  *leafCount,
                   Node  *parent,
                   uint   splitParameterMax,
                   uint   splitValueMaxFactSize,
                   uint  *splitValueMaxFactPtr,
                   double splitValueMaxCont) {
  char factorFlag;
  char daughterFlag;
  uint i;
  char result;
  if (getTraceFlag() & FORK_DEF_TRACE) {
    Rprintf("\nforkAndUpdate() ENTRY ...\n");
  }
  result = forkNode(parent, 
                    splitParameterMax, 
                    splitValueMaxFactSize, 
                    splitValueMaxFactPtr, 
                    splitValueMaxCont);
  if (result == TRUE) {
    if (getTraceFlag() & FORK_DEF_TRACE) {
      Rprintf("\nForking On:  ");
      Rprintf("\n   LeafCnt   SpltParm  ");
      Rprintf("\n%10d %10d ", parent -> leafCount, parent -> splitParameter);
      if (parent -> splitValueFactSize > 0) {
        for (i = parent -> splitValueFactSize; i >= 1; i--) {
          Rprintf("%8x ", (parent -> splitValueFactPtr)[i]);
        }
        Rprintf("\n");
      }
      else {
        Rprintf(" %20.4f \n", parent -> splitValueCont);
      }
    }
    (*leafCount)++;
    factorFlag = FALSE;
    if (strcmp(_xType[splitParameterMax], "C") == 0) {
      factorFlag = TRUE;
      _totalMWCPCount += parent -> splitValueFactSize;
    }
    for (i = 1; i <= _observationSize; i++) {
      if (_nodeMembership[i] == parent) {
        daughterFlag = RIGHT;
        if (factorFlag == TRUE) {
          if (_observation[splitParameterMax][i] != 0) {
            daughterFlag = splitOnFactor((uint) _observation[splitParameterMax][i], _splitValueMaxFactPtr);
            if (getTraceFlag() & FORK_DEF_TRACE) {
              Rprintf("\nNode Membership:  %10d %16d", i, (uint) _observation[splitParameterMax][i]);
            }
          }
          else {
            Rprintf("\nRSF:  *** ERROR *** ");
            Rprintf("\nRSF:  Attempt to fork on NA value on (index, parameter):  (%10d, %10d)", i, splitParameterMax);
            Rprintf("\nRSF:  Please Contact Technical Support.");
            Rprintf("\nRSF:  The application will now exit.\n");
            exit(TRUE);
          }
        }
        else {
          if (!ISNA(_observation[splitParameterMax][i])) {
            if (getTraceFlag() & FORK_DEF_TRACE) {
              Rprintf("\nNode Membership:  %10d %16.4f, %16.4f", i, _observation[splitParameterMax][i], _splitValueMaxCont);
            }
            if (_observation[splitParameterMax][i] <= _splitValueMaxCont) {
              daughterFlag = LEFT;
            }
          }
          else {
            Rprintf("\nRSF:  *** ERROR *** ");
            Rprintf("\nRSF:  Attempt to fork on NA value on (index, parameter):  (%10d, %10d)", i, splitParameterMax);
            Rprintf("\nRSF:  Please Contact Technical Support.");
            Rprintf("\nRSF:  The application will now exit.\n");
            exit(TRUE);
          }
        }
        if (daughterFlag == LEFT) {
          _nodeMembership[i] = parent -> left;
          ((parent -> left) -> leafCount) = (parent -> leafCount);
          if (getTraceFlag() & FORK_DEF_TRACE) {
            Rprintf(" --> LEFT ");
          }
        }
        else {
          _nodeMembership[i] = parent -> right;
          ((parent -> right) -> leafCount) = *leafCount;
          if (getTraceFlag() & FORK_DEF_TRACE) {
            Rprintf(" --> RGHT ");
          }
        }
      }
    }
  }
  else {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  forkNode() failed.");
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    exit(TRUE);
  }
  if (getTraceFlag() & FORK_DEF_TRACE) {
    Rprintf("\nforkAndUpdate() EXIT ...\n");
  }
  return result;
}
char makeTree (char     multipleImputeFlag,
               uint     b,
               Node    *parent,
               uint     depth,
               uint    *maximumDepth) {
  char result;
  uint splitParameterMax;
  Node *reversePtr;
  uint i;
  if (getTraceFlag() & SPLT_DEF_TRACE) {
    Rprintf("\nmakeTree() ENTRY ...\n");
  }
  result = TRUE;
  if (getTraceFlag() & SPLT_DEF_TRACE) {
    Rprintf("\nMake tree leaf:  %10d at depth %10d \n", parent -> leafCount, parent -> depth);
  }
  parent -> depth = depth;
  if (multipleImputeFlag == FALSE) {
    if (_mRecordSize > 0) {
      result = (testNodeSize(parent) || (_leafCount_[b] == 1));
      if (result) {
        imputeNode(RSF_GROW,
                   TRUE,
                   b,
                   parent);
        if (_mTimeIndexFlag == TRUE) {
          updateTimeIndexArray(parent);
        }
      }
    }
  }
  if (result == TRUE) {
    result = getBestSplit(parent, & splitParameterMax);
    if (result == TRUE) {
      result = forkAndUpdate(_leafCount_ + b,
                             parent,
                             splitParameterMax,
                             _splitValueMaxFactSize, 
                             _splitValueMaxFactPtr, 
                             _splitValueMaxCont);
      if (result == TRUE) {
        if (getTraceFlag() & SPLT_DEF_TRACE) {
          Rprintf("\nNode Membership:  \n");
          for (i=1; i <=  _observationSize; i++) {
            Rprintf("%10d %10d \n", i, _nodeMembership[i] -> leafCount);
          }
        }
        if (getTraceFlag() & SPLT_DEF_TRACE) {
          Rprintf("\nmakeTree() LEFT:  \n");
        }
        makeTree (multipleImputeFlag,
                  b,
                  parent -> left,
                  depth + 1,
                  maximumDepth);
        if (getTraceFlag() & SPLT_DEF_TRACE) {
          Rprintf("\nmakeTree() RIGHT:  \n");
        }
        makeTree (multipleImputeFlag,
                  b,
                  parent -> right,
                  depth + 1,
                  maximumDepth);
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
      if (getTraceFlag() & SPLT_DEF_TRACE) {
        Rprintf("\ngetBestSplit() FAILED ...\n");
      }
    }
  }  
  else {
    parent -> splitFlag = FALSE;
    if (getTraceFlag() & SPLT_DEF_TRACE) {
      Rprintf("\ntestMinimumDeath() FAILED ...\n");
    }
  }
  if (!result) {
    if (_eventTypeSize > 1) {
      if (!(_opt & OPT_IMPU_ONLY)) {
        parent -> poe = uivector(1, _eventTypeSize);
        parent -> subSurvival = dmatrix(1, _eventTypeSize, 1, _sortedTimeInterestSize);
      }
    }
  }
  if (!result) {
    *maximumDepth = ((depth > *maximumDepth) ? parent -> depth : *maximumDepth);
    if (depth > 0) {
      parent -> splitDepth = uivector(1, parent -> depth);
      reversePtr = parent;
      for (i = 1; i <= depth; i++) {
        if ((reversePtr -> parent) == NULL) {
          Rprintf("\nRSF:  *** ERROR *** ");
          Rprintf("\nRSF:  Reverse parsing of tree failed in forkAndUpdate().");
          Rprintf("\nRSF:  Please Contact Technical Support.");
          Rprintf("\nRSF:  The application will now exit.\n");
          exit(TRUE);
        }
        (parent -> splitDepth)[depth - i + 1] = (reversePtr -> parent) -> splitParameter;
        reversePtr = reversePtr -> parent;
      }    
    }
  }
  if (getTraceFlag() & SPLT_DEF_TRACE) {
    Rprintf("\nmakeTree() EXIT ...\n");
  }
  return result;
}
char restoreTree(uint    b,
                 Node   *parent,
                 uint   *leafCount,
                 uint   *offset,
                 uint   *treeID,
                 uint   *nodeID,
                 uint   *parmID,
                 double *contPT,
                 uint   *mwcpSZ,
                 uint  **mwcpPtr,
                 uint    depth,
                 uint   *maximumDepth) {
  char result;
  Node *reversePtr;
  uint i;
  if (getTraceFlag() & SPLT_DEF_TRACE) {
    Rprintf("\nrestoreTree() ENTRY ...\n");
  }
  if (b != treeID[*offset]) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Invalid forest input record at line:  %10d", b);
    Rprintf("\nRSF:  Please Contact Technical Support.");
    Rprintf("\nRSF:  The application will now exit.\n");
    Rprintf("\nDiagnostic Trace of Tree Record:  \n");
    Rprintf("\n    treeID     nodeID     parmID       spltPT     mwcpSZ \n");
    Rprintf("%10d %10d %10d %12.4f %10d \n", treeID[*offset], nodeID[*offset], parmID[*offset], contPT[*offset], mwcpSZ[*offset]);
    exit(TRUE);
  }
  parent -> depth = depth;
  parent -> left  = NULL;
  parent -> right = NULL;
  for (i = 1; i <= parent -> xSize; i++) {
    parent -> permissibleSplit[i] = FALSE;
  }
  parent -> splitFlag = FALSE;
  parent -> mortality = NA_REAL;
  parent -> leafCount = nodeID[*offset];
  parent -> splitParameter = parmID[*offset];
  if ((parent -> splitParameter) != 0) {
    if (strcmp(_xType[parent -> splitParameter], "C") == 0) {
      if (getTraceFlag() & SPLT_DEF_TRACE) {
        Rprintf("  mwcpPT (reversed):  [%16x] 0x ", (*mwcpPtr) + 1);
      }
      parent -> splitValueFactSize = mwcpSZ[*offset];
      parent -> splitValueFactPtr = uivector(1, mwcpSZ[*offset]);    
      for (i = 1; i <= parent -> splitValueFactSize; i++) {
        (*mwcpPtr) ++;
        (parent -> splitValueFactPtr)[i] = **mwcpPtr;
        if (getTraceFlag() & SPLT_DEF_TRACE) {
          Rprintf("%8x ",  (parent -> splitValueFactPtr)[i]);
        }
      }
      parent -> splitValueCont = NA_REAL;
    }
    else {
      parent -> splitValueCont = contPT[*offset];
      parent -> splitValueFactSize = 0;
      parent -> splitValueFactPtr = NULL;
    }
  }
  else {
    parent -> splitValueCont     = NA_REAL;
    parent -> splitValueFactSize = 0;
    parent -> splitValueFactPtr  = NULL;
  }
  (*offset) ++;
  if ((parent -> splitParameter) != 0) {
    if (getTraceFlag() & FORK_DEF_TRACE) {
      getNodeInfo(parent);
    }
    result = TRUE;
    (*leafCount)++;
    parent -> left  = makeNode(parent -> xSize);
    setParent(parent ->  left, parent);
    ((parent -> left) -> leafCount) = (parent -> leafCount);
    restoreTree(b, 
                parent -> left, 
                leafCount, 
                offset, 
                treeID, 
                nodeID, 
                parmID,
                contPT,
                mwcpSZ,
                mwcpPtr,
                depth + 1,
                maximumDepth);
    parent -> right = makeNode(parent -> xSize);
    setParent(parent -> right, parent);
    ((parent -> right) -> leafCount) = *leafCount;
    restoreTree(b, 
                parent -> right, 
                leafCount, 
                offset, 
                treeID, 
                nodeID, 
                parmID,
                contPT,
                mwcpSZ,
                mwcpPtr,
                depth + 1,
                maximumDepth);
  }
  else {
    result = FALSE;
  }
  if (!result) {
    if (_eventTypeSize > 1) {
      if (!(_opt & OPT_IMPU_ONLY)) {
        parent -> poe = uivector(1, _eventTypeSize);
        parent -> subSurvival = dmatrix(1, _eventTypeSize, 1, _sortedTimeInterestSize);
      }
    }
  }
  if (!result) {
    *maximumDepth = ((depth > *maximumDepth) ? parent -> depth : *maximumDepth);
    if (depth > 0) {
      parent -> splitDepth = uivector(1, parent -> depth);
      reversePtr = parent;
      for (i = 1; i <= depth; i++) {
        if ((reversePtr -> parent) == NULL) {
          Rprintf("\nRSF:  *** ERROR *** ");
          Rprintf("\nRSF:  Reverse parsing of tree failed in restoreTree().");
          Rprintf("\nRSF:  Please Contact Technical Support.");
          Rprintf("\nRSF:  The application will now exit.\n");
          exit(TRUE);
        }
        (parent -> splitDepth)[depth - i + 1] = (reversePtr -> parent) -> splitParameter;
        reversePtr = reversePtr -> parent;
      }    
    }
  }
  if (getTraceFlag() & SPLT_DEF_TRACE) {
    Rprintf("\nrestoreTree() EXIT ...\n");
  }
  return result;
}
void saveTree(uint    b,
              Node   *parent,
              uint   *offset,
              uint   *treeID,
              uint   *nodeID,
              uint   *parmID,
              double *contPT,
              uint   *mwcpSZ,
              uint  **mwcpPtr) {
  uint i;
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\nsaveTree() ENTRY ...\n");
  }
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    getNodeInfo(parent);
  }
  treeID[*offset] = b;
  nodeID[*offset] = parent -> leafCount;
  parmID[*offset] = parent -> splitParameter;
  if ((parent -> splitParameter) != 0) {
    if (strcmp(_xType[parent -> splitParameter], "C") == 0) {
      if (getTraceFlag() & SUMM_HGH_TRACE) {
        Rprintf("  mwcpPT (reversed):  [%16x] 0x ", *(*mwcpPtr + 1));
      }
      mwcpSZ[*offset] = parent -> splitValueFactSize;
      for (i = 1; i <= mwcpSZ[*offset]; i++) {
        (*mwcpPtr) ++;
        **mwcpPtr = (parent -> splitValueFactPtr)[i];
        if (getTraceFlag() & SUMM_HGH_TRACE) {
          Rprintf("%8x ",  (parent -> splitValueFactPtr)[i]);
        }
      }
      contPT[*offset] = NA_REAL;
    }
    else {
      contPT[*offset] = parent -> splitValueCont;
      mwcpSZ[*offset] = 0;
    }
  }  
  else {
    contPT[*offset] = NA_REAL;
    mwcpSZ[*offset] = 0;
  }
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\n     index     treeID     nodeID     parmID       contPT     mwcpSZ \n");
    Rprintf("%10d %10d %10d %10d %12.4f %10d \n", 
            *offset, treeID[*offset], 
            nodeID[*offset], parmID[*offset], 
            contPT[*offset], mwcpSZ[*offset]);
  }
  (*offset) ++;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    saveTree(b, parent ->  left, offset, treeID, nodeID, parmID, contPT, mwcpSZ, mwcpPtr);
    saveTree(b, parent -> right, offset, treeID, nodeID, parmID, contPT, mwcpSZ, mwcpPtr);
  }
  if (getTraceFlag() & SUMM_HGH_TRACE) {
    Rprintf("\nsaveTree() EXIT ...\n");
  }
}
void freeTree(Node *parent) {
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    freeTree(parent -> left);
    freeTree(parent -> right);
  }
  freeNode(parent, _eventTypeSize, _sortedTimeInterestSize);
}
