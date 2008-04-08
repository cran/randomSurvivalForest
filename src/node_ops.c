//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.2.3
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
#include "node_ops.h"
extern uint getTraceFlag();
Node *makeNode() {
  Node *parent = (Node*) malloc((size_t)sizeof(Node));
  parent -> permissibleSplit = cvector(1, _xSize);
  return parent;
}
#define NR_END 1
#define FREE_ARG char*
Node ***nodePtrMatrix(ulong nrl, ulong nrh, ulong ncl, ulong nch) {
  ulong i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  Node ***m;
  m = (Node ***) malloc((size_t)((nrow+1)*sizeof(Node**)));
  if (!m) nrerror("allocation failure 1 in nodePtrMatrix()");
  m += 1;
  m -= nrl;
  m[nrl] = (Node **) malloc((size_t)((nrow*ncol+1)*sizeof(Node*)));
  if (!m[nrl]) nrerror("allocation failure 2 in nodePtrMatrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  return m;
}
Node **nodePtrVector(ulong nl, ulong nh) {
  Node **v;
  v = (Node **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(Node*)));     
  if (!v) nrerror("allocation failure in nodePtrVector()");
  return v-nl+NR_END;
}
void free_nodePtrMatrix(Node ***m,
                        ulong nrl,
                        ulong nrh,
                        ulong ncl,
                        ulong nch) {
  free((FREE_ARG) (m[nrl]+ncl-1));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_nodePtrVector(Node **v,
                        ulong nl,
                        ulong nh) {
  free((FREE_ARG) (v+nl-NR_END));
}
void free_Node(Node *parent) {
  free_cvector(parent -> permissibleSplit, 1, _xSize);
  free((FREE_ARG) parent);
}
#undef NR_END
#undef FREE_ARG
void getNodeInfo(Node *leaf) {
  uint i;
  Rprintf("\nNodeInfo:  ");
  Rprintf("\n   LeafCnt   SpltParm      SpltVal");
  Rprintf("\n%10d %10d %12.4f", leaf -> leafCount, leaf -> splitParameter, leaf -> splitValue);
  Rprintf("\nPermissible Splits \n");
  for (i=1; i <= _xSize; i++) {
    Rprintf("%10d \n", (leaf -> permissibleSplit)[i]);
  }
}
void setParent(Node *daughter, Node *parent) {
  daughter -> parent = parent;
}
void setLeftDaughter(Node *daughter, Node *parent) {
  parent -> left = daughter;
}
void setRightDaughter(Node *daughter, Node *parent) {
  parent -> right = daughter;
}
char forkNode (Node  *parent, 
               uint   splitParameter, 
               double splitValue) {
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nforkNode() ENTRY ...\n");
  }
  if (parent == NULL) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Inconsistent call to forkNode().  ");
    Rprintf("\nRSF:  The parent node is NULL.");
    return FALSE;
  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Inconsistent call to forkNode().  ");
    Rprintf("\nRSF:  The daughter nodes are NON-NULL.");
    return FALSE;
  }
  if (parent -> splitFlag == FALSE) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Inconsistent call to forkNode().  ");
    Rprintf("\nRSF:  The split flag is FALSE.");
    return FALSE;
  }
  if (_xSize < splitParameter) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Inconsistent call to forkNode().  ");
    Rprintf("\nRSF:  The split parameter index is out of range [1, _xSize].");
    return FALSE;
  }
  if ((parent -> permissibleSplit)[splitParameter] == FALSE) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Inconsistent call to forkNode().  ");
    Rprintf("\nRSF:  The split parameter is marked unsplittable.");
    return FALSE;
  }
  Node *left  = makeNode();
  Node *right = makeNode();
  parent -> splitParameter = splitParameter;
  parent -> splitValue = splitValue;
  setParent(left, parent);
  setParent(right, parent);
  setLeftDaughter(left, parent);
  setRightDaughter(right, parent);
  nrCopyVector(left  -> permissibleSplit, parent -> permissibleSplit, _xSize);
  nrCopyVector(right -> permissibleSplit, parent -> permissibleSplit, _xSize);
  left  -> left  = left  -> right = NULL;
  right -> left  = right -> right = NULL;
  left -> splitParameter  = 0;
  left -> splitValue      = NA_REAL;
  left -> splitFlag       = TRUE;
  left -> leafCount       = 0;
  right -> splitParameter  = 0;
  right -> splitValue      = NA_REAL;
  right -> splitFlag       = TRUE;
  right -> leafCount       = 0;
  parent -> splitFlag = FALSE;
  if (getTraceFlag() & DL2_TRACE) {
    Rprintf("\nParent Info:  "); getNodeInfo(parent);
    Rprintf("\nLeft Info:    "); getNodeInfo(left);
    Rprintf("\nRight Info:   "); getNodeInfo(right);
  }
  if (getTraceFlag() & DL1_TRACE) {
    Rprintf("\nforkNode() EXIT ...\n");
  }
  return TRUE;
}
