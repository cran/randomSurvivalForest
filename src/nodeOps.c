////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.4
////
////  Copyright 2013, Cleveland Clinic Foundation
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
////  Written by:
////    Hemant Ishwaran, Ph.D.
////    Director of Statistical Methodology
////    Professor, Division of Biostatistics
////    Clinical Research Building, Room 1058
////    1120 NW 14th Street
////    University of Miami, Miami FL 33136
////
////    email:  hemant.ishwaran@gmail.com
////    URL:    http://web.ccs.miami.edu/~hishwaran
////    --------------------------------------------------------------
////    Udaya B. Kogalur, Ph.D.
////    Adjunct Staff
////    Dept of Quantitative Health Sciences
////    Cleveland Clinic Foundation
////    
////    Kogalur & Company, Inc.
////    5425 Nestleway Drive, Suite L1
////    Clemmons, NC 27012
////
////    email:  commerce@kogalur.com
////    URL:    http://www.kogalur.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************

#include        <stdlib.h>
#include        "nrutil.h"
#include       "nodeOps.h"
#include   <R_ext/Print.h>
#define FORK_DEF_TRACE  0x000100  
extern unsigned int getTraceFlag();
#include <R_ext/Arith.h>
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
#define NR_END 1
#define FREE_ARG char*
Node ***nodePtrMatrix(unsigned long nrl, 
                      unsigned long nrh, 
                      unsigned long ncl, 
                      unsigned long nch) {
  unsigned long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
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
void free_nodePtrMatrix(Node ***m,
                        unsigned long nrl,
                        unsigned long nrh,
                        unsigned long ncl,
                        unsigned long nch) {
  free((FREE_ARG) (m[nrl]+ncl-1));
  free((FREE_ARG) (m+nrl-NR_END));
}
Node **nodePtrVector(unsigned long nl, unsigned long nh) {
  Node **v;
  v = (Node **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(Node*)));     
  if (!v) nrerror("allocation failure in nodePtrVector()");
  return v-nl+NR_END;
}
void free_nodePtrVector(Node **v,
                        unsigned long nl,
                        unsigned long nh) {
  free((FREE_ARG) (v+nl-NR_END));
}
Node *makeNode(unsigned int xSize) {
  unsigned int i;
  Node *parent = (Node*) malloc((size_t)sizeof(Node));
  parent -> xSize = xSize;
  parent -> permissibleSplit = cvector(1, xSize);
  for (i = 1; i <= xSize; i++) {
    (parent -> permissibleSplit)[i] = TRUE;
  }
  parent -> left               = NULL;  
  parent -> right              = NULL;  
  parent -> splitFlag            = TRUE;  
  parent -> mortality            = NA_REAL;  
  parent -> splitParameter       = 0;
  parent -> splitValueCont       = 0.0;  
  parent -> splitValueFactSize   = 0;
  parent -> splitValueFactPtr    = NULL;
  parent -> leafCount            = 0;
  parent -> depth                = 0;
  parent -> splitDepth           = NULL;
  parent -> poe                  = NULL;
  parent -> subSurvival          = NULL;
  parent -> eventCount           = 0;
  return parent;
}
void freeNode(Node        *parent, 
              unsigned int eventTypeSize, 
              unsigned int sortedTimeInterestSize) {
  free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
  if ((parent -> splitValueFactSize) > 0) {
    free_uivector(parent -> splitValueFactPtr, 1, parent -> splitValueFactSize);
  }
  if (((parent -> depth) > 0) && ((parent -> splitParameter) == 0)) {
    free_uivector(parent -> splitDepth, 1, parent -> depth);
  }
  if ((parent -> splitParameter) == 0) {
    if (parent -> poe != NULL) {
      free_uivector(parent -> poe, 1, eventTypeSize);
    }
    if (parent -> subSurvival != NULL) {
      free_dmatrix(parent -> subSurvival, 1, eventTypeSize, 1, sortedTimeInterestSize);
    }
  }
  free((FREE_ARG) parent);
}
#undef NR_END
#undef FREE_ARG
void getNodeInfo(Node *leaf) {
  unsigned int i;
  Rprintf("\nNodeInfo:  ");
  Rprintf("\n   LeafCnt   SpltParm  ");
  Rprintf("\n%10d %10d \n", leaf -> leafCount, leaf -> splitParameter);
  if (leaf -> splitValueFactSize > 0) {
    Rprintf("0x ");
    for (i = leaf -> splitValueFactSize; i >= 1; i--) {
      Rprintf("%8x ", (leaf -> splitValueFactPtr)[i]);
    }
  }
  else {
    Rprintf(" %12.4f \n", leaf -> splitValueCont);
  }
  Rprintf("\nPermissible Splits \n");
  for (i=1; i <= leaf -> xSize; i++) {
    Rprintf(" %10d", i);
  }
  Rprintf("\n");
  for (i=1; i <= leaf -> xSize; i++) {
    Rprintf(" %10d", (leaf -> permissibleSplit)[i]);
  }
  Rprintf("\n");
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
char forkNode(Node         *parent,
              unsigned int  splitParameter,
              unsigned int  splitValueMaxFactSize,
              unsigned int *splitValueMaxFactPtr,
              double        splitValueMaxCont) {
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
  if (parent -> xSize < splitParameter) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Inconsistent call to forkNode().  ");
    Rprintf("\nRSF:  The split parameter index is out of range [1, xSize].");
    return FALSE;
  }
  if ((parent -> permissibleSplit)[splitParameter] == FALSE) {
    Rprintf("\nRSF:  *** WARNING *** ");
    Rprintf("\nRSF:  Inconsistent call to forkNode().  ");
    Rprintf("\nRSF:  The split parameter is marked unsplittable.");
    return FALSE;
  }
  Node *left  = makeNode(parent -> xSize);
  Node *right = makeNode(parent -> xSize);
  parent -> splitParameter = splitParameter;
  parent -> splitValueCont = splitValueMaxCont;
  parent -> splitValueFactSize = splitValueMaxFactSize;
  parent -> splitValueFactPtr = splitValueMaxFactPtr;
  setParent(left, parent);
  setParent(right, parent);
  setLeftDaughter(left, parent);
  setRightDaughter(right, parent);
  nrCopyVector(left  -> permissibleSplit, parent -> permissibleSplit, left -> xSize);
  nrCopyVector(right -> permissibleSplit, parent -> permissibleSplit, right -> xSize);
  parent -> splitFlag = FALSE;
  return TRUE;
}
