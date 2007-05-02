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

#include "node.h"
char logRankApprox(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit
);
char logRankScore(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit,
  uint   **masterSplitOrder
);
char conserveEvents(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit
);
char logRank(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit
);
uint selectRandomCovariates(
  Node *parent,
  uint *covariateIndex);
char getBestSplit(
  Node    *parent,
  uint    *splitParameterMax,
  uint    *splitValueMax,
  uint     masterDeathTimeSize,
  uint     masterTimeSize,
  double **masterSplit,
  uint   **masterSplitOrder
);
char makeTree(
  Node    *parent,
  uint    *leafCount, 
  uint     masterDeathTimeSize,
  uint     masterTimesize,
  double **masterSplit,
  uint   **masterSplitOrder
);
char forkAndUpdate(
  uint  *leafCount,
  Node  *parent,
  uint   splitParameter,
  uint   splitValueIndex,
  double splitValue
);
char bootstrap(
  uint     mode,
  uint     b,
  Node    *rootPtr,
  uint    *oobSampleSize,
  double **masterSplit,
  uint    *masterSplitSize,
  uint   **masterSplitOrder
);
SEXP rsfGrow(
  SEXP traceFlag,
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
);
SEXP rsfPredict(
  SEXP traceFlag,
  SEXP mup,
  SEXP forestsize,
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
);
SEXP rsf(uint mode, uint traceFlag);
