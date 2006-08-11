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

#include "node.h"
#define LOG_RANK        1
#define CONSERVE_EVENTS 2
void freeTree(Node *parent);
Node* getMembership(
  Node *parent,
  uint  index
);
void getCumulativeHazardEstimate(
  double **cumulativeHazard,
  Node    *parent,
  double  *timeInterest,
  uint     sortedTimeInterestSize,
  Node   **nodeMembership,
  double  *masterTime
);
char logRank(
  Node  *parent,
  uint  *splitParameterMax,
  uint  *splitValueMax,
  Node **nodeMembership,
  uint   masterDeathTimeSize
);
char conserveEvents(
  Node  *parent,
  uint  *splitParameterMax,
  uint  *splitValueMax,
  Node **nodeMembership,
  uint   masterDeathTimeSize
);
char getBestSplit(
  Node  *parent,
  uint  *splitParameterMax,
  uint  *splitValueMax,
  Node **nodeMembership,
  uint   masterDeathTimeSize,
  uint  *splitRule
);
char makeTree(
  Node  *parent,
  Node **nodeMembership,
  uint  *leafCount, 
  uint   masterDeathTimeSize,
  uint  *splitRule
);
char forkAndUpdate(
  Node **nodeMembership,
  uint  *leafCount,
  Node  *parent,
  uint   splitParameter,
  uint   splitValue
);
void rsf(
  uint   *traceFlag,
  int    *memoryUseProtocol,
  int    *seed,
  uint   *splitRule,
  uint   *randomCovariateCount,
  uint   *bootstrapSampleCount,
  uint   *minimumDeathCount,
  uint   *observationSizePtr,
  double *time,
  uint   *status,
  uint   *timeInterestSize,
  double *timeInterest,
  uint   *xSizePtr,
  double *xData,
  double *ensembleEstimator,
  double *performanceMeasure,
  uint   *leafCount,
  uint   *proximity
);
