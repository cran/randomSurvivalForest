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

#ifndef RSFSPLITUTIL_H
#define RSFSPLITUTIL_H
#include "node.h"
void updateMaximumSplit(double  delta, 
                        uint    randomCovariate,
                        uint    jLong,
                        char    factorFlag,
                        uint    mwcpSizeAbsolute,
                        double *deltaMax,
                        uint   *splitParameterMax,
                        void   *permissibleSplitPtr);
uint stackAndSelectRandomCovariates(Node     *parent,
                                    uint      nodeSize,
                                    uint     *nodeIndex,
                                    uint    **covariateIndex,
                                    double ***permissibleSplit,
                                    uint    **permissibleSplitSize);
void unstackRandomCovariates(uint     nodeSize, 
                             uint    *covariateIndex,
                             double **permissibleSplit,
                             uint    *permissibleSplitSize);
uint getSelectableElement(uint    length,
                          char   *permissible,
                          double *weight);
void stackSplit(uint **localMembershipIndex, 
                uint **localDeathTimeCount, 
                uint **localDeathTimeIndex);
void unstackSplit(uint *localMembershipIndex, 
                  uint *localDeathTimeCount, 
                  uint *localDeathTimeIndex);
char getDeathCount(Node *parent, 
                   uint *localMembershipIndex, 
                   uint *localDeathTimeCount, 
                   uint *localDeathTimeIndex,
                   uint *localMembershipSize,
                   uint *localDeathTimeSize);
void stackSplitCompact(uint   deathTimeSize,
                       uint **nodeParentDeath,
                       uint **nodeParentAtRisk,
                       uint **nodeLeftDeath,
                       uint **nodeLeftAtRisk,
                       uint **nodeRightDeath,
                       uint **nodeRightAtRisk,
                       uint   nodeSize,
                       char **localSplitIndicator);
void unstackSplitCompact(uint  deathTimeSize,
                         uint *nodeParentDeath,
                         uint *nodeParentAtRisk,
                         uint *nodeLeftDeath,
                         uint *nodeLeftAtRisk,
                         uint *nodeRightDeath,
                         uint *nodeRightAtRisk,
                         uint  nodeSize,
                         char *localSplitIndicator);
void getAtRisk(uint *localMembershipIndex,
               uint *localDeathTimeCount,
               uint *localDeathTimeIndex,
               uint  localMembershipSize,
               uint  localDeathTimeSize,
               uint *nodeParentDeath,
               uint *nodeParentAtRisk);
uint stackAndConstructSplitVector (uint     localMembershipSize,
                                   uint     randomCovariateIndex,
                                   double  *permissibleSplit,
                                   uint     permissibleSplitSize,
                                   char    *factorFlag,
                                   char    *deterministicSplitFlag,
                                   uint    *mwcpSizeAbsolute,
                                   void   **permissibleSplitPtr);
void unstackSplitVector(uint   permissibleSplitSize,
                        uint   splitLength,
                        char   factorFlag,
                        char   deterministicSplitFlag,
                        void  *permissibleSplitPtr);
void virtuallySplitNode (uint  localMembershipSize,
                         char  factorFlag,
                         uint  mwcpSizeAbsolute,
                         uint  randomCovariate,
                         uint *localMembershipIndex,
                         void *permissibleSplitPtr,
                         uint  offset,
                         uint  localDeathTimeSize,
                         uint *localDeathTimeIndex,
                         uint *nodeParentAtRisk,
                         uint *nodeParentDeath,
                         uint *nodeLeftAtRisk,
                         uint *nodeLeftDeath,
                         uint *leftDeathTimeSize,
                         uint *nodeRightAtRisk,
                         uint *nodeRightDeath,
                         uint *rightDeathTimeSize,
                         char *localSplitIndicator);
void getReweightedRandomPair(uint relativefactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result);
void getRandomPair(uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result);
void createRandomBinaryPair(uint    relativeFactorSize, 
                            uint    absoluteFactorSize, 
                            uint    groupSize, 
                            double *absolutelevel, 
                            uint   *pair);
void convertRelToAbsBinaryPair(uint    relativeFactorSize, 
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel, 
                               uint   *pair);
char summarizeSplitResult(uint splitParameterMax, double deltaMax);
#endif
