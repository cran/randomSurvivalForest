//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.5.1
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

#ifndef RSFSPLITUTIL_H
#define RSFSPLITUTIL_H
#include "extern.h"
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
void stackSplitCompact(uint size,
                       uint **nodeParentDeath,
                       uint **nodeParentAtRisk,
                       uint **nodeLeftDeath,
                       uint **nodeLeftAtRisk,
                       uint **nodeRightDeath,
                       uint **nodeRightAtRisk);
void unstackSplitCompact(uint size,
                         uint *nodeParentDeath,
                         uint *nodeParentAtRisk,
                         uint *nodeLeftDeath,
                         uint *nodeLeftAtRisk,
                         uint *nodeRightDeath,
                         uint *nodeRightAtRisk);
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
                         uint *rightDeathTimeSize);
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
