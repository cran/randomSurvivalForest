////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.2
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
