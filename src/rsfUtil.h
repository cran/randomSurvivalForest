//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.2.1
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

#ifndef RSFUTIL_H
#define RSFUTIL_H
#include "extern.h"
uint updateTimeStamp(uint before);
void freeTree(Node *parent);
char restoreTree(uint    b,
                 Node   *parent,
                 uint   *leafCount,
                 uint   *offset,
                 uint   *treeID,
                 uint   *nodeID,
                 uint   *parmID,
                 double *spltPT);
void saveTree(uint    b,
              Node   *parent,
              uint   *offset,
              uint   *treeID,
              uint   *nodeID,
              uint   *parmID,
              double *spltPT);
char testNodeSize(Node *parent);
Node* getMembership(Node    *parent,
                    double **predictor,
                    uint     index);
Node *getTerminalNode(uint leaf);
Node* randomizeMembership(Node    *parent, 
                          double **predictor, 
                          uint     individual, 
                          uint     splitParameter);
double getConcordanceIndex(int     polarity,
                           uint    size, 
                           double *statusPtr, 
                           double *timePtr, 
                           double *predictedOutcome,
                           uint   *oobCount);
void getNelsonAalenEstimate(double **cumulativeHazard,
                            uint     treeID,
                            uint     sortedTimeInterestSize);
void updateEnsembleCHF(uint     mode, 
                       uint     sortedTimeInterestSize,
                       uint     treeID,
                       double **cumulativeHazard);
void getMeanSurvivalTime(double *meanSurvivalTime,
                         uint    treeID);
void updateEnsembleSurvivalTime(uint    mode, 
                                uint    treeID, 
                                double *meanSurvivalTime);
void getVariableImportance(char     mode,
                           uint     sortedTimeInterestSize,
                           uint     leafCount,
                           double **cumulativeHazard,
                           double  *meanSurvivalTime,
                           Node    *rootPtr,
                           uint     b);
void getVariablesUsed(Node *rootPtr, uint *varUsedPtr);
void updateEnsembleEvents (char multipleImputeFlag,
                           uint      mode,
                           uint      sortedTimeInterestSize,
                           Node     *rootPtr,
                           uint      b,
                           char    **dmRecordBootFlag,
                           double ***dmvImputation);
#endif
