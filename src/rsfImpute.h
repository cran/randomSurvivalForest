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

#ifndef RSFIMPUTE_H
#define RSFIMPUTE_H
#include "node.h"
void imputeInteraction (uint treeID, Node *parent);
char imputeNode (uint     type,
                 char     seedChainFlag,
                 uint     treeID, 
                 Node    *parent);
void imputeTree(uint mode, uint b, Node *parent, char rootFlag);
void imputeUpdateShadow (uint      mode, 
                         char      selectionFlag,
                         double ***dmvImputation,
                         double   *shadowStatus, 
                         double   *shadowTime, 
                         double  **shadowPredictor);
void imputeUpdateSummary (uint     mode, 
                          double  *statusPtr, 
                          double  *timePtr, 
                          double **predictorPtr, 
                          uint     treeID,
                          double ***dmvImputationPtr);
void imputeSummary(uint      mode,
                   char      selectionFlag,
                   char    **dmRecordBootFlag,
                   double ***dmvImputation);
void imputeConcordance(uint      mode,
                       uint      b,
                       char    **dmRecordBootFlag,
                       double ***dmvImputation,
                       double   *tempStatus,
                       double   *tempTime);
void imputeCommon(uint      mode,
                  uint      b,
                  char      selectionFlag,
                  char    **dmRecordBootFlag,
                  double ***dmvImputation,
                  char      predictorFlag);
void unImpute (uint mode);
void imputeMultipleTime (char selectionFlag);
double getMaximalValue(double *value, uint size);
double getMedianValue(double *value, uint size);
double getMeanValue(double *value, uint size);
double getSampleValue(double *value, uint size, char chainFlag);
uint getRecordMap(uint   *map, 
                  uint    size, 
                  double *status, 
                  double *time, 
                  double *data);
char getForestSign (uint mode, uint b);
void updateTimeIndexArray(Node *parent);
void updateEventTypeSubsets(double *summaryStatus, 
                            uint    mRecordSize,
                            int   **mvSign,
                            uint   *mRecordIndex);
#endif
