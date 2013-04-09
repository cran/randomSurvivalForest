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

#ifndef RSFIMPORTANCE_H
#define RSFIMPORTANCE_H
#include "node.h"
Node *getProxyMember(Node    *parent,
                     double **predictor,
                     uint     index);
Node *randomizeMembership(Node    *parent, 
                          double **predictor, 
                          uint     individual, 
                          uint     splitParameter);
void permute(uint n, uint *indx);
void getVariableImportance (uint      mode,
                            uint      leafCount,
                            Node     *rootPtr,
                            uint      b);
void getVimpRandom (uint      mode,
                    Node     *rootPtr,
                    double  **predictorPtr,
                    uint      obsSize,
                    uint      varSize,
                    char      selectionFlag);
void getVimpPermute(uint      mode,
                    Node     *rootPtr,
                    double  **predictorPtr,
                    uint      b,
                    uint      obsSize,
                    uint      varSize,
                    char      selectionFlag);
void finalizeVariableImportance(uint       mode,
                                uint       rejectedTreeCount, 
                                char     **dmRecordBootFlag,
                                double  ***dmvImputation);
#endif
