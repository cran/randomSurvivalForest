//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.0.0
//
//  Copyright 2007, Cleveland Clinic
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

#ifndef RSFIMPUTE_H
#define RSFIMPUTE_H
#include "extern.h"
char imputeNode(uint type, uint b, Node *parent);
void imputeTree(char mode, uint b, Node *parent);
void imputeUpdate(char mode, double **mvImputationPtr);
void imputeSummary(char      mode,
                   char    **mRecordBootFlag,
                   double ***mvImputation,
                   double   *sumStatus,
                   double   *sumTime,
                   double  **sumPredictor);
void imputeConcordance(char      mode,
                       uint      b,
                       char    **mRecordBootFlag,
                       double ***mvImputation,
                       double   *sumStatus,
                       double   *sumTime,
                       double   *status,
                       double   *time);
void imputeCommon(char      mode,
                  char      recordFlag,
                  uint      b,
                  char    **mRecordBootFlag,
                  double ***mvImputation,
                  double   *sumStatus,
                  double   *sumTime,
                  double  **sumPredictor);
void unImpute (char);
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
#endif
