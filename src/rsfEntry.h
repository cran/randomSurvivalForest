////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.0
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

#ifndef RSFENTRY_H
#define RSFENTRY_H
SEXP rsfGrow(SEXP traceFlag,
             SEXP opt,
             SEXP seedPtr,
             SEXP splitRule,
             SEXP splitRandomRule,  
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
             SEXP randomCovariateWeight,
             SEXP xType,
             SEXP reimputeSize);
SEXP rsfPredict(SEXP traceFlag,
                SEXP opt,
                SEXP seedPtr,  
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
                SEXP contPT,
                SEXP mwcpSZ,
                SEXP mwcpPT,
                SEXP seed,
                SEXP xType);
SEXP rsfInteraction(SEXP traceFlag,
                    SEXP opt,
                    SEXP seedPtr,  
                    SEXP forestSize, 
                    SEXP observationSize,
                    SEXP time,
                    SEXP status,
                    SEXP xSize,
                    SEXP xData,
                    SEXP timeInterestSize,
                    SEXP timeInterest,
                    SEXP treeID,
                    SEXP nodeID,
                    SEXP parmID,
                    SEXP contPT,
                    SEXP mwcpSZ,
                    SEXP mwcpPT,
                    SEXP seed,
                    SEXP xType,
                    SEXP intrPredictorSize,
                    SEXP intrPredictor,
                    SEXP fobservationSize,
                    SEXP intrObservation);
SEXP rsfImpute(SEXP traceFlag,
               SEXP seedPtr,  
               SEXP splitRule,  
               SEXP splitRandomRule,  
               SEXP randomCovariateCount,  
               SEXP forestSize,  
               SEXP minimumDeathCount,
               SEXP observationSize,
               SEXP time,
               SEXP status,
               SEXP xSize,
               SEXP xData,
               SEXP randomCovariateWeight,
               SEXP xType,
               SEXP imputeSize);
#endif
