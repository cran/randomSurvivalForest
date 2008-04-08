//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.2.3
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

extern int      *_seed_;
extern double   *_performance_;
extern uint     *_leafCount_;
extern uint     *_varUsed_;
extern uint     *_treeID_;
extern uint     *_nodeID_;
extern uint     *_parmID_;
extern double   *_spltPT_;
extern uint      _opt;
extern uint      _splitRule;
extern uint      _imputeSize;
extern uint      _forestSize;
extern uint      _minimumDeathCount;
extern uint      _randomCovariateCount;
extern double   *_randomCovariateWeight;
extern uint      _observationSize;
extern uint      _xSize;
extern double   *_time;
extern double   *_status;
extern double   *_xData;
extern uint      _fobservationSize;
extern double   *_ftime;
extern double   *_fstatus;
extern double   *_fxData;
extern uint      _timeInterestSize;
extern double   *_timeInterest;
extern SEXP      _sexp_xType;
extern uint      _intrPredictorSize;
extern uint     *_intrPredictor;
extern uint     *_intrObservation;
extern char    **_xType;
extern double  **_observation;
extern double  **_fobservation;
extern double   *_masterTime;
extern uint     *_masterTimeIndex;
extern uint      _masterTimeSize;
extern char     *_importanceFlag;
extern char      _mTimeIndexFlag; 
extern uint     *_mRecordMap;
extern uint     *_fmRecordMap;
extern uint      _mRecordSize;
extern uint      _fmRecordSize;
extern uint     *_mRecordIndex;
extern uint     *_fmRecordIndex;
extern uint      _mvSize;
extern uint      _fmvSize;
extern int     **_mvSign;
extern int     **_fmvSign;
extern int      *_mvIndex;
extern int      *_fmvIndex;
extern int     **_mvForestSign;
extern int     **_fmvForestSign;
extern double   *_mStatus;
extern double   *_fmStatus;
extern double   *_mTime;
extern double   *_fmTime;
extern double  **_oobEnsemblePtr;
extern double  **_fullEnsemblePtr;
extern double   *_ensembleRun;
extern double  **_vimpEnsembleRun;  
extern uint     *_oobEnsembleDen;
extern uint     *_fullEnsembleDen;
extern double   *_sImputeStatusPtr;
extern double   *_sImputeTimePtr;
extern double  **_sImputePredictorPtr;
extern double   *_sOOBImputeStatusPtr;
extern double   *_sOOBImputeTimePtr;
extern double  **_sOOBImputePredictorPtr;
extern uint   **_varUsedPtr;
extern int      *_seed1Ptr;
extern int      *_seed2Ptr;
extern Node    **_nodeMembership;
extern uint     *_bootMembershipIndex;
extern char     *_bootMembershipFlag;
extern uint     *_oobSampleSize;
extern Node    **_fnodeMembership;
extern char     *_fbootMembershipFlag;
extern uint     *_foobSampleSize;
