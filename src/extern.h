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

#ifndef EXTERNAL_H
#define EXTERNAL_H
#include          "node.h"
#include        "factor.h"
extern int      *_seed_;
extern double   *_performance_;
extern uint     *_leafCount_;
extern uint     *_varUsed_;
extern uint     *_treeID_;
extern uint     *_nodeID_;
extern uint     *_parmID_;
extern uint     *_mwcpSZ_;
extern double   *_contPT_;
extern uint     *_mwcpPT_;
extern uint      _opt;
extern uint      _splitRule;
extern uint      _splitRandomRule;
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
extern double   *rsf_ftime;
extern double   *_fstatus;
extern double   *_fxData;
extern uint      _timeInterestSize;
extern double   *_timeInterest;
extern SEXP      _sexp_xType;
extern uint      _intrPredictorSize;
extern uint     *_intrPredictor;
extern uint     *_intrIndividual;
extern uint     *_individualIndex;
extern uint     *_predictorIndex;
extern char    **_xType;
extern double  **_observation;
extern double  **_fobservation;
extern uint     *_eventType;
extern uint      _eventTypeSize;
extern uint     *_eventTypeIndex;
extern uint     *_eventTypeIndexSize;
extern double   *_masterTime;
extern uint     *_masterTimeIndex;
extern uint      _masterTimeSize;
extern char     *_importanceFlag;
extern uint      _sortedTimeInterestSize;
extern uint      _factorCount;
extern uint     *_factorMap;
extern uint     *_factorIndex;
extern uint     *_factorSize;
extern uint      _maxFactorLevel;
extern Factor  **_factorList;
extern uint      _mFactorSize;
extern uint      _fmFactorSize;
extern uint     *_mFactorIndex;
extern uint     *_fmFactorIndex;
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
extern double **_performancePtr;
extern double  ***_oobEnsemblePtr;
extern double  ***_fullEnsemblePtr;
extern uint      *_oobEnsembleDen;
extern uint      *_fullEnsembleDen;
extern uint     **_oobVimpInvalidDen;
extern double   **_importancePtr;
extern double   **_vimpMortality;
extern double ****_crVimpEnsemble;
extern double  ***_crVimpPOE;
extern uint      _mStatusSize;
extern uint     *_eIndividualSize;
extern uint     *_meIndividualSize;
extern uint    **_eIndividual;
extern double  ***_oobSubSurvivalPtr;
extern double  ***_fullSubSurvivalPtr;
extern double  ***_oobSubDistributionPtr;
extern double  ***_fullSubDistributionPtr;
extern double   **_oobPOEPtr;
extern double   **_fullPOEPtr;
extern double   *_sImputeStatusPtr;
extern double   *_sImputeTimePtr;
extern double  **_sImputePredictorPtr;
extern double   *_sOOBImputeStatusPtr;
extern double   *_sOOBImputeTimePtr;
extern double  **_sOOBImputePredictorPtr;
extern uint   **_varUsedPtr;
extern double **_splitDepthPtr;
extern uint     _totalMWCPCount;
extern int      *_seed1Ptr;
extern int      *_seed2Ptr;
extern Node    **_nodeMembership;
extern uint     *_bootMembershipIndex;
extern char     *_bootMembershipFlag;
extern uint     *_oobSampleSize;
extern char     *_genericMembershipFlag;
extern Node    **_fnodeMembership;
extern uint     *_foobSampleSize;
extern double   _splitValueMaxCont;
extern uint     _splitValueMaxFactSize;
extern uint    *_splitValueMaxFactPtr;
extern clock_t _benchTime;
extern clock_t _splitTime;
extern clock_t _hazrdTime;
extern clock_t _chazrdTime;
extern clock_t _ensblTime;
extern clock_t _censblTime;
extern clock_t _censblTimeSub;
extern clock_t _vimprTime;
extern clock_t _cindxTime;
#endif
