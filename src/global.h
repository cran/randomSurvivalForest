//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.5.0
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

#include <R_ext/Print.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "node.h"
#include "factor.h"
#ifndef NULL
#define NULL 0
#endif
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
#define ACTIVE   0x02
#define INACTIVE 0xFF
#define LEFT      0x01
#define RIGHT     0x00
#define SUMM_USR_TRACE  0x0001
#define SUMM_LOW_TRACE  0x0002
#define SUMM_MED_TRACE  0x0004
#define SUMM_HGH_TRACE  0x0008
#define SPLT_LOW_TRACE  0x0010
#define SPLT_MED_TRACE  0x0020
#define SPLT_HGH_TRACE  0x0040
#define FORK_DEF_TRACE  0x0080
#define MISS_LOW_TRACE  0x0100
#define MISS_MED_TRACE  0x0200
#define MISS_HGH_TRACE  0x0400
#define OUTP_DEF_TRACE  0x0800
#define NUMR_DEF_TRACE  0x1000
#define FACT_LOW_TRACE  0x2000
#define FACT_HGH_TRACE  0x4000
#define TIME_DEF_TRACE  0x8000
#define TURN_OFF_TRACE  0x0000
#define TURN_ON_TRACE   0x0001
#define EPSILON 1.0e-7
#define OPT_FENS      0x0001
#define OPT_OENS      0x0002
#define OPT_PERF      0x0004
#define OPT_PROX      0x0008
#define OPT_LEAF      0x0010
#define OPT_TREE      0x0020
#define OPT_SEED      0x0040
#define OPT_MISS      0x0080
#define OPT_OMIS      0x0100
#define OPT_VIMP_TYPE 0x0200
#define OPT_VIMP_JOIN 0x0400
#define OPT_VIMP      0x0800
#define OPT_VUSE_TYPE 0x1000
#define OPT_VUSE      0x2000
#define OPT_POUT_TYPE 0x4000
#define CENS_IDX -1
#define TIME_IDX -2
#define RSF_OUTP_ID   0  
#define RSF_STRG_ID   1  
#define RSF_FENS_ID   2  
#define RSF_OENS_ID   3  
#define RSF_PERF_ID   4  
#define RSF_PROX_ID   5  
#define RSF_LEAF_ID   6  
#define RSF_TREE_ID   7  
#define RSF_NODE_ID   8  
#define RSF_PARM_ID   9  
#define RSF_CONT_PT  10  
#define RSF_MWCP_SZ  11  
#define RSF_MWCP_PT  12  
#define RSF_SEED_ID  13  
#define RSF_VIMP_ID  14  
#define RSF_MISS_ID  15  
#define RSF_OMIS_ID  16  
#define RSF_VUSE_ID  17  
#define RSF_SEXP_CNT 18  
#define RSF_GROW 0x01
#define RSF_PRED 0x02
#define RSF_INTR 0x04
#define LOG_RANK        1
#define CONSERVE_EVENTS 2
#define LOG_RANK_SCORE  3
#define RANDOM_SPLIT    4
#define APROX 0
#define EXACT 1
#define SIZE_OF_INTEGER sizeof(uint)
#define MAX_EXACT_LEVEL SIZE_OF_INTEGER * 8
#define SAFE_FACTOR_SIZE 16
#define uint  unsigned int
#define ulong unsigned long
