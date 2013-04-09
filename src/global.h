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

#ifndef GLOBAL_H
#define GLOBAL_H
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
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
#define EPSILON 1.0e-7
#define OPT_FENS      0x000001
#define OPT_OENS      0x000002
#define OPT_PERF      0x000004
#define OPT_PROX      0x000008
#define OPT_LEAF      0x000010
#define OPT_TREE      0x000020
#define OPT_SEED      0x000040
#define OPT_MISS      0x000080
#define OPT_OMIS      0x000100
#define OPT_VIMP_TYPE 0x000200
#define OPT_VIMP_JOIN 0x000400
#define OPT_VIMP      0x000800
#define OPT_VUSE_TYPE 0x001000
#define OPT_VUSE      0x002000
#define OPT_POUT_TYPE 0x004000
#define OPT_SPLT_DPTH 0x008000
#define OPT_IMPU_ONLY 0x010000
#define OPT_VOUT_TYPE 0x020000
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
#define RSF_DPTH_ID  18  
#define RSF_FPOE_ID  19  
#define RSF_OPOE_ID  20  
#define RSF_SEXP_CNT 21  
#define RSF_GROW   0x01
#define RSF_PRED   0x02
#define RSF_INTR   0x04
#define LOG_RANK        1
#define CONSERVE_EVENTS 2
#define LOG_RANK_SCORE  3
#define RANDOM_SPLIT    4
#define SPLIT_RULE_CR   4
#define LOG_RANK_LAU_CR 5
#define LOG_RANK_CR     6
#define SUB_CUM_HAZ_CR  7
#define SPLIT_RULE_CNT  7
#define APROX 0
#define EXACT 1
#define SIZE_OF_INTEGER sizeof(uint)
#define MAX_EXACT_LEVEL SIZE_OF_INTEGER * 8
#define SAFE_FACTOR_SIZE 16
typedef unsigned int   uint;
typedef unsigned long  ulong; 
typedef unsigned char  uchar; 
#endif
