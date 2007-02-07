//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 2.0.0
//
//  Copyright 2006, Cleveland Clinic
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
#include <math.h>
#include <float.h>
#include <time.h>
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
#define ACTIVE 0x02
#define DL0_TRACE 0x01
#define DL1_TRACE 0x02
#define DL2_TRACE 0x04
#define DL3_TRACE 0x08
#define EPSILON 1.0e-7
#define MUP_PROX 0x01
#define MUP_TREE 0x02
#define MUP_PERF 0x04
#define MUP_VIMP 0x08
#define RSF_OUTP_ID   0  
#define RSF_STRG_ID   1  
#define RSF_FENS_ID   2  
#define RSF_OENS_ID   3  
#define RSF_PERF_ID   4  
#define RSF_LEAF_ID   5  
#define RSF_PROX_ID   6  
#define RSF_VIMP_ID   7  
#define RSF_TREE_ID   8  
#define RSF_NODE_ID   9  
#define RSF_PARM_ID  10  
#define RSF_SPLT_PT  11  
#define RSF_SEED_ID  12  
#define RSF_SEXP_CNT 13  
#define RSF_GROW 0x10
#define RSF_PRED 0x20
#define LOG_RANK        1
#define CONSERVE_EVENTS 2
#define LOG_RANK_SCORE  3
#define LOG_RANK_APPROX 4
#define uchar unsigned char
#define uint  unsigned int
#define ulong unsigned long
