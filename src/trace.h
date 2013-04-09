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

#ifndef TRACE_H
#define TRACE_H
#define SUMM_USR_TRACE  0x000001
#define SUMM_LOW_TRACE  0x000002
#define SUMM_MED_TRACE  0x000004
#define SUMM_HGH_TRACE  0x000008
#define SPLT_DEF_TRACE  0x000010
#define SPLT_LOW_TRACE  0x000020
#define SPLT_MED_TRACE  0x000040
#define SPLT_HGH_TRACE  0x000080
#define FORK_DEF_TRACE  0x000100  
#define MISS_LOW_TRACE  0x000200
#define MISS_MED_TRACE  0x000400
#define MISS_HGH_TRACE  0x000800
#define OUTP_DEF_TRACE  0x001000
#define NUMR_DEF_TRACE  0x002000  
#define FACT_LOW_TRACE  0x004000
#define FACT_HGH_TRACE  0x008000
#define ENSB_LOW_TRACE  0x010000
#define ENSB_HGH_TRACE  0x020000
#define TIME_DEF_TRACE  0x040000
#define TURN_OFF_TRACE  0x00000000
#define TURN_ON_TRACE   0x00000001
#define TRACE_MASK      0x00FFFFFF
#define TRACE_ITER              24
void setTraceFlag(unsigned int traceFlag);
unsigned int getTraceFlag();
void updateTraceFlag(char reset);
unsigned int updateTimeStamp(unsigned int before);
#endif
