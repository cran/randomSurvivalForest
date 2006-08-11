//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 1.0.0
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

#ifndef NRUTIL_H
#define NRUTIL_H
#include "global.h"
float ran1(int *idum);
void hpsort(
  double *ra,
  uint n
);
double **dmatrix(
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
double **pdvector(
  ulong nl,
  ulong nh
);
double *dvector(
  ulong nl,
  ulong nh
);
uint **uimatrix(
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
uint *uivector(
  ulong nl,
  ulong nh
);
uchar *ucvector(
  ulong nl,
  ulong nh
);
void nrerror(char error_text[]);
void free_dmatrix(
  double **m,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_pdvector(
  double **v,
  ulong nl,
  ulong nh
);
void free_dvector(
  double *v,
  ulong nl,
  ulong nh
);
void free_uimatrix(
  uint **m,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_uivector(
  uint *v,
  ulong nl,
  ulong nh
);
void free_ucvector(
  uchar *v,
  ulong nl,
  ulong nh
);
void copyMatrix(
  uint **new,
  uint **old,
  uint nrow,
  uint ncol
);
#endif
