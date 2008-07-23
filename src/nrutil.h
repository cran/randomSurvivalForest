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

#ifndef NRUTIL_H
#define NRUTIL_H
float ran1(int *idum);
float ran2(int *idum);
uint upower (uint x, uint n);
uint upower2 (uint n);
uint ulog2 (uint n);
void hpsort(
  double *ra,
  uint n
);
void hpsortui(
  uint *ra,
  uint n
);
void indexx(
  uint n,
  double *arr, 
  uint *indx
);
void permute(
  uint n,
  uint *indx
);
double ***dmatrix3(
  ulong n3l,
  ulong n3h,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
double **dmatrix(
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
int **imatrix(
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
uint **uimatrix(
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
ulong **ulmatrix(
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
char **cmatrix(
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
double *dvector(
  ulong nl,
  ulong nh
);
double **pdvector(
  ulong nl,
  ulong nh
);
int *ivector(
  ulong nl,
  ulong nh
);
int **pivector(
  ulong nl,
  ulong nh
);
uint *uivector(
  ulong nl,
  ulong nh
);
uint **puivector(
  ulong nl, 
  ulong nh
);
uint ***ppuivector(
  ulong nl,
  ulong nh
);
ulong *ulvector(
  ulong nl,
  ulong nh
);
ulong **pulvector(
  ulong nl,
  ulong nh
);
char *cvector(
  ulong nl,
  ulong nh
);
char **pcvector(
  ulong nl,
  ulong nh
);
void nrerror(char error_text[]);
void free_dmatrix3(
  double ***m,
  ulong n3l,
  ulong n3h,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_dmatrix(
  double **m,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_imatrix(
  int **m,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_uimatrix(
  uint **m,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_ulmatrix(
  ulong **m,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_cmatrix(
  char **m,
  ulong nrl,
  ulong nrh,
  ulong ncl,
  ulong nch
);
void free_dvector(
  double *v,
  ulong nl,
  ulong nh
);
void free_pdvector(
  double **v,
  ulong nl,
  ulong nh
);
void free_ivector(
  int *v,
  ulong nl,
  ulong nh
);
void free_uivector(
  uint *v,
  ulong nl,
  ulong nh
);
void free_pivector(
  int **v,
  ulong nl,
  ulong nh
);
void free_puivector(
  uint **v,
  ulong nl,
  ulong nh
);
void free_ppuivector(
  uint ***v,
  ulong nl,
  ulong nh
);
void free_ulvector(
  ulong *v,
  ulong nl,
  ulong nh
);
void free_pulvector(
  ulong **v, 
  ulong nl, 
  ulong nh
);
void free_cvector(
  char *v,
  ulong nl,
  ulong nh
);
void free_pcvector(
  char **v,
  ulong nl,
  ulong nh
);
void nrCopyMatrix(
  uint **new,
  uint **old,
  uint nrow,
  uint ncol
);
void nrCopyVector(
  char *new, 
  char *old, 
  uint ncol
);
#endif
