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

#ifndef NRUTIL_H
#define NRUTIL_H
float ran1(int *idum);
float ran2(int *idum);
unsigned int upower (unsigned int x, unsigned int n);
unsigned int upower2 (unsigned int n);
unsigned int ulog2 (unsigned int n);
void hpsort(
  double *ra,
  unsigned int n
);
void hpsortui(
  unsigned int *ra,
  unsigned int n
);
void indexx(
  unsigned int n,
  double *arr, 
  unsigned int *indx
);
double ****dmatrix4(unsigned long n4l,
                    unsigned long n4h,
                    unsigned long n3l, 
                    unsigned long n3h, 
                    unsigned long nrl, 
                    unsigned long nrh, 
                    unsigned long ncl, 
                    unsigned long nch);
double ***dmatrix3(
  unsigned long n3l,
  unsigned long n3h,
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
double **dmatrix(
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
int **imatrix(
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
unsigned int **uimatrix(
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
unsigned long **ulmatrix(
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
char **cmatrix(
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
double *dvector(
  unsigned long nl,
  unsigned long nh
);
double **pdvector(
  unsigned long nl,
  unsigned long nh
);
double ***ppdvector(
  unsigned long nl,
  unsigned long nh
);
int *ivector(
  unsigned long nl,
  unsigned long nh
);
int **pivector(
  unsigned long nl,
  unsigned long nh
);
unsigned int *uivector(
  unsigned long nl,
  unsigned long nh
);
unsigned int **puivector(
  unsigned long nl, 
  unsigned long nh
);
unsigned int ***ppuivector(
  unsigned long nl,
  unsigned long nh
);
unsigned long *ulvector(
  unsigned long nl,
  unsigned long nh
);
unsigned long **pulvector(
  unsigned long nl,
  unsigned long nh
);
char *cvector(
  unsigned long nl,
  unsigned long nh
);
char **pcvector(
  unsigned long nl,
  unsigned long nh
);
void nrerror(char error_text[]);
void free_dmatrix4(double    ****m,
                   unsigned long n4l,
                   unsigned long n4h,
                   unsigned long n3l,
                   unsigned long n3h,
                   unsigned long nrl,
                   unsigned long nrh,
                   unsigned long ncl,
                   unsigned long nch
                   );
void free_dmatrix3(
  double ***m,
  unsigned long n3l,
  unsigned long n3h,
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
void free_dmatrix(
  double **m,
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
void free_imatrix(
  int **m,
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
void free_uimatrix(
  unsigned int **m,
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
void free_ulmatrix(
  unsigned long **m,
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
void free_cmatrix(
  char **m,
  unsigned long nrl,
  unsigned long nrh,
  unsigned long ncl,
  unsigned long nch
);
void free_dvector(
  double *v,
  unsigned long nl,
  unsigned long nh
);
void free_pdvector(
  double **v,
  unsigned long nl,
  unsigned long nh
);
void free_ppdvector(
  double ***v,
  unsigned long nl,
  unsigned long nh
);
void free_ivector(
  int *v,
  unsigned long nl,
  unsigned long nh
);
void free_uivector(
  unsigned int *v,
  unsigned long nl,
  unsigned long nh
);
void free_pivector(
  int **v,
  unsigned long nl,
  unsigned long nh
);
void free_puivector(
  unsigned int **v,
  unsigned long nl,
  unsigned long nh
);
void free_ppuivector(
  unsigned int ***v,
  unsigned long nl,
  unsigned long nh
);
void free_ulvector(
  unsigned long *v,
  unsigned long nl,
  unsigned long nh
);
void free_pulvector(
  unsigned long **v, 
  unsigned long nl, 
  unsigned long nh
);
void free_cvector(
  char *v,
  unsigned long nl,
  unsigned long nh
);
void free_pcvector(
  char **v,
  unsigned long nl,
  unsigned long nh
);
void nrCopyMatrix(
  unsigned int **new,
  unsigned int **old,
  unsigned int nrow,
  unsigned int ncol
);
void nrCopyVector(
  char *new, 
  char *old, 
  unsigned int ncol
);
#endif
