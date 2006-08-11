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

#include "nrutil.h"
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float ran1(int *idum) {
  int j;
  int k;
  static int iy = 0;
  static int iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) {
      *idum = 1;
    }
    else {
      *idum = -(*idum);
    }
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum) / IQ;
      *idum = IA * (*idum - k * IQ) - IR * k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0) *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp = AM * iy) > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
void hpsort(double *ra, uint n) {
  uint i, ir, j, l;
  double rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
#define FREE_ARG char*
#define NR_END 1
double **dmatrix(ulong nrl, ulong nrh, ulong ncl, ulong nch) {
  ulong i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("Allocation Failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("Allocation Failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  return m;
}
double **pdvector(ulong nl, ulong nh) {
  double **v;
  v=(double **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double*)));
  if (!v) nrerror("Allocation Failure in pdvector()");
  return v-nl+NR_END;
}
uint **uimatrix(ulong nrl, ulong nrh, ulong ncl, ulong nch) {
  ulong i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  uint **m;
  m=(uint **) malloc((size_t)((nrow+NR_END)*sizeof(uint*)));
  if (!m) nrerror("Allocation Failure 1 in uimatrix()");
  m += NR_END;
  m -= nrl;
  m[nrl]=(uint *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(uint)));
  if (!m[nrl]) nrerror("Allocation Failure 2 in uimatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  return m;
}
double *dvector(ulong nl, ulong nh) {
  double *v;
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation Failure in dvector()");
  return v-nl+NR_END;
}
uint *uivector(ulong nl, ulong nh) {
  uint *v;
  v=(uint *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(uint)));
  if (!v) nrerror("Allocation Failure in uivector()");
  return v-nl+NR_END;
}
uchar *ucvector(ulong nl, ulong nh) {
  uchar *v;
  v=(uchar *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(uchar)));
  if (!v) nrerror("Allocation Failure in ucvector()");
  return v-nl+NR_END;
}
void nrerror(char error_text[]) {
  Rprintf("\nRSF:  *** ERROR *** ");
  Rprintf("\nRSF:  Numerical Recipes Run-Time Error:");
  Rprintf("\nRSF:  %s", error_text);
  Rprintf("\nRSF:  Please Contact Technical Support.");
  Rprintf("\nRSF   The application will now exit.\n");
  exit(FALSE);
}
void free_dmatrix(double **m,
                  ulong nrl,
                  ulong nrh,
                  ulong ncl,
                  ulong nch
                 ) {
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_pdvector(double **v, ulong nl, ulong nh) {
  free((FREE_ARG) (v+nl-NR_END));
}
void free_uimatrix(uint **m,
                   ulong nrl,
                   ulong nrh,
                   ulong ncl,
                   ulong nch
                  ) {
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_dvector(double *v, ulong nl, ulong nh) {
  free((FREE_ARG) (v+nl-NR_END));
}
void free_uivector(uint *v,ulong nl,ulong nh) {
  free((FREE_ARG) (v+nl-NR_END));
}
void free_ucvector(uchar *v,ulong nl,ulong nh) {
  free((FREE_ARG) (v+nl-NR_END));
}
#undef FREE_ARG
#undef NR_END
void copyMatrix(uint **new, uint **old, uint nrow, uint ncol) {
  uint i,j;
  for (i = 1; i <= nrow; i++) {  
    for (j = 1; j <= ncol; j++) {  
      new[i][j] = old[i][j];
    }
  }
}
