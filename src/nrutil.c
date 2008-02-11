//**********************************************************************
//**********************************************************************
//
//  RANDOM SURVIVAL FOREST 3.2.0
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

#include   "global.h"
#include   "nrutil.h"
extern int  *_seed2Ptr;
extern uint getTraceFlag();
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
float ran2(int *idum) {
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
void hpsortui(uint *ra, uint n) {
  uint i, ir, j, l;
  uint rra;
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
#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
void indexx(uint n, double *arr, uint *indx) {
  uint i, j, k, l=1;
  uint indxt, itemp, ir=n;
  uint *istack, jstack=0;
  double a;
  istack = uivector(1, NSTACK);
  for (j=1; j<=n; j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1; j<=ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i=j-1; i>=l; i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1] = indx[i];
        }
        indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l], indx[ir])
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1], indx[ir])
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        SWAP(indx[l], indx[l+1])
      }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i], indx[j])
      }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free_uivector(istack, 1, NSTACK);
}
#undef SWAP
#undef M
#undef NSTACK
void permute(uint n, uint *indx) {
  uint i,j,k;
  for (i=1; i<= n; i++) {
    indx[i] = 0;
  }
  for (i=n; i > 0; i--) {
    k = ceil(ran2(_seed2Ptr)*(i*1.0));
    for (j = 1; k > 0; j++) {
      if (indx[j] == 0) {
        k--;
      }
    }
    indx[j-1] = i;
  }
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nRequest for Permutation:  \n");
    Rprintf("     index   permute\n");
    for (i=1; i <= n; i++) {
      Rprintf("%10d %10d \n", i, indx[i]);
    }
  }
}
#define FREE_ARG char*
#define NR_END 1
double ***dmatrix3(ulong n3l, ulong n3h, ulong nrl, ulong nrh, ulong ncl, ulong nch) {
  ulong i,j,k,  n3 = n3h - n3l + 1, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double ***m;
  m=(double ***) malloc((size_t)((n3+NR_END)*sizeof(double**)));
  if (!m) nrerror("Allocation Failure 1 in dmatrix3()");
  m += NR_END;
  m -= nrl;
  m[n3l]=(double **) malloc((size_t)((n3*nrow+NR_END)*sizeof(double*)));
  if (!m[n3l]) nrerror("Allocation Failure 2 in dmatrix3()");
  m[n3l] += NR_END;
  m[n3l] -= nrl;
  m[n3l][nrl]=(double *) malloc((size_t)((n3*nrow*ncol+NR_END)*sizeof(double)));
  if (!m[n3l][nrl]) nrerror("Allocation Failure 3 in dmatrix3()");
  m[n3l][nrl] += NR_END;
  m[n3l][nrl] -= ncl;
  for(i=n3l+1;i<=n3h;i++) {
    m[i] = m[i-1] + nrow;
  }
  k = 0;
  for(i=n3l;i<=n3h;i++) {
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = m[n3l][nrl] + ncol * k;
      k++;
    }
  }
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ndmatrix3() alloc:  %20x  \n", m);
  }
  return m;
}
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
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ndmatrix() alloc:  %20x  \n", m);
  }
  return m;
}
int **imatrix(ulong nrl, ulong nrh, ulong ncl, ulong nch) {
  ulong i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("Allocation Failure 1 in imatrix()");
  m += NR_END;
  m -= nrl;
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) nrerror("Allocation Failure 2 in imatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nimatrix() alloc:  %20x  \n", m);
  }
  return m;
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
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nuimatrix() alloc:  %20x  \n", m);
  }
  return m;
}
char **cmatrix(ulong nrl, ulong nrh, ulong ncl, ulong nch) {
  ulong i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  char **m;
  m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char*)));
  if (!m) nrerror("Allocation Failure 1 in cmatrix()");
  m += NR_END;
  m -= nrl;
  m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
  if (!m[nrl]) nrerror("Allocation Failure 2 in cmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ncmatrix() alloc:  %20x  \n", m);
  }
  return m;
}
double *dvector(ulong nl, ulong nh) {
  double *v;
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation Failure in dvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ndvector() alloc:  %20x  \n", v);
  }
  return v;
}
double **pdvector(ulong nl, ulong nh) {
  double **v;
  v=(double **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double*)));
  if (!v) nrerror("Allocation Failure in pdvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npdvector() alloc:  %20x  \n", v);
  }
  return v;
}
int *ivector(ulong nl, ulong nh) {
  int *v;
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("Allocation Failure in uivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nivector() alloc:  %20x  \n", v);
  }
  return v;
}
int **pivector(ulong nl, ulong nh) {
  int **v;
  v=(int **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int*)));
  if (!v) nrerror("Allocation Failure in puivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npivector() alloc:  %20x  \n", v);
  }
  return v;
}
uint *uivector(ulong nl, ulong nh) {
  uint *v;
  v=(uint *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(uint)));
  if (!v) nrerror("Allocation Failure in uivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nuivector() alloc:  %20x  \n", v);
  }
  return v;
}
uint **puivector(ulong nl, ulong nh) {
  uint **v;
  v=(uint **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(uint*)));
  if (!v) nrerror("Allocation Failure in puivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npdvector() alloc:  %20x  \n", v);
  }
  return v;
}
char *cvector(ulong nl, ulong nh) {
  char *v;
  v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char)));
  if (!v) nrerror("Allocation Failure in cvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ncvector() alloc:  %20x  \n", v);
  }
  return v;
}
char **pcvector(ulong nl, ulong nh) {
  char **v;
  v=(char **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char*)));
  if (!v) nrerror("Allocation Failure in pcvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npcvector() alloc:  %20x  \n", v);
  }
  return v;
}
void nrerror(char error_text[]) {
  Rprintf("\nRSF:  *** ERROR *** ");
  Rprintf("\nRSF:  Numerical Recipes Run-Time Error:");
  Rprintf("\nRSF:  %s", error_text);
  Rprintf("\nRSF:  Please Contact Technical Support.");
  Rprintf("\nRSF:  The application will now exit.\n");
  exit(FALSE);
}
void free_dmatrix3(double ***m,
                  ulong n3l,
                  ulong n3h,
                  ulong nrl,
                  ulong nrh,
                  ulong ncl,
                  ulong nch
                 ) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ndmatrix3() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[n3l][nrl]+ncl-NR_END));
  free((FREE_ARG) (m[n3l]+nrl-NR_END));
  free((FREE_ARG) (m+n3l-NR_END));
}
void free_dmatrix(double **m,
                  ulong nrl,
                  ulong nrh,
                  ulong ncl,
                  ulong nch
                 ) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ndmatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_imatrix(int **m,
                  ulong nrl,
                  ulong nrh,
                  ulong ncl,
                  ulong nch
                  ) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nimatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_uimatrix(uint **m,
                   ulong nrl,
                   ulong nrh,
                   ulong ncl,
                   ulong nch
                  ) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nuimatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_cmatrix(char **m,
                  ulong nrl,
                  ulong nrh,
                  ulong ncl,
                  ulong nch
                  ) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ncmatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_dvector(double *v, ulong nl, ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ndvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_pdvector(double **v, ulong nl, ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npdvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_ivector(int *v,ulong nl,ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_pivector(int **v, ulong nl, ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_uivector(uint *v,ulong nl,ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\nuivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_puivector(uint **v, ulong nl, ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npuivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_cvector(char *v,ulong nl,ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\ncvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_pcvector(char **v, ulong nl, ulong nh) {
  if (getTraceFlag() & DL3_TRACE) {
    Rprintf("\npcvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
#undef FREE_ARG
#undef NR_END
void nrCopyMatrix(uint **new, uint **old, uint nrow, uint ncol) {
  uint i,j;
  for (i = 1; i <= nrow; i++) {  
    for (j = 1; j <= ncol; j++) {  
      new[i][j] = old[i][j];
    }
  }
}
