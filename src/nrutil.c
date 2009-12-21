////**********************************************************************
////**********************************************************************
////
////  RANDOM SURVIVAL FOREST 3.6.0
////
////  Copyright 2009, Cleveland Clinic Foundation
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
////  ----------------------------------------------------------------
////  Project Partially Funded By:
////    --------------------------------------------------------------
////    National Institutes of Health,  Grant HHSN268200800026C/0001
////
////    Michael S. Lauer, M.D., FACC, FAHA 
////    National Heart, Lung, and Blood Institute
////    6701 Rockledge Dr, Room 10122
////    Bethesda, MD 20892
////
////    email:  lauerm@nhlbi.nih.gov
////
////    --------------------------------------------------------------
////    Case Western Reserve University/Cleveland Clinic  
////    CTSA Grant:  UL1 RR024989, National Center for
////    Research Resources (NCRR), NIH
////
////    --------------------------------------------------------------
////    Dept of Defense Era of Hope Scholar Award, Grant W81XWH0910339
////    Andy Minn, M.D., Ph.D.
////    Department of Radiation and Cellular Oncology, and
////    Ludwig Center for Metastasis Research
////    The University of Chicago, Jules F. Knapp Center, 
////    924 East 57th Street, Room R318
////    Chicago, IL 60637
//// 
////    email:  aminn@radonc.uchicago.edu
////
////    --------------------------------------------------------------
////    Bryan Lau, Ph.D.
////    Department of Medicine, Johns Hopkins School of Medicine,
////    Baltimore, Maryland 21287
////
////    email:  blau1@jhmi.edu
////
////  ----------------------------------------------------------------
////  Written by:
////    --------------------------------------------------------------
////    Hemant Ishwaran, Ph.D.
////    Dept of Quantitative Health Sciences/Wb4
////    Cleveland Clinic Foundation
////    9500 Euclid Avenue
////    Cleveland, OH 44195
////
////    email:  hemant.ishwaran@gmail.com
////    phone:  216-444-9932
////    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
////
////    --------------------------------------------------------------
////    Udaya B. Kogalur, Ph.D.
////    Dept of Quantitative Health Sciences/Wb4
////    Cleveland Clinic Foundation
////    
////    Kogalur Shear Corporation
////    5425 Nestleway Drive, Suite L1
////    Clemmons, NC 27012
////
////    email:  ubk2101@columbia.edu
////    phone:  919-824-9825
////    URL:    www.kogalur-shear.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************

#include      <stdlib.h>
#include      "nrutil.h"
#include <R_ext/Print.h>
#define NUMR_DEF_TRACE  0x002000  
extern unsigned int getTraceFlag();
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
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
unsigned int upower (unsigned int x, unsigned int n) {
  unsigned int p;
  if ((x >= 2) & (n > (sizeof(unsigned int) * 8))) {
    nrerror("Overflow in upower(), exponent too large.");
  }
  for (p = 1; n > 0; --n) {
    p = p * x;
  }
  return p;
}
unsigned int upower2 (unsigned int n) {
  unsigned int p;
  if (n > (sizeof(unsigned int) * 8)) {
    nrerror("Overflow in upower2(), exponent too large.");
  }
  p = 1 << n;
  return p;
}
unsigned int ulog2 (unsigned int n) {
  unsigned int p;
  p = 0;
  while (n > 1) {
    n = n >> 1;
    p++;
  }
  return p;
}
void hpsort(double *ra, unsigned int n) {
  unsigned int i, ir, j, l;
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
void hpsortui(unsigned int *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  unsigned int rra;
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
void indexx(unsigned int n, double *arr, unsigned int *indx) {
  unsigned int i, j, k, l=1;
  unsigned int indxt, itemp, ir=n;
  unsigned int *istack, jstack=0;
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
#define FREE_ARG char*
#define NR_END 1
double ****dmatrix4(unsigned long n4l,
                    unsigned long n4h,
                    unsigned long n3l, 
                    unsigned long n3h, 
                    unsigned long nrl, 
                    unsigned long nrh, 
                    unsigned long ncl, 
                    unsigned long nch) {
  unsigned long i, j, k, p;
  unsigned long n4   = n4h - n4l + 1;
  unsigned long n3   = n3h - n3l + 1; 
  unsigned long nrow = nrh - nrl + 1; 
  unsigned long ncol = nch - ncl + 1;
  double ****m;
  m=(double ****) malloc((size_t)((n4+NR_END)*sizeof(double***)));
  if (!m) nrerror("Allocation Failure 1 in dmatrix4()");
  m += NR_END;
  m -= n4l;
  m[n4l]=(double ***) malloc((size_t)((n4*n3+NR_END)*sizeof(double**)));
  if (!m[n4l]) nrerror("Allocation Failure 2 in dmatrix4()");
  m[n4l] += NR_END;
  m[n4l] -= n3l;
  m[n4l][n3l]=(double **) malloc((size_t)((n4*n3*nrow+NR_END)*sizeof(double*)));
  if (!m[n4l][n3l]) nrerror("Allocation Failure 3 in dmatrix4()");
  m[n4l][n3l] += NR_END;
  m[n4l][n3l] -= nrl;
  m[n4l][n3l][nrl]=(double *) malloc((size_t)((n4*n3*nrow*ncol+NR_END)*sizeof(double)));
  if (!m[n4l][n3l][nrl]) nrerror("Allocation Failure 4 in dmatrix4()");
  m[n4l][n3l][nrl] += NR_END;
  m[n4l][n3l][nrl] -= ncl;
  for(i=n4l+1;i<=n4h;i++) {
    m[i] = m[i-1] + n3;
  }
  p = 0;
  for(i=n4l;i<=n4h;i++) {
    for(j=n3l;j<=n3h;j++) {
      m[i][j] = m[n4l][n3l] + nrow * p;
      p++;
    }
  }
  p = 0;
  for(i=n4l;i<=n4h;i++) {
    for(j=n3l;j<=n3h;j++) {
      for(k=nrl;k<=nrh;k++) {
        m[i][j][k] = m[n4l][n3l][nrl] + ncol * p;
        p++;
      }
    }
  }
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndmatrix4() alloc:  %20x  \n", m);
  }
  return m;
}
double ***dmatrix3(unsigned long n3l, 
                   unsigned long n3h, 
                   unsigned long nrl, 
                   unsigned long nrh, 
                   unsigned long ncl, 
                   unsigned long nch) {
  unsigned long i, j, k;  
  unsigned long n3   = n3h - n3l + 1;
  unsigned long nrow = nrh - nrl + 1; 
  unsigned long ncol = nch - ncl + 1;
  double ***m;
  m=(double ***) malloc((size_t)((n3+NR_END)*sizeof(double**)));
  if (!m) nrerror("Allocation Failure 1 in dmatrix3()");
  m += NR_END;
  m -= n3l;
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
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndmatrix3() alloc:  %20x  \n", m);
  }
  return m;
}
double **dmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  unsigned long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
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
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndmatrix() alloc:  %20x  \n", m);
  }
  return m;
}
int **imatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  unsigned long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nimatrix() alloc:  %20x  \n", m);
  }
  return m;
}
unsigned int **uimatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  unsigned long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  unsigned int **m;
  m=(unsigned int **) malloc((size_t)((nrow+NR_END)*sizeof(unsigned int*)));
  if (!m) nrerror("Allocation Failure 1 in uimatrix()");
  m += NR_END;
  m -= nrl;
  m[nrl]=(unsigned int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned int)));
  if (!m[nrl]) nrerror("Allocation Failure 2 in uimatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nuimatrix() alloc:  %20x  \n", m);
  }
  return m;
}
unsigned long **ulmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  unsigned long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  unsigned long **m;
  m=(unsigned long **) malloc((size_t)((nrow+NR_END)*sizeof(unsigned long*)));
  if (!m) nrerror("Allocation Failure 1 in ulmatrix()");
  m += NR_END;
  m -= nrl;
  m[nrl]=(unsigned long *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned long)));
  if (!m[nrl]) nrerror("Allocation Failure 2 in ulmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nulmatrix() alloc:  %20x  \n", m);
  }
  return m;
}
char **cmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  unsigned long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ncmatrix() alloc:  %20x  \n", m);
  }
  return m;
}
double *dvector(unsigned long nl, unsigned long nh) {
  double *v;
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation Failure in dvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndvector() alloc:  %20x  \n", v);
  }
  return v;
}
double **pdvector(unsigned long nl, unsigned long nh) {
  double **v;
  v=(double **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double*)));
  if (!v) nrerror("Allocation Failure in pdvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npdvector() alloc:  %20x  \n", v);
  }
  return v;
}
double ***ppdvector(unsigned long nl, unsigned long nh) {
  double ***v;
  v=(double ***) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double**)));
  if (!v) nrerror("Allocation Failure in ppdvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nppdvector() alloc:  %20x  \n", v);
  }
  return v;
}
int *ivector(unsigned long nl, unsigned long nh) {
  int *v;
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("Allocation Failure in ivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nivector() alloc:  %20x  \n", v);
  }
  return v;
}
int **pivector(unsigned long nl, unsigned long nh) {
  int **v;
  v=(int **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int*)));
  if (!v) nrerror("Allocation Failure in puivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npivector() alloc:  %20x  \n", v);
  }
  return v;
}
unsigned int *uivector(unsigned long nl, unsigned long nh) {
  unsigned int *v;
  v=(unsigned int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned int)));
  if (!v) nrerror("Allocation Failure in uivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nuivector() alloc:  %20x  \n", v);
  }
  return v;
}
unsigned int **puivector(unsigned long nl, unsigned long nh) {
  unsigned int **v;
  v=(unsigned int **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned int*)));
  if (!v) nrerror("Allocation Failure in puivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npuivector() alloc:  %20x  \n", v);
  }
  return v;
}
unsigned int ***ppuivector(unsigned long nl, unsigned long nh) {
  unsigned int ***v;
  v=(unsigned int ***) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned int**)));
  if (!v) nrerror("Allocation Failure in ppuivector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nppuivector() alloc:  %20x  \n", v);
  }
  return v;
}
unsigned long *ulvector(unsigned long nl, unsigned long nh) {
  unsigned long *v;
  v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned long)));
  if (!v) nrerror("Allocation Failure in ulvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nulvector() alloc:  %20x  \n", v);
  }
  return v;
}
unsigned long **pulvector(unsigned long nl, unsigned long nh) {
  unsigned long **v;
  v=(unsigned long **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned long*)));
  if (!v) nrerror("Allocation Failure in pulvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npulvector() alloc:  %20x  \n", v);
  }
  return v;
}
char *cvector(unsigned long nl, unsigned long nh) {
  char *v;
  v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char)));
  if (!v) nrerror("Allocation Failure in cvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ncvector() alloc:  %20x  \n", v);
  }
  return v;
}
char **pcvector(unsigned long nl, unsigned long nh) {
  char **v;
  v=(char **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char*)));
  if (!v) nrerror("Allocation Failure in pcvector()");
  v = v-nl+NR_END;
  if (getTraceFlag() & NUMR_DEF_TRACE) {
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
  exit(TRUE);
}
void free_dmatrix4(double    ****m,
                   unsigned long n4l,
                   unsigned long n4h,
                   unsigned long n3l,
                   unsigned long n3h,
                   unsigned long nrl,
                   unsigned long nrh,
                   unsigned long ncl,
                   unsigned long nch
                   ) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndmatrix4() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[n4l][n3l][nrl]+ncl-NR_END));
  free((FREE_ARG) (m[n4l][n3l]+nrl-NR_END));
  free((FREE_ARG) (m[n4l]+n3l-NR_END));
  free((FREE_ARG) (m+n4l-NR_END));
}
void free_dmatrix3(double ***m,
                  unsigned long n3l,
                  unsigned long n3h,
                  unsigned long nrl,
                  unsigned long nrh,
                  unsigned long ncl,
                  unsigned long nch
                 ) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndmatrix3() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[n3l][nrl]+ncl-NR_END));
  free((FREE_ARG) (m[n3l]+nrl-NR_END));
  free((FREE_ARG) (m+n3l-NR_END));
}
void free_dmatrix(double **m,
                  unsigned long nrl,
                  unsigned long nrh,
                  unsigned long ncl,
                  unsigned long nch
                 ) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndmatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_imatrix(int **m,
                  unsigned long nrl,
                  unsigned long nrh,
                  unsigned long ncl,
                  unsigned long nch
                  ) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nimatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_uimatrix(unsigned int **m,
                   unsigned long nrl,
                   unsigned long nrh,
                   unsigned long ncl,
                   unsigned long nch
                  ) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nuimatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_ulmatrix(unsigned long **m,
                   unsigned long nrl,
                   unsigned long nrh,
                   unsigned long ncl,
                   unsigned long nch
                  ) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nulmatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_cmatrix(char **m,
                  unsigned long nrl,
                  unsigned long nrh,
                  unsigned long ncl,
                  unsigned long nch
                  ) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ncmatrix() de-allocating:  %20x  \n", m);
  }
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
void free_dvector(double *v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ndvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_pdvector(double **v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npdvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_ppdvector(double ***v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nppdvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_ivector(int *v,unsigned long nl,unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_pivector(int **v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_uivector(unsigned int *v,unsigned long nl,unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nuivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_puivector(unsigned int **v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npuivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_ppuivector(unsigned int ***v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nppuivector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_ulvector(unsigned long *v,unsigned long nl,unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\nulvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_pulvector(unsigned long **v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npulvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_cvector(char *v,unsigned long nl,unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\ncvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
void free_pcvector(char **v, unsigned long nl, unsigned long nh) {
  if (getTraceFlag() & NUMR_DEF_TRACE) {
    Rprintf("\npcvector() de-allocating:  %20x  \n", v);
  }
  free((FREE_ARG) (v+nl-NR_END));
}
#undef FREE_ARG
#undef NR_END
void nrCopyMatrix(unsigned int **new, unsigned int **old, unsigned int nrow, unsigned int ncol) {
  unsigned int i,j;
  for (i = 1; i <= nrow; i++) {  
    for (j = 1; j <= ncol; j++) {  
      new[i][j] = old[i][j];
    }
  }
}
void nrCopyVector(char *new, char *old, unsigned int ncol) {
  unsigned int j;
  for (j = 1; j <= ncol; j++) {  
    new[j] = old[j];
  }
}
void testEndianness() {
  unsigned int     test = 0x12345678;
  unsigned int *testPtr = & test;
  Rprintf("\n %2x %2x %2x %2x \n", 
          *((char *) testPtr), 
          *((char *) testPtr + 1), 
          *((char *) testPtr + 2), 
          *((char *) testPtr + 3));
}
