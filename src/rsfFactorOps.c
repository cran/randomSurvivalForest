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

#include   "global.h"
#include   "nrutil.h"
#include   "rsfFactorOps.h"
extern uint  getTraceFlag();
#define FREE_ARG char*
#define NR_END 1
Factor **factorPtrVector(ulong nl, ulong nh) {
  Factor **v;
  v = (Factor **) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(Factor*)));     
  if (!v) nrerror("allocation failure in factorPtrVector()");
  return v-nl+NR_END;
}
void free_factorPtrVector(Factor **v,
                          ulong    nl,
                          ulong    nh) {
  free((FREE_ARG) (v+nl-NR_END));
}
Factor *makeFactor(uint r, char bookFlag) {
  uint i;
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nmakeFactor(%2x) ENTRY ...\n", bookFlag);
  }
  if (r < 2) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Minimum number of factor levels violated. ");
    Rprintf("\nRSF:  Requested %10d, Minimum Allowed %10d. ", r, 2);
    Rprintf("\nRSF:  The application will now exit. \n");
    exit(FALSE);
  }
  Factor *f = (Factor*) malloc((size_t)sizeof(Factor));
  f -> r = r;
  f -> cardinalGroupCount = (uint) floor(r/2);
  f -> mwcpSize = (r >> (3 + ulog2(SIZE_OF_INTEGER))) + ((r & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
  if (r <= MAX_EXACT_LEVEL) {
    f -> cardinalGroupSize = uivector(1, (f -> cardinalGroupCount) + 1);
    f -> complementaryPairCount =  ((uint*) (f -> cardinalGroupSize)) + (f -> cardinalGroupCount) + 1;
  }
  else {
    f -> cardinalGroupSize = dvector(1, (f -> cardinalGroupCount) + 1);
    f -> complementaryPairCount =  ((double*) (f -> cardinalGroupSize)) + (f -> cardinalGroupCount) + 1;
  }
  if (r <= MAX_EXACT_LEVEL) {
    *((uint*) f -> complementaryPairCount) = upower2(r-1) - 1;
  }
  else {
    *((double*) f -> complementaryPairCount) = pow(2, r-1) - 1;
  }
  for (i=1; i <= f -> cardinalGroupCount; i++) {
    if (r <= MAX_EXACT_LEVEL) {
      nChooseK(r, i, EXACT, ((uint*) f -> cardinalGroupSize) + i);
    }
    else {
      nChooseK(r, i, APROX, ((double*) f -> cardinalGroupSize) + i);
    }
    f -> cardinalGroupBinary = NULL;
  }
  if (!((f -> r) & 0x01)) {
    if (r <= MAX_EXACT_LEVEL) {
      ((uint*) f -> cardinalGroupSize)[f -> cardinalGroupCount] = ((uint*) f -> cardinalGroupSize)[f -> cardinalGroupCount] >> 1;
    }
    else {
      ((double*) f -> cardinalGroupSize)[f -> cardinalGroupCount] = ((double*) f -> cardinalGroupSize)[f -> cardinalGroupCount] / 2;
    }
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    if (r <= MAX_EXACT_LEVEL) {
      Rprintf("\n    Levels   GrpCount  PairCount\n");
      Rprintf("%10d %10d %12d \n", f -> r, f -> cardinalGroupCount, *((uint*) f -> complementaryPairCount));
      Rprintf("\n     Group       Size \n");
      for (i=1; i <= f -> cardinalGroupCount; i++) {
        Rprintf("%10d %12d \n", i, ((uint*) f -> cardinalGroupSize)[i]);
      }
    }
    else {
      Rprintf("\n    Levels   GrpCount  PairCount\n");
      Rprintf("%10d %10d %24.0f \n", f -> r, f -> cardinalGroupCount, *((double*) f -> complementaryPairCount));
      Rprintf("\n     Group       Size \n");
      for (i=1; i <= f -> cardinalGroupCount; i++) {
        Rprintf("%10d %24.0f \n", i, ((double*) f -> cardinalGroupSize)[i]);
      }
    }
  }
  if (bookFlag && (r <= MAX_EXACT_LEVEL)) {
    bookFactor(f);
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nmakeFactor() EXIT ...\n");
  }
  return f;
}
void free_Factor(Factor *f) {
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nfree_Factor(%4d) ENTRY ...\n", f -> r);
  }
  unBookFactor(f);
  if (f -> r <= MAX_EXACT_LEVEL) {
    free_uivector(f -> cardinalGroupSize, 1, (f -> cardinalGroupCount) + 1);
  }
  else {
    free_dvector(f -> cardinalGroupSize, 1, (f -> cardinalGroupCount) + 1);
  }
  free((FREE_ARG) f);
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nfree_Factor() EXIT ...\n");
  }
}
char bookFactor(Factor *f) {
  uint i, j;
  uint row;
  char result;
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nbookFactor(%4d) ENTRY ...\n", f -> r);
  }
  if (f -> r > MAX_EXACT_LEVEL) {
    Rprintf("\nRSF:  *** ERROR *** ");
    Rprintf("\nRSF:  Maximum number of factor levels violated in bookFactor(). ");
    Rprintf("\nRSF:  Requested %10d, Maximum Allowed %10d. ", f -> r, MAX_EXACT_LEVEL);
    Rprintf("\nRSF:  The application will now exit. \n");
    exit(FALSE);
  }
  if (f -> cardinalGroupBinary == NULL) {
    uint *leftLevel = uivector(1, f -> r);
    f -> cardinalGroupBinary = puivector(1, f -> cardinalGroupCount);
    for (i=1; i <= f -> cardinalGroupCount; i++) {
      (f -> cardinalGroupBinary)[i] = uivector(1, ((uint*) f -> cardinalGroupSize)[i]);
      row = 0;
      for (j = 1; j <= i; j++) {
        leftLevel[j] = 0;
      }
      if (getTraceFlag() & FACT_HGH_TRACE) {
        Rprintf("\nBooking Group:  %10d", i);
        Rprintf("\n     Index       Binary");
      }
      bookPair(f -> r , i, 1, &row, leftLevel, f);
      if (getTraceFlag() & FACT_HGH_TRACE) {
        Rprintf("\n");
      }
    }
    free_uivector(leftLevel, 1, f -> r);
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nbookFactor(%2x) EXIT ...\n", result);
  }
  return result;
}
char unBookFactor(Factor *f) {
  char result;
  uint i;
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nunBookFactor() ENTRY ...\n");
  }
  if (f -> cardinalGroupBinary != NULL) {
    for (i = 1; i <= f -> cardinalGroupCount; i++) {
      free_uivector((f -> cardinalGroupBinary)[i], 1, ((uint*) f -> cardinalGroupSize)[i]);
    }
    free_puivector(f -> cardinalGroupBinary, 1, f -> cardinalGroupCount);
    f -> cardinalGroupBinary = NULL;
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nunBookFactor(%2x) EXIT ...\n", result);
  }
  return result;
}
void bookPair (uint   levelCount, 
               uint    groupIndex, 
               uint    setColumn, 
               uint   *setRow, 
               uint   *daughter, 
               Factor *f) {
  uint i;
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nbookPair() ENTRY ...\n");
  }
  daughter[setColumn] ++;
  if (setColumn < groupIndex) {
    setColumn ++; 
    daughter[setColumn] ++;
    while (daughter[setColumn] < daughter[setColumn-1]) {
      daughter[setColumn] ++;
    }
    bookPair(levelCount, groupIndex, setColumn, setRow, daughter, f);
    daughter[setColumn] = 0;
    setColumn --;
    if ((*setRow) < ((uint*) (f -> cardinalGroupSize))[groupIndex]) {
      if (daughter[setColumn] < levelCount - (groupIndex - setColumn)) {
        bookPair(levelCount, groupIndex, setColumn, setRow, daughter, f);
      }
    }
  }
  else {
    (*setRow)++;
    (f -> cardinalGroupBinary)[groupIndex][*setRow] = 0;
    for (i=1; i <=groupIndex; i++) {
      (f -> cardinalGroupBinary)[groupIndex][*setRow] += upower(2, daughter[i] - 1);
    }
    if (getTraceFlag() & FACT_HGH_TRACE) {
      Rprintf("\n%10d", *setRow);
      Rprintf(", hex = %16x", (ulong) (f -> cardinalGroupBinary)[groupIndex][*setRow]);
    }
    if ( (levelCount > 2) && (daughter[setColumn] < levelCount)) {
      bookPair(levelCount, groupIndex, setColumn, setRow, daughter, f);
    }
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nbookPair() EXIT ...\n");
  }
}
void nChooseK (uint n, uint r, char type, void *result) {
  if (type == EXACT) {
    uint total, multiplier, divisor, newMultiplier, newDivisor, k;
    total = 1;
    divisor = 1;
    multiplier = n;
    k = ((r < (n-r)) ? r : (n-r));
    while(divisor <= k) {
      newMultiplier = multiplier;
      newDivisor = divisor;
      reduceFraction(& total, & newDivisor);
      reduceFraction(& newMultiplier, & newDivisor);
      if (newMultiplier > (UINT_MAX / total)) {
        Rprintf("\nRSF:  *** ERROR *** ");
        Rprintf("\nRSF:  Arithmetic Overflow Encountered in nChooseK(n, k). ");
        Rprintf("\nRSF:  Incoming parameters are (%10d, %10d). ", n, r);
        Rprintf("\nRSF:  The application will now exit. \n");
        exit(FALSE);
      }
      total = (total * newMultiplier) / newDivisor;
      multiplier--;
      divisor++;
    }
    *((uint*) result) = total;
  }
  else {
    double total, multiplier, divisor, k;
    total = 1;
    divisor = 1;
    multiplier = (double) n;
    k = (double) ((r < (n-r)) ? r : (n-r));
    while(divisor <= k) {
      total = (total * multiplier) / divisor;
      multiplier--;
      divisor++;
    }
    *((double*) result) = total;
  }
}
char reduceFraction(uint *numerator, uint *denominator) {
  uint numRemain, denRemain;
  char result;
  uint i;
  i = 2;
  result = FALSE;
  while (i <= *denominator) {
    numRemain = *numerator % i;
    if (numRemain == 0) {
      denRemain = *denominator % i;
      if (denRemain == 0) {
        *numerator = *numerator / i;
        *denominator = *denominator / i;
        result = TRUE;
      }
    }
    i++;
  }
  return result;
}
char splitOnFactor(uint level, uint *mwcp) {
  char daughterFlag;
  uint mwcpLevelIdentifier = (level >> (3 + ulog2(SIZE_OF_INTEGER))) + ((level & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
  uint mwcpLevelWord = upower(2, level - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  daughterFlag = RIGHT;
  if (mwcpLevelWord & mwcp[mwcpLevelIdentifier]) {
    daughterFlag = LEFT;
  }
  if (getTraceFlag() & FACT_HGH_TRACE) {
    Rprintf("\nMWCP Info:  level= %8d, word= %8d : %8x %8x --> ", level, mwcpLevelIdentifier, mwcpLevelWord, mwcp[mwcpLevelIdentifier]);
    if (daughterFlag == LEFT) {
      Rprintf("LEFT");
    }
    else {
      Rprintf("RGHT");
    }
  }
  return daughterFlag;
}
#undef NR_END
#undef FREE_ARG
