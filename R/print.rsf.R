##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.2.0
##
##  Copyright 2008, Cleveland Clinic Foundation
##
##  This program is free software; you can redistribute it and/or
##  modify it under the terms of the GNU General Public License
##  as published by the Free Software Foundation; either version 2
##  of the License, or (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public
##  License along with this program; if not, write to the Free
##  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##  Boston, MA  02110-1301, USA.
##
##  Project funded by:
##    National Institutes of Health, HL072771-01
##
##    Michael S Lauer, MD, FACC, FAHA
##    Cleveland Clinic Lerner College of Medicine of CWRU
##    9500 Euclid Avenue
##    Cleveland, OH 44195
##
##    email:  lauerm@ccf.org
##    phone:   216-444-6798
##
##  Written by:
##    Hemant Ishwaran, Ph.D.
##    Dept of Quantitative Health Sciences/Wb4
##    Cleveland Clinic Foundation
##    9500 Euclid Avenue
##    Cleveland, OH 44195
##
##    email:  hemant.ishwaran@gmail.com
##    phone:  216-444-9932
##    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
##    --------------------------------------------------------------
##    Udaya B. Kogalur, Ph.D.
##    Kogalur Shear Corporation
##    5425 Nestleway Drive, Suite L1
##    Clemmons, NC 27012
##
##    email:  ubk2101@columbia.edu
##    phone:  919-824-9825
##    URL:    www.kogalur-shear.com
##
##**********************************************************************
##**********************************************************************

print.rsf <- function(x, ...) {

  ### set default printing if object is a forest
  if (sum(inherits(x, c("rsf", "forest"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }

  ### check that object is interpretable
  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")

  ### printing depends upon the object
  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    cat("\nCall:\n", deparse(x$call), "\n\n")
    cat("                         Sample size: ", x$n,                 "\n", sep="")
    cat("                    Number of deaths: ", x$ndead,             "\n", sep="")
    if (!is.null(x$imputedIndv)) {
      cat("                    Was data imputed: ", "yes",               "\n", sep="")
      cat("                         Missingness: ",
          round(100*length(x$imputedIndv)/x$n,2), "%\n", sep="")      
    }
    cat("                     Number of trees: ", x$ntree,             "\n",sep="")
    cat("          Minimum terminal node size: ", x$nodesize,          "\n", sep="")
    cat("       Average no. of terminal nodes: ", mean(x$leaf.count),  "\n", sep="")
    cat("No. of variables tried at each split: ", x$mtry,              "\n", sep="")
    cat("              Total no. of variables: ", length(x$predictorNames), "\n", sep="")  
    cat("                      Splitting rule: ", x$splitrule,         "\n", sep="")
    cat("              Estimate of error rate: ",
                             round(x$err.rate[x$ntree]*100, dig=2), "%\n\n", sep="")
  }
  else {
    cat("\nCall:\n", deparse(x$call), "\n\n")
    cat("  Sample size of test (predict) data: ", x$n,                 "\n", sep="")
    if (!is.null(x$ndead)) {
      cat("       Number of deaths in test data: ", x$ndead,             "\n", sep="")
    }
    if (!is.null(x$imputedData)) {
      cat("               Was test data imputed: ", "yes",               "\n", sep="")
      cat("                         Missingness: ",
          round(100*length(x$imputedIndv)/x$n,2), "%\n", sep="")      
    }
    cat("                Number of grow trees: ", x$ntree,             "\n",sep="")
    cat("  Average no. of grow terminal nodes: ", mean(x$leaf.count),  "\n", sep="")
    cat("         Total no. of grow variables: ", length(x$predictorNames), "\n", sep="")  
    if (!is.null(x$err.rate)) {
      if (sum(is.na(x$err.rate)) == 0) {
        cat("                     Test error rate: ",
            round(x$err.rate[x$ntree]*100, dig=2), "%\n\n", sep="")
      }
    }
  }
}
