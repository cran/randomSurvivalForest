####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.4
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  Written by:
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  commerce@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

print.rsf <- function(x, ...) {
  if (sum(inherits(x, c("rsf", "forest"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  if (!is.null(x$cens)) {
    events  <- na.omit(x$cens)[na.omit(x$cens) > 0]
    n.event <- length(unique(events))
    event.freq <- paste(tapply(events, events, length), collapse = ", ")
  }
  else {
    n.event <- 1
  }
  if (!is.null(x$err.rate)) {
    x$err.rate <- rbind(x$err.rate)
    err.rate <- round(100*x$err.rate[, ncol(x$err.rate)], 2) 
    if (!all(is.na(err.rate))) {
      err.rate <- paste(err.rate, "%", collapse=", ", sep = "")
    }
    else {
      err.rate <- NULL
    }
  }
  else {
    err.rate <- NULL
  }
  if (is.null(x$nsplit)) x$nsplit <- 0
  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")
  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    cat("\nCall:\n", deparse(x$call), "\n\n")
    cat("                         Sample size: ", x$n,                 "\n", sep="")
    if (n.event == 1) {
      cat("                    Number of deaths: ", x$ndead,             "\n", sep="")
    }
    else {
      cat("                    Number of events: ", event.freq,          "\n", sep="")
    }
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
    if (x$nsplit > 0 & x$splitrule != "random") {
      cat("                      Splitting rule: ", paste(x$splitrule,"*random*"),         "\n", sep="")
      cat("       Number of random split points: ", x$nsplit                   ,         "\n", sep="")
    }
    else {
      cat("                      Splitting rule: ", x$splitrule,         "\n", sep="")
    } 
    cat("              Estimate of error rate: ", err.rate,            "\n\n", sep="")
  }
  else {
    cat("\nCall:\n", deparse(x$call), "\n\n")
    cat("  Sample size of test (predict) data: ", x$n,                 "\n", sep="")
    if (!is.null(x$cens)) {
      if (n.event == 1) {
        cat("       Number of deaths in test data: ", x$ndead,             "\n", sep="")
      }
      else {
        cat("       Number of events in test data: ", event.freq,           "\n", sep="")
      }
    }
    if (!is.null(x$imputedData)) {
      cat("               Was test data imputed: ", "yes",               "\n", sep="")
      cat("                         Missingness: ",
          round(100*length(x$imputedIndv)/x$n,2), "%\n", sep="")      
    }
    cat("                Number of grow trees: ", x$ntree,             "\n",sep="")
    cat("  Average no. of grow terminal nodes: ", mean(x$leaf.count),  "\n", sep="")
    cat("         Total no. of grow variables: ", length(x$predictorNames), "\n", sep="")  
    if (!is.null(err.rate)) {
      cat("                     Test error rate: ", err.rate, "\n\n", sep="")
    }
  }
}
