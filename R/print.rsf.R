####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.2
####
####  Copyright 2009, Cleveland Clinic Foundation
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
####  ----------------------------------------------------------------
####  Project Partially Funded By:
####    --------------------------------------------------------------
####    National Institutes of Health,  Grant HHSN268200800026C/0001
####
####    Michael S. Lauer, M.D., FACC, FAHA 
####    National Heart, Lung, and Blood Institute
####    6701 Rockledge Dr, Room 10122
####    Bethesda, MD 20892
####
####    email:  lauerm@nhlbi.nih.gov
####
####    --------------------------------------------------------------
####    Case Western Reserve University/Cleveland Clinic  
####    CTSA Grant:  UL1 RR024989, National Center for
####    Research Resources (NCRR), NIH
####
####    --------------------------------------------------------------
####    Dept of Defense Era of Hope Scholar Award, Grant W81XWH0910339
####    Andy Minn, M.D., Ph.D.
####    Department of Radiation and Cellular Oncology, and
####    Ludwig Center for Metastasis Research
####    The University of Chicago, Jules F. Knapp Center, 
####    924 East 57th Street, Room R318
####    Chicago, IL 60637
#### 
####    email:  aminn@radonc.uchicago.edu
####
####    --------------------------------------------------------------
####    Bryan Lau, Ph.D.
####    Department of Medicine, Johns Hopkins School of Medicine,
####    Baltimore, Maryland 21287
####
####    email:  blau1@jhmi.edu
####
####  ----------------------------------------------------------------
####  Written by:
####    --------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Dept of Quantitative Health Sciences/Wb4
####    Cleveland Clinic Foundation
####    9500 Euclid Avenue
####    Cleveland, OH 44195
####
####    email:  hemant.ishwaran@gmail.com
####    phone:  216-444-9932
####    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
####
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Dept of Quantitative Health Sciences/Wb4
####    Cleveland Clinic Foundation
####    
####    Kogalur Shear Corporation
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk2101@columbia.edu
####    phone:  919-824-9825
####    URL:    www.kogalur-shear.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

print.rsf <- function(x, ...) {

  ### set default printing if object is a forest
  if (sum(inherits(x, c("rsf", "forest"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }

  #
  
  # number of event types
  if (!is.null(x$cens)) {
    events  <- na.omit(x$cens)[na.omit(x$cens) > 0]
    n.event <- length(unique(events))
    event.freq <- paste(tapply(events, events, length), collapse = ", ")
  }
  else {
    n.event <- 1
  }

  # error rates 
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
    
  
  ### ensure backward compatibility for nsplit
  if (is.null(x$nsplit)) x$nsplit <- 0

  ### check that object is interpretable
  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")

  ### printing depends upon the object
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
