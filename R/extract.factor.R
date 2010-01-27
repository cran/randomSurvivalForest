####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.1
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

rm.na.levels <- function(data, prednames=NULL) {
  factor.names <- names(data)[sapply(1:ncol(data), function(k) {is.factor(data[ , k])})]
  if (!is.null(prednames)) factor.names <- intersect(factor.names , prednames)
  if (length(factor.names) > 0) {
    levels.na.pt <- sapply(1:length(factor.names), function(k) {
        any(levels(data[ , names(data) == factor.names[k]]) == "NA")
    })
    if (any(levels.na.pt)) {
      factor.names <- factor.names[levels.na.pt]
      for (k in 1:length(factor.names)) {
        x <- data[ , names(data) == factor.names[k]]
        levels(x)[levels(x) == "NA"]  <- NA
        data[ , names(data) == factor.names[k]] <- x
      }
    }
  }
  data
}  
  

extract.factor <- function(data, prednames=NULL) {
  predictorType  <- xfactor.levels <- xfactor.order.levels <- NULL
  xfactor <- names(data)[sapply(1:ncol(data), function(k)
                  {is.factor(data[ , k]) && !is.ordered(data[ , k])})]
  xfactor.order <- names(data)[sapply(1:ncol(data), function(k) {is.ordered(data[ , k])})]
  if (!is.null(prednames)) xfactor <- intersect(xfactor , prednames)
  if (!is.null(prednames)) xfactor.order <- intersect(xfactor.order , prednames)
  if (length(xfactor) > 0) {
    xfactor.levels <- apply(cbind(1:(1+length(xfactor))), 1, function(k) {
        if (k <= length(xfactor)) {
          levels(data[ , names(data) == xfactor[k]])
        }
        else {
          NULL
        }
    })
    xfactor.levels <- xfactor.levels[-(1+length(xfactor))]
  }
  if (length(xfactor.order) > 0 ) {
    xfactor.order.levels <- sapply(1:(1+length(xfactor.order)), function(k) {
        if (k <= length(xfactor.order)) {
          levels(data[ , names(data) == xfactor.order[k]])
        }
        else {
          NULL
        }
    })
    xfactor.order.levels <- xfactor.order.levels[-(1+length(xfactor.order))]
  }
  if (!is.null(prednames)) {
    predictorType <- rep("R", length(prednames))
    if (length(xfactor) > 0) {
      predictorType[is.element(prednames, xfactor)] <- "C"
    }
    if (length(xfactor.order) > 0) {
      predictorType[is.element(prednames, xfactor.order)] <- "I"
    }
  }
  return(list(
    xfactor=xfactor,
    xfactor.order=xfactor.order,
    xfactor.levels=xfactor.levels,
    xfactor.order.levels=xfactor.order.levels,
    predictorType=predictorType
    ))              
}
