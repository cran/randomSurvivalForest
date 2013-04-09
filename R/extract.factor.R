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
