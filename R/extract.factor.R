##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.5.0
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

extract.factor <- function(data, prednames=NULL) {
  predictorType  <- xfactor.levels <- xfactor.order.levels <- NULL
  xfactor <- names(data)[apply(cbind(1:ncol(data)), 1, function(k)
                  {is.factor(data[ , k]) && !is.ordered(data[ , k])})]
  xfactor.order <- names(data)[apply(cbind(1:ncol(data)), 1, function(k) {is.ordered(data[ , k])})]
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
    xfactor.order.levels <- apply(cbind(1:(1+length(xfactor.order))), 1, function(k) {
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