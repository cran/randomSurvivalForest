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

impute.rsf <- function(formula,#
                       data = NULL,#
                       ntree = 1000,#
                       mtry = NULL,#
                       nodesize = NULL,#
                       splitrule = NULL,#
                       nsplit = 0,#
                       big.data = FALSE,#
                       nimpute = 1,#
                       predictorWt = NULL,#
                       seed = NULL,#
                       do.trace = FALSE,#
                       ...)
{
  importance <- c("randomsplit", "permute", "none")[3]
  na.action <- c("na.omit", "na.impute")[2]
  forest <- FALSE
  proximity <- FALSE
  varUsed <- NULL
  split.depth <- FALSE
  class(data) <- c("data.frame", "impute.only")
  object <- rsf(formula = formula, data = data, ntree = ntree, mtry = mtry, nodesize = nodesize,
                splitrule = splitrule, nsplit = nsplit, big.data = big.data,
                nimpute = nimpute, predictorWt = predictorWt, seed = seed, do.trace = do.trace,
                importance = importance, na.action = na.action, forest = forest,
                proximity = proximity, varUsed = varUsed, split.depth = split.depth,
                impute.only = TRUE)
  imputed.data <- cbind(cens = object$cens, time = object$time, object$predictors)
  if (!is.null(object$imputedIndv)) {
    imputed.data[object$imputedIndv, ] <- object$imputedData
  }
  colnames(imputed.data)[c(2,1)] <- all.vars(object$formula)[1:2]
  invisible(imputed.data) 
}
