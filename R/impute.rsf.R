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

########################################################################
# 
# impute data using rsf
#
#
########################################################################

impute.rsf <- function(formula,
                       data = NULL,
                       ntree = 1000,
                       mtry = NULL,
                       nodesize = NULL,
                       splitrule = NULL,
                       nsplit = 0,
                       big.data = FALSE,
                       nimpute = 1,
                       predictorWt = NULL,
                       seed = NULL,
                       do.trace = FALSE,
                       ...)
{

  ## set parameters accordingly
  ## special data.frame class used to flag impute.only in rsf.default
  importance <- c("randomsplit", "permute", "none")[3]
  na.action <- c("na.omit", "na.impute")[2]
  forest <- FALSE
  proximity <- FALSE
  varUsed <- NULL
  split.depth <- FALSE
  class(data) <- c("data.frame", "impute.only")
  
  ## rsf grow call
  object <- rsf(formula = formula, data = data, ntree = ntree, mtry = mtry, nodesize = nodesize,
                splitrule = splitrule, nsplit = nsplit, big.data = big.data,
                nimpute = nimpute, predictorWt = predictorWt, seed = seed, do.trace = do.trace,
                importance = importance, na.action = na.action, forest = forest,
                proximity = proximity, varUsed = varUsed, split.depth = split.depth,
                impute.only = TRUE)

  ## interlay missing and non-missing data
  imputed.data <- cbind(cens = object$cens, time = object$time, object$predictors)

  if (!is.null(object$imputedIndv)) {
    imputed.data[object$imputedIndv, ] <- object$imputedData
  }

  colnames(imputed.data)[c(2,1)] <- all.vars(object$formula)[1:2]

  ## return the goodies 
  invisible(imputed.data) 

}
