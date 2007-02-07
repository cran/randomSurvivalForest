##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 2.0.0
##
##  Copyright 2006, Cleveland Clinic
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

predict.rsf <- function(
    object = NULL,
    newdata = NULL,
    proximity = FALSE,
    do.trace = FALSE,
    ...)
{

    ### check that 'object' and 'newdata' are appropriate types
    if (is.null(object)) stop("Object is empty!")
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
        sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2  &
        sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
      if (is.null(object$forest)) 
        stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
      object <- object$forest
    }
    if (is.null(newdata)) stop("'newdata' is null.")
    if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.")
    
    ### check that newdata matches original training data
    ### add noise variable if necessary
    ### clean up newdata if necessary
    ### check to see if training data contains time and censoring info
    ### NOTE: automatic pass given if object is of class (rsf, partial)
    if (sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) == 2) {
      predictorsTestdata  <- as.matrix(newdata)
      rownames(predictorsTestdata) <- colnames(predictorsTestdata) <- NULL
      Time <- rep(0, dim(predictorsTestdata)[1])
      Cens <- rep(1, dim(predictorsTestdata)[1])
      performance <- F
    }
    else {
      fNames <- all.vars(object$formula)
      if (sum(is.element(fNames,"noise")) == 1) {
        newdata <- as.data.frame(cbind(newdata, noise = rnorm(dim(newdata)[1])))
      }
      if (sum(is.element(fNames[-(1:2)], names(newdata))) != length(fNames[-(1:2)])) {
        stop("'newdata' does not match original training data.")
      }
      if (sum(is.element(names(newdata), fNames[1:2])) == 2) {
        ftermLabels <- c(fNames[1:2], attr(terms(object$formula), "term.labels"))
      }
      else {
        ftermLabels <- attr(terms(object$formula), "term.labels")
      }
      newdata <- as.data.frame(
               model.matrix(as.formula(paste("~ -1 +",
               paste(ftermLabels, collapse="+"))), newdata))
      predictorsTestdata  <- as.matrix(newdata[,is.element(names(newdata), object$predictorNames)])
      rownames(predictorsTestdata) <- colnames(predictorsTestdata) <- NULL
      fNames.pt <- is.element(fNames[1:2], names(newdata))
      if (sum(fNames.pt) == 2) {
        Time <- newdata[,is.element(names(newdata), fNames[1])]
        Cens <- newdata[,is.element(names(newdata), fNames[2])]
        if (min(Cens) < 0 || max(Cens) > 1) {
          stop("censoring variable in newdata must be coded as 0 [censored] and 1 [death].")
        }
        performance <- T
      }
      else {
        Time <- rep(0, dim(predictorsTestdata)[1])
        Cens <- rep(1, dim(predictorsTestdata)[1])
        performance <- F
      }
    }
    
    ### work out individuals at risk (used later for mortality calculation)
    ### do.trace details
    tunique <- object$timeInterest
    ntree <- length(unique(object$nativeArray[,1]))
    Risk <- apply(cbind(1:length(tunique)),
                  1,
                  function(i, tau, tunq){sum(tau >= tunq[i])},
                  tau = object$Time, tunq = tunique)
    Risk <- Risk - c(Risk[-1],0)
    if (!is.logical(do.trace)) {
      if (do.trace >= 1){
        do.trace <- 256*round(do.trace) + 1
      }
      else {
        do.trace <- 0
      }
    }
    else {
      do.trace <- 1*(do.trace)
    }
 
    ###################################################################
    # Parameters passed the C function rsfPredict(...) are as follows:
    ###################################################################
    # 00 - C function name
    # 01 - trace output flag 
    #    - 0  = no trace output
    #    - !0 = various levels of trace output
    # 02 - memory useage protocol for return objects
    #    - any combination of the following are allowed
    #    - 0x00 = only the default objects are returned
    #    - 0x01 = return proximity information
    #    - 0x02 = N/A
    #    - 0x04 = return performance measure
    #    - 0x08 = N/A
    # 03 - number of bootstrap iterations
    #    - integer > 0
    # 04 - number of observations in GROW data set
    #    - integer > 1
    # 05 - vector of GROW observed times of death
    #    - vector of double values > 0 
    #    - optional, invalid values permitted
    # 06 - vector of GROW observed event types
    #    - 0 = censored
    #    - 1 = death
    #    - optional, invalid values permitted
    # 07 - number of GROW predictors 
    #    - integer > 0
    # 08 - [n x p] matrix of GROW predictor observations
    # 09 - number of observations in PRED data set
    #    - integer > 1
    # 10 - vector of PRED observed times of death
    #    - vector of double values > 0 
    # 11 - vector of PRED observed event types
    #    - 0 = censored
    #    - 1 = death
    # 12 - [fn x p] matrix of PRED predictor observations
    # 13 - number of time points of interest
    #    - integer > 0
    # 14 - vector of time points of interest
    #    - vector of double values
    # 15 - vector representing treeID
    # 16 - vector representing nodeID
    # 17 - vector representing parmID
    # 18 - vector representing spltPT
    # 19 - vector of GROW bootstrap random seeds
    #############################################################

    #############################################################
    # SEXP outputs (see native code for description):
    # Note that outputs depend on the MUP flags.
    #
    # fullEnsemble - default output
    # performance  - default output
    # leafCount    - default output
    # proximity    - optional
    #############################################################    

    nativeOutput <- .Call("rsfPredict",
        as.integer(do.trace),
        as.integer((4 * (if (performance) 1 else 0)) +
                   (if (proximity) 1 else 0)),
        as.integer(ntree),
        as.integer(dim(object$predictors)[1]),
        as.double(object$Time),
        as.integer(object$Cens),
        as.integer(dim(object$predictors)[2]),
        as.numeric(object$predictors),
        as.integer(dim(predictorsTestdata)[1]),
        as.double(Time),
        as.integer(Cens),
        as.numeric(predictorsTestdata),
        as.integer(length(object$timeInterest)),
        as.double(object$timeInterest),
        as.integer(object$nativeArray[,1]),
        as.integer(object$nativeArray[,2]),
        as.integer(object$nativeArray[,3]),
        as.double(object$nativeArray[,4]),
        as.integer(object$bootstrapSeed))

    mortality <- apply(matrix(nativeOutput$fullEnsemble, 
                       nrow = length(Cens),
                       byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt)},
                       wt = Risk)
        
    rsfOutput <- list(
        call = match.call(),
        forest = object,
        ntree = ntree,
        leaf.count = nativeOutput$leafCount,
        timeInterest = object$timeInterest,
        n = length(Cens),
        ndead = (if (performance) sum(Cens == 1) else NULL),
        Time = (if (performance) Time else NULL),
        Cens = (if (performance) Cens else NULL),
        predictorNames = object$predictorNames,
        predictors = predictorsTestdata,
        ensemble = matrix(nativeOutput$fullEnsemble, nrow = length(Cens), byrow = FALSE),
        mortality = (if (max(mortality) <= length(object$Cens)) mortality else
                        round(mortality*(length(object$Cens))/(1*(max(mortality) == 0)+max(mortality)))),
        err.rate = (if (performance) nativeOutput$performance else NULL),
        proximity = (if (proximity) nativeOutput$proximity else NULL)
    )
    class(rsfOutput) <- c("rsf", "predict")
    return(rsfOutput)

  }
