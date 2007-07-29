##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.0.0
##
##  Copyright 2007, Cleveland Clinic
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
    importance = TRUE,
    na.action = c("na.omit", "na.impute")[1],
    proximity = FALSE,
    seed = NULL,
    do.trace = FALSE,
    ...)
{

    ## check that 'object' and 'newdata' are appropriate types
    ## flag class (rsf, grow, bigdata)
    if (is.null(object)) stop("Object is empty!")
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
        sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2  &
        sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")
    big.data <- F
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
      if (is.null(object$forest)) 
        stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
      if (inherits(object, c("rsf", "grow", "bigdata"), TRUE) [3] == 3) big.data <- T
      object <- object$forest
    }
    if (is.null(newdata)) stop("'newdata' is null.")
    if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.")

    ## automatic pass given if object is of class (rsf, partial)
    ## send void outcomes to native library
    if (sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) == 2) {
      predictorsTest  <- as.matrix(newdata)
      rownames(predictorsTest) <- colnames(predictorsTest) <- NULL
      Time <- rep(-1, dim(predictorsTest)[1])
      Cens <- rep(-1, dim(predictorsTest)[1])
      performance <- F
      remove(newdata)
    }

    ## details for other rsf class objects
    ## check that newdata matches original training data
    ## add noise variable if necessary
    ## extract predictor and outcome names
    ## note special treatment for big.data
    if (!(sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) == 2)) {
      fNames <- all.vars(object$formula, max.names=1e7)
      if (sum(is.element(fNames,"noise")) == 1) {
        if (any(names(newdata) == "noise")) names(newdata)[names(newdata) == "noise"] <- "NOISE"
        newdata <- as.data.frame(cbind(newdata, noise = rnorm(dim(newdata)[1])))
      }
      if ((sum(is.element(fNames[-(1:2)], names(newdata))) != length(fNames[-(1:2)])) |
          (sum(is.element(names(newdata), fNames[-(1:2)])) != length(fNames[-(1:2)]))) {
        stop("'newdata' does not match original training data.")
      }
      if (sum(is.element(names(newdata), fNames[1:2])) == 2) {
        if (!big.data) {
          ftermLabels <- c(fNames[1:2], attr(terms(object$formula), "term.labels"))
        }
        else {
          ftermLabels <- fNames
        }
      }
      else {
        if (!big.data) {
          ftermLabels <- attr(terms(object$formula), "term.labels")
        }
        else {
          ftermLabels <- fNames[-c(1:2)]
        }
      }
      ## NA processing
      ## add bogus row to overcome model.matrix strangeness with all NA's 
      ## special treatment for big.data=T (coerce factors to real)
      ## remove records with all values NA
      ## removing NA's depletes the data: check if enough data left
      if (na.action != "na.omit" & na.action != "na.impute") na.action <- "na.omit"
      if (!big.data) { 
        old.na.action <- options()$na.action
        na.keep <- function(x){x}
        na.takeOut <- function(x){na.omit(x)}
        if (na.action == "na.omit") options(na.action=na.takeOut) else options(na.action=na.keep)
        newdata <- rbind(rep(-1, dim(newdata)[2]), newdata)
        newdata <- as.data.frame(
                     model.matrix(as.formula(paste("~ -1 +",
                     paste(ftermLabels, collapse="+"))), newdata))
        newdata <- newdata[-1, ]
        options(na.action=old.na.action)
        if (dim(newdata)[1] == 0)
            stop("No records in the NA-processed test data.  Consider imputing using na.impute.")
      }
      else {
        if (na.action == "na.omit") {
          newdata <- na.omit(newdata[,is.element(names(newdata), ftermLabels)])
          if (dim(newdata)[1] == 0)
            stop("No records in the NA-processed test data.  Consider imputing using na.impute.")
        }
        else {
         newdata <- newdata[,is.element(names(newdata), ftermLabels)]
        }
      }
      pt.row <- apply(newdata, 1, function(x){all(is.na(x))})
      newdata <- newdata[!pt.row, ]
      if (dim(newdata)[1] == 0)
        stop("No records in the NA-processed test data.  Consider imputing using na.impute.")
      ## extract the test predictor matrix
      ## time and censoring processing
      ## if no outcomes present send void outcomes to native library
      ## allow for all time and or censoring to be missing
      predictorsTest  <- as.matrix(newdata[,is.element(names(newdata), object$predictorNames)])
      rownames(predictorsTest) <- colnames(predictorsTest) <- NULL
      if (sum(is.element(fNames[1:2], names(newdata))) == 2) {
        Time <- newdata[,is.element(names(newdata), fNames[1])]
        Cens <- newdata[,is.element(names(newdata), fNames[2])]
        if (!all(is.na(Cens))) {
          if (!all(is.element(unique(na.omit(Cens)), c(0,1)))) {
              stop("censoring variable in 'newdata' can only be NA, 0 or 1.")
          }
        }
        if (!all(is.na(Time))) {
          if (!all(na.omit(Time) > 0)) stop("Time must be strictly positive.")
        }
        performance <- T
      }
      else {
        Time <- rep(-1, dim(predictorsTest)[1])
        Cens <- rep(-1, dim(predictorsTest)[1])
        performance <- F
      }
      remove(newdata)
    }
      
    ## determine how many records having missing data
    nMiss <- sum(is.na(Cens) | is.na(Time) | apply(predictorsTest, 1, function(x){any(is.na(x))}))
    
    ## set predictor types
    if (!big.data) {
      whole.number <- function(x) {
        n.whole.number <- max(10, round(0.25*length(x)))
        if (length(unique(x)) <= n.whole.number & all((na.omit(x) - floor(na.omit(x))) == 0)) "I" else "R"
      }
      predictorType <- apply(object$predictors, 2, whole.number)
    }
    else {
      predictorType <- rep("R", dim(object$predictors)[2])
    }

    ## work out individuals at risk (used later for mortality calculation)
    ## do.trace details
    tunique <- object$timeInterest
    ntree <- length(unique(object$nativeArray[,1]))
    Risk <- apply(cbind(1:length(tunique)),
                  1,
                  function(i, tau, tunq){sum(tau >= tunq[i])},
                  tau = object$time, tunq = tunique)
    Risk <- Risk - c(Risk[-1],0)

    ## seed details
    ## generate seed
    if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
    seed <- -round(abs(seed))

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

    ## #################################################################
    ## Parameters passed the C function rsfPredict(...) are as follows:
    ## #################################################################
    ## 00 - C function name
    ## 01 - trace output flag 
    ##    - 0  = no trace output
    ##    - !0 = various levels of trace output
    ## 02 - memory useage protocol for return objects
    ##    - any combination of the following are allowed
    ##    - 0x00 = only the default objects are returned
    ##    - 0x01 = return proximity information
    ##    - 0x02 = N/A
    ##    - 0x04 = return performance measure
    ##    - 0x08 = N/A
    ## 03 - random seed for repeatability
    ##    - integer < 0 
    ## 04 - number of bootstrap iterations
    ##    - integer > 0
    ## 05 - number of observations in GROW data set
    ##    - integer > 1
    ## 06 - vector of GROW observed times of death
    ##    - vector of double values > 0 
    ##    - optional, invalid values permitted
    ## 07 - vector of GROW observed event types
    ##    - 0 = censored
    ##    - 1 = death
    ##    - optional, invalid values permitted
    ## 08 - number of GROW predictors 
    ##    - integer > 0
    ## 09 - [n x p] matrix of GROW predictor observations
    ## 10 - number of observations in PRED data set
    ##    - integer > 1
    ## 11 - vector of PRED observed times of death
    ##    - vector of double values > 0 
    ## 12 - vector of PRED observed event types
    ##    - 0 = censored
    ##    - 1 = death
    ## 13 - [fn x p] matrix of PRED predictor observations
    ## 14 - number of time points of interest
    ##    - integer > 0
    ## 15 - vector of time points of interest
    ##    - vector of double values
    ## 16 - vector representing treeID
    ## 17 - vector representing nodeID
    ## 18 - vector representing parmID
    ## 19 - vector representing spltPT
    ## 20 - vector of GROW bootstrap random seeds
    ## 21 - vector of predictor types
    ##      "R" - real value
    ##      "I" - integer value 
    ##      "C" - categorical 
    ## ###########################################################

    ## ###########################################################
    ##  SEXP outputs (see native code for description):
    ## Note that outputs depend on the MUP flags.
    ##
    ## fullEnsemble - default output
    ## performance  - default output
    ## leafCount    - default output
    ## proximity    - optional
    ## ###########################################################    

    nativeOutput <- .Call("rsfPredict",
        as.integer(do.trace),
        as.integer((8 * (if (importance) 1 else 0)) +
                   (4 * (if (performance) 1 else 0)) +
                   (if (proximity) 1 else 0)),
        as.integer(seed),
        as.integer(ntree),
        as.integer(dim(object$predictors)[1]),
        as.double(object$time),
        as.double(object$cens),
        as.integer(dim(object$predictors)[2]),
        as.numeric(object$predictors),
        as.integer(dim(predictorsTest)[1]),
        as.double(Time),
        as.double(Cens),
        as.numeric(predictorsTest),
        as.integer(length(object$timeInterest)),
        as.double(object$timeInterest),
        as.integer(object$nativeArray[,1]),
        as.integer(object$nativeArray[,2]),
        as.integer(object$nativeArray[,3]),
        as.double(object$nativeArray[,4]),
        as.integer(object$seed),
        as.character(predictorType))

    ## check for error return condition in the native code.
    if(is.null(nativeOutput)) {
      stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
    }

    ## predict ensemble mortality
    mortality <- apply(matrix(nativeOutput$fullEnsemble, 
                       nrow = dim(predictorsTest)[1],
                       byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt)},
                       wt = Risk)

    ## check if there was missing data
    ## Assign imputed data
    if (nMiss > 0) {
      imputedData <- matrix(nativeOutput$imputation, nrow = nMiss, byrow = FALSE)
      imputedIndv <- imputedData[,1]
      imputedData <- as.matrix(imputedData[,-1])
      if (nMiss == 1) imputedData <- t(imputedData)
    }

    ## add column names to predictor matrix
    ## add column names to imputed data
    colnames(predictorsTest) <- object$predictorNames
    if (nMiss > 0) {
      colnames(imputedData) <- c(fNames[2], fNames[1], object$predictorNames)
    }

    rsfOutput <- list(
        call = match.call(),
        forest = object,
        ntree = ntree,
        leaf.count = nativeOutput$leafCount,
        timeInterest = object$timeInterest,
        n = length(Cens),
        ndead = (if (performance) sum(Cens == 1) else NULL),
        time = (if (performance) Time else NULL),
        cens = (if (performance) Cens else NULL),
        predictorNames = object$predictorNames,
        predictors = predictorsTest,
        ensemble = matrix(nativeOutput$fullEnsemble, nrow = length(Cens), byrow = FALSE),
        mortality = mortality,
        err.rate = (if (performance) nativeOutput$performance else NULL),
        importance = (if (importance) nativeOutput$importance-nativeOutput$performance[ntree]
                      else NULL),
        proximity = (if (proximity) nativeOutput$proximity else NULL),
        imputedIndv = (if (nMiss>0) imputedIndv else NULL),
        imputedData = (if (nMiss>0) imputedData else NULL)
    )
    class(rsfOutput) <- c("rsf", "predict")
    return(rsfOutput)

  }
