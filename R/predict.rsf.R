##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.5.1
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

predict.rsf <- function(
    object = NULL,
    test = NULL,
    importance = c("randomsplit", "permute", "none")[1],
    na.action = c("na.omit", "na.impute")[1],
    proximity = FALSE,
    seed = NULL,
    do.trace = FALSE,
    ...)
{

  ## Check that 'object' is of the appropriate type.
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2  &
      sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("rsf", "partial.rough"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'.")

  ## Acquire the forest and flag class (rsf, grow, bigdata)
  big.data <- F
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    if (inherits(object, c("rsf", "grow", "bigdata"), TRUE) [3] == 3) big.data <- T
    object <- object$forest
  }

  ## Check test data is a non-null data frame
  if (is.null(test)) stop("test data is null.")
  if (!is.data.frame(test)) stop("test data must be a data frame.")

  ## Ensure backward compatability with 3.0 user scripts.
  if (importance == TRUE)  importance <- "permute"
  if (importance == FALSE) importance <- "none"    

  ## Native code parameter for handling rough option
  predictedOutcomeBits <- 0

  ##
  ## Automatic pass given if object is of internal class (rsf, partial)
  ## Send void outcomes to the native code
  ##
  if (sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) == 2 |
      sum(inherits(object, c("rsf", "partial.rough"), TRUE) == c(1, 2)) == 2) {
    if (sum(inherits(object, c("rsf", "partial.rough"), TRUE) == c(1, 2)) == 2) {
      class(object) <- c("rsf", "partial")
      predictedOutcomeBits <- 2^14
    }
    test <- as.data.frame(data.matrix(test))
    predictorsTest <- as.matrix(test)
    rownames(predictorsTest) <- colnames(predictorsTest) <- NULL
    Time <- rep(-1, nrow(predictorsTest))
    Cens <- rep(-1, nrow(predictorsTest))
    performance <- F
    importance <- "none"
    remove(test)
  }
  
  ## Deterimine which grow variables are unordered/ordered factors
  ## Save factor levels of grow data for later reconversion
  ## Assign predictor type
  get.factor <- extract.factor(object$predictors)
  xfactor <- get.factor$xfactor
  xfactor.levels <- get.factor$xfactor.levels
  xfactor.order <- get.factor$xfactor.order
  xfactor.order.levels <- get.factor$xfactor.order.levels
  
  ##
  ## Details for all other rsf class objects:
  ##
  if (!(sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) == 2)) {

    ## Check if test data matches training (grow) data
    ## Check that factors are the same
    ## Confirm factor labels overlap
    ## Force test levels to equal grow factor levels
    ## (this last step is crucial to ensuring an immutable map)
    fNames <- all.vars(object$formula, max.names=1e7)
    test <- test[, is.element(names(test), fNames)]
    if (is.null(dim(test))) {
      test <- cbind(test)
      colnames(test) <- fNames[-c(1,2)]
    }
    get.factor <- extract.factor(test, fNames[-c(1,2)])
    xfactor.test <- get.factor$xfactor
    xfactor.order.test <- get.factor$xfactor.order
    
    if (!setequal(xfactor, xfactor.test))
      stop("(unordered) factors from test data do not match original training data")
    if (!setequal(xfactor.order, xfactor.order.test))
      stop("ordered factors from test data do not match original training data")
    if ((ncol(test) - sum(is.element(fNames[1:2], names(test)))) != (length(fNames) - 2)) {
      stop("test data does not match original training data")
    }
    if (length(xfactor) > 0) {
      for (k in 1:length(xfactor)) {
        xk.train <- object$predictors[ , names(object$predictors) == xfactor[k] ]
        xk.test.org <- xk.test <- xk.naomit.test <- test[ , names(test) == xfactor[k] ]
        if (na.action == "na.omit") {
          xk.naomit.test <- na.omit(test)[ , names(test) == xfactor[k] ]
        }
        if (length(na.omit(setdiff(unique(xk.naomit.test) , unique(xk.train)))) > 0)
          stop("(unordered) factor labels in test data differ from training data")
        levels(xk.test) <- xfactor.levels[[k]]
        for (l in 1:length(xk.test)) { xk.test[l] <- xk.test.org[l] }
        test[ , names(test) == xfactor[k] ] <- xk.test
      }
    }
    if (length(xfactor.order) > 0) {
      for (k in 1:length(xfactor.order)) {
        xk.train <- object$predictors[ , names(object$predictors) == xfactor.order[k] ]
        xk.test.org <- xk.test <- xk.naomit.test <- test[ , names(test) == xfactor.order[k] ]
        if (na.action == "na.omit") {
          xk.naomit.test <- na.omit(test)[ , names(test) == xfactor.order[k] ]
        }
        if (length(na.omit(setdiff(unique(xk.naomit.test) , unique(xk.train)))) > 0)
          stop("ordered factor labels in test data differ from training data")
        levels(xk.test) <- xfactor.order.levels[[k]]
        for (l in 1:length(xk.test)) { xk.test[l] <- xk.test.org[l] }
        test[ , names(test) == xfactor.order[k] ] <- xk.test
      }
    }

    ## Extract predictor and outcome names
    if (sum(is.element(names(test), fNames[1:2])) == 2) {
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

    ## Data conversion to numeric mode
    ## Add bogus row to overcome model.matrix strangeness with all NA's
    ## Special treatment for big.data
    if (na.action != "na.omit" & na.action != "na.impute") na.action <- "na.omit"
    if (!big.data) {
      test <- as.data.frame(data.matrix(test))
      old.na.action <- options()$na.action
      na.keep <- function(x){ x }
      na.takeOut <- function(x){ na.omit(x) }
      if (na.action == "na.omit") options(na.action=na.takeOut) else options(na.action=na.keep)
      all.na <- all(is.na(test))
      if (all.na) test <- rbind(0, test)
      test <- as.data.frame(model.matrix(as.formula(paste("~ -1 +",
                            paste(ftermLabels, collapse="+"))), test))
      if (all.na) test <- test[-1, ]
      options(na.action=old.na.action)
      if (nrow(test) == 0)
        stop("No records in the NA-processed test data.  Consider imputing using na.impute.")
    }
    else {#big.data
      test <- test[ , is.element(names(test), ftermLabels)]
      test <- as.data.frame(data.matrix(test))
      if (na.action == "na.omit") {
        test <- na.omit(test)
        if (nrow(test) == 0)
          stop("No records in the NA-processed test data.  Consider imputing using na.impute.")
      }
    }
    
    ## Extract test predictor matrix and sort columns
    ## Time and censoring processing
    ## If no outcomes present send void outcomes to native library
    ## Allow for all time and or censoring to be missing
    predictorsTest  <- as.matrix(test[,is.element(names(test), object$predictorNames)])
    if (length(object$predictorNames) > 1) {
      sort.col <- NULL
      for (k in 1:length(object$predictorNames)) {
        sort.col[k] <- which(colnames(predictorsTest) == object$predictorNames[k])
      }
      predictorsTest <- predictorsTest[ , sort.col]
    }    
    rownames(predictorsTest) <- colnames(predictorsTest) <- NULL

    if (sum(is.element(fNames[1:2], names(test))) == 2) {
      Time <- test[, is.element(names(test), fNames[1])]
      Cens <- test[, is.element(names(test), fNames[2])]
      if (!all(is.na(Cens))) {
        if (!all(is.element(unique(na.omit(Cens)), c(0,1)))) {
          stop("Censoring variable in test data must be NA, 0 or 1")
        }
      }
      if (!all(is.na(Time))) {
        if (!all(na.omit(Time) > 0)) stop("time must be strictly positive")
      }
      performance <- T
    }
    else {
      Time <- rep(-1, nrow(predictorsTest))
      Cens <- rep(-1, nrow(predictorsTest))
      performance <- F
      importance <- "none"
    }
    remove(test)

  }
  
  ## Data conversion for training predictors
  predictors <- as.matrix(data.matrix(object$predictors))
  rownames(predictors) <- colnames(predictors) <- NULL

  ## Set predictor type
  predictorType <- rep("R", ncol(predictorsTest))
  if (length(xfactor) > 0) {
    predictorType[is.element(object$predictorNames, xfactor)] <- "C"
  }
  if (length(xfactor.order) > 0) {
    predictorType[is.element(object$predictorNames, xfactor.order)] <- "I"
  }
  
  ## Work out individuals at risk (used later for mortality calculation)
  ## do.trace details
  tunique <- object$timeInterest
  ntree <- length(unique(object$nativeArray[,1]))
  Risk <- apply(cbind(1:length(tunique)),
                1,
                function(i, tau, tunq){sum(tau >= tunq[i])},
                tau = na.omit(object$time), tunq = tunique)
  Risk <- Risk - c(Risk[-1] , 0)

  ## seed details
  ## generate seed
  if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
  seed <- -round(abs(seed))

  if (!is.logical(do.trace)) {
    if (do.trace >= 1){
      do.trace <- 2^16 *round(do.trace) + 1
    }
    else {
      do.trace <- 0
    }
  }
  else {
    do.trace <- 1 * do.trace
  }

  ## Convert importance option into native code parameter.
  if (importance == "none") {
    importance.bits <- 0
  }
  else if (importance == "permute") {
    importance.bits <- 2048 + 0
  }
  else if (importance == "randomsplit") {
    importance.bits <- 2048 + 512
  }
  else {
    stop("Invalid choice for 'importance' option:  ", importance)
  }
  
  ####################################################################
  ## predict.rsf(...)
  ##
  ## Parameters passed:
  ## #################################################################
  ## 00 - C function call
  ##
  ## 01 - trace output flag
  ##    -  0 = no trace output
  ##    - !0 = various levels of trace output
  ##
  ## 02 - option protocol for output objects.
  ##    - only some of the following are relevant (*)
  ##      GROW PRED INTR
  ##       *         *    0x0001 = OENS
  ##       *    *         0x0002 = FENS
  ##       *    *    *    0x0004 = PERF
  ##       **   *         0x0008 = PROX
  ##       *    *    *    0x0010 = LEAF
  ##       **             0x0020 = TREE \ part of 
  ##       **             0x0040 = SEED / the forest 
  ##       ***  ***       0x0080 = MISS
  ##       ***       ***  0x0100 = OMIS
  ##       **   **   *    0x0200 = \  VIMP_TYPE
  ##       **   **   *    0x0400 =  | VIMP_JOIN
  ##                **    0x4000 =  | VIMP_APRX
  ##       **   **   *    0x0800 = /  VIMP
  ##       **             0x1000 = \  VUSE_TYPE
  ##       **             0x2000 = /  VUSE
  ##       *    *    *    0x4000 = POUT_TYPE
  ##
  ##       (*)   default  output
  ##       (**)  optional output
  ##       (***) default  output 
  ##             - dependent on data and potentially suppressed
  ##
  ## 03 - random seed for repeatability, integer < 0 
  ## 04 - number of bootstrap iterations, integer > 0
  ## 05 - number of observations in GROW data set, integer > 1
  ## 06 - vector of GROW observed times of death, doubles > 0 
  ##    - optional, invalid (negative) values permitted
  ## 07 - vector of GROW observed event types
  ##       0 = censored
  ##       1 = death
  ##    - optional, invalid (negative) values permitted
  ## 08 - number of GROW predictors, integer > 0
  ## 09 - [p x n] matrix of GROW predictor observations
  ## 10 - number of observations in PRED data set, integer > 1
  ## 11 - vector of PRED observed times of death, doubles > 0 
  ## 12 - vector of PRED observed event types
  ##       0 = censored
  ##       1 = death
  ## 13 - [fn x p] matrix of PRED predictor observations
  ## 14 - number of time points of interest, integer > 0
  ## 15 - vector of time points of interest, doubles > 0
  ## 16 - vector representing treeID
  ## 17 - vector representing nodeID
  ## 18 - vector representing parmID
  ## 19 - vector representing contPT
  ## 20 - vector representing mwcpSZ
  ## 21 - vector representing mwcpPT  
  ## 22 - GROW bootstrap random seed
  ## 23 - vector of predictor types
  ##      "R" - real value
  ##      "I" - integer value 
  ##      "C" - categorical 
  ## ###########################################################

  ## ###########################################################
  ##  SEXP outputs:  See parameter 2 above and note "*".
  ## ###########################################################    

  nativeOutput <- .Call("rsfPredict",
                        as.integer(do.trace),
                        as.integer(predictedOutcomeBits +
                                   importance.bits +
                                   (8 * (if (proximity) 1 else 0)) +
                                   (4 * (if (performance) 1 else 0))),
                        as.integer(seed),
                        as.integer(ntree),
                        as.integer(nrow(predictors)),
                        as.double(object$time),
                        as.double(object$cens),
                        as.integer(ncol(predictors)),
                        as.numeric(predictors),
                        as.integer(nrow(predictorsTest)),
                        as.double(Time),
                        as.double(Cens),
                        as.numeric(predictorsTest),
                        as.integer(length(object$timeInterest)),
                        as.double(object$timeInterest),
                        as.integer((object$nativeArray)$treeID),
                        as.integer((object$nativeArray)$nodeID),
                        as.integer((object$nativeArray)$parmID),
                        as.double((object$nativeArray)$contPT),
                        as.integer((object$nativeArray)$mwcpSZ),                        
                        as.integer(object$nativeFactorArray),                        
                        as.integer(object$seed),
                        as.character(predictorType))

  ## remove data converted predictors
  ## check for error return condition in the native code
  remove(predictors)
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }

  ## predict ensemble mortality (check if rough estimation in place)
  if (predictedOutcomeBits == 0) {
    mortality <- apply(matrix(nativeOutput$fullEnsemble, 
                              nrow = nrow(predictorsTest),
                              byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt , na.rm = TRUE)},
                       wt = Risk)
  }
  else {
    mortality <- nativeOutput$fullEnsemble
  }

  ## Check if there was missing data, and assign data if necessary.
  nMiss <- sum(is.na(Cens) | is.na(Time) | apply(predictorsTest, 1, function(x){any(is.na(x))}))
  if (nMiss > 0) {
    imputedData <- matrix(nativeOutput$imputation, nrow = nMiss, byrow = FALSE)
    imputedIndv <- imputedData[,1]
    imputedData <- as.matrix(imputedData[,-1])
    if (nMiss == 1) imputedData <- t(imputedData)
  }

  ## Add column names to test predictor matrix 
  ## Add names to importance values
  ## Add column names to imputed data
  predictorsTest <- as.data.frame(predictorsTest)
  colnames(predictorsTest) <- object$predictorNames
  if (importance != "none" & performance) {
    VIMP <- nativeOutput$importance-nativeOutput$performance[ntree]
    names(VIMP) <- object$predictorNames
  }
  else {
    VIMP <- NULL
  }
  if (nMiss > 0) {
    imputedData <- as.data.frame(imputedData)
    colnames(imputedData) <- c(fNames[2], fNames[1], object$predictorNames)
  }
  
  ## Map test predictor factors back to original values
  if (length(xfactor) > 0) {
    for (k in 1:length(xfactor)) {
      ptk <- (colnames(predictorsTest) == xfactor[k])
      xk.org <- xk <- factor(xfactor.levels[[k]][predictorsTest[ , ptk ]])
      levels(xk) <- xfactor.levels[[k]]
      for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      predictorsTest[ , ptk] <- xk
    }
  }
  if (length(xfactor.order) > 0) {
    for (k in 1:length(xfactor.order)) {
      ptk <- (colnames(predictorsTest) == xfactor.order[k])
      xk.org <- xk <- factor(xfactor.order.levels[[k]][predictorsTest[ , ptk ]], ordered = TRUE)
      levels(xk) <- xfactor.order.levels[[k]]
      for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      predictorsTest[ , ptk] <- xk
    }
  }

  ## Map imputed data factors back to original values
  if (length(xfactor) > 0 & nMiss > 0) {
    for (k in 1:length(xfactor)) {
      ptk <- (colnames(imputedData) == xfactor[k])
      xk.org <- xk <- factor(xfactor.levels[[k]][imputedData[ , ptk ]])
      levels(xk) <- xfactor.levels[[k]]
      for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      imputedData[ , ptk] <- xk
    }
  }
  if (length(xfactor.order) > 0 & nMiss > 0) {
    for (k in 1:length(xfactor.order)) {
      ptk <- (colnames(imputedData) == xfactor.order[k])
      xk.org <- xk <- factor(xfactor.order.levels[[k]][imputedData[ , ptk ]], ordered = TRUE)
      levels(xk) <- xfactor.order.levels[[k]]
      for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      imputedData[ , ptk] <- xk
    }
  }

  rsfOutput <- list(
                    call = match.call(),
                    forest = object,
                    ntree = ntree,
                    leaf.count = nativeOutput$leafCount,
                    timeInterest = object$timeInterest,
                    n = length(Cens),
                    ndead = (if (performance) sum(Cens == 1 , na.rm = T) else NULL),
                    time = (if (performance) Time else NULL),
                    cens = (if (performance) Cens else NULL),
                    predictorNames = object$predictorNames,
                    predictors = predictorsTest,
                    ensemble = if (predictedOutcomeBits == 0) 
                    matrix(nativeOutput$fullEnsemble, nrow = length(Cens), byrow = FALSE)
                    else NULL,
                    mortality = mortality,
                    err.rate = (if (performance) nativeOutput$performance else NULL),
                    importance = VIMP,
                    proximity = (if (proximity) nativeOutput$proximity else NULL),
                    imputedIndv = (if (nMiss>0) imputedIndv else NULL),
                    imputedData = (if (nMiss>0) imputedData else NULL)
                    )
  class(rsfOutput) <- c("rsf", "predict")
  return(rsfOutput)

}
