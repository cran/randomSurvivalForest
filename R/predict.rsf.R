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

predict.rsf <- function(
    object = NULL,#
    test = NULL,#
    importance = c("randomsplit", "permute", "none")[1],#
    na.action = c("na.omit", "na.impute")[1],#
    outcome = c("train", "test")[1],#
    proximity = FALSE,#
    split.depth = FALSE,#
    seed = NULL,#
    do.trace = FALSE,#
    ...)
{
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2  &
      sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("rsf", "partial.rough"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'")
  big.data <- F
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    if (inherits(object, c("rsf", "grow", "bigdata"), TRUE) [3] == 3) big.data <- T
    object <- object$forest
  }
  if (is.null(test)) stop("test data is null.")
  if (!is.data.frame(test)) stop("test data must be a data frame.")
  if (length(which(outcome == c("train", "test"))) != 1) stop("invalid outcome specified:", outcome) 
  if (importance == TRUE)  importance <- "permute"
  if (importance == FALSE) importance <- "none"    
  predictedOutcomeBits <- 0
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
    performance <- FALSE
    importance <- "none"
    remove(test)
  }
  get.factor <- extract.factor(object$predictors)
  xfactor <- get.factor$xfactor
  xfactor.levels <- get.factor$xfactor.levels
  xfactor.order <- get.factor$xfactor.order
  xfactor.order.levels <- get.factor$xfactor.order.levels
  if (!(sum(inherits(object, c("rsf", "partial"), TRUE) == c(1, 2)) == 2)) {
    fNames <- all.vars(object$formula, max.names=1e7)
    test <- test[, is.element(names(test), fNames)]
    if (is.null(dim(test))) {
      test <- cbind(test)
      colnames(test) <- fNames[-c(1,2)]
    }
    test <- rm.na.levels(test, fNames[-c(1,2)])
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
    else { 
      test <- test[ , is.element(names(test), ftermLabels)]
      test <- as.data.frame(data.matrix(test))
      if (na.action == "na.omit") {
        test <- na.omit(test)
        if (nrow(test) == 0)
          stop("No records in the NA-processed test data.  Consider imputing using na.impute.")
      }
    }
    predictorsTest  <- as.matrix(test[,is.element(names(test), object$predictorNames)])
    if (length(object$predictorNames) > 1) {
      sort.col <- NULL
      for (k in 1:length(object$predictorNames)) {
        sort.col[k] <- which(colnames(predictorsTest) == object$predictorNames[k])
      }
      if (nrow(predictorsTest) > 1) {
        predictorsTest <- predictorsTest[ , sort.col]
      }
      else {
        predictorsTest <- rbind(predictorsTest[ , sort.col])
      }
    }    
    rownames(predictorsTest) <- colnames(predictorsTest) <- NULL
    if (sum(is.element(fNames[1:2], names(test))) == 2) {
      Time <- test[, is.element(names(test), fNames[1])]
      Cens <- test[, is.element(names(test), fNames[2])]
      eventType <- unique(na.omit(Cens))
      if (sum(eventType >= 0) != length(eventType)) {
        stop("censoring variable must be coded as NA, 0, or greater than 0")    
      }
      if (!all(is.na(Time))) {
        if (!all(na.omit(Time) > 0)) stop("time must be strictly positive")
      }
      performance <- TRUE
      if (length(setdiff(na.omit(Cens), na.omit(object$cens))) > 1) {
        stop("event types in test data do not match training data")
      }
    }
    else {
      Time <- rep(-1, nrow(predictorsTest))
      Cens <- rep(-1, nrow(predictorsTest))
      performance <- FALSE
      importance <- "none"
    }
    remove(test)
  }
  n.event <- length(unique(na.omit(object$cens)[na.omit(object$cens) > 0]))
  if (n.event > 1) n.event <- n.event + 1
  if (outcome == "train") {
    predictors <- as.matrix(data.matrix(object$predictors))
    time <- object$time
    cens <- object$cens
  }
  else {
    if (!performance) stop("survival data missing in the test set: outcome='test' cannot be used")
    predictors <- predictorsTest
    time <- Time
    cens <- Cens
    newobject <- object
    newobject$time <- time
    newobject$cens <- cens
    newobject$predictors <- predictors
    newobject$outcome <- "test"
    vimp.out <- vimp(newobject, joint = FALSE)
    rm(newobject)
    VIMP <- vimp.out$importance
    ERR <- vimp.out$err.rate
    if (n.event > 1) ERR <- cbind(vimp.out$err.rate)
  }
  rownames(predictors) <- colnames(predictors) <- NULL
  predictorType <- rep("R", ncol(predictorsTest))
  if (length(xfactor) > 0) {
    predictorType[is.element(object$predictorNames, xfactor)] <- "C"
  }
  if (length(xfactor.order) > 0) {
    predictorType[is.element(object$predictorNames, xfactor.order)] <- "I"
  }
  tunique <- object$timeInterest
  ntree <- length(unique(object$nativeArray[,1]))
  Risk <- apply(cbind(1:length(tunique)),
                1,
                function(i, tau, tunq){sum(tau >= tunq[i])},
                tau = na.omit(time), tunq = tunique)
  Risk <- Risk - c(Risk[-1] , 0)
  if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
  seed <- -round(abs(seed))
  if (!is.logical(do.trace)) {
    if (do.trace >= 1){
      do.trace <- 2^24 *round(do.trace) + 1
    }
    else {
      do.trace <- 0
    }
  }
  else {
    do.trace <- 1 * do.trace
  }
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
  if (outcome == "test") {
    importance.bits <- 0
  }
  nativeOutput <- .Call("rsfPredict",
                        as.integer(do.trace),
                        as.integer(predictedOutcomeBits +
                                   importance.bits +
                                   (8 * (if (proximity) 1 else 0)) +
                                   (4 * (if (performance) 1 else 0)) + 
                                   (2^15 * (if (split.depth) 1 else 0)) +
                                   (2^17 * (if (outcome == "test") 1 else 0))),
                        as.integer(seed),
                        as.integer(ntree),
                        as.integer(nrow(predictors)),
                        as.double(time),
                        as.double(cens),
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
  remove(predictors)
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }
  nMiss <- sum(is.na(Cens) | is.na(Time) | apply(predictorsTest, 1, function(x){any(is.na(x))}))
  if (nMiss > 0) {
    imputedData <- matrix(nativeOutput$imputation, nrow = nMiss, byrow = FALSE)
    imputedIndv <- imputedData[,1]
    imputedData <- as.matrix(imputedData[,-1])
    if (nMiss == 1) imputedData <- t(imputedData)
  }
  predictorsTest <- as.data.frame(predictorsTest)
  colnames(predictorsTest) <- object$predictorNames
  if (nMiss > 0) {
    imputedData <- as.data.frame(imputedData)
    colnames(imputedData) <- c(fNames[2], fNames[1], object$predictorNames)
  }
  if (length(xfactor) > 0) {
    for (k in 1:length(xfactor)) {
      ptk <- (colnames(predictorsTest) == xfactor[k])
      factor.k <- xfactor.levels[[k]][predictorsTest[ , ptk ]]
      labels.k <- xfactor.levels[[k]][sort(unique(predictorsTest[ , ptk ]))]
      xk.org <- xk <- factor(factor.k, labels = labels.k, levels = labels.k, exclude = NULL)  
      if (length(setdiff(xfactor.levels[[k]], labels.k)) > 0) {
        levels(xk) <- xfactor.levels[[k]]
        for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      }
      predictorsTest[ , ptk] <- xk
    }
  }
  if (length(xfactor.order) > 0) {
    for (k in 1:length(xfactor.order)) {
      ptk <- (colnames(predictorsTest) == xfactor.order[k])
      factor.k <- xfactor.order.levels[[k]][predictorsTest[ , ptk ]]
      labels.k <- xfactor.order.levels[[k]][sort(unique(predictorsTest[ , ptk ]))]
      xk.org <- xk <- factor(factor.k, labels = labels.k, levels = labels.k, exclude = NULL, ordered = TRUE)
      if (length(setdiff(xfactor.order.levels[[k]], labels.k)) > 0) {
        levels(xk) <- xfactor.order.levels[[k]]
        for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      }
      predictorsTest[ , ptk] <- xk
    }
  }
  if (length(xfactor) > 0 & nMiss > 0) {
    for (k in 1:length(xfactor)) {
      ptk <- (colnames(imputedData) == xfactor[k])
      factor.k <- xfactor.levels[[k]][imputedData[ , ptk ]]
      labels.k <- xfactor.levels[[k]][sort(unique(imputedData[ , ptk ]))]
      xk.org <- xk <- factor(factor.k, labels = labels.k, levels = labels.k, exclude = NULL)  
      if (length(setdiff(xfactor.levels[[k]], labels.k)) > 0) {
        levels(xk) <- xfactor.levels[[k]]
        for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      }
      imputedData[ , ptk] <- xk
    }
  }
  if (length(xfactor.order) > 0 & nMiss > 0) {
    for (k in 1:length(xfactor.order)) {
      ptk <- (colnames(imputedData) == xfactor.order[k])
      factor.k <- xfactor.order.levels[[k]][imputedData[ , ptk ]]
      labels.k <- xfactor.order.levels[[k]][sort(unique(imputedData[ , ptk ]))]
      xk.org <- xk <- factor(factor.k, labels = labels.k, levels = labels.k, exclude = NULL, ordered = TRUE)
      if (length(setdiff(xfactor.order.levels[[k]], labels.k)) > 0) {
        levels(xk) <- xfactor.order.levels[[k]]
        for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      }
      imputedData[ , ptk] <- xk
    }
  }
  if (n.event > 1) {    
    ens.names <- list(NULL, NULL, c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")))
    err.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), NULL)
    vimp.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), object$predictorNames)
    poe.names <- list(paste("event.", 1:(n.event - 1), sep = ""), NULL)
  }
  else {
    ens.names <- list(NULL, NULL, NULL)
    err.names <- list(NULL, NULL)
    vimp.names <- list(NULL, object$predictorNames)
  }
  if ((predictedOutcomeBits == 0)) {
    mortality <- apply(rbind(array(nativeOutput$fullEnsemble, c(length(Cens), length(object$timeInterest),
      n.event))[,,1]), 1, function(x, wt) {sum(x*wt , na.rm = TRUE)}, wt = Risk)
  }
  else {
    mortality <- array(nativeOutput$fullEnsemble, c(length(Cens), length(object$timeInterest),
      n.event))[,1,1]
  }
  if (outcome == "train") {
    if (performance) {
      ERR <- matrix(nativeOutput$performance, ncol=ntree, byrow=TRUE,
                    dimnames=err.names)[1:n.event, ]
    }
    else {
     ERR <- NULL
   }
  }
  make.vec <- function(x) {
    if (!is.null(dim(x)) && nrow(x) == 1) {
      x.names <- colnames(x)
      x <- c(x)
      names(x) <- x.names
      x
    }
    else {
      x
    }
  }
  if (outcome == "train") {
    if (importance != "none" & performance) {
      VIMP <- make.vec(matrix(nativeOutput$importance, ncol=length(object$predictorNames),
                    byrow=TRUE, dimnames = vimp.names)[1:n.event,, drop = FALSE] -
        matrix(nativeOutput$performance, ncol=ntree, byrow=TRUE,
                    dimnames = err.names)[1:n.event, ntree])
    }
    else {
      VIMP <- NULL
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
    time =  (if (performance) Time else NULL),
    cens =  (if (performance) Cens else NULL),
    predictorNames = object$predictorNames,
    predictors = predictorsTest,
    ensemble = (if (predictedOutcomeBits == 0) 
      array(nativeOutput$fullEnsemble, c(length(Cens), length(object$timeInterest),
      n.event), dimnames = ens.names)[,,1:n.event] else NULL),
    poe = (if (!is.null(nativeOutput$fullPOE)) matrix(nativeOutput$fullPOE,
      ncol=length(Cens), byrow=TRUE, dimnames=poe.names) else NULL),
    mortality = mortality,
    err.rate = ERR,
    importance = VIMP,
    proximity = nativeOutput$proximity,
    imputedIndv = (if (nMiss>0) imputedIndv else NULL),
    imputedData = (if (nMiss>0) imputedData else NULL),
    splitDepth  = (if (split.depth) matrix(nativeOutput$splitDepth,
                    nrow = length(Cens), byrow=FALSE) else NULL)
  )
  class(rsfOutput) <- c("rsf", "predict")
  return(rsfOutput)
}
