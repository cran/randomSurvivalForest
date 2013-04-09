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

rsf <- function(
    formula,#
    data = NULL,#
    ntree = 1000,#
    mtry = NULL,#
    nodesize = NULL,#
    splitrule = NULL,#
    nsplit = 0,#
    importance = c("randomsplit", "permute", "none")[1],#
    big.data = FALSE,#
    na.action = c("na.omit", "na.impute")[1],#
    nimpute = 1,#
    predictorWt = NULL,#
    forest = TRUE,#
    proximity = FALSE,#
    varUsed = NULL,#
    split.depth = FALSE,#
    seed = NULL,#
    do.trace = FALSE,#
    ...)
{
  warning("\nThis package has been deprecated. Please upgrade to the new package 'randomForestSRC'.\n")
  if (!inherits(formula, "formula")) stop("'formula' is not a formula object.")
  if (is.null(data)) stop("'data' is missing.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE 
  mf[[1]] <- as.name("model.frame")
  fNames <- all.vars(formula, max.names=1e7)
  if (!( (all.names(formula)[2] == "Survrsf" | all.names(formula)[2] == "Surv")
       & sum(is.element(names(data), fNames[1:2]))==2)) 
    stop("Formula is specified incorrectly.")
  if (sum(inherits(data, c("data.frame", "impute.only"), TRUE) == c(1, 2)) == 2) {
    impute.only <- TRUE
  }
  else {
    impute.only <- FALSE
  }
  if (importance == TRUE)  importance <- "permute"
  if (importance == FALSE) importance <- "none"    
  pt.Time <- is.na(data[,is.element(names(data), fNames[1])])
  pt.Cens <- is.na(data[,is.element(names(data), fNames[2])])
  if (all(pt.Time) | all(pt.Cens))
    stop("All records missing on survival and (or) censoring.  Analysis not meaningful.")
  pt.col <- apply(data, 2, function(x){all(is.na(x))})
  data <- data[, !pt.col]    
  if (nrow(data) <= 2)
    stop("No predictors left in the NA-processed data.  Analysis not meaningful.")
  pt.row <- apply(data, 1, function(x){all(is.na(x))})
  data <- data[!pt.row, ]    
  if (ncol(data) <= 1)
    stop("Less than 2 records in the NA-processed data.  Analysis not meaningful.")
  if (fNames[3] != ".") {
    if (!big.data) {
      predTempNames <- attr(terms(formula), "term.labels")
    }
    else {
      predTempNames <- fNames[-c(1:2)]
    }
  }
  else {
    predTempNames <- names(data)[!is.element(names(data), fNames[1:2])]
    formula <- as.formula(paste(paste("Surv(",fNames[1],",",fNames[2],") ~"),
                                paste(predTempNames, collapse="+")))
  }
  data <- rm.na.levels(data,
        prednames=(if (fNames[3] != ".") fNames[-c(1,2)] else predTempNames))
  get.factor <- extract.factor(data,
        prednames=(if (fNames[3] != ".") fNames[-c(1,2)] else predTempNames))
  xfactor <- get.factor$xfactor
  xfactor.levels <- get.factor$xfactor.levels
  xfactor.order <- get.factor$xfactor.order
  xfactor.order.levels <- get.factor$xfactor.order.levels
  if (na.action != "na.omit" & na.action != "na.impute") na.action <- "na.omit"
  if (!big.data) {
    data <- as.data.frame(data.matrix(data))
    old.na.action <- options()$na.action
    na.keep <- function(x){ x }
    na.takeOut <- function(x){ na.omit(x) }
    if (na.action == "na.omit") options(na.action=na.takeOut) else options(na.action=na.keep)
    data <- as.data.frame(model.matrix(as.formula(paste("~ -1 +",
                 paste(c(fNames[1:2], predTempNames), collapse="+"))), data))
    options(na.action=old.na.action)
  }
  else { 
    data <- data[ , is.element(names(data), c(fNames[1:2], predTempNames))]
    data <- as.data.frame(data.matrix(data))
    if (na.action == "na.omit") {
      data <- na.omit(data)
    }
  }
  predictorNames <- names(data)[!is.element(names(data), fNames[1:2])]
  predictors <- as.matrix(data[,is.element(names(data), predictorNames)])
  rownames(predictors) <- colnames(predictors) <- NULL
  n <- nrow(predictors)
  ncov <- ncol(predictors)
  if (n <= 1)
    stop("Less than two records in the NA-processed data.  Analysis not meaningful.")
  predictorType <- rep("R", ncov)
  if (length(xfactor) > 0) {
    predictorType[is.element(predictorNames, xfactor)] <- "C"
  }
  if (length(xfactor.order) > 0) {
    predictorType[is.element(predictorNames, xfactor.order)] <- "I"
  }
  nsplit <- round(nsplit)
  if (nsplit < 0) stop("Invalid nsplit value specified.")    
  if (is.null(predictorWt)) {
    predictorWt <- rep(1/ncov, ncov)
  }
  else {
    if (any(predictorWt < 0) | length(predictorWt) != ncov | all(predictorWt == 0)) {
      predictorWt <- rep(1/ncov, ncov)
    }
    else {
      predictorWt <-predictorWt/sum(predictorWt)
    }
  }
  ntree <- round(ntree)
  if (ntree < 1) stop("Invalid choice of 'ntree'.  Cannot be less than 1.")
  nimpute <- round(nimpute)
  if (nimpute < 1) stop("Invalid choice of 'nimpute'.  Cannot be less than 1.")
  if (!is.null(mtry)) {
    mtry <- round(mtry)
    if (mtry < 1 | mtry > ncov) mtry <- max(1, min(mtry, ncov))
  }
  else {
    mtry <- max(floor(sqrt(ncov)), 1)
  }
  if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
  seed <- -round(abs(seed))
  Time <- data[,is.element(names(data), fNames[1])]
  Cens <- data[,is.element(names(data), fNames[2])]
  nMiss <- sum(is.na(Cens) | is.na(Time) | apply(predictors, 1, function(x){any(is.na(x))}))
  if (all(na.omit(Cens) == 0)) {
    stop("No deaths in data!")
  }
  eventType <- unique(na.omit(Cens))
  if (sum(eventType >= 0) != length(eventType)) {
    stop("Censoring variable must be coded as NA, 0, or greater than 0.")    
  }
  eventType <- eventType[eventType > 0]
  splitrule.names <- c("logrank", "conserve", "logrankscore", "random", "logrankCR")
  if (is.null(splitrule)) {
    if (length(eventType) ==  1) {
      splitrule.idx <- which(splitrule.names == "logrank")
    }
    else {
      splitrule.idx <- which(splitrule.names == "logrankCR")
    }
    splitrule <- splitrule.names[splitrule.idx]
  }
  else {
    splitrule.idx <- which(splitrule.names == splitrule)
    if (length(splitrule.idx) != 1) {
      stop("Invalid split rule specified:  ", splitrule)
    if (length(eventType) ==  1 & splitrule.idx == 5)
      stop("Cannot specify logrankCR splitting for right-censored data")
    }
    if (splitrule == "random" & nsplit == 0) {
      nsplit <- 1
    }
  }
  splitrule.nsplit <- nsplit
  if (!all(na.omit(Time) > 0)) {
    stop("Time must be strictly positive.")
  }
  remove(data)
  nonMissingOutcome <- which(!is.na(Cens) & !is.na(Time))
  nonMissingDeathFlag <- (Cens[nonMissingOutcome] != 0)    
  timeInterest <- sort(unique(Time[nonMissingOutcome[nonMissingDeathFlag]]))
  N <- length(timeInterest)
  if (N <= 1) stop("Less than 2 unique event times.  Analysis not meaningful.")
  if (!is.null(nodesize)) {
    nodesize <- min(round(nodesize), N) 
    if (nodesize < 1) stop("Invalid choice of 'nodesize'. Cannot be less than 1.")
  }
  else {
    nodesize <- min(c(round(0.632*sum(na.omit(Cens) != 0)), c(3, 6)[1+1*(length(eventType) > 1)], N))
  }
  if (impute.only) {
    if (nMiss == 0) stop("data has no missing values, using 'impute' makes no sense")
  }
  else {
    Risk <- apply(cbind(1:length(timeInterest)),
                1,
                function(i, tau, tunq) {sum(tau >= tunq[i])},
                tau = Time[nonMissingOutcome], tunq = timeInterest)
    Risk <- Risk - c(Risk[-1] , 0)
    impute.only <- FALSE
  }   
  if (!is.logical(do.trace)) {
    if (do.trace >= 1) {
      do.trace <- 2^24 * round(do.trace) + 1
    }
    else {
      do.trace <- 0
    }
  }
  else {
    do.trace <- 1 * do.trace
  }
  if (importance == "none") {
    importanceBits <- 0
  }
  else if (importance == "permute") {
    importanceBits <- 2^11 + 0
  }
  else if (importance == "randomsplit") {
    importanceBits <- 2^11 + 2^9
  }
  else {
    stop("Invalid choice for 'importance' option:  ", importance)
  }
  if (is.null(varUsed)) {
    varUsedBits <- 0
  }
  else if (varUsed == "all.trees") {
    varUsedBits <- 2^13 + 0
  }
  else if (varUsed == "by.tree") {
    varUsedBits <- 2^13 + 2^12
  }
  else {
    stop("Invalid choice for 'varUsed' option:  ", varUsed)
  }
  if (impute.only == TRUE) {
    impute.only <- 2^16
  }
  else {
    impute.only <- 0
  }
  nativeOutput <- .Call("rsfGrow",
                        as.integer(do.trace),
                        as.integer(impute.only +
                                   varUsedBits +
                                   importanceBits +
                                   (32 * (if (forest) 1 else 0)) +
                                   (8  * (if (proximity) 1 else 0)) +
                                   (2^15 * (if (split.depth) 1 else 0))),
                        as.integer(seed),
                        as.integer(splitrule.idx), 
                        as.integer(splitrule.nsplit),
                        as.integer(mtry),
                        as.integer(ntree),
                        as.integer(nodesize),
                        as.integer(n),
                        as.double(Time),
                        as.double(Cens),
                        as.integer(ncov),
                        as.numeric(predictors),
                        as.integer(length(timeInterest)),
                        as.double(timeInterest),
                        as.double(predictorWt),
                        as.character(predictorType),
                        as.integer(nimpute))
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }
  if (nMiss > 0) {
    imputedData <- matrix(nativeOutput$imputation, nrow = nMiss, byrow = FALSE)
    imputedIndv <- imputedData[,1]
    imputedData <- as.matrix(imputedData[,-1])
    if (nMiss == 1) imputedData <- t(imputedData)
    if (nimpute > 1) {
      if (nimpute == 2) {
        imputedOOBData <- matrix(nativeOutput$oobImputation, nrow = nMiss, byrow = FALSE)
        imputedOOBData <- as.matrix(imputedOOBData[,-1])
        if (nMiss == 1) imputedOOBData <- t(imputedOOBData)
        Cens[imputedIndv] <- imputedOOBData[, 1]
        Time[imputedIndv] <- imputedOOBData[, 2]
        predictors[imputedIndv, ] <- imputedOOBData[, -1:-2]
      }
      else {
        Cens[imputedIndv] <- imputedData[, 1]
        Time[imputedIndv] <- imputedData[, 2]
        predictors[imputedIndv, ] <- imputedData[, -1:-2]
      }
      imputedIndv    <- NULL
      imputedData    <- NULL
      imputedOOBData <- NULL
    }  
    else {
      colnames(imputedData) <- c(fNames[2], fNames[1], predictorNames)
      imputedData=as.data.frame(imputedData)
    }
  }
  predictors <- as.data.frame(predictors)
  colnames(predictors) <- predictorNames
  if (length(xfactor) > 0) {
    for (k in 1:length(xfactor)) {      
      ptk <- (colnames(predictors) == xfactor[k])
      factor.k <- xfactor.levels[[k]][predictors[ , ptk ]]
      labels.k <- xfactor.levels[[k]][sort(unique(predictors[ , ptk ]))]
      xk.org <- xk <- factor(factor.k, labels = labels.k, levels = labels.k, exclude = NULL)  
      if (length(setdiff(xfactor.levels[[k]], labels.k)) > 0) {
        levels(xk) <- xfactor.levels[[k]]
        for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      }
      predictors[ , ptk] <- xk
    }
  }
  if (length(xfactor.order) > 0) {
    for (k in 1:length(xfactor.order)) {
      ptk <- (colnames(predictors) == xfactor.order[k])
      factor.k <- xfactor.order.levels[[k]][predictors[ , ptk ]]
      labels.k <- xfactor.order.levels[[k]][sort(unique(predictors[ , ptk ]))]
      xk.org <- xk <- factor(factor.k, labels = labels.k, levels = labels.k, exclude = NULL, ordered = TRUE)
      if (length(setdiff(xfactor.order.levels[[k]], labels.k)) > 0) {
        levels(xk) <- xfactor.order.levels[[k]]
        for (l in 1:length(xk)) { xk[l]  <- xk.org[l] }
      }
      predictors[ , ptk] <- xk
    }
  }
  if (length(xfactor) > 0 & nMiss > 0 & nimpute < 2) {
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
  if (length(xfactor.order) > 0 & nMiss > 0 & nimpute < 2) {
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
  if (!is.null(varUsed)) {
    if (varUsed == "all.trees") {
      varUsed <- nativeOutput$varUsed
      names(varUsed) <- predictorNames
    }
    else {
      varUsed <- matrix(nativeOutput$varUsed, nrow = ntree, byrow = TRUE)
      colnames(varUsed) <- predictorNames
    }
  }
  if (forest) {
    nativeArray <- as.data.frame(cbind(nativeOutput$treeID,
                                       nativeOutput$nodeID,
                                       nativeOutput$parmID,
                                       nativeOutput$contPT,
                                       nativeOutput$mwcpSZ))
    names(nativeArray) <- c("treeID", "nodeID", "parmID", "contPT", "mwcpSZ")
    nativeFactorArray <- nativeOutput$mwcpPT
    forest <- list(nativeArray = nativeArray,
                   nativeFactorArray = nativeFactorArray,
                   timeInterest = timeInterest, 
                   predictorNames = predictorNames,
                   seed = nativeOutput$seed,
                   predictors = predictors,
                   formula = formula,
                   time = Time,
                   cens = Cens)
    class(forest) <- c("rsf", "forest")
  }
  else {
    forest <- NULL
  }
  if ((length(eventType) > 1) & (splitrule.idx >= 4)) {
    n.event <- length(eventType) + 1
    ens.names <- list(NULL, NULL, c("CHF", paste("condCHF.", 1:(n.event-1), sep = "")))
    err.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), NULL)
    vimp.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), predictorNames)
    poe.names <- list(paste("event.", 1:(n.event - 1), sep = ""), NULL)
  }
  else {
    Cens[Cens != 0 & !is.na(Cens)] <-  min(Cens[Cens != 0 & !is.na(Cens)])
    n.event <- 1
    ens.names <- list(NULL, NULL, NULL)
    vimp.names <- list(NULL, predictorNames)
    err.names <- list(NULL, NULL)
  }
  if (!impute.only) {
    mortality <- apply(array(nativeOutput$fullEnsemble, c(n, length(timeInterest), n.event))[,,1],
                     1, function(x, wt) {sum(x*wt)}, wt = Risk)
    oob.mortality <- apply(array(nativeOutput$oobEnsemble, c(n, length(timeInterest), n.event))[,,1],
                     1,function(x, wt) {sum(x*wt)}, wt = Risk)
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
  if (importance != "none") {
      VIMP <- make.vec(matrix(nativeOutput$importance, ncol = ncov, 
                    byrow = TRUE, dimnames = vimp.names)[1:n.event,, drop = FALSE] -
      matrix(nativeOutput$performance, ncol = ntree, byrow = TRUE,
                    dimnames = err.names)[1:n.event, ntree])
  }
  else {
    VIMP <- NULL
  }
  if (impute.only) {
    nativeOutput$fullEnsemble <- nativeOutput$oobEnsemble <- 
    nativeOutput$fullPOE <- nativeOutput$oobPOE <- mortality <- oob.mortality <-
    nativeOutput$performance <- nativeOutput$leafCount <- VIMP <- forest <-
    nativeOutput$proximity <- varUsed <-  NULL
  }
  rsfOutput <- list(
    call = match.call(),
    formula = formula,
    n = n,
    ndead = sum(na.omit(Cens) != 0),
    ntree = ntree,
    nimpute = nimpute,                      
    mtry = mtry,
    nodesize = nodesize,
    splitrule = splitrule,
    nsplit = nsplit,
    time = Time,
    cens = Cens,
    timeInterest = timeInterest,
    predictorNames = predictorNames,
    predictorWt = predictorWt,
    predictors = predictors,
    ensemble = (if (!is.null(nativeOutput$fullEnsemble))
       array(nativeOutput$fullEnsemble, c(n, length(timeInterest), n.event), dimnames=ens.names)[,,1:n.event]
       else NULL),
    oob.ensemble = (if (!is.null(nativeOutput$oobEnsemble))
        array(nativeOutput$oobEnsemble,  c(n, length(timeInterest), n.event), dimnames=ens.names)[,,1:n.event]
       else NULL),
    poe = (if (!is.null(nativeOutput$fullPOE))
           matrix(nativeOutput$fullPOE, ncol=n, byrow=TRUE, dimnames=poe.names) else NULL),
    oob.poe = (if (!is.null(nativeOutput$oobPOE))
           matrix(nativeOutput$oobPOE, ncol=n, byrow=TRUE, dimnames=poe.names) else NULL),
    mortality = mortality,
    oob.mortality = oob.mortality,
    err.rate = (if (!is.null(nativeOutput$performance))
                matrix(nativeOutput$performance, ncol=ntree, byrow=TRUE, dimnames=err.names)[1:n.event, ]
                else NULL),
    leaf.count = nativeOutput$leafCount,
    importance = VIMP,
    forest = forest,
    proximity = nativeOutput$proximity,
    varUsed = varUsed,
    imputedIndv = (if (nMiss > 0) imputedIndv else NULL),
    imputedData = (if (nMiss > 0) imputedData else NULL),
    splitDepth  = (if (split.depth) matrix(nativeOutput$splitDepth, nrow=n, byrow=FALSE) else NULL)
  )
  if (!big.data) {
    class(rsfOutput) <- c("rsf", "grow")
  }
  else {
    class(rsfOutput) <- c("rsf", "grow", "bigdata")
  }
  return(rsfOutput)
}
