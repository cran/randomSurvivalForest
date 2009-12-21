####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.0
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
# Primary R function for Random Surival Forests
# ----------------------------------------------------------------------
#  Prediction and variable selection for right censored survival and
#  competing risk data using Random Survival Forests (RSF) (Ishwaran,
#  Kogalur, Blackstone and Lauer, 2008).  A random forest (Breiman, 2001)
#  of survival trees is grown and used to estimate an ensemble cumulative
#  hazard function (CHF) and ensemble conditional cumulative hazard function
#  (CCHF) (the latter in the case of competing risks).  Trees can be
#  grown using different survival tree splitting rules.  An
#  out-of-bag estimate of Harrell's concordance index (Harrell,
#  1982) is provided for assessing prediction accuracy.  Variable
#  importance (VIMP) for single, as well as grouped variables, can be
#  used to filter variables and to assess variable predictiveness.
#  Missing data (x-variables, survival times, censoring indicators) can
#  be imputed on both training and test data. Note this is the default
#  generic method for the package.
########################################################################

rsf.default <- function(
    formula,
    data = NULL,
    ntree = 1000,
    mtry = NULL,
    nodesize = NULL,
    splitrule = NULL,
    nsplit = 0,
    importance = c("randomsplit", "permute", "none")[1],
    big.data = FALSE,
    na.action = c("na.omit", "na.impute")[1],
    nimpute = 1,
    predictorWt = NULL,
    forest = FALSE,
    proximity = FALSE,
    varUsed = NULL,
    split.depth = FALSE,
    seed = NULL,
    do.trace = FALSE,
    ...)
{

  ## Preliminary checks for formula and data
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

  ## Ensure backward compatability with 3.0 user scripts
  if (importance == TRUE)  importance <- "permute"
  if (importance == FALSE) importance <- "none"    

  ## Preliminary processing of missing data:
  ## If all survival times or all censoring values missing: stop
  ## Remove all predictors will missing values in all entries
  ## Next remove all records with missing values in all entries
  ## At each step check if enough data remains
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

  ## Get predictor names and formula
  ## If factor has an NA in its levels, remove it
  ## Identify unordered/ordered factors 
  ## Save factor levels for immutable map
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
    formula <- as.formula(paste(paste("Survrsf(",fNames[1],",",fNames[2],") ~"),
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

  ## Data conversion to numeric mode
  ## Get predictors
  ## Be mindful of NA's
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
  else {#big.data
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

  ## Set predictor type
  predictorType <- rep("R", ncov)
  if (length(xfactor) > 0) {
    predictorType[is.element(predictorNames, xfactor)] <- "C"
  }
  if (length(xfactor.order) > 0) {
    predictorType[is.element(predictorNames, xfactor.order)] <- "I"
  }
  ## Set predictor weight
  ## Check that nsplit is set correctly 
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

  ## Final checks on option parameters
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

  ## Time and censoring checks:
  ## Check to see if there are deaths.  
  ## Check that censoring and time are properly coded.
  Time <- data[,is.element(names(data), fNames[1])]
  Cens <- data[,is.element(names(data), fNames[2])]
  nMiss <- sum(is.na(Cens) | is.na(Time) | apply(predictors, 1, function(x){any(is.na(x))}))
  if (all(na.omit(Cens) == 0)) {
    stop("No deaths in data!")
  }

  ## Extract the unique event types.
  eventType <- unique(na.omit(Cens))
  ## Ensure they are all greater than or equal to zero.
  if (sum(eventType >= 0) != length(eventType)) {
    stop("Censoring variable must be coded as NA, 0, or greater than 0.")    
  }
  ## Discard the censored state, if it exists.
  eventType <- eventType[eventType > 0]

  
  ## Ensure coherency of the split rule.  If it has not been specified, and
  ## the data set does not contain competing risks (CR), then the split rule is
  ## logrank.  If it has not been specified, and the data set does contain CR, then
  ## the split rule is logrankCR.  This allows the user to specify a non-CR analysis
  ## even if the data contains CR.
  
  splitrule.names <- c("logrank", "conserve", "logrankscore", "random", "logrankCR")
  if (is.null(splitrule)) {
    ## No split rule specified
    if (length(eventType) ==  1) {
      splitrule.idx <- which(splitrule.names == "logrank")
    }
    else {
      splitrule.idx <- which(splitrule.names == "logrankCR")
    }
    splitrule <- splitrule.names[splitrule.idx]
  }
  else {
    ## User specified split rule
    splitrule.idx <- which(splitrule.names == splitrule)
    if (length(splitrule.idx) != 1) {
      stop("Invalid split rule specified:  ", splitrule)
    if (length(eventType) ==  1 & splitrule.idx == 5)
      stop("Cannot specify logrankCR splitting for right-censored data")
    }
    ## the native code allows nsplit > 0, but is dynamically resized
    if (splitrule == "random" & nsplit == 0) {
      nsplit <- 1
    }
  }

  ## Initialize nsplit, noting that it may be been overridden.
  splitrule.nsplit <- nsplit


  if (!all(na.omit(Time) > 0)) {
    stop("Time must be strictly positive.")
  }

  ## Don't need the data anymore
  remove(data)
  
  ## Set grid of time points
  ## Final checks for event time consistency
  ## Set nodesize
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

  ## Check status of impute mode (internal use only).  Check all is ok.
  ## Otherwise ... work out individuals at risk (used later for mortality calculation)
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

  ## Convert trace into native code parameter
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

  ## Convert importance option into native code parameter
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


  ## Convert varUsed option into native code parameter
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

  ## Convert impute.only option into native code parameter.
  ## It will always be initialized.
  if (impute.only == TRUE) {
    impute.only <- 2^16
  }
  else {
    impute.only <- 0
  }

####################################################################
  ## rsf.default(...)
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
  ##    - note that the native code is capable of some optional
  ##      outputs not implemented on the R-side by (**).
  ##
  ##      GROW PRED INTR IMPUTE
  ##       *         *          0x000001 = OENS
  ##       *    *               0x000002 = FENS
  ##       *    *    *          0x000004 = PERF
  ##       **   **              0x000008 = PROX
  ##       *    *    *    (**)  0x000010 = LEAF
  ##       **                   0x000020 = TREE \ part of 
  ##       **                   0x000040 = SEED / the forest 
  ##       **   **   **    *    0x000080 = MISS
  ##       ***       **    *    0x000100 = OMIS
  ##       **   **   *          0x000200 = \  VIMP_TYPE
  ##       **   **   *          0x000400 =  | VIMP_JOIN
  ##       **   **   *          0x000800 = /  VIMP
  ##       **  (**) (**)        0x001000 = \  VUSE_TYPE
  ##       **  (**) (**)        0x002000 = /  VUSE
  ##      (**) (**)  *          0x004000 = POUT_TYPE (non SEXP)
  ##       **   **   (**)       0x008000 = SPLT_DPTH
  ##       *         *                     OPOE (tracks OENS)
  ##       *    *                          FPOE (tracks FENS)
  ##
  ##       *   default  output
  ##       **  optional output
  ##      (**) optional output not implemented on the R-side
  ##
  ## 03 - random seed for repeatability, integer < 0 
  ## 04 - split rule to be used 
  ##       1 = Log-Rank
  ##       2 = Conservation of Events
  ##       3 = Log-Rank Score
  ##       4 = Random Split
  ##       5 = Log-Rank CR
  ## 05 - number of random split points, non-negative integer
  ## 06 - number of covariates to be randomly selected for
  ##      growing tree, integer > 0
  ## 07 - number of trees in the forest, integer > 0
  ## 08 - minimum number of deaths allowed in a node, integer > 0
  ## 09 - number of observations in data set, integer > 1
  ## 10 - vector of observed event times, doubles > 0
  ## 11 - vector of observed event types
  ##         0  = censored
  ##       > 0  = events
  ## 12 - number of predictors, integer > 0
  ## 13 - [p x n] matrix of predictor observations
  ## 14 - number of time points of interest, integer > 0
  ## 15 - vector of time points of interest, doubles
  ## 16 - vector of random covariate weights, doubles
  ## 17 - vector of predictor types
  ##       "R" - real value
  ##       "I" - integer value 
  ##       "C" - categorical 
  ## 18 - number of impute iterations, integer > 0
  ## 
  ## #################################################################

  ## ###########################################################
  ##  SEXP outputs:  See parameter 2 above and note "*".
  ## ###########################################################    

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

  ## check for error return condition in the native code
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }
    
  ## Check if there was missing data, and assign imputed data if so.
  if (nMiss > 0) {
    imputedData <- matrix(nativeOutput$imputation, nrow = nMiss, byrow = FALSE)
    imputedIndv <- imputedData[,1]
    imputedData <- as.matrix(imputedData[,-1])
    if (nMiss == 1) imputedData <- t(imputedData)
    
    ## Fill NA's in original GROW data with multiply imputed values.
    ## This will now serve as the forest data set and will enable
    ## recovery of terminal node membership with the head of the
    ## seed chain.
    if (nimpute > 1) {
      if (nimpute == 2) {
        ## The forest was grown using overlaid OOB summary values.
        imputedOOBData <- matrix(nativeOutput$oobImputation, nrow = nMiss, byrow = FALSE)
        imputedOOBData <- as.matrix(imputedOOBData[,-1])
        if (nMiss == 1) imputedOOBData <- t(imputedOOBData)
        Cens[imputedIndv] <- imputedOOBData[, 1]
        Time[imputedIndv] <- imputedOOBData[, 2]
        predictors[imputedIndv, ] <- imputedOOBData[, -1:-2]
      }
      else {
        ## The forest was grown using overlaid full (all) summary values.         
        Cens[imputedIndv] <- imputedData[, 1]
        Time[imputedIndv] <- imputedData[, 2]
        predictors[imputedIndv, ] <- imputedData[, -1:-2]
      }
      ## Remove the imputed data outputs.
      imputedIndv    <- NULL
      imputedData    <- NULL
      imputedOOBData <- NULL
      
    }  
    else {
      ## Add column names to the imputed data outputs in the absence
      ## of multiple imputation.
      ## names(Cens) = fNames[2]
      ## names(Time) = fNames[1]    
      colnames(imputedData) <- c(fNames[2], fNames[1], predictorNames)
      imputedData=as.data.frame(imputedData)
    }
  }

  ## Add column names to predictor matrix 
  predictors <- as.data.frame(predictors)
  colnames(predictors) <- predictorNames

  ## Map predictor factors back to original values
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
  
  ## Map imputed data factors back to original values
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

  ## Get varUsed information
  ## Add names for nicety
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

  ## Define the forest.
  if (forest) {
    nativeArray <- as.data.frame(cbind(nativeOutput$treeID,
                                       nativeOutput$nodeID,
                                       nativeOutput$parmID,
                                       nativeOutput$contPT,
                                       nativeOutput$mwcpSZ))
    names(nativeArray) <- c("treeID", "nodeID", "parmID", "contPT", "mwcpSZ")

    ## This can be NULL if there are no factor splits.
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

  # number of event types
  # pretty names
  if ((length(eventType) > 1) & (splitrule.idx >= 4)) {
    ## Split rule can be logrankCR or random
    n.event <- length(eventType) + 1
    ens.names <- list(NULL, NULL, c("CHF", paste("condCHF.", 1:(n.event-1), sep = "")))
    err.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), NULL)
    vimp.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), predictorNames)
    poe.names <- list(paste("event.", 1:(n.event - 1), sep = ""), NULL)
  }
  else {
    # map all non-zero, non-na censoring values to a unique integer value
    Cens[Cens != 0 & !is.na(Cens)] <-  min(Cens[Cens != 0 & !is.na(Cens)])
    n.event <- 1
    ens.names <- list(NULL, NULL, NULL)
    vimp.names <- list(NULL, predictorNames)
    err.names <- list(NULL, NULL)
  }

  ## ensemble mortality (in-bag and OOB)
  if (!impute.only) {
    mortality <- apply(array(nativeOutput$fullEnsemble, c(n, length(timeInterest), n.event))[,,1],
                     1, function(x, wt) {sum(x*wt)}, wt = Risk)
    oob.mortality <- apply(array(nativeOutput$oobEnsemble, c(n, length(timeInterest), n.event))[,,1],
                     1,function(x, wt) {sum(x*wt)}, wt = Risk)
  }
  
  ## VIMP
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

  ## assign NULL values under impute mode (some redundancy)
  if (impute.only) {
    nativeOutput$fullEnsemble <- nativeOutput$oobEnsemble <- 
    nativeOutput$fullPOE <- nativeOutput$oobPOE <- mortality <- oob.mortality <-
    nativeOutput$performance <- nativeOutput$leafCount <- VIMP <- forest <-
    nativeOutput$proximity <- varUsed <-  NULL
  }

  ## make the output object
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

