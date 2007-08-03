##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.0.1
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

########################################################################
# Primary R function for Random Surival Forests
# ----------------------------------------------------------------------
# Random Survival Forests (RSF) for right censored survival data (Ishwaran,
# Kogalur, Blackstone and Lauer, 2007).  RSF is an extension of
# Breiman's Random Forests (Breiman, 2001) to survival analysis
# settings.  Algorithm uses a binary recursive tree growing procedure
# with different splitting rules for growing an ensemble cumulative
# hazard function.  An out-of-bag (OOB) estimate of Harrell's concordance
# index (Harrell, 1982) is provided for assessing prediction.
# Importance values for variables can be computed.
# Prediction on test data is also available.  Missing data (x-variables,
# survival times, censoring indicators) can be imputed on both training
# and test data. Note this is the default generic method for the
# package.
########################################################################

rsf.default <- function(
    formula,
    data = NULL,
    ntree = 1000,
    mtry = NULL,
    nodesize = NULL,
    splitrule = c("logrank", "conserve", "logrankscore", "logrankapprox")[1],
    importance = TRUE,
    big.data = FALSE,
    na.action = c("na.omit", "na.impute")[1],
    predictorWt = NULL,
    forest = FALSE,
    proximity = FALSE,
    seed = NULL,
    ntime = NULL,
    add.noise = FALSE,
    do.trace = FALSE,
    ...)
{
    ## preliminary checks for formula and data
    if (!inherits(formula, "formula")) stop("'formula' is not a formula object.")
    if (is.null(data)) stop("'data' is missing.")
    if (!is.data.frame(data)) stop("'data' must be a data frame.")
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula","data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE 
    mf[[1]] <- as.name("model.frame")
    fNames <- all.vars(formula, max.names=1e7)
    if (!(all.names(formula)[2]=="Survrsf" &
        sum(is.element(names(data), fNames[1:2]))==2)) 
      stop("Outcome is not a random survival object.  Use 'Survrsf' for the formula.")

    ## details for handling NA's
    ## if all survival times or all censoring values missing stop
    ## first remove all predictors will missing values in all entries
    ## next remove all records with missing values in all entries
    ## at each step check if enough data is left.
    pt.Time <- is.na(data[,is.element(names(data), fNames[1])])
    pt.Cens <- is.na(data[,is.element(names(data), fNames[2])])
    if (all(pt.Time) | all(pt.Cens))
      stop("All records missing on survival and (or) censoring.  Analysis not meaningful.")
    pt.col <- apply(data, 2, function(x){all(is.na(x))})
    data <- data[, !pt.col]    
    if (dim(data)[2] <= 2)
      stop("No predictors left in the NA-processed data.  Analysis not meaningful.")
    pt.row <- apply(data, 1, function(x){all(is.na(x))})
    data <- data[!pt.row, ]    
    if (dim(data)[1] <= 1)
      stop("Less than 2 records in the NA-processed data.  Analysis not meaningful.")

    ## Get predictor names and formula 
    ## Special treatment for big.data=T
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
       formula = as.formula(paste(paste("Survrsf(",fNames[1],",",fNames[2],") ~"),
                            paste(predTempNames, collapse="+")))
    }

    ## More NA processing
    ## Special treatment for big.data=T (coerce factors to real)
    ## Removing NA's depletes the data.  Check if enough data is left.
    if (na.action != "na.omit" & na.action != "na.impute") na.action <- "na.omit"
    if (!big.data) {
      old.na.action <- options()$na.action
      na.keep <- function(x){x}
      na.takeOut <- function(x){na.omit(x)}
      if (na.action == "na.omit") options(na.action=na.takeOut) else options(na.action=na.keep)
      data <- as.data.frame(
            model.matrix(as.formula(paste("~ -1 +",
            paste(c(fNames[1:2], predTempNames), collapse="+"))), data))
      options(na.action=old.na.action)
    }
    else {
      if (na.action == "na.omit") {
         data <- na.omit(data[,is.element(names(data), c(fNames[1:2], predTempNames))])
      }
      else {
        data <- data[,is.element(names(data), c(fNames[1:2], predTempNames))]
      }
    }
    predictorNames <- names(data)[!is.element(names(data), fNames[1:2])]
    predictors <- as.matrix(data[,is.element(names(data), predictorNames)])
    n <- dim(predictors)[1]
    if (n <= 1)
      stop("Less than 2 records in the NA-processed data.  Analysis not meaningful.")
    if (big.data) {###coerce factors to NA
      predictors <- matrix(as.real(predictors), nrow = n, byrow = F)
      coerced.pt <- apply(predictors, 2, function(x){all(is.na(x))})
      if (sum(coerced.pt) == dim(predictors)[2])
        stop("All variables are factors.  Not allowed with big.data=T.")
      if (sum(coerced.pt) > 0){
        predictors <- predictors[,!coerced.pt]
        predictorNames <- predictorNames[!coerced.pt]
      }
    }

    ## add noise variable if requested (and adjust formula and predictorWt)
    if (add.noise || (length(predictorNames) == 1 && length(unique(predictors)) <= 10)) {
      noise <- rnorm(dim(predictors)[1])
      predictors <- cbind(predictors, noise)
      if (any(predictorNames == "noise")) predictorNames[predictorNames == "noise"] <- "NOISE"
      predictorNames <- c(predictorNames, "noise")
      formula <- as.formula(paste(paste("Survrsf(",fNames[1],",",fNames[2],") ~"),
                            paste(c(predTempNames, "noise"), collapse="+")))
      if (!is.null(predictorWt)) 
        predictorWt <- c(predictorWt, min(predictorWt))
    }
    rownames(predictors) <- colnames(predictors) <- NULL
    ncov <- length(predictorNames)
    Time <- data[,is.element(names(data), fNames[1])]
    Cens <- data[,is.element(names(data), fNames[2])]


    ## determine how many records having missing data
    nMiss <- sum(is.na(Cens) | is.na(Time) | apply(predictors, 1, function(x){any(is.na(x))}))
    
    ## set predictor types
    if (!big.data) {
      whole.number <- function(x) {
        n.whole.number <- max(10, round(0.25*length(x)))
        if (length(unique(x)) <= n.whole.number & all((na.omit(x) - floor(na.omit(x))) == 0)) "I" else "R"
      }
      predictorType <- apply(predictors, 2, whole.number)
    }
    else {
      predictorType <- rep("R", ncov)
    }

    ## check to see if there are deaths
    ## check that censoring and time are properly coded
    if (all(na.omit(Cens) == 0)) {
      stop("No deaths in data: consider converting 0's to 1's.")
    }
    if (!all(is.element(unique(na.omit(Cens)), c(0,1)))) {
      stop("Censoring variable must be coded as NA, 0 or 1.")
    }
    if (!all(na.omit(Time) > 0)) {
      stop("Time must be strictly positive.")
    }
  
    ## check that option parameters are correctly specified
    ntree <- round(ntree)
    if (ntree < 1) stop("Invalid choice of 'ntree'.  Cannot be less than 1.")
    if (!is.null(mtry)) {
      mtry <- round(mtry)
      if (mtry < 1 | mtry > ncov) mtry <- max(1, min(mtry, ncov))
    }
    else {
      mtry <- max(floor(sqrt(ncov)), 1)
    }
    splitrule.names = c("logrank", "conserve", "logrankscore", "logrankapprox")
    splitrule.idx = which(splitrule.names == splitrule)
    if (length(splitrule.idx) != 1)
      stop("Invalid split rule specified:  ", splitrule)

    ## seed details
    ## generate seed
    if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
    seed <- -round(abs(seed))

    ## cap the number of unique points at a maximum of ntime
    ## set the nodesize
    nonMissingOutcome <- which(!is.na(Cens) & !is.na(Time))
    nonMissingDeathFlag <- (Cens[nonMissingOutcome] == 1)    
    timeInterest <- sort(unique(Time[nonMissingOutcome[nonMissingDeathFlag]]))
    N <- length(timeInterest)
    if (N<=1) stop("Less than 2 unique event times.  Analysis not meaningful.")
    if (!is.null(nodesize)) {
      nodesize <- min(round(nodesize), N) 
      if (nodesize < 1) stop("Invalid choice of 'nodesize'. Cannot be less than 1.")
    }
    else {
        nodesize <- min(c(round(0.632*sum(na.omit(Cens) == 1)), 3, N))
    }
    if (is.null(ntime) || ntime <= 1 || ntime > N) ntime <- N
    ntime <- as.integer(ntime)
    if (ntime > 5 & ntime < N)
        timeInterest <- timeInterest[unique(as.integer(seq(1, N, length = ntime)))]
    N <- length(timeInterest)
    
    ## work out individuals at risk (used later for mortality calculation)
    Risk <- apply(cbind(1:length(timeInterest)),
                  1,
                  function(i, tau, tunq) {sum(tau >= tunq[i])},
                  tau = Time[nonMissingOutcome], tunq = timeInterest)
    Risk <- Risk - c(Risk[-1],0)

    ## predictorWt details    
    if (is.null(predictorWt)) {
      predictorWt <- rep(1/ncov, ncov)
    }
    else {
      if (any(predictorWt < 0) | length(predictorWt) != ncov | all(predictorWt == 0)) {
        predictorWt <- rep(1/ncov, ncov)
      }
      else {
        predictorWt[predictorWt == 0] <- min(1e-5, 1/ncov^2)
        predictorWt <-predictorWt/sum(predictorWt)
      }
    }
    ## do.trace details
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
    ## Parameters passed the C function rsf(...) are as follows:
    ## #################################################################
    ## 00 - C function name
    ## 01 - trace output flag
    ##    -  0 = no trace output
    ##    - !0 = various levels of trace output
    ## 02 - memory useage protocol (MUP) for return objects
    ##    - any combination of the following are allowed
    ##    - 0x00 = only the default objects are returned
    ##    - 0x01 = return proximity information
    ##    - 0x02 = return resulting forest
    ##    - 0x04 = N/A
    ##    - 0x08 = return variable importance
    ## 03 - random seed for repeatability
    ##    - integer < 0 
    ## 04 - split rule to be used 
    ##    - 1 = Log-Rank
    ##    - 2 = Conservation of Events
    ##    - 3 = Log-Rank Score
    ##    - 4 = Approximate Log-Rank
    ## 05 - number of covariates to be randomly selected for
    ##      growing tree
    ##    - integer > 0
    ## 06 - number of bootstrap iterations
    ##    - integer > 0
    ## 07 - minimum number of deaths allowed in a node
    ##    - integer > 0
    ## 08 - number of observations in data set
    ##    - integer > 1
    ## 09 - vector of observed times of death
    ##    - vector of double values
    ## 10 - vector of observed event types
    ##    - 0 = censored
    ##    - 1 = death
    ## 11 - number of predictors 
    ##    - integer > 0
    ## 12 - [p x n] matrix of predictor observations
    ## 13 - number of time points of interest
    ##    - integer > 0
    ## 14 - vector of time points of interest
    ##    - vector of double values
    ## 15 - random covariate weight
    ##    - vector of double values
    ## 16 - vector of predictor types
    ##      "R" - real value
    ##      "I" - integer value 
    ##      "C" - categorical 
    ## 
    ## #################################################################

    ##############################################################
    ## SEXP outputs (see native code for description):
    ## Note that outputs depend on the MUP flags.
    ##
    ## fullEnsemble - default output
    ## oobEnsemble  - default output
    ## performance  - default output
    ## leafCount    - default output
    ## proximity    - optional
    ## importance   - optional
    ## treeID       \\
    ## nodeID        \\
    ## parmID         || - forest optional
    ## spltPT        ##     in PRED mode
    ## seed         ##
    ## imputedData  - output is dependent on missing data
    ##############################################################    

    nativeOutput <- .Call("rsfGrow",
        as.integer(do.trace),
        as.integer((8 * (if (importance) 1 else 0)) +
                   (2 * (if (forest) 1 else 0)) +
                   (if (proximity) 1 else 0)),
        as.integer(seed),
        as.integer(splitrule.idx), 
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
        as.character(predictorType))

    ## check for error return condition in the native code
    if(is.null(nativeOutput)) {
      stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
    }
    
    ## ensemble mortality (in-bag and OOB) 
    mortality <- apply(matrix(nativeOutput$fullEnsemble, 
                       nrow = n,
                       byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt)},
                       wt = Risk)
    oob.mortality <- apply(matrix(nativeOutput$oobEnsemble, 
                       nrow = n,
                       byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt)},
                       wt = Risk)
    
    ## check if there was missing data
    ## assign imputed data 
    if (nMiss > 0) {
      imputedData <- matrix(nativeOutput$imputation, nrow = nMiss, byrow = FALSE)
      imputedIndv <- imputedData[,1]
      imputedData <- as.matrix(imputedData[,-1])
      if (nMiss == 1) imputedData <- t(imputedData)
    }

    ## define the forest
    if (forest) {
      nativeArray <- as.data.frame(cbind(nativeOutput$treeID,
                                        nativeOutput$nodeID,
                                        nativeOutput$parmID,
                                        nativeOutput$spltPT))
      names(nativeArray) <- c("treeID", "nodeID", "parmID", "spltPT")
      forest <- list(nativeArray = nativeArray, 
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

    ## add column names to predictor matrix
    ## add column names to imputed data
    colnames(predictors) <- predictorNames
    if (nMiss > 0) {
      colnames(imputedData) <- c(fNames[2], fNames[1], predictorNames)
    }
    
    ## create the output object
    rsfOutput <- list(
        call = match.call(),
        formula = formula,
        n = n,
        ndead = sum(na.omit(Cens) == 1),
        ntree = ntree,
        mtry = mtry,
        nodesize = nodesize,
        splitrule = splitrule,
        time = Time,
        cens = Cens,
        timeInterest = timeInterest,
        predictorNames = predictorNames,
        predictorWt = predictorWt,
        predictors = predictors,
        ensemble = matrix(nativeOutput$fullEnsemble, nrow = n, byrow = FALSE),
        oob.ensemble = matrix(nativeOutput$oobEnsemble, nrow = n, byrow = FALSE),
        mortality = mortality,
        oob.mortality = oob.mortality,
        err.rate = nativeOutput$performance,
        leaf.count = nativeOutput$leafCount,
        importance = (if (importance) nativeOutput$importance-nativeOutput$performance[ntree]
                      else NULL),
        forest = forest,
        proximity = (if (proximity) nativeOutput$proximity else NULL),
        imputedIndv = (if (nMiss>0) imputedIndv else NULL),
        imputedData = (if (nMiss>0) imputedData else NULL)
    )

    if (!big.data) {
      class(rsfOutput) <- c("rsf", "grow")
    }
    else {
      class(rsfOutput) <- c("rsf", "grow", "bigdata")
    }

    return(rsfOutput)
    
  }
