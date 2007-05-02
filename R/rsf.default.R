##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 2.1.0
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

#################################################################
# Primary R function for Random Survival Forests
# ---------------------------------------------------------------
# Description:
#  Ishwaran and Kogalur's Random Survival Forests algorithm for right
#  censored survival data (Ishwaran and Kogalur, 2006).  This is a direct
#  extension of Breiman's Random Forests method (Breiman, 2001) to
#  survival analysis settings.  Algorithm uses a binary recursive tree
#  growing procedure with different splitting rules for growing an
#  ensemble cumulative hazard function.  An out-of-bag (OOB) estimate of
#  Harrell's concordance index (Harrell, 1982) is provided for assessing
#  prediction.  Importance values for predictors can also be computed.
#  Prediction on test data is also available.  Note that this is the
#  default generic method for the package.
#################################################################

rsf.default <- function(
    formula,
    data = NULL,
    ntree = 1000,
    mtry = NULL,
    nodesize = NULL,
    splitrule = c("logrank", "conserve", "logrankscore", "logrankapprox")[1],
    importance = TRUE,
    big.data = FALSE,
    predictorWt = NULL,
    forest = FALSE,
    do.trace = FALSE,
    proximity = FALSE,
    seed = NULL,
    ntime = NULL,
    add.noise = FALSE,
    ...)
{
    ### preliminary checks for formula and data
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

    ### process data for native code
    ### remove any NA's and issue warning
    ### special treatement for big.data=T
    ### add noise variable if requested (and adjust formula and predictorWt)
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
    if (!big.data) {
      data <- as.data.frame(
            model.matrix(as.formula(paste("~ -1 +",
            paste(c(fNames[1:2], predTempNames), collapse="+"))), data))
    }
    else {
       data <- na.omit(data[,is.element(names(data), c(fNames[1:2], predTempNames))])
    } 
    predictorNames <- names(data)[!is.element(names(data), fNames[1:2])]
    predictors <- as.matrix(data[,is.element(names(data), predictorNames)])
    if (dim(predictors)[1] <= 1)
      stop("Less than one observation in training data.  Analysis is not meaningful.")
    if (add.noise || (length(predictorNames) == 1 && length(unique(predictors)) <= 10)) {
      noise <- rnorm(dim(predictors)[1])
      predictors <- cbind(predictors, noise)
      predictorNames <- c(predictorNames, "noise")
      formula <- as.formula(paste(paste("Survrsf(",fNames[1],",",fNames[2],") ~"),
                            paste(c(predTempNames, "noise"), collapse="+")))
      if (!is.null(predictorWt)) 
        predictorWt <- c(predictorWt, min(predictorWt))
    }
    rownames(predictors) <- colnames(predictors) <- NULL
    ncov <- length(predictorNames)
    Cens <- data[,is.element(names(data), fNames[2])]
    Time <- data[,is.element(names(data), fNames[1])]

    ### internal checks
    if (min(Cens) < 0 | max(Cens) > 1) 
        stop("Censoring variable must be coded as 0 [censored] and 1 [death].")
    if (sum(Cens == 0) == length(Cens)) 
        stop("No deaths in data: convert 0's to 1's and try again.")

  
    ### prepare variables for native code
    ### check that option parameters are correctly specified
    ### cap the number of unique points at a maximum of ntime
    ### work out individuals at risk (used later for mortality calculation)
    ### do.trace details
    ntree <- round(ntree)
    if (ntree < 1) stop("Invalid choice of 'ntree'.  Cannot be less than 1.")
    if (!is.null(mtry)) {
      mtry <- round(mtry)
      if (mtry < 1 || mtry > ncov) mtry <- max(1, min(mtry, ncov))
    }
    else {
      mtry <- max(floor(sqrt(ncov)), 1)
    }
    if (!is.null(nodesize)) {
      nodesize <- round(nodesize)
      if (nodesize < 1) stop("Invalid choice of 'nodesize'. Cannot be less than 1.")
    }
    else {
      nodesize <- min(sum(Cens == 1), 3)
    }
    splitrule.names = c("logrank", "conserve", "logrankscore", "logrankapprox")
    splitrule.idx = which(splitrule.names == splitrule)
    if (length(splitrule.idx) != 1)
      stop("Invalid split rule specified:  ", splitrule)
    if (is.null(seed)) seed <- rnorm(1)*10000
    if (abs(seed) < 1) seed <- -1 else seed <- as.integer(-1 * abs(seed))
    timeInterest <- sort(unique(Time[Cens == 1]))
    N <- length(timeInterest)
    if (is.null(ntime) || ntime <= 1 || ntime > N) ntime <- N
    ntime <- as.integer(ntime)
    if (ntime > 5 & ntime < N)
        timeInterest <- timeInterest[unique(as.integer(seq(1, N, length = ntime)))]
    N <- length(timeInterest)
    Risk <- apply(cbind(1:length(timeInterest)),
                  1,
                  function(i, tau, tunq) {sum(tau >= tunq[i])},
                  tau = Time, tunq = timeInterest)
    Risk <- Risk - c(Risk[-1],0)

    if (is.null(predictorWt)) {
      predictorWt <- rep(1.0/dim(predictors)[2], dim(predictors)[2])
    }
    else {
      if (any(predictorWt <= 0) | length(predictorWt) != dim(predictors)[2]) {
        predictorWt <- rep(1.0/dim(predictors)[2], dim(predictors)[2])
      }
      else {
        predictorWt <-predictorWt/sum(predictorWt)
      }
    }
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
      
    #############################################################
    # Parameters passed the C function rsf(...) are as follows:
    #############################################################
    # 00 - C function name
    # 01 - trace output flag
    #    -  0 = no trace output
    #    - !0 = various levels of trace output
    # 02 - memory useage protocol (MUP) for return objects
    #    - any combination of the following are allowed
    #    - 0x00 = only the default objects are returned
    #    - 0x01 = return proximity information
    #    - 0x02 = return resulting forest
    #    - 0x04 = N/A
    #    - 0x08 = return variable importance
    # 03 - random seed for repeatability
    #    - integer < 0 
    # 04 - split rule to be used 
    #    - 1 = Log-Rank
    #    - 2 = Conservation of Events
    #    - 3 = Log-Rank Score
    #    - 4 = Approximate Log-Rank
    # 05 - number of covariates to be randomly selected for
    #      growing tree
    #    - integer > 0
    # 06 - number of bootstrap iterations
    #    - integer > 0
    # 07 - minimum number of deaths allowed in a node
    #    - integer > 0
    # 08 - number of observations in data set
    #    - integer > 1
    # 09 - vector of observed times of death
    #    - vector of double values
    # 10 - vector of observed event types
    #    - 0 = censored
    #    - 1 = death
    # 11 - number of predictors 
    #    - integer > 0
    # 12 - [n x p] matrix of predictor observations
    # 13 - number of time points of interest
    #    - integer > 0
    # 14 - vector of time points of interest
    #    - vector of double values
    # 15 - random covariate weight
    #    - vector of double values
    #############################################################

    #############################################################
    # SEXP outputs (see native code for description):
    # Note that outputs depend on the MUP flags.
    #
    # fullEnsemble - default output
    # oobEnsemble  - default output
    # performance  - default output
    # leafCount    - default output
    # proximity    - optional
    # importance   - optional
    # treeID       \\
    # nodeID        \\
    # parmID         || - forest optional
    # spltPT        ##
    # seed         ##
    #############################################################    

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
        as.integer(dim(predictors)[1]),
        as.double(Time),
        as.integer(Cens),
        as.integer(dim(predictors)[2]),
        as.numeric(predictors),
        as.integer(length(timeInterest)),
        as.double(timeInterest),
        as.double(predictorWt)
    )
    mortality <- apply(matrix(nativeOutput$fullEnsemble, 
                       nrow = length(Cens),
                       byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt)},
                       wt = Risk)
    oob.mortality <- apply(matrix(nativeOutput$oobEnsemble, 
                       nrow = length(Cens),
                       byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt)},
                       wt = Risk)
    
    if (forest) {
      nativeArray <- as.data.frame(cbind(nativeOutput$treeID,
                                        nativeOutput$nodeID,
                                        nativeOutput$parmID,
                                        nativeOutput$spltPT))
      names(nativeArray) <- c("treeID", "nodeID", "parmID", "spltPT")
      forest <- list(nativeArray = nativeArray, 
                     timeInterest = timeInterest, 
                     predictorNames = predictorNames,
                     bootstrapSeed = nativeOutput$seed,
                     predictors = predictors,
                     formula = formula,
                     Time = Time,
                     Cens = Cens)
       class(forest) <- c("rsf", "forest")
    }
    else
      forest <- NULL
  
    
    rsfOutput <- list(
        call = match.call(),
        formula = formula,
        n = length(Cens),
        ndead = sum(Cens == 1),
        ntree = ntree,
        mtry = mtry,
        nodesize = nodesize,
        splitrule = splitrule,
        Time = Time,
        Cens = Cens,
        timeInterest = timeInterest,
        predictorNames = predictorNames,
        predictorWt = predictorWt,
        predictors = predictors,
        ensemble = matrix(nativeOutput$fullEnsemble, nrow = length(Cens), byrow = FALSE),
        oob.ensemble = matrix(nativeOutput$oobEnsemble, nrow = length(Cens), byrow = FALSE),
        mortality = (if (max(mortality) <= length(Cens)) mortality else
                     round(mortality*(length(Cens))/(1*(max(mortality) == 0)+
                     max(mortality)))),
        oob.mortality = (if (max(oob.mortality,na.rm=T) <= length(Cens)) oob.mortality else
                        round(oob.mortality*(length(Cens))/(1*(max(oob.mortality,na.rm=T) == 0)+
                        max(oob.mortality,na.rm=T)))),
        err.rate = nativeOutput$performance,
        leaf.count = nativeOutput$leafCount,
        importance = (if (importance) nativeOutput$importance-nativeOutput$performance[ntree]
                      else NULL),
        forest = forest,
        proximity = (if (proximity) nativeOutput$proximity else NULL)
    )
    if (!big.data) {
      class(rsfOutput) <- c("rsf", "grow")
    }
    else {
      class(rsfOutput) <- c("rsf", "grow", "bigdata")
    }
    return(rsfOutput)
    
  }

