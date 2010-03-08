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

vimp <- function(object = NULL,
                 predictorNames = NULL,
                 subset = NULL,
                 joint = TRUE,
                 rough = FALSE,                            
                 importance = c("randomsplit", "permute", "none")[1],
                 seed = NULL,
                 do.trace = FALSE,
                 ...) {

  ## Check that 'object' is of the appropriate type.
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'.")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    object <- object$forest
  }

  ## Map interaction names to columns of GROW predictor matrix
  ## Ensure interaction names are coherent
  predictorNames <- unique(predictorNames)
  if (is.null(predictorNames)) {
    predictorNames <- object$predictorNames
  }
  else {
    predictorNames <- intersect(predictorNames, object$predictorNames)
  }
  if (length(predictorNames) == 0)
    stop("predictor names do not match object predictor matrix")
  predictorNames <- apply(cbind(1:length(predictorNames)), 1, function(i) {
                   which(is.element(object$predictorNames, predictorNames[i]))})

  ## The user has not specified a subset of the GROW data.
  ## Assume the entire data set is to be used.
  if (is.null(subset)) {
    subset <- 1:nrow(object$predictors)
  }
  else {
    subset <- unique(subset[subset >= 1
                         & subset <= nrow(object$predictors)])
    if (length(subset) == 0) stop("'subset' not set properly.")
  }
  
  ## Get information from the grow object
  ## Identify unordered/ordered factors
  ## Determine predictor types
  ## Data conversion to numeric
  ntree <- length(unique(object$nativeArray[,1]))
  get.factor <- extract.factor(object$predictors, object$predictorNames)
  xfactor <- get.factor$xfactor
  xfactor.order <- get.factor$xfactor.order
  predictorType <- get.factor$predictorType
  object$predictors <- as.matrix(data.matrix(object$predictors))
  rownames(object$predictors) <- colnames(object$predictors) <- NULL
  
  ## Generate the random seed
  if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
  seed <- -round(abs(seed))
  
  ## Set the trace level
  if (!is.logical(do.trace)) {
    if (do.trace >= 1){
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
  if (joint == TRUE) {
    importanceBits <- importanceBits + 2^10
  }
  else if (joint == FALSE) {
    ## Do nothing.  All is well.
  }
  else {
    stop("Invalid choice for 'joint' option:  ", joint)
  }
  if (rough == TRUE) {
    ## Use the mean survival time as predicted outcome
    predictedOutcomeBits <- 2^14
  }
  else if (rough == FALSE) {
    ## Use the CHF as predicted outcome
    predictedOutcomeBits <- 0    
  }
  else {
    stop("Invalid choice for 'rough' option:  ", rough)
  }


  ## convert outcome (if it is present) into opt bit
  if (is.null(object$outcome)) {
    outcomeBits <- 0
  }
  else {
    if (object$outcome == "test") {
      outcomeBits <- 2^17
    }
    else {
      outcomeBits <- 0
    }
  }

  ########################################################################
  ## rsfInteraction(...)
  ##
  ## Parameters passed:
  ## #####################################################################
  ## 00 - C function call
  ##
  ## 01 - trace output flag 
  ##        0 = no trace output
  ##       !0 = various levels of trace output
  ##
  ## 02 - option protocol for output objects.
  ##    - see primary GROW call for details
  ##
  ## 03 - random seed for repeatability, integer < 0 
  ## 04 - number of trees in the forest, integer > 0
  ## 05 - number of individuals in GROW data set, integer > 1
  ## 06 - vector of GROW observed times of death, doubles > 0 
  ## 07 - vector of GROW observed event types
  ##        0 = censored
  ##        1 = death
  ## 08 - number of GROW predictors, integer > 0
  ## 09 - [p x n] matrix of GROW predictor observations
  ## 10 - number of time points of interest, integer > 0
  ## 11 - vector of time points of interest, doubles
  ## 12 - vector representing treeID
  ## 13 - vector representing nodeID
  ## 14 - vector representing parmID
  ## 15 - vector representing contPT
  ## 16 - vector representing mwcpSZ
  ## 17 - vector representing mwcpPT  
  ## 18 - head of seed chain, integer < 0
  ## 19 - vector of predictor types
  ##        "R" - real value
  ##        "I" - integer value 
  ##        "C" - categorical
  ## 20 - number of predictors in interaction, integer > 0
  ## 21 - vector of predictors in interaction, integers > 0
  ## 22 - number of individuals in the (potentially proper)
  ##      subsetted GROW data set, integer > 1
  ## 23 - index of individuals in the subsetted GROW data set,
  ##      integers > 1
  ## ###########################################################

  ## ###########################################################
  ##  SEXP outputs:  See parameter 2 above and note "*".
  ## ###########################################################    

  nativeOutput <- .Call("rsfInteraction",
                        as.integer(do.trace),
                        as.integer(predictedOutcomeBits +
                                   importanceBits +
                                   outcomeBits),
                        as.integer(seed),
                        as.integer(ntree),
                        as.integer(nrow(object$predictors)),
                        as.double(object$time),
                        as.double(object$cens),
                        as.integer(ncol(object$predictors)),
                        as.numeric(object$predictors),
                        as.integer(length(object$timeInterest)),
                        as.double(object$timeInterest),
                        as.integer((object$nativeArray)$treeID),
                        as.integer((object$nativeArray)$nodeID),
                        as.integer((object$nativeArray)$parmID),
                        as.double((object$nativeArray)$contPT),
                        as.integer((object$nativeArray)$mwcpSZ),                        
                        as.integer(object$nativeFactorArray),                        
                        as.integer(object$seed),
                        as.character(predictorType),
                        as.integer(length(predictorNames)),
                        as.integer(predictorNames),
                        as.integer(length(subset)),
                        as.integer(subset))
  
  ## Check for error return condition in the native code.
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }

  # number of event types
  # pretty names
  n.event <- length(unique(na.omit(object$cens)[na.omit(object$cens) > 0]))
  if (n.event > 1) n.event <- n.event + 1
  if (joint) j.names <- "joint" else j.names <- object$predictorNames[predictorNames]
  if (n.event > 1) {    
    err.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), NULL)
    vimp.names <- list(c("CHF", paste("condCHF.", 1:(n.event - 1), sep = "")), j.names)
  }
  else {
    err.names <- list(NULL, NULL)
    vimp.names <- list(NULL, j.names)
  }

  # error rates
  err.rate <- matrix(nativeOutput$performance, ncol=ntree, byrow=T,
                      dimnames=err.names)[1:n.event, ntree, drop = FALSE]
  colnames(err.rate) <- "all.vars"
  if(!is.null(nativeOutput$importance)) {
    err.perturb.rate <- matrix(nativeOutput$importance, ncol=max(1, length(j.names)), 
                             byrow = TRUE, dimnames = vimp.names)[1:n.event,, drop = FALSE]
    importance <- t(sweep(t(err.perturb.rate), 2, err.rate))
  }
  else {
    err.perturb.rate <- importance <- NULL
  }

  ## for pretty output
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
  # make.vec <- function(x) {if (!is.null(dim(x)) && nrow(x) == 1) c(x) else x}
  
  rsfOutput <- list(
                    err.rate = make.vec(err.rate),
                    err.perturb.rate = make.vec(err.perturb.rate),
                    importance = make.vec(importance)
                    )

 return(rsfOutput)

}
