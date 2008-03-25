##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.2.2
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

interaction.rsf <- function(
                            object = NULL,
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

  ## Map interaction names to columns of predictor matrix
  ## Map interaction terms and the (potentially proper) subset of
  ## the GROW data to the native code requirements.
  predictorNames <- unique(predictorNames)
  predictorNames <- apply(cbind(1:length(predictorNames)), 1, function(i){
                   which(is.element(colnames(object$predictors), predictorNames[i]))})
  if (length(predictorNames)==0 || is.null(predictorNames)) {
    ## The user has not specified any variable names.  Exit immediately. 
    stop("Variable names must be specified.")
  }
  else {
    if (is.null(subset)) {
      ## The user has not specified a subset of the GROW data.
      ## Assume the entire data set is to be used.
      subset <- 1:dim(object$predictors)[1]
    }
  }
  subset <- unique(subset[subset >= 1
                         & subset <= dim(object$predictors)[1]])
  if (length(subset) == 0) stop("'subset' not set properly.")
  
  ## Get names
  fNames <- all.vars(object$formula, max.names=1e7)
  ntree <- length(unique(object$nativeArray[,1]))

  ## Set the predictor types.
  if (inherits(object, c("rsf", "grow", "bigdata"), TRUE) [3] == 3) big.data <- T else big.data <- F
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

  ## Generate the random seed.
  if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
  seed <- -round(abs(seed))
  
  ## Set the trace level.
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

  ## Convert importance option into native code parameter.
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
    ## Use the mean survival time as predicted outcome.
    predictedOutcomeBits <- 2^14
  }
  else if (rough == FALSE) {
    ## Use the CHF as predicted outcome.
    predictedOutcomeBits <- 0    
  }
  else {
    stop("Invalid choice for 'rough' option:  ", rough)
  }
    
  
  ## #####################################################################
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
  ## 15 - vector representing spltPT
  ## 16 - head of seed chain, integer < 0
  ## 17 - vector of predictor types
  ##        "R" - real value
  ##        "I" - integer value 
  ##        "C" - categorical
  ## 18 - number of predictors in interaction, integer > 0
  ## 19 - vector of predictors in interaction, integers > 0
  ## 20 - number of individuals in the (potentially proper)
  ##      subsetted GROW data set, integer > 1
  ## 21 - index of individuals in the subsetted GROW data set,
  ##      integers > 1
  ## ###########################################################

  ## ###########################################################
  ##  SEXP outputs:  See parameter 2 above and note "*".
  ## ###########################################################    

  nativeOutput <- .Call("rsfInteraction",
                        as.integer(do.trace),
                        as.integer(predictedOutcomeBits +
                                   importanceBits),
                        as.integer(seed),
                        as.integer(ntree),
                        as.integer(dim(object$predictors)[1]),
                        as.double(object$time),
                        as.double(object$cens),
                        as.integer(dim(object$predictors)[2]),
                        as.numeric(object$predictors),
                        as.integer(length(object$timeInterest)),
                        as.double(object$timeInterest),
                        as.integer(object$nativeArray[,1]),
                        as.integer(object$nativeArray[,2]),
                        as.integer(object$nativeArray[,3]),
                        as.double(object$nativeArray[,4]),
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

  rsfOutput <- list(
                    err.rate = nativeOutput$performance,
                    importance = nativeOutput$importance-nativeOutput$performance[ntree]
                    )

 return(rsfOutput)

}
