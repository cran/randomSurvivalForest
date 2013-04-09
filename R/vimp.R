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

vimp <- function(object = NULL,#
                 predictorNames = NULL,#
                 subset = NULL,#
                 joint = TRUE,#
                 rough = FALSE,#                            
                 importance = c("randomsplit", "permute", "none")[1],#
                 seed = NULL,#
                 do.trace = FALSE,#
                 ...)
{
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'.")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    object <- object$forest
  }
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
  if (is.null(subset)) {
    subset <- 1:nrow(object$predictors)
  }
  else {
    subset <- unique(subset[subset >= 1
                         & subset <= nrow(object$predictors)])
    if (length(subset) == 0) stop("'subset' not set properly.")
  }
  ntree <- length(unique(object$nativeArray[,1]))
  get.factor <- extract.factor(object$predictors, object$predictorNames)
  xfactor <- get.factor$xfactor
  xfactor.order <- get.factor$xfactor.order
  predictorType <- get.factor$predictorType
  object$predictors <- as.matrix(data.matrix(object$predictors))
  rownames(object$predictors) <- colnames(object$predictors) <- NULL
  if (is.null(seed) || abs(seed)<1) seed <- runif(1,1,1e6)
  seed <- -round(abs(seed))
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
  }
  else {
    stop("Invalid choice for 'joint' option:  ", joint)
  }
  if (rough == TRUE) {
    predictedOutcomeBits <- 2^14
  }
  else if (rough == FALSE) {
    predictedOutcomeBits <- 0    
  }
  else {
    stop("Invalid choice for 'rough' option:  ", rough)
  }
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
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }
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
  rsfOutput <- list(
                    err.rate = make.vec(err.rate),
                    err.perturb.rate = make.vec(err.perturb.rate),
                    importance = make.vec(importance)
                    )
 return(rsfOutput)
}
