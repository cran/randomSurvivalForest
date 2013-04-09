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

max.subtree.rsf <- function(object, max.order=2, sub.order=FALSE, ...) {
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2)
  stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    object <- object$forest
  }
  nativeArray <- object$nativeArray
  if (is.null(nativeArray)) {
    stop("RSF nativeArray content is NULL.  Please ensure the object is valid.")
  }
  predictorNames <- object$predictorNames
  if (is.null(predictorNames)) {
    stop("RSF predictorNames content is NULL.  Please ensure the object is valid.")
  }
  if (is.null(object$predictors)) {
    stop("RSF predictors content is NULL.  Please ensure the object is valid.")
  }
  max.order <- floor(max.order)
  if (max.order < 0) {
    stop("RSF 'max.order' requested for distance order statistic must be an integer greater than zero (0).")
  }
  MAX.DEPTH <- 100
  numTree <- length(as.vector(unique(nativeArray$treeID)))
  numParm <- length(predictorNames)
  subtree <- vector("list", 8)
  names(subtree) <- c("count",
                      "order",
                      "meanSum",
                      "depth",
                      "terminalDepthSum",
                      "subOrder",
                      "subOrderDiag",
                      "nodesAtDepth")
  forestMeanSum     <- rep(0, numParm)
  if (max.order > 0) {
    orderSum    <- matrix(0, nrow=numParm, ncol=max.order)
  }
  else {
    order.tree    <- matrix(NA, nrow=numParm, ncol=numTree)
  }
  recursiveObject <- list(offset = 1,
                          subtree = subtree,
                          diagnostic = 0,
                          diagnostic2 = 0)
  subtreeCountSum  <- rep(0, numParm)
  subOrderSum <- matrix(0, nrow=numParm, ncol=numParm)
  terminalDepth <- rep(0, numTree)
  nodesAtDepthMatrix <- matrix(NA, nrow = MAX.DEPTH, ncol = numTree)
  offsetMark <- 1
  stumpCnt <- 0
  for (b in 1:numTree) {
    recursiveObject$subtree$nodesAtDepth <- rep(NA, MAX.DEPTH)
    recursiveObject$subtree$meanSum <- rep(NA, numParm)
    if (max.order > 0) {
      recursiveObject$subtree$order <- matrix(NA, nrow=numParm, ncol=max.order)
    }
    else {
      recursiveObject$subtree$order <- rep(NA, numParm)
    }
    if (sub.order ==TRUE) {
      recursiveObject$subtree$subOrder <- matrix(1.0, nrow=numParm, ncol=numParm)
      recursiveObject$subtree$subOrderDiag <- rep(NA, numParm)
    }
    recursiveObject$subtree$depth <- 0
    recursiveObject$subtree$terminalDepthSum <- 0
    recursiveObject$subtree$count <- rep(0, numParm)
    rootParmID <- nativeArray$parmID[recursiveObject$offset] 
    offsetMark <- recursiveObject$offset
    recursiveObject <- rsfParseTree(
      recursiveObject,
      max.order,
      sub.order,
      nativeArray,
      b,
      distance=0,
      subtreeFlag=rep(FALSE, numParm))
    if (rootParmID != 0) {
      index <- which(recursiveObject$subtree$count == 0)
      recursiveObject$subtree$meanSum[index] <- recursiveObject$subtree$depth
      forestMeanSum <- forestMeanSum + recursiveObject$subtree$meanSum
      if (max.order > 0) {
        index <- which(is.na(recursiveObject$subtree$order))
        recursiveObject$subtree$order[index] <- recursiveObject$subtree$depth
        orderSum   <- orderSum + recursiveObject$subtree$order
      }
      else {
        index <- which(is.na(recursiveObject$subtree$order))
        recursiveObject$subtree$order[index] <- recursiveObject$subtree$depth
        order.tree[ , b] <- recursiveObject$subtree$order
      }
      subtreeCountSum <- subtreeCountSum + (recursiveObject$subtree$count / ((recursiveObject$offset - offsetMark + 1) / 4))
      terminalDepth[b] <- recursiveObject$subtree$terminalDepthSum / ((recursiveObject$offset - offsetMark + 1) / 2)
      if (sub.order == TRUE) {
        index <- which(recursiveObject$subtree$count > 0)
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$subOrderDiag[index]
        index <- which(recursiveObject$subtree$count == 0)
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$depth
        diag(recursiveObject$subtree$subOrder) <-  diag(recursiveObject$subtree$subOrder) / recursiveObject$subtree$depth
        subOrderSum <- subOrderSum + recursiveObject$subtree$subOrder
      }
      nodesAtDepthMatrix[, b] <- recursiveObject$subtree$nodesAtDepth
    }
    else {
      stumpCnt <- stumpCnt + 1
      nodesAtDepthMatrix[, b] <- NA      
    }
  }  
  nameVector <- c("mean",
                  "order",
                  "count",
                  "terminal",
                  "nodesAtDepth",
                  "subOrder")
  result <- vector("list", length(nameVector))
  names(result) <- nameVector
  if(numTree != stumpCnt) {
    result$terminal <- terminalDepth
    result$mean <- forestMeanSum / (numTree - stumpCnt)
    names(result$mean) <- predictorNames
    if (max.order > 0) {
      result$order <- orderSum / (numTree - stumpCnt)
      rownames(result$order) <- predictorNames
    }
    else {
      result$order <- order.tree
      rownames(result$order) <- predictorNames
    }
    result$count <- subtreeCountSum / (numTree - stumpCnt)
    names(result$count) <- predictorNames
    result$nodesAtDepth <- nodesAtDepthMatrix
    if (sub.order == TRUE) {      
      result$subOrder   <- subOrderSum / (numTree - stumpCnt)
      rownames(result$subOrder) <- predictorNames
      colnames(result$subOrder) <- predictorNames
    }
  }
  threshold <- ExactThreshold(result)
  result <- c(result, threshold=threshold)
  return (result)
}
rsfParseTree <- function(recursiveObject,
                         max.order,
                         sub.order,
                         nativeArray,
                         b,
                         distance,
                         subtreeFlag) {
  recursiveObject$diagnostic <- recursiveObject$diagnostic + 1
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  if (distance > 0) {
    if (distance <= length(recursiveObject$subtree$nodesAtDepth)) {    
      if (is.na(recursiveObject$subtree$nodesAtDepth[distance])) {
        recursiveObject$subtree$nodesAtDepth[distance] <- 1
      }
      else {
        recursiveObject$subtree$nodesAtDepth[distance] <- recursiveObject$subtree$nodesAtDepth[distance] + 1
      }
    }
  }
  splitParameter <- nativeArray$parmID[recursiveObject$offset]
  if (splitParameter == 0) {
    terminalFlag <- TRUE
  }
  else if (splitParameter != 0) {
    terminalFlag <- FALSE
  }
  if (!terminalFlag) {
    if (subtreeFlag[splitParameter] == FALSE) {
      recursiveObject$subtree$count[splitParameter] <- recursiveObject$subtree$count[splitParameter] + 1
      if (is.na(recursiveObject$subtree$meanSum[splitParameter])) {
        recursiveObject$subtree$meanSum[splitParameter] <- distance
      }
      else {
        recursiveObject$subtree$meanSum[splitParameter] <- recursiveObject$subtree$meanSum[splitParameter] + distance
      }
      if (max.order > 0) {
        orderVector <- c(recursiveObject$subtree$order[splitParameter, ], distance, NA)
        index <- which.max(is.na(orderVector))
        orderVector[index] <- distance
        sortedVector <- sort(orderVector[1:index])
        if (index <= max.order) {
          orderVector <- c(sortedVector, rep(NA, max.order-index))
        }
        else {
          orderVector <- sortedVector[1:max.order]
        }
        recursiveObject$subtree$order[splitParameter, ] <- orderVector
      }
      else {
        if (is.na(recursiveObject$subtree$order[splitParameter])) {
          recursiveObject$subtree$order[splitParameter] <- distance
        }
        else {
          recursiveObject$subtree$order[splitParameter] <- min(recursiveObject$order[splitParameter], distance)
        }
      }
      subtreeFlag[splitParameter] <- TRUE
      if (sub.order == TRUE) {
        if (is.na(recursiveObject$subtree$subOrderDiag[splitParameter])) {
          recursiveObject$subtree$subOrderDiag[splitParameter] <- distance
        }
        else {
          recursiveObject$subtree$subOrderDiag[splitParameter] <- min(recursiveObject$subtree$subOrderDiag[splitParameter], distance)
        }
        recursive2Object <- list(offset = recursiveObject$offset,
          depth = 0,
          minimumVector = rep(NA, dim(recursiveObject$subtree$subOrder)[2]),
          diagnostic = recursiveObject$diagnostic2)
        subtree2Flag <- rep(FALSE, dim(recursiveObject$subtree$subOrder)[2])        
        subtree2Flag[splitParameter] <- TRUE
        recursive2Object <- rsfParse2Tree(recursive2Object,
          nativeArray,
          b,
          distance=0,
          subtreeFlag=subtree2Flag)
        recursiveObject$diagnostic2 <- recursiveObject$diagnostic2 + recursive2Object$diagnostic
        recursive2Object$minimumVector[splitParameter] <- recursive2Object$depth
        recursive2Object$minimumVector[which(is.na(recursive2Object$minimumVector))] <- recursive2Object$depth
        recursive2Object$minimumVector <- recursive2Object$minimumVector / recursive2Object$depth
        recursiveObject$subtree$subOrder[splitParameter, ] <- pmin(recursiveObject$subtree$subOrder[splitParameter, ], recursive2Object$minimumVector)
      }  
    }  
  }  
  recursiveObject$subtree$depth <- max(recursiveObject$subtree$depth, distance)
  recursiveObject$offset <- recursiveObject$offset + 1
  if (terminalFlag == FALSE) {
    distance <- distance + 1
    recursiveObject <- rsfParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)
    recursiveObject <- rsfParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)
  }
  else {
    recursiveObject$subtree$terminalDepthSum <- recursiveObject$subtree$terminalDepthSum + distance
  }
  return (recursiveObject)
}
rsfParse2Tree <- function(recursiveObject,
                          nativeArray,
                          b,
                          distance,
                          subtreeFlag) {
  recursiveObject$diagnostic = recursiveObject$diagnostic + 1
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  splitParameter = nativeArray$parmID[recursiveObject$offset]
  if (splitParameter == 0) {
    terminalFlag = TRUE
  }
  else if (splitParameter != 0) {
    terminalFlag = FALSE
  }
  if (splitParameter != 0) {
    if (subtreeFlag[splitParameter] == FALSE) {
      if (is.na(recursiveObject$minimumVector[splitParameter])) {
        recursiveObject$minimumVector[splitParameter] = distance
      }
      else {
        recursiveObject$minimumVector[splitParameter] = min(recursiveObject$minimumVector[splitParameter], distance)
      }
      subtreeFlag[splitParameter] = TRUE
    }
  }
  recursiveObject$depth = max(recursiveObject$depth, distance)
  distance = distance + 1
  recursiveObject$offset = recursiveObject$offset + 1
  if (terminalFlag == FALSE) {
    recursiveObject = rsfParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)
    recursiveObject = rsfParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)
  }
  return (recursiveObject)
}
  maxDepthProb <- function(p, D, l) {
    if (!is.null(l)) Ld <- 0
    prob <- rep(0, D+1)
    for (d in 0:(D-1)) {
      if (is.null(l)) {
        Ld <- 2^d-1
        ld <- 2^d
      }
      else{
       ld <- l[d+1]
       if (d > 0) Ld <- Ld + l[d] 
      }
      prob.d.1 <- Ld*log(1-1/p)
      prob.d.2 <- ld*(log(1-1/p))
      prob[d+1] <- exp(prob.d.1)*(1-exp(prob.d.2))
    }
    prob[D+1] = 1 - sum(prob[1:D])
    if (prob[D+1] < 0) {
      prob[D+1] <- 0
      prob <- prob/sum(prob)
    }
      prob
  }
  maxDepthStat <- function(pseq, D=NULL, l=NULL) {
    mn <- std <- rep(0, length(pseq))
    if (is.null(D) & is.null(l)) stop("set D or l")
    if (!is.null(l)) {
      D <- length(l)
    }
    D.support <- (0:D)
      for (j in 1:length(pseq)) {
        prob <- maxDepthProb(pseq[j], D=D, l=l)
        mn[j] <- sum(D.support*prob)
        std[j] <- sqrt(sum((D.support^2)*prob) - mn[j]^2)
      }
    return(list(mean=mn, std=std))
  }
  ExactThreshold <- function(v) {
    if (is.null(v$mean)) return(NULL)
    n.at.d <- round(c(1, c(na.omit(apply(v$nodesAtDepth, 1, mean, na.rm=TRUE)))))
    avg.depth <- round(mean(apply(v$nodesAtDepth, 2, function(x){sum(!is.na(x))}), na.rm=TRUE))
    l <- n.at.d[1:max(avg.depth - 1, 1)]
    if (length(v$mean) == 1) {
      return(0)
    }
    else {
      return(maxDepthStat(length(v$mean), l=l)$mean)
    }
  }
max.subtree <-  max.subtree.rsf
