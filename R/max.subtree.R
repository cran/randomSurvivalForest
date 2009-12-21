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

#########################################################################
##
##  Maximal Subtree and Supporting Functions
##
#########################################################################

### User callable function:

max.subtree <- function(object,
                        max.order=2,
                        sub.order=FALSE,
                        ...) {

  ## Incoming parameter checks.  All are fatal.
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2)
  stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'")

  ## Acquire the forest 
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

  ## Maximum depth monitored for nodes at depth counts. 
  MAX.DEPTH <- 100

  ## Count the number of trees in the forest.
  numTree <- length(as.vector(unique(nativeArray$treeID)))

  ## Count the number of parameters in the data set.
  numParm <- length(predictorNames)

  ## Create the (local) subtree object for recursion.  This is NOT the
  ## output object, but is closely related to it.

  subtree <- vector("list", 8)
  names(subtree) <- c("count",
                      "order",
                      "meanSum",
                      "depth",
                      "terminalDepthSum",
                      "subOrder",
                      "subOrderDiag",
                      "nodesAtDepth")

  ## Initialize the (local) subtree object for recursion.  This is NOT
  ## the output object, but is closely related to it.  This is
  ## SPECIFIC to a single tree.  Note that the maximal distance from
  ## the root node of a tree or sub-tree is zero-based.  This means
  ## that the root node represents a distance of zero (0).
  ##
  ## $count   - of length [numParm].
  ##          - integer count of the number of maximal subtrees for each
  ##            predictor in each tree.
  ##          - this will be zero (0) if the predictor is not used in
  ##            the forest.
  ## $order   - if max.order > 0, matrix of dim [numParm] x [max.order]
  ##            representing the order statistic for the maximal
  ##            subtree distance for each predictor
  ##          - if max.order == 0, matrix of [numParm] x [numTree]
  ##            representing the minimum maximal subtree distance for
  ##            each predictor
  ##          - rows represent the distribution for the order
  ##            statistic for a predictor
  ##          - thus [2][1] gives the minimum distance for predictor 2.
  ##          - elements will have the value (maximum + 1), the
  ##            penalty, if the predictor is not used in the tree, or
  ##            the statistic does not exist
  ## $meanSum - vector of length [numParm]
  ##          - the pre-mean mimimum maximal subtree distance for each
  ##            predictor in each tree.
  ##          - if the predictor is not used in the tree, this will be
  ##            the maximum depth for the tree, the penalty.
  ## $depth   - integer depth (defined as the maximum terminal node
  ##            distance from the root node) possible in the tree
  ##            topology.  This will be always greater than zero (0)
  ##            for non-trivial trees
  ## $terminalDepthSum
  ##          - integer sum of the depths of the terminal
  ##            nodes of the tree
  ##          - used in calculating the average depth of the
  ##            nodes for tree
  ## $subOrder
  ##          - matrix of [numParm] x [numParm] representing the
  ##            normalized minimum maximal w-subtree distance for
  ##            parameter v's maximal subtree.
  ##          - thus if v=2, w=1, consider a maximal 2-subtree with
  ##            the "root" node of this subtree associated with
  ##            relative depth = 0.  Then consider a 1-subtree of the
  ##            2-subtree.  Then, [2][1] gives the normalized minimum
  ##            maximal subtree distance of all 1-subtrees within the
  ##            2-subtree.  If the maximal w-subtree does not exist,
  ##            the distance is set to be one (1).  The normalized
  ##            distance is calculated as follows.  First, the
  ##            maximum depth of the terminal nodes of the v-subtree
  ##            is determined.  Then the relative depth of the
  ##            w-subtree (from the root node) of the v-subtree is
  ##            determined.  This latter quantity is divided by the
  ##            prior maximum terminal node depth quantity to give a
  ##            number between (0,1].
  ##          - the diagonal [i][i] is the normalized minimum maximal
  ##            v-subtree distance for the parameter i.  If it is not
  ##            split on in the tree, it is assigned the value one
  ##            (1).
  ##          - the resulting matrix is NA-free.
  ## $subOrderDiag
  ##          - interim vector of length [numParm] representing the 
  ##            minimum maximal v-subtree distance for the tree that
  ##            will constitute the diagonal of $subOrder.
  ##          - variables that are not split on are assigned the
  ##            penalty value of the maximum depth.
  ## $nodesAtDepth
  ##          - number of nodes at each depth, where [1] = 2 for
  ##            a non-trivial tree, [2] = 4 for a balanced tree.
  ##            In general, for a balanced tree, at depth d, the
  ##            count is 2^d.
  ##

  ## Pre-mean maximal subtree distance by parameter over the forest.
  forestMeanSum     <- rep(0, numParm)

  if (max.order > 0) {
    ## Order statistics by parameter of maximal subtree distance. 
    orderSum    <- matrix(0, nrow=numParm, ncol=max.order)
  }
  else {
    ## This is by tree over the forest.
    ## Minimum maximal subree distance by parameter by tree over the forest.
    order.tree    <- matrix(NA, nrow=numParm, ncol=numTree)
  }

  ## Create the recursive output object.  This would be unnecessary 
  ## if it was possible to declare global variables in a package.
  recursiveObject <- list(offset = 1,
                          subtree = subtree,
                          diagnostic = 0,
                          diagnostic2 = 0)

  ## Normalized maximal subtree count sum over the forest.  Linked to $count.
  subtreeCountSum  <- rep(0, numParm)

  ## Sub-order interim sum.
  subOrderSum <- matrix(0, nrow=numParm, ncol=numParm)

  ## Average depth of each tree.  Linked to $terminalDepthSum 
  terminalDepth <- rep(0, numTree)

  ## Final list of count of nodes at depth by tree.
  nodesAtDepthMatrix <- matrix(NA, nrow = MAX.DEPTH, ncol = numTree)

  offsetMark <- 1
  stumpCnt <- 0
  
  ## Loop through all trees.
  for (b in 1:numTree) {

    ## Global dependencies:  (predictorNames, forest)

    ## Reset the nodes at depth count.
    recursiveObject$subtree$nodesAtDepth <- rep(NA, MAX.DEPTH)
    
    ## Reset the mean maximal subtree distance.
    ## This will be NA if the predictor is not used in this tree.
    recursiveObject$subtree$meanSum <- rep(NA, numParm)
      
    ## Reset the order statistic matrix representing maximal distance.
    if (max.order > 0) {
      ## This is tree specific.  This will be NA if the predictor is
      ## not used in this tree.
      recursiveObject$subtree$order <- matrix(NA, nrow=numParm, ncol=max.order)
    }
    else {
      ## This is tree specific.  This represents the column in
      ## order.tree corresponding to this tree.  We don't want to drag
      ## the entire data structure order.tree up and down the stack,
      ## so we subset it.
      recursiveObject$subtree$order <- rep(NA, numParm)
    }


    if (sub.order ==TRUE) {
      ## Reset the sub-order statistic matrix representing the minimum
      ## maximal w-subtree distance for v != w.  
      recursiveObject$subtree$subOrder <- matrix(1.0, nrow=numParm, ncol=numParm)

      ## Reset the subOrder diagonals to NA.
      recursiveObject$subtree$subOrderDiag <- rep(NA, numParm)
    }
  
    
    ## Reset the maximum maximal subtree distance.
    recursiveObject$subtree$depth <- 0

    ## Reset the terminal node depth sums.
    recursiveObject$subtree$terminalDepthSum <- 0
    
    ## Reset the count of the maximal subtrees in the tree.
    recursiveObject$subtree$count <- rep(0, numParm)
    
    ## Identify the root node split parameter.
    rootParmID <- nativeArray$parmID[recursiveObject$offset] 

    ## Save the previous value of the offset.  This is used in
    ## determining the number of terminal nodes in the tree.
    offsetMark <- recursiveObject$offset

    
    ## ----------------------------
    ##  PRIMARY RECURSION PRIMARY
    ## ----------------------------
    
    ## Recursively parse the tree in the primary protocol.
    recursiveObject <- rsfParseTree(
      recursiveObject,
      max.order,
      sub.order,
      nativeArray,
      b,
      distance=0,
      subtreeFlag=rep(FALSE, numParm))

    ## Check that the current tree is not a stump.
    if (rootParmID != 0) {

      ## Update the forest sum ...

      ## Determine which predictors _were_not_ split on, in the current tree.
      index <- which(recursiveObject$subtree$count == 0)

      ## Penalize these unused predictors by making their mean distance the worst.
      recursiveObject$subtree$meanSum[index] <- recursiveObject$subtree$depth
      forestMeanSum <- forestMeanSum + recursiveObject$subtree$meanSum
      
      if (max.order > 0) {
        ## This is tree specific.  Here we view the matrix as a vector.
        index <- which(is.na(recursiveObject$subtree$order))
        ## Penalize these order statistics making their values the worst.      
        recursiveObject$subtree$order[index] <- recursiveObject$subtree$depth

        orderSum   <- orderSum + recursiveObject$subtree$order
        
      }
      else {
        ## This is tree specific.  Here we view modify the element corresponding to this tree only.
        index <- which(is.na(recursiveObject$subtree$order))
        ## Penalize these order statistics making their values the worst.      
        recursiveObject$subtree$order[index] <- recursiveObject$subtree$depth

        ## Copy the order statistic for this tree to the persistent matrix.
        order.tree[ , b] <- recursiveObject$subtree$order
        
      }

      ## Explanation of demominator: The total number of nodes N
      ## (internal and external) in a tree is given by the number of
      ## records in a tree.  In this case, offset - offsetMark = N,
      ## due to the nature of the pointer incrementation.  The number
      ## of external (terminal) nodes N(TERM) is given by (N + 1) / 2.
      ## The maximum number of subtrees for the tree is given by
      ## N(TERM) / 2 which is the denominator.

      ## We normalize the actual number of maximal subtrees found in
      ## the tree by the maximum number maximal subtrees possible.
      subtreeCountSum <- subtreeCountSum + (recursiveObject$subtree$count / ((recursiveObject$offset - offsetMark + 1) / 4))

      ## Get the average depth of the terminal nodes of this tree.
      ## The denominator is just the number of terminal nodes.
      terminalDepth[b] <- recursiveObject$subtree$terminalDepthSum / ((recursiveObject$offset - offsetMark + 1) / 2)
      
      if (sub.order == TRUE) {
        ## Determine which predictors _were_ split on, in the current tree.
        index <- which(recursiveObject$subtree$count > 0)

        ## Over-ride those diagonal elements with the minimum depth encountered.
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$subOrderDiag[index]

        ## Determine which predictors _were_not_ split on, in the
        ## current tree.  This can also be determined by examining the
        ## diagonal vector, and noting those elements with the value
        ## NA
        index <- which(recursiveObject$subtree$count == 0)
        
        ## Over-ride those missing diagonal elements with penalty value.
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$depth

        ## Normalize the depths.
        diag(recursiveObject$subtree$subOrder) <-  diag(recursiveObject$subtree$subOrder) / recursiveObject$subtree$depth
        
        ## Add the result to the forest sum.
        subOrderSum <- subOrderSum + recursiveObject$subtree$subOrder
        
      }

      ## Update the nodes at depth count.
      nodesAtDepthMatrix[, b] <- recursiveObject$subtree$nodesAtDepth
      
    }
    else {
      ## The tree is a stump.  It will be excluded.
      stumpCnt <- stumpCnt + 1

      ## The nodes at depth count will already be NA.
      nodesAtDepthMatrix[, b] <- NA      
    }
    
  }  ## for (b in 1:subtree) ...


  ## --------------------
  ## Output preparation  
  ## --------------------
  
  ## Prepare the ensemble values for the output object.

  ## $mean - vector of length [numParm]
  ##       - mean minimum maximal subtree distance by parameter over the
  ##         forest.
  ## $order - if max.order > 0, matrix of dim [numParm] x [max.order]
  ##            representing the order statistic for the maximal
  ##            subtree distance for each predictor over the forest
  ##          - if max.order == 0, matrix of [numParm] x [numTree]
  ##            representing the minimum maximal subtree distance for
  ##            each predictor over the forest
  ## $count - vector of length [numParm]
  ##        - normalized maximal subtree count by parameter over the
  ##          forest
  ##        - normalization is by the number of topologically possible
  ##          subtrees in each tree
  ## $subOrder - matrix of [numParm] x [numParm] representing the
  ##             normalized minimum maximal w-subtree distance for
  ##             parameter v's maximal subtree averaged over all trees.
  ##           - see the tree specific description for details.
  ##           - output if sub.order=TRUE
  ## $terminal - vector of length [numTree]
  ##           - average terminal node depth of each tree.
  ## $nodesAtDepth - matrix of [MAX.DEPTH] x [numTree]
  ##               - number of nodes at each depth, where [1] = 2 for
  ##                 a non-trivial tree, [2] = 4 for a balanced tree.
  ##                 In general, for a balanced tree, at depth d, the
  ##                 count is [d] = 2^d.

  
  nameVector <- c("mean",
                  "order",
                  "count",
                  "terminal",
                  "nodesAtDepth",
                  "subOrder")
  
  result <- vector("list", length(nameVector))
  names(result) <- nameVector

  ## Precautionary check for all stumps.
  if(numTree != stumpCnt) {

    ## Attach the average terminal node depth to the output object.
    result$terminal <- terminalDepth
    
    ## Determine the forest average maximal subtree distance.
    result$mean <- forestMeanSum / (numTree - stumpCnt)
    names(result$mean) <- predictorNames
    
    if (max.order > 0) {
      ## Summarize the order statistics containing the maximal subtree distances. 
      result$order <- orderSum / (numTree - stumpCnt)
      rownames(result$order) <- predictorNames
    }
    else {
      result$order <- order.tree
      rownames(result$order) <- predictorNames
    }
    
    ## Determine the forest average subtree count.
    result$count <- subtreeCountSum / (numTree - stumpCnt)
    names(result$count) <- predictorNames

    ## Copy the matrix for the number of nodes at each depth.
    result$nodesAtDepth <- nodesAtDepthMatrix
    

    if (sub.order == TRUE) {      
      ## Calculate the forest average of the minimum maximal w-subtree distances.
      result$subOrder   <- subOrderSum / (numTree - stumpCnt)
      rownames(result$subOrder) <- predictorNames
      colnames(result$subOrder) <- predictorNames
    }
  }

  ## print (recursiveObject$diagnostic)
  ## print (recursiveObject$diagnostic2)


  ## Minimal depth for thresholding weak variables
  threshold <- ExactThreshold(result)
  result <- c(result, threshold=threshold)

  ##exit; return
  return (result)

}


#########################################################################
##
##  PRIMARY RECURSION
##
#########################################################################

### Recursive function to determine first order maximal distances.
rsfParseTree <- function(recursiveObject,
                         max.order,
                         sub.order,
                         nativeArray,
                         b,
                         distance,
                         subtreeFlag) {
                          

  ## Diagnostic count of calls.
  recursiveObject$diagnostic <- recursiveObject$diagnostic + 1

  ## Weak consistency check to ensure that the iteration matches the treeID in the nativeArray record.
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }

  ## Check if we are at the root node.  Otherwise increment the node at depth count.
  if (distance > 0) {
    ## We are not at the root node.
    if (distance <= length(recursiveObject$subtree$nodesAtDepth)) {    
      if (is.na(recursiveObject$subtree$nodesAtDepth[distance])) {
        ## This is the first split at this depth, so initialize the nodes at depth count.
        recursiveObject$subtree$nodesAtDepth[distance] <- 1
      }
      else {
        ## Increment the count of nodes at this depth.
        recursiveObject$subtree$nodesAtDepth[distance] <- recursiveObject$subtree$nodesAtDepth[distance] + 1
      }
    }
  }
  
    

  ## Read the current nativeArray split parameter.
  splitParameter <- nativeArray$parmID[recursiveObject$offset]

  ## Determine whether this is a terminal node.  
  if (splitParameter == 0) {
    terminalFlag <- TRUE
  }
  else if (splitParameter != 0) {
    terminalFlag <- FALSE
  }

  if (!terminalFlag) {

    ## Update the maximal subtree information for the parameter if it
    ## has not already been encountered in the tree.

    if (subtreeFlag[splitParameter] == FALSE) {


      ## Increment the subtree count for the parameter.
      recursiveObject$subtree$count[splitParameter] <- recursiveObject$subtree$count[splitParameter] + 1


      ## Update the (pre) mean distance.
      if (is.na(recursiveObject$subtree$meanSum[splitParameter])) {
        recursiveObject$subtree$meanSum[splitParameter] <- distance
      }
      else {
        recursiveObject$subtree$meanSum[splitParameter] <- recursiveObject$subtree$meanSum[splitParameter] + distance
      }

      if (max.order > 0) {
        ## Update distance encountered for the order statistic.

        ## The last element, NA, in the temporary vector is a filler to handle which.max() silliness. 
        orderVector <- c(recursiveObject$subtree$order[splitParameter, ], distance, NA)
        ## Identify the first NA occurence.
        index <- which.max(is.na(orderVector))

        ## The first (relevant) NA is filled with the current distance.
        orderVector[index] <- distance
        ## The vector is then sorted.
        sortedVector <- sort(orderVector[1:index])
        ## Restore the length of the order vector, if necessary, by apppending NA's.
        if (index <= max.order) {
          orderVector <- c(sortedVector, rep(NA, max.order-index))
        }
        else {
          orderVector <- sortedVector[1:max.order]
        }

        ## Update the order statistic.
        recursiveObject$subtree$order[splitParameter, ] <- orderVector

      }
      else {

        ## Update the minimum distance encountered for the parameter.
        if (is.na(recursiveObject$subtree$order[splitParameter])) {
          recursiveObject$subtree$order[splitParameter] <- distance
        }
        else {
          recursiveObject$subtree$order[splitParameter] <- min(recursiveObject$order[splitParameter], distance)
        }
      }

      ## Indicate that the parameter has been split on.
      subtreeFlag[splitParameter] <- TRUE

      if (sub.order == TRUE) {
        
        ## Update the diagonal element with the minimum distance encountered for the parameter.
        if (is.na(recursiveObject$subtree$subOrderDiag[splitParameter])) {
          recursiveObject$subtree$subOrderDiag[splitParameter] <- distance
        }
        else {
          recursiveObject$subtree$subOrderDiag[splitParameter] <- min(recursiveObject$subtree$subOrderDiag[splitParameter], distance)
        }

        ## Create the recursive2 list object.  Note that the relative depth starts at zero (0) for the subtree.  Also note
        ## that the object returns a vector containing the unnormalized minimum maximal w-subtree distance.  
        recursive2Object <- list(offset = recursiveObject$offset,
          depth = 0,
          minimumVector = rep(NA, dim(recursiveObject$subtree$subOrder)[2]),
          diagnostic = recursiveObject$diagnostic2)

        ## Initialize the split flags for the sub-recursion.  Ignore
        ## the diagonal element by setting the split flag.
        subtree2Flag <- rep(FALSE, dim(recursiveObject$subtree$subOrder)[2])        
        subtree2Flag[splitParameter] <- TRUE
        
        ## Do the sub-recursion.
        recursive2Object <- rsfParse2Tree(recursive2Object,
          nativeArray,
          b,
          distance=0,
          subtreeFlag=subtree2Flag)
      
        recursiveObject$diagnostic2 <- recursiveObject$diagnostic2 + recursive2Object$diagnostic

        ## Set the element corresponding to the diagonal element temporarily to the penalized value, to avoid being manipulated.
        recursive2Object$minimumVector[splitParameter] <- recursive2Object$depth

        ## Determine which w-subtrees in this v-subtree do not exist and set them to the penalized value.
        recursive2Object$minimumVector[which(is.na(recursive2Object$minimumVector))] <- recursive2Object$depth
        
        ## Normalize the relative distances by the maximum terminal node depth encountered for this v-subtree.
        recursive2Object$minimumVector <- recursive2Object$minimumVector / recursive2Object$depth

        ## Update the minimum w-subtree distance encountered in the tree.
        recursiveObject$subtree$subOrder[splitParameter, ] <- pmin(recursiveObject$subtree$subOrder[splitParameter, ], recursive2Object$minimumVector)
        
      }  ## if (sub.order == TRUE) ...
      
    }  ## if (subtreeFlag[splitParameter] == FALSE) ...

  }  ## if (splitParameter != 0) ...


  ## Update the maximum depth encountered for this tree.
  recursiveObject$subtree$depth <- max(recursiveObject$subtree$depth, distance)
    
  ## Increment the offset.
  recursiveObject$offset <- recursiveObject$offset + 1


  ## Parse left and then right, if this is not a terminal node.
  if (terminalFlag == FALSE) {

    ## Increment the (parsed) tree distance. 
    distance <- distance + 1

    ## Parse left:
    recursiveObject <- rsfParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)

    ## Parse right:
    recursiveObject <- rsfParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)

  }
  else {
    ## Update the terminal node depth sum. 
    recursiveObject$subtree$terminalDepthSum <- recursiveObject$subtree$terminalDepthSum + distance
  }
    

  return (recursiveObject)

}


#########################################################################
##
##  SECONDARY RECURSION
##
#########################################################################
  
## Determine the second order maximal subtree distance.
## Recursive Object Definition:
##
## recursiveObject =
##   list(offset,
##        depth = 0,
##        minimumVector = rep(NA, numParm),
##        diagnostic)
##
## Inputs:
##   distance
##     - distance (depth) from the root of the (primary)
##       v-subtree.  Thus, this always starts at zero (0).
##     
##
## Outputs:  
##   recursiveObject$minimumDistance
##     - vector of length [numParm]
##     - second order unnormalized minimum maximum w-subtree distances
##       relative to the subtree defined by the initial root call.
##   recursiveObject$depth
##     - maximum terminal node depth, relative to the v-subtree.
##
## Algorithm Notes:  
##
## We use the terminology (v,w) for a w-maximal subtree within
## (relative to) a v-maximal subtree.
## 
## STEP-1:
## 
## For a given v, and a given v-max subtree, we search for all w-max
## subtrees within this v-max subtree, where v != w.  We single out
## the w-subtree with minimum depth -- that is, the w-split closest to
## the v-split.  Thus, when we talk about depths of w, these are
## relative depths, relative to the v-max subtree.  So the depth at
## root of the v-max subtree is defined to be zero (0).
## 
## For a given v, and a given v-max subtree, we also determine the
## maximum terminal node depth.  Let's call this the penalty depth.
## Again, this is relative to the v-max subtree.
## 
## If a particular w-max subtree does not exist for a given v and
## v-max subtree, it is assigned the penalty depth.
## 
## Now, there are p-1 minimum w-max subtree depths, given by analyzing
## all w != v.  These distances are normalized (divided by) by the
## penalty depth for the v-max subtree, so that all p-1 values are in
## (0,1].
## 
## 
## STEP-2:
## Consider v as fixed.  We do STEP-1 for all v-max subtrees.  Say
## there are count(v) of these.  We then take the MINIMUM normalized depth across the count(v)
## vectors of length p-1 that give the minimum w-max subtree depths by
## dividing by count(v).
## 
## STEP-3:
## Currently, if the v-max subtree does not exist, the row
## corresponding to v in the $subOrder matrix, namely [v, ] will be
## one (1) except for the element corresponding to the diagonal [v,v].
## It will be zero (0).
## 
## STEP-4:
## 
## All of the above are repeated for each tree in the forest such that
## the [p]x[p] matrix for each tree is normalized (divided by) the the
## number of valid (non-stumped) trees in the forest.
## 
## There is the rare possibility that a predictor v will never be
## split upon within the forest.  The entire row for this predictor,
## given by [v, ], will be NA.

rsfParse2Tree <- function(recursiveObject,
                          nativeArray,
                          b,
                          distance,
                          subtreeFlag) {

  ## Diagnostic count of calls.
  recursiveObject$diagnostic = recursiveObject$diagnostic + 1

  ## Weak consistency check to ensure that the iteration matches the treeID in the nativeArray record.
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  
  ## Read the current nativeArray split parameter.
  splitParameter = nativeArray$parmID[recursiveObject$offset]

  ## Determine whether this is a terminal node.  
  if (splitParameter == 0) {
    terminalFlag = TRUE
  }
  else if (splitParameter != 0) {
    terminalFlag = FALSE
  }

  if (splitParameter != 0) {

    ## Check that the parameter has not already been encountered along this path.
    if (subtreeFlag[splitParameter] == FALSE) {

      ## Update the minimum maximal subtree distance encountered for the parameter.
      if (is.na(recursiveObject$minimumVector[splitParameter])) {
        recursiveObject$minimumVector[splitParameter] = distance
      }
      else {
        recursiveObject$minimumVector[splitParameter] = min(recursiveObject$minimumVector[splitParameter], distance)
      }
     
      ## Indicate that the parameter has been split on
      subtreeFlag[splitParameter] = TRUE
      
    }


  }

  ## Update the (relative) maximum depth encountered for this tree.
  recursiveObject$depth = max(recursiveObject$depth, distance)
    
  ## Increment the (parsed) tree distance. 
  distance = distance + 1

  ## Increment the offset.
  recursiveObject$offset = recursiveObject$offset + 1
  
  ## Parse left and then right, if this is not a terminal node.
  if (terminalFlag == FALSE) {

    ## Parse left:
    recursiveObject = rsfParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)

    ## Parse right:
    recursiveObject = rsfParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)

  }

  return (recursiveObject)

}



#########################################################################
##
##  Supporting Functions
##
#########################################################################

## Minimal depth threshold functions

  maxDepthProb <- function(p, D, l) {
    # density for minimal depth
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
    # compute mean and std of minimal depth under null
    # mean currently used only
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
    #n.at.d = forest averaged node counts at depth d (rounded)
    #avg.d  = forest averaged tree depth (rounded)
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

