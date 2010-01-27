####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.1
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

rsf2pmml <- function(object, ...) {

  ## Ensure the forest object is coherent.
  rsfForest <- checkForestObject(object)

  ## Initialize the local variables extracted from the forest object.
  nativeArray <- rsfForest$nativeArray
  timeInterest <- rsfForest$timeInterest
  formula <- rsfForest$formula
  forestSeed <- rsfForest$seed
  predictorNames <- rsfForest$predictorNames

  ## Extract the predictor types.
  get.factor <- extract.factor(rsfForest$predictors, predictorNames)
  predictorTypes <- get.factor$predictorType

  ## This may be null in the absence of factors.
  nativeFactorArray <- rsfForest$nativeFactorArray
        
  ## Count the number of trees in the forest.
  numTrees <- length(as.vector(unique(nativeArray$treeID)))

  ## Define the root elements of the PMML file.  This is a quick work-around for an issue
  ## with this version of the XML package, and the inablility to add namespace information
  ## and attributes concurrently. 
  rootString <- getRootString()

  ## Define the document and the root node.
  pmmlDoc <- xmlTreeParse(rootString, asText=TRUE)
  pmmlRoot <- xmlRoot(pmmlDoc)

  ## Define the MiningBuildTask node for the document.
  miningBldTskNode <- xmlNode("MiningBuildTask")
  
  ## Define the Extension node for the document.
  mbtExtentionNode <- xmlNode("Extension")

  ## Add the formula to the Extension node.
  mbtExtentionNode <- append.XMLNode(mbtExtentionNode, xmlNode("X-RSF-Formula", attrs=c(name=formula)))

  ## Add the bootstrap seed to the Extension node.
  mbtExtentionNode <- append.XMLNode(mbtExtentionNode, xmlNode("X-RSF-ForestSeed", attrs=c(value=forestSeed)))

  ## Add the times of interest to the Extension node.
  mbtExtentionNode <- append.XMLNode(mbtExtentionNode, 
    xmlNode("X-RSF-TimesOfInterest", 
            xmlNode("Array", 
                    attrs=c(type="real", n=length(timeInterest)), 
                    paste(timeInterest, collapse="  \n  "))))
  
  ## Add the Extension node to the MiningBuildTask.
  miningBldTskNode <- append.XMLNode(miningBldTskNode, mbtExtentionNode)

  ## Add the MiningBuildTask to the root node.
  pmmlRoot <- append.XMLNode(pmmlRoot, miningBldTskNode)
  
  ## Add the DataDictionary to the root node.
  pmmlRoot <- append.XMLNode(pmmlRoot, getDataDictNode(predictorNames=predictorNames, predictorTypes=predictorTypes))

  ## Create a dummy XML node object to insert into the recursive output object.
  internalNode <- xmlNode("Null")

  ## Define the variables for the offsets and leaf count in the recursive output object.
  offset <- mwcpOffset <- leafCount <- 1

  ## Create the recursive output object.  This would be unnecessary 
  ## if it was possible to declare global variables in a package.
  recursiveObject <- list(internalNode = internalNode, offset = offset, mwcpOffset = mwcpOffset, leafCount = leafCount)

  ## Loop through all trees in the forest and extract the data.
  for (b in 1:numTrees) {
    treeModelNode <- xmlNode("TreeModel", attrs=c(modelName=b, functionName="prediction", algorithmName="rsf", splitCharacteristic="binary"))

    miningSchemaNode <- xmlNode("MiningSchema")
    
    for (k in 1:length(predictorNames)) {
      miningSchemaNode <- append.XMLNode(miningSchemaNode, xmlNode("MiningField", attrs=c(name=predictorNames[k])))
    }

    treeModelNode <- append.XMLNode(treeModelNode, miningSchemaNode)

    ## Global dependencies:  (predictorNames, forest)

    ## Initialize the root node.  This differs from the rest of the internal nodes in the PMML structure.
    treeRoot <- xmlNode("Node", attrs=c(score=0, id=1))
    treeRoot <- append.XMLNode(treeRoot, xmlNode("True"))
    
    rootParmID <- nativeArray$parmID[recursiveObject$offset] 
    rootContPT <- nativeArray$contPT[recursiveObject$offset]
    rootMwcpSZ <- nativeArray$mwcpSZ[recursiveObject$offset]
    
    ## Initialize the multi-word complementary pair if necessary.
    if (rootMwcpSZ == 0) {
      rootMwcpPT <- NULL
    }
    else {
      rootMwcpPT <- nativeFactorArray[recursiveObject$mwcpOffset:(recursiveObject$mwcpOffset + rootMwcpSZ - 1)]
      recursiveObject$mwcpOffset <- recursiveObject$mwcpOffset + rootMwcpSZ 
    }

    recursiveObject$offset <- recursiveObject$offset + 1
    recursiveObject$leafCount <- 1
    
    ## Check that the current tree is not a stump (root node only with no branches)
    if (rootParmID != 0) {

      ## The tree must be created in two phases.  First, the root left daughter branches
      ## are created.  Second, the root right daughter branches are created.  This is
      ## due to the root node having a slightly different structure using the PMML
      ## protocol.  The root node actually has no split information.  The split information
      ## is encoded into the daughter nodes.  Thus, instead of making a check for the root
      ## node in the recursive routine, we call the recursive routine twice.

      ## Create the left daughter nodes.  Note that the object node content is irrelevant as input.
      recursiveObject$internalNode <- NULL
      recursiveObject <- rsfMakeTree(recursiveObject, nativeArray, nativeFactorArray, predictorNames, predictorTypes, b, -1, rootParmID, rootContPT, rootMwcpSZ, rootMwcpPT)

      treeRoot <- append.XMLNode(treeRoot, recursiveObject$internalNode)

      recursiveObject$leafCount <- recursiveObject$leafCount + 1

      ## Creat the right daughter nodes.  Note that the object node content is irrelevant as input.
      recursiveObject$internalNode <- NULL
      recursiveObject <- rsfMakeTree(recursiveObject, nativeArray, nativeFactorArray, predictorNames, predictorTypes, b, +1, rootParmID, rootContPT, rootMwcpSZ, rootMwcpPT)

      treeRoot <- append.XMLNode(treeRoot, recursiveObject$internalNode)

    }
    
    ## Add the current tree to the PMML data structure.
    treeModelNode <- append.XMLNode(treeModelNode, treeRoot)
    pmmlRoot <- append.XMLNode(pmmlRoot, treeModelNode)
  }

  return (pmmlRoot)
}

rsfMakeTree <- function(recursiveObject,
                        nativeArray,
                        nativeFactorArray,
                        predictorNames,
                        predictorTypes,
                        b,
                        daughter,
                        splitParameter,
                        contValue,
                        mwcpSZ,
                        mwcpPT) {

  ## Node information encoded in a PMML TreeModel follows a slightly different protocol
  ## than that encoded in our RSF matrix representation.   The RSF representation
  ## is linear in nature   Each record containing node information must encode the split
  ## information, particularly the split parameter and split point, in the record itself.
  ## In contrast, the PMML TreeModel indicates a split by the presence of daughters in 
  ## the node.  The split parameter and split point are encoded by a SimplePredicate tag
  ## in the daughters.  Thus every node except the root node contains the SimplePredicate tag.
  ## In creating a PMML tree from an RSF tree, the recursive algorithm
  ## requires a "look back" to the previous record in the RSF tree to determine the split
  ## parameter and value.  This is accomplished via the parameters passed by the parent 
  ## call to this routine. 

  ## Weak consistency check to ensure that the iteration matches the treeID in the nativeArray record.
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }

  ## Read the current nativeArray and nativeFactorArray records, and save the "look back" information.
  fwdSplitParameter <- nativeArray$parmID[recursiveObject$offset]
  fwdContValue <- nativeArray$contPT[recursiveObject$offset]
  fwdMwcpSZ <-  nativeArray$mwcpSZ[recursiveObject$offset]
  fwdMwcpOffset <- recursiveObject$mwcpOffset
  
  ## Initialize the multi-word complementary pair if necessary.
  if (fwdMwcpSZ == 0) {
    fwdMwcpPT <- NULL
  }
  else {
    ## Note that words with the value 0x80000000 are interpreted as NA in R.
    fwdMwcpPT <- nativeFactorArray[fwdMwcpOffset:(fwdMwcpOffset + fwdMwcpSZ - 1)]
    recursiveObject$mwcpOffset <- recursiveObject$mwcpOffset + fwdMwcpSZ 
  }

  ## Determine whether this is a terminal node.  Create the node that will be returned on this call.  
  if (fwdSplitParameter == 0) {
    rsfNode <- xmlNode("Node", attrs=c(score=1, id=nativeArray$nodeID[recursiveObject$offset])) 
    terminalFlag <- TRUE
  }
  else if (fwdSplitParameter != 0) {
    rsfNode <- xmlNode("Node", attrs=c(score=0, id=nativeArray$nodeID[recursiveObject$offset])) 
    terminalFlag <- FALSE
  }
  
  ## Determine whether this the left of right daughter.
  if (daughter == -1) {
    parseString <- "lessOrEqual"
  }
  else if (daughter == +1) {
    parseString <- "greaterThan"
  }
  else {
    ## Ensure that the function call is coherent to aid in debugging.
    stop("Invalid parse direction encountered during recursion.  Please contact Technical Support.")
  }
  
  ## Define the split information to this node via the look back.
  simplePredicateNode <- xmlNode("SimplePredicate", attrs=c(field=predictorNames[splitParameter], operator=parseString, value=contValue))

  ## In the case of factors, add the multi-word complementary pair information.
  
  if (predictorTypes[splitParameter] == "C") {

    ## Note that words with the value 0x80000000 are interpreted as NA in R.
    mwcpNode <- xmlNode("Extension",
      xmlNode("X-RSF-MWCP", 
              xmlNode("Array", 
                      attrs=c(type="integer",n=mwcpSZ), 
                      paste(mwcpPT, collapse="  \n  "))))

    simplePredicateNode <- append.XMLNode(simplePredicateNode, mwcpNode)

  }

  rsfNode <- append.XMLNode(rsfNode, simplePredicateNode)

  ## Increment the offset, always.
  recursiveObject$offset <- recursiveObject$offset + 1

  ## Parse left and then right, if this is not a terminal node.
  if (terminalFlag == FALSE) {

    ## Parse left:
    ## Do not increment the leafCount.  Internally increment the offset, always.
    ## Note that the object node content is irrelevant as input.
    recursiveObject$internalNode <- NULL
    recursiveObject <- rsfMakeTree(recursiveObject, nativeArray, nativeFactorArray, predictorNames, predictorTypes, b, daughter = -1, fwdSplitParameter, fwdContValue, fwdMwcpSZ, fwdMwcpPT)

    rsfNode <- append.XMLNode(rsfNode, recursiveObject$internalNode)

    ## Parse right:
    ## Increment the leafCount.  Internally increment the offset, always.
    ## Note that the object node content is irrelevant as input.
    recursiveObject$leafCount <- recursiveObject$leafCount + 1
    recursiveObject$internalNode <- NULL
    recursiveObject <- rsfMakeTree(recursiveObject, nativeArray, nativeFactorArray, predictorNames, predictorTypes, b, daughter = +1, fwdSplitParameter, fwdContValue, fwdMwcpSZ, fwdMwcpPT)

    rsfNode <- append.XMLNode(rsfNode, recursiveObject$internalNode)

  }

  ## Modify the recursive object with the new internal node structure.
  recursiveObject$internalNode <- rsfNode

  return (recursiveObject)

}

## ## ## ##########################
## ## ## ## Internal functions:  ##
## ## ## ##########################

## Coherency check of the forest object.  Failure stops execution, hence there is no return value.
checkForestObject <- function(object) {

  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'")
  }

  if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) {
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    }
    rsfForest <- object$forest
  }
  else {
    rsfForest <- object
  }

  if (is.null(rsfForest$nativeArray)) {
    stop("RSF nativeArray content is NULL.  Please ensure the object is valid.")
  }
  
  if (is.null(rsfForest$timeInterest)) {
    stop("RSF timeInterest content is NULL.  Please ensure the object is valid.")
  }

  if (is.null(rsfForest$formula)) {
    stop("RSF formula content is NULL.  Please ensure the object is valid.")
  }

  if (is.null(rsfForest$seed)) {
    stop("RSF forestSeed content is NULL.  Please ensure the object is valid.")
  }

  if (is.null(rsfForest$predictorNames)) {
    stop("RSF predictorNames content is NULL.  Please ensure the object is valid.")
  }

  if (is.null(rsfForest$predictors)) {
    stop("RSF predictors content is NULL.  Please ensure the object is valid.")
  }

  return (rsfForest)
}


##  The root string is accessed here by rsf2pmml() and by rsf2rfz().
getRootString <- function() {
  rootString <- 
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
         <PMML version=\"3.1\" xmlns=\"http://www.dmg.org/PMML-3_1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">
           <Header copyright=\"Copyright 2008, Cleveland Clinic\" description=\"Random Survival Forest Tree Model\">
              <Application name=\"Random Survival Forest\" version=\"3.0\"/>
           </Header>
         </PMML>
       "
  return (rootString)
}

## Form the data dictionary node.  This is also accessed by rsf2pmml() and by rsf2rfz().
getDataDictNode <-  function(predictorNames, predictorTypes) {

  ## Define the DataDictionary node for the document.
  dataDictNode <- xmlNode("DataDictionary", attrs=c(numberOfFields=length(predictorNames)))

  ## Add the predictor names to the DataDictionary.
  for (k in 1:length(predictorNames)) {
    if (predictorTypes[k] == "C") {
      dataDictNode <- append.XMLNode(dataDictNode, xmlNode("DataField", attrs=c(name=predictorNames[k], optype="categorical", dataType="string")))
    }
    if (predictorTypes[k] == "I") {
      dataDictNode <- append.XMLNode(dataDictNode, xmlNode("DataField", attrs=c(name=predictorNames[k], optype="ordinal", dataType="integer")))
    }
    if (predictorTypes[k] == "R") {
      dataDictNode <- append.XMLNode(dataDictNode, xmlNode("DataField", attrs=c(name=predictorNames[k], optype="continuous", dataType="double")))
    }
  }

  return (dataDictNode)
  
}
