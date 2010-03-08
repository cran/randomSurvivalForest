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

pmml2rsf <- function(pmmlRoot, ...) {

  ## Ensure that the root node is not null.
  if (is.null(pmmlRoot)) {
    stop("PMML content is NULL.  Please ensure the data is valid.")
  }

  ## Ensure that the root node is of type PMML.
  if (xmlName(pmmlRoot) != "PMML") {
    stop("XML content malformed or not of type PMML.  Please ensure the data is valid.")
  }
  
  ## Point to the MiningBuildTask node, there should only be one.
  rsfMiningBldTsk <- pmmlRoot[names(pmmlRoot) == "MiningBuildTask"]
  ## Stop if there is not one and only one.
  if(length(rsfMiningBldTsk) != 1) {
    stop("MiningBuildTask not found or malformed in PMML content.  Please ensure the data is valid.")
  }
  
  ## Point to the DataDictionary node, there should only be one.
  rsfDataDict <- pmmlRoot[names(pmmlRoot) == "DataDictionary"]
  ## Stop if there is not one and only one.
  if(length(rsfDataDict) != 1) {
    stop("DataDictionary not found or malformed in PMML content.  Please ensure the data is valid.")
  }

  ## Point to the TreeModel nodes, there will be one for each tree.
  rsfTreeModel <- pmmlRoot[names(pmmlRoot) == "TreeModel"]

  ## Count the number of trees in the forest.
  numTrees = length(rsfTreeModel)

  ## Stop if there are no trees in the forest.
  if(numTrees < 1) {
    stop("TreeModel not found in PMML content.  Please ensure the data is valid.")
  }

  ## Extract the MiningBuildTask node.
  miningBldTskNode = rsfMiningBldTsk[1]$MiningBuildTask
  
  ## Extract the DataDictionary node.
  dataDictNode = rsfDataDict[1]$DataDictionary

  ## Extract the numberOfFields attribute.
  numFields = as.integer(xmlGetAttr(dataDictNode, "numberOfFields"))

  ## Point to the DataField nodes, there will be one for each predictor.
  rsfDataField = dataDictNode[names(dataDictNode) == "DataField" ]

  ## Count the number of predictors.
  numPredictors = length(rsfDataField)
  ## Stop if there are no predictors.
  if(numPredictors < 1) {
    stop("DataField not found in PMML content.  Please ensure the data is valid.")
  }

  ## Check that the number of fields attribute is consistent with the actual number of predictors listed.
  if(numPredictors != numFields) {
    stop("DataField count inconsistent in PMML content.  Please ensure the data is valid.")
  }

  ## Create an empty vector of strings for the predictor names.
  predictorNames = rep(NA, numPredictors)

  ##  Create a default vector of strings for the predictor types.
  predictorType = rep("R", numPredictors)

  ## Extract the predictor names and predictor types.
  for (i in 1:numPredictors) {
    dataFieldNode = rsfDataField[i]$DataField
    predictorNames[i] = xmlGetAttr(dataFieldNode, "name")
    if (xmlGetAttr(dataFieldNode, "optype") == "categorical") {
      predictorType[i] = "C"
    }
    if (xmlGetAttr(dataFieldNode, "optype") == "ordinal") {
      predictorType[i] = "I"
    }
  }

  ## Point to all Extension nodes, we need only the one containing the RSF extensions.
  rsfMiningBldTskExtn = miningBldTskNode[names(miningBldTskNode) == "Extension"]

  ## Count the number of Extension nodes.
  numExtn = length(rsfMiningBldTskExtn)

  extnFlag = FALSE
  ## Loop through all Extension nodes, searching for the RSF extensions.
  for (i in 1:numExtn) {
    if (extnFlag == FALSE) {
      ## Determine if this extension is the RSF extension.
      mbteNode = rsfMiningBldTskExtn[i]$Extension

      rsfFormula = mbteNode[names(mbteNode) == "X-RSF-Formula"]            
      rsfForestSeed = mbteNode[names(mbteNode) == "X-RSF-ForestSeed"]
      rsfTimeInt = mbteNode[names(mbteNode) == "X-RSF-TimesOfInterest"]


      if (length(rsfFormula) == 1) {
        formulaNode = rsfFormula[1]$'X-RSF-Formula'
        extnFlag = TRUE
      }
      else {
        stop("Formula extension not found in PMML content.  Please ensure the data is valid.")    
      }

      formula = as.formula(xmlGetAttr(formulaNode, "name"))

      if (length(rsfForestSeed) == 1) {
        forestSeedNode = rsfForestSeed[1]$'X-RSF-ForestSeed'
        extnFlag = TRUE
      }
      else {
        stop("ForestSeed extension not found in PMML content.  Please ensure the data is valid.")    
      }

      forestSeed = as.integer(xmlGetAttr(forestSeedNode, "value"))
      
      if (length(rsfTimeInt) == 1) {
        timeIntNode = rsfTimeInt[1]$'X-RSF-TimesOfInterest'
        ## Point to the Array node containing the time points of interest.
        rsfTIArray = timeIntNode[names(timeIntNode) == "Array"]
        ## Ensure there is one and only one Array node.
        if (length(rsfTIArray) == 1) {
          tiaNode = rsfTIArray[1]$Array
          ## Count the number of time points of interest.
          numTimeInt = as.integer(xmlGetAttr(tiaNode, "n"))
          if (numTimeInt > 0) {
            ## Extract the raw string of numbers.
            timeIntStr = xmlValue(tiaNode)
            ## Strip spaces and new line characters, convert to double vector.
            timeInterest = as.double(strsplit(gsub(" ", "", timeIntStr), "\n")[[1]])

            ## Check that the resulting dimension of the timeInterest vector is as defined in the XML file.
            if (length(timeInterest) != numTimeInt) {
              stop("Array dimension inconsistent in PMML content.  Please ensure that the data is valid.")
            }
            extnFlag = TRUE
          }
          else {
            stop("TimesOfInterest array dimension invalid in PMML content.  Please ensure the data is valid.")    
          }
        }
        else {
          stop("TimesOfInterest array not found in PMML content.  Please ensure the data is valid.")    
        }
      }
      else {
        stop("TimesOfInterest extension not found in PMML content.  Please ensure the data is valid.")    
      }
    }
  }


  ## Create the data frame that will contain the nativeArray and nativeFactorArray
  ## records representing the forest.
  nativeArrayNames = c("treeID", "nodeID", "parmID", "contPT", "mwcpSZ")
  nativeArray = as.data.frame(rbind(rep(0, length(nativeArrayNames))))
  names(nativeArray) = nativeArrayNames

  nativeFactorArray = as.vector(0)

  ## Define the variables for the offsets and leaf count in the recursive output object.
  offset = mwcpOffset = leafCount = 1

  ## Create the recursive output object.  This would be unnecessary 
  ## if it was possible to declare global variables in a package.
  recursiveObject = list(nativeArray  = nativeArray, nativeFactorArray = nativeFactorArray, offset = offset, mwcpOffset = mwcpOffset, leafCount = leafCount)

  ## Loop through all trees in the forest and extract the data.
  for (b in 1:numTrees) {
    treeModelNode = rsfTreeModel[b]$TreeModel

    if (xmlGetAttr(treeModelNode, "splitCharacteristic") != "binary") {
      stop("TreeModel splitCharacteristic in tree ", b, " not binary or malformed in PMML content.  Please ensure the data is valid.")
    }

    ## Point to the MiningSchema node, there should only be one.
    rsfMiningSchema = treeModelNode[names(treeModelNode) == "MiningSchema" ]
    ## Stop if there is not one and only one.
    if(length(rsfMiningSchema) != 1) {
      stop("MiningSchema not found or malformed in PMML content.  Please ensure the data is valid.")
    }

    ## Extract the MiningSchema node.
    miningSchemaNode = rsfMiningSchema[1]$MiningSchema

    ## Point to the MiningField nodes, there will be one for each covariate.
    rsfMiningField  = miningSchemaNode[names(miningSchemaNode) == "MiningField" ]

    ## Count the number of covariates used in this tree.
    numCovariates = length(rsfMiningField)
    ## Stop if there are no predictors.
    if(numCovariates < 1) {
      stop("MiningField in tree ", b, " not found in PMML content.  Please ensure the data is valid.")
    }

    ## Create an empty vector of strings for the names of the covariates.
    covariates = rep(NA, numCovariates)

    ## Extract the covariate names and populate the vector of covariate names.
    for (j in 1:numCovariates) {
      miningFieldNode = rsfMiningField[j]$MiningField
      covariates[j] = xmlGetAttr(miningFieldNode, "name")
    }


    ## Consistency check to ensure that the covariates are a subset of the predictors.
    if (sum(is.element(covariates, predictorNames)) != numCovariates) {
      stop("Covariates in tree ", b, " not a subset of predictors.  Please ensure the data is valid.")
    }

    ## Point to the root node for the tree, there should be only one.
    rsfNodePtr = treeModelNode[names(treeModelNode) == "Node"]
    if(length(rsfNodePtr) != 1) {
      stop("Node in tree ", b, " not found or malformed in PMML content.  Please ensure the data is valid.")
    }
    ## Extract the root node for the tree.
    rsfNode = rsfNodePtr[1]$Node

    ## Re-initialize the leafCount for this tree.
    recursiveObject$leafCount = 1

    recursiveObject = pmmlRestoreTree(recursiveObject=recursiveObject, treeID=b, parent=rsfNode, predictorNames=predictorNames, predictorType=predictorType)
  }

  ## Delete the first and empty record in the nativeArray and nativeFactorArray.  This is an artifact of initialization.
  recursiveObject$nativeArray = recursiveObject$nativeArray[-1,]
  recursiveObject$nativeFactorArray = recursiveObject$nativeFactorArray[-1]

  outputObject = list (nativeArray = recursiveObject$nativeArray, nativeFactorArray = recursiveObject$nativeFactorArray, timeInterest = timeInterest, predictorNames = predictorNames, formula = formula, seed = forestSeed)
  class(outputObject) = "rsfForest"
  return (outputObject)
}

pmmlRestoreTree <- function (recursiveObject, treeID, parent, predictorNames, predictorType) {

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

  ## Increment the nativeArray record pointer and add an empty record to the nativeArray in preparation
  ## for the upcomming recursive calls for the left and right daughters.
  recursiveObject$offset =  recursiveObject$offset + 1

  recursiveObject$nativeArray = rbind(recursiveObject$nativeArray, rep(0, length(dim(recursiveObject$nativeArray)[1])))

  ## Note that the nativeArray record pointer points to an empty record on entry.  On exit
  ## and before each recursive call, this pointer is incremented and an empty 
  ## record is added to the nativeArray.  

  ## Initialize the tree identifier and node identifier.
  recursiveObject$nativeArray$treeID[recursiveObject$offset] = treeID
  recursiveObject$nativeArray$nodeID[recursiveObject$offset] = recursiveObject$leafCount

  ## Mark the split parameter and continuous split point as zero and NA for now.
  recursiveObject$nativeArray$parmID[recursiveObject$offset] = 0
  recursiveObject$nativeArray$contPT[recursiveObject$offset] = NA
  ## Set the multi-word complementary pair size to zero for now.
  recursiveObject$nativeArray$mwcpSZ[recursiveObject$offset] = 0

  ## Determine whether this is a terminal node by testing for the presence of daughter nodes.
  ## We could, optionally, test for the presense of a SimplePredicate tag, but the former
  ## method is implemented.

  daughterPtrs = parent[names(parent) == "Node"]

  if (length(daughterPtrs) == 2) {
    ## There are daughters in the node. Extract the daughter nodes.

    daughters = list (daughterPtrs[1]$Node, daughterPtrs[2]$Node)

    ## Determine which daughters are the left and right nodes.  We do not assume any
    ## ordering.  We determine left and right by the censoring information in each
    ## daughter node.  If the SimplePredicate tag has attribute lessOrEqual, it is
    ## the left daughter.  If the attribute is greaterThan, it is the right daughter.

    ## We perform the following consistency check.  We determine that attribute does
    ## confirm a binary split with a unique value as the split point.

    splitInfoPtrs = list (daughters[[1]][names(daughters[[1]]) == "SimplePredicate"], daughters[[2]][names(daughters[[2]]) == "SimplePredicate"])

    if ((length(splitInfoPtrs[[1]]) != 1) | (length(splitInfoPtrs[[2]]) != 1)) {
      stop("SimplePredicate in tree ", treeID, " not found or malformed in PMML content.  Please ensure the data is valid.")          
    }

    ## Extract the split information from both daughters and continue with the consistency check.
    splitInfo = list (splitInfoPtrs[[1]]$SimplePredicate,  splitInfoPtrs[[2]]$SimplePredicate)

    splitParm = list (xmlGetAttr(splitInfo[[1]], "field"), xmlGetAttr(splitInfo[[2]], "field"))

    if (splitParm[[1]] != splitParm[[2]]) {
      stop("SimplePredicate field attribute in tree ", treeID, " malformed in PMML content.  Please ensure the data is valid.")          
    }

    contSplitValue = list (xmlGetAttr(splitInfo[[1]], "value"), xmlGetAttr(splitInfo[[2]], "value"))

    if ((contSplitValue[[1]] != "NA") & (contSplitValue[[2]] != "NA")) {
      if (as.double(contSplitValue[[1]]) != as.double(contSplitValue[[2]])) {
        stop("SimplePredicate value attribute in tree ", treeID, " malformed in PMML content.  Please ensure the data is valid.")          
      }
      contPT = as.double(contSplitValue[[1]])
    }
    else {
      if (contSplitValue[[1]] != contSplitValue[[2]]) {
        stop("SimplePredicate value attribute in tree ", treeID, " malformed in PMML content.  Please ensure the data is valid.")          
      }
      contPT = NA
    }
    
    splitOperator = list (xmlGetAttr(splitInfo[[1]], "operator"), xmlGetAttr(splitInfo[[2]], "operator"))

    leftIndex = rightIndex = 0
    ## Ensure that the split is coherent and binary, and determine which daughter is
    ## the left and right so that the appropriate node is passed to the recursive call.
    if ((splitOperator[[1]] == "lessOrEqual") & (splitOperator[[2]] == "greaterThan")) {
      leftIndex = 1
      rightIndex = 2
    }
    else if  ( (splitOperator[[2]] == "lessOrEqual") & (splitOperator[[1]] == "greaterThan")) {
      leftIndex = 2
      rightIndex = 1
    }
    else {
      stop("SimplePredicate operator attribute in tree ", treeID, " malformed in PMML content.  Please ensure the data is valid.")          
    }

    ## Ensure that the split parameter is a valid parameter.
    parmID = which(predictorNames == splitParm[[1]])
    if (length(parmID) != 1) {
      stop("SimplePredicate field attribute in tree ", treeID, " malformed in PMML content.  Please ensure the data is valid.")          
    }

    ## Initialize the split parameter and continuous split point.
    ## The split point will be NA in the case of categorical factors.
    recursiveObject$nativeArray$parmID[recursiveObject$offset] = parmID
    recursiveObject$nativeArray$contPT[recursiveObject$offset] = contPT

    ## The default multi-word complementary pair size is zero.
    recursiveObject$nativeArray$mwcpSZ[recursiveObject$offset] = 0

    ## Check for a categorical split.
    if (predictorType[parmID] == "C") {
      ## Extract the multi-word complementary pair and its size.

      ## Point to all Extension nodes in the SimplePredicate node.
      ## We need only the one containing the RSF MWCP extension.
      rsfSPExtn = splitInfo[[1]][names(splitInfo[[1]]) == "Extension"]

      ## Count the number of Extension nodes.
      numExtn = length(rsfSPExtn)

      extnFlag = FALSE
      ## Loop through all Extension nodes, searching for the RSF MWCP extension.
      for (i in 1:numExtn) {
        if (extnFlag == FALSE) {
          ## Determine if this extension is the RSF MWCP extension.
          spNode = rsfSPExtn[i]$Extension

          rsfMwcp = spNode[names(spNode) == "X-RSF-MWCP"]            

          if (length(rsfMwcp) == 1) {
            mwcpNode = rsfMwcp[1]$'X-RSF-MWCP'
            ## Point to the Array node.
            rsfMwcpArray = mwcpNode[names(mwcpNode) == "Array"]
            ## Ensure there is one and only one Array node.
            if (length(rsfMwcpArray) == 1) {
              mwcpNode = rsfMwcpArray[1]$Array
              ## Count the number of words in the multi-word complementary pair.
              numMwcp = as.integer(xmlGetAttr(mwcpNode, "n"))

              if (numMwcp > 0) {
                ## Extract the raw string of numbers.
                mwcpStr = xmlValue(mwcpNode)
                ## Strip spaces and new line characters, convert to an integer vector.
                mwcpStrArray = strsplit(gsub(" ", "", mwcpStr), "\n")[[1]]
                ## Note that words with the value 0x80000000 are interpreted as NA in R.
                mwcpNA = (mwcpStrArray == "NA")
                mwcpPT = rep(NA, length(mwcpStrArray))
                mwcpPT[!mwcpNA] = as.integer(mwcpStrArray[!mwcpNA]) 
                
                ## Check that the resulting dimension of the vector is as defined in the XML file.
                if (length(mwcpPT) != numMwcp) {
                  stop("Array dimension inconsistent in PMML content.  Please ensure that the data is valid.")
                }

                ## Increment the nativeFactorArray record pointer and add the MWCP to the nativeFactorArray.
                recursiveObject$mwcpOffset =  recursiveObject$mwcpOffset + 1
                recursiveObject$nativeFactorArray = c(recursiveObject$nativeFactorArray, mwcpPT)

                ## Save the size of the multi-word complementary pair.
                recursiveObject$nativeArray$mwcpSZ[recursiveObject$offset] = numMwcp
                extnFlag = TRUE
              }
              else {
                stop("MWCP array dimension invalid in PMML content.  Please ensure the data is valid.")    
              }
            }
            else {
              stop("MWCP array not found in PMML content.  Please ensure the data is valid.")    
            }
          }
          else {
            stop("MWCP extension not found in PMML content.  Please ensure the data is valid.")    
          }
        }
      }

    }

    ## The left daughter will assume the parent's leaf count.  The right daughter will
    ## increment the parent's leaf count by one.  This information is encoded in the
    ## PMML Node identifier, but will be assumed redundant and NOT checked, since forests
    ## produced  by other applications may not encode this information in a similar manner.
    ## In addition, the PMML Node score encodes whether or not the node is a terminal.  
    ## This information will also be considered redundant.

    recursiveObject = pmmlRestoreTree(recursiveObject, treeID, daughters[[leftIndex]], predictorNames, predictorType)

    ## Increment the leaf count for the right daughter only.
    recursiveObject$leafCount =  recursiveObject$leafCount + 1

    recursiveObject = pmmlRestoreTree(recursiveObject, treeID, daughters[[rightIndex]], predictorNames, predictorType)
  }
  else {
    ## This is a terminal node.  Recurse, recurse, recurse.
  }

  return (recursiveObject)
}

