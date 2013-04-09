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

rsf2rfz <- function(object, forestName = NULL, ...) {
  rsfForest <- checkForestObject(object)
  if (is.null(forestName)) {
    stop("RSF forest name is NULL.  Please provide a valid name for the forest .rfz file.")
  }
  if (nchar(forestName) > 4) {
    if (substr(forestName, nchar(forestName)-3, nchar(forestName)) == ".rfz") {
      forestName <- substr(forestName, 1, nchar(forestName)-4)
    }
  }
  nativeArray <- rsfForest$nativeArray
  timeInterest <- rsfForest$timeInterest
  formula <- rsfForest$formula
  forestSeed <- rsfForest$seed
  predictorNames <- rsfForest$predictorNames
  get.factor <- extract.factor(rsfForest$predictors, predictorNames)
  predictorTypes <- get.factor$predictorType
  nativeFactorArray <- rsfForest$nativeFactorArray
  numTrees <- length(as.vector(unique(nativeArray$treeID)))
  rootString <- getRootString()
  pmmlDoc <- xmlTreeParse(rootString, asText=TRUE)
  pmmlRoot <- xmlRoot(pmmlDoc)
  pmmlRoot <- append.XMLNode(pmmlRoot, getDataDictNode(predictorNames=predictorNames, predictorTypes=predictorTypes))
  write.table(nativeArray, 
              paste(forestName, ".txt", sep=""), quote = FALSE)
  write.table(nativeFactorArray, 
                paste(forestName, ".factor.txt", sep=""), col.names=FALSE, quote = FALSE)
  xmlFile <- file(paste(forestName, ".xml", sep=""), open="w")
  saveXML(pmmlRoot, xmlFile)
  close(xmlFile)
  zipCommand <- paste("zip", sep=" ",
    paste(forestName, ".rfz", sep=""),
    paste(forestName, ".txt", sep=""),
    paste(forestName, ".factor.txt", sep=""),
    paste(forestName, ".xml", sep="")) 
  system(command = zipCommand)
  unlink(paste(forestName, ".txt", sep=""))
  unlink(paste(forestName, ".factor.txt", sep=""))
  unlink(paste(forestName, ".xml", sep=""))
}
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
getDataDictNode <-  function(predictorNames, predictorTypes) {
  dataDictNode <- xmlNode("DataDictionary", attrs=c(numberOfFields=length(predictorNames)))
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
