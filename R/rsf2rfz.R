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

rsf2rfz <- function(object, forestName = NULL, ...) {

  ## Ensure the forest object is coherent.
  rsfForest <- checkForestObject(object)

  if (is.null(forestName)) {
    stop("RSF forest name is NULL.  Please provide a valid name for the forest .rfz file.")
  }

  ## If the user has already provided the .rfz extension, remove it.
  if (nchar(forestName) > 4) {
    if (substr(forestName, nchar(forestName)-3, nchar(forestName)) == ".rfz") {
      forestName <- substr(forestName, 1, nchar(forestName)-4)
    }
  }
  
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

  ## Add the DataDictionary to the root node.
  pmmlRoot <- append.XMLNode(pmmlRoot, getDataDictNode(predictorNames=predictorNames, predictorTypes=predictorTypes))

  ## Write the native array information.
  write.table(nativeArray, 
              paste(forestName, ".txt", sep=""), quote = FALSE)

  ## Write the native factor array information if it exists.
  ## *** WARNING ***  *** WARNING ***  *** WARNING *** 
  ## In 32-bit and 64-bit systems, the integer value 0x8000000 is
  ## interpreted as NA by R as it output from the native code SEXP
  ## object.  Thus write.table will contain NA's.  These need to be
  ## handled specially in the JUNG code that parses the .rfz file.
  ## *** WARNING ***  *** WARNING ***  *** WARNING ***   
  write.table(nativeFactorArray, 
                paste(forestName, ".factor.txt", sep=""), col.names=FALSE, quote = FALSE)


  ## Write the predictor names and types.
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


