\name{pmml_to_rsf}
\alias{pmml_to_rsf}
\title{Restore Random Survival Forest From PMML}
\description{
    \code{pmml_to_rsf} implements the Predictive Model Markup
    Language specification for a \pkg{randomSurvivalForest} forest
    object.  In particular, this function gives the user the ability to
    restore the geometry of a forest from a PMML XML document.
}
\usage{
    pmml_to_rsf(pmmlRoot, ...)
}
\arguments{
    \item{pmmlRoot}{The top-level \dQuote{XMLNode} object, or
      equivalently the root node, resulting from parsing an XML
      document.  This node must be of type PMML.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  The Predictive Model Markup Language is an XML based language which
  provides a way for applications to define statistical and data mining
  models and to share models between PMML compliant applications.  More
  information about PMML and the Data Mining Group can be found at
  http://www.dmg.org.  

  Use of PMML and \code{pmml_to_rsf} requires the \pkg{XML} package.  Be
  aware that XML is a very verbose data format.  Reasonably sized trees
  and data sets can lead to extremely large text files.  XML, while
  achieving interoperability, is not an efficient data storage mechanism
  in this case.

  It is anticipated that \code{pmml_to_rsf} will be used to import the
  geometry of a forest from other PMML compliant applications.  In
  addition, the user may wish to restore the geometry of a forest that
  was previously saved using \code{rsf_to_pmml}.
}    
\value{
     A \pkg{randomSurvivalForest} \code{forest} object.  See note below.
}
\note{
  One cautionary note is in order.  The PMML representation of the
  forest object is incomplete, in that the object needs to be massaged
  in order for prediction to be possible.  This will be clear in the
  examples.  This deficiency will be addressed in future releases of
  this package.  However, it was felt that the current functionality was
  important enough and mature enough to warrant release in this version
  of the product.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran and Udaya B. Kogalur (2006).  Random Survival
    Forests.  \emph{Cleveland Clinic Technical Report}.

    http://www.dmg.org
}	

\seealso{
  \code{xmlTreeParse},
  \code{xmlRoot},
  \code{saveXML},
  \code{rsf},
  \code{predict.rsf},
  \code{rsf_to_pmml}.
}
\examples{
  # Example 1:  Growing a forest, saving it as a PMML document,
  # restoring the forest from the PMML document, and using this forest to
  # perform prediction.

\dontrun{
  library("XML")

  data(veteran, package = "randomSurvivalForest")
  veteran.out <- rsf(Survrsf(time, status)~.,
        data = veteran,
        ntree = 5,
        forest = TRUE)
  veteran.forest <- veteran.out$forest
  veteran.pmml <- rsf_to_pmml(veteran.forest)

    # Save the document to disk.
    saveXML(veteran.pmml, "veteran.forest.xml")
    veteran.pmml <- xmlRoot(xmlTreeParse("veteran.forest.xml"))
  
  partial.forest <- pmml_to_rsf(veteran.pmml)

  # The PMML forest object must be massaged before it can be used
  # for prediction as follows:
  veteran.restored.forest <- list(
                      nativeArray=partial.forest$nativeArray, 
                      timeInterest=partial.forest$timeInterest, 
                      predictorNames=partial.forest$predictorNames,
                      bootstrapSeed=partial.forest$bootstrapSeed,
                      formula=partialForest$formula,
                      predictors=veteran.forest$predictors,
                      Time=veteran.forest$Time,
                      Cens=veteran.forest$Cens)

  # The actual time, censoring and prediction values of the data set
  # used to grow the forest are not contained in the PMML
  # representation of the forest.  If the user has access to the original
  # datafile that was used to grow the forest, this information can be
  # easily recovered.  The names corresponding to the time, censoring and
  # prediction data are all retained in the PMML representation of the forest.

  class(veteran.restored.forest) <- c("rsf", "forest")
  veteran.restored.out <- predict.rsf(veteran.restored.forest, veteran)
}

}
\keyword{survival}
\keyword{tree}
