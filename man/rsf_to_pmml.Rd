\name{rsf_to_pmml}
\alias{rsf_to_pmml}
\title{Save Random Survival Forest As PMML}
\description{
    \code{rsf_to_pmml} implements the Predictive Model Markup
    Language specification for a \pkg{randomSurvivalForest} forest
    object.  In particular, this function gives the user the ability to save the geometry of a
    forest as a PMML XML document.
}
\usage{
    rsf_to_pmml(rsfForest, ...)
}
\arguments{
    \item{rsfForest}{The \code{forest} object contained in an object of class
      \pkg{randomSurvivalForest}, as that contained in the object
      returned by the function \code{rsf} with the parameter
      \dQuote{forest=TRUE}.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  The Predictive Model Markup Language is an XML based language which
  provides a way for applications to define statistical and data mining
  models and to share models between PMML compliant applications.  More
  information about PMML and the Data Mining Group can be found at
  http://www.dmg.org.  

  Use of PMML and \code{rsf_to_pmml} requires the \pkg{XML} package.  Be
  aware that XML is a very verbose data format.  Reasonably sized trees
  and data sets can lead to extremely large text files.  XML, while
  achieving interoperability, is not an efficient data storage mechanism
  in this case.

  It is anticipated that \code{rsf_to_pmml} will be used to export the
  geometry of the forest to other PMML compliant applications, including
  graphics packages that are capable of printing binary trees.  In
  addition, the user may wish to save the geometry of the forest for later
  retrieval and prediction on new data sets using \code{rsf_to_pmml}
  together with \code{pmml_to_rsf}.
}    
\value{
  An object of class \code{XMLNode} as that defined by the \pkg{XML}
  package.  This represents the top level, or root node, of the XML
  document and is of type PMML.
}
\note{
  One cautionary note is in order.  The PMML representation of the
  \pkg{randomSurvivalForest} forest object is incomplete, in that the
  object needs to be massaged in order for prediction to be possible.
  This will be clear in the examples.  This deficiency will be addressed
  in future releases of this package.  However, it was felt that the
  current functionality was important enough and mature enough to
  warrant release in this version of the product.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran, Udaya B. Kogalur, Eugene H. Blackstone and Michael S. Lauer
    (2007).  Random Survival Forests. \emph{Cleveland Clinic Technical
    Report}.

    http://www.dmg.org
}	

\seealso{
  \code{xmlTreeParse},
  \code{xmlRoot},
  \code{saveXML},
  \code{rsf},
  \code{predict.rsf},
  \code{pmml_to_rsf}.
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
