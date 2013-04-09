\name{rsf2rfz}
\alias{rsf2rfz}
\title{Save Random Survival Forest as a .rfz Compressed File}
\description{
    \code{rsf2rfz} saves a \pkg{randomSurvivalForest} forest object as an \code{.rfz}
    compressed file that is readable by the RSF Java plugin that is
    capable of visualizing the trees in the forest.
}
\usage{
    rsf2rfz(object, forestName = NULL, ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
	forest)}.}
    \item{forestName}{The desired prefix name for forest as a string.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  An \code{.rfz} compressed file is actually a \code{.zip} file
  consisting of three files. The first is an ASCII file of type
  \code{.txt} containing the \code{$nativeArray} component of the
  forest.  The second is an ASCII file of type \code{.factor.txt}
  containing the \code{$nativefactorArray} component of the forest.  The
  third is an ASCII file of type \code{.xml} containing the \code{PMML
  DataDictionary} component.

  PMML or the Predictive Model Markup Language is an XML based language
  which provides a way for applications to define statistical and data
  mining models and to share models between PMML compliant applications.
  More information about PMML and the Data Mining Group can be found at
  http://www.dmg.org.

  This function, \code{rsf2rfz}, is used to import the geometry of the
  forest to a Java plugin that is capable of visualizing the trees in
  the forest.

  The geometry of the forest is saved as a file called
  \code{forestName.rfz} in the users working directory.  This file can
  then be read by the \pkg{randomSurvivalForest} Java plugin.

  Contact the authors on downloading the Java plugin.
}    
\value{
  None.
}
\note{
  Contact the authors on downloading the Java plugin.
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    http://www.dmg.org
}	

\seealso{
  \code{rsf}.
}
\examples{
\dontrun{
# Example 1:  Growing a forest, saving it as a \emph{.rfz} file ready
# for import into the Java plugin.


library("XML")

data(veteran, package = "randomSurvivalForest")
veteran.out <- rsf(Surv(time, status)~., data = veteran, ntree = 5)
veteran.forest <- veteran.out$forest
rsf2rfz(veteran.forest, forestName="veteran")
}
}
\keyword{survival}
\keyword{tree}
