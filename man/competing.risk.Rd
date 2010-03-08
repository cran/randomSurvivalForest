\name{competing.risk}
\alias{competing.risk}
\alias{competing.risk}
\title{Summary values for Competing Risks}
\description{
  Plot the ensemble cumulative incidence function (CIF), ensemble
  subsurvival function, ensemble conditional survival function, and
  ensemble unconditional survival function from a RSF competing risk
  analysis.
}
\usage{
    competing.risk(x, plot = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{plot}{Logical.  If TRUE, proximity is plotted.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  From top to bottom, left to right are plots of the: (a) ensemble CIF;
  (b) ensemble subsurvival function; (c) ensemble conditional survival
  functions; and (d) ensemble (unconditional) survival function.
  See Gray (1988) for motivation regarding the CIF and subsurvival
  functions.
  
  For right-censored data, only the survival function is plotted.

  Whenever possible, out-of-bag (OOB) values are plotted.
}
\value{
   Invisibly, the ensemble CIF, ensemble subsurvival function and
   conditional mortality for each event type.  Whenever possible,
   out-of-bag (OOB) values are returned.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
    Gray R.J. (1988).  A class of k-sample tests for comparing the
    cumulative incidence of a competing risk, \emph{Ann. Statist.},
    16: 1141-1154.

    Ishwaran H., Kogalur U.B. and Lau B.M. (2010).  Random survival
    forests for competing risks.
}
\seealso{
  \code{find.interaction},
  \code{predict.rsf},
  \code{rsf},
  \code{wihs}.
}
\examples{
\dontrun{
  data(follic, package = "randomSurvivalForest")
  follic.out <- rsf(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
  competing.risk(follic.out)
}
}
\keyword{file}
