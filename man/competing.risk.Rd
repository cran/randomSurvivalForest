\name{competing.risk}
\alias{competing.risk}
\alias{competing.risk}
\title{Summary Plots for Competing Risks}
\description{
  Plot the ensemble cumulative incidence function (CIF), ensemble
  subsurvival function, ensemble conditional survival function, and
  ensemble unconditional survival function from a random survival
  forests competing risk analysis (Ishwaran et al., 2010).
}
\usage{
    competing.risk(x, plot = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{plot}{Should curves be plotted?}
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
   OOB values are returned.
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    Gray R.J. (1988).  A class of k-sample tests for comparing the
    cumulative incidence of a competing risk, \emph{Ann. Statist.},
    16:1141-1154.

    Ishwaran H., Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2010).  Random survival forests for competing risks.
}
\seealso{
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
