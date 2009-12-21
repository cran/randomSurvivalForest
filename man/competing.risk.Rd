\name{competing.risk}
\alias{competing.risk}
\alias{competing.risk}
\title{Plots Summarizing Competings Risks Analysis}
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
  
  For right-censored data only the survival function is plotted. 
}
\value{
   Invisibly, the ensemble CIF and ensemble subsurvival function for
   each event type.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    R. J. Gray (1988).  A class of k-sample tests for comparing the
    cumulative incidence of a competing risk, \emph{Ann. Statist.},
    16: 1141-1154.

    H. Ishwaran, and Udaya B. Kogalur (2009).  Competing risks: a forest
    approach. \emph{Cleveland Clinic Technical Report}.
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
