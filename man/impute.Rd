\name{impute.rsf}
\alias{impute.rsf}
\title{Impute Only Mode}
\description{
  Imputation for right censored survival and competing risk data.  A
  random survival forest is grown and used to impute missing data.  No
  ensemble estimates or error rates are calculated.  This is a fast way
  to impute data.
}
\usage{
impute.rsf(formula, data = NULL, ntree = 1000, mtry = NULL,
    nodesize = NULL, splitrule = NULL, nsplit = 0, big.data = FALSE,
    nimpute = 1, predictorWt = NULL, seed = NULL, do.trace = FALSE,
    ...)
}
\arguments{
    \item{formula}{A symbolic description of the model to be fit.}
    \item{data}{Data frame containing the data to be imputed.}
    \item{ntree}{Number of trees to grow.}
    \item{mtry}{Number of variables randomly sampled at each split.}
    \item{nodesize}{Minimum terminal node size.}
    \item{splitrule}{Splitting rule used to grow trees.}
    \item{nsplit}{Non-negative integer value used to specify random splitting.}
    \item{big.data}{Set this value to TRUE for large data.}
    \item{nimpute}{Number of iterations of missing data algorithm.}
    \item{predictorWt}{Weights for selecting variables for splitting on.}
    \item{seed}{Seed for random number generator.}
    \item{do.trace}{Should trace output be enabled?} 
    \item{...}{Further arguments passed to or from other methods.}
  }
\details{
  Grows a RSF and uses this to impute missing data.  All external
  calculations such as ensemble calculations, error rates, etc. are
  turned off.  Use this function if your only interest is imputing the
  data.

  All options are the same as for \code{rsf}.
}   
\value{
  Invisibly, the data frame containing the orginal data with imputed
  data overlayed.
}
\author{
  Hemant Ishwaran \email{hemant.ishwaran@gmail.com}

  Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App. Statist.}, 2:841-860.
}
\seealso{
  \code{rsf}.
}
\examples{
\dontrun{
data(pbc, package = "randomSurvivalForest")
imputed.data <- impute.rsf(Surv(days, status) ~ ., data = pbc, nsplit = 3)
}
}
\keyword{survival}
\keyword{tree}
