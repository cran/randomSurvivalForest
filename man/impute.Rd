\name{impute.rsf}
\alias{impute.rsf}
\title{Random Survival Forest Impute Only Mode}
\description{
  Imputation mode for right censored survival and competing risk data
  using Random Survival Forests (RSF) (Ishwaran, Kogalur, Blackstone and
  Lauer, 2008).  A random forest (Breiman, 2001) of survival trees is
  grown and used to impute missing data.  No ensemble estimates or error
  rates are calculated.
}
\usage{
    impute.rsf(formula,
      data = NULL,
      ntree = 1000,
      mtry = NULL,
      nodesize = NULL,
      splitrule = NULL,
      nsplit = 0,
      big.data = FALSE,
      nimpute = 1,
      predictorWt = NULL,
      seed = NULL,
      do.trace = FALSE,
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
  Grows a RSF and uses this to impute data.  All external calculations
  such as ensemble calculations, error rates, etc. are turned off.  Use
  this function if your only interest is imputing the data.

  All options are the same as \code{rsf} and the user should consult the
  help file for \code{rsf} for details.
}   
\value{
  Invisibly, the data frame containing the orginal data with imputed
  data overlayed.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\seealso{
  \code{rsf}.
}
\examples{
\dontrun{
#------------------------------------------------------------
# Example 1:  Veteran's Administration lung cancer data
# Randomized trial of two treatment regimens for lung cancer
# See Kalbfleisch & Prentice

data(pbc, package = "randomSurvivalForest")
imputed.data <- impute.rsf(Survrsf(days, status) ~ ., data = pbc, nsplit = 3)
}
}
\keyword{survival}
\keyword{tree}
