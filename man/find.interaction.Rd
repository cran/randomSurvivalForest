\name{find.interaction}
\alias{find.interaction}
\title{Find Interactions Between Pairs of Variables}
\description{
  Test for pairwise interactions between variables by comparing
  pairwise importance values to additive individual importance
  values.
}
\usage{
    find.interaction(object,
                  predictorNames = NULL,
                  sorted = TRUE,
                  npred = NULL,
                  subset = NULL, 
                  nrep = 1,
                  rough = FALSE,
                  importance = c("randomsplit", "permute")[1],
                  seed = NULL,
                  do.trace = FALSE,
                  ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
      forest)}.  Note: \code{forest}=TRUE must be used in the
      original \code{rsf} call.}
    \item{predictorNames}{Character vector of variable names to be
      considered.  Default is to use all variables.}
    \item{sorted}{Should variables be sorted by importance values?  Only
      applies when \code{predictorNames}=NULL.}
    \item{npred}{Use the first npred variables as ordered by VIMP (only 
      applies when \code{predictorNames}=NULL).  Default uses all variables.}
    \item{subset}{Indices indicating which rows of the predictor matrix
      to be used (note: this applies to the \emph{object} predictor
      matrix, \code{predictors}).  Default is to use all rows.}
    \item{nrep}{Number of Monte Carlo replicates.}
    \item{rough}{Logical value indicating whether fast approximation
      should be used.  Default is FALSE.}
    \item{importance}{Method used to compute variable importance (VIMP).}
    \item{seed}{Seed for random number generator.  Must be a negative
      integer (the R wrapper handles incorrectly set seed values).}
    \item{do.trace}{Logical. Should trace output be enabled? Default is
      FALSE. Integer values can also be passed.  A positive value
      causes output to be printed each \code{do.trace} iteration.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
     Using a previously grown forest, identify pairwise interactions for
     all pairs of variables from a specified list.  Two variables are
     paired and their paired VIMP calculated (refered to as 'Paired'
     importance).  The VIMP for each separate variable is also
     calculated.  The sum of these two values is refered to as
     'Additive' importance.  A large positive or negative difference
     between 'Paired' and 'Additive' indicates an association worth
     pursuing if the VIMP's for each variable are reasonably large
     (Ishwaran, 2007).

     Depending on the size of the data, computations might be slow.  In
     such cases, users should consider setting \code{npred} to a smaller
     number, or restricting the analysis to a subset of the data.
     
     If \code{nrep} is greater than 1, the analysis is repeated
     \code{nrep} times and results averaged over the replications.

     Note that \code{find.interaction} calls the lower level function
     \code{interaction.rsf}.  For programming purposes, users may
     consider doing likewise.
}
\value{
   Invisibly, the interaction table.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran (2007).  Variable importance in binary regression
    trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

}
\seealso{
  \code{interaction.rsf}.
}
\examples{
#------------------------------------------------------------------------
# Explore relationship between the top two predictors from veteran data.

data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Survrsf(time,status)~., veteran, ntree = 1000, forest = TRUE)
find.interaction(v.out, npred = 2, nrep=1)

\dontrun{
#------------------------------------------------------------------------
# All pairwise interactions: PBC data.
# Use fast approximation to speed up computations.

data(pbc, package = "randomSurvivalForest") 
rsf.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000, forest = TRUE)
find.interaction(rsf.out, nrep=3, rough=T)
}

}
\keyword{file}
