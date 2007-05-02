\name{find.interaction}
\alias{find.interaction}
\title{Find Interactions Between Pairs of Predictors}
\description{
  Test for pairwise interactions between predictors by comparing
  pairwise importance values to additive individual importance
  values.
}
\usage{
    find.interaction(x,
                  sort = TRUE,
                  predictorNames = NULL,
                  n.pred = NULL,
                  n.rep = 1,
                  ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)}.}
    \item{sort}{Should predictors be sorted by importance values (only
        applies if importance values are available and
        \code{predictorNames}=NULL)?  Default is TRUE.}
    \item{predictorNames}{Character vector of predictor names.  Only these
        predictors will be considered.  Default is all.}
    \item{n.pred}{Number of predictors to be plotted (only applies when
	\code{predictorNames}=NULL).  Default is all.}
    \item{n.rep}{Number of Monte Carlo replicates.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
     Pairwise interactions between the first \code{n.pred} predictors
     are investigated (if \code{sort}=TRUE analysis is based on the
     first \code{n.pred} as ordered by importance values).  Predictors
     are paired and then removed from the model.  The increase in error
     rate is computed and is defined to be the 'Paired' importance value
     of the two predictors.  The increase in error rate if each
     predictor is removed individually is also computed and the sum of
     these two values yields an 'Additive' importance value.  The
     difference between the 'Paired' and 'Additive' values can be used
     to identify potential interactions.  If the difference is positive,
     and relatively large, the two predictors may have a tree-type
     interaction.

     Computations can be fairly heavy if the data is large, thus users
     should consider setting \code{n.pred} to a fairly small value to
     guage computational times.

     If \code{n.rep} is greater than 1, the analysis is replicated
     \code{n.rep} times.  Results reported are averaged values in this
     case.
}
\value{
   Invisibly, the interaction table.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran and Udaya B. Kogalur (2006).  Random Survival
    Forests.  \emph{Cleveland Clinic Technical Report}.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time,status)~., veteran, ntree = 1000)
  find.interaction(v.out, n.pred = 3, n.rep=5)
}
\keyword{file}
