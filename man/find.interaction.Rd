\name{find.interaction}
\alias{find.interaction}
\title{Find Interactions Between Pairs of Variables}
\description{
  Test for pairwise interactions between variables by comparing
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
    \item{sort}{Should variables be sorted by importance values (only
        applies if importance values are available and
        \code{predictorNames}=NULL)?  Default is TRUE.}
    \item{predictorNames}{Character vector of variable names.  Only these
        variables will be considered.  Default is all.}
    \item{n.pred}{Number of variables to be plotted (only applies when
	\code{predictorNames}=NULL).  Default is all.}
    \item{n.rep}{Number of Monte Carlo replicates.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
     Pairwise interactions between the first \code{n.pred} variables
     are investigated (if \code{sort}=TRUE analysis is based on the
     first \code{n.pred} as ordered by importance values).  Variables
     are paired and then removed from the model.  The increase in error
     rate is computed and is defined to be the 'Paired' importance value
     of the two variables.  The increase in error rate if each
     variable is removed individually is also computed and the sum of
     these two values yields an 'Additive' importance value.  Large
     positive or negative differences between 'Paired' and 'Additive'
     values indicate potentially interesting associations.  

     Computations can be fairly heavy if the data is large, thus users
     should consider setting \code{n.pred} to a fairly small value to
     gauge computational times.

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
    H. Ishwaran (2007).  Variable importance in binary
    trees.  \emph{Cleveland Clinic Technical Report}.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time,status)~., veteran, ntree = 1000)
  find.interaction(v.out, n.pred = 2, n.rep=1)
}
\keyword{file}
