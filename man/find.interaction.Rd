\name{find.interaction}
\alias{find.interaction}
\title{Find Interactions Between Pairs of Variables}
\description{
  Find pairwise interactions between variables.
}
\usage{
   find.interaction(object, predictorNames = NULL,
          method = c("maxsubtree", "vimp")[1], sorted = TRUE, 
          npred = NULL, subset = NULL, nrep = 1, rough = FALSE,
          importance = c("randomsplit", "permute")[1],
          seed = NULL, do.trace = FALSE, ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
      forest)}.}
    \item{predictorNames}{Character vector of names of target x-variables.
      Default is to use all variables.}
    \item{method}{Method of analysis: maximal subtree or VIMP.  See details
      below.}
    \item{sorted}{Should variables be sorted?}
    \item{npred}{Use the first npred ordered variables. Default is to use all variables.}
    \item{subset}{Indices indicating which rows of the predictor matrix
      to be used (note: this applies to the \emph{object} predictor
      matrix, \code{predictors}).  Default is to use all rows.}
    \item{nrep}{Number of Monte Carlo replicates.  Applies only when
      \code{method}="\code{vimp}".}
    \item{rough}{Should fast approximation be used? Applies only when
      \code{method}="\code{vimp}".}
    \item{importance}{Type of variable importance (VIMP).  Applies only
      when \code{method}="\code{vimp}".}
    \item{seed}{Seed (negative integer) for random number generator.}
    \item{do.trace}{Logical.  Should trace output be enabled? Integer
      values can also be passed.  A positive value causes output to be
      printed each \code{do.trace} iteration.  Applies only when
      \code{method}="\code{vimp}".}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
     Using a previously grown forest, identify pairwise interactions for
     all pairs of variables from a specified list.  There are two
     distinct approaches specified by the \code{method} option.

     If \code{method}="\code{maxsubtree}", then a maximal subtree
     analysis is used.  In this case, a matrix is returned where entries
     [i][i] are the normalized minimal depth of variable [i] relative to
     the root node (normalized w.r.t. the size of the tree) and entries
     [i][j] indicate the normalized minimal depth of a variable [j]
     w.r.t. the maximal subtree for variable [i] (normalized w.r.t. the
     size of [i]'s maximal subtree).  Smaller [i][i] entries indicate
     predictive variables.  Small [i][j] entries having small [i][i]
     entries are a sign of an interaction between variable i and j
     (note: the user should scan rows, not columns, for small entries).
     See Ishwaran et al. (2010) for more details.
     
     If \code{method}="\code{vimp}", then a joint-VIMP approach is used.
     Two variables are paired and their paired VIMP calculated (refered
     to as Paired importance).  The VIMP for each separate variable is
     also calculated.  The sum of these two values is refered to as
     Additive importance.  A large positive or negative difference
     between Paired and Additive indicates an association worth pursuing
     if the VIMP's for each variable are reasonably large.  See Ishwaran
     (2007) for more details.

     Computations might be slow depending upon the size of the data and
     the forest.  In such cases, consider setting \code{npred} to a
     smaller number, or using \code{rough}=TRUE if
     \code{method}="\code{vimp}".  If \code{method}="\code{maxsubtree}",
     consider using a smaller number of trees in the original grow call.

     If \code{nrep} is greater than 1, the analysis is repeated
     \code{nrep} times and results averaged over the replications
     (applies only when \code{method}="\code{vimp}").

    For competing risk data, maximal subtree analyses correspond to
    unconditional values (i.e., they are non-event specific).  Setting
    \code{method}="\code{vimp}", however, yields pairwise interactions
    for both event and non-event specific settings.  } \value{
    Invisibly, the interaction table (a list for competing risk data) or
    the maximal subtree matrix.
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    Ishwaran H. (2007).  Variable importance in binary regression
    trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

    Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and 
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.
}
\seealso{
  \code{max.subtree},
  \code{vimp}.
}
\examples{
\dontrun{
#------------------------------------------------------------------------
# Maximal subtree approach, top 8 predictors (PBC data).

data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Surv(days,status) ~ ., pbc, nsplit = 10)
find.interaction(pbc.out, npred = 8)

#------------------------------------------------------------------------
# VIMP approach (PBC data). 
# Use fast approximation to speed up computations.

data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Surv(days,status) ~ ., pbc, nsplit = 10)
find.interaction(pbc.out, method = "vimp", nrep=3, rough=T)

#------------------------------------------------------------------------
# Competing risks (WIHS data).

data(wihs, package = "randomSurvivalForest")
wihs.out <- rsf(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 200)
find.interaction(wihs.out, method = "vimp")
}
}
\keyword{survival}
\keyword{tree}
