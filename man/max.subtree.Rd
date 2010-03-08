\name{max.subtree}
\alias{max.subtree}
\title{Acquire Maximal Subtree Information}
\description{
    Extract maximal subtree information from a forest.  Used for
    variable selection and identifying interactions between variables.
}
\usage{
    max.subtree(object, max.order = 2, sub.order = FALSE, ...)
}
\arguments{
   \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
	forest)}. Requires \code{forest}=TRUE in the
        original \code{rsf} call.}
    \item{max.order}{Non-negative integer specifying the target number
      of order depths.  Default is to return the first and second order
      depths.  Used to identify predictive variables.  See details
      below.}
    \item{sub.order}{Set this value to TRUE to return the
      minimal depth of each variable relative to another variable.
      Used to identify interrelationship between variables.  See
      details below.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
    The maximal subtree for a variable \code{x} is the largest subtree
    whose root node splits on \code{x}.  Thus, all parent nodes of
    \code{x}'s maximal subtree have nodes that split on variables other
    than \code{x}.  The largest maximal subtree possible is the root
    node.  In general, however, there can be more than one maximal
    subtree for a variable. A maximal subtree may also not exist if
    there are no splits on the variable.  More details can be found in
    Ishwaran et al. (2009).
  
    The minimal depth of a maximal subtree measures predictiveness of a
    variable \code{x}.  It equals the shortest distance (the depth) from
    the root node to the parent node of the maximal subtree (zero is the
    smallest value possible).  The smaller the minimal depth, the more
    impact \code{x} has on prediction.  The second order depth is the
    shortest distance from the root node to the second node split using
    \code{x}.  To specify the target order depth, use the
    \code{max.order} option (e.g., setting \code{max.order}=2 returns
    the first and second order depths).

    Set sub.order=TRUE to obtain the minimal depth of a variable
    relative to another variable.  This returns a \code{p}x\code{p}
    matrix, where \code{p} is the number of variables, and entries
    [i][j] are the normalized relative minimal depth of a variable [j]
    within the maximal subtree for variable [i], where normalization
    adjusts for the size of [i]'s maximal subtree.  Entry [i][i] is the
    normalized minimal depth of i relative to the root node.  The matrix
    should be read by looking across rows (not down columns) and
    identifies interrelationship between variables.  Small [i][j]
    entries indicate interactions.  See \code{find.interaction} for
    further details.

    For competing risk data, analysis corresponds to unconditional
    values (i.e., they are non-event specific).
}
\value{
     A list with the following components:
     \item{mean}{Minimal depth averaged over a tree and forest for
       each variable.  Vector of size \code{p}.}
     \item{order}{Order depths for a given variable up to \code{max.order}
       averaged over a tree and the forest.  Matrix of dimension
       \code{p}x\code{max.order}.  If \code{max.order}=0, a matrix of
       \code{p}x\code{ntree} is returned containing the minimum maximal
       subtree distance for each variable by tree.}
     \item{count}{Averaged number of  maximal subtrees, normalized by
       the size of a tree, for each variable. Vector of size \code{p}.}
     \item{terminal}{Average terminal depth of each tree.  Vector of
       size \code{ntree}.}
     \item{nodesAtDepth}{Number of nodes per depth per tree.  
       Matrix of dimension \code{maxDepth}x\code{ntree}.}
     \item{subOrder}{Average minimal depth of a variable relative to another
       variable.  Matrix of dimension \code{p}x\code{p}.  Can be NULL.}
     \item{threshold}{Threshold used to select variables.  Variables
        whose minimal depth exceeds this value are considered to be weak
        variables.  See Ishwaran et al. (2009).}
}
\author{
    Hemant Ishwaran  \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
    Ishwaran H., Kogalur U.B>, Gorodeski E.Z, Minn A.J. and 
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc. (in press)}.
}
\seealso{
  \code{find.interaction},
  \code{varSel}.
}
\examples{
\dontrun{
#------------------------------------------------------------------------
# First and second order depths for all variables

data(veteran, package = "randomSurvivalForest")
veteran.out <- rsf(Survrsf(time, status) ~ . , data = veteran, forest = TRUE)
v <- max.subtree(veteran.out)

# first and second order depths
print(round(v$order, 3))

# weak variables have minimal depth greater than the following threshold
print(v$threshold)

}
}

\keyword{survival}
\keyword{tree}
