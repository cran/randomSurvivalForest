\name{interaction.rsf}
\alias{interaction.rsf}
\title{VIMP for Single or Grouped Variables}
\description{
  Calculate variable importance (VIMP) for a single variable or group of
  variables.
}
\usage{
    interaction.rsf(object,
                  predictorNames = NULL,
                  subset = NULL,
                  joint = TRUE,
                  rough = FALSE,
                  importance = c("randomsplit", "permute", "none")[1],
                  seed = NULL,
                  do.trace = FALSE,
                  ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
      forest)}.  Note: \code{forest}=TRUE must be used in the
      original \code{rsf} call.}
    \item{predictorNames}{Character vector of x-variable names to be
      considered.  If NULL (the default) all variables are used.  Only
      x-variables listed in the object predictor matrix will be used.}
    \item{subset}{Indices indicating which rows of the predictor
      matrix to be used (note: this applies to the \emph{object}
      predictor matrix, \code{predictors}).  Default is to use all rows.}
    \item{joint}{Should joint-VIMP or individual VIMP be calculated?
      See details below.}
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
    Using a previously grown forest, and restricting the data to that
    indicated by \code{subset}, calculate the VIMP for variables listed
    in \code{predictorNames}.  If \code{joint}=TRUE, a joint-VIMP is
    calculated.  The joint-VIMP is the importance for the group of
    variables, when the group is perturbed simultaneously.  If
    \code{joint}=FALSE, the VIMP for each variable considered separately
    is calculated.

    Depending upon the option \code{importance}, VIMP is calculated
    either by random daugther assignment, by random permutation of
    the variable(s), or none (no perturbing).
}
\value{
    A list with the following components:
    \item{err.rate}{OOB error rate for the (unperturbed) ensemble
      restricted to the subsetted data.}
    \item{err.perturb.rate}{OOB error rate for the perturbed ensemble
      restricted to the subsetted data.   Either a vector
      or a single number depending upon the option \code{joint}.}
    \item{importance}{Variable importance (VIMP).  Either a vector
      or a single number depending upon the option \code{joint}.}
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
  \code{find.interaction}.
}
\examples{
#------------------------------------------------------------------------
# Example of paired-VIMP. 
# Veteran data.

data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Survrsf(time,status)~., veteran, ntree = 1000, forest = TRUE)
interaction.rsf(v.out, c("karno","celltype"))$importance

\dontrun{
#------------------------------------------------------------------------
# Individual VIMP for data restricted to events only.
# PBC data.

data(pbc, package = "randomSurvivalForest") 
rsf.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000, forest = TRUE)
o.r <- rev(order(rsf.out$importance))
VIMP <- rsf.out$importance[o.r]
VIMP.events <- rep(0, length(VIMP))
names(VIMP.events) <- names(VIMP) 
events <- which(rsf.out$cens == 1)
VIMP.events <-
 interaction.rsf(rsf.out, names(VIMP), events, joint = FALSE)$importance
VIMP.all <- as.data.frame(cbind(VIMP.events = VIMP.events, VIMP = VIMP))
print(round(VIMP.all, 3))

#------------------------------------------------------------------------
# Estimate variability of VIMP in two ways (PBC data):
# (i)  Monte Carlo:  Estimates variability of the procedure
# (ii) Bootstrap:    Estimates statistical variability

data(pbc, package = "randomSurvivalForest") 
rsf.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000, forest = TRUE)
o.r <- rev(order(rsf.out$importance))
VIMP <- rsf.out$importance[o.r]
subset.index <- 1:nrow(rsf.out$predictors)
VIMP.mc <- VIMP.boot <- NULL
for (k in 1:100) {
  cat("iteration:", k , "\n")
  VIMP.mc <-
    cbind(VIMP.mc, interaction.rsf(rsf.out, names(VIMP), joint = FALSE)$importance)
  VIMP.boot <-
    cbind(VIMP.boot, interaction.rsf(rsf.out, names(VIMP),
      subset = sample(subset.index, replace = TRUE), joint = FALSE)$importance)
}
VIMP.mc <- as.data.frame(cbind(VIMP.mean = apply(VIMP.mc, 1, mean),
                    VIMP.sd = apply(VIMP.mc, 1, sd)))
VIMP.boot <- as.data.frame(cbind(VIMP.mean = apply(VIMP.boot, 1, mean),
                    VIMP.sd = apply(VIMP.boot, 1, sd)))
rownames(VIMP.mc) <- rownames(VIMP.boot) <- names(VIMP)
print(round(VIMP.mc, 3))
print(round(VIMP.boot, 3))
}
}
\keyword{file}

