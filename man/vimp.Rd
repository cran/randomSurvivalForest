\name{vimp}
\alias{vimp}
\title{VIMP for Single or Grouped Variables}
\description{
  Calculate variable importance (VIMP) for a single variable or group of
  variables.
}
\usage{
    vimp(object,
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
      forest)}. Requires \code{forest}=TRUE in the original \code{rsf}
      call.}
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
    \item{importance}{Type of VIMP.}
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
    calculated.  The joint-VIMP is the importance for a group of
    variables, when the group is perturbed simultaneously.  If
    \code{joint}=FALSE, the VIMP for each variable considered separately
    is calculated.

    Depending upon the option \code{importance}, VIMP is calculated
    either by random daugther assignment, by random permutation of
    the variable(s), or without perturbation (none).

    For competing risk data, VIMP and error rates are given for the
    ensemble CHF and the conditional CHF (CCHF) for each event type.
}
\value{
    A list with the following components:
    \item{err.rate}{OOB error rate for the (unperturbed) ensemble
      restricted to the subsetted data.}
    \item{err.perturb.rate}{OOB error rate for the perturbed ensemble
      restricted to the subsetted data.   Dimension depends upon the
      option \code{joint}.}
    \item{importance}{Variable importance (VIMP). Dimension depends
      upon the option \code{joint}.} 
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
  \code{vimp}.
}
\examples{
#------------------------------------------------------------------------
# Example of paired-VIMP. 
# Veteran data.

data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Survrsf(time,status)~., veteran, ntree = 1000, forest = TRUE)
vimp(v.out, c("karno","celltype"))$importance

\dontrun{
#------------------------------------------------------------------------
# Individual VIMP for data restricted to events only.
# PBC data.

data(pbc, package = "randomSurvivalForest") 
rsf.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000,
               forest = TRUE, nsplit = 3)
o.r <- rev(order(rsf.out$importance))
imp <- rsf.out$importance[o.r]
imp.events <- rep(0, length(imp))
events <- which(rsf.out$cens == 1)
imp.events <-
 vimp(rsf.out, names(imp), events, joint = FALSE)$importance
imp.all <- as.data.frame(cbind(imp.events = imp.events, imp = imp))
print(round(imp.all, 3))

#------------------------------------------------------------------------
# Estimate variability of VIMP in two ways (PBC data):
# (i)  Monte Carlo:  Estimates variability of the procedure
# (ii) Bootstrap:    Estimates statistical variability

data(pbc, package = "randomSurvivalForest") 
rsf.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000,
               nsplit = 3, forest = TRUE)
o.r <- rev(order(rsf.out$importance))
imp.names <- names(rsf.out$importance[o.r])
subset.index <- 1:nrow(rsf.out$predictors)
imp.mc <- imp.boot <- NULL
for (k in 1:100) {
  cat("iteration:", k , "\n")
  imp.mc <-
    cbind(imp.mc, vimp(rsf.out, imp.names, joint = FALSE)$importance)
  imp.boot <-
    cbind(imp.boot, vimp(rsf.out, imp.names,
    subset = sample(subset.index, replace = TRUE), joint = FALSE)$importance)
}
imp.mc <- as.data.frame(cbind(imp.mean = apply(imp.mc, 1, mean),
                    imp.sd = apply(imp.mc, 1, sd)))
imp.boot <- as.data.frame(cbind(imp.mean = apply(imp.boot, 1, mean),
                    imp.sd = apply(imp.boot, 1, sd)))
print(round(imp.mc, 3))
print(round(imp.boot, 3))
}
}
\keyword{file}

