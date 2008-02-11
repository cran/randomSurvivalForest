\name{predict.rsf}
\alias{predict.rsf}
\title{Random Survival Forest Prediction}
\description{
    Prediction on new data using Random Survival Forests.  
}
\usage{
  predict.rsf(object = NULL,
              newdata = NULL,
              importance = c("randomsplit", "permute", "none")[1],
              na.action = c("na.omit", "na.impute")[1],
              proximity = FALSE,
              seed = NULL,
              do.trace = FALSE,
              ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
	forest)}. Note that \code{forest}=TRUE must be used in the
        original \code{rsf} call for prediction to work.}
    \item{newdata}{Data frame containing test data. Missing values
        allowed.}
    \item{importance}{Method used to compute variable importance (VIMP).  Only
        applies when test data contains outcomes.}
    \item{na.action}{Action to be taken if the data contain NA's.  Possible
        values are \code{na.omit}, which removes the entire record if
        even one of its entries is NA, and \code{na.impute}, which
        imputes the test data.  See details below.}
    \item{proximity}{Logical. Should proximity measure between
        test observations be calculated? Can be large. Default is FALSE.}
    \item{seed}{Seed for random number generator.  Must be a negative
        integer (the R wrapper handles incorrectly set seed values).}
    \item{do.trace}{Logical. Should trace output be enabled? Default is
        FALSE. Integer values can also be passed.  A positive value
	causes output to be printed each \code{do.trace} iteration.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  \code{predict.rsf} takes a test data set, drops it down the forest
   grown from the training data, and computes an ensemble cumulative
   hazard function (CHF).  CHF's are calculated for all individuals in
   the test data set at the unique death time points of the original
   grow (training) data.  The error rate and VIMP for variables are
   computed on the test data if outcome information is available.
   Setting the option \code{na.action} to \code{na.impute} imputes
   missing test data (x-variables or outcomes) using the forest grown on
   the training data (Ishwaran et al. 2007).  Only training data is used
   in imputing test data to avoid biasing error rates.
}    
\value{
     An object of class \code{(rsf, predict)}, which is a list with the
     following components:
     \item{call}{The original grow call to \code{rsf}.}
     \item{forest}{The grow forest.}
     \item{ntree}{Number of trees in the grow forest.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
         grow forest.  Vector of length \code{ntree}.}
     \item{timeInterest}{Sorted unique event times from grow (training) data.
         Ensemble values given for these time points only.}
     \item{n}{Sample size of test data (depends upon NA's, see \code{na.action}).}
     \item{ndead}{Number of deaths in test data (can be NULL).}
     \item{time}{Vector of survival times from test data (can be NULL).}
     \item{cens}{Vector of censoring indicators from test data (can be NULL).}
     \item{predictorNames}{Character vector of variable names.}
     \item{predictors}{Test data matrix of x-variables used for prediction.}
     \item{ensemble}{Matrix containing the ensemble CHF for the test data. Each
         row corresponds to a test data individual's CHF evaluated at
	 each of the time points in \code{timeInterest}.}
     \item{mortality}{Vector containing ensemble mortality for each
         individual in the test data.  Ensemble mortality values should
         be interpreted in terms of total number of training deaths.}
     \item{err.rate}{Vector of length \code{ntree} containing error
         rate of the test data.  Can be NULL.}
     \item{importance}{VIMP of each variable in the test data. Can be NULL.}
     \item{proximity}{If \code{proximity}=TRUE, a matrix recording
	 proximity of the inputs from test data is computed.  Value
	 returned is a vector of the lower diagonal of the matrix.  Use
	 \code{plot.proximity()} to extract this information.}
     \item{imputedIndv}{Vector of indices of records in test data with
       missing values.  Can be NULL.}
     \item{imputedData}{Matrix of imputed test data.  First two columns
        are censoring and survival time, respectively. The remaining
        columns are the x-variables.  Row i contains imputed outcomes
        and x-variables for row \code{imputedIndv}[i] of
        \code{predictors}.  Can be NULL.}
}
\note{
    The key deliverable is the matrix \code{ensemble} which contains the
    ensemble CHF for each individual in the test data evaluated at a set
    of distinct time points.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    L. Breiman (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    H. Ishwaran, U.B. Kogalur, E.H. Blackstone and M.S. Lauer
    (2007).  Random survival forests, \emph{Cleveland Clinic Technical
    Report}.

    H. Ishwaran, U.B. Kogalur (2007).  Random survival forests for R,
    \emph{Rnews}, 7/2:25-31.
}	
\seealso{
  \code{rsf},
  \code{print.rsf},
  \code{plot.ensemble},
  \code{plot.variable},
  \code{plot.error},
  \code{plot.proximity},
  \code{pmml_to_rsf},
  \code{rsf_to_pmml}.
}
\examples{
data(veteran, package = "randomSurvivalForest")
veteran.out <- rsf(Survrsf(time, status)~., forest = TRUE, data = veteran)
baseForest <- veteran.out$forest
veteran.pred1 <- predict.rsf(baseForest, veteran, proximity = FALSE)
veteran.pred2 <- predict.rsf(veteran.out, veteran, proximity = TRUE)
}

\keyword{survival}
\keyword{tree}
