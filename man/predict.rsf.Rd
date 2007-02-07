\name{predict.rsf}
\alias{predict.rsf}
\title{Random Survival Forest Prediction}
\description{
    Prediction on new data using Random Survival Forests.  
}
\usage{
  predict.rsf(object = NULL,
              newdata = NULL,
              proximity = FALSE, 
              do.trace = FALSE,
              ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
	forest)}. Note that \code{forest}=TRUE must be used in the
        original \code{rsf} call for prediction to work.}
    \item{newdata}{Data frame containing test data. Missing NA values
        are not encouraged. These are dealt with by removing the entire
        record if even one of its entries is NA.}
    \item{proximity}{Logical. Should proximity measure between
        test observations be calculated? Can be huge. Default is FALSE.}
    \item{do.trace}{Logical. Should trace output be enabled? Default is
        FALSE. Integer values can also be passed.  A positive value
	causes output to be printed each \code{do.trace} iteration.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  \code{predict.rsf} takes a test data set, drops it down the forest grown from
   the training data, and then computes ensemble cumulative hazard
   functions.  CHF's are predicted for all individuals in the test data
   set at the unique death time points of the original grow (training)
   data.
}    
\value{
     An object of class \code{(rsf, predict)}, which is a list with the
     following components:
     \item{call}{The original grow call to \code{rsf}.}
     \item{forest}{The grow forest.}
     \item{ntree}{Number of trees in grow forest.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
         grow forest.  Vector of length \code{ntree}.}
     \item{timeInterest}{Sorted unique event times from grow data.
         Ensemble values given for these time points only.}
     \item{n}{Sample size of test data.}
     \item{ndead}{Number of deaths in test data (can be NULL).}
     \item{Time}{Vector recording survival times from test data (can be NULL).}
     \item{Cens}{Vector recording censoring information from test data (can be NULL).}
     \item{predictorNames}{Character vector of predictor names.}
     \item{predictors}{Test data matrix of predictors used for prediction.}
     \item{ensemble}{Matrix of the ensemble cumulative hazard function
         for the test data. Each row corresponds to a test data
         individual's CHF evaluated at each of the time points in \code{timeInterest}.}
     \item{mortality}{A vector representing the estimated ensemble mortality 
         for each individual in the test data.  Ensemble mortality values should
         be interpreted in terms of total number of deaths.}
     \item{err.rate}{Vector of length \code{ntree} containing error
         rate of the test data.  Can be NULL.}
     \item{proximity}{If \code{proximity}=TRUE, a matrix recording
	 proximity of the inputs from test data is computed.  Value
	 returned is a vector of the lower diagonal of the matrix.  Use
	 \code{plot.proximity()} to extract this information.}
}
\note{
    The key deliverable is the matrix \code{ensemble} which contains the
    estimated ensemble cumulative hazard function for each individual in 
    the test data evaluated at a set of distinct time points.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran and Udaya B. Kogalur (2006).  Random Survival
    Forests.  \emph{Cleveland Clinic Technical Report}.

    L. Breiman (2001). Random forests, \emph{Machine Learning}, 45:5-32.
}	

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
