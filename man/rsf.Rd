\name{rsf}
\alias{rsf}
\title{Random Survival Forests Primary R Function}
\description{
    \pkg{rsf} implements Ishwaran and Kogalur's Random Survival
    Forests algorithm for right censored survival data.  The algorithm
    uses a binary recursive tree growing procedure with different
    splitting rules for growing an ensemble cumulative hazard
    function.  An out-of-bag (OOB) estimate of Harrell's concordance
    index is provided for assessing prediction.
}
\usage{
rsf(formula = formula(data),
    data = sys.parent(),
    ntree = 100,
    mtry = NULL,
    nodesize = NULL,
    ntime = NULL,
    splitrule = c("conserve", "logrank"),
    proximity = FALSE, 
    seed = NULL,
    ntime.upper = 1000,
    trace = FALSE,
    ...)
}
\arguments{
    \item{formula}{A symbolic description of the model that is
        to be fit.  The details of model specification are given below.}
    \item{data}{Data frame containing the predictors (variables) in the
        model. Missing NA values are not encouraged.  These are dealt with
        by removing the entire record if even one of its entries is NA.}
    \item{ntree}{Number of trees to grow.  This should not be set
        to a number too small, in order to ensure that every input
        row gets predicted at least a few times.}
    \item{mtry}{Number of predictors randomly sampled as candidates
        at each split.  The default is sqrt(p) where p
        is the number of predictors.}
    \item{nodesize}{Minimum number of deaths (must have unique event
        times) in a terminal node.  Default equals
        min(3,\code{ndead}).  Larger values cause
        smaller trees to be grown.  The algorithm thus takes less time.}
    \item{ntime}{Maximum number of distinct time points
        considered for evaluating ensemble.  Default equals number of distinct events.}
    \item{splitrule}{Splitting rule used for splitting nodes in growing
        the survival tree.  Possible values are \dQuote{logrank} and \dQuote{conserve}.
        Default value is \dQuote{conserve}.  See details below.}
    \item{proximity}{Logical.  Should proximity measure between
        observations be calculated?  Creates an \code{n}x\code{n}
        matrix (which can be huge).  Default is FALSE.}
    \item{seed}{Seed for random number generator.  Must be a negative
        integer (the R wrapper handles incorrectly set seed values).}
    \item{ntime.upper}{Upper limit allowable for \code{ntime}.}
    \item{trace}{Logical. Should trace output be enabled? Default is FALSE.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  A typical formula has the form \code{Survrsf(time, censoring) ~ terms}
  where \code{time} is survival time and \code{censoring} is a binary
  censoring indicator.  Note that censoring must be coded as 0=censored
  and 1=death (event).  Unlike standard formula calls in R, a
  specification of the form \code{first:second}, used to indicate terms
  obtained by taking the interactions of all terms in \code{first} with
  all terms in \code{second}, will be interpreted by \pkg{rsf} to simply
  mean \code{first+second}.  This is because interactions are implicit
  in growing random forests.  Attempts to coerce terms to include higher
  order terms using functions from other packages like \code{bs()} or
  \code{poly()} are also fruitless.  A formula such as
  \code{Survrsf(time, censoring) ~ bs(term1, df=5) + term2} is interpreted
  exactly the same as \code{Survrsf(time, censoring) ~ term1 + term2}.

  Predictors which are encoded as factors will be coerced into dummy
  variables.  These dummy variables will be automatically labelled using
  the original predictor name.  For example, if marital status is a
  predictor named \dQuote{marital} encoded as a factor with levels
  \dQuote{S}, \dQuote{M} and \dQuote{D}, two new dummy variables will be
  created labeled \dQuote{maritalM} and \dQuote{maritalS}.
  
  The \dQuote{logrank} splitting rule grows trees by splitting nodes by
  maximization of the log-rank test statistic.  The default splitting
  rule is \dQuote{conserve} which splits nodes by finding daughters
  closest to the conservation of events principle (the cumulative hazard
  summing to total number of deaths; see references below).  The
  conservation of events splitting rule is less susceptible to an
  edge-effect bias as sometimes seen with the log-rank splitting rule.
  Both rules often yield roughly the same prediction error performance,
  although \dQuote{logrank} is sometimes slightly better if proportional
  hazards is not too severely violated.  The \dQuote{conserve} splitting
  rule is faster, especially with larger data sets.

  For \emph{very} large data sets, discretizing continuous predictors
  can greatly speed up computational times.  Discretization does not
  have to be overly granular for gains to be seen.
}    
\value{
     An object of class \code{randomSurvivalForest}, which is a list with the
     following components:
     \item{call}{The original call to \pkg{rsf}.}
     \item{formula}{The original formula used in calling \pkg{rsf}.}
     \item{terms}{The terms used in formula.}
     \item{n}{Sample size of the data.}
     \item{ndead}{Number of deaths.}
     \item{ntree}{Number of trees grown.}
     \item{mtry}{Number of predictors sampled for splitting at each node.}
     \item{nodesize}{Minimum size of terminal nodes.}
     \item{splitrule}{Splitting rule used.}
     \item{Time}{Vector of length \code{n} recording survival times.}
     \item{Died}{Vector of length \code{n} recording censoring
       information (0=censored, 1=death).}
     \item{Time.unique}{Sorted unique event times.  Ensemble values are given
       for these time points only.}
     \item{prednames}{A character vector of predictor names used to
           grow the forest.}
     \item{predictors}{Matrix of predictors used to grow the forest.}
     \item{ensemble}{A matrix of the ensemble cumulative hazard
       function with each row corresponding to an individual's CHF
       evaluated at each of the time points in \code{Time.unique}.}
     \item{mortality}{A vector of length \code{n}, with each value
       representing the estimated ensemble mortality for an
       individual in the data.  Ensemble mortality values should
       be interpreted in terms of total number of deaths.}
     \item{err.rate}{Vector of length \code{ntree} containing OOB error
       rates for the ensemble, with the b-th element being the error
       rate for the ensemble formed using the first b trees.  Error
       rate is 1-C, where C is Harrell's concordance index.  Error
       rates are between 0 and 1, with 0.5 representing the benchmark
       value of a procedure based on random guessing.  A value of 0
       is perfect.}
     \item{leaf.count}{Number of terminal nodes for each tree in the forest.  Vector
       of length \code{ntree}.}
     \item{proximity}{If \code{proximity}=TRUE, a matrix
       of dimension \code{n}x\code{n} of proximity measures among the input is
	  calculated (based on the frequency that pairs of data points
	  are in the same terminal nodes). Value returned is a vector of
	  the lower diagonal of the matrix.  See the call
	  \code{plot.proximity()} for how to extract information.}
}
\note{
    The key deliverable is the matrix \code{ensemble} which contains the
    estimated ensemble cumulative hazard function for each individual
    evaluated at a set of distinct time points.  The vector
    \code{mortality} is a weighted sum over the columns of
    \code{ensemble}, weighted by the number of individuals at risk at
    the different time points.  Entry \code{i} of the vector represents
    the estimated total mortality of individual \code{i} in terms of
    total number of deaths.  In other words, if \code{i} has a mortality
    value of 100, then if all individuals had the same predictor as
    \code{i}, there would be on average 100 deaths in the dataset.
    
    Different R wrappers are provided with the package to aid in
    interpreting the ensemble.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran and Udaya B. Kogalur (2006).  Random Survival
    Forests.  \emph{Cleveland Clinic Technical Report}.

    L. Breiman (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    F.E. Harrell et al. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247:2543-2546.

    D.C. Naftel, E.H. Blackstone and M.E. Turner (1985).
    Conservation of events, unpublished notes.
}	

}
\seealso{
  \code{print.rsf},
  \code{plot.ensemble},
  \code{plot.variable},
  \code{plot.error},
  \code{plot.proximity}
  \code{Survrsf}.
}
\examples{
# Example 1:  Minimal argument call using included sample data set.
data(veteran, package = "randomSurvivalForest")
veteran.out <- rsf(Survrsf(time, status)~., data = veteran)
print.rsf(veteran.out)
plot.error(veteran.out)

# Example 2:  Full argument call using included sample data set.
veteran.f <- as.formula(Survrsf(time, status)~.)
data(veteran)
ntree <- 100
mtry <- 2
nodesize <- 3
ntime <- NULL
splitrule <- "logrank"
proximity <- FALSE
seed <- -1
ntime.upper <- 1000
trace <- TRUE
veteran2.out <- rsf(veteran.f, veteran, ntree, mtry, nodesize,
                    ntime, splitrule, proximity, seed, ntime.upper, trace)
print.rsf(veteran2.out)

# Example 3:  Veteran's Administration lung cancer trial from
#  Kalbfleisch & Prentice.  Randomized trial of two treatment
#  regimens for lung cancer.  Assumes "survival"
#  library is loaded.  This is a standard survival analysis
#  data set, included in the package.

if (library("survival", logical.return = TRUE))
{
        data(veteran, package = "randomSurvivalForest")
        veteran3.out <- rsf(Survrsf(time, status)~ celltype,
                       veteran,
                       ntree = 1000)
        plot.ensemble(veteran3.out)
        par(mfrow = c(1,1))
        plot(survfit(Surv(time, status) ~ celltype, data = veteran))
}


# Example 4:  Primary biliary cirrhosis (PBC) of the liver.
#  Data found in Appendix D.1 of Fleming and Harrington, Counting
#  Processes and Survival Analysis, Wiley, 1991 (only differences
#  are that age is in days and sex and stage variables are not
#  missing for observations 313-418).  This is a sample data
#  set included in the package.

data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000)
print.rsf(pbc.out)
plot.variable(pbc.out)  ### par(ask = TRUE) might be useful here


# Example 5:  Compare Cox regression to Random Survival Forests
#  using the PBC data.  Compute OOB estimate of Harrell's
#  concordance index (B = 100).  Assumes "Hmisc" and "survival"
#  libraries are loaded. 

if (library("survival",logical.return = TRUE) 
    & library("Hmisc",logical.return = TRUE))
{
    data(pbc, package = "randomSurvivalForest")
    pbc2.out <- rsf(Survrsf(days,status)~., pbc, mtry = 2, ntree = 1000)
    B <- 100 
    cox.err <- 0 
    cox.f <- as.formula(Surv(days,status)~.)  
    pbc.data <- pbc[apply(is.na(pbc), 1, sum) == 0,] ##remove NA's 
    cat("Out-of-bag Cox Analysis ...", "\n")
    for (b in 1:B){
        cat("Cox bootstrap", b, "\n") 
        bag.sample <- sample(1:dim(pbc.data)[1],
                             dim(pbc.data)[1],
                             replace = TRUE) 
        oob.sample <- setdiff(1:dim(pbc.data)[1],bag.sample)
        train <- pbc.data[bag.sample,]
        test <- pbc.data[oob.sample,]
        cox.out <- coxph(cox.f, train)
        cox.predict <- predict(cox.out, test)
        cox.err <- cox.err +
	     rcorr.cens(cox.predict, 
                  Surv(pbc.data$days[oob.sample],
                  pbc.data$status[oob.sample])
                 )[1]
     }
     cat("Error rates:", "\n")
     cat("Random Survival Forests:",
         pbc2.out$err.rate[pbc2.out$ntree],
         "\n")
     cat("         Cox Regression:", cox.err[1]/B, "\n")
}


# Example 6:  Using an external data set.
    \dontrun{
        file.in <- "other.data"
        other.data <- read.table(file.in, header = TRUE)

        rsf.f <- as.formula(Survrsf(time, status)~.)
        rsf.out <- rsf(formula = rsf.f, data = other.data) 
    }
}

\keyword{survival}
\keyword{tree}
