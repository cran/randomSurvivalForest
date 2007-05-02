\name{rsf.default}
\alias{rsf.default}
\alias{rsf}
\alias{randomSurvivalForest}
\title{Random Survival Forest Entry Point}
\description{
  Ishwaran and Kogalur's Random Survival Forests algorithm for right
  censored survival data (Ishwaran and Kogalur, 2006).  This is a direct
  extension of Breiman's Random Forests method (Breiman, 2001) to
  survival analysis settings.  Algorithm uses a binary recursive tree
  growing procedure with different splitting rules for growing an
  ensemble cumulative hazard function.  An out-of-bag (OOB) estimate of
  Harrell's concordance index (Harrell, 1982) is provided for assessing
  prediction.  Importance values for predictors can also be computed.
  Prediction on test data is also available.  Note that this is the
  default generic method for the package.
}
\usage{
\method{rsf}{default}(formula,
    data = NULL,
    ntree = 1000,
    mtry = NULL,
    nodesize = NULL,
    splitrule = c("logrank", "conserve", "logrankscore", "logrankapprox")[1],
    importance = TRUE,
    big.data = FALSE,
    predictorWt = NULL,
    forest = FALSE,
    do.trace = FALSE,
    proximity = FALSE, 
    seed = NULL,
    ntime = NULL,
    add.noise = FALSE,
    ...)
}
\arguments{
    \item{formula}{A symbolic description of the model that is
        to be fit.  The details of model specification are given below.}
    \item{data}{Data frame containing the predictors (variables) in the
        model. Missing NA values are not encouraged.  These are dealt with
        by removing the entire record if even one of its entries is NA
	(applies to those entries specifically called in \code{formula}).}
    \item{ntree}{Number of trees to grow.  This should not be set
        to a number too small, in order to ensure that every input
        row gets predicted at least a few times.}
      \item{mtry}{Number of predictors randomly sampled at each split.
	  The default is sqrt(\code{p}) where \code{p} equals the number of
	  predictors.}
    \item{nodesize}{Minimum number of deaths (must have unique event
        times) in a terminal node.  Default equals
        min(3,\code{ndead}).  Larger values cause smaller trees to be
	grown. The algorithm thus takes less time.}
    \item{splitrule}{Splitting rule used for splitting nodes in growing
        the survival tree.  Possible values are \dQuote{logrank},
        \dQuote{conserve}, \dQuote{logrankscore} and \dQuote{logrankapprox}.
        Default value is \dQuote{logrank}.  See details below.}
    \item{importance}{Logical.  Should importance of predictors be
      estimated?  Default is TRUE.}
    \item{big.data}{Logical.  Set this value to true when the number of predictors
      \code{p} is \emph{very} large, or the data is very large.  See details below.}
    \item{predictorWt}{Vector of positive weights where entry \code{k}
        represents the likehood of selecting predictor 
        \code{k} as a candidate for splitting.  Default is to use uniform
        weights.  Vector must be of dimension \code{p}, where \code{p}
        equals the number of predictors.}
    \item{forest}{Logical. Should the forest object be returned?  Used
        for prediction on new data.  Default is FALSE.}
    \item{do.trace}{Logical. Should trace output be enabled? Default is
        FALSE. Integer values can also be passed.  A positive value
	causes output to be printed each \code{do.trace} iteration.}
    \item{proximity}{Logical. Should proximity measure between
        observations be calculated?  Creates an \code{n}x\code{n}
        matrix (which can be huge).  Default is FALSE.}
    \item{seed}{Seed for random number generator.  Must be a negative
        integer (the R wrapper handles incorrectly set seed values).}
    \item{ntime}{Maximum number of desired distinct time points
        considered for evaluating ensemble.  Default equals number of
	distinct events.}
    \item{add.noise}{Logical. Should noise variable be added?  Default
        is FALSE.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{  
  The default rule, the \dQuote{logrank} splitting rule, grows trees by
  splitting nodes by maximization of the log-rank test statistic
  (Leblanc and Crowley, 1993).  The \dQuote{conserve} splitting rule
  splits nodes by finding daughters closest to the conservation of
  events principle (the cumulative hazard summing to total number of
  deaths; see Naftel, Blackstone and Turner, 1985).  The conservation of
  events splitting rule is less susceptible to an edge-effect bias as
  sometimes seen with the log-rank splitting rule.  The
  \dQuote{logrankscore} splitting rule uses a standardized log-rank
  statistic (Hothorn and Lausen, 2003).  The \dQuote{logrankapprox}
  splitting rule splits nodes using an approximation to the log-rank
  test (suggested by Michael Leblanc; also see Cox and Oakes, page 105).

  All four rules often yield roughly the same prediction error
  performance, but users are encouraged to try all methods in any given
  example.  The \dQuote{logrankapprox} splitting rule is almost always
  fastest, especially with large data sets.  After that,
  \dQuote{conserve} is often second fastest.  For \emph{very} large data
  sets, discretizing continuous predictors and/or the observed survival
  times can greatly speed up computational times.  Discretization does
  not have to be overly granular for substantial gains to be seen.

  A typical formula has the form \code{Survrsf(time, censoring) ~
  terms}, where \code{time} is survival time and \code{censoring} is a
  binary censoring indicator.  Note that censoring must be coded as
  0=censored and 1=death (event) and \code{time} must be strictly
  positive.

  Predictors which are encoded as factors will be coerced into dummy
  variables.  These dummy variables will be automatically labelled using
  the original predictor name.  For example, if marital status is a
  predictor named \dQuote{marital} encoded as a factor with levels
  \dQuote{S}, \dQuote{M} and \dQuote{D}, two new dummy variables will be
  created labeled \dQuote{maritalM} and \dQuote{maritalS}.

  For \emph{very} large data sets, or data with a \emph{large} number of
  predictors, users should consider setting the logical flag
  \code{big.data} to TRUE.  This bypasses the large overhead needed by R
  in creating design matrices and parsing formula.  Users should be
  aware, however, that predictors are not processed and are interpreted
  \emph{as is} when this option is turned on.  Think of the data frame
  as containing time and censoring information and the rest of the data
  as the pre-processed design matrix when this option is on.  Side
  effects are that factors are not transformed to dummy values (in fact
  they are coerced to NAs), and transformations used in the formula
  (such as logs etc.) are ignored.
  
}    
\value{
     An object of class \code{(rsf, grow)}, which is a list with the
     following components:
     \item{call}{The original call to \code{rsf}.}
     \item{formula}{The formula used in the call.}
     \item{n}{Sample size of the data.}
     \item{ndead}{Number of deaths.}
     \item{ntree}{Number of trees grown.}
     \item{mtry}{Number of predictors randomly selected for splitting at each node.}
     \item{nodesize}{Minimum size of terminal nodes.}
     \item{splitrule}{Splitting rule used.}
     \item{Time}{Vector of length \code{n} recording survival times.}
     \item{Cens}{Vector of length \code{n} recording censoring
         information (0=censored, 1=death).}
     \item{timeInterest}{Sorted unique event times.  Ensemble values are given
         for these time points only.}
     \item{predictorNames}{A character vector of the predictor names used in
         growing the forest.}
     \item{predictorWt}{Vector of positive weights used for randomly
         sampling predictors for splitting.}
     \item{predictors}{Matrix of predictors used to grow the forest.}
     \item{ensemble}{A matrix of the ensemble cumulative hazard
         function with each row corresponding to an individual's CHF
         evaluated at each of the time points in \code{timeInterest}.}
     \item{oob.ensemble}{Same as \code{ensemble}, but based on the OOB
         cumulative hazard function.}
     \item{mortality}{A vector of length \code{n}, with each value
         representing the estimated ensemble mortality for an
         individual in the data.  Ensemble mortality values should
         be interpreted in terms of total number of deaths.}
     \item{oob.mortality}{Same as \code{mortality}, but based on \code{oob.ensemble}.}
     \item{err.rate}{Vector of length \code{ntree} containing OOB error
         rates for the ensemble, with the b-th element being the error
         rate for the ensemble formed using the first b trees. Error
         rate is 1-C, where C is Harrell's concordance index. Error
         rates are between 0 and 1, with 0.5 representing the benchmark
         value of a procedure based on random guessing. A value of 0
         is perfect.}
     \item{leaf.count}{Number of terminal nodes for each tree in the forest.  Vector
         of length \code{ntree}.}
     \item{importance}{Importance measure of each predictor.  For each
         predictor this is the difference
         in the OOB error rate when the predictor is randomly permuted
         compared to the OOB error rate without any permutation (i.e. the
         B-th component of \code{err.rate}).  Large positive values
         indicate informative variables, whereas small values, or negative
         values, indicate predictors unlikely to be informative.} 
     \item{forest}{If \code{forest}=TRUE, the forest object is returned.
         This object can then be used for prediction with new data sets.}
     \item{proximity}{If \code{proximity}=TRUE, a matrix
         of dimension \code{n}x\code{n} of proximity measures among the input is
         calculated (based on the frequency that pairs of data points
         are in the same terminal nodes). Value returned is a vector of
         the lower diagonal of the matrix.  Use 
         \code{plot.proximity()} to extract this information.}
}
\note{
    The key deliverable is the matrix \code{ensemble} which contains the
    estimated ensemble cumulative hazard function for each individual
    evaluated at a set of distinct time points (an OOB ensemble,
    \code{oob.ensemble}, is also returned).  The vector \code{mortality}
    (likewise \code{oob.mortality}) is a weighted sum over the columns of
    \code{ensemble}, weighted by the number of individuals at risk at
    the different time points.  Entry \code{i} of the vector represents
    the estimated total mortality of individual \code{i} in terms of
    total number of deaths.  In other words, if \code{i} has a mortality
    value of 100, then if all individuals had the same predictors as
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

    M. LeBlanc and J. Crowley (1993).  Survival trees by goodness of split,
    \emph{J. Amer. Stat. Assoc.}, 88:457-467.
    
    D.C. Naftel, E.H. Blackstone and M.E. Turner (1985). Conservation of
    events, unpublished notes.

    T. Hothorn and B. Lausen (2003). On the exact distribution of maximally selected
    rank statistics, \emph{Computational Statistics & Data Analysis}, 43:121-137.    

    D.R. Cox and D. Oakes (1988).  \emph{Analysis of Survival Data},
    Chapman and Hall.
    
    A. Liaw and M. Wiener (2002). Classification and regression by randomForest,
    \emph{R News}, 2:18-22.    
}
\seealso{
  \code{plot.ensemble},
  \code{plot.variable},
  \code{plot.error},
  \code{plot.proximity},
  \code{predict.rsf},
  \code{print.rsf},
  \code{find.interaction},
  \code{pmml_to_rsf},
  \code{rsf_to_pmml},
  \code{Survrsf}.
}
\examples{
# Example 1:  Veteran's Administration lung cancer trial from
# Kalbfleisch & Prentice.  Randomized trial of two treatment
# regimens for lung cancer.  Minimal argument call.  Print
# results, then plot error rate and importance values.

data(veteran, package = "randomSurvivalForest")
veteran.out <- rsf(Survrsf(time, status)~., data = veteran)
print.rsf(veteran.out)
plot.error(veteran.out)

# Example 2:  Richer argument call.
# Note that forest option is set to true to illustrate 
# how one might use 'rsf' for prediction (see 'rsf.predict'
# for more details).

data(veteran, package = "randomSurvivalForest")
veteran.f <- as.formula(Survrsf(time, status)~.)
ntree <- 500
mtry <- 2
nodesize <- 3
splitrule <- "logrank"
proximity <- TRUE
forest <- TRUE
seed <- -1
ntime <- NULL
do.trace <- 25
veteran2.out <- rsf(veteran.f, veteran, ntree, mtry, nodesize, splitrule,
                    proximity = proximity,
                    forest = forest,
                    seed = seed,
                    ntime = ntime,
                    do.trace = do.trace)
print.rsf(veteran2.out)
plot.proximity(veteran2.out)

# Take a peek at the forest ...
head(veteran2.out$forest$nativeArray)

# Partial plot of mortality for sorted predictors.
# Set partial to FALSE for a faster plot.
plot.variable(veteran2.out, 3, partial = TRUE)


# Example 3:  Veteran data again. Look specifically at
# Karnofsky performance score.  Compare to Kaplan-Meier.
# Assumes "survival" library is loaded.

if (library("survival", logical.return = TRUE))
{
        data(veteran, package = "randomSurvivalForest")
        veteran3.out <- rsf(Survrsf(time, status)~karno,
                       veteran,
                       ntree = 1000)
        plot.ensemble(veteran3.out)
        par(mfrow = c(1,1))
        plot(survfit(Surv(time, status)~karno, data = veteran))
}

# Example 4:  Primary biliary cirrhosis (PBC) of the liver.
# Data found in Appendix D.1 of Fleming and Harrington, Counting
# Processes and Survival Analysis, Wiley, 1991 (only differences
# are that age is in days and sex and stage variables are not
# missing for observations 313-418).  

data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000)
print(pbc.out)

# Example 5:  Compare Cox regression to Random Survival Forests
# for PBC data.  Compute OOB estimate of Harrell's concordance 
# index for Cox regression using B = 100 bootstrap draws.
# Assumes "Hmisc" and "survival" libraries are loaded. 

if (library("survival", logical.return = TRUE) 
    & library("Hmisc", logical.return = TRUE))
{
    data(pbc, package = "randomSurvivalForest")
    pbc2.out <- rsf(Survrsf(days,status)~., pbc, mtry = 2, ntree = 1000)
    B <- 100 
    cox.err <- rep(NA, B) 
    cox.f <- as.formula(Surv(days,status)~.)  
    pbc.data <- pbc[apply(is.na(pbc), 1, sum) == 0,] ##remove NA's 
    cat("Out-of-bag Cox Analysis ...", "\n")
    for (b in 1:B) {
        cat("Cox bootstrap", b, "\n") 
        bag.sample <- sample(1:dim(pbc.data)[1],
                             dim(pbc.data)[1],
                             replace = TRUE) 
        oob.sample <- setdiff(1:dim(pbc.data)[1], bag.sample)
        train <- pbc.data[bag.sample,]
        test <- pbc.data[oob.sample,]
        cox.out <- coxph(cox.f, train)
        cox.out <- tryCatch({coxph(cox.f, train)}, error=function(ex){NULL})
        if (is.list(cox.out)) {
          cox.predict <- predict(cox.out, test)
          cox.err[b] <- rcorr.cens(cox.predict, 
                 Surv(pbc.data$days[oob.sample], pbc.data$status[oob.sample]))[1]
        }
     }
     cat("Error rates:", "\n")
     cat("Random Survival Forests:", pbc2.out$err.rate[pbc2.out$ntree], "\n")
     cat("         Cox Regression:", mean(cox.err, na.rm = TRUE), "\n")
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
