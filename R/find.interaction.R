##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.0.1
##
##  Copyright 2007, Cleveland Clinic
##
##  This program is free software; you can redistribute it and/or
##  modify it under the terms of the GNU General Public License
##  as published by the Free Software Foundation; either version 2
##  of the License, or (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public
##  License along with this program; if not, write to the Free
##  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##  Boston, MA  02110-1301, USA.
##
##  Project funded by:
##    National Institutes of Health, HL072771-01
##
##    Michael S Lauer, MD, FACC, FAHA
##    Cleveland Clinic Lerner College of Medicine of CWRU
##    9500 Euclid Avenue
##    Cleveland, OH 44195
##
##    email:  lauerm@ccf.org
##    phone:   216-444-6798
##
##  Written by:
##    Hemant Ishwaran, Ph.D.
##    Dept of Quantitative Health Sciences/Wb4
##    Cleveland Clinic Foundation
##    9500 Euclid Avenue
##    Cleveland, OH 44195
##
##    email:  hemant.ishwaran@gmail.com
##    phone:  216-444-9932
##    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
##    --------------------------------------------------------------
##    Udaya B. Kogalur, Ph.D.
##    Kogalur Shear Corporation
##    5425 Nestleway Drive, Suite L1
##    Clemmons, NC 27012
##
##    email:  ubk2101@columbia.edu
##    phone:  919-824-9825
##    URL:    www.kogalur-shear.com
##
##**********************************************************************
##**********************************************************************

find.interaction <- function (
    x, 
    sort = TRUE,
    predictorNames = NULL,
    n.pred = NULL,
    n.rep  = 1,
    ...) {

    ### check that object is interpretable
    if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2)
      stop("Function only works for objects of class `(rsf, grow)'.")

    ### extract data, predictor names, splitrule, ntree, and formula for subsequent calls 
    ### special treatment needed if user passes predictorNames
    ### should predictors be sorted by importance?
    ### determine number of variables to be paired-up
    predictors <- x$predictors
    if (!is.null(x$imputedIndv)) predictors[x$imputedIndv, ] <- x$imputedData
    cov.names <- x$predictorNames
    splitrule <- x$splitrule
    ntree <- x$ntree
    n.interact <- n.cov <- length(cov.names)
    fNames <- all.vars(x$formula, max.names=1e7)
    rsf.f.org <- paste("Survrsf(",fNames[1], ",", fNames[2], ") ~", sep="")
    newData <- as.data.frame(cbind(x$time, x$cens, predictors))
    colnames(newData) <- c(fNames[1:2], cov.names)
    if (!is.null(predictorNames)) {
      if (sum(is.element(cov.names, predictorNames)) == 0) {
           cat("Coefficient list does not match available predictors:","\n")
           print(cov.names)
           stop()
      }
      predictorNames <- unique(predictorNames[is.element(predictorNames, cov.names)])
      if (length(predictorNames) == 1)
        stop("Pairwise comparisons require more than one candidate predictor")
      o.pt <- 1:(length(predictorNames))
      for (k in 1:length(predictorNames)) {
        o.pt[k] <- (1:n.cov)[cov.names == predictorNames[k]]
      }
      o.pt <- 2 + c(o.pt, setdiff(1:n.cov, o.pt))
      newData <- newData[, c(1, 2, o.pt)]
      cov.names <- names(newData)[-c(1, 2)]
      sort <- F
      n.pred <- NULL
      n.interact <- length(predictorNames)
    }
    if (sort) {
      if (!is.null(x$importance)) {
        o.r <- rev(order(x$importance))
        newData <- newData[, c(1,2,(3:(2+n.cov))[o.r])]
        cov.names <- cov.names[o.r]
      }
    }
    if (!is.null(n.pred)) {
      n.interact <- min(n.cov, max(round(n.pred), 1))
    }
    if (n.interact == 1) stop("Pairwise comparisons require more than one candidate predictor")

    ### interaction loop: call RSF multiple times to work out various error rates
    interact.imp <- rownames.interact.imp <- NULL
    rsf.f <- as.formula(paste(rsf.f.org,"."))
    for (k in 1:(n.interact-1)) {
      n.2.cov <- n.interact-k
      imp <- rep(0,1+n.2.cov)
      imp.2 <- rep(0,n.2.cov)
      for (l in (k+1):n.interact) {
        cat("Pairing",cov.names[k],"with",cov.names[l],"\n")
        rsf.k.f <- as.formula(paste(rsf.f.org,paste(cov.names[-k],collapse="+")))
        rsf.l.f <- as.formula(paste(rsf.f.org,paste(cov.names[-l],collapse="+")))
        rsf.2.f <- as.formula(paste(rsf.f.org,paste(cov.names[-c(k,l)],collapse="+")))
        for (m in 1:n.rep) {
          seed <- -1*sample(1:1e5,1)
          rsfOutput    <-  rsf(rsf.f,newData,ntree,splitrule=splitrule,seed=seed)
          rsfOutput.k  <-  rsf(rsf.k.f,newData,ntree,splitrule=splitrule,seed=seed)
          rsfOutput.l  <-  rsf(rsf.l.f,newData,ntree,splitrule=splitrule,seed=seed)
          rsfOutput.2  <-  rsf(rsf.2.f,newData,ntree,splitrule=splitrule,seed=seed)
          imp[1] <- imp[1]+rsfOutput.k$err[ntree]-rsfOutput$err[ntree]
          imp[l-k+1] <- imp[l-k+1]+rsfOutput.l$err[ntree]-rsfOutput$err[ntree]
          imp.2[l-k] <- imp.2[l-k]+rsfOutput.2$err[ntree]-rsfOutput$err[ntree]
        }
      }
      imp[1] <- imp[1]/n.2.cov
      imp <- imp/n.rep
      imp.2 <- imp.2/n.rep
      interact.imp <- rbind(interact.imp,cbind(imp.2,(imp[1]+imp)[-1],imp.2-(imp[1]+imp)[-1]))
      rownames.interact.imp <- c(rownames.interact.imp,
                                 paste(cov.names[k],":",cov.names[(k+1):n.interact],
                                 sep=""))
    }
    colnames(interact.imp) <- c("Paired","Additive","Difference")
    rownames(interact.imp) <- rownames.interact.imp

    ### output table
    cat("\n")
    cat("                    No. of variables: ", n.cov,               "\n", sep="")
    cat("                   Variables sorted?: ", sort,                "\n", sep="")
    cat("   No. of variables used for pairing: ", n.interact,          "\n", sep="")
    cat("    Total no. of paired interactions: ", dim(interact.imp)[1],"\n", sep="")
    cat("            Monte Carlo replications: ", n.rep,               "\n", sep="")
    cat("                      Splitting rule: ", x$splitrule,         "\n", sep="")
    cat("                     Number of trees: ", x$ntree,             "\n",sep="")
    cat("\n")
    print(round(interact.imp,4))

    ### return the goodies
    invisible(list(interaction.table=interact.imp))
  }
