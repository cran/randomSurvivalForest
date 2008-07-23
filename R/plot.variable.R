##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.5.0
##
##  Copyright 2008, Cleveland Clinic Foundation
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

plot.variable <- function(
    x, 
    plots.per.page = 4,
    granule = 5,
    sorted = TRUE,
    type = c("mort", "rel.freq", "surv", "time")[1],
    partial = FALSE,
    predictorNames = NULL,
    npred = NULL,                           
    npts = 25,
    subset = NULL, 
    ...) {

    ### don't want to use x for object 
    ### check that object is interpretable
    object <- x
    rm(x)
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
        sum(inherits(object, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
      stop("Function only works for objects of class `(rsf, grow)', '(rsf, predict)'.")
    if (type != "mort" & type != "rel.freq" & type != "surv" & type != "time")
      stop("Invalid choice for 'type:  " + type)
    if (!partial & type == "time") stop("Type 'time' can only be used for partial plots.")
    
    ### subset the data?
    if (!is.null(subset) & length(unique(subset)) != 0) {
      subset <- subset[subset >=1 & subset <= dim(object$predictors)[1]]
      subset <- unique(subset)
      if (length(subset) == 0) stop("'subset' not set properly.")
      if (length(subset) > 1) {
        object$predictors <- object$predictors[subset, ]
      }
      else {
        object$predictors <- t(as.matrix(object$predictors[subset, ]))
      }
      object$mortality <- object$mortality[subset]
      if (sum(subset) > 1) {
        object$ensemble <- object$ensemble[subset, ]
      }
      else {
        object$ensemble <- t(as.matrix(object$ensemble[subset, ]))
      }
    }
   
    ### get predictor matrix (use imputed values if available)
    ### extract predictor names to be plotted
    ### should predictors be sorted by importance?
    predictors <- object$predictors
    if (!is.null(object$imputedIndv)) predictors[object$imputedIndv, ] <- object$imputedData[, -c(1:2)]
    if (is.null(predictorNames)) {
      cov.names <- object$predictorNames
    }
    else {
      cov.names <- predictorNames
      if (sum(is.element(object$predictorNames, cov.names)) == 0){
           cat("Coefficient list does not match available predictors:","\n")
           print(object$predictorNames)
           stop()
      }
      cov.names <- unique(cov.names[is.element(cov.names, object$predictorNames)])
      npred <- NULL
    }
    n.cov <- length(cov.names)
    if (sorted) {
      if (!is.null(object$importance)) {
        n.cov <- length(cov.names)
        cov.imp <- rep(0, n.cov)
        for (k in 1:n.cov) {
          cov.imp[k] <- object$importance[object$predictorNames == cov.names[k]]
        }
        cov.names <- cov.names[rev(order(cov.imp))]
      }
    }
    if (!is.null(npred)) {
      npred <- max(round(npred), 1)
      n.cov <- min(length(cov.names), npred)
      cov.names <- cov.names[1:n.cov]
    }
    n <- dim(predictors)[1]
    

    ## Save par settings
    old.par <- par(no.readonly = TRUE)


    ## nice y-label
    if (type == "mort") {
      ylabel <- "mortality"
    }
    else if (type == "rel.freq") {
      ylabel <- "standardized mortality"
    }
    else if (type == "surv") {
      ylabel <- "predicted median survival"
    }
    else {
      ylabel <- "predicted survival time"
    }
    
    ##--------------------------------------------------------------------------------
    ## Marginal plots
    ##--------------------------------------------------------------------------------
    if (!partial) {
      if (n > 500) cex <- 0.5 else cex <- 0.75
      plots.per.page <- max(round(min(plots.per.page,n.cov)), 1)
      granule <- max(round(granule),1)
      par(mfrow = c(min(plots.per.page, ceiling(n.cov/plots.per.page)), plots.per.page))
      if (type == "mort") {
        yhat <- object$mortality
      }
      else if (type == "rel.freq") {
        yhat <- object$mortality/max(n, na.omit(object$mortality))
      }
      else {
        yhat <-
          100*exp(-object$ensemble[ , max(which(object$timeInterest<=na.omit(median(object$time))))])
      }
      for (k in 1:n.cov) {
        x <- predictors[, object$predictorNames == cov.names[k]]
        x.uniq <- unique(x)
        n.x <- length(x.uniq)
        if (!is.factor(x) & n.x > granule) {
            plot(x,
                 yhat,
                 xlab = cov.names[k],
                 ylab = ylabel,
                 type = "n",
                 cex.lab = 1.5)
            points(x[object$cens == 1], yhat[object$cens == 1],
                   pch = 16, col = 4, cex = cex)
            points(x[object$cens==0], yhat[object$cens == 0],
                   pch = 16, cex = cex)
            lines(lowess(x[!is.na(x)], yhat[!is.na(x)]), col = 2, lwd=3)
        }
        else {
          if (is.factor(x)) x <- factor(x, exclude = NULL)          
          boxplot(yhat ~ x, na.action = "na.omit",
                    xlab = cov.names[k],
                    ylab = ylabel,
                    notch = TRUE,
                    outline = FALSE,
                    data = predictors,
                    col = "bisque",
                    names = rep("", n.x),
                    xaxt = "n",
                    pars = list(cex.lab = 1.5))
          at.pretty <- unique(round(pretty(1:n.x, min(30, n.x))))
          at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
          axis(1,
               at = at.pretty,
               labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
               tick = TRUE)
        }
      }
    }
    ##--------------------------------------------------------------------------------
    ## Partial plots
    ##--------------------------------------------------------------------------------
    else {
      if (is.null(object$forest)) {
        stop("Forest is empty!  Re-run rsf (grow) analysis with forest set to 'TRUE'.")
      }
      plots.per.page <- max(round(min(plots.per.page,n.cov)), 1)
      granule <- max(round(granule),1)
      par(mfrow = c(min(plots.per.page, ceiling(n.cov/plots.per.page)), plots.per.page))
      baseForest <- object$forest
      if (type == "time") {
        class(baseForest) <- c("rsf", "partial.rough")
      }
      else { 
        class(baseForest) <- c("rsf", "partial")
      }
      if (npts < 1) npts <- 1 else npts <- round(npts)
      for (k in 1:n.cov) {
        x <- predictors[, object$predictorNames == cov.names[k]]
        if (is.factor(x)) x <- factor(x, exclude = NULL)          
        n.x <- length(unique(x))
        if (!is.factor(x) & n.x > npts) {
          x.uniq <- sort(unique(x))[unique(as.integer(seq(1, n.x, length = min(npts, n.x))))]
        }
        else {
           x.uniq <- sort(unique(x))
        }
        n.x <- length(x.uniq)
        if (n.x > 25) cex <- 0.5 else cex <- 0.75
        yhat <- yhat.se <- NULL
        newdata.x <- predictors
        for (l in 1:n.x) {
          newdata.x[, object$predictorNames == cov.names[k]] <- rep(x.uniq[l], n)
          if (type == "mort" | type == "rel.freq" | type == "time") {
            pred.temp <- predict.rsf(baseForest, newdata.x)$mortality
          }
          else if (type == "surv") {
            pred.temp <-
             100*exp(-predict.rsf(baseForest, newdata.x)$ensemble[,
                    max(which(object$timeInterest<=median(na.omit(object$time))))])
          }
          if (!is.factor(x) & (n.x > granule | n.cov == 1)) {
            yhat <- c(yhat, mean(pred.temp , na.rm = TRUE))
            yhat.se <- c(yhat.se, sd(pred.temp/sqrt(n) , na.rm = TRUE))
          }
          else {
            mean.temp <- mean(pred.temp , na.rm = TRUE)
            pred.temp <- mean.temp + (pred.temp-mean.temp)/sqrt(n)
            yhat <- c(yhat, pred.temp)
          }
        }        
        if (type == "rel.freq") nAdj <- max(n, yhat, na.rm = TRUE)
        if (!is.factor(x) & (n.x > granule | n.cov == 1)) {
          if (type == "rel.freq") {
            yhat <- yhat/nAdj
            yhat.se <- yhat.se/n
          }
          plot(c(min(x), x.uniq, max(x), x.uniq, x.uniq),
               c(NA, yhat, NA, yhat+2*yhat.se, yhat-2*yhat.se),
               xlab = cov.names[k],
               ylab = ylabel,
               type = "n",
               cex.lab = 1.5)
          points(x.uniq, yhat, pch = 16, cex = cex, col = 2)
          if (!is.na(yhat.se) && any(yhat.se > 0)) {
            lines(lowess(x.uniq, yhat+2*yhat.se), lty = 3, col = 2)
            lines(lowess(x.uniq, yhat-2*yhat.se), lty = 3, col = 2)
          }
          lines(lowess(x.uniq, yhat), lty = 2, lwd=2)
          rug(x, ticksize=0.03)
        }
        else {
          if (type != "rel.freq") {
            y.se <- 2
          }
          else {
            yhat <- yhat/nAdj
            y.se <- 2/n
          }          
          bxp.call <- boxplot(yhat ~ rep(x.uniq, rep(n, n.x)), range = 2, plot = F)
          boxplot(yhat ~ rep(x.uniq, rep(n, n.x)),
                  xlab = cov.names[k],
                  ylab = ylabel,
                  notch = TRUE,
                  outline = FALSE,
                  range = 2,
                  ylim = c(min(bxp.call$stats[1,], na.rm=TRUE)
                      - y.se,max(bxp.call$stats[5,], na.rm=TRUE) + y.se),
                  data = predictors,
                  col = "bisque",
                  names = rep("",n.x),
                  xaxt = "n",
                  pars = list(cex.lab = 1.5))
          at.pretty <- unique(round(pretty(1:n.x, min(30,n.x))))
          at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
          axis(1,
               at = at.pretty,
               labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
               tick = TRUE)
        }
      }
    }

    ## Restore par settings
    par(old.par)


  }
