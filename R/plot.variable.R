##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 2.0.0
##
##  Copyright 2006, Cleveland Clinic
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

plot.variable <- function (
    x, 
    plots.per.page = 4,
    granule = 5,
    sort = TRUE,
    partial = FALSE,
    predictorNames = NULL,
    n.pred = NULL,                           
    n.pts = 25,
    ...) {

    ### check that object is interpretable
    if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
        sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
      stop("Function only works for objects of class `(rsf, grow)', '(rsf, predict)'.")

    ### extract predictor names to be plotted
    ### should predictors be sorted by importance?
    predictors <- x$predictors
    if (is.null(predictorNames)) {
      cov.names <- x$predictorNames
    }
    else {
      cov.names <- predictorNames
      if (sum(is.element(x$predictorNames, cov.names)) == 0){
           cat("Coefficient list does not match available predictors:","\n")
           print(x$predictorNames)
           stop()
         }
      cov.names <- cov.names[is.element(cov.names, x$predictorNames)]
    }
    if (sort) {
      if (!is.null(x$importance)) {
        n.cov <- length(cov.names)
        cov.imp <- rep(0, n.cov)
        for (k in 1:n.cov) {
          cov.imp[k] <- x$importance[x$predictorNames == cov.names[k]]
        }
        cov.names <- cov.names[rev(order(cov.imp))]
      }
    }
    if (!is.null(n.pred)) {
      n.pred <- max(round(n.pred), 1)
      n.cov <- min(length(cov.names), n.pred)
      cov.names <- cov.names[1:n.cov]
    }
    n <- dim(predictors)[1]
    mortality.ensemble <- x$mortality

      
    ### plots (marginal, partial)
    old.par <- par(no.readonly = TRUE)
    if (!partial) {
      if (n > 500) cex <- 0.5 else cex <- 0.75
      plots.per.page <- max(round(min(plots.per.page,n.cov)), 1)
      granule <- max(round(granule),1)
      par(mfrow = c(min(plots.per.page, ceiling(n.cov/plots.per.page)), plots.per.page))
      for (k in 1:n.cov) {
        y <- predictors[, x$predictorNames == cov.names[k]]
        y.uniq <- unique(y)
        n.y <- length(y.uniq)
        if (n.y > granule) {
            plot(y,
                 mortality.ensemble,
                 xlab=cov.names[k],
                 ylab = "",
                 type = "n",
                 cex.lab = 1.5)
            points(y[x$Cens == 1], mortality.ensemble[x$Cens == 1],
                   pch = 16, col = 4, cex = cex)
            points(y[x$Cens==0], mortality.ensemble[x$Cens == 0],
                   pch = 16, cex = cex)
            lines(lowess(y, mortality.ensemble), col = 2)
        }
        else {
            boxplot(mortality.ensemble ~ y,
                    xlab = cov.names[k],
                    notch = T,
                    outline = F,
                    data = predictors,
                    col = "bisque",
                    names = rep("", n.y),
                    xaxt = "n",
                    pars = list(cex.lab = 1.5))
            at.pretty <- unique(round(pretty(1:n.y, min(30, n.y))))
            at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.y]
            axis(1,
                 at = at.pretty,
                 labels = format(sort(y.uniq)[at.pretty], trim = T, digits = 4),
                 tick = T)
        }
      }
    }
    else {
      if (is.null(x$forest)) {
        stop("Forest is empty!  Re-run rsf (grow) analysis with forest set to 'TRUE'.")
      }
      plots.per.page <- max(round(min(plots.per.page,n.cov)), 1)
      granule <- max(round(granule),1)
      par(mfrow = c(min(plots.per.page, ceiling(n.cov/plots.per.page)), plots.per.page))
      baseForest <- x$forest
      class(baseForest) <- c("rsf", "partial")
      if (n.pts < 1) n.pts <- 1 else n.pts <- round(n.pts)
      for (k in 1:n.cov) {
        y <- predictors[, x$predictorNames == cov.names[k]]
        n.y <- length(unique(y))
        if (n.y > n.pts) {
          y.uniq <- sort(unique(y))[unique(as.integer(seq(1, n.y, length = min(n.pts, n.y))))]
        }
        else {
           y.uniq <- sort(unique(y))
        }
        n.y <- length(y.uniq)
        if (n.y > 25) cex <- 0.5 else cex <- 0.75
        mortality.y <- mortality.y.se <- NULL
        newdata.y <- as.data.frame(predictors)
        colnames(newdata.y) <- x$predictorNames
        for (l in 1:n.y) {
           newdata.y[, x$predictorNames == cov.names[k]] <- rep(y.uniq[l], n)
           if (n.y > granule | n.cov == 1) {
             mortality.y <- c(mortality.y, mean(predict.rsf(baseForest, newdata.y)$mortality))
             mortality.y.se <- c(mortality.y.se,
                                 sd(predict.rsf(baseForest, newdata.y)$mortality)/sqrt(n))
           }
           else {
             pred.temp <- predict.rsf(baseForest, newdata.y)$mortality
             mean.temp <- mean(pred.temp)
             pred.temp <- mean.temp + (pred.temp-mean.temp)/sqrt(n)
             mortality.y <- c(mortality.y, pred.temp)
           }
        }
        if (n.y > granule | n.cov == 1) {
            plot(c(min(y), y.uniq,max(y), y.uniq, y.uniq),
                 c(NA, mortality.y, NA, mortality.y+2*mortality.y.se, mortality.y-2*mortality.y.se),
                 xlab=cov.names[k],
                 ylab = "",
                 type = "n",
                 cex.lab = 1.5)
            points(y.uniq, mortality.y, pch = 16, cex = cex, col = 2)
            if (any(mortality.y.se > 0)) {
              lines(lowess(y.uniq, mortality.y+2*mortality.y.se), lty = 3, col = 2)
              lines(lowess(y.uniq, mortality.y-2*mortality.y.se), lty = 3, col = 2)
            }
            lines(lowess(y.uniq, mortality.y), lty = 2)
            rug(y, ticksize=0.03)
        }
        else {
            bxp.call <- boxplot(mortality.y ~ rep(y.uniq, rep(n, n.y)), range = 2, plot = F)
            boxplot(mortality.y ~ rep(y.uniq, rep(n, n.y)),
                    xlab = cov.names[k],
                    notch = T,
                    outline = F,
                    range = 2,
                    ylim = c(min(bxp.call$stats[1,])-2,max(bxp.call$stats[5,])+2),
                    data = predictors,
                    col = "bisque",
                    names = rep("",n.y),
                    xaxt = "n",
                    pars = list(cex.lab = 1.5))
            at.pretty <- unique(round(pretty(1:n.y, min(30,n.y))))
            at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.y]
            axis(1,
                 at = at.pretty,
                 labels = format(sort(y.uniq)[at.pretty], trim = T, digits = 4),
                 tick = T)
        }
      }
    }
    par(old.par)
  }
