##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.5.1
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

plot.ensemble <- function (x, plots.one.page = TRUE, ...) {

    ### check that object is interpretable
    if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
      stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")
    if (sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) == 2) {
      rsfPred <- TRUE
    }
    else {
      rsfPred <- FALSE
    }
    
    ### null case can occur for '(rsf, predict)' objects, so check 
    if (is.null(x$ndead)) return()

    ### use imputed missing time or censoring indicators
    if (!is.null(x$imputedIndv)) {
      x$cens[x$imputedIndv]=x$imputedData[,1]
      x$time[x$imputedIndv]=x$imputedData[,2]
    }
    
    ### no point in producing plots if sample size too small or
    ### not enough deaths
    if (x$n < 5 | sum(x$cens) < 2) return()

    # survival curves
    survival.ensemble <- t(exp(-x$ensemble))
    survival.mean.ensemble <- apply(survival.ensemble,1,mean)
    if (!rsfPred) {
      allTime.unique <- sort(unique(x$time[x$cens == 1]))
      unique.pt <- is.element(allTime.unique, x$timeInterest)
      Y <- apply(cbind(1:length(allTime.unique)),
                 1,
                 function(j, tau, t.unq) {sum(tau >= t.unq[j])},
                 tau = x$time,
                 t.unq = allTime.unique)
      d <- apply(cbind(1:length(allTime.unique)),
                 1,
                 function(j, d, tau, t.unq) {sum(tau == t.unq[j])},
                 tau = x$time[x$cens == 1],
                 t.unq = allTime.unique)
      r <- d/(Y+1*(Y == 0))
      survival.aalen <- exp(-cumsum(r))[unique.pt]
    }

    # Brier score
    # . all distinct time points
    # . stratified by mortality percentiles
    brier.score <- matrix(NA, length(x$timeInterest), 4)
    mort.perc <- c(min(x$mortality)-1e-5, quantile(x$mortality, (1:4)/4))
    brier.wt <- c(apply(cbind(1:x$n),
                   1,
                   function(i, tau, event, t.unq) {
                     1*(tau[i] > t.unq) + 1*(tau[i] <= t.unq & event[i] == 1)
                   },
                   tau =  x$time, event = x$cens,
                   t.unq = x$timeInterest))
    delta <- t(apply(cbind(1:x$n),
                   1,
                   function(i, tau, event, t.unq) {1*(tau[i] > t.unq)},
                   tau =  x$time, event = x$cens,
                   t.unq = x$timeInterest)  - survival.ensemble)^2
    for (k in 1:4){
      mort.pt <- (x$mortality > mort.perc[k]) & (x$mortality <= mort.perc[k+1])
      brier.score[,k] <- apply(as.matrix(delta[mort.pt,]), 2, function(x) {
                         sum(brier.wt[mort.pt]*x)/(1*(sum(brier.wt[mort.pt])==0)
                                 +sum(brier.wt[mort.pt]))})
    }                  

    # plots
    old.par <- par(no.readonly = TRUE)
    if (plots.one.page) par(mfrow = c(2,2)) else par(mfrow=c(1,1))
    par(cex = 1.0)
    if (x$n > 500) {
        r.pt <- sample(1:x$n, 500, replace = FALSE)
        matplot(x$timeInterest,
                survival.ensemble[,r.pt],
                xlab = "Time",
                ylab = "Survival",
                type = "l",
                col = 1, 
                lty = 3)
    }
    else {
        matplot(x$timeInterest,
                survival.ensemble,
                xlab = "Time",
                ylab = "Survival",
                type = "l",
                col = 1,
                lty = 3)
    }
    if (!rsfPred) lines(x$timeInterest, survival.aalen, lty = 1, col = 3, lwd = 3)
    lines(x$timeInterest, survival.mean.ensemble, lty = 1, col = 2, lwd = 3)
    rug(x$timeInterest, ticksize=-0.03)
    if (plots.one.page) title("Ensemble Survival", cex.main = 1.25) 
    if (!rsfPred) {
      plot(survival.aalen,
           survival.mean.ensemble,
           xlab = "Nelson-Aalen Survival",
           ylab = "Ensemble Survival",
           type = "l")
      abline(0, 1, col = 2, lty = 2)
      if (plots.one.page) title("Survival", cex.main = 1.25)
    }
    matplot(x$timeInterest, brier.score,
         xlab = "Time",
         ylab = "Score",
         type = "l",
         col  = 1,
         lty  = 1:4)
    lines(x$timeInterest, apply(delta, 2, mean), col=2, lwd=3)
    point.x=round(length(x$timeInterest)*c(3,4)/4)
    text(x$timeInterest[point.x],brier.score[point.x,1],"0-25",col=4)
    text(x$timeInterest[point.x],brier.score[point.x,2],"25-50",col=4)
    text(x$timeInterest[point.x],brier.score[point.x,3],"50-75",col=4)
    text(x$timeInterest[point.x],brier.score[point.x,4],"75-100",col=4)
    rug(x$timeInterest, ticksize=0.03)
    if (plots.one.page) title("Brier Score",cex.main = 1.25)
    plot(x$time, x$mortality, xlab = "Time", ylab = "Ensemble Mortality", type = "n")
    if (plots.one.page) title("Mortality vs Time", cex.main = 1.25)
    if (x$n > 500) cex <- 0.5 else cex <- 0.75
    points(x$time[x$cens == 1], x$mortality[x$cens == 1], pch = 16, col = 4, cex = cex)
    points(x$time[x$cens == 0], x$mortality[x$cens == 0], pch = 16, cex = cex)
    if (sum(x$cens == 1) > 1)
      points(supsmu(x$time[x$cens == 1][order(x$time[x$cens == 1])],
                  x$mortality[x$cens == 1][order(x$time[x$cens == 1])]),
             type = "l",
             lty = 3,
             col = 4,
             cex = cex)
    if (sum(x$cens == 0) > 1)
      points(supsmu(x$time[x$cens == 0][order(x$time[x$cens == 0])],
                  x$mortality[x$cens == 0][order(x$time[x$cens == 0])]),
             type = "l",
             lty = 3,
             cex = cex)
    rug(x$timeInterest, ticksize=-0.03)
    par(old.par)      
}
