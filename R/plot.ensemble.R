####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.4
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  Written by:
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  commerce@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

plot.ensemble.rsf <- function (x, plots.one.page = TRUE, ...) {
    if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
      stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")
    if (sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) == 2) {
      rsfPred <- TRUE
    }
    else {
      rsfPred <- FALSE
    }
    if (is.null(x$ndead)) return()
    if (!is.null(x$imputedIndv)) {
      x$cens[x$imputedIndv]=x$imputedData[,1]
      x$time[x$imputedIndv]=x$imputedData[,2]
    }
    if (x$n < 5 | sum(x$cens) < 2) return()
    if (rsfPred) {
      mort    <- x$mortality
      ensb    <- x$ensemble
      y.lab   <- "Mortality"
      title.1 <- "Ensemble Survival"
      title.2 <- "Mortality vs Time"
    }
    else {
      mort    <- x$oob.mortality
      ensb    <- x$oob.ensemble
      y.lab   <- "OOB Mortality"
      title.1 <- "OOB Ensemble Survival"
      title.2 <- "OOB Mortality vs Time"
    }
    n.event  <- length(unique(na.omit(x$cens)[na.omit(x$cens) > 0]))
    if (n.event > 1) ensb <- ensb[,,1]
    surv.ensb <- t(exp(-ensb))
    surv.mean.ensb <- apply(surv.ensb, 1, mean, na.rm = TRUE)
    if (!rsfPred) {
      Y <- sapply(1:length(x$timeInterest),
                 function(j, tau, t.unq) {sum(tau >= t.unq[j])},
                 tau = x$time,
                 t.unq = x$timeInterest)
      d <- sapply(1:length(x$timeInterest),
                 function(j, d, tau, t.unq) {sum(tau == t.unq[j])},
                 tau = x$time[x$cens != 0],
                 t.unq = x$timeInterest)
      r <- d/(Y+1*(Y == 0))
      surv.aalen <- exp(-cumsum(r))
      sIndex <- function(Ju,Ev) { sapply(1:length(Ev), function(j) {sum(Ju <= Ev[j])}) }
      censTime <- sort(unique(x$time[x$cens == 0]))
      censTime.pt <- sIndex(censTime, x$timeInterest)
      Y <- sapply(1:length(censTime),
                 function(j, tau, t.unq) {sum(tau >= t.unq[j])},
                 tau = x$time,
                 t.unq = censTime)
      d <- sapply(1:length(censTime),
                 function(j, d, tau, t.unq) {sum(tau == t.unq[j])},
                 tau = x$time[x$cens == 0],
                 t.unq = censTime)
      r <- d/(Y+1*(Y == 0))
      cens.aalen <- c(1, exp(-cumsum(r)))[1+censTime.pt]
    }
    if (!rsfPred) {
      brier.wt <- t(apply(cbind(1:x$n),
                   1,
                   function(i, tau, event, t.unq) {
                     pt <- sIndex(t.unq, tau[i])
                     c1 <- 1*(tau[i] <= t.unq & event[i] != 0)/c(1,cens.aalen)[1+pt]
                     c2 <- 1*(tau[i] > t.unq)/cens.aalen
                     (c1 + c2)
                   },
                   tau =  x$time, event = x$cens,
                   t.unq = x$timeInterest))
      delta <- t(apply(cbind(1:x$n),
                   1,
                   function(i, tau, event, t.unq) {1*(tau[i] > t.unq)},
                   tau =  x$time, event = x$cens,
                   t.unq = x$timeInterest)  -  surv.ensb)^2
      brier.score <- matrix(NA, length(x$timeInterest), 4)
      mort.perc   <- c(min(mort, na.rm = TRUE) - 1e-5, quantile(mort, (1:4)/4, na.rm = TRUE))
      for (k in 1:4){
        mort.pt <- (mort > mort.perc[k]) & (mort <= mort.perc[k+1])
        brier.score[, k] <- sapply(1:length(x$timeInterest),
                                  function(j) {mean(brier.wt[mort.pt, j]*delta[mort.pt, j], na.rm = TRUE)})
      }
      brier.score <- cbind(brier.score, sapply(1:length(x$timeInterest),
                                  function(j) {mean(brier.wt[, j]*delta[, j], na.rm = TRUE)}))
    }
    old.par <- par(no.readonly = TRUE)
    if (plots.one.page) {
      if (!rsfPred) par(mfrow = c(2,2)) else par(mfrow = c(1,2))
    }
    else {
      par(mfrow=c(1,1))
    }
    par(cex = 1.0)
    if (x$n > 500) {
        r.pt <- sample(1:x$n, 500, replace = FALSE)
        matplot(x$timeInterest,
                surv.ensb[,r.pt],
                xlab = "Time",
                ylab = title.1,
                type = "l",
                col = 1, 
                lty = 3)
    }
    else {
        matplot(x$timeInterest,
                surv.ensb,
                xlab = "Time",
                ylab = title.1,
                type = "l",
                col = 1,
                lty = 3)
    }
    if (!rsfPred) lines(x$timeInterest, surv.aalen, lty = 1, col = 3, lwd = 3)
    lines(x$timeInterest, surv.mean.ensb, lty = 1, col = 2, lwd = 3)
    rug(x$timeInterest, ticksize=-0.03)
    if (plots.one.page) title(title.1, cex.main = 1.25) 
    if (!rsfPred) {
      plot(surv.aalen,
           surv.mean.ensb,
           xlab = "Nelson-Aalen Survival",
           ylab = title.1,
           type = "l")
      abline(0, 1, col = 2, lty = 2)
      if (plots.one.page) title("Survival", cex.main = 1.25)
    }
    if (!rsfPred) {
      matplot(x$timeInterest, brier.score,
         xlab = "Time",
         ylab = "Score",
         type = "l",
         lwd  = c(rep(1, 4), 2),
         col  = c(rep(1, 4), 2),
         lty  = c(1:4, 1))
      point.x=round(length(x$timeInterest)*c(3,4)/4)
      text(x$timeInterest[point.x],brier.score[point.x,1],"0-25",col=4)
      text(x$timeInterest[point.x],brier.score[point.x,2],"25-50",col=4)
      text(x$timeInterest[point.x],brier.score[point.x,3],"50-75",col=4)
      text(x$timeInterest[point.x],brier.score[point.x,4],"75-100",col=4)
      rug(x$timeInterest, ticksize=0.03)
      if (plots.one.page) title("Brier Score",cex.main = 1.25)
    }
    plot(x$time, mort, xlab = "Time", ylab = y.lab, type = "n")
    if (plots.one.page) title(title.2, cex.main = 1.25)
    if (x$n > 500) cex <- 0.5 else cex <- 0.75
    points(x$time[x$cens == 1], mort[x$cens == 1], pch = 16, col = 4, cex = cex)
    points(x$time[x$cens == 0], mort[x$cens == 0], pch = 16, cex = cex)
    if (sum(x$cens == 1) > 1)
      points(supsmu(x$time[x$cens == 1][order(x$time[x$cens == 1])],
                  mort[x$cens == 1][order(x$time[x$cens == 1])]),
             type = "l",
             lty = 3,
             col = 4,
             cex = cex)
    if (sum(x$cens == 0) > 1)
      points(supsmu(x$time[x$cens == 0][order(x$time[x$cens == 0])],
                  mort[x$cens == 0][order(x$time[x$cens == 0])]),
             type = "l",
             lty = 3,
             cex = cex)
    rug(x$timeInterest, ticksize=-0.03)
    par(old.par)      
}
plot.ensemble <- plot.ensemble.rsf
