##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 1.0.0
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

plot.ensemble <- function(x, main = deparse(substitute(x)), ...) {
    # internal checks
    if (!inherits(x,"randomSurvivalForest"))
        stop("This function only works for objects of class `randomSurvivalForest'")

    # survival curves
    survival.ensemble <- t(exp(-x$ensemble))
    survival.mean.ensemble <- apply(survival.ensemble,1,mean)
    allTime.unique <- sort(unique(x$Time[x$Died == 1]))
    unique.pt <- is.element(allTime.unique, x$Time.unique)
    Y <- apply(cbind(1:length(allTime.unique)),
               1,
               function(j, tau, t.unq) {sum(tau >= t.unq[j])},
               tau = x$Time,
               t.unq = allTime.unique)
    d <- apply(cbind(1:length(allTime.unique)),
               1,
               function(j, d, tau, t.unq) {sum(tau == t.unq[j])},
               tau = x$Time[x$Died == 1],
               t.unq = allTime.unique)
    r <- d/(Y+1*(Y == 0))
    survival.aalen <- exp(-cumsum(r))[unique.pt]

    # Brier score
    # for all distinct time points
    # stratified by mortality percentiles
    brier.score <- NULL
    mort.perc <- c(min(x$mortality)-1e-5, quantile(x$mortality, (1:4)/4))
    delta <- t(apply(cbind(1:x$n),
                   1,
                   function(i, tau, t.unq) {1*(tau[i] > t.unq)},
                   tau =  x$Time,
                   t.unq = x$Time.unique)  - survival.ensemble)^2
    for (k in 1:4){
      mort.pt <- (x$mortality > mort.perc[k]) & (x$mortality <= mort.perc[k+1]) 
      brier.score <- cbind(brier.score,apply(delta[mort.pt,], 2, mean))
    }                  

    # plots
    old.cex <- par()$cex
    par(mfrow = c(2,2))
    par(cex = 1.0)
    if (x$n > 500) {
        r.pt <- sample(1:x$n, 500, replace = FALSE)
        matplot(x$Time.unique,
                survival.ensemble[,r.pt],
                xlab = "Time",
                ylab = "Survival",
                type = "l",
                col = 1, 
                lty = 3)
    }
    else {
        matplot(x$Time.unique,
                survival.ensemble,
                xlab = "Time",
                ylab = "Survival",
                type = "l",
                col = 1,
                lty = 3)
    }
    lines(x$Time.unique, survival.aalen, lty = 1, col = 3, lwd = 3)
    lines(x$Time.unique, survival.mean.ensemble, lty = 1, col = 2, lwd = 3)
    title("Ensemble Survival", cex.main = 1.25)
    plot(survival.aalen,
         survival.mean.ensemble,
         xlab = "Nelson-Aalen Survival",
         ylab = "Ensemble Survival",
         type = "l")
    abline(0, 1, col = 2, lty = 2)
    title("Survival", cex.main = 1.25)
    matplot(x$Time.unique, brier.score,
         xlab = "Time",
         ylab = "Score",
         type = "l",
         col  = 1,
         lty  = 1:4)
    lines(x$Time.unique, apply(delta, 2, mean), col=3, lwd=1)
    point.x=round(length(x$Time.unique)*c(3,4)/4)
    text(x$Time.unique[point.x],brier.score[point.x,1],"0-25",col=4)
    text(x$Time.unique[point.x],brier.score[point.x,2],"25-50",col=4)
    text(x$Time.unique[point.x],brier.score[point.x,3],"50-75",col=4)
    text(x$Time.unique[point.x],brier.score[point.x,4],"75-100",col=4)
    title("Brier Score",cex.main = 1.25)
    plot(x$Time, x$mortality, xlab = "Time", ylab = "Ensemble Mortality", type = "n")
    title("Mortality vs Time", cex.main = 1.25)
    if (x$n > 500) cex <- 0.5 else cex <- 0.75
    points(x$Time[x$Died == 1], x$mortality[x$Died == 1], pch = 16, col = 4, cex = cex)
    points(x$Time[x$Died == 0], x$mortality[x$Died == 0], pch = 16, cex = cex)
    points(supsmu(x$Time[x$Died == 1][order(x$Time[x$Died == 1])],
                  x$mortality[x$Died == 1][order(x$Time[x$Died == 1])]),
           type = "l",
           lwd = 2,
           col = 4,
           cex = cex)
    if (sum(x$Died == 0) > 1)
    points(supsmu(x$Time[x$Died == 0][order(x$Time[x$Died == 0])],
                  x$mortality[x$Died == 0][order(x$Time[x$Died == 0])]),
           type = "l",
           lwd = 2,
           cex = cex)
    par(cex = old.cex)
}
