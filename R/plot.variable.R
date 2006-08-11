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

plot.variable <- function(x, 
    plots.per.page = 4,
    granule = 25,
    prednames = NULL,
    main = deparse(substitute(x)),
    ...)
{
    if (!inherits(x,"randomSurvivalForest"))
        stop("This function only works for objects of class `randomSurvivalForest'")
    predictors <- x$predictors
    if (is.null(prednames)) {
      cov.names = x$prednames
    }
    else {
      cov.names = prednames
      if (sum(is.element(x$prednames, cov.names)) == 0){
           cat("coefficient list does not match available predictors:","\n")
           print(x$prednames)
           stop()
         }
      cov.names = cov.names[is.element(cov.names, x$prednames)]
    }
    n.cov <- length(cov.names)
    mortality.ensemble <- x$mortality
    old.cex <- par()$cex
    if (x$n > 500) cex <- 0.5 else cex <- 0.75
    plots.per.page <- max(round(min(plots.per.page,n.cov)),1)
    granule <- max(round(granule),1)
    par(mfrow = c(min(plots.per.page, ceiling(n.cov/plots.per.page)), plots.per.page))
    for (k in 1:n.cov) {
        y <- predictors[, x$prednames == cov.names[k]]
        n.y <- length(unique(y))
        if (n.y > granule) {
            plot(y, mortality.ensemble, xlab=cov.names[k], ylab = "", type = "n", cex.lab = 1.5)
            points(y[x$Died == 1], mortality.ensemble[x$Died == 1], pch = 16, col = 4, cex = cex)
            points(y[x$Died==0], mortality.ensemble[x$Died == 0], pch = 16, cex = cex)
            lines(lowess(y, mortality.ensemble), col = 2)
        }
        else {
            boxplot(mortality.ensemble~y,
                    xlab = cov.names[k],
                    notch = T,
                    outline = F,
                    data = predictors,
                    col = "bisque",
                    names = rep("",n.y),
                    xaxt = "n",
                    pars = list(cex.lab = 1.5))
           x.at <- unique(round(pretty(1:n.y, min(30,n.y))))
           x.at <- x.at[x.at >= 1 & x.at <= n.y]
           axis(1, at = x.at, labels = sort(unique(y))[x.at], tick = T)
        }
    }
    par(cex = old.cex)
}
