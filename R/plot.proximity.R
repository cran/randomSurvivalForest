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

plot.proximity.rsf <- function (x, plot = TRUE, ...) {
   if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'")
  if (is.null(x$proximity))
    stop("Proximity is NULL.  Re-run analysis with proximity set to 'TRUE'.")
  if (length(x$proximity) == 1) return()
  proximity <- matrix(0, x$n, x$n)
  count <- 0
  for (k in 1:x$n){
    proximity[k,1:k] <- x$proximity[(count+1):(count+k)]
    proximity[1:k,k] <- proximity[k,1:k]
    count <- count+k
  }
  proximity <- proximity/diag(proximity)
  if (plot) {
    loc <- cmdscale(sqrt(1-proximity))
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(1,1))
    plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="cmdscale(Ensemble Mortality)")
    if (sum(x$cens == 1) > 1 & sum(x$cens == 0) > 1) {
      text(loc[,1], -loc[,2],
           pmax(1, ceiling(100*(x$mortality -
                      min(x$mortality, na.rm = TRUE))/diff(range(x$mortality, na.rm = TRUE)))),
           col = 3*x$cens+1,
           cex = 0.8)
    }
    else {
      text(loc[,1], -loc[,2],
           pmax(1, ceiling(100*(x$mortality -
                      min(x$mortality, na.rm = TRUE))/diff(range(x$mortality, na.rm = TRUE)))),
           cex = 0.8)
    }
    par(old.par)
  }
  invisible(proximity)
}
plot.proximity <- plot.proximity.rsf
