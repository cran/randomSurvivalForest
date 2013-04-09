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

plot.error.rsf <- function (x, sorted = TRUE, ...) {
  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")
    if (is.null(x$importance)) x$importance <- NA
    if (all(is.na(x$err.rate)) & all(is.na(x$importance))) return()
    if (all(is.na(x$importance))) {
      if (x$ntree > 1 & !all(is.na(x$err.rate))) {
        old.par <- par(no.readonly = TRUE)
        err <- rbind(x$err.rate)  
        par(mfrow = c(1,1))
        plot.err(err)
        par(old.par)      
      }
    }
    else {
      old.par <- par(no.readonly = TRUE)
      err <- rbind(x$err.rate)  
      imp <- t(rbind(x$importance))
      n.pred <- nrow(imp)
      if (sorted) pred.order <- order(imp[, 1]) else pred.order <- n.pred:1
      if (ncol(imp) == 1) max.pred <- 100 else max.pred <- 80/ncol(imp)
      if (n.pred > max.pred) {
        dotchart.labels <- rep("", n.pred)
        pretty.pt <- pretty(1:n.pred, n = max.pred)
        pretty.pt <- pretty.pt[pretty.pt > 0 & pretty.pt < n.pred]
        dotchart.labels[pretty.pt] <- x$predictorNames[pred.order][pretty.pt]
      }
      else {
        dotchart.labels <- x$predictorNames[pred.order]
      }
      if (x$ntree > 1 & !all(is.na(x$err.rate))) par(mfrow = c(1,2)) else par(mfrow = c(1,1))
      if (x$ntree > 1 & !all(is.na(x$err.rate))) {
        plot.err(err)
      }
      if (ncol(imp) > 1) { 
        colnames(imp) <- c("CHF", paste("cond CHF", 1:(ncol(imp)-1), "          "))
        imp.out <- imp[rev(pred.order),, drop = FALSE]
        colnames(imp.out) <- c("CHF", paste("condCHF.", 1:(ncol(imp)-1), sep=""))
        dotChart(imp[pred.order,, drop = FALSE], dotchart.labels)
      }
      if (ncol(imp) == 1) {
        dotChart(imp[pred.order, ], dotchart.labels)
        if (!is.null(x$predictorWt) & length(unique(x$predictorWt)) > 1 ) {
          if (length(unique(x$predictorWt)) == 1) x$predictorWt <- 1
          imp.out=as.data.frame(cbind(imp,imp/max(abs(imp),na.rm=T) ,x$predictorWt),
                            row.names=x$predictorNames)[rev(pred.order),]
          if (nrow(imp.out) == 1) imp.out[1 , 2] <- 1
          colnames(imp.out) <- c("Importance","Relative Imp","predictorWt")
        }
        else {
          imp.out=as.data.frame(cbind(imp,imp/max(abs(imp),na.rm=T)),
                            row.names=x$predictorNames)[rev(pred.order),]
          if (nrow(imp.out) == 1) imp.out[1 , 2] <- 1
          colnames(imp.out) <- c("Importance","Relative Imp")
        }
      }
      cat("\n")
      print(round(imp.out[1:min(n.pred, max.pred),, drop = FALSE],4), justify="right", print.gap=3)
      par(old.par)      
    } 
}
plot.err <- function(err) {
  matplot(1:ncol(err), t(err),
          xlab = "Number of Trees",
          ylab = "Error Rate",
          type = c("p", "l")[1 + 1 * (ncol(err) > 1)], pch = 16, lty = 1,
          lwd = c(3, 1.5)[1 + 1 * (nrow(err) > 10)])
  if (nrow(err) > 1) {
    legend("topright",
           legend = c("CHF", paste("cond CHF   ", 1:(nrow(err)-1), " ")),
           col = 1:nrow(err), lty = 1, xjust = 1,
           lwd = c(3, 1.5)[1 + 1 * (nrow(err) > 10)],
           cex = c(1, 0.75)[1 + 1 * (nrow(err) > 10)])
  }
}
dotChart <- function(x, labels = NULL) {
    if (!is.null(dim(x))) {
      ncol  <- ncol(x)
      x.dot <- NULL
      for (k in ncol(x):1) {x.dot <- c(x.dot, x[, k])}
      gcolor <- 1:ncol
    }
    else {
      x.dot <- x
      gcolor <- par("fg")
    }
    y.dot <- dot.chart.main(x, labels = labels, cex=0.90, pch="", lwd = 2, lcolor = "white", gcolor = gcolor)
    segments(rep(max(0, min(x.dot, na.rm = TRUE)) - 1e-6, length(y.dot)),
             y.dot, x.dot, y.dot, col=c(2,4)[1 + 1*(x.dot>0)], lwd = 4)
    if (min(x.dot, na.rm = TRUE) < 0) abline(v=0, lwd = 2, lty = 2, col = 1)
    mtext("Variable Importance", side = 1, line = 3, cex = 1)
  }
dot.chart.main <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pch = 21, gpch = 21, bg = par("bg"), color = par("fg"), gcolor = par("fg"), 
    lcolor = "gray", xlim = range(x[is.finite(x)]), main = NULL, 
    xlab = NULL, ylab = NULL, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    if (!is.numeric(x)) 
        stop("'x' must be a numeric vector or matrix")
    n <- length(x)
    if (is.matrix(x)) {
        if (is.null(labels)) 
            labels <- rownames(x)
        if (is.null(labels)) 
            labels <- as.character(1:nrow(x))
        labels <- rep(labels, length.out = n)
        if (is.null(groups)) 
            groups <- col(x, as.factor = TRUE)
        glabels <- levels(groups)
    }
    else {
        if (is.null(labels)) 
            labels <- names(x)
        glabels <- if (!is.null(groups)) 
            levels(groups)
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2] <- nmai[4] + max(linch + goffset, ginch) + 0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- 1:n
        y <- o
        ylim <- c(0, n + 1)
    }
    else {
        o <- sort.list(as.numeric(groups), decreasing = TRUE)
        x <- x[o]
        groups <- groups[o]
        color <- rep(color, length.out = length(groups))[o]
        lcolor <- rep(lcolor, length.out = length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
            col = color, las = 2, cex = cex, ...)
    }
    abline(h = y, lty = "dotted", col = lcolor)
    points(x, y, pch = pch, col = color, bg = bg)
    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                ...)
        }
    }
    axis(1)
    box()
    title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible(y)
  }
plot.error <- plot.error.rsf
