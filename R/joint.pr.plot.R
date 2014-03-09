###########################################################################
# joint.pr.plot                                                           #
#                                                                         #
# The purpose of the joint.pr.plot function is to plot a joint            #
# probability region. The joint.pr.plot function uses modified forms of   #
# the ellipse and dataEllipse functions from the car package.             #
###########################################################################

joint.pr.plot <- function(x, y, quantiles=c(0.25,0.50,0.75,0.95))
     {
     ### Function from car package for ellipses
     ellipse <- function(center, shape, radius, log="", center.pch=19,
          center.cex=1.5, segments=51, draw=TRUE, add=draw, xlab="",
          ylab="", col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3,
          ...)
          {
          trans.colors <- function(col, alpha=0.5, names=NULL) {
               # this function by Michael Friendly
               nc <- length(col)
               na <- length(alpha)
               # make lengths conform, filling out to the longest
               if(nc != na) {
                    col <- rep(col, length.out=max(nc,na))
                    alpha <- rep(alpha, length.out=max(nc,na))}
               clr <-rbind(col2rgb(col)/255, alpha=alpha)
               col <- rgb(clr[1,], clr[2,], clr[3,], clr[4,], names=names)
               return(col)}
          logged <- function(axis=c("x", "y")){
               axis <- match.arg(axis)
               0 != length(grep(axis, log))}
          if(!(is.vector(center) && 2==length(center)))
               stop("center must be a vector of length 2")
          if(!(is.matrix(shape) && all(2==dim(shape))))
               stop("shape must be a 2 by 2 matrix")
          if(max(abs(shape - t(shape)))/max(abs(shape)) > 1e-10)
               stop("shape must be a symmetric matrix")
          angles <- (0:segments)*2*pi/segments 
          unit.circle <- cbind(cos(angles), sin(angles)) 
          Q <- chol(shape, pivot=TRUE)
          order <- order(attr(Q, "pivot"))
          ellipse <- t(center + radius*t( unit.circle %*% Q[,order]))
          colnames(ellipse) <- c("x", "y")
          if(logged("x")) ellipse[, "x"] <- exp(ellipse[, "x"])
          if(logged("y")) ellipse[, "y"] <- exp(ellipse[, "y"])
          fill.col <- trans.colors(col, fill.alpha)
          if(draw) {
               if(add) {
                    lines(ellipse, col=col, lwd=lwd, ...) 
                    if(fill) polygon(ellipse, col=fill.col, border=NA)
               }
               else {
                    plot(ellipse, type="n", xlab=xlab, ylab=ylab, ...) 
                    lines(ellipse, col=col, lwd=lwd, ... )
                    if(fill) polygon(ellipse, col=fill.col, border=NA)
               }
               if((center.pch != FALSE) && (!is.null(center.pch)))
                    points(center[1], center[2], pch=center.pch,
                         cex=center.cex, col=col)}
          invisible(ellipse)
          }
     ### Function from car package for ellipses
     dataEllipse <- function(x, y, weights, log="", quantiles=c(0.5, 0.95),
          center.pch=19, center.cex=1.5, draw=TRUE, plot.points=draw,
          add=!plot.points, segments=51, xlab=deparse(substitute(x)),
          ylab=deparse(substitute(y)), col=palette()[1:2], lwd=1,
          fill=FALSE, fill.alpha=0.3, ...)
          {
          if(length(col) == 1) col <- rep(col, 2)
          if(missing(y)) {
               if(is.matrix(x) && ncol(x) == 2) {
                    if(missing(xlab)) xlab <- colnames(x)[1]
                    if(missing(ylab)) ylab <- colnames(x)[2]
                    y <- x[,2]
                    x <- x[,1]
                    }
              else stop("x and y must be vectors, or x must be a 2 column matrix.")
              }
          else if(!(is.vector(x) && is.vector(y) && length(x) == length(y)))
               stop("x and y must be vectors of the same length")
          if(missing(weights)) weights <- rep(1, length(x))
          if(length(weights) != length(x))
               stop("weights must be of the same length as x and y")
          if(draw) {
               if(!add) {
                    mycol <- rgb(0, 100, 0, 50, maxColorValue=255)
                    plot(x, y, type="n", col=mycol, pch=16, xlab=xlab,
                         ylab=ylab, ...)
               if(plot.points) points(x, y, col=mycol, pch=16, ...)}}
          dfn <- 2
          dfd <- length(x) - 1
          v <- cov.wt(cbind(x, y), wt=weights)
          shape <- v$cov
          center <- v$center
          result <- vector("list", length=length(quantiles))
          names(result) <- quantiles
          for (i in seq(along=quantiles)) {
               level <- quantiles[i]
               radius <- sqrt(dfn * qf(level, dfn, dfd))
               result[[i]] <- ellipse(center, shape, radius, log=log,
                    center.pch=center.pch, center.cex=center.cex,
                    segments=segments, col=col[1], lwd=lwd, fill=fill,
                    fill.alpha=fill.alpha, draw=draw, ...)}
          invisible(if(length(quantiles) == 1) result[[1]] else result)
          }
     dataEllipse(x, y, quantiles=quantiles)
     }

#End
