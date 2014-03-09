###########################################################################
# joint.density.plot                                                      #
#                                                                         #
# The purpose of the joint.density.plot function is to produce a joint    #
# density plot from samples of two marginal posterior distributions. This #
# function is derived from the kde2d and bandwidth.nrd functions from the #
# MASS package.                                                           #
###########################################################################

joint.density.plot <- function(x, y, Title=NULL, contour=TRUE, color=FALSE,
     Trace=NULL)
     {
     ### Initial Checks
     xname <- deparse(substitute(x))
     yname <- deparse(substitute(y))
     x <- as.vector(x)
     y <- as.vector(y)
     if(!identical(length(y), length(x)))
          stop("vectors x and y must be the same length.")
     if(any(!is.finite(x)))
          stop("x must have finite values.")
     if(any(!is.finite(y)))
          stop("y must have finite values.")
     ### Two-Dimensional Kernel Density Estimates
     kde2d <- function(x, y, h, n=25, lims=c(range(x), range(y)))
          {
          nx <- length(x)
          if(any(!is.finite(lims)))
               stop("x and y must have finite values.")
          n <- rep(n, length.out=2L)
          gx <- seq.int(lims[1L], lims[2L], length.out=n[1L])
          gy <- seq.int(lims[3L], lims[4L], length.out=n[2L])
          h <- if(missing(h)) c(bandwidth.nrd(x), bandwidth.nrd(y))
          else rep(h, length.out=2L)
          h <- h / 4 # for S's bandwidth scale
          ax <- outer(gx, x, "-" ) / h[1L]
          ay <- outer(gy, y, "-" ) / h[2L]
          z <- tcrossprod(matrix(dnorm(ax), , nx),
               matrix(dnorm(ay), , nx)) / (nx * h[1L] * h[2L])
          list(x=gx, y=gy, z=z)
          }
     bandwidth.nrd <- function(x)
          {
          r <- quantile(x, c(0.25, 0.75))
          h <- (r[2L] - r[1L]) / 1.34
          4 * 1.06 * min(sqrt(var(x)), h) * length(x) ^ (-1/5)
          }
     dd <- kde2d(x,y)
     if(color == FALSE) {
     plot(x, y, cex=0.1, main=Title, xlab=xname, ylab=yname, col="gray")}
     else if(color == TRUE) {
          crp <- colorRampPalette(c("black","red","yellow","white"), space="rgb")
          image(dd, main=Title, xlab=xname, ylab=yname, col=crp(200))
          }
     if(contour == TRUE) {contour(dd, nlevels=10, add=TRUE)}
     if(!is.null(Trace)) {
          if(length(Trace) != 2) stop("Trace requires 2 elements.")
          if(Trace[1] >= Trace[2])
             stop("Trace[1] not smaller than Trace[2].")
          if(Trace[1] < 1) stop("Trace[1] < 1.")
          if(Trace[2] > length(x)) stop("Trace[2] > length(x).")
          lines(x[Trace[1]:Trace[2]], y[Trace[1]:Trace[2]], col="green")
          points(x[Trace[1]], y[Trace[1]], cex=0.5, col="green")
          }
     }

#End
