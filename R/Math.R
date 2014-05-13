###########################################################################
# Math                                                                    #
#                                                                         #
# This is a collection of functions to facilitate math.                   #
###########################################################################

GaussHermiteQuadRule <- function(N)
     {
     N <- abs(round(N))
     i <- seq(1, N-1, by=1)
     d <- sqrt(i/2)
     Diag <- function(x, k=0) {
          if(!is.numeric(x) && !is.complex(x))
               stop("Argument 'x' must be a real or complex vector or matrix.")
          if(!is.numeric(k) || k != round(k))
               stop("Argument 'k' must be an integer.")
          if(is.matrix(x)) {
               n <- nrow(x); m <- ncol(x)
               if(k >= m || -k >= n) {
                    y <- matrix(0, nrow=0, ncol=0)
               } else {
                    y <- x[col(x) == row(x) + k]
                    }
          } else {
               if(is.vector(x)) {
                    n <- length(x)
                    m <- n + abs(k)
                    y <- matrix(0, nrow=m, ncol=m)
                    y[col(y) == row(y) + k] <- x
               } else {
                    stop("Argument 'x' must be a real or complex vector or matrix.")
                    }
               }
          return(y)
          }
     E <- eigen(Diag(d, 1) + Diag(d, -1), symmetric=TRUE)
     L <- E$values
     V <- E$vectors
     inds <- order(L)
     x <- L[inds]
     V <- t(V[, inds])
     w <- sqrt(pi) * V[, 1]^2
     out <- list(nodes=x, weights=w)
     class(out) <- "gausshermitequadrule"
     return(out)
     }
Hermite <- function(x, N, prob=TRUE)
     {
     N <- abs(round(N))
     isBadLength <- (length(N) != 1) && (length(x) != length(N)) &&
          (length(x) != 1)
     if(isBadLength == TRUE)
     stop(paste("Argument 'n' must be either a vector of same length",
          "as argument 'x',\n  a single integer or 'x' must be a ",
          "single value!", sep=""))
     H <- function(x, N)
          {
          if(N <= 1) return(switch(N + 1, 1, x))
          else return(x * Recall(x, N - 1) - (N - 1) * Recall(x, N - 2))
          }
     scale <- 1
     if(prob == FALSE) {
          x <- sqrt(2) * x
          scale <- 2^(N / 2)}
     return(scale * mapply(H, x, N))
     }
logadd <- function(x, add=TRUE) 
     {
     x <- as.vector(x)
     x <- sort(x[is.finite(x)], decreasing=TRUE)
     x <- c(x[1], x[which(x != x[1])])
     if(length(x) == 1) return(x)
     n <- length(x)
     if(add == TRUE)
          z <- x[1] + log(1 + sum(exp(x[-1] - x[1])))
     else 
          z <- x[1] + sum(log(1 - exp(x[-1] - x[1])))
     return(z)
     }
partial <- function(Model, parm, Data, Interval=1e-6, Method="simple")
     {
     f <- Model(parm, Data)[["LP"]]
     n <- length(parm)
     if(Method == "simple") {
          if(n == 1)
               return({Model(parm + Interval, Data)[["LP"]] - f} / Interval)
          df <- rep(NA, n)
          for (i in 1:n) {
               dx <- parm
               dx[i] <- dx[i] + Interval
               df[i] <- {Model(dx, Data)[["LP"]] - f} / Interval}
          df[which(!is.finite(df))] <- 0
          return(df)
          }
     else if(Method == "Richardson") {
          zero.tol <- sqrt(.Machine$double.eps / 7e-7)
          d <- 0.0001
          r <- 4
          v <- 2
          a <- matrix(NA, r, n)
          h <- abs(d*parm) + Interval*{abs(parm) < zero.tol}
          for (k in 1:r) {
               if(n == 1)
                    a[k,] <- {Model(parm + h, Data)[["LP"]] -
                         Model(parm - h, Data)[["LP"]]} / (2*h)
               else for (i in 1:n) {
                    if((k != 1) && {abs(a[(k-1),i]) < 1e-20}) a[k,i] <- 0
	            else a[k,i] <- (Model(parm + h*(i == seq(n)), Data)[["LP"]] - 
	                 Model(parm - h*(i == seq(n)), Data)[["LP"]]) / (2*h[i])}
               a[k,which(!is.finite(a[k,]))] <- 0
               h <- h / v}
          for (m in 1:(r - 1))
               a <- (a[2:(r+1-m),,drop=FALSE]*(4^m)-a[1:(r-m),,drop=FALSE])/(4^m-1)
          return(c(a))
          }
     else stop("The", Method, "method is unknown.")
     }

#End
