###########################################################################
# Matrices                                                                #
#                                                                         #
# These are utility functions for matrices.                               #
###########################################################################

as.indicator.matrix <- function(x)
     {
     n <- length(x)
     x <- as.factor(x)
     X <- matrix(0, n, length(levels(x)))
     X[(1:n) + n*(unclass(x)-1)] <- 1
     dimnames(X) <- list(names(x), levels(x))
     return(X)
     }
as.inverse <- function(x)
     {
     if(!is.matrix(x)) x <- matrix(x)
     if(!is.square.matrix(x)) stop("x must be a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x must be a symmetric matrix.")
     tol <- .Machine$double.eps
     options(show.error.messages=FALSE)
     xinv <- try(solve(x))
     if(inherits(xinv, "try-error")) {
          k <- nrow(x)
          eigs <- eigen(x, symmetric=TRUE)
          if(min(eigs$values) < tol) {
               tolmat <- diag(k)
               for (i in 1:k)
                    if(eigs$values[i] < tol) tolmat[i,i] <- 1/tol
                    else tolmat[i,i] <- 1/eigs$values[i]
               }
          else tolmat <- diag(1/eigs$values, nrow=length(eigs$values))
          xinv <- eigs$vectors %*% tolmat %*% t(eigs$vectors)
          }
     options(show.error.messages=TRUE)
     xinv <- as.symmetric.matrix(xinv)
     return(xinv)
     }
as.parm.matrix <- function(x, k, parm, Data, a=-Inf, b=Inf, restrict=FALSE,
     chol=FALSE)
     {
     X <- matrix(0, k, k)
     if(restrict == TRUE) {
          X[upper.tri(X, diag=TRUE)] <- c(1,
               parm[grep(deparse(substitute(x)),
               Data[["parm.names"]])])}
     else {
          X[upper.tri(X, diag=TRUE)] <- parm[grep(deparse(substitute(x)),
               Data[["parm.names"]])]}
     if(chol == TRUE) {
          if(a != -Inf | b != Inf) {
               x <- as.vector(X[upper.tri(X, diag=TRUE)])
               x.num <- which(x < a)
               x[x.num] <- a
               x.num <- which(x > b)
               x[x.num] <- b
               X[upper.tri(X, diag=TRUE)] <- x
               diag(X) <- abs(diag(X))
               }
          X[lower.tri(X)] <- 0
          return(X)
          }
     X[lower.tri(X)] <- t(X)[lower.tri(X)]
     if(a != -Inf | b != Inf) {
          x <- as.vector(X[upper.tri(X, diag=TRUE)])
          x.num <- which(x < a)
          x[x.num] <- a
          x.num <- which(x > b)
          x[x.num] <- b
          X[upper.tri(X, diag=TRUE)] <- x
          X[lower.tri(X)] <- t(X)[lower.tri(X)]
          }
     if(!is.symmetric.matrix(X)) X <- as.symmetric.matrix(X)
     if(!exists("LDEnv")) LDEnv <- new.env()
     if(restrict == FALSE) {
          if(is.positive.definite(X)) {
               assign("LaplacesDemonMatrix", as.vector(X[upper.tri(X,
                    diag=TRUE)]), envir=LDEnv)}
          else {
               if(exists("LaplacesDemonMatrix", envir=LDEnv)) {
                    X[upper.tri(X,
                         diag=TRUE)] <- as.vector(get("LaplacesDemonMatrix",
                         envir=LDEnv))
                    X[lower.tri(X)] <- t(X)[lower.tri(X)]}
               else {X <- diag(k)}}
          }
     if(restrict == TRUE) {
          if(is.positive.definite(X)) {
               assign("LaplacesDemonMatrix", as.vector(X[upper.tri(X,
                    diag=TRUE)][-1]), envir=LDEnv)}
          else {
               if(exists("LaplacesDemonMatrix", envir=LDEnv)) {
                    X[upper.tri(X, diag=TRUE)] <- c(1,
                         as.vector(get("LaplacesDemonMatrix",
                         envir=LDEnv)))
                    X[lower.tri(X)] <- t(X)[lower.tri(X)]
                    if(!is.symmetric.matrix(X)) X <- as.symmetric.matrix(X)
                    }
               else {X <- diag(k)}}
          }
     return(X)
     }
as.positive.definite <- function(x)
     {
     eig.tol <- 1e-06
     conv.tol <- 1e-07
     posd.tol <- 1e-08
     iter <- 0; maxit <- 100
     n <- ncol(x)
     D_S <- x
     D_S[] <- 0
     X <- x
     converged <- FALSE
     conv <- Inf
     while (iter < maxit && !converged) {
          Y <- X
          R <- Y - D_S
          e <- eigen(R, symmetric=TRUE)
          Q <- e$vectors
          d <- e$values
          p <- d > eig.tol * d[1]
          if(!any(p))
               stop("Matrix seems negative semi-definite.")
          Q <- Q[, p, drop=FALSE]
          X <- tcrossprod(Q * rep(d[p], each=nrow(Q)), Q)
          D_S <- X - R
          conv <- norm(Y - X, "I") / norm(Y, "I")
          iter <- iter + 1
          converged <- (conv <= conv.tol)
          }
     if(!converged) {
          warning("as.positive.definite did not converge in ", iter,
               " iterations.")}
     e <- eigen(X, symmetric=TRUE)
     d <- e$values
     Eps <- posd.tol * abs(d[1])
     if(d[n] < Eps) {
          d[d < Eps] <- Eps
          Q <- e$vectors
          o.diag <- diag(X)
          X <- Q %*% (d * t(Q))
          D <- sqrt(pmax(Eps, o.diag)/diag(X))
          X[] <- D * X * rep(D, each=n)}
     X <- as.symmetric.matrix(X)
     return(X)
     }
as.positive.semidefinite <- function(x)
     {
     if(!is.matrix(x)) x <- matrix(x)
     if(!is.square.matrix(x)) stop("x must be a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x must be a symmetric matrix.")
     iter <- 0; maxit <- 100
     converged <- FALSE
     while (iter < maxit && !converged) {
          iter <- iter + 1
          out <- eigen(x=x, symmetric=TRUE)
          mGamma <- t(out$vectors)
          vLambda <- out$values
          vLambda[vLambda < 0] <- 0
          x <- t(mGamma) %*% diag(vLambda) %*% mGamma
          x <- as.symmetric.matrix(x)
          if(is.positive.semidefinite(x)) converged <- TRUE
          }
     if(converged == FALSE) {
          warning("as.positive.semidefinite did not converge in ", iter,
               " iterations.")}
     return(x)
     }
as.symmetric.matrix <- function(x, k=NULL)
     {
     if(is.vector(x)) {
          if(any(!is.finite(x))) stop("x must have finite values.")
          if(is.null(k)) k <- (-1 + sqrt(1 + 8 * length(x))) / 2
          symm <- matrix(0, k, k)
          symm[lower.tri(symm, diag=TRUE)] <- x
          symm2 <- symm
          symm2[upper.tri(symm2, diag=TRUE)] <- 0
          symm <- symm + t(symm2)
          }
     else if(is.matrix(x)) {
          if(!is.square.matrix(x)) stop("x must be a square matrix.")
          if(any(!is.finite(diag(x))))
               stop("The diagonal of x must have finite values.")
          symm <- x
          x.lower.fin <- FALSE; x.upper.fin <- FALSE
          if(all(is.finite(x[lower.tri(x, diag=TRUE)])))
               x.lower.fin <- TRUE
          if(all(is.finite(x[upper.tri(x, diag=TRUE)])))
               x.upper.fin <- TRUE
          if(x.lower.fin) symm[upper.tri(x)] <- t(x)[upper.tri(x)]
          else if(x.upper.fin) symm[lower.tri(x)] <- t(x)[lower.tri(x)]
          else {
               new.up <- x[upper.tri(x)]
               new.low <- x[lower.tri(x)]
               new.up[which(!is.finite(new.up))] <- t(x)[lower.tri(x)][which(!is.finite(new.up))]
               new.low[which(!is.finite(new.low))] <- t(x)[upper.tri(x)][which(!is.finite(new.low))]
               if(any(!is.finite(c(new.up, new.low))))
                    stop("Off-diagonals in x must have finite values.")
               else {
                    symm[upper.tri(symm)] <- new.up
                    symm[lower.tri(symm)] <- new.low
                    }
               }
     }
     else stop("x must be a vector or matrix.")
     return(symm)
     }
.colVars <- function(X)
     {
     N <- nrow(X)
     Y <- X - matrix(colMeans(X), N, ncol(X), byrow=TRUE)
     Z <- colMeans(Y*Y)*N/{N-1}
     return(Z)
     }
Cov2Cor <- function(Sigma)
     {
     if(missing(Sigma)) stop("Sigma is a required argument.")
     if(any(!is.finite(Sigma))) stop("Sigma must have finite values.")
     if(is.matrix(Sigma)) {
          if(!is.positive.definite(Sigma))
               stop("Sigma is not positive-definite.")
          x <- 1 / sqrt(diag(Sigma))
          R <- x * t(x * Sigma)}
     else if(is.vector(Sigma)) {
          k <- as.integer(sqrt(length(Sigma)))
          Sigma <- matrix(Sigma, k, k)
          x <- 1 / sqrt(diag(Sigma))
          R <- as.vector(x * t(x * Sigma))}
     return(R)
     }
CovEstim <- function(Model, parm, Data, Method="Hessian")
     {
     if(Method == "Hessian") {
          VarCov <- try(-as.inverse(Hessian(Model, parm, Data)),
               silent=TRUE)
          if(!inherits(VarCov, "try-error"))
               diag(VarCov)[which(diag(VarCov) <= 0)] <- .Machine$double.eps
          else {
               cat("\nWARNING: Failure to solve matrix inversion of ",
                    "Approx. Hessian.\n", sep="")
               cat("NOTE: Identity matrix is supplied instead.\n")
                    VarCov <- diag(length(parm))}
          }
     else if(Method == "Identity") VarCov <- diag(length(parm))
     else if(Method == "OPG") {
          if(is.null(Data[["X"]])) stop("X is required in the data.")
          y <- TRUE
          if(is.null(Data[["y"]])) {
               y <- FALSE
               if(is.null(Data[["Y"]]))
                    stop("y or Y is required in the data.")}
          if(y == TRUE) {
               if(length(Data[["y"]]) != nrow(Data[["X"]]))
                    stop("length of y differs from rows in X.")
               }
          else {
               if(nrow(Data[["Y"]]) != nrow(Data[["X"]]))
                    stop("The number of rows differs in y and X.")}
          LIV <- length(parm)
          VarCov <- matrix(0, LIV, LIV)
          for (i in 1:nrow(Data[["X"]])) {
               Data.temp <- Data
               Data.temp$X <- Data.temp$X[i,,drop=FALSE]
               if(y == TRUE) Data.temp$y <- Data.temp$y[i]
               else Data.temp$Y <- Data.temp$Y[i,]
               g <- partial(Model, parm, Data.temp)
               VarCov <- VarCov + tcrossprod(g,g)}
          VarCov <- as.inverse(as.symmetric.matrix(VarCov))
          }
     else if(Method == "Sandwich") {
          B <- as.inverse(Hessian(Model, parm, Data))
          if(is.null(Data[["X"]])) stop("X is required in the data.")
          y <- TRUE
          if(is.null(Data[["y"]])) {
               y <- FALSE
               if(is.null(Data[["Y"]]))
                    stop("y or Y is required in the data.")}
          if(y == TRUE) {
               if(length(Data[["y"]]) != nrow(Data[["X"]]))
                    stop("length of y differs from rows in X.")
               }
          else {
               if(nrow(Data[["Y"]]) != nrow(Data[["X"]]))
                    stop("The number of rows differs in y and X.")}
          LIV <- length(parm)
          M <- matrix(0, LIV, LIV)
          n <- nrow(Data[["X"]])
          for (i in 1:n) {
               Data.temp <- Data
               Data.temp$X <- Data.temp$X[i,,drop=FALSE]
               if(y == TRUE) Data.temp$y <- Data.temp$y[i]
               else Data.temp$Y <- Data.temp$Y[i,]
               g <- partial(Model, parm, Data.temp)
               M <- M + tcrossprod(g,g)}
          M <- as.symmetric.matrix(M)
          VarCov <- B %*% M %*% B #Bread, Meat, Bread
          }
     else cat("\nWARNING: CovEst Method is unrecognized.")
     return(VarCov)
     }
GaussHermiteCubeRule <- function(N, dims, rule)
     {
     if(missing(rule)) Q <- GaussHermiteQuadRule(N)
     else Q <- rule
     if(dims == 1) return(Q)
     patterns_eq <- function(N, dims)
          {
          I <- matrix(1:N)
          for (i in 2:dims) {
               nf <- dim(I)[1]
               nc <- dim(I)[2]
               I2 <- cbind(kronecker(matrix(I[1, ], 1, nc),
                    matrix(1, N, 1)), (1:N))
               for (j in 2:nf)
                    I2 <- rbind(I2, cbind(kronecker(matrix(I[j, ], 1, nc),
                         matrix(1, N, 1)), (1:N)))
               I <- I2}
          return(I)
          }
     I <- patterns_eq(N, dims)
     n <- dim(I)[1]
     X2 <- matrix(0, n, dims)
     A2 <- matrix(1, n, 1)
     for (i in 1:n)
          for (j in 1:dims) {
               X2[i, j] <- Q$nodes[I[i, j]]
               A2[i, 1] <- A2[i, 1] * Q$weights[I[i, j]]}
     Max <- Q$weights[1] * Q$weights[round((N + 1)/2)]/15
     keep <- (A2 > Max)
     n2 <- sum(keep)
     X <- matrix(0, n2, dims)
     A <- matrix(1, n2, 1)
     k <- 0
     for (i in 1:n)
          if(keep[i]) {
               k <- k + 1
               X[k, ] <- X2[i, ]
               A[k, ] <- A2[i, ]}
     out <- list(nodes=X, weights=as.vector(A))
     class(out) <- "gausshermitecuberule"
     return(out)
     }
Hessian <- function(Model, parm, Data, Interval=1e-6, Method="Richardson")
     {
     if(Method == "simple") {
          parm.len <- length(parm)
          eps <- Interval * parm
          H <- matrix(0, parm.len, parm.len)
          for (i in 1:parm.len) {
               for (j in i:parm.len) {
                    x1 <- x2 <- x3 <- x4 <- parm
                    x1[i] <- x1[i] + eps[i]
                    x1[j] <- x1[j] + eps[j]
                    x2[i] <- x2[i] + eps[i]
                    x2[j] <- x2[j] - eps[j]
                    x3[i] <- x3[i] - eps[i]
                    x3[j] <- x3[j] + eps[j]
                    x4[i] <- x4[i] - eps[i]
                    x4[j] <- x4[j] - eps[j]
                    H[i, j] <- {Model(x1, Data)[["LP"]] -
                         Model(x2, Data)[["LP"]] - Model(x3, Data)[["LP"]] +
                         Model(x4, Data)[["LP"]]} / {4 * eps[i] * eps[j]}
                    }
               }
          H[lower.tri(H)] <- t(H)[lower.tri(H)]
          return(H)
          }
     else if(Method != "Richardson") stop("Method is unknown.")
     genD <- function(Model, parm, Data, Interval)
          {
          d <- 0.0001
          r <- 4
          v <- 2
          zero.tol <- sqrt(.Machine$double.eps / 7e-7)
          f0 <- Model(parm, Data)[["LP"]]
          p <- length(parm)
          h0 <- abs(d*parm) + Interval*(abs(parm) < zero.tol)
          D <- matrix(0, length(f0), (p*(p + 3)) / 2)
          Daprox <- matrix(0, length(f0), r)
          Hdiag  <- matrix(0, length(f0), p)
          Haprox <- matrix(0, length(f0), r)
          for (i in 1:p) {
               h <- h0
               for (k in 1:r) {
                    f1 <- Model(parm + (i == (1:p))*h, Data)[["LP"]]
                    f2 <- Model(parm - (i == (1:p))*h, Data)[["LP"]]
                    Daprox[,k] <- (f1 - f2)  / (2*h[i])
                    Haprox[,k] <- (f1 - 2*f0 + f2) / h[i]^2
                    h <- h / v
                    NULL}
               for (m in 1:(r - 1))
                    for (k in 1:(r - m)) {
                         Daprox[,k] <- {Daprox[,k+1]*(4^m) -
                              Daprox[,k]} / (4^m - 1)
                         Haprox[,k] <- {Haprox[,k+1]*(4^m) -
                              Haprox[,k]} / (4^m - 1)
                         NULL}
               D[,i] <- Daprox[,1]
               Hdiag[,i] <- Haprox[,1]
               NULL}
          u <- p
          for (i in 1:p) {
               for (j in 1:i) {
                    u <- u + 1
                    if(i == j) {D[,u] <- Hdiag[,i]; NULL}
                    else {
                         h <- h0
                         for (k in 1:r) {
                              f1 <- Model(parm + (i == (1:p))*h +
                                   (j == (1:p))*h, Data)[["LP"]]
                              f2 <- Model(parm - (i == (1:p))*h -
                                   (j == (1:p))*h, Data)[["LP"]]
                              Daprox[,k] <- {f1 - 2*f0 + f2 -
                                   Hdiag[,i]*h[i]^2 - 
                                   Hdiag[,j]*h[j]^2} / (2*h[i]*h[j])
                              h <- h / v}
                         for (m in 1:(r - 1))
                              for (k in 1:(r-m)) {
                                   Daprox[,k] <- {Daprox[,k+1]*(4^m) -
                                        Daprox[,k]} / (4^m - 1); NULL}
                         D[,u] <- Daprox[,1]
                         NULL}}}
          invisible(D)
          }
     D <- genD(Model, parm, Data, Interval)
     if(1 != nrow(D)) stop("BUG! should not get here.")
     H <- diag(NA, length(parm))
     u <- length(parm)
     for (i in 1:length(parm)) {
          for (j in 1:i) {
               u <- u + 1
               H[i,j] <- D[,u]}}
     H <- H + t(H)
     diag(H) <- diag(H) / 2
     return(H)
     }
is.positive.definite <- function(x)
     {
     if(!is.matrix(x)) stop("x is not a matrix.")
     if(!is.square.matrix(x)) stop("x is not a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x is not a symmetric matrix.")
     ### Deprecated Method 1
     #pd <- TRUE
     #ed <- eigen(x, symmetric=TRUE)
     #ev <- ed$values
     #if(!all(ev >= -1e-06 * abs(ev[1]))) pd <- FALSE
     ### Deprecated Method 2
     #eval <- eigen(x, only.values=TRUE)$values
     #if(any(is.complex(eval))) eval <- rep(0, max(dim(x)))
     #tol <- max(dim(x)) * max(abs(eval)) * .Machine$double.eps
     #if(all(eval > tol)) pd <- TRUE
     #else pd <- FALSE
     ### Currently Active, Method 3
     eigs <- eigen(x, symmetric=TRUE)$values
     if(any(is.complex(eigs))) return(FALSE)
     if(all(eigs > 0)) pd <- TRUE
     else pd <- FALSE
     return(pd)
     }
is.positive.semidefinite <- function(x)
     {
     if(!is.matrix(x)) stop("x is not a matrix.")
     if(!is.square.matrix(x)) stop("x is not a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x is not a symmetric matrix.")
     eigs <- eigen(x, symmetric=TRUE)$values
     if(any(is.complex(eigs))) return(FALSE)
     if(all(eigs >= 0)) pd <- TRUE
     else pd <- FALSE
     return(pd)
     }
is.square.matrix <- function(x) {return(nrow(x) == ncol(x))}
is.symmetric.matrix <- function(x) {return(sum(x == t(x)) == (nrow(x)^2))}
Jacobian <- function(Model, parm, Data, Interval=1e-6, Method="simple")
     {
     f <- Model(parm, Data)[[1]]
     n <- length(parm)
     if(Method == "simple") {
          df <-matrix(NA, length(f), n)
          for (i in 1:n) {
               dx <- parm
               dx[i] <- dx[i] + Interval
               df[,i] <- {Model(dx, Data)[[1]] - f} / Interval}
          return(df)
          }
     else if(Method == "Richardson") {
          d <- 0.0001
          zero.tol <- sqrt(.Machine$double.eps / 7e-7)
          r <- 4
          v <- 2
          a <- array(NA, c(length(f), r, n))
          h <- abs(d*parm) + Interval*(abs(parm) < zero.tol)
          for (k in 1:r) {
               for (i in 1:n) {
                    a[,k,i] <- {Model(parm + h*(i == seq(n)), Data)[[1]] -
	                 Model(parm - h*(i == seq(n)), Data)[[1]]} / (2*h[i])}
               h <- h / v}
          for (m in 1:(r - 1))
               a <- {a[,2:(r+1-m),,drop=FALSE]*(4^m) -
                    a[,1:(r-m),,drop=FALSE]} / (4^m-1)
          return(array(a, dim(a)[c(1,3)]))
          }
     else stop("The", Method, "is unknown.")
     }
logdet <- function(x)
     {
     return(2*sum(log(diag(chol(x)))))
     }
lower.triangle <- function(x, diag=FALSE)
     {
     return(x[lower.tri(x, diag=diag)])
     }
read.matrix <- function(file, header=FALSE, sep=",", nrow=0, samples=0,
     size=0, na.rm=FALSE) {
     if(nrow <= 0) nrow <- length(count.fields(file))
     con <- file(file, open="r")
     on.exit(close(con))
     if(size <= 0) size <- nrow
     if(samples <= 0) samples <- nrow
     use <- sort(sample(nrow, samples))
     now <- strsplit(readLines(con, 1), sep)[[1]]
     ncol <- length(now)
     if(header == TRUE) {
          col.names <- now
          read <- 1
          skip <- 1
          }
     else {
          col.names <- paste("X[,", 1:ncol, "]", sep="")
          read <- 0
          skip <- 0
          }
     seek(con, 0)
     X <- matrix(0, nrow=samples, ncol=ncol)
     rownames(X) <- use
     left <- nrow
     got <- 1
     while (left > 0) {
          now <- matrix(scan(file=con, sep=sep, skip=skip, n=size*ncol,
               quiet=TRUE), ncol=ncol, byrow=TRUE)
          print(dim(now))
          begin <- read + 1
          end <- read + size
          want <- (begin:end)[begin:end %in% use] - read
          if(length(want) > 0) {
               nowdat <- now[want,]
               newgot <- got + length(want) - 1
               X[got:newgot,] <- nowdat
               got <- newgot + 1}
          read <- read + size 
          left <- left - size 
          }
     colnames(X) <- col.names
     if(na.rm == TRUE) {
          num.mis <- sum(is.na(X))
          if(num.mis > 0) {
               cat("\n", num.mis, "missing value(s) found.")
               cat("\n", sum(complete.cases(X)),
                    "row(s) found with missing values.")}}
     return(X)
     }
.rowVars <- function(X)
     {
     N <- ncol(X)
     Y <- X - matrix(rowMeans(X), nrow(X), N)
     Z <- rowMeans(Y*Y)*N/{N-1}
     return(Z)
     }
SparseGrid <- function(J, K)
     {
     ### Initial Checks
     J <- max(abs(round(J)), 1)
     K <- max(abs(round(K)), 1)
     type <- "GQN" #Gauss-Hermite
     GQN <- function(level)
          {
          switch(level, {
               n <- c(0)
               w <- c(1)
          }, {
               n <- c(1)
               w <- c(0.5)
          }, {
               n <- c(0, 1.73205080756888)
               w <- c(0.666666666666667, 0.166666666666667)
          }, {
               n <- c(0.741963784302726, 2.3344142183389800)
               w <- c(0.454124145231931, 0.0458758547680685)
          }, {
               n <- c(0, 1.35562617997427, 2.85697001387281)
               w <- c(0.533333333333333, 0.222075922005613,
                      0.0112574113277207)
          }, {
               n <- c(0.616706590192594, 1.88917587775371,
                      3.32425743355212)
               w <- c(0.408828469556029, 0.0886157460419145,
                      0.00255578440205624)
          }, {
               n <- c(0, 1.15440539473997, 2.36675941073454,
                      3.75043971772574)
               w <- c(0.4571428571428580, 0.240123178605013000,
                      0.0307571239675865, 0.000548268855972219)
          }, {
               n <- c(0.539079811351375, 1.63651904243511,
                      2.802485861287540, 4.14454718612589)
               w <- c(0.373012257679077, 0.117239907661759,
                      0.00963522012078826, 0.000112614538375368)
          }, {
               n <- c(0, 1.02325566378913, 2.07684797867783,
                      3.205429002856470, 4.51274586339978)
               w <- c(0.406349206349207, 0.244097502894939,
                      0.049916406765218, 0.00278914132123177,
                      2.23458440077466e-05)
          }, {
               n <- c(0.484935707515498, 1.46598909439116,
                      2.484325841638950, 3.58182348355193,
                      4.859462828332310)
               w <- c(0.3446423349320190, 0.13548370298026700,
                      0.0191115805007703, 0.00075807093431222,
                      4.31065263071831e-06)
          }, {
               n <- c(0, 0.928868997381064, 1.87603502015485,
                      2.86512316064364, 3.93616660712998,
                      5.18800122437487)
               w <- c(0.369408369408370000, 0.24224029987397000,
                      0.066138746071057600, 0.00672028523553727,
                      0.000195671930271223, 8.1218497902149e-07)
          }, {
               n <- c(0.444403001944139, 1.34037519715162,
                      2.259464451000800, 3.22370982877010,
                      4.271825847932280, 5.50090170446775)
               w <- c(0.3216643615128300, 0.14696704804533000,
                      0.0291166879123641, 0.00220338068753318,
                      4.83718492259061e-05, 1.49992716763716e-07)
          }, {
               n <- c(0, 0.85667949351945, 1.72541837958824,
                      2.62068997343221, 3.56344438028163,
                      4.59139844893652, 5.80016725238650)
               w <- c(0.340992340992341, 0.237871522964136,
                      0.0791689558604501, 0.0117705605059965,
                      0.000681236350442926, 1.15265965273339e-05,
                      2.7226276428059e-08)
          }, {
               n <- c(0.412590457954602, 1.24268895548546,
                      2.088344745701940, 2.96303657983867,
                      3.886924575059770, 4.89693639734556,
                      6.08740954690129)
               w <- c(0.3026346268130190, 0.15408333984251400,
                      0.0386501088242534, 0.00442891910694741,
                      0.000200339553760744, 2.66099134406763e-06,
                      4.86816125774839e-09)
          }, {
               n <- c(0, 0.799129068324548, 1.60671006902873,
                      2.43243682700976, 3.28908242439877,
                      4.19620771126902, 5.19009359130478,
                      6.36394788882984)
               w <- c(0.3182595182595180, 0.2324622936097320,
                      0.0894177953998444, 0.0173657744921376,
                      0.00156735750354996, 5.64214640518902e-05,
                      5.9754195979206e-07, 8.58964989963318e-10)
          }, {
               n <- c(0.386760604500557, 1.16382910055496,
                      1.951980345716330, 2.76024504763070,
                      3.600873624171550, 4.49295530252001,
                      5.472225705949340, 6.63087819839313)
               w <- c(0.286568521238012, 0.158338372750949,
                      0.0472847523540141, 0.00726693760118474,
                      0.00052598492657391, 1.53000321624873e-05,
                      1.30947321628682e-07, 1.49781472316183e-10)
          }, {
               n <- c(0, 0.751842600703896, 1.50988330779674,
                      2.28101944025299, 3.07379717532819,
                      3.90006571719801, 4.77853158962998,
                      5.74446007865941, 6.88912243989533)
               w <- c(0.299538370126608, 0.226706308468979,
                      0.0974063711627181, 0.0230866570257112,
                      0.00285894606228465, 0.000168491431551339,
                      4.01267944797987e-06, 2.80801611793058e-08,
                      2.58431491937492e-11)
          }, {
               n <- c(0.365245755507698, 1.0983955180915,
                      1.83977992150865, 2.59583368891124,
                      3.37473653577809, 4.1880202316294,
                      5.05407268544274, 6.0077459113596,
                      7.13946484914648)
               w <- c(0.2727832346542880, 0.160685303893513,
                      0.0548966324802227, 0.0105165177519414,
                      0.00106548479629165, 5.1798961441162e-05,
                      1.02155239763698e-06, 5.90548847883655e-09,
                      4.41658876935871e-12)
          }, {
               n <- c(0, 0.71208504404238, 1.42887667607837,
                      2.15550276131694, 2.89805127651575,
                      3.66441654745064, 4.46587262683103,
                      5.32053637733604, 6.26289115651325,
                      7.38257902403043)
               w <- c(0.28377319275152100, 0.220941712199144000,
                      0.10360365727614400, 0.028666691030118500,
                      0.00450723542034204, 0.000378502109414268,
                      1.53511459546667e-05, 2.53222003209287e-07,
                      1.22037084844748e-09, 7.48283005405723e-13)
          }, {
               n <- c(0.346964157081356, 1.04294534880275,
                      1.74524732081413, 2.45866361117237,
                      3.18901481655339, 3.94396735065732,
                      4.73458133404606, 5.5787388058932,
                      6.51059015701366, 7.61904854167976)
               w <- c(0.260793063449555, 0.161739333984,
                      0.061506372063976, 0.013997837447101,
                      0.00183010313108049, 0.000128826279961929,
                      4.40212109023086e-06, 6.12749025998296e-08,
                      2.48206236231518e-10, 1.25780067243793e-13)
          }, {
               n <- c(0, 0.678045692440644, 1.35976582321123,
                      2.04910246825716, 2.75059298105237,
                      3.46984669047538, 4.21434398168842,
                      4.99496394478203, 5.82938200730447,
                      6.75144471871746, 7.84938289511382)
               w <- c(0.270260183572877, 0.21533371569506,
                      0.108392285626419, 0.0339527297865428,
                      0.00643969705140878, 0.000708047795481537,
                      4.21923474255159e-05, 1.22535483614825e-06,
                      1.45066128449307e-08, 4.97536860412175e-11,
                      2.09899121956567e-14)
          }, {
               n <- c(0.331179315715274, 0.995162422271216,
                      1.66412483911791, 2.34175999628771,
                      3.03240422783168, 3.74149635026652,
                      4.47636197731087, 5.24772443371443,
                      6.07307495112290, 6.98598042401882,
                      8.07402998402171)
               w <- c(0.250243596586935, 0.161906293413675,
                      0.0671963114288899, 0.0175690728808058,
                      0.00280876104757721, 0.000262283303255964,
                      1.33459771268087e-05, 3.319853749814e-07,
                      3.36651415945821e-09, 9.84137898234601e-12,
                      3.47946064787714e-15)
          }, {
               n <- c(0, 0.648471153534496, 1.29987646830398,
                      1.95732755293342, 2.62432363405918,
                      3.30504002175297, 4.00477532173330,
                      4.73072419745147, 5.49347398647179,
                      6.31034985444840, 7.21465943505186,
                      8.29338602741735)
               w <- c(0.258509740808839, 0.209959669577543,
                      0.112073382602621, 0.0388671837034809,
                      0.00857967839146566, 0.00116762863749786,
                      9.3408186090313e-05, 4.08997724499215e-06,
                      8.77506248386172e-08, 7.67088886239991e-10,
                      1.92293531156779e-12, 5.73238316780209e-16)
          }, {
               n <- c(0.317370096629452, 0.953421922932109,
                      1.593480429816420, 2.240467851691750,
                      2.89772864322331, 3.56930676407356,
                      4.26038360501991, 4.97804137463912,
                      5.73274717525120, 6.54167500509863,
                      7.43789066602166, 8.50780351919526)
               w <- c(0.240870115546641, 0.161459512867,
                      0.0720693640171784, 0.021126344408967,
                      0.00397660892918131, 0.000464718718779398,
                      3.2095005652746e-05, 1.21765974544258e-06,
                      2.26746167348047e-08, 1.71866492796487e-10,
                      3.71497415276242e-13, 9.39019368904192e-17)
          }, {
               n <- c(0, 0.622462279186076, 1.24731197561679,
                      1.87705836994784, 2.51447330395221,
                      3.16277567938819, 3.82590056997249,
                      4.50892992296729, 5.21884809364428,
                      5.96601469060670, 6.76746496380972,
                      7.65603795539308, 8.71759767839959)
               w <- c(0.248169351176485, 0.20485102565034,
                      0.114880924303952, 0.043379970167645,
                      0.0108567559914623, 0.0017578504052638,
                      0.000177766906926527, 1.06721949052025e-05,
                      3.5301525602455e-07, 5.73802386889938e-09,
                      3.79115000047719e-11, 7.10210303700393e-14,
                      1.53003899799868e-17)
          })
          return(list(nodes=n, weights=w))
          }
     SparseGridGetSeq <- function(J, norm)
          {
          seq.vec <- rep(0, J)
          a <- norm - J
          seq.vec[1] <- a
          fs <- matrix(seq.vec, nrow=1, ncol=length(seq.vec))
          cnt <- 1
          while (seq.vec[J] < a) {
               if(cnt == J) {
                    for (i in seq(cnt - 1, 1, -1)) {
                         cnt <- i
                         if(seq.vec[i] != 0) {
                              break
                              }
                         }
                    }
               seq.vec[cnt] <- seq.vec[cnt] - 1
               cnt <- cnt + 1
               seq.vec[cnt] <- a - sum(seq.vec[1:(cnt - 1)])
               if(cnt < J) {
                    seq.vec[(cnt + 1):J] <- rep(0, J -
                         cnt)
                    }
               fs <- rbind(fs, seq.vec)
               }
          fs <- fs + 1
          return(fs)
          }
     SparseGridKronProd <- function(n1D, w1D)
          {
          nodes <- matrix(n1D[[1]], nrow=length(n1D[[1]]), ncol=1)
          weights <- w1D[[1]]
          if(length(n1D) > 1) {
               for (j in 2:length(n1D)) {
                    newnodes <- n1D[[j]]
                    nodes <- cbind(kronecker(nodes, rep(1, length(newnodes))),
                         kronecker(rep(1, nrow(nodes)), newnodes))
                    weights <- kronecker(weights, w1D[[j]])
                    }
               }
          return(list(nodes=nodes, weights=weights))
          }
     sortrows <- function (A, index.return=FALSE)
          {
          if(!is.matrix(A)) {
               stop("SparseGrid:::sortrows expects a matrix as argument A.")
               }
          A.nrow <- nrow(A)
          A.ncol <- ncol(A)
          if(index.return == TRUE) indices <- 1:nrow(A)
          for (col.cnt in seq(ncol(A), 1, -1)) {
               tmp.indices <- order(A[, col.cnt])
               A <- A[tmp.indices, , drop=FALSE]
               if(index.return == TRUE) indices <- indices[tmp.indices]
               }
          if(index.return == TRUE) {
               res <- list(x=matrix(A, nrow=A.nrow, ncol=A.ncol),
                    ix=indices)
               }
          else res <- matrix(A, nrow=A.nrow, ncol=A.ncol)
          return(res)
     }
     tryCatch({
          n1D <- vector(mode="list", length=K)
          w1D <- vector(mode="list", length=K)
          R1D <- rep(0, K)
          for (level in 1:K) {
               res <- eval(call(type, level))
               nodes <- res$nodes
               weights <- res$weights
               R1D[level] <- length(weights)
               n1D[[level]] <- nodes
               w1D[[level]] <- weights
               }
          }, error=function(e) {cat("Error evaluating the 1D rule\n")})
     minq <- max(0, K - J)
     maxq <- K - 1
     nodes <- matrix(0, nrow=0, ncol=J)
     weights <- numeric(0)
     for (q.cnt in minq:maxq) {
          r <- length(weights)
          bq <- (-1)^(maxq - q.cnt) * choose(J - 1, J +
               q.cnt - K)
          indices.mat <- SparseGridGetSeq(J, J + q.cnt)
          Rq <- sapply(1:nrow(indices.mat), function(row.cnt) {
               prod(R1D[indices.mat[row.cnt, ]])})
          Rq.sum <- sum(Rq)
          nodes <- rbind(nodes, matrix(0, nrow=Rq.sum, ncol=J))
          weights <- c(weights, rep(0, Rq.sum))
          for (j in 1:nrow(indices.mat)) {
               midx <- indices.mat[j, ]
               res <- SparseGridKronProd(n1D[midx], w1D[midx])
               nodes[(r + 1):(r + Rq[j]), ] <- res$nodes
               weights[(r + 1):(r + Rq[j])] <- bq * res$weights
               r <- r + Rq[j]
               }
          nodes.sorted <- sortrows(nodes, index.return=TRUE)
          nodes <- nodes.sorted$x
          weights <- weights[nodes.sorted$ix]
          keep <- 1
          lastkeep <- 1
          if(nrow(nodes) > 1) {
               for (j in 2:nrow(nodes)) {
                    if(all(nodes[j, ] == nodes[j - 1, ])) {
                         weights[lastkeep] <- weights[lastkeep] + weights[j]
                         }
                    else {
                         lastkeep <- j
                         keep <- c(keep, j)
                         }
                    }
               }
          nodes <- matrix(nodes[keep, ], nrow=length(keep), ncol=J)
          weights <- weights[keep]
          }
     nr <- length(weights)
     m <- n1D[[1]]
     for (j in 1:J) {
          keep <- rep(0, nr)
          numnew <- 0
          for (r in 1:nr) {
               if(nodes[r, j] != m) {
                    numnew <- numnew + 1
                    keep[numnew] <- r
                    }
               }
          if(numnew > 0) {
               nodes <- rbind(nodes, matrix(nodes[keep[1:numnew],
                    ], nrow=numnew, ncol=J))
               nodes[(nr + 1):(nr + numnew), j] <- 2 * m - nodes[(nr +
                    1):(nr + numnew), j]
               weights <- c(weights, weights[keep[1:numnew]])
               nr <- nr + numnew
               }
          }
     nodes.sorted <- sortrows(nodes, index.return=TRUE)
     nodes <- nodes.sorted$x
     weights <- weights[nodes.sorted$ix]
     weights <- abs(weights)/sum(abs(weights))
     return(list(nodes=nodes, weights=weights))
     }
TransitionMatrix <- function(theta.y=NULL, y.theta=NULL, p.theta=NULL)
     {
     if(!is.null(theta.y)) {
          theta.y <- as.vector(theta.y)
          N <- length(theta.y)
          theta.y <- as.matrix(table(theta.y[-N], theta.y[-1]))
          theta.y <- theta.y / rowSums(theta.y)
          return(theta.y)
          }
     else {
          N <- length(y.theta)
          y.theta <- as.matrix(table(From=y.theta[-N], To=y.theta[-1]))
          y.theta <- y.theta / rowSums(y.theta)
          if(!is.null(p.theta)) {
               if(!is.matrix(p.theta)) p.theta <- as.matrix(p.theta)
               if(!identical(dim(y.theta),dim(p.theta)))
                    stop("Dimensions of p.theta differ from the matrix of y.theta.")
               theta.y <- y.theta * p.theta
               theta.y <- theta.y / rowSums(theta.y)
               }
          else p.theta <- y.theta
          return(p.theta)
          }
     }
tr <- function(x) {return(sum(diag(x)))}
upper.triangle <- function(x, diag=FALSE)
     {
     return(x[upper.tri(x, diag=diag)])
     }

#End
