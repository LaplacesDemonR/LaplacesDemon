###########################################################################
# Heidelberger.Diagnostic                                                 #
#                                                                         #
# The purpose of the Heidelberger.Diagnostic function is to perform the   #
# Heidelberger and Welch MCMC convergence diagnostic on Markov chains.    #
###########################################################################

Heidelberger.Diagnostic <- function(x, eps=0.1, pvalue=0.05) 
     {
     if(missing(x)) stop("x is a required argument")
     if(!identical(class(x), "demonoid"))
          stop("x must be an object of class demonoid.")
     if(all(is.na(x$Posterior2))) x <- x$Posterior1
     else x <- x$Posterior2
     HW.mat0 <- matrix(0, ncol=6, nrow=ncol(x))
     dimnames(HW.mat0) <- list(colnames(x),
          c("stest", "start", "pvalue", "htest", "mean", "halfwidth"))
     HW.mat <- HW.mat0
     spectrum0 <- function(x, max.freq=0.5, order=1, max.length=200)
          {
          x <- as.matrix(x)
          if(!is.null(max.length) && nrow(x) > max.length) {
               batch.size <- ceiling(nrow(x) / max.length)
               x <- aggregate(ts(x, frequency=batch.size), nfreq=1,
                    FUN=mean)
               }
          else batch.size <- 1
          out <- do.spectrum0(x, max.freq=max.freq, order=order)
          out$spec <- out$spec * batch.size
          return(out)
          }
     do.spectrum0 <- function(x, max.freq=0.5, order=1)
          {
          fmla <- switch(order+1, spec ~ one, spec ~ f1, spec ~ f1 + f2)
          if(is.null(fmla)) stop("invalid order")
          N <- nrow(x)
          Nfreq <- floor(N/2)
          freq <- seq(from=1/N, by=1/N, length=Nfreq)
          f1 <- sqrt(3) * (4 * freq - 1)
          f2 <- sqrt(5) * (24 * freq^2 - 12 * freq + 1)
          v0 <- numeric(ncol(x))
          for (i in 1:ncol(x)) {
               y <- x[,i]
               v <- var(y, na.rm=TRUE)
               if(!is.finite(v)) v <- 0
               if(v == 0) v0[i] <- 0
               else {
                    yfft <- fft(y)
                    spec <- Re(yfft * Conj(yfft)) / N
                    spec.data <- data.frame(one=rep(1, Nfreq), f1=f1,
                         f2=f2, spec=spec[1 + (1:Nfreq)],
                         inset=I(freq <= max.freq))
                    glm.out <- try(glm(fmla, family=Gamma(link="log"),
                         data=spec.data), silent=TRUE)
                    if(!inherits(glm.out, "try-error"))
                         v0[i] <- predict(glm.out, type="response",
                              newdata=data.frame(spec=0, one=1,
                              f1=-sqrt(3), f2=sqrt(5)))
                    else v0[i] <- 0}}
          return(list(spec=v0))
          }
     pcramer <- function(q, eps=1.0e-5)
          {
          log.eps <- log(eps)
          y <- matrix(0, nrow=4, ncol=length(q))
          for (k in 0:3) {
               z <- gamma(k + 0.5) * sqrt(4*k + 1) /
                    (gamma(k+1) * pi^(3/2) * sqrt(q))
               u <- (4*k + 1)^2/(16*q)
               y[k+1,] <- ifelse(u > -log.eps, 0,
                    z * exp(-u) * besselK(x=u, nu=1/4))}
          return(colSums(y))
          }
     ### Heidelberger and Welch Diagnostic
     for (j in 1:ncol(x)) {
          start.vec <- seq(from=1, to=nrow(x)/2, by=nrow(x)/10)
          Y <- x[, j, drop=TRUE]
          n1 <- length(Y)
          ### Schruben's test for convergence, applied sequentially
          S0 <- spectrum0(Y[(n1/2):n1])$spec
          converged <- FALSE
          for (i in seq(along=start.vec)) {
               Y <- Y[start.vec[i]:length(Y)]
               n <- length(Y)
               ybar <- mean(Y)
               B <- cumsum(Y) - ybar * (1:n)
               Bsq <- (B * B) / (n * S0)
               I <- sum(Bsq) / n
               if(converged <- !is.na(I) && pcramer(I) < 1 - pvalue)
                    break}
          ### Recalculate S0 using section of chain that passed convergence test
          S0ci <- spectrum0(Y)$spec
          halfwidth <- 1.96 * sqrt(S0ci/n)
          passed.hw <- !is.na(halfwidth) & (abs(halfwidth/ybar) <= eps)
          if(!converged || is.na(I) || is.na(halfwidth)) {
               nstart <- NA
               passed.hw <- NA
               halfwidth <- NA
               ybar <- NA
               }
          else nstart <- start(Y)[1]
          HW.mat[j, ] <- c(converged, nstart, 1 - pcramer(I), 
               passed.hw, ybar, halfwidth)}
     class(HW.mat) <- "heidelberger"
     return(HW.mat)
     }

#End
