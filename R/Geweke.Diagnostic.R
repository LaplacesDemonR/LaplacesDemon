###########################################################################
# Geweke.Diagnostic                                                       #
#                                                                         #
# The purpose of the Geweke.Diagnostic function is to estimate            #
# stationarity in samples according to Geweke's diagnostic. Although the  #
# code is slightly different, it is essentially the same as the           #
# geweke.diag function in the coda package.                               #
###########################################################################

Geweke.Diagnostic <- function(x) 
     {
     x <- as.matrix(x)
     if(nrow(x) < 100) return(rep(NA, ncol(x)))
     frac1 <- 0.1; frac2 <- 0.5
     startx <- 1; endx <- nrow(x)
     xstart <- c(startx, endx - frac2 * {endx - startx})
     xend <- c(startx + frac1 * {endx - startx}, endx)
     y.variance <- y.mean <- vector("list", 2)
     for (i in 1:2) {
          y <- x[xstart[i]:xend[i],]
          y.mean[[i]] <- colMeans(as.matrix(y))
          yy <- as.matrix(y)
          y <- as.matrix(y)
          max.freq <- 0.5; order <- 1; max.length <- 200
          if(nrow(yy) > max.length) {
               batch.size <- ceiling(nrow(yy) / max.length)
               yy <- aggregate(ts(yy, frequency=batch.size), nfreq=1, 
                    FUN=mean)}
          else {batch.size <- 1}
          yy <- as.matrix(yy)
          fmla <- switch(order + 1,
               spec ~ one,
               spec ~ f1,
               spec ~ f1 + f2)
          if(is.null(fmla)) stop("Invalid order.")
          N <- nrow(yy)
          Nfreq <- floor(N/2)
          freq <- seq(from=1/N, by=1/N, length=Nfreq)
          f1 <- sqrt(3) * {4 * freq - 1}
          f2 <- sqrt(5) * {24 * freq * freq - 12 * freq + 1}
          v0 <- numeric(ncol(yy))
          for (j in 1:ncol(yy)) {
               zz <- yy[,j]
               if(var(zz) == 0) v0[j] <- 0
               else {
                    yfft <- fft(zz)
                    spec <- Re(yfft * Conj(yfft)) / N
                    spec.data <- data.frame(one=rep(1, Nfreq), f1=f1,
                         f2=f2, spec=spec[1 + {1:Nfreq}],
                         inset=I(freq <= max.freq))
                    glm.out <- try(glm(fmla, family=Gamma(link="log"),
                         data=spec.data), silent=TRUE)
                    if(!inherits(glm.out, "try-error")) 
                         v0[j] <- predict(glm.out, type="response",
                              newdata=data.frame(spec=0, one=1,
                              f1=-sqrt(3), f2=sqrt(5)))
                    }
               }
          spec <- list(spec=v0)
          spec$spec <- spec$spec * batch.size
          y.variance[[i]] <- spec$spec / nrow(y)
          }
     z <- {y.mean[[1]] - y.mean[[2]]} /
          sqrt(y.variance[[1]] + y.variance[[2]])
     return(z)
     }

#End
