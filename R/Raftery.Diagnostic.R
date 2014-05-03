###########################################################################
# Raftery.Diagnostic                                                      #
#                                                                         #
# The purpose of the Raftery.Diagnostic function is to perform MCMC       #
# diagnostics on an object of class demonoid.                             #
###########################################################################

Raftery.Diagnostic <- function(x, q=0.025, r=0.005, s=0.95, eps=0.001)
     {
     if(missing(x)) stop("x is a required argument")
     if(!identical(class(x), "demonoid"))
          stop("x must be an object of class demonoid.")
     if(all(is.na(x$Posterior2))) post <- x$Posterior1
     else post <- x$Posterior2
     Thinning <- x$Thinning
     resmatrix <- matrix(nrow=ncol(post), ncol=4,
          dimnames=list(colnames(post),
          c("M", "N", "Nmin", "I")))
     phi <- qnorm(0.5 * (1 + s))
     nmin <- as.integer(ceiling((q * (1 - q) * phi^2) / r^2))
     if(nmin > nrow(post)) resmatrix <- c("Error", nmin)
     else for (i in 1:ncol(post)) {
          quant <- quantile(post[, i, drop=TRUE], probs=q)
          dichot <- post[, i, drop=TRUE] <= quant
          kthin <- 0
          bic <- 1
          while (bic >= 0) {
               kthin <- kthin + Thinning
               testres <- as.vector(Thin(dichot, By=kthin))
               testres <- factor(testres, levels=c(FALSE,TRUE))
               newdim <- length(testres)
               testtran <- table(testres[1:(newdim - 2)],
                    testres[2:(newdim - 1)], testres[3:newdim])
               testtran <- array(as.double(testtran), dim=dim(testtran))
               g2 <- 0
               for (i1 in 1:2) {
                     for (i2 in 1:2) {
                          for (i3 in 1:2) {
                               if(testtran[i1, i2, i3] != 0) {
                                    fitted <- (sum(testtran[i1, i2, 1:2]) * 
                                         sum(testtran[1:2, i2, i3])) /
                                         (sum(testtran[1:2, i2, 1:2]))
                                    g2 <- g2 + testtran[i1, i2, i3] *
                                         log(testtran[i1, i2, i3]/fitted) * 2
                                    }}}}
               bic <- g2 - log(newdim - 2) * 2
               }
          finaltran <- table(testres[1:(newdim - 1)], testres[2:newdim])
          alpha <- finaltran[1, 2]/(finaltran[1, 1] + finaltran[1, 2])
          beta <- finaltran[2, 1]/(finaltran[2, 1] + finaltran[2, 2])
          tempburn <- log((eps * (alpha + beta))/max(alpha, 
               beta))/(log(abs(1 - alpha - beta)))
          nburn <- as.integer(ceiling(tempburn) * kthin)
          tempprec <- ((2 - alpha - beta) * alpha * beta * phi^2) /
               (((alpha + beta)^3) * r^2)
          nkeep <- as.integer(ceiling(tempprec) * kthin)
          iratio <- (nburn + nkeep) / nmin
          resmatrix[i, 1] <- nburn
          resmatrix[i, 2] <- nkeep + nburn
          resmatrix[i, 3] <- nmin
          resmatrix[i, 4] <- signif(iratio, digits=3)
          }
     y <- list(params=c(q=q, r=r, s=s), resmatrix=resmatrix)
     class(y) <- "raftery"
     return(y)
     }

#End
