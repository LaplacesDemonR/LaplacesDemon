###########################################################################
# Gelman.Diagnostic                                                       #
#                                                                         #
# The purpose of the Gelman.Diagnostic function is to perform the         #
# Gelman-Rubin MCMC diagnostic on a set of posterior samples. This        #
# function is similar to the gelman.diag function in the coda package,    #
# but has been modified to work with objects of class demonoid.           #
###########################################################################

Gelman.Diagnostic <- function(x, confidence=0.95, transform=FALSE)
     {
     Nchain <- length(x)
     if(Nchain < 2) stop("More than one chain is required.")
     if(!all(sapply(x, class) == "demonoid"))
        stop("At least one item in list x is not of class demonoid.")
     Burn <- Ntot <- Niter <- Nvar <- rep(0, Nchain)
     for (i in 1:Nchain) {
          Burn[i] <- x[[i]]$Rec.BurnIn.Thinned
          Ntot[i] <- nrow(x[[i]]$Posterior1)
          if(Burn[i] >= Ntot[i]) Niter[i] <- Ntot[i]
          if(Burn[i] < Ntot[i]) Niter[i] <- Ntot[i] - Burn[i]
          Nvar[i] <- x[[i]]$Parameters}
     if(length(unique(Ntot)) != 1)
          stop("Total number of iterations differs with demonoid objects.")
     Ntot <- Ntot[1]
     Burn <- max(Burn)
     Niter <- min(Niter)
     if(length(unique(Nvar)) != 1)
          stop("Total number of parameters differs with demonoid objects.")
     else Nvar <- Nvar[1]
     xnames <- colnames(x[[1]]$Posterior1)
     if(transform == TRUE) {
          Gelman.Transform <- function(x, Nvar, Nchain)
               {
               for (i in 1:Nchain) {
                    for (j in 1:Nvar) {
                         if(min(x[[i]][,j]) > 0) {
                              if(max(x[[i]][,j]) < 1) {
                                   x[[i]][,j] <- log(x[[i]][,j] /
                                   (1 - x[[i]][,j]))}
                              else x[[i]][,j] <- log(x[[i]][,j])}
                         }
                    }
               return(x)
               }
          }
     ## Multivariate (upper case)
     if(Burn < Ntot) {
          for (i in 1:Nchain) {
               x[[i]]$Posterior1 <- x[[i]]$Posterior1[Burn:Ntot,]}}
     else warning("Non-stationary samples were used.")
     temp <- unlist(x, recursive=FALSE)
     x <- temp[names(temp) == "Posterior1"]
     if(transform == TRUE) x <- Gelman.Transform(x, Nvar, Nchain)
     S2 <- array(sapply(x, var, simplify=TRUE), dim=c(Nvar,Nvar,Nchain))
     W <- apply(S2, c(1,2), mean)
     if(Nvar > 1){
          xbar <- matrix(sapply(x, apply, 2, mean, simplify=TRUE),
               nrow=Nvar, ncol=Nchain)
          }
     else {
          xbar <- matrix(sapply(x, mean, simplify=TRUE), nrow=Nvar,
               ncol=Nchain)}
     B <- Niter * var(t(xbar))
     if(Nvar > 1) {
          CW <- chol(W)
          emax <- eigen(backsolve(CW, t(backsolve(CW, B,
               transpose=TRUE)), transpose=TRUE),
               symmetric=TRUE, only.values=TRUE)$values[1]
          mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
          }
     else mpsrf <- NULL
     ## Univariate (lower case)
     w <- diag(W)
     b <- diag(B)
     s2 <- matrix(apply(S2, 3, diag), nrow=Nvar, ncol=Nchain)
     muhat <- rowMeans(xbar)
     var.w <- .rowVars(s2) / Nchain
     var.b <- (2 * b^2) / (Nchain - 1)
     cov.wb <- (Niter / Nchain) * diag(var(t(s2), t(xbar^2)) -
          2 * muhat * var(t(s2), t(xbar)))
     V <- (Niter - 1) * w / Niter + (1 + 1/Nchain) * b / Niter
     var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 *
          var.b + 2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb) / Niter^2
     df.V <- (2 * V^2) / var.V
     df.adj <- (df.V + 3) / (df.V + 1)
     B.df <- Nchain - 1
     W.df <- (2 * w^2) / var.w
     R2.fixed <- (Niter - 1) / Niter
     R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
     R2.estimate <- R2.fixed + R2.random
     R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * R2.random
     psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
     dimnames(psrf) <- list(xnames, c("Point Est.", "Upper C.I."))
     out <- list(PSRF=psrf, MPSRF=mpsrf)
     return(out)
     }

#End


















