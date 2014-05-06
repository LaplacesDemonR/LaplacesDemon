###########################################################################
# MISS                                                                    #
#                                                                         #
# The MISS function performs multiple imputation via sequential sampling. #
###########################################################################

MISS <- function(X, Iterations=100, Algorithm="GS", Fit=NULL, verbose=TRUE)
     {
     ### Initial Checks
     if(missing(X)) stop("X is a required argument.")
     if(!is.matrix(X)) stop("X is not a matrix.")
     N <- nrow(X)
     J <- ncol(X)
     if(J < 2) stop("At least 2 columns in X are required.")
     Nmiss <- apply(X, 2, function(x) sum(is.na(x)))
     if(sum(Nmiss) == 0) stop("There are no missing values in X.")
     cat("\nNumber of Missing Values by Variable:\n")
     print(Nmiss)
     cat("\n")
     if(Algorithm != "GS") stop("Algorithm unknown.")
     for (i in 1:N)
          if(sum(is.na(X[i,])) == J) stop("All missing row found.")
     ### Parameters and Variable Type
     if(is.null(Fit)) {
          Type <- rep(1, J)
          parm <- list()
          for (j in 1:J) {
               uniq <- unique(X[complete.cases(X[,j]),j])
               if(length(uniq) == 1) stop("Constant found.")
               else if(length(uniq) == 2) {
                    if(all(c(0,1) == uniq[order(uniq)])) {
                         ### Binary Logit or Robit
                         Type[j] <- 2
                         parm[[j]] <- list(z=rep(0, sum(!is.na(X[,j]))),
                              beta=rep(0, J),
                              lambda=rep(1, sum(!is.na(X[,j]))))
                         }
                    else parm[[j]] <- list(beta=rep(0,J), sigma=0)
                    }
               else {
                    ### Linear Regression
                    if(Algorithm == "ESS")
                         parm[[j]] <- list(beta=rep(0,J), gamma=0, sigma=0)
                    else parm[[j]] <- list(beta=rep(0,J), sigma=0)
                    }
               }
          rm(uniq)
          }
     else {
          if(!identical(class(Fit), "miss"))
               stop("Fit is not an object of class miss.")
          Algorithm <- Fit$Algorithm
          parm <- Fit$parm
          Type <- Fit$Type}
     ### Observed indicator matrix O
     varnames <- colnames(X)
     O <- matrix(TRUE, N, J)
     O[which(is.na(X))] <- FALSE
     ### Initial Missing Values
     if(is.null(Fit)) {
          for (j in 1:J)
               if(Type[j] == 1)
                    X[which(is.na(X[,j])),j] <- mean(X[,j], na.rm=TRUE)
               else X[which(is.na(X[,j])),j] <- 1
          }
     else {
          if(sum(Nmiss) != length(Fit$Imp[,1]))
               stop("Length of Initial.Missings differs from number missing.")
          X[which(is.na(X))] <- Fit$Imp[,ncol(Fit$Imp)]}
     Imp <- matrix(0, sum(Nmiss), Iterations)
     ### Multiple Imputation Samplers
     GibbsLinReg <- function(y, obs, X) {
          X <- cbind(1, as.matrix(X))
          Xobs <- X[obs,]
          yobs <- y[obs]
          XtX <- t(Xobs) %*% Xobs
          ridge <- 1e-5
          penalty <- ridge * diag(XtX)
          if(length(penalty) == 1) penalty <- matrix(penalty)
          v <- as.inverse(as.symmetric.matrix(XtX + diag(penalty)))
          coef <- t(yobs %*% Xobs %*% v)
          resid <- yobs - Xobs %*% coef
          df <- max(sum(obs) - ncol(Xobs), 1)
          sigma <- sqrt(sum((resid)^2) / rchisq(1, df))
          beta <- coef + {t(chol({v + t(v)} / 2)) %*%
               rnorm(ncol(Xobs))} * sigma
          imp <- X[!obs,] %*% beta + rnorm(sum(!obs)) * sigma
          out <- list(imp=imp, beta=beta, sigma=sigma)
          return(out)}
     GibbsRobit <- function(y, obs, X, z, beta, lambda) {
          X <- cbind(1, as.matrix(X))
          Xobs <- X[obs,]
          yobs <- y[obs]
          n <- length(yobs)
          nu <- 8
          mu <- Xobs %*% beta
          z[yobs==0] <- qnorm(runif(n, 0, pnorm(0, mu, sqrt(1/lambda))),
               mu, sqrt(1/lambda))[yobs==0]
          z[yobs==1] <- qnorm(runif(n, pnorm(0, mu, sqrt(1/lambda)),1),
               mu, sqrt(1/lambda))[yobs==1]
          W <- diag(lambda)
          vbeta <- as.inverse(as.symmetric.matrix(t(Xobs) %*% W %*% Xobs))
          betahat <- vbeta %*% {t(Xobs) %*% W %*% z}
          beta <- c(rmvn(1, t(betahat), vbeta))
          lambda <- rgamma(n, {nu + 1} / 2,
               scale=2/(nu + {z - Xobs %*% beta}^2))
          mu <- X[!obs,] %*% beta
          eta <- invlogit(mu)
          imp <- rbern(length(eta), eta)
          out <- list(imp=imp, z=z, beta=beta, lambda=lambda)
          }
     ### Main Loop
     cat("\nImputation begins...\n")
     for (i in 1:Iterations) {
          cat("\nIteration:", i, "  ")
          imp <- NULL
          for (j in 1:J) {
               if(Nmiss[j] > 0) {
                    if(verbose == TRUE) cat(" V", j, " ", sep="")
                    if(Algorithm == "GS" & Type[j] == 1) {
                         out <- GibbsLinReg(y=X[,j], obs=O[,j], X=X[,-j])
                         parm[[j]]$beta <- out$beta
                         parm[[j]]$sigma <- out$sigma
                         }
                    else if(Algorithm == "GS" & Type[j] == 2) {
                         out <- GibbsRobit(y=X[,j], obs=O[,j],
                              X=X[,-j],
                              z=parm[[j]]$z, beta=parm[[j]]$beta,
                              lambda=parm[[j]]$lambda)
                         parm[[j]]$z <- out$z
                         parm[[j]]$beta <- out$beta
                         parm[[j]]$sigma <- out$sigma
                         }
                    X[,j][!O[,j]] <- out$imp
                    imp <- c(imp, out$imp)}
               if(j == J) Imp[,i] <- imp}}
     cat("\n\nEstimating Posterior Modes...")
     PostMode <- apply(Imp, 1, Mode)
     cat("\nFinished.\n")
     out <- list(Algorithm=Algorithm, Imp=Imp, parm=parm,
          PostMode=PostMode, Type=Type)
     class(out) <- "miss"
     return(out)
     }

#End
