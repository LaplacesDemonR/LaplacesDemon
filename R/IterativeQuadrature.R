###########################################################################
# IterativeQuadrature                                                     #
#                                                                         #
# The purpose of the IterativeQuadrature function is to perform iterative #
# quadrature on a Bayesian model.                                         #
###########################################################################

IterativeQuadrature <- function(Model, parm, Data, Covar=NULL,
     Iterations=100, Algorithm="CAGH", Specs=NULL, Samples=1000, sir=TRUE,
     Stop.Tolerance=c(1e-5,1e-15), CPUs=1, Type="PSOCK")
     {
     cat("\nIterativeQuadrature was called on ", date(), "\n", sep="")
     time1 <- proc.time()
     IQ.call <- match.call()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n")
     if(missing(Model)) stop("Model is a required argument.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data)) stop("Data is a required argument.")
     if(missing(parm)) {
          cat("Initial values were not supplied, and\n")
          cat("have been set to zero prior to LaplaceApproximation().\n")
          parm <- rep(0, length(Data[["parm.names"]]))}
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n")}}}}
     Iterations <- min(max(round(Iterations), 2), 1000000)
     "%!in%" <- function(x,table) return(match(x, table, nomatch=0) == 0)
     if(Algorithm %!in% c("AGH","AGHSG","CAGH"))
          stop("Algorithm is unknown.")
     if(Algorithm == "AGH") {
          Algorithm <- "Adaptive Gauss-Hermite"
          if(missing(Specs)) stop("The Specs argument is required.")
          if(length(Specs) != 4) stop("The Specs argument is incorrect.")
          if(Specs[["N"]] == Specs[[1]])
               N <- max(abs(round(Specs[[1]])), 2)
          if(Specs[["Nmax"]] == Specs[[2]])
               Nmax <- max(abs(round(Specs[[2]])), N)
          Packages <- NULL
          if(!is.null(Specs[["Packages"]]))
               Packages <- Specs[["Packages"]]
          Dyn.libs <- NULL
          if(!is.null(Specs[["Dyn.libs"]]))
               Dyn.libs <- Specs[["Dyn.libs"]]
          }
     else if(Algorithm == "AGHSG") {
          Algorithm <- "Adaptive Gauss-Hermite Sparse Grid"
          if(missing(Specs)) stop("The Specs argument is required.")
          if(length(Specs) != 4) stop("The Specs argument is incorrect.")
          if(Specs[["K"]] == Specs[[1]])
               K <- max(abs(round(Specs[[1]])), 2)
          if(Specs[["Kmax"]] == Specs[[2]])
               Kmax <- max(abs(round(Specs[[2]])), K)
          Packages <- NULL
          if(!is.null(Specs[["Packages"]]))
               Packages <- Specs[["Packages"]]
          Dyn.libs <- NULL
          if(!is.null(Specs[["Dyn.libs"]]))
               Dyn.libs <- Specs[["Dyn.libs"]]
          }
     else if(Algorithm == "CAGH") {
          Algorithm <- "Componentwise Adaptive Gauss-Hermite"
          if(missing(Specs)) stop("The Specs argument is required.")
          if(length(Specs) != 4) stop("The Specs argument is incorrect.")
          if(Specs[["N"]] == Specs[[1]])
               N <- max(abs(round(Specs[[1]])), 2)
          if(Specs[["Nmax"]] == Specs[[2]])
               Nmax <- max(abs(round(Specs[[2]])), N)
          Packages <- NULL
          if(!is.null(Specs[["Packages"]]))
               Packages <- Specs[["Packages"]]
          Dyn.libs <- NULL
          if(!is.null(Specs[["Dyn.libs"]]))
               Dyn.libs <- Specs[["Dyn.libs"]]
          }
     if(any(Stop.Tolerance <= 0)) Stop.Tolerance <- c(1e-5,1e-15)
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out,1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of apply functions\n")
          cat(     "were found in the Model specification. Iteration speed will\n")
          cat("     increase if apply functions are vectorized in R or coded\n")
          cat("     in a faster language such as C++ via the Rcpp package.\n")}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if for loops are vectorized in R or coded in a\n")
          cat("     faster language such as C++ via the Rcpp package.\n")}
     ###########################  Preparation  ############################
     m.old <- Model(parm, Data)
     if(!is.list(m.old)) stop("Model must return a list.")
     if(length(m.old) != 5) stop("Model must return five components.")
     if(any(names(m.old) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.")
     if(length(m.old[["LP"]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(parm), length(m.old[["parm"]])))
          stop("The number of initial values and parameters differs.")
     if(!is.finite(m.old[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n")
          if(!is.null(Data[["PGF"]]))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data, PGF=FALSE)
          m.old <- Model(Initial.Values, Data)
          }
     if(!is.finite(m.old[["LP"]])) stop("The posterior is non-finite.")
     if(!is.finite(m.old[["Dev"]])) stop("The deviance is non-finite.")
     parm <- m.old[["parm"]]
     LIV <- length(parm)
     ScaleF <- 2.381204 * 2.381204 / LIV
     if(is.null(Covar)) Covar <- diag(LIV) * ScaleF
     if(is.vector(Covar))
          if(identical(length(Covar), LIV)) Covar <- diag(LIV) * Covar
     if(!is.matrix(Covar)) stop("Covar must be a matrix.")
     if(!is.square.matrix(Covar)) stop("Covar must be a square matrix.")
     if(!is.symmetric.matrix(Covar)) Covar <- as.symmetric.matrix(Covar)
     ####################  Begin Iterative Quadrature  ####################
     cat("Algorithm:", Algorithm, "\n")
     cat("\nIterativeQuadrature is beginning to update...\n")
     if(Algorithm == "Adaptive Gauss-Hermite") {
          IQ <- .iqagh(Model, parm, Data, Covar, Iterations,
               Stop.Tolerance, CPUs, m.old, N, Nmax, Packages, Dyn.libs)
          }
     else if(Algorithm == "Adaptive Gauss-Hermite Sparse Grid") {
          IQ <- .iqaghsg(Model, parm, Data, Covar, Iterations,
               Stop.Tolerance, CPUs, m.old, K, Kmax, Packages, Dyn.libs)
          }
     else if(Algorithm == "Componentwise Adaptive Gauss-Hermite") {
          IQ <- .iqcagh(Model, parm, Data, Covar, Iterations,
               Stop.Tolerance, CPUs, m.old, N, Nmax, Packages, Dyn.libs)
          }
     VarCov <- IQ$Covar
     Dev <- as.vector(IQ$Dev)
     iter <- IQ$iter
     LPw <- IQ$LPw
     M <- IQ$M
     mu <- IQ$mu
     N <- IQ$N
     parm.len <- ncol(mu)
     parm.new <- mu[nrow(mu),]
     parm.old <- IQ$parm.old
     tol <- IQ$tol
     Z <- IQ$Z
     rm(IQ)
     if(iter == 1) stop("IterativeQuadrature stopped at iteration 1.")
     converged <- FALSE
     if({tol[1] <= Stop.Tolerance[1]} & {tol[2] <= Stop.Tolerance[2]})
          converged <- TRUE
     ### Column names to samples
     if(ncol(mu) == length(Data[["parm.names"]]))
          colnames(mu) <- Data[["parm.names"]]
     rownames(mu) <- 1:nrow(mu)
     cat("\n\n")
     #################  Sampling Importance Resampling  ##################
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Sampling from Posterior with Sampling Importance Resampling\n")
          posterior <- SIR(Model, Data, mu=parm.new, Sigma=VarCov,
               n=Samples, CPUs=CPUs, Type=Type)
          Mon <- matrix(0, nrow(posterior), length(Data[["mon.names"]]))
          dev <- rep(0, nrow(posterior))
          for (i in 1:nrow(posterior)) {
               mod <- Model(posterior[i,], Data)
               dev[i] <- mod[["Dev"]]
               Mon[i,] <- mod[["Monitor"]]
               }
          colnames(Mon) <- Data[["mon.names"]]}
     else {
          if({sir == TRUE} & {converged == FALSE})
               cat("Posterior samples are not drawn due to Converge=FALSE\n")
          posterior <- NA; Mon <- NA}
     #####################  Summary, Point-Estimate  ######################
     cat("Creating Summary from Point-Estimates\n")
     Summ1 <- matrix(NA, parm.len, 4, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","LB","UB")))
     Summ1[,1] <- parm.new
     Summ1[,2] <- sqrt(diag(VarCov))
     Summ1[,3] <- parm.new - 2*Summ1[,2]
     Summ1[,4] <- parm.new + 2*Summ1[,2]
     ###################  Summary, Posterior Samples  ####################
     Summ2 <- NA
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Creating Summary from Posterior Samples\n")
          Summ2 <- matrix(NA, ncol(posterior), 7,
               dimnames=list(Data[["parm.names"]],
                    c("Mean","SD","MCSE","ESS","LB","Median","UB")))
          Summ2[,1] <- colMeans(posterior)
          Summ2[,2] <- sqrt(.colVars(posterior))
          Summ2[,3] <- Summ2[,2] / sqrt(nrow(posterior))
          Summ2[,4] <- rep(nrow(posterior), ncol(posterior))
          Summ2[,5] <- apply(posterior, 2, quantile, c(0.025))
          Summ2[,6] <- apply(posterior, 2, quantile, c(0.500))
          Summ2[,7] <- apply(posterior, 2, quantile, c(0.975))
          Deviance <- rep(0, 7)
          Deviance[1] <- mean(dev)
          Deviance[2] <- sd(dev)
          Deviance[3] <- sd(dev) / sqrt(nrow(posterior))
          Deviance[4] <- nrow(posterior)
          Deviance[5] <- as.numeric(quantile(dev, probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(dev, probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(dev, probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:ncol(Mon)) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[,j])
               Monitor[2] <- sd(as.vector(Mon[,j]))
               Monitor[3] <- sd(as.vector(Mon[,j])) / sqrt(nrow(Mon))
               Monitor[4] <- nrow(Mon)
               Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
                    na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
                    na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
                    na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data[["mon.names"]][j]
               }
          }
     ###############  Logarithm of the Marginal Likelihood  ###############
     LML <- list(LML=NA, VarCov=VarCov)
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          lml <- LML(theta=posterior, LL=(dev*(-1/2)), method="NSIS")
          LML[[1]] <- lml[[1]]}
     else if({sir == FALSE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(Model, Data, Modes=parm.new, Covar=VarCov,
               method="LME")}
     colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]
     time2 <- proc.time()
     #############################  Output  ##############################
     cat("Creating Output\n")
     IQ <- list(Algorithm=Algorithm,
          Call=IQ.call,
          Converged=converged,
          Covar=VarCov,
          Deviance=as.vector(Dev),
          History=mu,
          Initial.Values=parm,
          Iterations=iter,
          LML=LML[[1]],
          LP.Final=as.vector(Model(parm.new, Data)[["LP"]]),
          LP.Initial=m.old[["LP"]],
          LPw=LPw,
          M=M,
          Minutes=round(as.vector(time2[3] - time1[3]) / 60, 2),
          Monitor=Mon,
          N=N,
          Posterior=posterior,
          Summary1=Summ1,
          Summary2=Summ2,
          Tolerance.Final=tol,
          Tolerance.Stop=Stop.Tolerance,
          Z=Z)
     class(IQ) <- "iterquad"
     cat("\nIterativeQuadrature is finished.\n")
     return(IQ)
     }
.iqagh <- function(Model, parm, Data, Covar, Iterations, Stop.Tolerance,
     CPUs, m.old, N, Nmax, Packages, Dyn.libs)
     {
     if(CPUs == 1) {
          Dev <- matrix(m.old[["Dev"]],1,1)
          LPworst <- m.old[["LP"]]
          parm.len <- length(parm)
          if(parm.len > 10)
               warning("The number of parameters may be too large.")
          mu <- rbind(parm)
          parm.old <- parm
          colnames(mu) <- Data[["parm.names"]]
          sigma <- sqrt(diag(Covar))
          rule <- GaussHermiteCubeRule(N=N, dims=parm.len)
          X <- rule$nodes
          W <- rule$weights
          n <- length(W)
          count <- 0
          expand <- FALSE
          IE <- 1
          tol <- c(1,1)
          ### Begin Iterative Quadrature
          options(warn=-1)
          cat("\nIteration: 1,   Univariate Nodes: ", N,
               ",   Multivariate Nodes: ", n, sep="")
          for (iter in 2:Iterations) {
               ### New Quadrature Rule
               if({expand == TRUE} & {N < Nmax}) {
                    N <- min(N+1, Nmax)
                    rule <- GaussHermiteCubeRule(N=N, dims=parm.len)
                    X <- rule$nodes
                    W <- rule$weights
                    n <- length(W)}
               adjust <- FALSE
               expand <- FALSE
               cat("\nIteration: ", iter, ",   Univariate Nodes: ", N,
                    ",   Multivariate Nodes: ", n, sep="")
               LP <- rep(LPworst, n)
               LPw <- M <- Z <- matrix(0, n, parm.len)
               mu <- rbind(mu, 0)
               ### Evaluate at the Nodes
               for (i in 1:n) {
                    Z[i,] <- mu[iter-1,] + sqrt(2)*sigma*X[i,]
                    mod <- Model(Z[i,], Data)
                    if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
                         mod[["Monitor"]]))))
                         LP[i] <- mod[["LP"]]
                    Z[i,] <- mod[["parm"]]}
               ### Weights
               LPw <- matrix(exp(LP - logadd(LP)), n, parm.len)
               LPw <- LPw / matrix(colSums(LPw), n, parm.len)
               LPw[which(!is.finite(LPw))] <- 0
               if(all(LPw == 0)) {
                    expand <- TRUE
                    LPw <- rep(.Machine$double.eps, n)}
               ### Adapt mu and Sigma
               mu[iter,] <- colSums(LPw * Z)
               temp <- cov.wt(Z, wt=LPw[,1])$cov
               if(all(is.finite(temp))) Covar <- as.symmetric.matrix(temp)
               else expand <- TRUE
               diag(Covar) <- abs(diag(Covar))
               diag(Covar)[which(diag(Covar) < .Machine$double.eps)] <- .Machine$double.eps
               M <- W * exp(X) * sqrt(2) *
                    matrix(sigma, n, parm.len, byrow=TRUE)
               sigma <- sqrt(diag(Covar))
               mod <- Model(mu[iter,], Data)
               Dev <- rbind(Dev, mod[["Dev"]])
               ### Accept Only Improved mu
               m.old <- Model(mu[iter-1,], Data)
               if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
                    mod[["Monitor"]])))) {
                    if(m.old[["LP"]] >= mod[["LP"]]) {
                         expand <- TRUE
                         mu[iter,] <- mu[iter-1,]}
                    }
               else {
                    expand <- TRUE
                    mu[iter,] <- mu[iter-1,]}
               ### Integration Error and Tolerance
               Mw <- M / matrix(colSums(M), n, parm.len, byrow=TRUE)
               IE <- c(IE, mean(abs(Mw - LPw)))
               tol[1] <- sqrt(sum({mu[iter,] - mu[iter-1,]}^2))
               tol[2] <- abs(IE[iter] - IE[iter-1])
               if(tol[1] <= Stop.Tolerance[1]) {
                    expand <- TRUE
                    if(tol[2] <= Stop.Tolerance[2]) {
                         if(count > 0) break
                         count <- count + 1}
                    }
               else count <- 0
               }
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n")
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n")
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...")
          cat("\n##################################################\n")
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n")
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
               ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())
          Dev <- matrix(m.old[["Dev"]],1,1)
          LPworst <- m.old[["LP"]]
          parm.len <- length(parm)
          if(parm.len > 10)
               warning("The number of parameters may be too large.")
          mu <- rbind(parm)
          parm.old <- parm
          colnames(mu) <- Data[["parm.names"]]
          sigma <- sqrt(diag(Covar))
          rule <- GaussHermiteCubeRule(N=N, dims=parm.len)
          X <- rule$nodes
          W <- rule$weights
          n <- length(W)
          count <- 0
          expand <- FALSE
          IE <- 1
          tol <- c(1,1)
          ### Begin Iterative Quadrature
          options(warn=-1)
          cat("\nIteration: 1,   Univariate Nodes: ", N,
               ",   Multivariate Nodes: ", n, sep="")
          for (iter in 2:Iterations) {
               ### New Quadrature Rule
               if({expand == TRUE} & {N < Nmax}) {
                    N <- min(N+1, Nmax)
                    rule <- GaussHermiteCubeRule(N=N, dims=parm.len)
                    X <- rule$nodes
                    W <- rule$weights
                    n <- length(W)}
               adjust <- FALSE
               expand <- FALSE
               cat("\nIteration: ", iter, ",   Univariate Nodes: ", N,
                    ",   Multivariate Nodes: ", n, sep="")
               LP <- rep(LPworst, n)
               LPw <- M <- Z <- matrix(0, n, parm.len)
               mu <- rbind(mu, 0)
               ### Evaluate at the Nodes
               Z <- matrix(mu[iter-1,], n, parm.len, byrow=TRUE) +
                    sqrt(2) * matrix(sigma, n, parm.len, byrow=TRUE) * X
               mod <- parLapply(cl, 1:n, function(x) Model(Z[x,], Data))
               LP <- as.vector(unlist(lapply(mod, function(x) x[["LP"]])))
               Z <- matrix(as.vector(unlist(lapply(mod,
                    function(x) x[["parm"]]))), n, parm.len, byrow=TRUE)
               LP[which(!is.finite(LP))] <- LPworst
               ### Weights
               LPw <- matrix(exp(LP - logadd(LP)), n, parm.len)
               LPw <- LPw / matrix(colSums(LPw), n, parm.len)
               LPw[which(!is.finite(LPw))] <- 0
               if(all(LPw == 0)) {
                    expand <- TRUE
                    LPw <- rep(.Machine$double.eps, n)}
               ### Adapt mu and Sigma
               mu[iter,] <- colSums(LPw * Z)
               temp <- cov.wt(Z, wt=LPw[,1])$cov
               if(all(is.finite(temp))) Covar <- as.symmetric.matrix(temp)
               else expand <- TRUE
               diag(Covar) <- abs(diag(Covar))
               diag(Covar)[which(diag(Covar) < .Machine$double.eps)] <- .Machine$double.eps
               M <- W * exp(X) * sqrt(2) *
                    matrix(sigma, n, parm.len, byrow=TRUE)
               sigma <- sqrt(diag(Covar))
               mod <- Model(mu[iter,], Data)
               Dev <- rbind(Dev, mod[["Dev"]])
               ### Accept Only Improved mu
               m.old <- Model(mu[iter-1,], Data)
               if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
                    mod[["Monitor"]])))) {
                    if(m.old[["LP"]] >= mod[["LP"]]) {
                         expand <- TRUE
                         mu[iter,] <- mu[iter-1,]}
                    }
               else {
                    expand <- TRUE
                    mu[iter,] <- mu[iter-1,]}
               ### Integration Error and Tolerance
               Mw <- M / matrix(colSums(M), n, parm.len, byrow=TRUE)
               IE <- c(IE, mean(abs(Mw - LPw)))
               tol[1] <- sqrt(sum({mu[iter,] - mu[iter-1,]}^2))
               tol[2] <- abs(IE[iter] - IE[iter-1])
               if(tol[1] <= Stop.Tolerance[1]) {
                    expand <- TRUE
                    if(tol[2] <= Stop.Tolerance[2]) {
                         if(count > 0) break
                         count <- count + 1}
                    }
               else count <- 0
               }
          }
     options(warn=0)
     ### Output
     IQ <- list(Covar=Covar, Dev=Dev, iter=iter, LPw=LPw, M=M,
          mu=mu, N=n, parm.old=parm.old, tol=tol, Z=Z)
     return(IQ)
     }
.iqaghsg <- function(Model, parm, Data, Covar, Iterations, Stop.Tolerance,
     CPUs, m.old, K, Kmax, Packages, Dyn.libs)
     {
     if(CPUs == 1) {
          Dev <- matrix(m.old[["Dev"]],1,1)
          LPworst <- m.old[["LP"]]
          parm.len <- length(parm)
          if(parm.len > 10)
               warning("The number of parameters may be too large.")
          mu <- rbind(parm)
          parm.old <- parm
          colnames(mu) <- Data[["parm.names"]]
          sigma <- sqrt(diag(Covar))
          rule <- SparseGrid(J=parm.len, K=K)
          X <- rule$nodes
          W <- rule$weights
          N <- length(W)
          count <- 0
          expand <- FALSE
          IE <- 1
          tol <- c(1,1)
          ### Begin Iterative Quadrature
          options(warn=-1)
          cat("\nIteration: 1,   K: ", K, ",   Nodes: ", N, sep="")
          for (iter in 2:Iterations) {
               ### New Quadrature Rule
               if({expand == TRUE} & {K < Kmax}) {
                    K <- min(K+1, Kmax)
                    rule <- SparseGrid(J=parm.len, K=K)
                    X <- rule$nodes
                    W <- rule$weights
                    N <- length(W)}
               expand <- FALSE
               cat("\nIteration: ", iter, ",   K: ", K, ",   Nodes: ", N,
                    sep="")
               LP <- rep(LPworst, N)
               LPw <- M <- Z <- matrix(0, N, parm.len)
               mu <- rbind(mu, 0)
               ### Evaluate at the Nodes
               for (i in 1:N) {
                    Z[i,] <- mu[iter-1,] + sqrt(2)*sigma*X[i,]
                    mod <- Model(Z[i,], Data)
                    if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
                         mod[["Monitor"]]))))
                         LP[i] <- mod[["LP"]]
                    Z[i,] <- mod[["parm"]]}
               ### Weights
               LPw <- matrix(exp(LP - logadd(LP)), N, parm.len)
               LPw <- LPw / matrix(colSums(LPw), N, parm.len)
               LPw[which(!is.finite(LPw))] <- 0
               if(all(LPw == 0)) {
                    expand <- TRUE
                    LPw <- rep(.Machine$double.eps, N)}
               ### Adapt mu and Sigma
               mu[iter,] <- colSums(LPw * Z)
               temp <- cov.wt(Z, wt=LPw[,1])$cov
               if(all(is.finite(temp))) Covar <- as.symmetric.matrix(temp)
               else expand <- TRUE
               diag(Covar) <- abs(diag(Covar))
               diag(Covar)[which(diag(Covar) < .Machine$double.eps)] <- .Machine$double.eps
               M <- W * exp(X) * sqrt(2) *
                    matrix(sigma, N, parm.len, byrow=TRUE)
               sigma <- sqrt(diag(Covar))
               mod <- Model(mu[iter,], Data)
               Dev <- rbind(Dev, mod[["Dev"]])
               ### Accept Only Improved mu
               m.old <- Model(mu[iter-1,], Data)
               if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
                    mod[["Monitor"]])))) {
                    if(m.old[["LP"]] >= mod[["LP"]]) {
                         expand <- TRUE
                         mu[iter,] <- mu[iter-1,]}
                    }
               else {
                    expand <- TRUE
                    mu[iter,] <- mu[iter-1,]}
               ### Integration Error and Tolerance
               Mw <- M / matrix(colSums(M), N, parm.len, byrow=TRUE)
               IE <- c(IE, mean(abs(Mw - LPw)))
               tol[1] <- sqrt(sum({mu[iter,] - mu[iter-1,]}^2))
               tol[2] <- abs(IE[iter] - IE[iter-1])
               if(tol[1] <= Stop.Tolerance[1]) {
                    expand <- TRUE
                    if(tol[2] <= Stop.Tolerance[2]) {
                         if(count > 0) break
                         count <- count + 1}
                    }
               else count <- 0
               }
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n")
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n")
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...")
          cat("\n##################################################\n")
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n")
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
               ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())
          Dev <- matrix(m.old[["Dev"]],1,1)
          LPworst <- m.old[["LP"]]
          parm.len <- length(parm)
          if(parm.len > 10)
               warning("The number of parameters may be too large.")
          mu <- rbind(parm)
          parm.old <- parm
          colnames(mu) <- Data[["parm.names"]]
          sigma <- sqrt(diag(Covar))
          rule <- SparseGrid(J=parm.len, K=K)
          X <- rule$nodes
          W <- rule$weights
          N <- length(W)
          count <- 0
          expand <- FALSE
          IE <- 1
          tol <- c(1,1)
          ### Begin Iterative Quadrature
          options(warn=-1)
          cat("\nIteration: 1,   K: ", K, ",   Nodes: ", N, sep="")
          for (iter in 2:Iterations) {
               ### New Quadrature Rule
               if({expand == TRUE} & {K < Kmax}) {
                    K <- min(K+1, Kmax)
                    rule <- SparseGrid(J=parm.len, K=K)
                    X <- rule$nodes
                    W <- rule$weights
                    N <- length(W)}
               expand <- FALSE
               cat("\nIteration: ", iter, ",   K: ", K, ",   Nodes: ", N,
                    sep="")
               LP <- rep(LPworst, N)
               LPw <- M <- Z <- matrix(0, N, parm.len)
               mu <- rbind(mu, 0)
               ### Evaluate at the Nodes
               Z <- matrix(mu[iter-1,], N, parm.len, byrow=TRUE) +
                    sqrt(2) * matrix(sigma, N, parm.len, byrow=TRUE) * X
               mod <- parLapply(cl, 1:N, function(x) Model(Z[x,], Data))
               LP <- as.vector(unlist(lapply(mod, function(x) x[["LP"]])))
               Z <- matrix(as.vector(unlist(lapply(mod,
                    function(x) x[["parm"]]))), N, parm.len, byrow=TRUE)
               LP[which(!is.finite(LP))] <- LPworst
               ### Weights
               LPw <- matrix(exp(LP - logadd(LP)), N, parm.len)
               LPw <- LPw / matrix(colSums(LPw), N, parm.len)
               LPw[which(!is.finite(LPw))] <- 0
               if(all(LPw == 0)) {
                    expand <- TRUE
                    LPw <- rep(.Machine$double.eps, N)}
               ### Adapt mu and Sigma
               mu[iter,] <- colSums(LPw * Z)
               temp <- cov.wt(Z, wt=LPw[,1])$cov
               if(all(is.finite(temp))) Covar <- as.symmetric.matrix(temp)
               else expand <- TRUE
               diag(Covar) <- abs(diag(Covar))
               diag(Covar)[which(diag(Covar) < .Machine$double.eps)] <- .Machine$double.eps
               M <- W * exp(X) * sqrt(2) *
                    matrix(sigma, N, parm.len, byrow=TRUE)
               sigma <- sqrt(diag(Covar))
               mod <- Model(mu[iter,], Data)
               Dev <- rbind(Dev, mod[["Dev"]])
               ### Accept Only Improved mu
               m.old <- Model(mu[iter-1,], Data)
               if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
                    mod[["Monitor"]])))) {
                    if(m.old[["LP"]] >= mod[["LP"]]) {
                         expand <- TRUE
                         mu[iter,] <- mu[iter-1,]}
                    }
               else {
                    expand <- TRUE
                    mu[iter,] <- mu[iter-1,]}
               ### Integration Error and Tolerance
               Mw <- M / matrix(colSums(M), N, parm.len, byrow=TRUE)
               IE <- c(IE, mean(abs(Mw - LPw)))
               tol[1] <- sqrt(sum({mu[iter,] - mu[iter-1,]}^2))
               tol[2] <- abs(IE[iter] - IE[iter-1])
               if(tol[1] <= Stop.Tolerance[1]) {
                    expand <- TRUE
                    if(tol[2] <= Stop.Tolerance[2]) {
                         if(count > 0) break
                         count <- count + 1}
                    }
               else count <- 0
               }
          }
     options(warn=0)
     ### Output
     IQ <- list(Covar=Covar, Dev=Dev, iter=iter, LPw=LPw, M=M,
          mu=mu, N=N, parm.old=parm.old, tol=tol, Z=Z)
     return(IQ)
     }
.iqcagh <- function(Model, parm, Data, Covar, Iterations, Stop.Tolerance,
     CPUs, m.old, N, Nmax, Packages, Dyn.libs)
     {
     if(CPUs == 1) {
          Dev <- matrix(m.old[["Dev"]],1,1)
          LPworst <- m.old[["LP"]]
          parm.len <- length(parm)
          mu <- rbind(parm)
          parm.old <- parm
          colnames(mu) <- Data[["parm.names"]]
          sigma <- sqrt(diag(Covar))
          LP <- rep(LPworst, N)
          LPw <- M <- Z <- matrix(0, N, parm.len)
          rule <- GaussHermiteQuadRule(N)
          count <- 0
          expand <- matrix(0, 4, parm.len)
          IE <- 1
          tol <- c(1,1)
          ### Begin Iterative Quadrature
          options(warn=-1)
          cat("\nIteration: 1,   Nodes: ", N, ",   LP:",
               round(m.old[["LP"]],3), sep="")
          for (iter in 2:Iterations) {
               cat("\nIteration: ", iter, ",   Nodes: ", N,
                    ",   LP:", round(m.old[["LP"]],3), sep="")
               mu <- rbind(mu, 0)
               ### New Quadarature Rule
               if({any(rowSums(expand) == parm.len)} & {N < Nmax}) {
                    N <- min(N + 1, Nmax)
                    LPw <- M <- Z <- matrix(0, N, parm.len)
                    LP <- rep(LPworst, N)
                    rule <- GaussHermiteQuadRule(N)}
               expand <- matrix(0, 4, parm.len)
               X <- matrix(rule$nodes, N, parm.len)
               W <- matrix(rule$weights, N, parm.len)
               for (j in sample(parm.len)) {
                    ### Evaluate at the Nodes
                    Z[,j] <- m.old[["parm"]][j] + sqrt(2)*sigma[j]*X[,j]
                    for (i in 1:N) {
                         prop <- m.old[["parm"]]
                         prop[j] <- Z[i,j]
                         m.new <- Model(prop, Data)
                         if(all(is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                              m.new[["Monitor"]]))))
                              LP[i] <- m.new[["LP"]]
                         Z[i,j] <- m.new[["parm"]][j]}
                    ### Correct for Constraints
                    X[,j] <- (Z[,j] - m.old[["parm"]][j]) / (sqrt(2)*sigma[j])
                    W[,j] <- (2^(N-1) * factorial(N) * sqrt(pi)) /
                         (N^2 * Hermite(X[,j], N-1, prob=FALSE)^2)
                    W[,j][which(!is.finite(W[,j]))] <- 0
                    if(all(W[,j] == 0)) W[,j] <- 1/N
                    ### Weights
                    LPw[,j] <- exp(LP - logadd(LP))
                    LPw[,j] <- LPw[,j] / sum(LPw[,j])
                    if(any(!is.finite(LPw[,j]))) {
                         LPw[,j] <- rep(0, N)
                         expand[1,j] <- 1}
                    ### Adapt mu and sigma
                    if(all(LPw[,j] == 0)) {
                         LPw[,j] <- .Machine$double.eps
                         mu[iter,j] <- m.old[["parm"]][j] + rnorm(1,0,1e-5)
                         expand[2,j] <- 1
                         }
                    else {
                         LPmax <- which.max(LP)[1]
                         if({LPmax == 1} | {LPmax == N})
                              mu[iter,j] <- Z[LPmax,j]
                         else mu[iter,j] <- sum(LPw[,j] * Z[,j])
                         M[,j] <- W[,j] * exp(X[,j]) * sqrt(2) * sigma[j]
                         sigma[j] <- min(max(sqrt(sum(LPw[,j] *
                              {sqrt(2)*sigma[j]*X[,j]}^2)),
                              .Machine$double.eps), 1e5)}
                    ### Accept Increases Only
                    m.new <- Model(mu[iter,], Data)
                    mu[iter,] <- m.new[["parm"]]
                    if(all(is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                         m.new[["Monitor"]])))) {
                         if(m.old[["LP"]] >= m.new[["LP"]]) {
                              expand[3,j] <- 1
                              mu[iter,] <- m.old[["parm"]]
                              }
                         else m.old <- m.new
                         }
                    else {
                         expand[3,j] <- 1
                         mu[iter,] <- m.old[["parm"]]
                         }
                    }
               Dev <- rbind(Dev, m.old[["Dev"]])
               ### Integration Error, Parameter Change, and Tolerance
               Mw <- M / matrix(colSums(M), N, parm.len, byrow=TRUE)
               IE <- c(IE, mean(abs(Mw - LPw)))
               tol[1] <- sqrt(sum({mu[iter,] - mu[iter-1,]}^2))
               tol[2] <- abs(IE[iter] - IE[iter-1])
               if(tol[1] <= Stop.Tolerance[1]) {
                    expand[4,] <- 1
                    if(tol[2] <= Stop.Tolerance[2]) {
                         if(count > 0) break
                         count <- count + 1}
                    }
               else count <- 0
               }
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n")
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n")
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...")
          cat("\n##################################################\n")
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n")
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
               ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())
          Dev <- matrix(m.old[["Dev"]],1,1)
          LPworst <- m.old[["LP"]]
          parm.len <- length(parm)
          mu <- rbind(parm)
          parm.old <- parm
          colnames(mu) <- Data[["parm.names"]]
          sigma <- sqrt(diag(Covar))
          LPbest <- rep(-Inf, parm.len)
          LP <- rep(LPworst, N)
          LPw <- M <- Z <- matrix(0, N, parm.len)
          rule <- GaussHermiteQuadRule(N)
          count <- 0
          expand <- matrix(0, 4, parm.len)
          IE <- 1
          tol <- c(1,1)
          ### Begin Iterative Quadrature
          options(warn=-1)
          cat("\nIteration: 1,   Nodes: ", N,
                    ",   LP:", round(m.old[["LP"]],3), sep="")
          for (iter in 2:Iterations) {
               cat("\nIteration: ", iter, ",   Nodes: ", N,
                    ",   LP:", round(m.old[["LP"]],3), sep="")
               mu <- rbind(mu, 0)
               ### New Quadarature Rule
               if({any(rowSums(expand) == parm.len)} & {N < Nmax}) {
                    N <- min(N + 1, Nmax)
                    LPw <- M <- Z <- matrix(0, N, parm.len)
                    LP <- rep(LPworst, N)
                    rule <- GaussHermiteQuadRule(N)}
               expand <- matrix(0, 4, parm.len)
               X <- matrix(rule$nodes, N, parm.len)
               W <- matrix(rule$weights, N, parm.len)
               for (j in sample(parm.len)) {
                    ### Evaluate at the Nodes
                    Z[,j] <- m.old[["parm"]][j] + sqrt(2)*sigma[j]*X[,j]
                    prop <- matrix(m.old[["parm"]], N, parm.len, byrow=TRUE)
                    prop[,j] <- Z[,j]
                    m.new <- parLapply(cl, 1:N, function(x) Model(prop[x,], Data))
                    LP <- as.vector(unlist(lapply(m.new, function(x) x[["LP"]])))
                    Z[,j] <- as.vector(unlist(lapply(m.new,
                         function(x) x[["parm"]][j])))
                    LP[which(!is.finite(LP))] <- LPworst
                    ### Correct for Constraints
                    X[,j] <- (Z[,j] - m.old[["parm"]][j]) / (sqrt(2)*sigma[j])
                    W[,j] <- (2^(N-1) * factorial(N) * sqrt(pi)) /
                         (N^2 * Hermite(X[,j], N-1, prob=FALSE)^2)
                    W[,j][which(!is.finite(W[,j]))] <- 0
                    if(all(W[,j] == 0)) W[,j] <- 1/N
                    ### Weights
                    LPw[,j] <- exp(LP - logadd(LP))
                    LPw[,j] <- LPw[,j] / sum(LPw[,j])
                    if(any(!is.finite(LPw[,j]))) {
                         LPw[,j] <- rep(0, N)
                         expand[1,j] <- 1}
                    ### Adapt mu and sigma
                    if(all(LPw[,j] == 0)) {
                         LPw[,j] <- .Machine$double.eps
                         mu[iter,j] <- m.old[["parm"]][j] + rnorm(1,0,1e-5)
                         expand[2,j] <- 1
                         }
                    else {
                         LPmax <- which.max(LP)[1]
                         if({LPmax == 1} | {LPmax == N})
                              mu[iter,j] <- Z[LPmax,j]
                         else mu[iter,j] <- sum(LPw[,j] * Z[,j])
                         M[,j] <- W[,j] * exp(X[,j]) * sqrt(2) * sigma[j]
                         sigma[j] <- min(max(sqrt(sum(LPw[,j] *
                              {sqrt(2)*sigma[j]*X[,j]}^2)),
                              .Machine$double.eps), 1e5)}
                    ### Accept Increases Only
                    m.new <- Model(mu[iter,], Data)
                    mu[iter,] <- m.new[["parm"]]
                    if(all(is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                         m.new[["Monitor"]])))) {
                         if(m.old[["LP"]] >= m.new[["LP"]]) {
                              expand[3,j] <- 1
                              mu[iter,] <- m.old[["parm"]]
                              }
                         else m.old <- m.new
                         }
                    else {
                         expand[3,j] <- 1
                         mu[iter,] <- m.old[["parm"]]
                         }
                    }
               Dev <- rbind(Dev, m.old[["Dev"]])
               ### Integration Error, Parameter Change, and Tolerance
               Mw <- M / matrix(colSums(M), N, parm.len, byrow=TRUE)
               IE <- c(IE, mean(abs(Mw - LPw)))
               tol[1] <- sqrt(sum({mu[iter,] - mu[iter-1,]}^2))
               tol[2] <- abs(IE[iter] - IE[iter-1])
               if(tol[1] <= Stop.Tolerance[1]) {
                    expand[4,] <- 1
                    if(tol[2] <= Stop.Tolerance[2]) {
                         if(count > 0) break
                         count <- count + 1}
                    }
               else count <- 0
               }
          }
     options(warn=0)
     ### Output
     IQ <- list(Covar=diag(sigma^2), Dev=Dev, iter=iter, LPw=LPw, M=M,
          mu=mu, N=N, parm.old=parm.old, tol=tol, Z=Z)
     return(IQ)
     }

#End
