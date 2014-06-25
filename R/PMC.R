###########################################################################
# PMC                                                                     #
#                                                                         #
# The purpose of the PMC function is update a model with Population       #
# Monte Carlo.                                                            #
###########################################################################

PMC <- function(Model, Data, Initial.Values, Covar=NULL, Iterations=10,
     Thinning=1, alpha=NULL, M=1, N=1000, nu=9, CPUs=1, Type="PSOCK")
     {
     cat("\nPMC was called on ", date(), "\n", sep="")
     time1 <- proc.time()
     pmc.call <- match.call()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n")
     if(missing(Model)) stop("A function must be entered for Model.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data))
          stop("A list containing data must be entered for Data.")
     if(is.null(Data[["mon.names"]])) stop("In Data, mon.names is NULL.")
     if(is.null(Data[["parm.names"]])) stop("In Data, parm.names is NULL.")
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n")}}}}
     M <- max(round(abs(M)), 1)
     N <- max(round(abs(N)), 1)
     if(missing(Initial.Values)) {
          cat("WARNING: Initial Values were not supplied.\n")
          Initial.Values <- matrix(0, M, length(Data[["parm.names"]]))}
     if(is.vector(Initial.Values)) {
          if(!identical(length(Initial.Values),
               length(Data[["parm.names"]]))) {
               cat("WARNING: The length of Initial Values differed from",
                    "Data$parm.names.\n")
          Initial.Values <- matrix(0, M, length(Data[["parm.names"]]))}}
     else {
          if(!identical(ncol(Initial.Values),
               length(Data[["parm.names"]]))) {
               cat("WARNING: Columns in Initial Values differed from",
                    "Data$parm.names.\n")
          Initial.Values <- matrix(0, M, length(Data[["parm.names"]]))}}
     if(any(!is.finite(Initial.Values))) {
          cat("WARNING: Initial Values contain non-finite values.\n")
          Initial.Values <- matrix(0, M, length(Data[["parm.names"]]))}
     if(is.vector(Initial.Values)) {
          LIV <- length(Initial.Values)
          Initial.Values <- matrix(Initial.Values, M, LIV, byrow=TRUE)}
     else if(nrow(Initial.Values) != M)
               stop("nrow(Initial.Values != M.")
     else LIV <- ncol(Initial.Values)
     ### Covar
     ScaleF <- 2.381204 * 2.381204 / LIV
     if(is.null(Covar))
          Covar <- array(diag(LIV)*ScaleF, dim=c(LIV,LIV,Iterations,M))
     else Covar <- array(Covar, dim=c(LIV, LIV, Iterations, M))
     ### Iterations and Thinning
     Iterations <- max(round(abs(Iterations)), 1)
     Thinning <- min(max(round(abs(Thinning)), 1), N*M/2)
     ### alpha
     if(is.null(alpha)) alpha <- matrix(rep(1/M, M), M, Iterations)
     else if(any(!is.finite(alpha))) {
          stop("\nWARNING: alpha had non-finite values. Creating alpha...\n.")
          alpha <- matrix(rep(1/M, M), M, Iterations)}
     if(is.vector(alpha)) alpha <- matrix(alpha, M, Iterations, byrow=TRUE)
     if(nrow(alpha) != M) stop("nrow(alpha) != M.")
     else if(ncol(alpha) != Iterations) stop("ncol(alpha) != Iterations.")
     if(any(colSums(alpha) != 1)) 
          alpha <- alpha / matrix(colSums(alpha), M, Iterations,
               byrow=TRUE)
     ### mu
     mu <- array(0, dim=c(Iterations, LIV, M))
     for (m in 1:M) {
          mu[,,m] <- matrix(Initial.Values[m,], Iterations, LIV,
               byrow=TRUE)}
     ### Miscellaneous
     nu <- max(abs(nu), 2.01)
     post <- array(0, dim=c(N, LIV, Iterations, M))
     LH <- LP <- array(0, dim=c(N, Iterations, M))
     LW <- matrix(0, N, Iterations)
     essn <- perp <- rep(0, Iterations)
     ##########################  Test the Model  ##########################
     for (m in 1:M) {
          M0 <- Model(as.vector(Initial.Values[m,]), Data)
          if(!is.list(M0)) stop("Model must return a list.")
          if(length(M0) != 5) stop("Model must return five components.")
          if(length(M0[["LP"]]) > 1) stop("Multiple joint posteriors exist!")
          if(!identical(length(M0[["Monitor"]]),
               length(Data[["mon.names"]])))
               stop("Length of mon.names differs from length of monitors.")
          if(any(!is.finite(c(M0[["LP"]],M0[["Dev"]],M0[["parm"]]))))
               stop("Model produces non-finite results.")}
     ### Looking for apply functions and for loops
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out, 1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of apply functions\n")
          cat("     were found in the Model specification. Sampling speed will\n")
          cat("     increase if apply functions are vectorized in R or coded\n")
          cat("     in a faster language such as C++ via the Rcpp package.\n")}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Sampling speed will\n")
          cat("     increase if for loops are vectorized in R or coded in a\n")
          cat("     faster language such as C++ via the Rcpp package.\n")}
     ######################  Laplace Approximation  #######################
     ### Sample Size of Data
     if(!is.null(Data[["n"]])) if(length(Data[["n"]]) == 1) NN <- Data[["n"]] 
     if(!is.null(Data[["N"]])) if(length(Data[["N"]]) == 1) NN <- Data[["N"]] 
     if(!is.null(Data[["y"]])) NN <- nrow(matrix(Data[["y"]]))
     if(!is.null(Data[["Y"]])) NN <- nrow(matrix(Data[["Y"]]))
     if(is.null(NN)) stop("Sample size of Data not found in n, N, y, or Y.")
     if({sum(abs(Initial.Values[1,]) == 0) == ncol(Initial.Values)} &
          {NN >= 5*ncol(Initial.Values)}) {
          cat("\nLaplace Approximation will be used on initial values.\n")
          Fit.LA <- LaplaceApproximation(Model, Initial.Values[1,], Data,
               Method="SPG", CovEst="Hessian", sir=FALSE)
          Covar <- array(Fit.LA$Covar, dim=c(LIV, LIV, Iterations, M))
          Initial.Values <- matrix(Fit.LA$Summary1[1:ncol(Initial.Values),1],
               M, LIV)}
     ### Covar symmetry and positive-definiteness
     for (m in 1:M) {
          if(!is.symmetric.matrix(Covar[,,1,m]))
               Covar[,,1,m] <- as.symmetric.matrix(Covar[,,1,m])
          if(!is.positive.definite(Covar[,,1,m]))
               Covar[,,1,m] <- as.positive.definite(Covar[,,1,m])}
     ###################  Prepare for Parallelization  ####################
     CPUs <- abs(round(CPUs))
     if(CPUs > 1) {
          detectedCores <- max(detectCores(),
               as.integer(Sys.getenv("NSLOTS")), na.rm=TRUE)
          cat("\n\nCPUs Detected:", detectedCores, "\n")
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n")
               CPUs <- detectedCores}}
     ############################  Begin PMC  #############################
     cat("\nPMC is beginning to update...\n")
     for (iter in 1:Iterations) {
          cat("Iteration: ", iter, sep="")
          bad <- FALSE
          ### Importance Sampling
          for (m in 1:M) {
               S <- nu/(nu-2) * Covar[,,iter,m]
               if(!is.symmetric.matrix(S))
                    S <- as.symmetric.matrix(S)
               if(!is.positive.definite(S))
                    S <- as.positive.definite(S)
               post[,,iter,m] <- rmvt(N, mu[iter,,m], S, nu)
               if(sum(!is.finite(post[,,iter,m]) > 0))
                    stop("Bad draws from importance distribution.")
               ### Non-Parallel Processing
               if(CPUs == 1) {
                    for (i in 1:N) {
                         mod <- Model(post[i,,iter,m], Data)
                         if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
                              mod[["parm"]]))))
                              M0 <- mod
                         LP[i,iter,m] <- M0[["LP"]]
                         post[i,,iter,m] <- M0[["parm"]]}
                    }
               else { ### Parallel Processing
                    cl <- makeCluster(CPUs, Type)
                    varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
                         ls(envir=parent.env(environment()))))
                    clusterExport(cl, varlist=varlist, envir=environment())
                    clusterSetRNGStream(cl)
                    mod <- parLapply(cl, 1:nrow(post),
                         function(x) Model(post[x,,iter,m], Data))
                    stopCluster(cl)
                    LP[,iter,m] <- unlist(lapply(mod,
                         function(x) x[["LP"]]))[1:nrow(post)]
                    post[,,iter,m] <- matrix(unlist(lapply(mod,
                         function(x) x[["parm"]])), dim(post)[1],
                         dim(post)[2], byrow=TRUE)}
               ### Proposal Sampling Distribution H(theta) ~ MVT
               LH[,iter,m] <- dmvt(post[,,iter,m], mu[iter,,m],
                    S, nu, log=TRUE)
               if(any(!is.finite(LH[,iter,m]))) {
                    if(iter > 1) LH[,iter,m] <- LH[,iter-1,m]
                    else stop("Sampling density H has non-finite values.")}}
          ### Weights
          LW[,iter] <- apply(LP[,iter,] -
               matrix(apply(matrix(alpha[,iter], N, M, byrow=TRUE) +
                    LH[,iter,], 1, logadd), N, M, byrow=TRUE), 1, logadd)
          if(all(!is.finite(LW[,iter]))) {
               bad <- TRUE; LW[,iter] <- 0}
          else if(any(!is.finite(LW[,iter]))) {
               w <- which(!is.finite(LW[,iter]))
               LW[w,iter] <- 1e-800}
          ### Normalize Weights
          LW[,iter] <- LW[,iter] - logadd(LW[,iter])
          if(any(!is.finite(LW[,iter]))){
               bad <- TRUE
               if(iter == 1) {
                    cat(", WARNING: Bad weights, setting to 1/N.\n")
                    LW[,iter] <- rep(log(1/N), N)}
               else {
                    cat(", WARNING: Bad weights, using last iteration.\n")
                    alpha[,iter] <- alpha[,iter-1]
                    for (m in 1:M) Covar[,,iter,m] <- Covar[,,iter-1,m]
                    LH[,iter,] <- LH[,iter-1,]
                    LP[,iter,] <- LP[,iter-1,]
                    post[,,iter,] <- post[,,iter-1,]
                    LW[,iter] <- LW[,iter-1]}}
          ### Convergence
          essn[iter] <- (1 / sum(exp(LW[,iter])^2)) / N
          perp[iter] <- exp(-sum(exp(LW[,iter]) * LW[,iter])) / N
          if(all(LW[,iter] == log(1/N))) essn[iter] <- perp[iter] <- 0
          if(bad == FALSE)
               cat(",   ESSN: ", round(essn[iter],5),
                    ",   Perplexity: ", round(perp[iter],5), "\n", sep="")
          ### Adaptation
          LpM <- matrix(log(alpha[,iter]), N, M, byrow=TRUE) + LH[,iter,]
          LpM <- LpM - matrix(apply(LpM, 1, logadd), N, M)
          if(iter < Iterations) {
               alpha[,iter+1] <- exp(apply(LW[,iter] + LpM, 2, logadd))
               if(any(!is.finite(alpha[,iter+1])))
                    alpha[,iter+1] <- ifelse(!is.finite(alpha[,iter+1]),
                         alpha[,iter], alpha[,iter+1])
                    alpha[which(alpha[,iter+1] < 0.001),iter+1] <- 0.001
                    alpha[,iter+1] <- alpha[,iter+1] / sum(alpha[,iter+1])
               for (m in 1:M) {
                    mu[iter+1,,m] <- colSums(exp(LW[,iter]) *
                         post[,,iter,m] * exp(LpM[,m])) / alpha[m,iter+1]
                    if(any(!is.finite(mu[iter+1,,m])))
                         mu[iter+1,,m] <- mu[iter,,m]}
               for (m in 1:M) {
                    Covar[,,iter+1,m] <- 0
                    for (i in 1:N) {
                         Covar[,,iter+1,m] <- Covar[,,iter+1,m] +
                         (exp(LW[i,iter] + LpM[i,m]) *
                         ((post[i,,iter,m] - mu[iter+1,,m]) %*%
                         t(post[i,,iter,m] - mu[iter+1,,m]))) /
                         alpha[m,iter+1]}
                    if(any(!is.finite(diag(Covar[,,iter+1,m]))))
                         Covar[,,iter+1,m] <- Covar[,,iter,m]
                    if(any(diag(Covar[,,iter+1,m] < 1e-100)))
                         Covar[,,iter+1,m] <- Covar[,,iter,m]
                    Covar[,,iter+1,m] <- as.symmetric.matrix(Covar[,,iter+1,m])
                    if(!is.positive.definite(Covar[,,iter+1,m]))
                         Covar[,,iter+1,m] <- Covar[,,iter,m]
                    if(bad == TRUE) Covar[,,iter+1,m] <- Covar[,,iter,m]}}
          }
     ### Combine Samples from Mixture Components
     Posterior2 <- post[,,Iterations,1]
     colnames(Posterior2) <- Data[["parm.names"]]
     if(M > 1) {
          for (m in 2:M) {
               if(alpha[m,Iterations] >= 0.002)
                    Posterior2 <- rbind(Posterior2,
                         post[,,Iterations,m])}}
     Posterior2 <- Thin(Posterior2, By=Thinning)
     ### Final Sampling
     cat("Final Sampling\n")
     Dev <- rep(0, nrow(Posterior2))
     Mon <- matrix(0, nrow(Posterior2), length(Data[["mon.names"]]))
     colnames(Mon) <- Data[["mon.names"]]
     for (i in 1:nrow(Posterior2)) {
          temp <- Model(Posterior2[i,], Data)
          Dev[i] <- temp[["Dev"]]
          Mon[i,] <- temp[["Monitor"]]}
     ### Posterior Summary Table
     cat("Creating Summaries\n")
     Summ <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ[,1] <- colMeans(Posterior2)
     Summ[,2] <- sqrt(.colVars(Posterior2))
     Summ[,3] <- 0
     Summ[,4] <- ESS(Posterior2)
     Summ[,5] <- apply(Posterior2, 2, quantile, c(0.025), na.rm=TRUE)
     Summ[,6] <- apply(Posterior2, 2, quantile, c(0.500), na.rm=TRUE)
     Summ[,7] <- apply(Posterior2, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(Posterior2)) {
          temp <- try(MCSE(Posterior2[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ[i,3] <- temp
          else Summ[i,3] <- MCSE(Posterior2[,i], method="sample.variance")}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)))
     if(inherits(temp, "try-error"))
          temp <- MCSE(as.vector(Dev), method="sample.variance")
     Deviance[3] <- temp
     Deviance[4] <- ESS(Dev)
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ <- rbind(Summ, Deviance)
     for (j in 1:ncol(Mon)) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(Mon[,j]), silent=TRUE)
          if(inherits(temp, "try-error")) 
               temp <- MCSE(Mon[,j], method="sample.variance")
          Monitor[3] <- temp
          Monitor[4] <- ESS(Mon[,j])
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.5,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ <- rbind(Summ, Monitor)
          rownames(Summ)[nrow(Summ)] <- Data[["mon.names"]][j]}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     cat("Estimating Log of the Marginal Likelihood\n")
     LML <- LML(theta=Posterior2,
          LL=as.vector(Dev)*(-1/2), method="NSIS")
     time2 <- proc.time()
     ### Compile Output
     cat("Creating Output\n")
     pmc.out <- (list(alpha=alpha,
          Call=pmc.call,
          Covar=Covar,
          Deviance=Dev,
          DIC=c(mean(as.vector(Dev)),
               var(as.vector(Dev))/2,
               mean(as.vector(Dev)) + var(as.vector(Dev))/2),
          ESSN=essn,
          Initial.Values=Initial.Values,
          Iterations=Iterations,
          LML=LML[[1]],
          M=M,
          Minutes=round(as.vector(time2[3] - time1[3]) / 60,2),
          Model=Model,
          N=N,
          nu=nu,
          Mu=mu,
          Monitor=Mon,
          Parameters=LIV,
          Perplexity=perp,
          Posterior1=post,
          Posterior2=Posterior2,
          Summary=Summ,
          Thinned.Samples=nrow(Posterior2),
          Thinning=Thinning,
          W=exp(LW)))
     class(pmc.out) <- "pmc"
     return(pmc.out)
     }

#End
