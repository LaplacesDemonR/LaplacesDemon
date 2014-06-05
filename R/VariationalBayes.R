###########################################################################
# VariationalBayes                                                        #
#                                                                         #
# The purpose of the VariationalBayes function is to estimate a model     #
# with a Variational Bayes algorithm.                                     #
###########################################################################

VariationalBayes <- function(Model, parm, Data, Covar=NULL,
     Interval=1.0E-6, Iterations=1000, Method="Salimans2", Samples=1000,
     sir=TRUE, Stop.Tolerance=1.0E-5, CPUs=1, Type="PSOCK")
     {
     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     VB.call <- match.call()
     if(missing(Model)) stop("Model is a required argument.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data)) stop("Data is a required argument.")
     if(missing(parm)) {
          cat("Initial values were not supplied, and\n")
          cat("have been set to zero prior to VariationalBayes().\n")
          parm <- rep(0, length(Data[["parm.names"]]))}
     LIV <- length(as.vector(parm))
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n")}}}}
     if({Interval <= 0} | {Interval > 1}) Interval <- 1.0E-6
     Iterations <- min(max(round(Iterations), 10), 1000000)
     "%!in%" <- function(x,table) return(match(x, table, nomatch=0) == 0)
     if(Method %!in% c("Salimans2"))
          stop("Method is unknown.")
     if(Stop.Tolerance <= 0) Stop.Tolerance <- 1.0E-5
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
     ### Sample Size of Data
     if(!is.null(Data[["n"]])) if(length(Data[["n"]]) == 1) N <- Data[["n"]] 
     if(!is.null(Data[["N"]])) if(length(Data[["N"]]) == 1) N <- Data[["N"]] 
     if(!is.null(Data[["y"]])) N <- nrow(matrix(Data[["y"]]))
     if(!is.null(Data[["Y"]])) N <- nrow(matrix(Data[["Y"]]))
     if(!is.null(N)) cat("Sample Size: ", N, "\n")
     else stop("Sample size of Data not found in n, N, y, or Y.")
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
     if(!identical(Model(m.old[["parm"]], Data)[["LP"]], m.old[["LP"]])) {
          cat("WARNING: LP differs when initial values are held constant.\n")
          cat("     Derivatives may be problematic if used.\n")}
     #####################  Begin Variational Bayes  #####################
     cat("Variational Bayes begins...\n")
     if(Method == "Salimans2") {
          VB <- .vbsalimans2(Model, parm, Data, Covar, Iterations,
               Interval, Stop.Tolerance, m.old)
          }
     Dev <- as.vector(VB$Dev)
     iter <- VB$iter
     parm.len <- LIV
     parm.new <- VB$parm.new
     parm.old <- VB$parm.old
     post <- VB$post
     Step.Size <- VB$Step.Size
     tol.new <- VB$tol.new
     VarCov <- VB$VarCov
     rm(VB)
     if(iter == 1) stop("VariationalBayes stopped at iteration 1.")
     if(tol.new <= Stop.Tolerance) converged <- TRUE
     else converged <- FALSE
     ### Column names to samples
     if(dim(post)[2] == length(Data[["parm.names"]]))
          dimnames(post) <- list(1:dim(post)[1], Data[["parm.names"]], 1:2)
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
     VB <- list(Call=VB.call,
          Converged=converged,
          Covar=VarCov,
          Deviance=as.vector(Dev),
          History=post,
          Initial.Values=parm,
          Iterations=iter,
          LML=LML[[1]],
          LP.Final=as.vector(Model(parm.new, Data)[["LP"]]),
          LP.Initial=m.old[["LP"]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60, 2),
          Monitor=Mon,
          Posterior=posterior,
          Step.Size.Final=Step.Size,
          Step.Size.Initial=Step.Size,
          Summary1=Summ1,
          Summary2=Summ2,
          Tolerance.Final=tol.new,
          Tolerance.Stop=Stop.Tolerance)
     class(VB) <- "vb"
     cat("Variational Bayes is finished.\n\n")
     return(VB)
     }
.vbsalimans2 <- function(Model, parm, Data, Covar, Iterations, Interval,
     Stop.Tolerance, m.old)
     {
     m.new <- m.old
     LIV <- length(parm)
     m <- parm
     if(is.null(Covar)) V <- diag(LIV) #Variance
     else {
          V <- as.positive.definite(Covar)
          if(nrow(V) != ncol(V)) stop("Covar is not square.")
          if(nrow(V) != LIV) V <- diag(LIV)
          else {
               V <- Covar
               diag(V) <- abs(diag(V))}}
     rm(Covar)
     z <- m #Guess of the mean
     P <- as.inverse(V) #Precision
     a <- rep(0, LIV)
     zbar <- rep(0, LIV)
     Pbar <- matrix(0, LIV, LIV)
     abar <- rep(0, LIV)
     w <- 1 / sqrt(Iterations) #step-size
     half1 <- Iterations / 2
     half2 <- 2 / Iterations
     post <- array(0, dim=c(Iterations, LIV, 2))
     post[1,,1] <- m
     post[1,,2] <- diag(V)
     Dev <- matrix(m.old[["Dev"]],1,1)
     ### Stochastic Approximation
     for (iter in 1:Iterations) {
          mbar <- mbar.last <- m #mean
          Vbar <- Vbar.last <- V #variance
          m.old <- m.new #model
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### Draw a sample from the approximating distribution q
          xstar <- m
          while(identical(xstar, m)) {
               xstar <- try(m +
                    as.vector(matrix(rnorm(LIV),1,LIV) %*% chol(V)),
                    silent=TRUE)
               if(inherits(xstar, "try-error"))
                    xstar <- rnorm(LIV,m,abs(diag(V)))
               m.temp <- try(Model(xstar, Data), silent=TRUE)
               if(inherits(xstar, "try-error")) xstar <- m
               else xstar <- m.temp[["parm"]]
               if(any(!is.finite(c(m.temp[["LP"]], m.temp[["Dev"]],
                    m.temp[["Monitor"]]))))
                    xstar <- m}
          ### Gradient and Hessian
          g <- partial(Model, xstar, Data, Interval=Interval)
          H <- Hessian(Model, xstar, Data, Interval=Interval)
          ### Stochastic Approx.
          a <- (1 - w)*a + w*g
          P <- (1 - w)*P - w*H
          z <- (1 - w)*z + w*xstar
          if(any(!is.finite(P))) P <- diag(LIV)
          if(!is.symmetric.matrix(P))
               P <- as.symmetric.matrix(P)
          diag(P) <- abs(diag(P))
          if(!is.positive.definite(P))
               P <- as.positive.definite(P)
          V <- as.inverse(P)
          m <- as.vector(V %*% a + z)
          ### Evaluate Proposal
          m.new <- try(Model(m, Data), silent=TRUE)
          if(inherits(m.new, "try-error")) m.new <- m.old
          else if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
               m.new[["Monitor"]]))))
               m.new <- m.old
          else if(log(runif(1)) >= (m.new[["LP"]] - m.old[["LP"]]))
               m.new <- m.old
          m <- m.new[["parm"]]
          ### Storage
          post[iter,,1] <- m
          post[iter,,2] <- diag(V)
          Dev <- rbind(Dev, m.new[["Dev"]])
          ### Do averaging if over half-way
          if(iter > half1) {
               abar <- abar + half2*g
               Pbar <- Pbar - half2*H
               zbar <- zbar + half2*xstar
               if(any(!is.finite(Pbar))) Pbar <- diag(LIV)
               if(!is.symmetric.matrix(Pbar))
                    Pbar <- as.symmetric.matrix(Pbar)
               diag(Pbar) <- abs(diag(Pbar))
               if(!is.positive.definite(Pbar))
                    Pbar <- as.positive.definite(Pbar)
               Vbar <- as.inverse(Pbar)
               mbar <- as.vector(Vbar %*% abar + zbar)
               mbar <- Model(mbar, Data)[["parm"]]
               ### Tolerance
               tol.new <- sum(sqrt(sum({mbar - mbar.last} *
                    {mbar - mbar.last})),
                    sqrt(sum({diag(Vbar) - diag(Vbar.last)} *
                    {diag(Vbar) - diag(Vbar.last)})))
               if(tol.new <= Stop.Tolerance) {
                    post <- post[1:iter,,]
                    break}
               }
          }
     Dev <- Dev[-1,]
     #mbar and Vbar should be returned, but not if wildly different...
     LB <- mbar - 3*diag(Vbar)
     UB <- mbar + 3*diag(Vbar)
     if(any((m < LB) | (m > UB))) {
          mbar <- m
          Vbar <- diag(LIV) * diag(V)}
     ### Output
     VB <- list(Dev=Dev, iter=iter, parm.new=mbar, parm.old=parm,
          post=post, Step.Size=w, tol.new=tol.new, VarCov=Vbar)
     return(VB)
     }

#End
