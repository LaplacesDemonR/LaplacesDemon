###########################################################################
# LaplaceApproximation                                                    #
#                                                                         #
# The purpose of the LaplaceApproximation function is to maximize the     #
# logarithm of the unnormalized joint posterior distribution of a         #
# Bayesian model with one of many optimization algorithms.                #
###########################################################################

LaplaceApproximation <- function(Model, parm, Data, Interval=1.0E-6,
     Iterations=100, Method="SPG", Samples=1000, CovEst="Hessian",
     sir=TRUE, Stop.Tolerance=1.0E-5, CPUs=1, Type="PSOCK")
     {
     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     LA.call <- match.call()
     if(missing(Model)) stop("Model is a required argument.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data)) stop("Data is a required argument.")
     if(missing(parm)) {
          cat("Initial values were not supplied, and\n")
          cat("have been set to zero prior to LaplaceApproximation().\n")
          parm <- rep(0, length(Data[["parm.names"]]))}
     if(is.null(Data[["mon.names"]]))
          stop("In Data, mon.names is NULL.")
     if(is.null(Data[["parm.names"]]))
          stop("In Data, parm.names is NULL.")
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
     if(Method %!in% c("AGA","BFGS","BHHH","CG","DFP","HAR","HJ","LBFGS",
          "LM","NM","NR","PSO","Rprop","SGD","SOMA","SPG","SR1","TR"))
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
          cat("     increase if apply functions are vectorized in R or coded\n")}
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
     if(Method == "SGD") if(!is.null(Data[["Nr"]])) N <- Data[["Nr"]]
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
     ####################  Begin Laplace Approximation  ###################
     cat("Laplace Approximation begins...\n")
     if(Method == "AGA") {
          LA <- .laaga(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "BFGS") {
          LA <- .labfgs(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "BHHH") {
          LA <- .labhhh(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "CG") {
          LA <- .lacg(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
          }
     else if(Method == "DFP") {
          LA <- .ladfp(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "HAR") {
          LA <- .lahar(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
          }
     else if(Method == "HJ") {
          LA <- .lahj(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
               m.old)
          }
     else if(Method == "LBFGS") {
          LA <- .lalbfgs(Model, parm, Data, Iterations, Stop.Tolerance,
               m.old)
          }
     else if(Method == "LM") {
          LA <- .lalm(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
          }
     else if(Method == "NM") {
          LA <- .lanm(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
          }
     else if(Method == "NR") {
          LA <- .lanr(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "PSO") {
          LA <- .lapso(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
          }
     else if(Method == "Rprop") {
          LA <- .larprop(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "SGD") {
          LA <- .lasgd(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
          }
     else if(Method == "SOMA") {
          LA <- .lasoma(Model, parm, Data,
               options=list(maxiter=Iterations, tol=Stop.Tolerance))
          }
     else if(Method == "SPG") {
          LA <- .laspg(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "SR1") {
          LA <- .lasr1(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)
          }
     else if(Method == "TR") {
          LA <- .latr(Model, parm, Data, Iterations, m.old)
          }
     Dev <- as.vector(LA$Dev)
     if(is.null(LA$H)) H <- FALSE
     else H <- LA$H
     iter <- LA$iter
     parm.len <- LA$parm.len
     parm.new <- LA$parm.new
     parm.old <- LA$parm.old
     post <- LA$post
     Step.Size <- LA$Step.Size
     tol.new <- LA$tol.new
     rm(LA)
     if(iter == 1) stop("LaplaceApproximation stopped at iteration 1.")
     if(tol.new <= Stop.Tolerance) converged <- TRUE
     else converged <- FALSE
     ### Column names to samples
     if(ncol(post) == length(Data[["parm.names"]]))
          colnames(post) <- Data[["parm.names"]]
     rownames(post) <- 1:nrow(post)
     ########################  Covariance Matirx  #########################
     cat("Estimating the Covariance Matrix\n")
     if(all(H == FALSE)) {
          VarCov <- CovEstim(Model, parm.new, Data, Method=CovEst)
          }
     else {
          VarCov <- -as.inverse(as.symmetric.matrix(H))
          diag(VarCov) <- abs(diag(VarCov))
          }
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
          c("Mode","SD","LB","UB")))
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
                    c("Mode","SD","MCSE","ESS","LB","Median","UB")))
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
     LA <- list(Call=LA.call,
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
          Step.Size.Initial=1,
          Summary1=Summ1,
          Summary2=Summ2,
          Tolerance.Final=tol.new,
          Tolerance.Stop=Stop.Tolerance)
     class(LA) <- "laplace"
     cat("Laplace Approximation is finished.\n\n")
     return(LA)
     }
.laaga <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old)
     {
     alpha.star <- 0.234
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.len <- length(as.vector(parm))
     parm.new <- parm.old <- m.old[["parm"]]
     names(parm.new) <- names(parm.old) <- Data[["parm.names"]]
     tol.new <- 1
     Step.Size  <- 1 / parm.len
     post <- matrix(parm.new, 1, parm.len)
     for (iter in 1:Iterations) {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Approximate Truncated Gradient
          approx.grad <- partial(Model, parm.old, Data, Interval)
          approx.grad <- interval(approx.grad, -1000, 1000, reflect=FALSE)
          ### Proposal
          parm.new <- parm.old + Step.Size * approx.grad
          if(any(!is.finite(parm.new))) parm.new <- parm.old
          m.new <- Model(parm.new, Data)
          tol.new <- max(sqrt(sum(approx.grad^2)),
               sqrt(sum({parm.new - parm.old}^2)))
          if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
               m.new[["Monitor"]])))) {
               m.new <- m.old
               parm.new <- parm.old}
          ### Accept/Reject and Adapt
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               parm.new <- m.new[["parm"]]
               Step.Size <- Step.Size + (Step.Size / (alpha.star *
                    (1 - alpha.star))) * (1 - alpha.star) / iter
               }
          else {
               m.new <- m.old
               parm.new <- parm.old
               Step.Size <- abs(Step.Size - (Step.Size / (alpha.star *
                    (1 - alpha.star))) * alpha.star / iter)
               }
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
.labfgs <- function(Model, parm, Data, Interval, Iterations,
     Stop.Tolerance, m.old) {
     m.new <- m.old
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.old <- parm
     parm.len <- length(as.vector(parm))
     post <- matrix(m.old[["parm"]], 1, parm.len)
     tol.new <- 1
     keepgoing <- TRUE
     g.old <- g.new <- rep(0, parm.len)
     B <- diag(parm.len) #Approximate Hessian
     options(warn=-1)
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Gradient and Direction p
          g.old <- g.new
          g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
          p <- as.vector(tcrossprod(g.new, -B))
          p[which(!is.finite(p))] <- 0
          ### Step-size Line Search
          Step.Size <- 0.8
          changed <- TRUE
          while(m.new[["LP"]] <= m.old[["LP"]] & changed == TRUE) {
               Step.Size <- Step.Size * 0.2
               s <- Step.Size*p
               prop <- m.old[["parm"]] + s
               changed <- !identical(m.old[["parm"]], prop)
               m.new <- Model(prop, Data)
               if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                    m.new[["Monitor"]]))))
                    m.new <- m.old
               }
          ### BFGS Update to Approximate Hessian B
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               g <- g.new - g.old
               CC <- sum(s*g) #Curvature condition
               if(CC > 0) {
                    y <- as.vector(crossprod(B, g))
                    D <- as.double(1 + crossprod(g, y)/CC)
                    B <- B - (tcrossprod(s, y) + tcrossprod(y, s) -
                         D * tcrossprod(s, s))/CC}
               if(any(!is.finite(B))) B <- diag(parm.len)
               }
          ### Storage
          post <- rbind(post, m.old[["parm"]])
          Dev <- rbind(Dev, m.old[["Dev"]])
          ### Tolerance
          tol.new <- sqrt(sum(s*s))
          if(keepgoing == FALSE) tol.new <- 0
          if(tol.new <= Stop.Tolerance) break
          }
     options(warn=0)
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len,
          parm.new=m.old[["parm"]], parm.old=parm.old, post=post,
          Step.Size=Step.Size, tol.new=tol.new)
     return(LA)
     }
.labhhh <- function(Model, parm, Data, Interval, Iterations=100,
     Stop.Tolerance, m.old)
     {
     ### Check data for X and y or Y
     if(is.null(Data[["X"]])) stop("X is required in the data.")
     y <- TRUE
     if(is.null(Data[["y"]])) {
          y <- FALSE
          if(is.null(Data[["Y"]])) stop("y or Y is required in the data.")}
     if(y == TRUE) {
          if(length(Data[["y"]]) != nrow(Data[["X"]]))
               stop("length of y differs from rows in X.")
          }
     else {
          if(nrow(Data[["Y"]]) != nrow(Data[["X"]]))
               stop("The number of rows differs in y and X.")}
     m.new <- m.old
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.old <- parm
     parm.len <- length(parm)
     post <- matrix(parm, 1, parm.len)
     tol.new <- 1
     options(warn=-1)
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Gradient p, OPG A from gradient g, and direction delta
          p <- partial(Model, m.old[["parm"]], Data, Interval)
          A <- matrix(0, parm.len, parm.len)
          for (i in 1:nrow(Data[["X"]])) {
               Data.temp <- Data
               Data.temp$X <- Data.temp$X[i,,drop=FALSE]
               if(y == TRUE) Data.temp$y <- Data.temp$y[i]
               else Data.temp$Y <- Data.temp$Y[i,]
               g <- partial(Model, m.old[["parm"]], Data.temp, Interval)
               A <- A + tcrossprod(g,g)}
          A <- -as.inverse(as.symmetric.matrix(A))
          delta <- as.vector(tcrossprod(p, -A))
          ### Step-size Line Search
          Step.Size <- 0.8
          changed <- TRUE
          while(m.new[["LP"]] <= m.old[["LP"]] & changed == TRUE) {
               Step.Size <- Step.Size * 0.2
               s <- Step.Size*delta
               prop <- m.old[["parm"]] + s
               changed <- !identical(m.old[["parm"]], prop)
               m.new <- Model(prop, Data)
               if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                    m.new[["Monitor"]]))))
                    m.new <- m.old
               }
          if(m.new[["LP"]] > m.old[["LP"]]) m.old <- m.new
          ### Storage
          post <- rbind(post, m.old[["parm"]])
          Dev <- rbind(Dev, m.old[["Dev"]])
          ### Tolerance
          tol.new <- sqrt(sum(delta*delta))
          if(tol.new < Stop.Tolerance) break
          }
     options(warn=0)
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len,
          parm.new=m.old[["parm"]], parm.old=parm.old, post=post,
          Step.Size=Step.Size, tol.new=tol.new)
     return(LA)
     }
.lacg <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     m.best <- m.new <- m.old
     parm.len <- length(parm)
     tol <- Stop.Tolerance
     tol.new <- 1
     Dev <- matrix(m.old[["Dev"]],1,1)
     post <- matrix(parm, 1, parm.len)
     eps <- 1e-7
     bvec <- parm.old <- parm
     n <- length(bvec)
     ig <- 0 # count gradient evaluations
     stepredn <- 0.15 # Step reduction in line search
     acctol <- 1e-04 # acceptable point tolerance
     reltest <- 100 # relative equality test
     accpoint <- as.logical(FALSE)  # so far do not have an acceptable point
     cyclimit <- min(2.5 * n, 10 + sqrt(n))
     fmin <- f <- m.old[["LP"]] * -1
     bvec <- m.old[["parm"]]
     keepgoing <- TRUE
     oldstep <- steplength <- 0.8
     fdiff <- NA # initially no decrease
     cycle <- 0 # !! cycle loop counter
     options(warn=-1)
     while (keepgoing == TRUE) {
          t <- as.vector(rep(0, n))  # zero step vector
          c <- t  # zero 'last' gradient
          while ({keepgoing == TRUE} && {cycle < cyclimit}) {
               cycle <- cycle + 1
               parm <- bvec
               ig <- ig + 1
               post <- rbind(post, m.best[["parm"]])
               Dev <- rbind(Dev, m.best[["Dev"]])
               ### Print Status
               if(ig %% round(Iterations / 10) == 0)
                    cat("Iteration: ", ig, " of ", Iterations,
                         ",   LP: ", round(m.old[["LP"]],1), "\n")
               if(ig > Iterations) {
                    ig <- Iterations
                    LA <- list(Dev=Dev, iter=ig, parm.len=parm.len,
                         parm.new=parm, parm.old=parm.old, post=post,
                         Step.Size=steplength, tol.new=tol.new)
                    return(LA)}
               g <- partial(Model, bvec, Data) * -1
               g1 <- sum(g * (g - c))
               g2 <- sum(t * (g - c))
               gradsqr <- sum(g * g)
               c <- g
               g3 <- 1
               if(gradsqr > tol * (abs(fmin) + reltest)) {
                    if(g2 > 0) {
                         betaDY <- gradsqr / g2
                         betaHS <- g1 / g2
                         g3 <- max(0, min(betaHS, betaDY))}
                    }
               else {
                    keepgoing <- FALSE
                    tol.new <- gradsqr
                    break}
               if(g3 == 0 || cycle >= cyclimit) {
                    fdiff <- NA
                    cycle <- 0
                    break
                    }
               else {
                    t <- t * g3 - g
                    gradproj <- sum(t * g)
                    ### Line search
                    OKpoint <- FALSE
                    steplength <- oldstep * 1.5
                    f <- fmin
                    changed <- TRUE
                    while ({f >= fmin} && {changed == TRUE}) {
                         bvec <- parm + steplength * t
                         changed <- !identical((bvec + reltest),
                              (parm + reltest))
                         if(changed == TRUE) {
                              m.old <- m.new
                              m.new <- Model(bvec, Data)
                              if(any(!is.finite(c(m.new[["LP"]],
                                   m.new[["Dev"]], m.new[["Monitor"]])))) {
                                   f <- fmin + 1
                                   m.new <- m.old
                                   }
                              else {
                                   if(m.new[["LP"]] > m.best[["LP"]])
                                        m.best <- m.new
                                   f <- m.new[["LP"]] * -1
                                   tol.new <- max(sqrt(sum(g*g)),
                                        sqrt(sum({m.new[["parm"]] -
                                        m.old[["parm"]]}^2)))
                                   }
                              bvec <- m.new[["parm"]]
                              if(f < fmin) f1 <- f
                              else {
                                   savestep <-steplength
                                   steplength <- steplength * stepredn
                                   if(steplength >= savestep) changed <-FALSE}}}
                    changed1 <- changed
                    if(changed1 == TRUE) {
                         newstep <- 2 * (f - fmin - gradproj * steplength)
                         if(newstep > 0)
                              newstep <- -(gradproj * steplength *
                                   steplength / newstep)
                         bvec <- parm + newstep * t
                         changed <- !identical((bvec + reltest), 
                              (parm + reltest))
                         if(changed == TRUE) {
                              m.old <- m.new
                              m.new <- Model(bvec, Data)
                              if(any(!is.finite(c(m.new[["LP"]],
                                   m.new[["Dev"]], m.new[["Monitor"]])))) {
                                   f <- fmin + 1
                                   m.new <- m.old
                                   }
                              else {
                                   if(m.new[["LP"]] > m.best[["LP"]])
                                        m.best <- m.new
                                   f <- m.new[["LP"]] * -1
                                   tol.new <- max(sqrt(sum(g*g)),
                                        sqrt(sum({m.new[["parm"]] -
                                        m.old[["parm"]]}^2)))
                                   }
                              bvec <- m.new[["parm"]]}
                         if(f < min(fmin, f1)) {
                              OKpoint <- TRUE
                              accpoint <- (f <= fmin + gradproj *
                                   newstep * acctol)
                              fdiff <- (fmin - f)
                              fmin <- f
                              oldstep <- newstep
                              }
                         else {
                              if(f1 < fmin) {
                                   bvec <- parm + steplength * t
                                   accpoint <- (f1 <= fmin + gradproj * 
                                        steplength * acctol)
                                   OKpoint <- TRUE
                                   fdiff <- (fmin - f1)
                                   fmin <- f1
                                   oldstep <- steplength
                                   }
                              else {
                                   fdiff <- NA
                                   accpoint <- FALSE}
                              }
                         if(accpoint == FALSE) {
                              keepgoing <- FALSE
                              break}
                         }
                    else {
                         if(cycle == 1) {
                              keekpgoing <- FALSE
                              break}}}
               }
          if(oldstep < acctol) oldstep <- acctol
          if(oldstep > 1) oldstep <- 1
          }
     options(warn=0)
     Dev <- Dev[-1,]; post <- post[-1,]
     if({ig < Iterations} & {tol.new > Stop.Tolerance})
          tol.new <- Stop.Tolerance
     ### Output
     LA <- list(Dev=Dev, iter=ig, parm.len=parm.len,
          parm.new=m.best[["parm"]], parm.old=parm.old, post=post,
          Step.Size=steplength, tol.new=tol.new)
     return(LA)
     }
.ladfp <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old) {
     m.new <- m.old
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.old <- parm
     parm.len <- length(as.vector(parm))
     post <- matrix(m.old[["parm"]], 1, parm.len)
     tol.new <- 1
     keepgoing <- TRUE
     g.old <- g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
     B <- Iden <- diag(parm.len) #Approximate Hessian
     options(warn=-1)
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Gradient and Direction p
          g.old <- g.new
          p <- as.vector(tcrossprod(g.old, -B))
          p[which(!is.finite(p))] <- 0
          ### Step-size Line Search
          Step.Size <- 0.8
          changed <- TRUE
          while(m.new[["LP"]] <= m.old[["LP"]] & changed == TRUE) {
               Step.Size <- Step.Size * 0.2
               s <- Step.Size*p
               prop <- m.old[["parm"]] + s
               changed <- !identical(m.old[["parm"]], prop)
               m.new <- Model(prop, Data)
               if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                    m.new[["Monitor"]]))))
                    m.new <- m.old
               }
          g.new <- -1*partial(Model, m.new[["parm"]], Data, Interval)
          ### DFP Update to Approximate Hessian B
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               g <- g.new - g.old
               if(sum(s*g) > 0) {
                    denom <- as.vector(t(g) %*% s)
                    B <- (Iden - ((g %*% t(s)) / denom)) %*% B %*%
                         (Iden - ((s %*% t(g)) / denom)) +
                         ((g %*% t(g)) / denom)
                    if(any(!is.finite(B))) B <- diag(parm.len)}}
          ### Storage
          post <- rbind(post, m.old[["parm"]])
          Dev <- rbind(Dev, m.old[["Dev"]])
          ### Tolerance
          tol.new <- sqrt(sum(g.new*g.new))
          if(keepgoing == FALSE) tol.new <- 0
          if(tol.new <= Stop.Tolerance) break
          }
     options(warn=0)
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len,
          parm.new=m.old[["parm"]], parm.old=parm.old, post=post,
          Step.Size=Step.Size, tol.new=tol.new)
     return(LA)
     }
.lahar <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     alpha.star <- 0.234
     tau <- 1
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.len <- length(as.vector(parm))
     parm.new <- parm.old <- parm
     names(parm.new) <- names(parm.old) <- Data[["parm.names"]]
     tol.new <- 1
     post <- matrix(parm.new, 1, parm.len)
     for (iter in 1:Iterations) {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Propose new parameter values
          theta <- rnorm(parm.len)
          d <- theta / sqrt(sum(theta*theta))
          Step.Size <- runif(1,0,tau)
          parm.new <- parm.old + Step.Size * d
          m.new <- Model(parm.new, Data)
          tol.new <- sqrt(sum({m.new[["parm"]] - parm.old}^2))
          if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
               m.new[["Monitor"]]))))
               m.new <- m.old
          ### Accept/Reject and Adapt
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               parm.new <- m.new[["parm"]]
               tau <- tau + (tau / (alpha.star *
                    (1 - alpha.star))) * (1 - alpha.star) / iter
               }
          else {
               m.new <- m.old
               parm.new <- parm.old
               tau <- abs(tau - (tau / (alpha.star *
                    (1 - alpha.star))) * alpha.star / iter)
               }
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
.lahj <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     n <- length(parm)
     post <- matrix(parm, 1, n)
     if(n == 1)
          stop("For univariate functions use some different method.")
     tol.new <- 1
     f <- function(x, Data) {
          fun <- Model(x, Data)
          fun[["LP"]] <- fun[["LP"]] * -1
          return(fun)}
     steps  <- 2^c(-(0:(Iterations-1)))
     dir <- diag(1, n, n)
     x <- parm.old <- parm
     fx <- f(x, Data)
     fcount <- 1
     ###  Search with a single scale
     .lahjsearch <- function(xb, f, h, dir, fcount, Data)
          {
          x  <- xb
          xc <- x
          sf <- 0
          finc <- 0
          hje  <- .lahjexplore(xb, xc, f, h, dir, Data=Data)
          x    <- hje$x
          fx   <- hje$fx
          sf   <- hje$sf
          finc <- finc + hje$numf
          ### Pattern move
          while (sf == 1) {
               d  <- x - xb
               xb <- x
               xc <- x + d
               fb <- fx
               hje  <- .lahjexplore(xb, xc, f, h, dir, fb, Data)
               x    <- hje$x
               fx   <- hje$fx
               sf   <- hje$sf
               finc <- finc + hje$numf
               if(sf == 0) {  # pattern move failed
                    hje  <- .lahjexplore(xb, xb, f, h, dir, fb, Data)
                    x    <- hje$x
                    fx   <- hje$fx
                    sf   <- hje$sf
                    finc <- finc + hje$numf}
               }
          return(list(x=x, fx=fx, sf=sf, finc=finc))
          }
     ###  Exploratory move
     .lahjexplore <- function(xb, xc, f, h, dir, fbold, Data)
          {
          n <- length(xb)
          x <- xb
          if(missing(fbold)) {
               fb <- f(x, Data)
               x <- fb[["parm"]]
               numf <- 1
               }
          else {
               fb <- fbold
               numf <- 0}
          fx <- fb
          xt <- xc
          sf <- 0 # do we find a better point ?
          dirh <- h * dir
          fbold <- fx
          for (k in sample.int(n, n)) { # resample orthogonal directions
               p <- xt + dirh[, k]
               ft <- f(p, Data)
               p <- ft[["parm"]]
               numf <- numf + 1
               if(ft[["LP"]] >= fb[["LP"]]) {
                    p <- xt - dirh[, k]
                    ft <- f(p, Data)
                    p <- ft[["parm"]]
                    numf <- numf + 1
                    }
               else {
                    sf <- 1
                    xt <- p
                    fb <- ft}
               }
          if(sf == 1) {
               x <- xt
               fx <- fb}
          return(list(x=x, fx=fx, sf=sf, numf=numf))
          }
     ### Start the main loop
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(fx[["LP"]],1), "\n")
          hjs <- .lahjsearch(x, f, steps[iter], dir, fcount, Data)
          x <- hjs$x
          fx <- hjs$fx
          sf <- hjs$sf
          fcount <- fcount + hjs$finc
          post <- rbind(post, x)
          Dev <- rbind(Dev, fx[["Dev"]])
          tol.new <- sqrt(sum({fx[["parm"]] - parm.old}^2))
          if(tol.new <= Stop.Tolerance) break
          parm.old <- x
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     LA <- list(Dev=Dev, iter=iter, parm.len=n, parm.new=x,
          parm.old=parm, post=post, Step.Size=steps[iter],
          tol.new=tol.new)
     return(LA)
     }
.lalbfgs <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.len <- length(as.vector(parm))
     parm.new <- parm.old <- m.old[["parm"]]
     names(parm.new) <- names(parm.old) <- Data[["parm.names"]]
     tol.new <- 1
     post <- matrix(parm.new, 1, parm.len)
     ModelWrapper <- function(parm.new) {
          out <- Model(parm.new, Data)[["LP"]]
          return(out)
          }
     for (iter in 1:Iterations) {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### LBFGS
          Fit <- optim(par=parm.new, fn=ModelWrapper,
               method="L-BFGS-B", control=list(fnscale=-1, maxit=1))
          m.new <- Model(Fit$par, Data)
          tol.new <- sqrt(sum({m.new[["parm"]] - parm.old}^2))
          if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
               m.new[["Monitor"]]))) | {m.new[["LP"]] < m.old[["LP"]]})
               m.new <- m.old
          m.old <- m.new
          parm.new <- m.new[["parm"]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=0,
          tol.new=tol.new)
     return(LA)
     }
.lalm <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     Norm <- function(x, p=2) {
          stopifnot(is.numeric(x) || is.complex(x), is.numeric(p),
               length(p) == 1)
          if(p > -Inf && p < Inf) return(sum(abs(x)^p)^(1/p))
          else if(p ==  Inf) return(max(abs(x)))
          else if(p == -Inf) return(min(abs(x)))
          else return(NULL)
          }
     Dev <- matrix(m.old[["Dev"]],1,1)
     tau <- 1e-3 ### Damping constant
     tolg <- 1e-8
     n <- length(parm)
     m <- 1
     x <- xnew <- parm
     mod <- modbest <- m.old
     r <- r.adj <- mod[["LP"]]
     if(r.adj > 0) r.adj <- r.adj * -1
     dold <- dnew <- mod[["Dev"]]
     x <- parm <- mod[["parm"]]
     post <- matrix(parm, 1, n)
     f <- 0.5 * r * r
     J <- rbind(partial(Model, x, Data))
     g <- t(J) %*% r.adj
     ng <- max(abs(g))
     A <- t(J) %*% J
     mu <- tau * max(diag(A)) ### Damping parameter
     nu <- 2
     nh <- 0
     iter <- 1
     while (iter < Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(mod[["LP"]],1), "\n")
          iter <- iter + 1
          R <- chol(A + mu*diag(n))
          h <- c(-t(g) %*% chol2inv(R))
          bad <- !is.finite(h)
          h[bad] <- 1e-10 * J[1,bad] #Small gradient if ill-conditioned
          nh <- Norm(h)
          xnew <- x + h
          h <- xnew - x
          dL <- sum(h*(mu*h - g)) / 2
          mod <- Model(xnew, Data)
          if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
               mod[["Monitor"]])))) {
               if(mod[["LP"]] > modbest[["LP"]]) modbest <- mod
               rn <- mod[["LP"]]
               xnew <- mod[["parm"]]
               }
          else {
               rn <- r
               xnew <- x}
          fn <- 0.5 * rn * rn
          Jn <- rbind(partial(Model, xnew, Data))
          df <- f - fn
          if(rn > 0) df <- fn - f
          if(dL > 0 && df > 0) {
               tol.new <- sqrt(sum({xnew - x}^2))
               x <- xnew
               f <- fn
               J <- Jn
               r <- r.adj <- rn
               if(r.adj > 0) r.adj <- r.adj * -1
               A <- t(J) %*% J
               g <- t(J) %*% r.adj
               ng <- Norm(g, Inf)
               mu <- mu * max(1/3, 1 - (2*df/dL - 1)^3)
               nu <- 2}
          else {mu <- mu*nu
                nu <- 2*nu}
          post <- rbind(post, modbest[["parm"]])
          Dev <- rbind(Dev, modbest[["Dev"]])
          if(ng <= Stop.Tolerance) {
               tol.new <- ng
               break
               }
          else if(nh <= Stop.Tolerance) {
               tol.new <- nh
               break}
          }
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=n, parm.new=modbest[["parm"]],
          parm.old=parm, post=post, Step.Size=nh, tol.new=tol.new)
     return(LA)
     }
.lanm <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     n <- length(as.vector(parm))
     parm.old <- m.old[["parm"]]
     Step.Size <- 1
     if(!is.numeric(parm) || length(parm) < 2)
          stop("Parameters must be a numeric vector of length > 1.")
     # simplex vertices around parm
     V <- t(1/(2*n) * cbind(diag(n), rep(-1, n)) + parm)
     P <- Q <- c()
     # Function values at vertices
     Y <- numeric(n+1)
     for (j in 1:(n+1)) {
          Mo <- Model(V[j, ], Data)
          Y[j] <- Mo[["LP"]] * -1
          V[j,] <- Mo[["parm"]]}
     ho <- lo <- which.min(Y)
     li <- hi <- which.max(Y)
     for (j in 1:(n+1)) {
          if(j != lo && j != hi && Y[j] <= Y[li]) li <- j
          if(j != hi && j != lo && Y[j] >= Y[ho]) ho <- j}
     iter <- 0
     while(Y[hi] > Y[lo] + Stop.Tolerance && iter < Iterations) {
          S <- numeric(n)
          for (j in 1:(n+1)) S <- S + V[j,1:n]
          M <- (S - V[hi,1:n])/n
          R <- 2*M - V[hi,1:n]
          Mo <- Model(R, Data)
          yR <- Mo[["LP"]] * -1
          R <- Mo[["parm"]]
          if(yR < Y[ho]) {
               if(Y[li] < yR) {
                    V[hi,1:n] <- R
                    Y[hi] <- yR
                    }
               else {
                    E <- 2*R - M
                    Mo <- Model(E, Data)
                    yE <- Mo[["LP"]] * -1
                    E <- Mo[["parm"]]
                    if(yE < Y[li]) {
                         V[hi,1:n] <- E
                         Y[hi] <- yE
                         }
                    else {
                         V[hi,1:n] <- R
                         Y[hi] <- yR
                         }
                    }
               }
          else {
               if(yR < Y[hi]) {
                    V[hi,1:n] <- R
                    Y[hi] <- yR}
               C <- (V[hi,1:n] + M) / 2
               Mo <- Model(C, Data)
               yC <- Mo[["LP"]] * -1
               C <- Mo[["parm"]]
               C2 <- (M + R) / 2
               Mo <- Model(C2, Data)
               yC2 <- Mo[["LP"]] * -1
               C2 <- Mo[["parm"]]
               if(yC2 < yC) {
                    C <- C2
                    yC <- yC2}
               if(yC < Y[hi]) {
                    V[hi,1:n] <- C
                    Y[hi] <- yC
                    }
               else {
                    for (j in 1:(n+1)) {
                         if(j != lo) {
                              V[j,1:n] <- (V[j,1:n] + V[lo,1:n]) / 2
                              Z <- V[j,1:n]
                              Mo <- Model(Z, Data)
                              Y[j] <- Mo[["LP"]] * -1
                              Z <- Mo[["parm"]]
                              V[j,1:n] <- Z}}}}
          ho <- lo <- which.min(Y)
          li <- hi <- which.max(Y)
          for (j in 1:(n+1)) {
               if(j != lo && j != hi && Y[j] <= Y[li]) li <- j
               if(j != hi && j != lo && Y[j] >= Y[ho]) ho <- j}
          iter <- iter + 1
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          P <- rbind(P, V[lo, ])
          Q <- c(Q, Y[lo])
          Dev <- rbind(Dev, Model(V[lo,], Data)[["Dev"]])
          }
     snorm <- 0
     for (j in 1:(n+1)) {
          s <- abs(V[j] - V[lo])
          if(s >= snorm) snorm <- s}
     V0 <- V[lo, 1:n]
     y0 <- Y[lo]
     dV <- snorm
     dy <- abs(Y[hi] - Y[lo])
     Dev <- Dev[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=n, parm.new=V0,
          parm.old=parm.old, post=P, Step.Size=Step.Size,
          tol.new=Y[hi] - Y[lo])
     return(LA)
     }
.lanr <- function(Model, parm, Data, Interval, Iterations=100, Stop.Tolerance,
     m.old)
     {
     m.new <- m.old
     Dev <- matrix(m.old[["Dev"]],1,1)
     Step.Size <- 1 / length(parm)
     parm.old <- parm
     parm.len <- length(parm)
     converged <- FALSE
     post <- matrix(parm, 1, parm.len)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          check1 <- TRUE; check2 <- FALSE
          m.old <- Model(parm, Data)
          p <- partial(Model, m.old[["parm"]], Data, Interval)
          H <- Hessian(Model, m.old[["parm"]], Data)
          if(all(is.finite(H))) {
               Sigma <- -as.inverse(H)
               delta <- as.vector(tcrossprod(p, Sigma))
               }
          else check1 <- FALSE
          if(check1 == TRUE) {
               if(any(!is.finite(delta))) check1 <- FALSE
               if(check1 == TRUE) {
                    temp1 <- temp2 <- temp3 <- parm
                    Step.Size1 <- Step.Size / 2
                    Step.Size3 <- Step.Size * 2
                    temp1 <- temp1 + Step.Size1 * delta
                    temp2 <- temp2 + Step.Size * delta
                    temp3 <- temp3 + Step.Size3 * delta
                    Mo1 <- Model(temp1, Data)
                    Mo2 <- Model(temp2, Data)
                    Mo3 <- Model(temp3, Data)
                    check2 <- FALSE
                    if({Mo1[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                         Mo3[["LP"]])} & Mo1[["LP"]] > m.old[["LP"]]) {
                         Step.Size <- Step.Size1
                         parm <- parm + Step.Size * delta
                         m.old <- m.new <- Mo1
                         }
                    else if({Mo2[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                         Mo3[["LP"]])} & Mo2[["LP"]] > m.old[["LP"]]) {
                         parm <- parm + Step.Size * delta
                         m.old <- m.new <- Mo2
                         }
                    else if({Mo3[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                         Mo3[["LP"]])} & Mo3[["LP"]] > m.old[["LP"]]) {
                         Step.Size <- Step.Size3
                         parm <- parm + Step.Size * delta
                         m.old <- m.new <- Mo3
                         }
                    else check2 <- TRUE}}
          if({check1 == FALSE} | {check2 == TRUE}) {
               delta <- p
               temp1 <- temp2 <- temp3 <- parm
               Step.Size1 <- Step.Size / 2
               Step.Size3 <- Step.Size * 2
               temp1 <- temp1 + Step.Size1 * delta
               temp2 <- temp2 + Step.Size * delta
               temp3 <- temp3 + Step.Size3 * delta
               Mo1 <- Model(temp1, Data)
               Mo2 <- Model(temp2, Data)
               Mo3 <- Model(temp3, Data)
               if({Mo1[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                    Mo3[["LP"]])} & Mo1[["LP"]] > m.old[["LP"]]) {
                    Step.Size <- Step.Size1
                    parm <- parm + Step.Size * delta
                    m.old <- m.new <- Mo1
                    }
               else if({Mo2[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                    Mo3[["LP"]])} & Mo2[["LP"]] > m.old[["LP"]]) {
                    parm <- parm + Step.Size * delta
                    f <- Mo2
                    }
               else if({Mo3[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                    Mo3[["LP"]])} & Mo3[["LP"]] > m.old[["LP"]]) {
                    Step.Size <- Step.Size3
                    parm <- parm + Step.Size * delta
                    m.old <- m.new <- Mo3
                    }
               else { #Jitter in case of failure
                    Step.Size <- Step.Size / 2
                    parm <- parm + Step.Size * rnorm(length(parm))}}
          post <- rbind(post, parm)
          Dev <- rbind(Dev, m.new[["Dev"]])
          Step.Size <- max(Step.Size, .Machine$double.eps)
          tol.new <- sqrt(sum(delta^2))
          if(tol.new < Stop.Tolerance) {converged <- TRUE; break}
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     LA <- list(Dev=Dev, H=H, iter=iter, parm.len=parm.len, parm.new=parm,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
.lapso <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     LP <- NA
     parm.len <- length(parm)
     parm.new <- parm.old <- parm
     tol.new <- 1
     post <- matrix(parm.new, 1, parm.len)
     p.s <- floor(10 + 2*sqrt(parm.len)) ## Swarm size
     k <- 3 ### Exponent for calculating number of informants
     p.p <- 1-(1-1/p.s)^k ### Average % of informants
     p.w0 <- p.w1 <- 1 / (2*log(2)) ### Exploitation constants
     p.c.p <- 0.5 + log(2) ### Local exploration constant
     p.c.g <- 0.5 + log(2) ### Global exploration constant
     p.randorder <- TRUE ### Randomize Particle Order
     X <- V <- matrix(parm, parm.len, p.s)
     if(!is.null(Data[["PGF"]])) {
          for (i in 2:ncol(X)) X[,i] <- GIV(Model, Data, PGF=TRUE)
          for (i in 1:ncol(V)) V[,i] <- GIV(Model, Data, PGF=TRUE)
          }
     else {
          for (i in 2:ncol(X)) X[,i] <- GIV(Model, Data, PGF=FALSE)
          for (i in 1:ncol(V)) V[,i] <- GIV(Model, Data, PGF=FALSE)}
     V <- (V - X) / 2 ### Velocity
     mods <- apply(X, 2, function(x) Model(x, Data))
     f.x <- sapply(mods, with, LP) * -1
     f.d <- sapply(mods, with, Dev)
     X <- matrix(sapply(mods, with, parm), nrow(X), ncol(X))
     P <- X
     f.p <- f.x
     P.improved <- rep(FALSE, p.s)
     i.best <- which.min(f.p)
     error <- f.p[i.best]
     init.links <- TRUE
     iter <- 1
     while (iter < Iterations && tol.new > Stop.Tolerance) {
          iter <- iter + 1
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(mod[["LP"]],1), "\n")
          if(p.p != 1 && init.links == TRUE) {
               links <- matrix(runif(p.s*p.s, 0, 1) <= p.p, p.s, p.s)
               diag(links) <- TRUE}
          if(p.randorder == TRUE) index <- sample(p.s)
          else index <- 1:p.s
          for (i in index) {
               if(p.p == 1) j <- i.best
               else j <- which(links[,i])[which.min(f.p[links[,i]])]
               temp <- p.w0 + (p.w1 - p.w0)*(iter / Iterations)
               V[,i] <- temp*V[,i]
               V[,i] <- V[,i] +
                    runif(parm.len, 0, p.c.p)*(P[,i] - X[,i])
               if(i != j)
                    V[,i] <- V[,i] +
                         runif(parm.len, 0, p.c.g)*(P[,j] - X[,i])
               X[,i] <- X[,i] + V[,i]
               mod <- Model(X[,i], Data)
               f.x[i] <- mod[["LP"]] * -1
               f.d[i] <- mod[["Dev"]]
               X[,i] <- mod[["parm"]]
               if(f.x[i] < f.p[i]) { ### Improvement
                    P[,i] <- X[,i]
                    f.p[i] <- f.x[i]
                    if(f.p[i] < f.p[i.best]) i.best <- i}
               }
          post <- rbind(post, P[,i.best])
          Dev <- rbind(Dev, f.d[i.best])
          tol.new <- mean(abs(V[,i.best]))
          init.links <- f.p[i.best] == error
          error <- f.p[i.best]
          }
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=P[,i.best],
          parm.old=parm.old, post=post, Step.Size=mean(abs(V[,i.best])),
          tol.new=tol.new)
     return(LA)
     }
.larprop <- function(Model, parm, Data, Interval, Iterations,
     Stop.Tolerance, m.old)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.len <- length(as.vector(parm))
     approx.grad.old <- approx.grad.new <- rep(0, parm.len)
     parm.old <- m.old[["parm"]]
     parm.new <- parm.old - 0.1 #First step
     names(parm.new) <- names(parm.old) <- Data[["parm.names"]]
     tol.new <- 1
     post <- matrix(parm.new, 1, parm.len)
     Step.Size <- rep(0.0125, parm.len)
     for (iter in 1:Iterations) {
          approx.grad.old <- approx.grad.new
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Approximate Truncated Gradient
          approx.grad.new <- partial(Model, parm.old, Data, Interval)
          approx.grad.new <- interval(approx.grad.new, -1000, 1000,
               reflect=FALSE)
          ### Determine if Gradients Changed Sign
          change <- (approx.grad.old >= 0) == (approx.grad.new >= 0)
          ### Adjust weight (step size) based on sign change
          Step.Size[change == TRUE] <- Step.Size[change == TRUE] * 0.5
          Step.Size[change == FALSE] <- Step.Size[change == FALSE] * 1.2
          Step.Size <- interval(Step.Size, 0.0001, 50, reflect=FALSE)
          ### Propose new state based on weighted approximate gradient
          parm.new <- parm.old + Step.Size * approx.grad.new
          m.new <- Model(parm.new, Data)
          tol.new <- max(sqrt(sum(approx.grad.new^2)),
               sqrt(sum({m.new[["parm"]] - parm.old}^2)))
          if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
               m.new[["Monitor"]]))) | {m.new[["LP"]] < m.old[["LP"]]}) {
               p.order <- sample(1:length(parm.new))
               parm.temp <- parm.old
               for (i in 1:length(p.order)) {
                    parm.temp[p.order[i]] <- parm.new[p.order[i]]
                    m.new <- Model(parm.temp, Data)
                    if(m.new[["LP"]] < m.old[["LP"]])
                         parm.temp[p.order[i]] <- parm.old[p.order[i]]}
               m.new <- Model(parm.temp, Data)}
          parm.new <- m.new[["parm"]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
.lasgd <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     m.new <- m.old
     parm.len <- length(as.vector(parm))
     parm.new <- parm.old <- parm
     names(parm.new) <- names(parm.old) <- Data[["parm.names"]]
     tol.new <- 1
     post <- matrix(parm.new, 1, parm.len)
     #Read SGD specifications
     if(is.null(Data[["file"]]))
          stop("SGD requires Data$file, which is missing.")
     if(is.null(Data[["Nr"]])) stop("SGD requires Data$Nr, which is missing.")
     if(is.null(Data[["Nc"]])) stop("SGD requires Data$Nc, which is missing.")
     if(is.null(Data[["size"]]))
          stop("SGD requires Data$size, which is missing.")
     if(is.null(Data[["epsilon"]]))
          stop("SGD requires Data$epsilon, which is missing.")
     con <- file(Data[["file"]], open="r")
     on.exit(close(con))
     for (iter in 1:Iterations) {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Sample Data
          seek(con, 0)
          skip.rows <- sample.int(Data[["Nr"]] - Data[["size"]], size=1)
          Data[["X"]] <- matrix(scan(file=con, sep=",", skip=skip.rows,
               nlines=Data[["size"]], quiet=TRUE), Data[["size"]],
               Data[["Nc"]], byrow=TRUE)
          ### Propose new parameter values
          g <- partial(Model, m.old[["parm"]], Data)
          parm.new <- m.new[["parm"]] + {Data[["epsilon"]] / 2} * g
          m.new <- Model(parm.new, Data)
          tol.new <- max(sqrt(sum(g^2)),
               sqrt(sum({m.new[["parm"]] - parm.old}^2)))
          if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
               m.new[["Monitor"]]))))
               m.new <- m.old
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Data[["epsilon"]],
          tol.new=tol.new)
     return(LA)
     }
.lasoma <- function(Model, parm, Data, bounds, options=list())
     {
     ### Initial Checks
     if(missing(bounds)) {
          bounds <- list(min=rep(-.Machine$double.xmax,
               length(parm)),
               max=rep(.Machine$double.xmax, length(parm)))
          }
     if(!(all(c("min","max") %in% names(bounds))))
          stop("Bounds list must contain \"min\" and \"max\" vector elements")
     if(length(bounds$min) != length(bounds$max))
          stop("Bounds are not of equal length.")
     ### Option Defaults
     defaultOptions <- list(pathLength=3, stepLength=0.11,
          perturbationChance=0.1, tol=1e-3, maxiter=1000,
          populationSize=10)
     defaultsNeeded <- setdiff(names(defaultOptions), names(options))
     spuriousOptions <- setdiff(names(options), names(defaultOptions))
     options[defaultsNeeded] <- defaultOptions[defaultsNeeded]
     if(length(spuriousOptions) > 0)
           warning("The following specified options are unused: ",
                paste(spuriousOptions,collapse=", "))
     ### Setup
     nParams <- length(bounds$min)
     nParamsTotal <- nParams * options$populationSize
     steps <- seq(0, options$pathLength, options$stepLength)
     nSteps <- length(steps)
     steps <- rep(steps, each=nParamsTotal)
     post <- matrix(parm, 1, nParams)
     ### Create Population
     population <- matrix(parm, nrow=nParams, ncol=options$populationSize)
     for (i in 2:ncol(population)) {
          if(!is.null(Data[["PGF"]]))
               population[,i] <- GIV(Model, Data, PGF=TRUE)
          else population[,i] <- GIV(Model, Data)}
     ### Calculate initial LP and Dev per individual in population
     temp <- apply(population, 2, function(x) Model(x, Data))
     population <- sapply(temp, function(l) l$parm)
     LP <- sapply(temp, function(l) l$LP)
     Dev <- sapply(temp, function(l) l$Dev)
     iteration <- 0
     DevHistory <- LPHistory <- numeric(0)
     if((max(LP) - min(LP)) < options$tol) {
          population <- matrix(runif(nParamsTotal,-10,10), nrow=nParams,
               ncol=options$populationSize)}
     if(all(!is.finite(LP))) stop("All individuals have non-finite LPs.")
     ### Evolution Begins
     repeat {
          ### Find the current leader
          leaderIndex <- which.max(LP)
          leaderDev <- Dev[leaderIndex]
          leaderLP <- LP[leaderIndex]
          ### Check termination criteria
          if(iteration == options$maxiter) break
          LPdiff <- max(LP) - min(LP)
          if(!is.finite(LPdiff)) LPdiff <- 1
          if(LPdiff < options$tol) break
          DevHistory <- c(DevHistory, leaderDev)
          LPHistory <- c(LPHistory, leaderLP)
          ### Find the migration direction for each individual
          directionsFromLeader <- apply(population, 2, "-",
               population[,leaderIndex])
          ### Establish which parameters will be changed
          toPerturb <- runif(nParamsTotal) < options$perturbationChance
          # Second line here has a minus, so directions are away from leader
          populationSteps <- array(rep(population, nSteps),
               dim=c(nParams, options$populationSize, nSteps))
          populationSteps <- populationSteps -
               steps * rep(directionsFromLeader * toPerturb, nSteps)
          ### Replace out-of-bounds parameters with random valid values
          outOfBounds <- which(populationSteps < bounds$min |
               populationSteps > bounds$max)
          randomSteps <- array(runif(nParamsTotal*nSteps),
               dim=c(nParams, options$populationSize, nSteps))
          randomSteps <- randomSteps * (bounds$max-bounds$min) + bounds$min
          populationSteps[outOfBounds] <- randomSteps[outOfBounds]
          ### Values over potential locations
          temp <- apply(populationSteps, 2:3, function(x) Model(x, Data))
          population <- sapply(temp, function(l) l$parm)
          LP <- sapply(temp, function(l) l$LP)
          LP <- matrix(LP, length(LP) / nSteps, nSteps)
          Dev <- sapply(temp, function(l) l$Dev)
          Dev <- matrix(Dev, length(Dev) / nSteps, nSteps)
          individualBestLocs <- apply(LP, 1, which.max)
          ### Migrate each individual to its best new location, and update LP
          indexingMatrix <- cbind(seq_len(options$populationSize),
               individualBestLocs)
          population <- t(apply(populationSteps, 1, "[", indexingMatrix))
          LP <- LP[indexingMatrix]
          Dev <- Dev[indexingMatrix]
          post <- rbind(post, population[,leaderIndex])
          iteration <- iteration + 1
          if(iteration %% round(options$maxiter / 10) == 0)
               cat("Iteration: ", iteration, " of ", options$maxiter, "\n")
          }
    ### Output
    LA <- list(Dev=DevHistory, 
         iter=iteration,
         parm.len=nParams,
         parm.new=population[,leaderIndex],
         parm.old=parm,
         post=post[-1,],
         Step.Size=options$stepLength,
         tol.new=signif(LPdiff,3))
    return(LA)
    }
.laspg <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old)
     {
     parm.old <- parm
     parm.len <- length(parm)
     Dev <- matrix(m.old[["Dev"]],1,1)
     post <- matrix(parm, 1, parm.len)
     M	   <- 10
     ftol     <- 1e-10
     maxfeval <- 10000
     quiet <- FALSE
     feval <- 1
     f0 <- fbest <- f <- m.old[["LP"]] * -1
     nmls <- function(p, f, d, gtd, lastfv, feval, Model, maxfeval, Data)
          {
          # Non-monotone line search of Grippo with safe-guarded quadratic
          # interpolation
          gamma <- 1.e-04
          fmax <- max(lastfv)
          alpha <- 1
          pnew <- p + alpha*d
          m.new <- try(Model(pnew, Data), silent=TRUE)
          feval <- feval + 1
          if(class(m.new) == "try-error" | !is.finite(m.new[["LP"]]))
               return(list(p=NA, f=NA, feval=NA, lsflag=1))
          fnew <- m.new[["LP"]] * -1
          pnew <- m.new[["parm"]]
          while(fnew > fmax + gamma*alpha*gtd) {
    	       if(alpha <= 0.1) alpha <- alpha/2
    	       else {
    	            atemp <- -(gtd*alpha^2) / (2*(fnew - f - alpha*gtd))
    	            if(atemp < 0.1 | atemp > 0.9*alpha) atemp <- alpha/2
    	            alpha <- atemp
    	            }
               pnew <- p + alpha*d
               m.new <- try(Model(pnew, Data), silent=TRUE)
               feval <- feval + 1
               if(class(m.new) == "try-error" | !is.finite(m.new[["LP"]]))
                    return(list(p=NA, f=NA, feval=NA, lsflag=1))
               fnew <- m.new[["LP"]] * -1
               pnew <- m.new[["parm"]]
               if(feval > maxfeval)
                    return(list(p=NA, f=NA, feval=NA, lsflag=2))
               }
          return(list(p=pnew, f=fnew, feval=feval, m.new=m.new, lsflag=0))
          }
     ### Initialization
     lmin <- 1e-30
     lmax <- 1e30
     iter <-  geval <- 0
     lastfv <- rep(-1e99, M)
     fbest <- NA
     fchg <- Inf
     if(any(!is.finite(parm))) stop("Failure in initial guess!")
     pbest <- parm
     g <- partial(Model, parm, Data, Interval) * -1
     geval <- geval + 1
     lastfv[1] <- fbest <- f
     pg <- parm - g
     if(any(is.nan(pg))) stop("Failure in initial projection!")
     pg <- pg - parm
     pg2n <- sqrt(sum(pg*pg))
     pginfn <- max(abs(pg))
     gbest <- pg2n
     if(pginfn != 0) lambda <- min(lmax, max(lmin, 1/pginfn))
     ###  Main iterative loop
     lsflag <- NULL
     while ({iter <= Iterations} &
          ({pginfn > Stop.Tolerance} | {fchg > ftol})) {
          iter <- iter + 1
          d <- parm - lambda * g
          d <- d - parm
          gtd <- sum(g * d)
          if(is.infinite(gtd)) {
               lsflag <- 4
               break}
          nmls.ans <- nmls(parm, f, d, gtd, lastfv, feval , Model,
               maxfeval, Data)
          lsflag <- nmls.ans$lsflag
          if(lsflag != 0) break
          fchg <- abs(f - nmls.ans$f)
          f     <- nmls.ans$f
          pnew  <- nmls.ans$p
          feval <- nmls.ans$feval
          m.new <- nmls.ans$m.new
          lastfv[(iter %% M) + 1] <- f
          gnew <- try(partial(Model, pnew, Data, Interval) * -1,
               silent=TRUE)
          geval <- geval + 1
          if(class(gnew) == "try-error" | any(is.nan(gnew))) {
               lsflag <- 3
               break}
          s <- pnew - parm
          y <- gnew - g
          sts <- sum(s*s)
          yty <- sum(y*y)
          sty <- sum(s*y)
          if(sts == 0 | yty == 0) lambda <- lmax
          else lambda <- min(lmax, max(lmin, sqrt(sts/yty)))
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.new[["LP"]],1), "\n")
          post <- rbind(post, pnew)
          Dev <- rbind(Dev, m.new[["Dev"]])
          parm <- pnew
          g   <- gnew
          pg <- parm - g - parm
          pg2n <- sqrt(sum(pg*pg))
          pginfn <- max(abs(pg))
          f.rep <- f * -1
          if(f < fbest) {
               fbest <- f
               pbest <- pnew
               gbest <- pginfn}
          }
     ### Output
     if(iter < Iterations) pginfn <- max(fchg, pginfn)
     if(is.null(lsflag)) lsflag <- 0
     if(lsflag == 0) parm <- pbest
     else {
          parm <- pbest
          pginfn <- gbest
          }
     m.new <- Model(parm, Data)
     Dev[nrow(Dev),] <- m.new[["Dev"]]
     post[nrow(post),] <- m.new[["parm"]]
     Dev <- Dev[-c(1:2),]; post <- post[-c(1:2),]
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm,
          parm.old=parm.old, post=post, Step.Size=1,
          tol.new=pginfn)
     return(LA)
     }
.lasr1 <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old) {
     m.new <- m.old
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.old <- parm
     parm.len <- length(as.vector(parm))
     post <- matrix(m.old[["parm"]], 1, parm.len)
     tol.new <- 1
     keepgoing <- TRUE
     g.old <- g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
     B <- diag(parm.len) #Approximate Hessian
     options(warn=-1)
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(m.old[["LP"]],1), "\n")
          ### Gradient and Direction p
          g.old <- g.new
          p <- as.vector(tcrossprod(g.old, -as.inverse(B)))
          p[which(!is.finite(p))] <- 0
          ### Step-size Line Search
          Step.Size <- 0.8
          changed <- TRUE
          while(m.new[["LP"]] <= m.old[["LP"]] & changed == TRUE) {
               Step.Size <- Step.Size * 0.2
               s <- Step.Size*p
               prop <- m.old[["parm"]] + s
               changed <- !identical(m.old[["parm"]], prop)
               m.new <- Model(prop, Data)
               if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                    m.new[["Monitor"]]))))
                    m.new <- m.old
               }
          g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
          ### SR1 Update to B, the Secant Approximation of the Hessian
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               g <- g.new - g.old
               if(sum(s*g) > 0) {
                    part1 <- g - B %*% s
                    if(abs(t(s) %*% part1) >=
                         1e-8*sqrt(sum(s*s))*sqrt(sum(part1*part1)))
                         B <- B + (part1 %*% t(part1)) /
                              as.vector(t(part1) %*% s)
                    if(any(!is.finite(B))) B <- diag(parm.len)}
               else B <- diag(parm.len)}
          ### Storage
          post <- rbind(post, m.old[["parm"]])
          Dev <- rbind(Dev, m.old[["Dev"]])
          ### Tolerance
          tol.new <- sqrt(sum(g.new*g.new))
          if(keepgoing == FALSE) tol.new <- 0
          if(tol.new <= Stop.Tolerance) break
          }
     options(warn=0)
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len,
          parm.new=m.old[["parm"]], parm.old=parm.old, post=post,
          Step.Size=Step.Size, tol.new=tol.new)
     return(LA)
     }
.latr <- function(Model, parm, Data, Iterations, m.old)
     {
     fterm <- sqrt(.Machine$double.eps)
     mterm <- sqrt(.Machine$double.eps)
     Norm <- function(x) return(sqrt(sum(x^2)))
     parm.len <- length(parm)
     post <- matrix(parm, 1, parm.len)
     Dev <- matrix(m.old[["Dev"]],1,1)
     r <- rinit <- 1
     rmax <- 1e10
     theta <- parm.old <- parm
     Mo <- m.old
     LP <- Mo[["LP"]]
     theta <- Mo[["parm"]]
     grad <- partial(Model, theta, Data)
     H <- Hessian(Model, theta, Data)
     accept <- TRUE
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% round(Iterations / 10) == 0)
               cat("Iteration: ", iter, " of ", Iterations,
                    ",   LP: ", round(Mo[["LP"]],1), "\n")
          if(accept == TRUE) {
               B <- H
               g <- grad
               f <- LP
               out.value.save <- f
                    B <- - B
                    g <- - g
                    f <- - f
               eout <- eigen(B, symmetric=TRUE)
               gq <- as.numeric(t(eout$vectors) %*% g)}
          is.newton <- FALSE
          if(all(eout$values > 0)) {
               ptry <- as.numeric(- eout$vectors %*% (gq / eout$values))
               if(Norm(ptry) <= r) is.newton <- TRUE}
          if(is.newton == FALSE) {
               lambda.min <- min(eout$values)
               beta <- eout$values - lambda.min
               imin <- beta == 0
               C1 <- sum((gq / beta)[!imin]^2)
               C2 <- sum(gq[imin]^2)
               C3 <- sum(gq^2)
               if(C2 > 0 || C1 > r^2) {
                    is.easy <- TRUE
                    is.hard <- (C2 == 0)
                    beta.dn <- sqrt(C2) / r
                    beta.up <- sqrt(C3) / r
                    beta.adj <- function(x) {
                         if(x == 0) {
                              if(C2 > 0) return(- 1 / r)
                              else return(sqrt(1 / C1) - 1 / r)}
                         return(sqrt(1 / sum((gq /
                              {beta + x})^2)) - 1 / r)}
                    if(beta.adj(beta.up) <= 0) uout <- list(root=beta.up)
                    else if(beta.adj(beta.dn) >= 0) uout <- list(root=beta.dn)
                    else uout <- uniroot(beta.adj, c(beta.dn, beta.up))
                    wtry <- gq / {beta + uout$root}
                    ptry <- as.numeric(-eout$vectors %*% wtry)
                    }
               else {
                    is.hard <- TRUE
                    is.easy <- FALSE
                    wtry <- gq / beta
                    wtry[imin] <- 0
                    ptry <- as.numeric(- eout$vectors %*% wtry)
                    utry <- sqrt(r^2 - sum(ptry^2))
                    if(utry > 0) {
                         vtry <- eout$vectors[ , imin, drop=FALSE]
                         vtry <- vtry[ , 1]
                         ptry <- ptry + utry * vtry}
                    }
               }
          preddiff <- sum(ptry * {g + as.numeric(B %*% ptry) / 2})
          theta.try <- theta + ptry
          Mo <- Model(theta.try, Data)
          LP <- Mo[["LP"]]
          theta.try <- Mo[["parm"]]
          grad <- partial(Model, theta.try, Data)
          H <- Hessian(Model, theta.try, Data)
          ftry <- LP
          ftry <- ftry * -1
          rho <- {ftry - f} / preddiff
          if(ftry < Inf) {
               is.terminate <- {abs(ftry - f) < fterm} || {
                    abs(preddiff) < mterm}
               }
          else {
               is.terminate <- FALSE
               rho <- -Inf}
          if(is.terminate == TRUE) {
               if(ftry < f) {
                    accept <- TRUE
                    theta <- theta.try
                    post <- rbind(post, theta)
                    Dev <- rbind(Dev, Mo[["Dev"]])}
               }
          else {
               if(rho < 1 / 4) {
                    accept <- FALSE
                    r <- r / 4
                    }
               else {
                    accept <- TRUE
                    theta <- theta.try
                    post <- rbind(post, theta)
                    Dev <- rbind(Dev, Mo[["Dev"]])
                    if({rho > 3 / 4} && {is.newton == FALSE})
                    r <- min(2 * r, rmax)}}
          if(is.terminate == TRUE) break
          }
     Mo <- Model(theta, Data)
     theta <- Mo[["parm"]]
     LA <- list(Dev=Dev, H=H, iter=iter, parm.len=parm.len,
          parm.new=theta, parm.old=parm.old, post=post,
          Step.Size=abs(ftry-f), tol.new=abs(preddiff))
     return(LA)
     }
#End
