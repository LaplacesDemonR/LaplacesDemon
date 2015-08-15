###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of the LaplacesDemon function is to use MCMC on the         #
# logarithm of the unnormalized joint posterior density of a Bayesian     #
# model.                                                                  #
###########################################################################

LaplacesDemon <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=10000, Status=100, Thinning=10, Algorithm="MWG",
     Specs=list(B=NULL), Debug=list(DB.chol=FALSE, DB.eigen=FALSE,
     DB.MCSE=FALSE, DB.Model=TRUE), LogFile="", ...)
     {
     cat("\nLaplace's Demon was called on ", date(), "\n", sep="",
          file=LogFile, append=TRUE)
     time1 <- proc.time()
     LDcall <- match.call()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n", file=LogFile, append=TRUE)
     if(missing(Model))
          stop("A function must be entered for Model.", file=LogFile,
                append=TRUE)
     if(!is.function(Model))
          stop("Model must be a function.", file=LogFile, append=TRUE)
     if(missing(Data))
          stop("A list containing data must be entered for Data.",
                file=LogFile, append=TRUE)
     if(is.null(Data[["mon.names"]]))
          stop("In Data, mon.names is NULL.", file=LogFile, append=TRUE)
     if(is.null(Data[["parm.names"]]))
          stop("In Data, parm.names is NULL.", file=LogFile, append=TRUE)
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n", file=LogFile,
                              append=TRUE)}}}}
     if(missing(Initial.Values)) {
          cat("WARNING: Initial Values were not supplied.\n", file=LogFile,
               append=TRUE)
          Initial.Values <- rep(0, length(Data[["parm.names"]]))}
     if(!identical(length(Initial.Values), length(Data[["parm.names"]]))) {
          cat("WARNING: The length of Initial Values differed from",
               "Data$parm.names.\n", file=LogFile, append=TRUE)
          Initial.Values <- rep(0, length(Data[["parm.names"]]))}
     if(any(!is.finite(Initial.Values))) {
          cat("WARNING: Initial Values contain non-finite values.\n",
               file=LogFile, append=TRUE)
          Initial.Values <- rep(0, length(Data[["parm.names"]]))}
     Iterations <- round(abs(Iterations))
     if(Iterations < 11) {
          Iterations <- 11
          cat("'Iterations' has been changed to ", Iterations, ".\n",
               sep="", file=LogFile, append=TRUE)}
     Status <- round(abs(Status))
     if({Status < 1} || {Status > Iterations}) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n",
               sep="", file=LogFile, append=TRUE)}
     Thinning <- round(abs(Thinning))
     if({Thinning < 1} || {Thinning > Iterations}) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n",
               sep="", file=LogFile, append=TRUE)}
     if(Algorithm %in% c("ADMG","AFSS","AGG","AHMC","AIES","AM","AMM",
          "AMWG","CHARM","DEMC","DRAM","DRM","ESS","Experimental","GG",
          "Gibbs","HARM","HMC","HMCDA","IM","INCA","MALA","MCMCMC","MTM",
          "MWG","NUTS","OHSS","pCN","RAM","Refractive","RDMH","RJ","RSS",
          "RWM","SAMWG","SGLD","Slice","SMWG","THMC","twalk","UESS",
          "USAMWG","USMWG")) {
          if(Algorithm == "ADMG") {
               Algorithm <- "Adaptive Directional Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(n=0, Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("n","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["n"]] <- abs(round(Specs[["n"]]))
               Specs[["Periodicity"]] <- max(abs(round(Specs[["Periodicity"]])),
                    length(Initial.Values))
               }
          else if(Algorithm == "AFSS") {
               Algorithm <- "Automated Factor Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(A=Inf, B=NULL, m=Inf, n=0, w=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","B","m","n","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- min(round(abs(Specs[["A"]])), Iterations)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               if(!identical(length(Specs[["m"]]), length(Initial.Values)))
                    Specs[["m"]] <- rep(Specs[["m"]],
                         length(Initial.Values))
               Specs[["n"]] <- abs(round(Specs[["n"]]))
               Specs[["w"]] <- abs(Specs[["w"]])
               if(!identical(length(Specs[["w"]]), length(Initial.Values)))
                    Specs[["w"]] <- rep(Specs[["w"]],
                         length(Initial.Values))
               }
          else if(Algorithm == "AGG") {
               Algorithm <- "Adaptive Griddy-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Grid","dparm","smax","CPUs","Packages","Dyn.libs")
                    %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.list(Specs[["Grid"]])) {
                    if(length(Specs[["Grid"]]) != length(Initial.Values)) {
                         Specs[["Grid"]] <- list(NULL)
                         for (i in 1:length(Initial.Values))
                              Specs[["Grid"]][[i]] <- GaussHermiteQuadRule(3)$nodes
                         cat("\nGrid was misspecified and changed to default.\n",
                              file=LogFile, append=TRUE)}
                    }
               else {
                    temp <- as.vector(Specs[["Grid"]])
                    Specs[["Grid"]] <- list(NULL)
                    for (i in 1:length(Initial.Values))
                         Specs[["Grid"]][[i]] <- temp}
               if(!is.null(Specs[["dparm"]])) {
                    Specs[["dparm"]] <- unique(interval(round(Specs[["dparm"]]), 1,
                         length(Initial.Values)))
                    Specs[["dparm"]] <- Specs[["dparm"]][order(Specs[["dparm"]])]}
               else Specs[["dparm"]] <- 0
               Specs[["smax"]] <- abs(Specs[["smax"]])
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "AHMC") {
               Algorithm <- "Adaptive Hamiltonian Monte Carlo"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(epsilon=rep(1/length(Initial.Values),
                         length(Initial.Values)), L=2, m=rep(1,
                         length(Initial.Values)), Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","L","m","Periodicity") %in% names(Specs)))
               if(!identical(names(Specs),
                    c("epsilon","L","m","Periodicity")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["epsilon"]]))
                    Specs[["epsilon"]] <- rep(1/length(Initial.Values),
                         length(Initial.Values))
               Specs[["epsilon"]] <- as.vector(abs(Specs[["epsilon"]]))
               if(length(Specs[["epsilon"]]) != length(Initial.Values)) {
                    cat("\nLength of epsilon is incorrect.\n",
                         file=LogFile, append=TRUE)
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]][1],
                         length(Initial.Values))}
               Specs[["L"]] <- abs(round(Specs[["L"]]))
               if(Specs[["L"]] < 1) {
                    cat("\nL has been increased to its minimum: 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["L"]] <- 1}
               if(is.null(Specs[["m"]]))
                    Specs[["m"]] <- diag(length(Initial.Values))
               }
          else if(Algorithm == "AIES") {
               Algorithm <- "Affine-Invariant Ensemble Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                          append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Nc","Z","beta","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["Nc"]] <- max(abs(round(Specs[["Nc"]])), 3)
               if(!is.null(Specs[["Z"]])) {
                    if(is.matrix(Specs[["Z"]])) {
                         if(ncol(Specs[["Z"]]) != length(Initial.Values))
                              stop("Z has the wrong number of columns.",
                                   file=LogFile, append=TRUE)
                         if(nrow(Specs[["Z"]]) != Specs[["Nc"]])
                              stop("Z has the wrong number of rows.",
                                   file=LogFile, append=TRUE)}}
               if(length(Specs[["beta"]]) != 1) {
                    cat("\nLength of beta is wrong. Changed to 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["beta"]] <- as.vector(Specs[["beta"]])[1]}
               if(Specs[["beta"]] <= 1) {
                    cat("\nbeta must be > 1. Changed to 2.\n",
                         file=LogFile, append=TRUE)
                         Specs[["beta"]] <- 2}
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               if(Specs[["CPUs"]] > 1 & Specs[["Nc"]] %% 2 != 0)
                    stop("For CPUs > 1, Nc must be even.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "AM") {
               Algorithm <- "Adaptive Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2),
                         Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "AMM") {
               Algorithm <- "Adaptive-Mixture Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2), B=NULL,
                         n=0, Periodicity=1, w=0.05)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","B","n","Periodicity","w") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.null(Specs[["B"]])) {
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Specs[["B"]])) {
                              Covar[[b]] <- diag(length(Specs[["B"]][[b]]))}}}
               Specs[["n"]] <- round(abs(Specs[["n"]]))
               Specs[["w"]] <- abs(Specs[["w"]])
               if(Specs[["w"]] <= 0 || Specs[["w"]] >= 1) {
                    Specs[["w"]] <- 0.05
                    cat("\nw was misspecified and changed to 0.05.\n",
                         file=LogFile, append=TRUE)}
               }
          else if(Algorithm == "AMWG") {
               Algorithm <- "Adaptive Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(B=NULL, n=0, Periodicity=50)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B","n","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "CHARM") {
               Algorithm <- "Componentwise Hit-And-Run Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(alpha.star=NA)
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("alpha.star") %in% names(Specs)))
                         stop("The Specs argument is incorrect",
                              file=LogFile, append=TRUE)
                    Specs[["alpha.star"]] <- abs(as.vector(Specs[["alpha.star"]])[1])
                    if(Specs[["alpha.star"]] <= 0 | Specs[["alpha.star"]] >= 1) {
                         cat("\nalpha.star not in (0,1), set to 0.44.\n",
                              file=LogFile, append=TRUE)
                         alpha.star <- 0.44}}
               }
          else if(Algorithm == "DEMC") {
               Algorithm <- "Differential Evolution Markov Chain"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Nc","Z","gamma","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["Nc"]] <- max(abs(round(Specs[["Nc"]])), 3)
               if(!is.null(Specs[["Z"]])) {
                    if(is.matrix(Specs[["Z"]])) {
                         if(ncol(Specs[["Z"]]) != length(Initial.Values))
                              stop("Z has the wrong number of columns.",
                                   file=LogFile, append=TRUE)
                         if(nrow(Specs[["Z"]]) != (floor(Iterations/Thinning)+1)) {
                              Z.temp <- Specs[["Z"]][nrow(Specs[["Z"]]),]
                              if(nrow(Specs[["Z"]]) < (floor(Iterations/Thinning)+1)) {
                                   Specs[["Z"]] <- rbind(Specs[["Z"]],
                                        Specs[["Z"]][1:(floor(Iterations/Thinning)+1-nrow(Specs[["Z"]])),])
                                   }
                              else if(nrow(Specs[["Z"]]) > (floor(Iterations/Thinning)+1))
                                   Specs[["Z"]] <- Specs[["Z"]][1:(floor(Iterations/Thinning)+1),]
                              Specs[["Z"]][1,] <- Z.temp
                              }
                         Specs[["Z"]] <- array(Specs[["Z"]], dim=c(floor(Iterations/Thinning)+1,
                              length(Initial.Values), Specs[["Nc"]]))}
                    if(dim(Specs[["Z"]])[1] != floor(Iterations/Thinning)+1)
                         stop("The first dimension of Z is incorrect.",
                              file=LogFile, append=TRUE)
                    if(dim(Specs[["Z"]])[2] != length(Initial.Values))
                         stop("The second dimension of Z is incorrect.",
                              file=LogFile, append=TRUE)
                    if(dim(Specs[["Z"]])[3] != Specs[["Nc"]])
                         stop("The third dimension of Z is incorrect.",
                              file=LogFile, append=TRUE)}
               if(is.null(Specs[["gamma"]]))
                    Specs[["gamma"]] <- 2.381204 /
                         sqrt(2*length(Initial.Values))
               else Specs[["gamma"]] <- abs(Specs[["gamma"]])
               Specs[["w"]] <- interval(Specs[["w"]], 0, 1)
               }
          else if(Algorithm == "DRAM") {
               Algorithm <- "Delayed Rejection Adaptive Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2),
                         Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "DRM") {
               Algorithm <- "Delayed Rejection Metropolis"
               Specs <- NULL
               }
          else if(Algorithm == "ESS") {
               Algorithm <- "Elliptical Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(B=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["B"]])) Specs[["B"]] <- list()
               }
          else if(Algorithm == "Experimental") {
               Specs=NULL
               }
          else if(Algorithm == "GG") {
               Algorithm <- "Griddy-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Grid","dparm","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.list(Specs[["Grid"]])) {
                    if(length(Specs[["Grid"]]) != length(Initial.Values)) {
                         Specs[["Grid"]] <- list(NULL)
                         for (i in 1:length(Initial.Values))
                              Specs[["Grid"]][[i]] <- seq(from=-0.1, to=0.1, len=5)
                         cat("\nGrid was misspecified and changed to default.\n",
                              file=LogFile, append=TRUE)}
                    }
               else {
                    temp <- as.vector(Specs[["Grid"]])
                    Specs[["Grid"]] <- list(NULL)
                    for (i in 1:length(Initial.Values))
                         Specs[["Grid"]][[i]] <- temp}
               if(!is.null(Specs[["dparm"]])) {
                    Specs[["dparm"]] <- unique(interval(round(Specs[["dparm"]]), 1,
                         length(Initial.Values)))
                    Specs[["dparm"]] <- Specs[["dparm"]][order(Specs[["dparm"]])]}
               else Specs[["dparm"]] <- 0
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "Gibbs") {
               Algorithm <- "Gibbs Sampler"
               if(missing(Specs) | is.null(Specs)) {
                    cat("\nSpecs missing or null, Algorithm changed to MWG.\n",
                         file=LogFile, append=TRUE)                         
                    Algorithm == "MWG"
                    Specs <- NULL
                    }
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("FC","MWG") %in% names(Specs)))
                         stop("The Specs argument is incorrect",
                              file=LogFile, append=TRUE)
                    if(!is.function(Specs[["FC"]]))
                         stop("FC must be a function.", file=LogFile,
                              append=TRUE)
                    FCtest <- try(Specs[["FC"]](Initial.Values, Data),
                         silent=TRUE)
                    if(inherits(FCtest, "try-error"))
                         stop("Error in FC.", file=LogFile, append=TRUE)
                    if(!is.vector(FCtest))
                         stop("FC must return a vector.", file=LogFile,
                              append=TRUE)
                    if(length(FCtest) != length(Initial.Values))
                         stop("Length of parameters to/from FC differs.",
                              file=LogFile, append=TRUE)
                    if(!is.null(Specs[["MWG"]]) &
                         !is.vector(Specs[["MWG"]]) &
                         !is.numeric(Specs[["MWG"]]))
                         stop("MWG must be a numeric vector.",
                              file=LogFile, append=TRUE)}               
               }
          else if(Algorithm == "HARM") {
               Algorithm <- "Hit-And-Run Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(alpha.star=NA, B=NULL)
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("alpha.star","B") %in% names(Specs)))
                         stop("The Specs argument is incorrect",
                              file=LogFile, append=TRUE)
                    Specs[["alpha.star"]] <- abs(as.vector(Specs[["alpha.star"]])[1])
                    if(Specs[["alpha.star"]] <= 0 | Specs[["alpha.star"]] >= 1) {
                         cat("\nalpha.star not in (0,1), set to 0.234.\n",
                              file=LogFile, append=TRUE)
                         alpha.star <- 0.234}
                    if(is.na(Specs[["alpha.star"]]) & !is.null(Specs[["B"]]))
                         alpha.star <- 0.234}
               if(is.null(Specs[["B"]])) Specs[["B"]] <- list()
               }
          else if(Algorithm == "HMC") {
               Algorithm <- "Hamiltonian Monte Carlo"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(epsilon=rep(1/length(Initial.Values),
                         length(Initial.Values)), L=2, m=rep(1,
                         length(Initial.Values)))
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","L","m") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["epsilon"]] <- abs(Specs[["epsilon"]])
               if(length(Specs[["epsilon"]]) != length(Initial.Values)) {
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]][1],
                         length(Initial.Values))}
               Specs[["L"]] <- abs(round(Specs[["L"]]))
               if(Specs[["L"]] < 1) {
                    cat("\nL has been increased to its minimum: 1.\n",
                         file=LogFile, append=TRUE)
                    L <- 1}
               if(is.null(Specs[["m"]]))
                    Specs[["m"]] <- diag(length(Initial.Values))
               }
          else if(Algorithm == "HMCDA") {
               Algorithm <- "Hamiltonian Monte Carlo with Dual-Averaging"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","delta","epsilon","Lmax","lambda") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- min(round(abs(Specs[["A"]])), Iterations)
               Specs[["delta"]] <- max(min(abs(Specs[["delta"]]), 1),
                    1/Iterations)
               if(!is.null(Specs[["epsilon"]]))
                    Specs[["epsilon"]] <- abs(Specs[["epsilon"]][1])
               Specs[["Lmax"]] <- abs(round(Specs[["Lmax"]]))
               Specs[["lambda"]] <- abs(Specs[["lambda"]])
               if(!is.null(Specs[["epsilon"]]))
                    if(Specs[["lambda"]] < Specs[["epsilon"]])
                         Specs[["lambda"]] <- Specs[["epsilon"]]
               }
          else if(Algorithm == "IM") {
               Algorithm <- "Independence Metropolis"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("mu") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["mu"]] <- as.vector(Specs[["mu"]])
               if(length(Specs[["mu"]]) != length(Initial.Values))
                    stop("length(mu) != length(Initial.Values).",
                         file=LogFile, append=TRUE)
               }
          else if(Algorithm == "INCA") {
               Algorithm <- "Interchain Adaptation"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2),
                         Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "MALA") {
               Algorithm <- "Metropolis-Adjusted Langevin Algorithm"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(A=1e7, alpha.star=0.574, delta=1,
                         epsilon=c(1e-6,1e-7))
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","alpha.star","gamma","delta","epsilon") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- abs(Specs[["A"]][1])
               Specs[["gamma"]] <- min(max(Specs[["gamma"]][1], 0),
                    Iterations)
               Specs[["delta"]] <- min(max(Specs[["delta"]][1], 1e-10),
                    1000)
               Specs[["epsilon"]] <- abs(Specs[["epsilon"]][1:2])
               }
          else if(Algorithm == "MCMCMC") {
               Algorithm <- "Metropolis-Coupled Markov Chain Monte Carlo"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(lambda=1, CPUs=1, Packages=NULL,
                         Dyn.libs=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("lambda","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["lambda"]] <- abs(Specs[["lambda"]])
               if(Specs[["CPUs"]] <= 1)
                    cat("\nCPUs must be at least 2. Attempting 2 CPUs...\n",
                         file=LogFile, append=TRUE)
               Specs[["CPUs"]] <- max(2, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "MTM") {
               Algorithm <- "Multiple-Try Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(K=4, CPUs=1, Packages=NULL, Dyn.libs=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("K","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["K"]] <- abs(round(Specs[["K"]]))
               if(Specs[["CPUs"]] < 1)
                    cat("\nCPUs must be at least 1.\n", file=LogFile,
                         append=TRUE)
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "MWG") {
               Algorithm <- "Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(B=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "NUTS") {
               Algorithm <- "No-U-Turn Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","delta","epsilon","Lmax") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- max(min(round(abs(Specs[["A"]])),
                    Iterations),1)
               Specs[["delta"]] <- max(min(abs(Specs[["delta"]]),
                    1), 1/Iterations)
               if(!is.null(Specs[["epsilon"]]))
                    Specs[["epsilon"]] <- abs(Specs[["epsilon"]][1])
               Specs[["Lmax"]] <- round(abs(Specs[["Lmax"]]))
               }
          else if(Algorithm == "OHSS") {
               Algorithm <- "Oblique Hyperrectangle Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(A=Iterations+1, n=0)
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("A", "n") %in% names(Specs)))
                          stop("The Specs argument is incorrect.",
                               file=LogFile, append=TRUE)
                    Specs[["A"]] <- round(abs(Specs[["A"]]))
                    Specs[["n"]] <- round(abs(Specs[["n"]]))}
               }
          else if(Algorithm == "pCN") {
               Algorithm <- "Preconditioned Crank-Nicolson"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(beta=0.01)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("beta") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["beta"]] <- max(min(Specs[["beta"]], 1), 0)
               }
          else if(Algorithm == "RAM") {
               Algorithm <- "Robust Adaptive Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(alpha.star=0.234, B=NULL, Dist="N",
                         gamma=0.66, n=0)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("alpha.star","B","Dist","gamma","n") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["alpha.star"]] <- Specs[["alpha.star"]][1]
               if(Specs[["alpha.star"]] <= 0 ||
                    Specs[["alpha.star"]] >= 1) {
                    cat("\nalpha.star not in (0,1). Changed to 0.234.\n",
                         file=LogFile, append=TRUE)
                    Specs[["alpha.star"]] <- 0.234}
               if(!is.null(Specs[["B"]])) {
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Specs[["B"]])) {
                              Covar[[b]] <- diag(length(Specs[["B"]][[b]]))}}}
               if(Specs[["Dist"]] != "t" & Specs[["Dist"]] != "N") {
                    cat("\nDist was not t or N, and changed to N.\n",
                         file=LogFile, append=TRUE)
                    Specs[["Dist"]] <- "N"}
               Specs[["gamma"]] <- Specs[["gamma"]][1]
               if(Specs[["gamma"]] <= 0.5 || Specs[["gamma"]] > 1) {
                    cat("\ngamma not in (0.5,1]. Changed to 0.66.\n",
                         file=LogFile, append=TRUE)
                    Specs[["gamma"]] <- 0.66}
               Specs[["n"]] <- abs(Specs[["n"]][1])
               }
          else if(Algorithm == "RDMH") {
               Algorithm <- "Random Dive Metropolis-Hastings"
               Specs <- NULL
               }
          else if(Algorithm == "Refractive") {
               Algorithm <- "Refractive Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=1, m=2, w=0.1, r=1.3)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","m","w","r") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               if(length(Specs[["m"]]) != 1)
                    Specs[["m"]] <- Specs[["m"]][1]
               if(Specs[["m"]] < 2) {
                    cat("\nm was misspecified, and is replaced with 2.\n",
                         file=LogFile, append=TRUE)
                    Specs[["m"]] <- 2}
               Specs[["w"]] <- abs(Specs[["w"]])
               if(length(Specs[["w"]]) != 1)
                    Specs[["w"]] <- Specs[["w"]][1]
               Specs[["r"]] <- abs(Specs[["r"]])
               if(length(Specs[["r"]]) != 1)
                    Specs[["r"]] <- Specs[["r"]][1]
               }
          else if(Algorithm == "RJ") {
               Algorithm <- "Reversible-Jump"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("bin.n","bin.p","parm.p","selectable","selected") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["bin.n"]] <- round(Specs[["bin.n"]])
               if(Specs[["bin.n"]] > length(Initial.Values))
                    Specs[["bin.n"]] <- length(Initial.Values)
               if(Specs[["bin.n"]] < 1) Specs[["bin.n"]] <- 1
               if(Specs[["bin.p"]] < 0 | Specs[["bin.p"]] > 1) {
                    Specs[["bin.p"]] <- interval(Specs[["bin.p"]],
                         0, 1, reflect=FALSE)
                    cat("\nbin.p must be in [0,1]. It's now",
                         round(Specs[["bin.p"]],5), "\n", file=LogFile,
                         append=TRUE)}
               Specs[["parm.p"]] <- as.vector(Specs[["parm.p"]])
               if(length(Specs[["parm.p"]]) != length(Initial.Values)) {
                    Specs[["parm.p"]] <- rep(Specs[["parm.p"]][1],
                         length(Initial.Values))
                    cat("\nparm.p now has the correct length, all equal to parm.p[1].\n",
                         file=LogFile, append=TRUE)}
               Specs[["selectable"]] <- as.vector(Specs[["selectable"]])
               if(length(Specs[["selectable"]]) != length(Initial.Values)) {
                    Specs[["selectable"]] <- rep(1, length(Initial.Values))
                    cat("\nselectable now has the correct length, all set to 1.\n",
                         file=LogFile, append=TRUE)}
               Specs[["selected"]] <- as.vector(Specs[["selected"]])
               if(length(Specs[["selected"]]) != length(Initial.Values)) {
                    Specs[["selected"]] <- rep(1, length(Initial.Values))
                    cat("\nselected now has the correct length, all set to 1.\n",
                         file=LogFile, append=TRUE)}
               }
          else if(Algorithm == "RSS") {
               Algorithm <- "Reflective Slice Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("m","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               if(length(Specs[["m"]]) != 1)
                    Specs[["m"]] <- Specs[["m"]][1]
               if(Specs[["m"]] < 1) {
                    cat("\nm was misspecified, and is replaced with 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["m"]] <- 1}
               Specs[["w"]] <- abs(Specs[["w"]])
               if(length(Specs[["w"]]) != length(Initial.Values))
                    Specs[["w"]] <- rep(Specs[["w"]],
                         len=length(Initial.Values))
               if(any(Specs[["w"]] <= 0)) {
                    cat("\nw was misspecified, and is replaced with 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["w"]][which(Specs[["w"]] <= 0)] <- 1}
               }
          else if(Algorithm == "RWM") {
               Algorithm <- "Random-Walk Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(B=list())
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "SAMWG") {
               Algorithm <- "Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          else if(Algorithm == "SGLD") {
               Algorithm <- "Stochastic Gradient Langevin Dynamics"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","file","Nr","Nc","size") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["Nr"]]  <- abs(round(Specs[["Nr"]]))
               Specs[["Nc"]] <- abs(round(Specs[["Nc"]]))
               Specs[["size"]] <- abs(round(Specs[["size"]]))
               if(Specs[["size"]] >= Specs[["Nr"]])
                    stop("size must be less than nr.")
               if(any(is.na(Specs[["epsilon"]])))
                    Specs[["epsilon"]] <- 1 / Specs[["Nr"]]
               if(length(Specs[["epsilon"]]) == 1)
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]],
                         Iterations)
               if(length(Specs[["epsilon"]]) > Iterations)
                    Specs[["epsilon"]] <- Specs[["epsilon"]][1:Iterations]
               }
          else if(Algorithm == "Slice") {
               Algorithm <- "Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(B=NULL, Bounds=c(-Inf,Inf), m=Inf,
                         Type="Continuous", w=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B","Bounds","m","Type","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["B"]])) {
                    B <- list()
                    B[[1]] <- 1:length(Initial.Values)
                    Bounds <- list()
                    Bounds[[1]] <- Specs[["Bounds"]]
                    m <- list()
                    m[[1]] <- Specs[["m"]]
                    Type <- list()
                    Type[[1]] <- Specs[["Type"]]
                    w <- list()
                    w[[1]] <- Specs[["w"]]
                    Specs[["B"]] <- B
                    Specs[["Bounds"]] <- Bounds
                    Specs[["m"]] <- m
                    Specs[["Type"]] <- Type
                    Specs[["w"]] <- w}
               if(!is.list(Specs[["B"]]))
                    stop("B must be a list.", file=LogFile, append=TRUE)
               if(!is.list(Specs[["Bounds"]])) {
                    Bounds <- list()
                    for (i in 1:length(Initial.Values))
                         Bounds[[i]] <- Specs[["Bounds"]]
                    Specs[["Bounds"]] <- Bounds}
               if(!is.list(Specs[["m"]])) {
                    Specs[["m"]] <- abs(Specs[["m"]][1])
                    m <- list()
                    for (i in 1:length(Initial.Values))
                         m[[i]] <- Specs[["m"]]
                    Specs[["m"]] <- m}
               if(!is.list(Specs[["Type"]])) {
                    Specs[["Type"]] <- Specs[["Type"]][1]
                    if(!Specs[["Type"]] %in% c("Continuous", "Nominal",
                         "Ordinal"))
                         Specs[["Type"]] <- "Continuous"
                    Type <- list()
                    for (i in 1:length(Initial.Values))
                         Type[[i]] <- Specs[["Type"]]
                    Specs[["Type"]] <- Type}
               if(!is.list(Specs[["w"]])) {
                    Specs[["w"]] <- abs(Specs[["w"]][1])
                    w <- list()
                    for (i in 1:length(Initial.Values))
                         w[[i]] <- Specs[["w"]]
                    Specs[["w"]] <- w}
               }
          else if(Algorithm == "SMWG") {
               Algorithm <- "Sequential Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          else if(Algorithm == "THMC") {
               Algorithm <- "Tempered Hamiltonian Monte Carlo"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","L","m", "Temperature") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["epsilon"]] <- as.vector(abs(Specs[["epsilon"]]))
               if(length(Specs[["epsilon"]]) != length(Initial.Values)) {
                    cat("\nLength of epsilon is incorrect.\n", file=LogFile,
                         append=TRUE)
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]][1],
                         length(Initial.Values))}
               Specs[["L"]] <- abs(round(Specs[["L"]]))
               if(Specs[["L"]] < 2) {
                    cat("\nL has been increased to its minimum: 2.\n",
                         file=LogFile, append=TRUE)
                    Specs[["L"]] <- 2}
               if(is.null(Specs[["m"]]))
                    Specs[["m"]] <- diag(length(Initial.Values))
               if(Specs[["Temperature"]] <= 0) {
                    cat("\nTemperature is incorrect, changed to 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["Temperature"]] <- 1}
               }
          else if(Algorithm == "twalk") {
               Algorithm <- "t-walk"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(SIV=NULL, n1=4, at=6, aw=1.5)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("SIV","n1","at","aw") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["SIV"]])) {
                    cat("\nGenerating SIV...\n", file=LogFile, append=TRUE)
                    if(!is.null(Data[["PGF"]]))
                         Specs[["SIV"]] <- GIV(Model, Data, PGF=TRUE)
                    else Specs[["SIV"]] <- GIV(Model, Data)}
               if(!identical(length(Specs[["SIV"]]),
                    length(Initial.Values))) {
                    cat("\nGenerating SIV due to length mismatch.\n",
                         file=LogFile, append=TRUE)
                    if(!is.null(Data[["PGF"]]))
                         Specs[["SIV"]] <- GIV(Model, Data, PGF=TRUE)
                    else Specs[["SIV"]] <- GIV(Model, Data)}
               Mo2 <- Model(Specs[["SIV"]], Data)
               if(!is.finite(Mo2[["LP"]]))
                    stop("SIV results in a non-finite posterior.",
                         file=LogFile, append=TRUE)
               if(!is.finite(Mo2[["Dev"]]))
                    stop("SIV results in a non-finite deviance.",
                         file=LogFile, append=TRUE)
               Specs[["SIV"]] <- Mo2[["parm"]]
               rm(Mo2)
               if(Specs[["n1"]] < 1) {
                    cat("\nn1 must be at least 1. Changed to 4.\n",
                         file=LogFile, append=TRUE)
                    Specs[["n1"]] <- 4}
               if(Specs[["at"]] <= 0) {
                    cat("\nat must be positive. Changed to 6.\n",
                         file=LogFile, append=TRUE)
                    Specs[["at"]] <- 6}
               if(Specs[["aw"]] <= 0) {
                    cat("\naw must be positive. Changed to 1.5.\n",
                         file=LogFile, append=TRUE)
                    Specs[["aw"]] <- 1.5}
               }
          else if(Algorithm == "UESS") {
               Algorithm = "Univariate Eigenvector Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(A=Inf, B=NULL, m=100, n=0)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","B","m","n") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- abs(round(Specs[["A"]]))
               if(!is.null(Specs[["B"]])) {
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Specs[["B"]])) {
                              Covar[[b]] <- diag(length(Specs[["B"]][[b]]))}}}
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               Specs[["n"]] <- abs(round(Specs[["n"]]))
               }
          else if(Algorithm == "USAMWG") {
               Algorithm <- "Updating Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn","Periodicity","Fit","Begin") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          else if(Algorithm == "USMWG") {
               Algorithm <- "Updating Sequential Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn","Fit","Begin") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          }
     else {cat("Unknown algorithm has been changed to Metropolis-within-Gibbs.\n",
                file=LogFile, append=TRUE)
          Algorithm <- "Metropolis-within-Gibbs"
          Specs <- NULL}
     if(!is.null(Specs[["Adaptive"]])) {
          Specs[["Adaptive"]] <- abs(Specs[["Adaptive"]])
          if({Specs[["Adaptive"]] < 1} |
               {Specs[["Adaptive"]] > Iterations})
               Specs[["Adaptive"]] <- Iterations + 1}
     if(!is.null(Specs[["B"]])) {
          if(length(Specs[["B"]]) > 0) {
               if(any(!is.finite(unlist(Specs[["B"]]))))
                    stop("Non-finite values in specification B.",
                         file=LogFile, append=TRUE)
               if(!identical(as.vector(as.numeric(unlist(Specs[["B"]]))),
                    round(abs(as.vector(as.numeric(unlist(Specs[["B"]])))))))
                    stop("Specification B must have only positive integers.",
                         file=LogFile, append=TRUE)
               if(!identical(length(unlist(Specs[["B"]])),
                    length(Initial.Values)))
                    stop("Non-integer values in specification B.",
                         file=LogFile, append=TRUE)}}
     if(!is.null(Specs[["Periodicity"]])) {
          Specs[["Periodicity"]] <- abs(Specs[["Periodicity"]])
          if({Specs[["Periodicity"]] < 1} |
               {Specs[["Periodicity"]] > Iterations})
               Specs[["Periodicity"]] <- Iterations + 1}
     Mo0 <- Model(Initial.Values, Data)
     if(!is.list(Mo0))
          stop("Model must return a list.", file=LogFile, append=TRUE)
     if(length(Mo0) != 5)
          stop("Model must return five components.", file=LogFile,
               append=TRUE)
     if(any(names(Mo0) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.",
               file=LogFile, append=TRUE)
     if(length(Mo0[["LP"]]) > 1)
          stop("Multiple joint posteriors exist!", file=LogFile,
               append=TRUE)
     if(!identical(length(Mo0[["Monitor"]]), length(Data[["mon.names"]])))
          stop("Length of mon.names differs from length of monitors.",
               file=LogFile, append=TRUE)
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
          cat("Suggestion:", acount, "possible instance(s) of apply functions\n",
               file=LogFile, append=TRUE)
          cat("     were found in the Model specification. Iteration speed will\n",
               file=LogFile, append=TRUE)
          cat("     increase if apply functions are vectorized in R or coded\n",
               file=LogFile, append=TRUE)
          cat("     in a faster language such as C++ via the Rcpp package.\n",
               file=LogFile, append=TRUE)}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of for loops\n",
               file=LogFile, append=TRUE)
          cat("     were found in the Model specification. Iteration speed will\n",
               file=LogFile, append=TRUE)
          cat("     increase if for loops are vectorized in R or coded in a\n",
               file=LogFile, append=TRUE)
          cat("     faster language such as C++ via the Rcpp package.\n",
               file=LogFile, append=TRUE)}
     rm(acount)
     if(!identical(Model(Mo0[["parm"]], Data)[["LP"]], Mo0[["LP"]])) {
          cat("WARNING: LP differs when initial values are held constant.\n",
               file=LogFile, append=TRUE)
          cat("     Derivatives may be problematic if used.\n",
               file=LogFile, append=TRUE)}
     #########################  Initial Settings  #########################
     Acceptance <- 0
     if(!is.finite(Mo0[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n",
               file=LogFile, append=TRUE)
          if(!is.null(Data[["PGF"]]))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data)
          Mo0 <- Model(Initial.Values, Data)
          }
     if(is.infinite(Mo0[["LP"]]))
          stop("The posterior is infinite!", file=LogFile, append=TRUE)
     if(is.nan(Mo0[["LP"]]))
          stop("The posterior is not a number!", file=LogFile, append=TRUE)
     if(is.na(Mo0[["Dev"]]))
          stop("The deviance is a missing value!", file=LogFile,
               append=TRUE)
     if(is.infinite(Mo0[["Dev"]]))
          stop("The deviance is infinite!", file=LogFile, append=TRUE)
     if(is.nan(Mo0[["Dev"]]))
          stop("The deviance is not a number!", file=LogFile, append=TRUE)
     if(any(is.na(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have a missing value!",
               file=LogFile, append=TRUE)
     if(any(is.infinite(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have an infinite value!",
               file=LogFile, append=TRUE)
     if(any(is.nan(Mo0[["Monitor"]])))
          stop("Monitored variable(s) include a value that is not a number!",
               file=LogFile, append=TRUE)
     if(Algorithm == "t-walk") {
          Mo0 <- Model(Initial.Values, Data)
          if(any(Mo0[["parm"]] == Specs[["SIV"]]))
              stop("Initial.Values and SIV not unique after model update.",
                   file=LogFile, append=TRUE)}
     ######################  Laplace Approximation  #######################
     ### Sample Size of Data
     if(!is.null(Data[["n"]])) if(length(Data[["n"]]) == 1) N <- Data[["n"]]
     if(!is.null(Data[["N"]])) if(length(Data[["N"]]) == 1) N <- Data[["N"]]
     if(!is.null(Data[["y"]])) N <- nrow(matrix(Data[["y"]]))
     if(!is.null(Data[["Y"]])) N <- nrow(matrix(Data[["Y"]]))
     if(is.null(N))
          stop("Sample size of Data not found in n, N, y, or Y.",
               file=LogFile, append=TRUE)
     if({all(Initial.Values == 0)} & {N >= 5*length(Initial.Values)}) {
          cat("\nLaplace Approximation will be used on initial values.\n",
               file=LogFile, append=TRUE)
          LIV <- length(Initial.Values)
          Fit.LA <- LaplaceApproximation(Model, Initial.Values, Data,
               Method="SPG", CovEst="Identity", sir=FALSE)
          Covar <- 2.381204 * 2.381204 / length(Initial.Values) *
               Fit.LA$Covar
          Initial.Values <- Fit.LA$Summary1[1:length(Initial.Values),1]
          cat("The covariance matrix from Laplace Approximation has been scaled\n",
               file=LogFile, append=TRUE)
          cat("for Laplace's Demon, and the posterior modes are now the initial\n",
               file=LogFile, append=TRUE)
          cat("values for Laplace's Demon.\n\n", file=LogFile, append=TRUE)}
     #########################  Prepare for MCMC  #########################
     Mo0 <- Model(Initial.Values, Data)
     Dev <- matrix(Mo0[["Dev"]], floor(Iterations/Thinning)+1, 1)
     Mon <- matrix(Mo0[["Monitor"]], floor(Iterations/Thinning)+1,
          length(Mo0[["Monitor"]]), byrow=TRUE)
     LIV <- length(Initial.Values)
     thinned <- matrix(Initial.Values, floor(Iterations/Thinning)+1,
          length(Initial.Values), byrow=TRUE)
     ScaleF <- 2.381204 * 2.381204 / LIV
     if(Algorithm %in% c("Adaptive Metropolis",
          "Adaptive-Mixture Metropolis",
          "Delayed Rejection Adaptive Metropolis",
          "Delayed Rejection Metropolis", "Interchain Adaptation",
          "Metropolis-Coupled Markov Chain Monte Carlo",
          "Random-Walk Metropolis")) {
          ### Algorithms that require both VarCov and tuning
          if(is.list(Covar) & Algorithm != "Adaptive-Mixture Metropolis" &
               Algorithm != "Random-Walk Metropolis") {
               Covar <- NULL}
          else if(is.matrix(Covar) & !is.list(Covar)) {
               diag(Covar)[which(diag(Covar) < 1e-100)] <- 1e-100
               tuning <- sqrt(diag(Covar))
               VarCov <- Covar
               }
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != LIV) tuning <- rep(ScaleF, LIV)
               tuning[which(tuning < 1e-100)] <- 1e-100
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning
               }
          else if(is.null(Covar)) {
               tuning <- rep(ScaleF, LIV)
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning
               }
          else if(is.list(Covar)) {
               tuning <- Covar
               for (i in 1:length(tuning)) {
                    tuning[[i]] <- sqrt(diag(tuning[[i]]))}
               VarCov <- Covar}
          if(is.matrix(VarCov) & !is.list(VarCov)) {
               DiagCovar <- matrix(diag(VarCov), 1, LIV)}
          else if(is.list(VarCov)) {
               DiagCovar <- matrix(1, 1, LIV)
               for (b in 1:length(Specs[["B"]])) {
                    DiagCovar[Specs[["B"]][[b]]] <- diag(VarCov[[b]])}}
          }
     else if(Algorithm %in% c("Adaptive Directional Metropolis-within-Gibbs",
          "Automated Factor Slice Sampler",
          "Elliptical Slice Sampler",
          "Independence Metropolis",
          "Metropolis-Adjusted Langevin Algorithm",
          "Oblique Hyperrectangle Slice Sampler",
          "Preconditioned Crank-Nicolson",
          "Robust Adaptive Metropolis",
          "Univariate Eigenvector Slice Sampler")) {
          ### Algorithms that require VarCov, but not tuning
          if(is.list(Covar)) VarCov <- Covar
          else if(is.matrix(Covar) & !is.list(Covar)) VarCov <- Covar
          else if(is.vector(Covar) & !is.list(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- abs(as.vector(Covar))
               diag(VarCov)[which(diag(VarCov) < 1e-100)] <- 1e-100
               }
          else if(is.null(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- rep(ScaleF, LIV)
               }
          else if(is.list(Covar)) VarCov <- Covar
          if(is.matrix(VarCov) & !is.list(VarCov))
               DiagCovar <- matrix(diag(VarCov), 1, LIV)
          else if(is.list(VarCov)) {
               DiagCovar <- matrix(1, 1, LIV)
               for (b in 1:length(Specs[["B"]])) {
                    DiagCovar[Specs[["B"]][[b]]] <- diag(VarCov[[b]])}}
          }
     else if(Algorithm %in% c("Adaptive Griddy-Gibbs",
          "Adaptive Metropolis-within-Gibbs",
          "Gibbs Sampler",
          "Metropolis-within-Gibbs",
          "Multiple-Try Metropolis",
          "Sequential Adaptive Metropolis-within-Gibbs",
          "Sequential Metropolis-within-Gibbs",
          "Updating Sequential Adaptive Metropolis-within-Gibbs",
          "Updating Sequential Metropolis-within-Gibbs")) {
          ### Algorithms that do not require VarCov, but require tuning
          if(is.list(Covar)) Covar <- NULL
          else if(is.matrix(Covar) & !is.list(Covar)) {
               tuning <- sqrt(diag(Covar))
               tuning[which(tuning < 1e-100)] <- 1e-100
               }
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != length(Initial.Values))
                    tuning <- rep(ScaleF, LIV)
               tuning[which(tuning < 1e-100)] <- 1e-100
               }
          else if(is.null(Covar)) {
               tuning <- rep(ScaleF, LIV)}
          VarCov <- NULL
          DiagCovar <- matrix(tuning, 1, LIV)
          }
     else {
          ### Algorithms that do not require VarCov or tuning
          VarCov <- NULL
          DiagCovar <- matrix(1, 1, LIV)
          }
     rm(Covar)
     ############################  Begin MCMC  ############################
     cat("Algorithm:", Algorithm, "\n", file=LogFile, append=TRUE)
     cat("\nLaplace's Demon is beginning to update...\n", file=LogFile,
          append=TRUE)
     options(warn=2)
     on.exit(options(warn=0))
     if(Algorithm == "Adaptive Directional Metropolis-within-Gibbs") {
          mcmc.out <- .mcmcadmg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive Griddy-Gibbs") {
          mcmc.out <- .mcmcagg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "Adaptive Hamiltonian Monte Carlo") {
          mcmc.out <- .mcmcahmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Affine-Invariant Ensemble Sampler") {
          mcmc.out <- .mcmcaies(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Adaptive Metropolis") {
          mcmc.out <- .mcmcam(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & !is.list(VarCov)) {
          mcmc.out <- .mcmcamm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & is.list(VarCov)) {
          mcmc.out <- .mcmcamm.b(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- .mcmcamwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "Automated Factor Slice Sampler") {
          mcmc.out <- .mcmcafss(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Componentwise Hit-And-Run Metropolis") {
          mcmc.out <- .mcmccharm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Delayed Rejection Adaptive Metropolis") {
          mcmc.out <- .mcmcdram(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Delayed Rejection Metropolis") {
          mcmc.out <- .mcmcdrm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Differential Evolution Markov Chain") {
          mcmc.out <- .mcmcdemc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Elliptical Slice Sampler") {
          mcmc.out <- .mcmcess(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Experimental") {
#          mcmc.out <- .mcmcexperimental(Model, Data, Iterations, Status,
#               Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0,
#               ScaleF, thinned, Debug, LogFile)}
          stop("Experimental function not found.", file=LogFile,
               append=TRUE)}
     else if(Algorithm == "Gibbs Sampler") {
          mcmc.out <- .mcmcgibbs(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "Griddy-Gibbs") {
          mcmc.out <- .mcmcgg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Hamiltonian Monte Carlo") {
          mcmc.out <- .mcmchmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging") {
          mcmc.out <- .mcmchmcda(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Hit-And-Run Metropolis") {
          mcmc.out <- .mcmcharm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Independence Metropolis") {
          mcmc.out <- .mcmcim(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Interchain Adaptation") {
          mcmc.out <- .mcmcinca(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Metropolis-Adjusted Langevin Algorithm") {
          mcmc.out <- .mcmcmala(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Metropolis-Coupled Markov Chain Monte Carlo") {
          mcmc.out <- .mcmcmcmcmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Multiple-Try Metropolis") {
          mcmc.out <- .mcmcmtm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned,
               tuning, Debug, LogFile)}
     else if(Algorithm == "Metropolis-within-Gibbs") {
          mcmc.out <- .mcmcmwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "No-U-Turn Sampler") {
          mcmc.out <- .mcmcnuts(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Oblique Hyperrectangle Slice Sampler") {
          mcmc.out <- .mcmcohss(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Preconditioned Crank-Nicolson") {
          mcmc.out <- .mcmcpcn(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Random Dive Metropolis-Hastings") {
          mcmc.out <- .mcmcrdmh(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Random-Walk Metropolis") {
          mcmc.out <- .mcmcrwm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Refractive Sampler") {
          mcmc.out <- .mcmcrefractive(Model, Data, Iterations, Status,
               Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Reflective Slice Sampler") {
          mcmc.out <- .mcmcrss(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned,
               Debug, LogFile)}
     else if(Algorithm == "Reversible-Jump") {
          mcmc.out <- .mcmcrj(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Robust Adaptive Metropolis") {
          mcmc.out <- .mcmcram(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- .mcmcsamwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else if(Algorithm == "Sequential Metropolis-within-Gibbs") {
          mcmc.out <- .mcmcsmwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else if(Algorithm == "Stochastic Gradient Langevin Dynamics") {
          mcmc.out <- .mcmcsgld(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Slice Sampler") {
          mcmc.out <- .mcmcslice(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Tempered Hamiltonian Monte Carlo") {
          mcmc.out <- .mcmcthmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "t-walk") {
          mcmc.out <- .mcmctwalk(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Univariate Eigenvector Slice Sampler") {
          mcmc.out <- .mcmcuess(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Updating Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- .mcmcusamwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else if(Algorithm == "Updating Sequential Metropolis-within-Gibbs") {
          mcmc.out <- .mcmcusmwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else stop("The algorithm is unrecognized.", file=LogFile, append=TRUE)
     options(warn=0)
     #########################  MCMC is Finished  #########################
     Acceptance <- mcmc.out$Acceptance
     Dev <- mcmc.out$Dev
     DiagCovar <- mcmc.out$DiagCovar
     Mon <- mcmc.out$Mon
     thinned <- mcmc.out$thinned
     VarCov <- mcmc.out$VarCov
     remove(mcmc.out)
     rownames(DiagCovar) <- NULL
     colnames(DiagCovar) <- Data[["parm.names"]]
     thinned <- matrix(thinned[-1,], nrow(thinned)-1, ncol(thinned))
     Dev <- matrix(Dev[-1,], nrow(Dev)-1, 1)
     Mon <- matrix(Mon[-1,], nrow(Mon)-1, ncol(Mon))
     if(is.matrix(VarCov) & !is.list(VarCov)) {
          colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]}
     else if(is.vector(VarCov) & !is.list(VarCov)) {
          names(VarCov) <- Data[["parm.names"]]}
     thinned.rows <- nrow(thinned)
     ### Warnings (After Updating)
     if(any(Acceptance == 0))
          cat("\nWARNING: All proposals were rejected.\n", file=LogFile,
               append=TRUE)
     ### Real Values
     thinned[which(!is.finite(thinned))] <- 0
     Dev[which(!is.finite(Dev))] <- 0
     Mon[which(!is.finite(Mon))] <- 0
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n", file=LogFile, append=TRUE)
     if(thinned.rows %% 10 == 0) thinned2 <- thinned
     if(thinned.rows %% 10 != 0) thinned2 <- thinned[1:(10*trunc(thinned.rows/10)),]
     HD <- BMK.Diagnostic(thinned2, batches=10)
     Ind <- 1 * (HD > 0.5)
     BurnIn <- thinned.rows
     batch.list <- seq(from=1, to=nrow(thinned2), by=floor(nrow(thinned2)/10))
     for (i in 1:9) {
          if(sum(Ind[,i:9]) == 0) {
               BurnIn <- batch.list[i] - 1
               break}}
     Stat.at <- BurnIn + 1
     rm(batch.list, HD, Ind, thinned2)
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n", file=LogFile, append=TRUE)
     acf.rows <- trunc(10*log10(thinned.rows))
     acf.temp <- matrix(1, acf.rows, LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV) {
          temp0 <- acf(thinned[,j], lag.max=acf.rows, plot=FALSE)
          if(length(temp0$acf[-1,1,1]) == acf.rows)
               acf.temp[,j] <- abs(temp0$acf[-1,1,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning}
     Rec.Thin[which(is.na(Rec.Thin))] <- nrow(acf.temp)
     ESS3 <- ESS(Mon)
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n", file=LogFile, append=TRUE)
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- sqrt(.colVars(thinned))
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     rm(ESS1)
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025), na.rm=TRUE)
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500), na.rm=TRUE)
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=!Debug[["DB.MCSE"]])
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else {
               if(Debug[["DB.MCSE"]] == TRUE)
                    cat("MCSE of", Data[["parm.names"]][i],
                         "failed in Summary1\n", file=LogFile, append=TRUE)
               Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")}}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)), silent=!Debug[["DB.MCSE"]])
     if(inherits(temp, "try-error")) {
          if(Debug[["DB.MCSE"]] == TRUE)
               cat("MCSE of deviance failed in Summary1\n", file=LogFile,
                    append=TRUE)
          temp <- MCSE(as.vector(Dev), method="sample.variance")}
     Deviance[3] <- temp
     Deviance[4] <- ESS(Dev)
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(as.vector(Mon[,j])), silent=!Debug[["DB.MCSE"]])
          if(inherits(temp, "try-error")) {
               if(Debug[["DB.MCSE"]] == TRUE)
                    cat("MCSE of", Data[["mon.names"]][j],
                         "failed in Summary1\n", file=LogFile, append=TRUE)
               temp <- MCSE(Mon[,j], method="sample.variance")}
          Monitor[3] <- temp
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data[["mon.names"]][j]}
     rm(ESS3)
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(Stat.at < thinned.rows) {
          ESS6 <- ESS(Mon[Stat.at:thinned.rows,])
          thinned2 <- matrix(thinned[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(thinned))
          Dev2 <- matrix(Dev[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Dev))
          Mon2 <- matrix(Mon[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Mon))
          Summ2[,1] <- colMeans(thinned2)
          Summ2[,2] <- sqrt(.colVars(thinned2))
          Summ2[,3] <- 0
          Summ2[,4] <- ESS(thinned[Stat.at:thinned.rows,])
          Summ2[,5] <- apply(thinned2, 2, quantile, c(0.025), na.rm=TRUE)
          Summ2[,6] <- apply(thinned2, 2, quantile, c(0.500), na.rm=TRUE)
          Summ2[,7] <- apply(thinned2, 2, quantile, c(0.975), na.rm=TRUE)
          for (i in 1:ncol(thinned2)) {
               temp <- try(MCSE(thinned2[,i]), silent=!Debug[["DB.MCSE"]])
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else {
                    if(Debug[["DB.MCSE"]] == TRUE)
                         cat("MCSE of", Data[["parm.names"]][i],
                              "failed in Summary2\n", file=LogFile,
                              append=TRUE)
                    Summ2[i,3] <- MCSE(thinned2[,i],
                         method="sample.variance")}}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- try(MCSE(as.vector(Dev2)), silent=!Debug[["DB.MCSE"]])
          if(inherits(temp, "try-error")) {
               if(Debug[["DB.MCSE"]] == TRUE)
                    cat("MCSE of deviance failed in Summary2\n",
                         file=LogFile, append=TRUE)
               temp <- MCSE(as.vector(Dev2), method="sample.variance")}
          Deviance[3] <- temp
          Deviance[4] <- ESS(Dev[Stat.at:thinned.rows,])
          Deviance[5] <- as.numeric(quantile(Dev2, probs=0.025,
               na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev2, probs=0.500,
               na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev2, probs=0.975,
               na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon2[,j])
               Monitor[2] <- sd(as.vector(Mon2[,j]))
               temp <- try(MCSE(as.vector(Mon2[,j])),
                    silent=!Debug[["DB.MCSE"]])
               if(inherits(temp, "try-error")) {
                    if(Debug[["DB.MCSE"]] == TRUE)
                         cat("MCSE of", Data[["mon.names"]][j],
                              "failed in Summary2\n", file=LogFile,
                              append=TRUE)
                    temp <- MCSE(Mon2[,j], method="sample.variance")}
               Monitor[3] <- temp
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon2[,j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon2[,j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon2[,j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data[["mon.names"]][j]}
          rm(ESS6)
          }
     ### Column names to samples
     if(identical(ncol(Mon), length(Data[["mon.names"]])))
          colnames(Mon) <- Data[["mon.names"]]
     if(identical(ncol(thinned), length(Data[["parm.names"]]))) {
          colnames(thinned) <- Data[["parm.names"]]}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     if(Algorithm %in% c("Adaptive Griddy-Gibbs",
          "Affine-Invariant Ensemble Sampler",
          "Automated Factor Slice Sampler",
          "Componentwise Hit-And-Run Metropolis",
          "Delayed Rejection Metropolis",
          "Elliptical Slice Sampler",
          "Gibbs Sampler",
          "Griddy-Gibbs",
          "Hamiltonian Monte Carlo",
          "Hit-And-Run Metropolis",
          "Independence Metropolis",
          "Metropolis-Adjusted Langevin Algorithm",
          "Metropolis-Coupled Markov Chain Monte Carlo",
          "Metropolis-within-Gibbs",
          "Multiple-Try Metropolis",
          "No-U-Turn Sampler",
          "Oblique Hyperrectangle Slice Sampler",
          "Preconditioned Crank-Nicolson",
          "Random Dive Metropolis-Hastings",
          "Random-Walk Metropolis",
          "Reflective Slice Sampler",
          "Refractive Sampler",
          "Reversible-Jump",
          "Sequential Metropolis-within-Gibbs",
          "Slice Sampler",
          "Stochastic Gradient Langevin Dynamics",
          "Tempered Hamiltonian Monte Carlo",
          "t-walk",
          "Univariate Eigenvector Slice Sampler") &
          {Stat.at < thinned.rows}) {
          cat("Estimating Log of the Marginal Likelihood\n", file=LogFile,
               append=TRUE)
          LML <- LML(theta=thinned2, LL=as.vector(Dev2)*(-1/2),
               method="NSIS")}
     time2 <- proc.time()
     ### Compile Output
     cat("Creating Output\n", file=LogFile, append=TRUE)
     LaplacesDemon.out <- list(Acceptance.Rate=round(Acceptance/Iterations,7),
          Algorithm=Algorithm,
          Call=LDcall,
          Covar=VarCov,
          CovarDHis=DiagCovar,
          Deviance=as.vector(Dev),
          DIC1=c(mean(as.vector(Dev)),
               var(as.vector(Dev))/2,
               mean(as.vector(Dev)) + var(as.vector(Dev))/2),
          DIC2=if(Stat.at < thinned.rows) {
               c(mean(as.vector(Dev2)),
               var(as.vector(Dev2))/2,
               mean(as.vector(Dev2)) +
               var(as.vector(Dev2))/2)}
               else rep(NA,3),
          Initial.Values=Initial.Values,
          Iterations=Iterations,
          LML=LML[[1]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60,2),
          Model=Model,
          Monitor=Mon,
          Parameters=LIV,
          Posterior1=thinned,
          Posterior2=if(Stat.at < thinned.rows) {
               thinned[Stat.at:thinned.rows,]}
               else thinned[thinned.rows,],
          Rec.BurnIn.Thinned=BurnIn,
          Rec.BurnIn.UnThinned=BurnIn*Thinning,
          Rec.Thinning=min(1000, max(Rec.Thin)),
          Specs=Specs,
          Status=Status,
          Summary1=Summ1,
          Summary2=Summ2,
          Thinned.Samples=thinned.rows,
          Thinning=Thinning)
     class(LaplacesDemon.out) <- "demonoid"
     cat("\nLaplace's Demon has finished.\n", file=LogFile, append=TRUE)
     return(LaplacesDemon.out)
     }
.mcmcadmg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     Debug, LogFile)
     {
     n <- Specs[["n"]]
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     AccRate <- rep(0, LIV)
     obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
     obs.scatter <- tcrossprod(Mo0[["parm"]])*n
     s <- svd(VarCov)
     U <- diag(s$u)
     tol <- LIV*max(s$d)*.Machine$double.eps
     problem <- any(s$d <= tol)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning)+1, LIV,
          byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Random-Scan Componentwise Estimation
          if(iter > 10) AccRate <- Acceptance / {iter - 1}
          if(problem == FALSE) lambda <- U*rnorm(LIV, 0,
               sqrt(0.01 + s$d*exp(2*s$d*(AccRate - 0.3))))
          else lambda <- rnorm(LIV, 0, sqrt(diag(VarCov)))
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- prop[j] + lambda[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed for",
                              Data[["parm.names"]][j], ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal for",
                              Data[["parm.names"]][j],
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         Acceptance[j] <- Acceptance[j] + 1}}}
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + Mo0[["parm"]]
          obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
          ### Adaptation
          if(iter %% Periodicity == 0) {
               VarCov <- obs.scatter/{n + iter} -
                    tcrossprod(obs.sum/{n + iter})
               diag(VarCov) <- diag(VarCov) + 1e-05
               s <- svd(VarCov)
               U <- diag(s$u)
               tol <- LIV*max(s$d)*.Machine$double.eps
               problem <- any(s$d <= tol)}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]
               DiagCovar[t.iter,] <- diag(VarCov)}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcafss <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     Debug, LogFile)
     {
     A <- Specs[["A"]]
     Block <- Specs[["B"]]
     m <- Specs[["m"]]
     n <- Specs[["n"]]
     w <- Specs[["w"]]
     B <- length(Block)
     targetRatio <- 0.5
     if(B == 0) {
          if(!is.symmetric.matrix(VarCov)) {
               cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
                    append=TRUE)
               VarCov <- as.symmetric.matrix(VarCov)}
          if(!is.positive.definite(VarCov)) {
               cat("\nNon-Positive-Definite Covar, correcting now...\n",
                    file=LogFile, append=TRUE)
               VarCov <- as.positive.definite(VarCov)}
          decomp.freq <- max(LIV * floor(Iterations / Thinning / 100), 10)
          cat("\nEigendecomposition will occur every", decomp.freq,
               "iterations.\n\n", file=LogFile, append=TRUE)
          factors <- eigen(VarCov)$vectors
          obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
          obs.scatter <- tcrossprod(Mo0[["parm"]])*n
          DiagCovar <- matrix(w, floor(Iterations/Thinning)+1, LIV,
               byrow=TRUE)
          nExpands <- nShrinks <- rep(0, LIV)
          IterPerAdapt <- 1
          nProposals <- 0
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Random-Scan Componentwise Estimation
               for (j in sample.int(LIV)) {
                    y.slice <- Mo0[["LP"]] - rexp(1)
                    upper <- runif(1,0,w[j])
                    lower <- upper - w[j]
                    ### Step Out
                    count <- 0
                    while (count <= m[j]) {
                         Mo1 <- try(Model(Mo0[["parm"]] +
                              lower*factors[,j], Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the lower",
                                        "bound failed for",
                                        Data[["parm.names"]][j],
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                              lower <- lower + w[j]
                              break}
                         else if(!is.finite(Mo1[["LP"]])) {
                              if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the lower",
                                        "bound for", Data[["parm.names"]][j],
                                        "resulted in a non-finite LP",
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                              lower <- lower + w[j]
                              break}
                         nExpands[j] <- nExpands[j] + 1
                         if(Mo1[["LP"]] <= y.slice) break
                         lower <- lower - w[j]
                         count <- count + 1
                         }
                    count <- 0
                    while (count <= m[j]) {
                         Mo1 <- try(Model(Mo0[["parm"]] +
                              upper*factors[,j], Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the upper",
                                        "bound failed for",
                                        Data[["parm.names"]][j],
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                              upper <- upper - w[j]
                              break}
                         else if(!is.finite(Mo1[["LP"]])) {
                              if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the upper",
                                        "bound for", Data[["parm.names"]][j],
                                        "resulted in a non-finite LP",
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                              upper <- upper - w[j]
                              break}
                         nExpands[j] <- nExpands[j] + 1
                         if(Mo1[["LP"]] <= y.slice) break
                         upper <- upper + w[j]
                         count <- count + 1
                         }
                    ### Rejection Sampling
                    repeat {
                         lower <- -abs(min(lower, upper))
                         upper <- abs(max(lower, upper))
                         prop <- runif(1, lower, upper)
                         Mo1 <- try(Model(Mo0[["parm"]] + prop *
                              factors[,j], Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Rejection sampling",
                                        "failed for",
                                        Data[["parm.names"]][j], "\n",
                                        file=LogFile, append=TRUE)
                              Mo1 <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                              Mo1[["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Rejection sampling for",
                                        Data[["parm.names"]][j],
                                        "resulted in non-finite value(s).\n",
                                        file=LogFile, append=TRUE)
                              Mo1 <- Mo0}
                         if(Mo1[["LP"]] >= y.slice) break
                         else if(abs(prop) < 1e-100) break
                         nShrinks[j] <- nShrinks[j] + 1
                         if(prop < 0) lower <- prop
                         else upper <- prop
                         }
                    Mo0 <- Mo1
                    }
               nProposals <- nProposals + 1
               obs.sum <- obs.sum + Mo0[["parm"]]
               obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
               ### Adaptation
               if({iter <= A} & {A - iter >= decomp.freq}) {
                    ### Tune Interval Widths
                    if(nProposals %% IterPerAdapt == 0) {
                         denom <- nExpands + nShrinks
                         for (j in 1:LIV) {
                              if(denom[j] > 0) {
                                   ratio <- nExpands[j] / denom[j]
                                   if(ratio == 0) ratio <- 1 / denom[j]
                                   multiplier <- ratio / targetRatio
                                   w[j] <- w[j]*multiplier
                                   }
                              }
                         nExpands <- nShrinks <- rep(0,LIV)
                         nProposals <- 0
                         IterPerAdapt <- IterPerAdapt * 2}
                    ### Tune Sampling Factors
                    if(iter %% decomp.freq == 0) {
                         VarCov <- obs.scatter/{n + iter} -
                              tcrossprod(obs.sum/{n + iter})
                         factors <- eigen(VarCov)$vectors
                         nExpands <- nShrinks <- rep(0,LIV)
                         IterPerAdapt <- 1
                         nProposals <- 0}
                    }
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]
                    DiagCovar[t.iter,] <- w}
               }
          }
     else {
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from number ",
                    "of blocks.", file=LogFile, append=TRUE)
          factors <- obs.sum <- obs.scatter <- list()
          decomp.freq <- rep(0, length(B))
          for (b in 1:B) {
               if(length(Block[[b]]) == 1)
                    stop("Single-parameter blocks are not allowed in AFSS.",
                         file=LogFile, append=TRUE)
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from block length.")
               if(!is.symmetric.matrix(VarCov[[b]])) {
                    cat("\nAsymmetric Covar block, correcting now...\n",
                         file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.symmetric.matrix(VarCov[[b]])}
               if(!is.positive.definite(VarCov[[b]])) {
                    cat("\nNon-Positive-Definite Covar block,",
                         "correcting now...\n", file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.positive.definite(VarCov[[b]])}
               decomp.freq[b] <- max(length(Block[[b]]) *
                    floor(Iterations / Thinning / 100), 10)
               factors[[b]] <-try(eigen(VarCov[[b]])$vectors,
                    silent=!Debug[["DB.eigen"]])
               if(inherits(factors[[b]], "try-error")) {
                    if(Debug[["DB.eigen"]] == TRUE)
                         cat("\nWARNING: Eigendecomposition of covariance",
                              "matrix failed for block", b, ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Eigendecomposition of an identity matrix",
                              "occurs instead.\n", file=LogFile, append=TRUE)
                    factors[[b]] <- diag(length(Block[[b]]))}
               obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n,
                    length(Block[[b]]), 1)
               obs.scatter[[b]] <- tcrossprod(Mo0[["parm"]][Block[[b]]])*n}
          if(all(decomp.freq == decomp.freq[1]))
               cat("\nEigendecomposition will occur every", decomp.freq[1],
                    "iterations.\n\n", file=LogFile, append=TRUE)
          else cat("\nEigendecomposition frequency varies by block,",
                    "and will occur between\n",
                    min(decomp.freq), "and", max(decomp.freq),
                    "iterations.\n\n", file=LogFile, append=TRUE)
          DiagCovar <- matrix(w, floor(Iterations/Thinning)+1, LIV,
               byrow=TRUE)
          nExpands <- nShrinks <- rep(0, LIV)
          IterPerAdapt <- rep(1, B)
          nProposals <- rep(0, B)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Blockwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Proceed by Block
               for (b in 1:B) {
                    ### Random-Scan Componentwise Estimation
                    for (j in sample(Block[[b]])) {
                         bj <- which(Block[[b]] == j)
                         y.slice <- Mo0[["LP"]] - rexp(1)
                         upper <- runif(1,0,w[j])
                         lower <- upper - w[j]
                         ### Step Out
                         count <- 0
                         while (count <= m[j]) {
                              prop <- Mo0[["parm"]]
                              prop[Block[[b]]] <- prop[Block[[b]]] +
                                   lower*factors[[b]][,bj]
                              Mo1 <- try(Model(prop, Data),
                                   silent=!Debug[["DB.Model"]])
                              if(inherits(Mo1, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the lower",
                                        "bound failed for",
                                        Data[["parm.names"]][j],
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                                   lower <- lower + w[j]
                                   break}
                              else if(!is.finite(Mo1[["LP"]])) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the lower",
                                        "bound for", Data[["parm.names"]][j],
                                        "resulted in a non-finite LP",
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                                   lower <- lower + w[j]
                                   break}
                              nExpands[j] <- nExpands[j] + 1
                              if(Mo1[["LP"]] <= y.slice) break
                              lower <- lower - w[j]
                              count <- count + 1
                              }
                         count <- 0
                         while (count <= m[j]) {
                              prop <- Mo0[["parm"]]
                              prop[Block[[b]]] <- prop[Block[[b]]] +
                                   upper*factors[[b]][,bj]
                              Mo1 <- try(Model(prop, Data),
                                   silent=!Debug[["DB.Model"]])
                              if(inherits(Mo1, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the upper",
                                        "bound failed for",
                                        Data[["parm.names"]][j],
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                                   upper <- upper - w[j]
                                   break}
                              else if(!is.finite(Mo1[["LP"]])) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                   cat("\nWARNING: Stepping out the upper",
                                        "bound for", Data[["parm.names"]][j],
                                        "resulted in a non-finite LP",
                                        "in step", count+1, ".\n",
                                        file=LogFile, append=TRUE)
                                   upper <- upper - w[j]
                                   break}
                              nExpands[j] <- nExpands[j] + 1
                              if(Mo1[["LP"]] <= y.slice) break
                              upper <- upper + w[j]
                              count <- count + 1
                              }
                         ### Rejection Sampling
                         repeat {
                              prop <- Mo0[["parm"]]
                              lower <- -abs(min(lower, upper))
                              upper <- abs(max(lower, upper))
                              u <- runif(1, lower, upper)
                              prop[Block[[b]]] <- prop[Block[[b]]] +
                                   u*factors[[b]][,bj]
                              Mo1 <- try(Model(prop, Data),
                                   silent=!Debug[["DB.Model"]])
                              if(inherits(Mo1, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                        cat("\nWARNING: Rejection sampling",
                                             "failed for",
                                             Data[["parm.names"]][j], "\n",
                                             file=LogFile, append=TRUE)
                                   Mo1 <- Mo0
                                   }
                              else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                   Mo1[["Monitor"]])))) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                        cat("\nWARNING: Rejection sampling for",
                                             Data[["parm.names"]][j],
                                             "resulted in non-finite value(s).\n",
                                             file=LogFile, append=TRUE)
                                   Mo1 <- Mo0}
                              if(Mo1[["LP"]] >= y.slice) break
                              else if(abs(u) < 1e-100) break
                              nShrinks[j] <- nShrinks[j] + 1
                              if(u < 0) lower <- u
                              else upper <- u
                              }
                         Mo0 <- Mo1
                         }
                    nProposals[b] <- nProposals[b] + 1
                    obs.sum[[b]] <- obs.sum[[b]] + Mo0[["parm"]][Block[[b]]]
                    obs.scatter[[b]] <- obs.scatter[[b]] +
                         tcrossprod(Mo0[["parm"]][Block[[b]]])
                    ### Adaptation
                    if({iter <= A} & {A - iter >= decomp.freq[b]}) {
                         ### Tune Interval Widths
                         if(nProposals[b] %% IterPerAdapt[b] == 0) {
                              for (j in Block[[b]]) {
                                   denom <- nExpands[j] + nShrinks[j]
                                   if(denom > 0) {
                                        ratio <- nExpands[j] / denom
                                        if(ratio == 0) ratio <- 1 / denom
                                        multiplier <- ratio / targetRatio
                                        w[j] <- w[j]*multiplier
                                        }
                                   }
                              nExpands[Block[[b]]] <- rep(0,length(Block[[b]]))
                              nShrinks[Block[[b]]] <- rep(0,length(Block[[b]]))
                              nProposals[b] <- 0
                              IterPerAdapt[b] <- IterPerAdapt[b] * 2}
                         ### Tune Sampling Factors
                         if(iter %% decomp.freq[b] == 0) {
                              VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
                                   tcrossprod(obs.sum[[b]]/{n + iter})
                              factors[[b]] <- eigen(VarCov[[b]])$vectors
                              nExpands[Block[[b]]] <- rep(0,length(Block[[b]]))
                              nShrinks[Block[[b]]] <- rep(0,length(Block[[b]]))
                              IterPerAdapt[b] <- 1
                              nProposals[b] <- 0}
                         }
                    ### Save Thinned Samples
                    if(iter %% Thinning == 0) {
                         t.iter <- floor(iter / Thinning) + 1
                         thinned[t.iter,] <- Mo0[["parm"]]
                         Dev[t.iter] <- Mo0[["Dev"]]
                         Mon[t.iter,] <- Mo0[["Monitor"]]
                         DiagCovar[t.iter,] <- w}
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcagg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     Debug, LogFile)
     {
     Grid <- Specs[["Grid"]]
     dparm <- Specs[["dparm"]]
     smax <- Specs[["smax"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     AGGCP <- function(Model, Data, j, Mo0, Grid, tuning, smax, Debug,
          LogFile)
          {
          G <- length(Grid[[j]])
          x <- Grid[[j]] * sqrt(2) * tuning[j]
          LP.grid <- rep(0, G)
          prop <- Mo0[["parm"]]
          theta <- prop[j] + x
          for (g in 1:G) {
               prop[j] <- theta[g]
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE)
                         cat("\nWARNING: Evaluating",
                              Data[["parm.names"]][j], "at",
                              round(prop[j],5), "failed.\n", file=LogFile,
                              append=TRUE)
                    Mo1 <- Mo0}
               LP.grid[g] <- Mo1[["LP"]]
               theta[g] <- Mo1[["parm"]][j]}
          if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
          LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
          LP.grid <- exp(LP.grid - logadd(LP.grid))
          LP.grid <- LP.grid / sum(LP.grid)
          s <- spline(theta, LP.grid, n=1000)
          s$y <- interval(s$y, 0, Inf, reflect=FALSE)
          if(length(which(s$y > 0)) == 0)
               prop[j] <- theta[which.max(LP.grid)[1]]
          else prop[j] <- sample(s$x, 1, prob=s$y)
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE)
                    cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                         "at", round(prop[j],5), "failed.\n",
                         file=LogFile, append=TRUE)
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE)
                    cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                         "at", round(prop[j],5),
                         "resulted in non-finite value(s).\n",
                         file=LogFile, append=TRUE)
               Mo1 <- Mo0
               }
          else tuning[j] <- min(max(sqrt(sum(LP.grid * x^2)),
                    1e-10), smax)
          Mo0 <- Mo1
          return(list(Mo0=Mo0, tuning=tuning))
          }
     AGGCPP <- function(Model, Data, j, Mo0, Grid, tuning, smax, Debug,
          LogFile, cl)
          {
          G <- length(Grid[[j]])
          x <- Grid[[j]] * sqrt(2) * tuning[j]
          LP.grid <- rep(0, G)
          LIV <- length(Mo0[["parm"]])
          prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
          prop[, j] <- prop[, j] + x
          Mo1 <- parLapply(cl, 1:G,
               function(x) Model(prop[x,], Data))
          LP.grid <- as.vector(unlist(lapply(Mo1,
               function(x) x[["LP"]])))
          prop <- matrix(as.vector(unlist(lapply(Mo1,
               function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
          theta <- prop[, j]
          if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
          LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
          LP.grid <- exp(LP.grid - logadd(LP.grid))
          LP.grid <- LP.grid / sum(LP.grid)
          s <- spline(theta, LP.grid, n=1000)
          s$y <- interval(s$y, 0, Inf, reflect=FALSE)
          prop <- Mo0[["parm"]]
          if(length(which(s$y > 0)) == 0)
               prop[j] <- theta[which.max(LP.grid)[1]]
          else prop[j] <- sample(s$x, 1, prob=s$y)
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE)
                    cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                         "at", round(prop[j],5), "failed.\n",
                         file=LogFile, append=TRUE)
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE)
                    cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                         "at", round(prop[j],5),
                         "resulted in non-finite value(s).\n",
                         file=LogFile, append=TRUE)
               Mo1 <- Mo0
               }
          else tuning[j] <- min(max(sqrt(sum(LP.grid * x^2)),
                    1e-10), smax)
          Mo0 <- Mo1
          return(list(Mo0=Mo0, tuning=tuning))
          }
     Acceptance <- matrix(0, 1, LIV)
     Grid.orig <- Grid
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     if(CPUs == 1) {
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- .mcmcggdp(Model, Data, j, Mo0, Grid, Debug,
                              LogFile)
                    else {
                         agg <- AGGCP(Model, Data, j, Mo0, Grid, tuning,
                              smax, Debug, LogFile)
                         Mo0 <- agg$Mo0
                         tuning[j] <- agg$tuning[j]}}
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]
                    DiagCovar <- rbind(DiagCovar, tuning)}}
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
               append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n",
                    file=LogFile, append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...",
               file=LogFile, append=TRUE)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
               ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- .mcmcggdpp(Model, Data, j, Mo0, Grid,
                              Debug, LogFile, cl)
                    else {
                         agg <- AGGCPP(Model, Data, j, Mo0, Grid,
                              tuning, smax, Debug, LogFile, cl)
                         Mo0 <- agg$Mo0
                         tuning[j] <- agg$tuning[j]}
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]
                    DiagCovar <- rbind(DiagCovar, tuning)}}}
     DiagCovar <- DiagCovar[-1,]
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar,
          Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcahmc <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     m <- Specs[["m"]]
     invm <- as.inverse(m)
     U <- chol(m)
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     DiagCovar <- matrix(epsilon, floor(Iterations/Periodicity), LIV,
          byrow=TRUE)
     gr0 <- partial(Model, post[1,], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Propose new parameter values
          prop <- post[iter,] <- Mo0[["parm"]]
          momentum0 <- as.vector(rnorm(LIV) %*% U)
          kinetic0 <- t(momentum0) %*% invm %*% momentum0 / 2
          momentum1 <- momentum0 + (epsilon / 2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + as.vector(epsilon %*% invm) * momentum1
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed in leapfrog", l,
                              ".\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal in leapfrog", l,
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1}
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed in leapfrog",
                                   l, ".\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),
                                   ")",sep=""), "\n", file=LogFile,
                                   append=TRUE)}
                         Mo1 <- Mo0.1
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in leapfrog",
                                   l, "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),
                                   ")", sep=""), "\n", file=LogFile,
                                   append=TRUE)}
                         Mo1 <- Mo0.1}}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon / 2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- t(momentum1) %*% invm %*% momentum1 / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               }
          ### Adaptation
          if(iter %% Periodicity == 0) {
               if(iter > 10) {
                    acceptances <- length(unique(post[(iter-9):iter,1]))
                    if(acceptances <= 1) epsilon <- epsilon * 0.8
                    else if(acceptances > 7) epsilon <- epsilon * 1.2}
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- epsilon}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=cov(thinned))
     return(out)
     }
.mcmcaies <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     Nc <- Specs[["Nc"]]
     Z <- Specs[["Z"]]
     beta <- Specs[["beta"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     Mo0 <- list(Mo0=Mo0)
     if(is.null(Z)) {
          Z <- matrix(Mo0[[1]][["parm"]], Nc, LIV, byrow=TRUE)
          for (i in 2:Nc) {
               if(!is.null(Data[["PGF"]])) {
                    Z[i,] <- GIV(Model, Data, PGF=TRUE)
                    }
               else Z[i,] <- GIV(Model, Data)
               }
          }
     for (i in 2:Nc) Mo0[[i]] <- Model(Z[i,], Data)
     if(CPUs == 1) {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile,
                         append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[[1]][["parm"]]
                    Dev[t.iter] <- Mo0[[1]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
               for (i in 1:Nc) {
                    ### Propose new parameter values with stretch move
                    z <- 1 / sqrt(runif(1, 1 / beta, beta))
                    s <- sample(c(1:Nc)[-i], 1)
                    prop <- Mo0[[s]][["parm"]] +
                         z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
                    if(i == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Multivariate,   LP: ",
                              round(Mo0[[1]][["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed in walker",
                                   i, ".\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),
                                   ")",sep=""), "\n", file=LogFile,
                                   append=TRUE)}
                         Mo1 <- Mo0[[i]]
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in walker", i,
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),
                                   ")",sep=""), "\n", file=LogFile,
                                   append=TRUE)}
                         Mo1 <- Mo0[[i]]}
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    else if(log.u < log.alpha) {
                         Mo0[[i]] <- Mo1
                         if(i == 1) {
                              Acceptance <- Acceptance + 1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[["parm"]]
                                   Dev[t.iter] <- Mo1[["Dev"]]
                                   Mon[t.iter,] <- Mo1[["Monitor"]]}
                              }
                         }
                    }
               }
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
               append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n", file=LogFile,
                    append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...",
               file=LogFile, append=TRUE)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
          ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())
          model.wrapper <- function(x, ...)
               {
               if(!is.null(Packages)) {
                    sapply(Packages,
                         function(x) library(x, character.only=TRUE,
                              quietly=TRUE))}
               if(!is.null(Dyn.libs)) {
                    sapply(Dyn.libs,
                         function(x) dyn.load(paste(wd, x, sep = "/")))
                    on.exit(sapply(Dyn.libs,
                         function(x) dyn.unload(paste(wd, x, sep = "/"))))}
               Model(prop[x,], Data)
               }
          prop <- Z
          batch1 <- 1:(Nc/2)
          batch2 <- batch1 + (Nc/2)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile,
                         append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[[1]][["parm"]]
                    Dev[t.iter] <- Mo0[[1]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
               for (i in 1:Nc) {
                    ### Propose new parameter values with stretch move
                    z <- 1 / sqrt(runif(1, 1 / beta, beta))
                    if(i <= (Nc/2)) s <- sample(batch2, 1)
                    else s <- sample(batch1, 1)
                    prop[i,] <- Mo0[[s]][["parm"]] +
                         z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
                    if(i == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Multivariate\n", file=LogFile,
                              append=TRUE)}
               ### Log-Posterior of the proposed state
               Mo1 <- clusterApply(cl, 1:Nc, model.wrapper,
                    Model, Data, prop)
               for (i in 1:Nc) {
                    if(any(!is.finite(c(Mo1[[i]][["LP"]],
                         Mo1[[i]][["Dev"]], Mo1[[i]][["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in walker", i,
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[i,],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1[[i]] <- Mo0[[i]]}
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[[i]][["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    else if(log.u < log.alpha) {
                         Mo0[[i]] <- Mo1[[i]]
                         if(i == 1) {
                              Acceptance <- Acceptance + 1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[[i]][["parm"]]
                                   Dev[t.iter] <- Mo1[[i]][["Dev"]]
                                   Mon[t.iter,] <- Mo1[[i]][["Monitor"]]}
                              }
                         }
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcam <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV,
          byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov),
               silent=!Debug[["DB.chol"]])
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               if(Debug[["DB.chol"]] == TRUE)
                    cat("\nWARNING: Cholesky decomposition failed for",
                         "proposal in iteration", iter, ".\n",
                         file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0}
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          else if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * cov(post[1:iter,])} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               ### Univariate Standard Deviations
               tuning <- sqrt(ScaleF * .colVars(post[1:iter,]) +
                    ScaleF * 1.0E-5)
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcamm <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Block <- Specs[["B"]]
     n <- Specs[["n"]]
     Periodicity <- Specs[["Periodicity"]]
     w <- Specs[["w"]]
     obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
     obs.scatter <- tcrossprod(Mo0[["parm"]])*n
     if(all(upper.triangle(VarCov) == 0)) prop.R <- NULL
     else prop.R <- ScaleF * chol(VarCov)
     tuning <- sqrt(0.0001 * ScaleF)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV,
          byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Propose new parameter values from a mixture
          if(is.null(prop.R) || runif(1) < w) {
               prop <- rnorm(LIV, Mo0[["parm"]], tuning)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Non-Adaptive Component,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)}
          else {
               prop <- Mo0[["parm"]] +
                    as.vector(rbind(rnorm(LIV)) %*% prop.R)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Adaptive Component,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else {
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1}}
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + Mo0[["parm"]]
          obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               VarCov <- obs.scatter/{n + iter} -
                    tcrossprod(obs.sum/{n + iter})
               diag(VarCov) <- diag(VarCov) + 1e-05
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               prop.R <- try(ScaleF * chol(VarCov),
                    silent=!Debug[["DB.chol"]])
               if(Debug[["DB.chol"]] == TRUE)
                    cat("\nWARNING: Cholesky decomposition failed for",
                         "proposal covariance in iteration", iter, ".\n",
                         file=LogFile, append=TRUE)
               if(!is.matrix(prop.R)) prop.R <- NULL}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcamm.b <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Block <- Specs[["B"]]
     n <- Specs[["n"]]
     Periodicity <- Specs[["Periodicity"]]
     w <- Specs[["w"]]
     B <- length(Block)
     if(!identical(length(VarCov), B))
          stop("Number of components in Covar differs from ",
               "number of blocks.", file=LogFile, append=TRUE)
     obs.scatter <- obs.sum <- prop.R <- list()
     for (b in 1:B) {
          if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
               stop("Diagonal of Covar[[",b,"]] differs from ",
                    "block length.", file=LogFile, append=TRUE)
          obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n,
               length(Block[[b]]), 1)
          obs.scatter[[b]] <- matrix(tcrossprod(Mo0[["parm"]][Block[[b]]])*n,
               length(Block[[b]]), length(Block[[b]]))
          if(all(upper.triangle(VarCov[[b]]) == 0)) prop.R[[b]] <- NA
          else prop.R[[b]] <- ScaleF * chol(VarCov[[b]])}
     tuning <- sqrt(0.0001 * ScaleF)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Proceed by Block
          for (b in 1:B) {
               ### Propose new parameter values from a mixture
               prop <- Mo0[["parm"]]
               if(any(is.na(prop.R[[b]])) || runif(1) < w) {
                    prop[Block[[b]]] <- rnorm(length(Block[[b]]),
                         Mo0[["parm"]][Block[[b]]], tuning)
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Blockwise,   LP: ",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)}
               else {
                    prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]] +
                         as.vector(rbind(rnorm(length(Block[[b]]))) %*%
                              prop.R[[b]])
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Blockwise,   LP: ",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)}
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop[Block[[b]]],
                              collapse=","),")",sep=""), "\n",
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal resulted in non-finite",
                              "value(s).\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop[Block[[b]]],
                              collapse=","),")",sep=""), "\n",
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance +
                              length(Block[[b]]) / LIV}}
               ### Update Sample and Scatter Sum
               obs.sum[[b]] <- obs.sum[[b]] + Mo0[["parm"]][Block[[b]]]
               obs.scatter[[b]] <- obs.scatter[[b]] +
                    tcrossprod(Mo0[["parm"]][Block[[b]]])
               ### Adapt the Proposal Variance
               if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
                    VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
                         tcrossprod(obs.sum[[b]]/{n + iter})
                    diag(VarCov[[b]]) <- diag(VarCov[[b]]) + 1e-05
                    if(b == 1) DiagCovar <- rbind(DiagCovar, rep(0,LIV))
                    DiagCovar[nrow(DiagCovar),Block[[b]]] <- diag(VarCov[[b]])
                    prop.R[[b]] <- try(ScaleF * chol(VarCov[[b]]),
                         silent=!Debug[["DB.chol"]])
                    if(Debug[["DB.chol"]] == TRUE)
                         cat("\nWARNING: Cholesky decomposition failed for",
                              "proposal covariance in iteration", iter,
                              ".\n", file=LogFile, append=TRUE)
                    if(!is.matrix(prop.R[[b]])) prop.R[[b]] <- NA}
               }
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcamwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     Debug, LogFile)
     {
     Block <- Specs[["B"]]
     n <- Specs[["n"]]
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- rep(0, LIV)
     B <- length(Block)
     DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV,
          byrow=TRUE)
     if(B == 0) {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Random-Scan Componentwise Estimation
               propdraw <- rnorm(LIV,0,tuning)
               for (j in sample.int(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + propdraw[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed for",
                                   Data[["parm.names"]][j], ".\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal for",
                                   Data[["parm.names"]][j],
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         u <- log(runif(1)) < {Mo1[["LP"]] - Mo0[["LP"]]}
                         if(u == TRUE) {
                              Mo0 <- Mo1
                              Acceptance[j] <- Acceptance[j] + 1}}}
               ### Adapt the Proposal Variance
               if(iter %% Periodicity == 0) {
                    size <- 1 / min(100, sqrt(n + iter))
                    Acceptance.Rate <- Acceptance / iter
                    log.tuning <- log(tuning)
                    tuning.num <- which(Acceptance.Rate > 0.44)
                    log.tuning[tuning.num] <- log.tuning[tuning.num] + size
                    tuning.num <- which(Acceptance.Rate <= 0.44)
                    log.tuning[tuning.num] <- log.tuning[tuning.num] - size
                    tuning <- exp(log.tuning)
                    a.iter <- floor(iter / Periodicity)
                    DiagCovar[a.iter,] <- tuning}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     else {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Random-Scan Componentwise Estimation
               propdraw <- rnorm(LIV,0,tuning)
               for (b in 1:B) {
                    for (j in sample(Block[[b]])) {
                         ### Propose new parameter values
                         prop <- Mo0[["parm"]]
                         prop[j] <- prop[j] + propdraw[j]
                         ### Log-Posterior of the proposed state
                         Mo1 <- try(Model(prop, Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal failed for",
                                        Data[["parm.names"]][j], ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop[j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                              Mo1[["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal for",
                                        Data[["parm.names"]][j],
                                        "resulted in non-finite value(s).\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop[j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0
                              }
                         else {
                              ### Accept/Reject
                              u <- log(runif(1)) < {Mo1[["LP"]] -
                                   Mo0[["LP"]]}
                              if(u == TRUE) {
                                   Mo0 <- Mo1
                                   Acceptance[j] <- Acceptance[j] + 1}}}}
               ### Adapt the Proposal Variance
               if(iter %% Periodicity == 0) {
                    size <- 1 / min(100, sqrt(n + iter))
                    Acceptance.Rate <- Acceptance / iter
                    log.tuning <- log(tuning)
                    tuning.num <- which(Acceptance.Rate > 0.44)
                    log.tuning[tuning.num] <- log.tuning[tuning.num] + size
                    tuning.num <- which(Acceptance.Rate <= 0.44)
                    log.tuning[tuning.num] <- log.tuning[tuning.num] - size
                    tuning <- exp(log.tuning)
                    a.iter <- floor(iter / Periodicity)
                    DiagCovar[a.iter,] <- tuning}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=mean(Acceptance),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
.mcmccharm <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     if(is.na(alpha.star)) {
          Acceptance <- matrix(0, 1, LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample.int(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed for",
                                   Data[["parm.names"]][j], ".\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal for",
                                   Data[["parm.names"]][j],
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                         if(u == TRUE) {
                              Mo0 <- Mo1
                              Acceptance[j] <- Acceptance[j] + 1}}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          ### Output
          out <- list(Acceptance=mean(as.vector(Acceptance)),
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=.colVars(thinned))
          return(out)
          }
     else {
          tau <- rep(1, LIV)
          Acceptance <- matrix(0, 1, LIV)
          DiagCovar <- matrix(tau, nrow(thinned), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample.int(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + tau[j]*lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed for",
                                   Data[["parm.names"]][j], ".\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal for",
                                   Data[["parm.names"]][j],
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                         if(u == TRUE) {
                              Mo0 <- Mo1
                              Acceptance[j] <- Acceptance[j] + 1
                              tau[j] <- tau[j] + (tau[j] / (alpha.star *
                              (1 - alpha.star))) * (1 - alpha.star) / iter
                              }
                         else {
                              tau[j] <- abs(tau[j] - (tau[j] / (alpha.star *
                              (1 - alpha.star))) * alpha.star / iter)}}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]
                    DiagCovar[t.iter,] <- tau}
               }
          ### Output
          out <- list(Acceptance=mean(as.vector(Acceptance)),
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=.colVars(thinned))
          return(out)
          }
     }
.mcmcdemc <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     Nc <- Specs[["Nc"]]
     Z <- Specs[["Z"]]
     gamma <- Specs[["gamma"]]
     w <- Specs[["w"]]
     const <- 2.381204 / sqrt(2)
     Mo0 <- list(Mo0=Mo0)
     if(is.null(Z)) {
          cat("\nGenerating Z...\n", file=LogFile, append=TRUE)
          Z <- array(0, dim=c(floor(Iterations/Thinning)+1, LIV, Nc))
          for (t in 1:dim(Z)[1]) {
               for (i in 1:Nc) {
                    if(t == 1 & i == 1) {
                         Z[t,,i] <- Mo0[[1]][["parm"]]
                         }
                    else {
                         if(!is.null(Data[["PGF"]])) {
                              Z[t,,i] <- GIV(Model, Data, PGF=TRUE)}
                         else Z[t,,i] <- GIV(Model, Data)
                         }
                    }
               }
          }
     else Z[1,,1] <- Mo0[[1]][["parm"]]
     for (i in 2:Nc) Mo0[[i]] <- Model(Z[1,,i], Data)
     for (iter in 1:Iterations) {
          ### Thinned Iteration
          t.iter <- floor(iter / Thinning) + 1
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               Z[t.iter,,] <- Z[t.iter-1,,]
               thinned[t.iter,] <- Mo0[[1]][["parm"]]
               Dev[t.iter] <- Mo0[[1]][["Dev"]]
               Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
          omega <- runif(1)
          for (i in 1:Nc) {
               r <- sample(dim(Z)[1], 2)
               s <- sample(c(1:Nc)[-i], 2)
               if(omega > w) {
                    ### Parallel Direction Move
                    prop <- Mo0[[i]][["parm"]] +
                         gamma*(Z[r[1],,s[1]] - Z[r[2],,s[2]]) +
                         runif(LIV, -0.001, 0.001)^LIV
                         }
               else {
                    ### Snooker Move
                    si <- sample(c(1:Nc)[-i], 1)
                    prop <- Mo0[[i]][["parm"]] + const*
                         ({Mo0[[si]][["parm"]] - Z[r[1],,s[1]]} -
                          {Mo0[[si]][["parm"]] - Z[r[2],,s[2]]})}
               if(i == 1 & iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP: ",
                         round(Mo0[[1]][["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed in chain", i,
                              ".\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0[[i]]
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal in chain", i,
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0[[i]]
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0[[i]] <- Mo1
                         Z[t.iter,,i] <- Mo1[["parm"]]
                         if(i == 1) {
                              Acceptance <- Acceptance + 1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[["parm"]]
                                   Dev[t.iter] <- Mo1[["Dev"]]
                                   Mon[t.iter,] <- Mo1[["Monitor"]]}
                              }
                         }
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=thinned,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcdram <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     DR <- 1
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV,
          byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov),
               silent=!Debug[["DB.chol"]])
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               if(Debug[["DB.chol"]] == TRUE)
                    cat("\nWARNING: Cholesky decomposition failed for",
                         "proposal 1 in iteration", iter, ".\n",
                         file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal 1 failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal 1 resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0}
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVNz <- try(rbind(rnorm(LIV)) %*%
                    chol(VarCov * 0.5), silent=!Debug[["DB.chol"]])
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    if(Debug[["DB.chol"]] == TRUE)
                         cat("\nWARNING: Cholesky decomposition failed for",
                              "proposal 2 in iteration", iter, ".\n",
                              file=LogFile, append=TRUE)
                    prop <- post[iter,]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo2, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal 2 failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo2 <- Mo0
                    }
               else if(any(!is.finite(c(Mo2[["LP"]], Mo2[["Dev"]],
                    Mo2[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal 2 resulted in non-finite",
                              "value(s).\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo2 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    options(warn=-1)
                    log.alpha.comp <- log(1 - exp(Mo1[["LP"]] -
                         Mo2[["LP"]]))
                    options(warn=0)
                    if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
                    log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                         {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] -
                         Mo0[["LP"]]))}
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo2
                         post[iter,] <- Mo2[["parm"]]
                         Acceptance <- Acceptance + 1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}
                         }
                    }
               }
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * cov(post[1:iter,])} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               ### Univariate Standard Deviations
               tuning <- sqrt(ScaleF * .colVars(post[1:iter,]) +
                    ScaleF * 1.0E-5)
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcdrm <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     DR <- 1
     U <- chol(VarCov)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% U, silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               prop <- Mo0[["parm"]]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal 1 failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal 1 resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0}
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1}
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVNz <- try(rbind(rnorm(LIV)) %*%
                    chol(VarCov * 0.5), silent=!Debug[["DB.chol"]])
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
               else {
                    if(Debug[["DB.chol"]] == TRUE)
                         cat("\nWARNING: Cholesky decomposition failed for",
                              "proposal 2 in iteration", iter, ".\n",
                              file=LogFile, append=TRUE)
                    prop <- Mo0[["parm"]]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo2, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal 2 failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo2 <- Mo0
                    }
               else if(any(!is.finite(c(Mo2[["LP"]], Mo2[["Dev"]],
                    Mo2[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal 2 resulted in non-finite",
                              "value(s).\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""),"\n", file=LogFile, append=TRUE)}
                    Mo2 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    options(warn=-1)
                    log.alpha.comp <- log(1 - exp(Mo1[["LP"]] -
                         Mo2[["LP"]]))
                    options(warn=0)
                    if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
                    log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                         {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] -
                         Mo0[["LP"]]))}
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo2
                         Acceptance <- Acceptance + 1}}
               }
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcess <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     Debug, LogFile)
     {
     Block <- Specs[["B"]]
     if(length(Block) == 0) {
          if(!is.symmetric.matrix(VarCov)) {
               cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
                    append=TRUE)
               VarCov <- as.symmetric.matrix(VarCov)}
          if(!is.positive.definite(VarCov)) {
               cat("\nNon-Positive-Definite Covar, correcting now...\n",
                    file=LogFile, append=TRUE)
               VarCov <- as.positive.definite(VarCov)}
          nu <- rnorm(LIV, 0, diag(VarCov))
          U <- chol(VarCov)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Propose new parameter values
               nu <- as.vector(rbind(rnorm(LIV)) %*% U)
               theta <- theta.max <- runif(1, 0, 2*pi)
               theta.min <- theta - 2*pi
               shrink <- TRUE
               log.u <- log(runif(1))
               ### Rejection Sampling
               while (shrink == TRUE) {
                    prop <- Mo0[["parm"]] * cos(theta) + nu*sin(theta)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Rejection sampling failed.\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Rejection sampling resulted",
                                   "in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0}
                    ### Accept/Reject
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         shrink <- FALSE
                         }
                    else {
                         if(theta < 0) theta.min <- theta
                         else theta.max <- theta
                         theta <- runif(1, theta.min, theta.max)}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     else {
          B <- length(Block)
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from ",
                    "number of blocks.", file=LogFile, append=TRUE)
          nu <- rep(NA, LIV)
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from ",
                         "block length.", file=LogFile, append=TRUE)
               if(!is.symmetric.matrix(VarCov[[b]])) {
                    cat("\nAsymmetric Covar block, correcting now...\n",
                         file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.symmetric.matrix(VarCov[[b]])}
               if(!is.positive.definite(VarCov[[b]])) {
                    cat("\nNon-Positive-Definite Covar block,",
                         "correcting now...\n", file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.positive.definite(VarCov[[b]])}
               nu[Block[[b]]] <- rnorm(length(Block[[b]]), 0,
                    diag(VarCov[[b]]))}
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Blockwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    blen <- length(Block[[b]])
                    nu[Block[[b]]] <- as.vector(rbind(rnorm(blen)) %*%
                         chol(VarCov[[b]]))
                    theta <- theta.max <- runif(1, 0, 2*pi)
                    theta.min <- theta - 2*pi
                    shrink <- TRUE
                    log.u <- log(runif(1))
                    ### Rejection Sampling
                    while (shrink == TRUE) {
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]] *
                              cos(theta) + nu[Block[[b]]]*sin(theta)
                         ### Log-Posterior of the proposed state
                         Mo1 <- try(Model(prop, Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Rejection sampling failed.\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                              Mo1[["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Rejection sampling resulted",
                                        "in non-finite value(s).\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0}
                         ### Accept/Reject
                         log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                         if(!is.finite(log.alpha)) log.alpha <- 0
                         if(log.u < log.alpha) {
                              Mo0 <- Mo1
                              shrink <- FALSE
                              }
                         else {
                              if(theta < 0) theta.min <- theta
                              else theta.max <- theta
                              theta <- runif(1, theta.min, theta.max)}}
                    }
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=cov(thinned))
     return(out)
     }
.mcmcgibbs <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     Debug, LogFile)
     {
     FC <- Specs[["FC"]]
     MWG <- Specs[["MWG"]]
     if(is.null(MWG)) {
          Acceptance <- Iterations
          MWGlen <- 0}
     else {
          MWGlen <- length(MWG)
          Acceptance <- matrix(0, 1, LIV)}
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Gibbs Sampling of Full Conditionals
          prop <- try(FC(Mo0[["parm"]], Data), silent=!Debug[["DB.Model"]])
          if(inherits(prop, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Gibbs proposal for full conditionals",
                         "failed.\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(Mo0[["parm"]], collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               prop <- Mo0[["parm"]]}
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Gibbs proposal failed.\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",sep=""),
                         "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0}
          Mo0 <- Mo1
          ### Metropolis-within-Gibbs
          if(MWGlen > 0) {
               ### Random-Scan Componentwise Estimation
               for (j in sample(MWG)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- rnorm(1, prop[j], tuning[j])
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: MWG proposal failed for",
                                   Data[["parm.names"]][j], ".\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: MWG proposal for",
                                   Data[["parm.names"]][j],
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                         if(u == TRUE) {
                              Mo0 <- Mo1
                              Acceptance[j] <- Acceptance[j] + 1}}}
               }
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     if(MWGlen > 0) Acceptance <- mean(as.vector(Acceptance[,MWG]))
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
.mcmcgg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     Grid <- Specs[["Grid"]]
     dparm <- Specs[["dparm"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     if(CPUs == 1) {
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- .mcmcggdp(Model, Data, j, Mo0, Grid, Debug,
                              LogFile)
                    else Mo0 <- .mcmcggcp(Model, Data, j, Mo0, Grid, Debug,
                              LogFile)
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]}}
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
               append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n",
                    file=LogFile, append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...",
               file=LogFile, append=TRUE)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
               ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- .mcmcggdpp(Model, Data, j, Mo0, Grid,
                              Debug, LogFile, cl)
                    else Mo0 <- .mcmcggcpp(Model, Data, j, Mo0, Grid,
                              Debug, LogFile, cl)
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]}}}
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar,
          Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
### Griddy-Gibbs Continuous Parameter (Non-Parallelized)
.mcmcggcp <- function(Model, Data, j, Mo0, Grid, Debug, LogFile)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     prop <- Mo0[["parm"]]
     theta <- prop[j] + Grid[[j]]
     for (g in 1:G) {
          prop[j] <- theta[g]
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE)
                    cat("\nWARNING: Evaluating",
                         Data[["parm.names"]][j], "at",
                         round(prop[j],5), "failed.\n", file=LogFile,
                         append=TRUE)
               Mo1 <- Mo0}
          LP.grid[g] <- Mo1[["LP"]]
          theta[g] <- Mo1[["parm"]][j]}
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     s <- spline(theta, LP.grid, n=1000)
     s$y <- interval(s$y, 0, Inf, reflect=FALSE)
     if(length(which(s$y > 0)) == 0)
          prop[j] <- theta[which.max(LP.grid)[1]]
     else prop[j] <- sample(s$x, 1, prob=s$y)
     Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
     if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5), "failed.\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0
          }
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5),
                    "resulted in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0}
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Continuous Parameter (Parallelized)
.mcmcggcpp <- function(Model, Data, j, Mo0, Grid, Debug, LogFile, cl)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     LIV <- length(Mo0[["parm"]])
     prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
     prop[, j] <- prop[, j] + Grid[[j]]
     Mo1 <- parLapply(cl, 1:G,
          function(x) Model(prop[x,], Data))
     LP.grid <- as.vector(unlist(lapply(Mo1,
          function(x) x[["LP"]])))
     prop <- matrix(as.vector(unlist(lapply(Mo1,
          function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
     theta <- prop[, j]
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     s <- spline(theta, LP.grid, n=1000)
     s$y <- interval(s$y, 0, Inf, reflect=FALSE)
     prop <- Mo0[["parm"]]
     if(length(which(s$y > 0)) == 0)
          prop[j] <- theta[which.max(LP.grid)[1]]
     else prop[j] <- sample(s$x, 1, prob=s$y)
     Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
     if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5), "failed.\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0
          }
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5),
                    "resulted in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0}
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Discrete Parameter (Non-Parallelized)
#where j is which parameter, and Grid are discrete values
.mcmcggdp <- function(Model, Data, j, Mo0, Grid, Debug, LogFile)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     prop <- Mo0[["parm"]]
     theta <- Grid[[j]]
     for (g in 1:G) {
          prop[j] <- theta[g]
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE)
                    cat("\nWARNING: Evaluating",
                         Data[["parm.names"]][j], "at",
                         round(prop[j],5), "failed.\n", file=LogFile,
                         append=TRUE)
               Mo1 <- Mo0}
          LP.grid[g] <- Mo1[["LP"]]
          theta[g] <- Mo1[["parm"]][j]}
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     prop[j] <- sample(theta, 1, prob=LP.grid)
     Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
     if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5), "failed.\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0
          }
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5),
                    "resulted in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0}
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Discrete Parameter (Parallelized)
.mcmcggdpp <- function(Model, Data, j, Mo0, Grid, Debug, LogFile, cl)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     LIV <- length(Mo0[["parm"]])
     prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
     prop[, j] <- prop[, j] + Grid[[j]]
     Mo1 <- parLapply(cl, 1:G,
          function(x) Model(prop[x,], Data))
     LP.grid <- as.vector(unlist(lapply(Mo1,
          function(x) x[["LP"]])))
     prop <- matrix(as.vector(unlist(lapply(Mo1,
          function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
     theta <- prop[, j]
     prop <- Mo0[["parm"]]
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     prop[j] <- sample(theta, 1, prob=LP.grid)
     Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
     if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5), "failed.\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0
          }
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE)
               cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
                    "at", round(prop[j],5),
                    "resulted in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
          Mo1 <- Mo0}
     Mo0 <- Mo1
     return(Mo0)
     }
.mcmcharm <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     Debug, LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     Block <- Specs[["B"]]
     if(is.na(alpha.star) & {length(Block) == 0}) {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile,
                         append=TRUE)
               ### Propose new parameter values
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1) * d
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal resulted in non-finite",
                              "value(s).\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + 1}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=.colVars(thinned))
          return(out)
          }
     else if(is.na(alpha.star) & {length(Block) > 0}) {
          B <- length(Block)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile,
                         append=TRUE)
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    theta <- rnorm(length(Block[[b]]))
                    d <- theta / sqrt(sum(theta*theta))
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- prop[Block[[b]]] + runif(1) * d
                    if({b == 1} & {iter %% Status == 0}) 
                         cat(",   Proposal: Blockwise,   LP: ",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed.\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal resulted in",
                                   "non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         log.u <- log(runif(1))
                         log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                         if(!is.finite(log.alpha)) log.alpha <- 0
                         if(log.u < log.alpha) {
                              Mo0 <- Mo1
                              Acceptance <- Acceptance +
                                   length(Block[[b]]) / LIV}}
                    }
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=.colVars(thinned))
          return(out)
          }
     else if(length(Block) == 0) {
          tau <- 1
          DiagCovar <- matrix(tau, nrow(thinned), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile,
                         append=TRUE)
               ### Propose new parameter values
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1,0,tau) * d
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal resulted in non-finite",
                              "value(s).\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + 1
                         tau <- tau + (tau / (alpha.star *
                              (1 - alpha.star))) * (1 - alpha.star) / iter
                         }
                    else {
                         tau <- abs(tau - (tau / (alpha.star *
                              (1 - alpha.star))) * alpha.star / iter)}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]
                    DiagCovar[t.iter,] <- tau}
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=.colVars(thinned))
          return(out)
          }
     else {
          B <- length(Block)
          tau <- rep(1,B)
          DiagCovar <- matrix(1, nrow(thinned), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile,
                         append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    theta <- rnorm(length(Block[[b]]))
                    d <- theta / sqrt(sum(theta*theta))
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- prop[Block[[b]]] +
                         runif(1,0,tau[b]) * d
                    if({b == 1} & {iter %% Status == 0}) 
                         cat(",   Proposal: Blockwise,   LP: ",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed.\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal resulted in",
                                   "non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         log.u <- log(runif(1))
                         log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                         if(!is.finite(log.alpha)) log.alpha <- 0
                         if(log.u < log.alpha) {
                              Mo0 <- Mo1
                              Acceptance <- Acceptance +
                                   length(Block[[b]]) / LIV
                              tau[b] <- tau[b] + (tau[b] / (alpha.star *
                                   (1 - alpha.star))) * (1 - alpha.star) / iter
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[["parm"]]
                                   Dev[t.iter] <- Mo1[["Dev"]]
                                   Mon[t.iter,] <- Mo1[["Monitor"]]
                                   DiagCovar[t.iter, Block[[b]]] <- tau[b]}
                              }
                         else {
                              tau[b] <- abs(tau[b] - (tau[b] / (alpha.star *
                                   (1 - alpha.star))) * alpha.star / iter)
                              if(iter %% Thinning == 0)
                                   DiagCovar[t.iter, Block[[b]]] <- tau[b]}
                         }
                    }
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=.colVars(thinned))
          return(out)
          }
     }
.mcmchmc <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     Debug, LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     m <- Specs[["m"]]
     invm <- as.inverse(m)
     U <- chol(m)
     gr0 <- partial(Model, Mo0[["parm"]], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum0 <- as.vector(rnorm(LIV) %*% U)
          kinetic0 <- t(momentum0) %*% invm %*% momentum0 / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + as.vector(epsilon %*% invm) * momentum1
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed in leapfrog", l,
                              ".\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal in leapfrog", l,
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1}
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed in leapfrog",
                                   l, ".\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0.1
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in leapfrog",
                                   l, "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0.1}}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- t(momentum1) %*% invm %*% momentum1 / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=matrix(epsilon, 1, LIV),
          Mon=Mon,
          thinned=thinned,
          VarCov=cov(thinned))
     return(out)
     }
.mcmchmcda <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     A <- Specs[["A"]]
     delta <- Specs[["delta"]]
     epsilon <- Specs[["epsilon"]]
     Lmax <- Specs[["Lmax"]]
     lambda <- Specs[["lambda"]]
     leapfrog <- function(theta, r, grad, epsilon, Model, Data, Mo0, Debug)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- try(Model(thetaprime, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed in leapfrog.\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(thetaprime, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal in leapfrog", 
                         "resulted in non-finite value(s).\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(thetaprime, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0}
          thetaprime <- Mo1[["parm"]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime,
               rprime=rprime,
               gradprime=gradprime,
               Mo1=Mo1)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data,
          LogFile)
          {
          cat("\nFinding a reasonable initial value for epsilon...",
               file=LogFile, append=TRUE)
          epsilon <- 0.001
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data, Mo0,
               Debug)
          if(!is.finite(leap$Mo1[["LP"]]))
               stop("LP is not finite in find.reasonable.epsilon().",
                    file=LogFile, append=TRUE)
          acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data,
                    Mo0, Debug)
               if(!is.finite(leap$Mo1[["LP"]]))
                    stop("LP is not finite in find.reasonable.epsilon().",
                         file=LogFile, append=TRUE)
               acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="",
               file=LogFile, append=TRUE)
          return(epsilon)
          }
     gr0 <- partial(Model, Mo0[["parm"]], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(Mo0[["parm"]], gr0, Mo0, Model,
               Data, LogFile)
     DiagCovar[1,] <- epsilon
     L <- max(1, round(lambda / epsilon))
     L <- min(L, Lmax)
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Begin HMCDA
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum1 <- momentum0 <- runif(LIV)
          joint <- Mo0[["LP"]] - 0.5 * as.vector(momentum0 %*% momentum0)
          L <- max(1, round(lambda / epsilon))
          L <- min(L, Lmax)
          gr1 <- gr0
          Mo0.1 <- Mo0
          ### Leapfrog Function
          for (l in 1:L) {
               momentum1 <- momentum1 + 0.5 * epsilon * gr1
               prop <- prop + epsilon * momentum1
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed in leapfrog", l,
                              ".\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""),"\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal in leapfrog", l, 
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1}
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed in leapfrog",
                                   l, ".\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0.1
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in leapfrog", l, 
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0.1}}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               momentum1 <- momentum1 + epsilon * gr1}
          ### Accept/Reject
          alpha <- min(1,
               exp(prop - 0.5 * as.vector(momentum1 %*% momentum1) - joint))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               gr0 <- gr1
               Acceptance <- Acceptance + 1}
          ### Adaptation
          if(iter > 1) {
               eta <- 1 / (iter - 1 + t0)
               Hbar <- (1 - eta) * Hbar + eta * (delta - alpha)
               if(iter <= A) {
                    epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
                    eta <- (iter-1)^-kappa
                    epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                         eta * log(epsilon))
                    DiagCovar <- rbind(DiagCovar, epsilon)}
               else epsilon <- epsilonbar}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcim <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     Debug, LogFile)
     {
     mu <- Specs[["mu"]]
     VarCov <- as.positive.definite(as.symmetric.matrix(VarCov * 1.1))
     Omega <- as.inverse(VarCov)
     U <- chol(VarCov)
     d <- eigen(VarCov, symmetric=TRUE)$values
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% U, silent=TRUE)
          if(!inherits(MVNz, "try-error")) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate,   LP: ",
                        round(Mo0[["LP"]],1), "\n", sep="",
                        file=LogFile, append=TRUE)
               prop <- as.vector(mu) + as.vector(MVNz)}
          else {prop <- as.vector(Mo0[["parm"]])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          ### Importance Densities (dmvn)
          ss <- prop - mu
          z <- rowSums({ss %*% Omega} * ss)
          d1 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
          ss <- Mo0[["parm"]] - mu
          z <- rowSums({ss %*% Omega} * ss)
          d0 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]] + d1 - d0
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=cov(thinned))
     return(out)
     }
.mcmcinca <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     con <- get("con")
     Chains <- get("Chains")
     DiagCovar <- matrix(0, floor(Iterations/Periodicity), LIV)
     ### Store all posteriors
     INCA_iter <- 1
     INCA_first <- TRUE
     tmpMean <- numeric(LIV)
     tmpCov <- matrix(0, LIV, LIV)
     tmpAlpha <- numeric(Periodicity)
     lambda <- ScaleF
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov),
               silent=!Debug[["DB.chol"]])
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               prop <- as.vector(post[iter,]) + as.vector(MVNz)}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               if(Debug[["DB.chol"]] == TRUE)
                    cat("\nWARNING: Cholesky decomposition failed for",
                         "proposal in iteration", iter, ".\n",
                         file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else {
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo1[["parm"]]
                    Acceptance <- Acceptance + 1}}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Save log.alpha
          if({iter %% Periodicity} == 0)
               tmpAlpha[Periodicity] <- min(1, exp(log.alpha))
          else tmpAlpha[(iter %% Periodicity)] <- min(1, exp(log.alpha))
          ### Shrinkage of Adaptive Proposal Variance
          if({iter < Adaptive} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               select_post <- cbind(post[(iter-Periodicity+1):iter,],
                    tmpAlpha)
               ### Ask for last posteriors to hpc_server
               tmp <- unserialize(con)
               ### Send new posteriors matrix to hpc_server      
               serialize(select_post, con)               
               if(is.matrix(tmp) && INCA_first == FALSE) {
                    for (i in 1:nrow(select_post)) {
                         tmpMean <- tmpMean + 1/(INCA_iter+1) *
                              (select_post[i, 1:LIV]-tmpMean)
                         tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
                              1/INCA_iter *
                              tcrossprod(select_post[i, 1:LIV]-tmpMean)
                         INCA_iter <- INCA_iter + 1}
                    for (i in 1:nrow(tmp)) {
                         tmpMean <- tmpMean + 1/(INCA_iter+1) *
                              (tmp[i, 1:LIV]-tmpMean)
                         tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
                              1/INCA_iter * 
                              tcrossprod(tmp[i, 1:LIV]-tmpMean)
                         INCA_iter <- INCA_iter + 1}
                    eta <- INCA_iter^-0.6
                    m1 <- median(select_post[, LIV+1])
                    m2 <- median(tmp[, LIV+1])
                    lambda <- exp(log(lambda) + eta * (m1 - 0.234))
                    lambda <- exp(log(lambda) + eta * (m2 - 0.234))}
               if(INCA_first == TRUE) {
                    for (i in 1:iter) {
                         tmpMean <- tmpMean + 1/(INCA_iter+1) *
                              (post[i, ]-tmpMean)
                         tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
                              1/INCA_iter * 
                              tcrossprod(post[i, ] - tmpMean)
                         INCA_iter <- INCA_iter + 1}
                    INCA_first <- FALSE}
               VarCov <- lambda * (tmpCov + 1e-9 * Iden.Mat)
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               ### Univariate Standard Deviations
               tuning <- sqrt(diag(VarCov))}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
                 Dev=Dev,
                 DiagCovar=DiagCovar,
                 Mon=Mon,
                 thinned=thinned,
                 VarCov=VarCov)
     return(out)
     }
.mcmcmala <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     VarCov, Debug, LogFile)
     {
     A <- Specs[["A"]]
     alpha.star <- Specs[["alpha.star"]]
     delta <- Specs[["delta"]]
     gamma.const <- Specs[["gamma"]]
     epsilon <- Specs[["epsilon"]]
     Gamm <- as.positive.definite(VarCov)
     mu <- Mo0[["parm"]]
     sigma2 <- 1 / (LIV*LIV)
     DiagCovar <- matrix(diag(Gamm), nrow(thinned), LIV)
     Iden <- diag(LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Propose new parameter values
          gr <- partial(Model, Mo0[["parm"]], Data)
          Dx <- {delta/max(delta, abs(gr))}*gr
          gamm <- min(gamma.const/iter, 1)
          Lambda <- Gamm + epsilon[2]*Iden
          U <- try(chol(sigma2*Lambda), silent=!Debug[["DB.chol"]])
          if(inherits(U, "try-error")) {
               if(Debug[["DB.chol"]] == TRUE)
                    cat("\nWARNING: Cholesky decomposition failed for",
                         "proposal in iteration", iter, ".\n",
                         file=LogFile, append=TRUE)
               U <- chol(as.positive.definite(sigma2*Lambda))}
          prop <- as.vector((Mo0[["parm"]] +
               {sigma2/2}*as.vector(Lambda %*% Dx)*Dx) +
               rbind(rnorm(LIV)) %*% U)
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else {
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1}}
          ### Adapt Gamma (first, since it uses mu[t] not [t+1])
          xmu <- Mo0[["parm"]] - mu
          Gamm.prop <- Gamm + gamm*{xmu %*% t(xmu) - Gamm}
          norm.Gamm <- norm(Gamm.prop, type="F")
          if(norm.Gamm <= A) Gamm <- Gamm.prop
          else if(!is.finite(norm.Gamm)) Gamm <- sigma2*Iden
          else Gamm <- {A/norm.Gamm}*Gamm.prop
          ### Adapt mu
          mu.prop <- mu + gamm*(Mo0[["parm"]] - mu)
          norm.mu <- sqrt(sum(mu.prop*mu.prop))
          if(norm.mu <= A) mu <- mu.prop
          else mu <- {A/norm.mu}*mu.prop
          ### Adapt sigma
          sigma2 <- interval(sqrt(sigma2) +
               gamm*(min(exp(log.alpha),1) - alpha.star),
               epsilon[1], A, reflect=FALSE)^2
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]
               DiagCovar[t.iter,] <- diag(Lambda)}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=Lambda)
     return(out)
     }
.mcmcmcmcmc <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     lambda <- Specs[["lambda"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     detectedCores <- detectCores()
     cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
          append=TRUE)
     if(CPUs > detectedCores) {
          cat("\nOnly", detectedCores, "will be used.\n",
               file=LogFile, append=TRUE)
          CPUs <- detectedCores}
     cat("\nLaplace's Demon is preparing environments for CPUs...",
          file=LogFile, append=TRUE)
     cat("\n##################################################\n",
          file=LogFile, append=TRUE)
     cl <- makeCluster(CPUs)
     cat("\n##################################################\n",
          file=LogFile, append=TRUE)
     on.exit(stopCluster(cl))
     varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
          ls(envir=parent.env(environment()))))
     clusterExport(cl, varlist=varlist, envir=environment())
     clusterSetRNGStream(cl)
     wd <- getwd()
     clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
          envir=environment())
     if(length(lambda) == 1) Temperature <- 1/(1 + lambda*(c(1:CPUs) - 1))
     else if(length(lambda) == LIV) Temperature <- lambda
     else Temperature <- 1/(1 + lambda[1]*(c(1:CPUs) - 1))
     coolest <- which.max(Temperature)[1]
     temp <- Mo0
     Mo0 <- list()
     for (i in 1:CPUs) Mo0[[i]] <- temp
     prop <- matrix(Mo0[[1]][["parm"]], CPUs, LIV, byrow=TRUE)
     Acceptance.swap <- 0
     U <- chol(VarCov)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[[coolest]][["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[[coolest]][["parm"]]
               Dev[t.iter] <- Mo0[[coolest]][["Dev"]]
               Mon[t.iter,] <- Mo0[[coolest]][["Monitor"]]}
          ### Propose new parameter values
          for (i in 1:CPUs)
               prop[i,] <- Mo0[[i]][["parm"]] + rbind(rnorm(LIV)) %*% U
          ### Log-Posterior of the proposed state
          Mo1 <- parLapply(cl, 1:CPUs, function(x)
               try(Model(prop[x,], Data), silent=!Debug[["DB.Model"]]))
          for (i in 1:CPUs) {
               if(inherits(Mo1[[i]], "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed in chain", i,
                              ".\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop[i,], collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1[[i]] <- Mo0[[i]]
                    }
               else if(any(!is.finite(c(Mo1[[i]][["LP"]], Mo1[[i]][["Dev"]],
                    Mo1[[i]][["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal in chain", i,
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop[i,], collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1[[i]] <- Mo0[[i]]
                    }
               }
          ### Accept/Reject
          for (i in 1:CPUs) {
               log.u <- log(runif(1))
               log.alpha <- (Mo1[[i]][["LP"]] - Mo0[[i]][["LP"]]) /
                    Temperature[i]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0[[i]] <- Mo1[[i]]
                    if(i == coolest) {
                         Acceptance <- Acceptance + 1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[[i]][["parm"]]
                              Dev[t.iter] <- Mo1[[i]][["Dev"]]
                              Mon[t.iter,] <- Mo1[[i]][["Monitor"]]}}}}
          ### Swap
          swap <- sample.int(CPUs, 2)
          log.u <- log(runif(1))
          log.alpha <- {(Mo0[[swap[1]]][["LP"]] - Mo0[[swap[2]]][["LP"]]) /
               Temperature[swap[2]]} +
               {(Mo0[[swap[2]]][["LP"]] - Mo0[[swap[1]]][["LP"]]) /
               Temperature[swap[1]]}
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Acceptance.swap <- Acceptance.swap + 1
               temp <- Mo0[[swap[2]]]
               Mo0[[swap[2]]] <- Mo0[[swap[1]]]
               Mo0[[swap[1]]] <- temp
               if({swap[1] == coolest} & {iter %% Thinning == 0}) {
                    thinned[t.iter,] <- Mo0[[swap[1]]][["parm"]]
                    Dev[t.iter] <- Mo0[[swap[1]]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[swap[1]]][["Monitor"]]}}
          }
     cat("\nSwap Acceptance Rate:",
          round(Acceptance.swap / Iterations, 5), "\n", file=LogFile,
          append=TRUE)
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcmtm <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, tuning, Debug,
     LogFile)
     {
     K <- Specs[["K"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     if(CPUs > 1) {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
               append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n",
                    file=LogFile, append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...",
               file=LogFile, append=TRUE)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n",
               file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
               ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())}
     Acceptance <- matrix(0, 1, LIV)
     Mo1 <- list()
     for (k in 1:K) Mo1[[k]] <- Mo0
     LW <- LP <- rep(0, K)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Random-Scan Componentwise Estimation
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop1 <- matrix(Mo0[["parm"]], K, LIV, byrow=TRUE)
               prop1[,j] <- rnorm(K, prop1[,j], tuning[j])
               ### Log-Posterior of the proposed states
               if(CPUs == 1) {
                    ### Non-parallel
                    for (k in 1:K) {
                         Mo1[[k]] <- try(Model(prop1[k,], Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1[[k]], "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "failed.\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop1[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]],
                              Mo1[[k]][["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "resulted in non-finite",
                                        "value(s).\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop1[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0}
                         LP[k] <- LW[k] <- Mo1[[k]][["LP"]]
                         prop1[k,] <- Mo1[[k]][["parm"]]}
                    }
               else {
                    ### Parallel
                    Mo1 <- parLapply(cl, 1:K, function(x)
                         try(Model(prop1[x,], Data),
                              silent=!Debug[["DB.Model"]]))
                    for (k in 1:K) {
                         if(inherits(Mo1[[k]], "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "failed.\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop1[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]],
                              Mo1[[k]][["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "resulted in non-finite",
                                        "value(s).\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop1[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0}
                         LP[k] <- LW[k] <- Mo1[[k]][["LP"]]
                         prop1[k,] <- Mo1[[k]][["parm"]]}
                    }
               ### Normalize Weights
               w <- exp(LW - logadd(LW))
               if(all(w == 0)) w <- rep(1/K, K)
               ### Sample a Proposal
               prop5 <- Mo0[["parm"]]
               prop2 <- sample(prop1[,j], size=1, prob=w)
               prop5[j] <- prop2
               ### Create Reference Set
               Mo2 <- try(Model(prop5, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo2, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop5[j],5),
                              file=LogFile, append=TRUE)}
                    Mo2 <- Mo0}
               prop3 <- c(rnorm(K-1, Mo2[["parm"]][j], tuning[j]),
                    Mo2[["parm"]][j])
               prop4 <- prop1
               prop4[,j] <- prop3
               ### Calculate Acceptance Probability
               numerator <- logadd(LP)
               denom <- rep(0, K)
               if(CPUs == 1) {
                    ### Non-parallel
                    for (k in 1:K) {
                         Mo1[[k]] <- try(Model(prop4[k,], Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1[[k]], "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "failed.\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop4[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]],
                              Mo1[[k]][["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "resulted in non-finite",
                                        "value(s).\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop4[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0}
                         denom[k] <- Mo1[[k]][["LP"]]}
                    }
               else {
                    ### Parallel
                    Mo1 <- parLapply(cl, 1:K, function(x)
                         try(Model(prop4[x,], Data),
                              silent=!Debug[["DB.Model"]]))
                    for (k in 1:K) {
                         if(inherits(Mo1[[k]], "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "failed.\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop4[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]],
                              Mo1[[k]][["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal ", k,
                                        "resulted in non-finite",
                                        "value(s).\n", file=LogFile,
                                        append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop4[k,j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1[[k]] <- Mo0}
                         denom[k] <- Mo1[[k]][["LP"]]}}
               denom <- logadd(denom)
               ### Accept/Reject
               u <- log(runif(1)) < (numerator - denom)
               if(u == TRUE) {
                    Mo0 <- Mo2
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcmwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     Debug, LogFile)
     {
     Block <- Specs[["B"]]
     B <- length(Block)
     Acceptance <- matrix(0, 1, LIV)
     if(B == 0) {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Random-Scan Componentwise Estimation
               propdraw <- rnorm(LIV,0,tuning)
               for (j in sample.int(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + propdraw[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed for",
                                   Data[["parm.names"]][j], ".\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal for",
                                   Data[["parm.names"]][j],
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter,
                                   "Current:", round(Mo0[["parm"]][j]),
                                   "Proposed:", round(prop[j],5),
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         u <- log(runif(1)) < {Mo1[["LP"]] - Mo0[["LP"]]}
                         if(u == TRUE) {
                              Mo0 <- Mo1
                              Acceptance[j] <- Acceptance[j] + 1}}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     else {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Componentwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Random-Scan Componentwise Estimation
               propdraw <- rnorm(LIV,0,tuning)
               ### Proceed by Block
               for (b in 1:B) {
                    for (j in sample(Block[[b]])) {
                         ### Propose new parameter values
                         prop <- Mo0[["parm"]]
                         prop[j] <- prop[j] + propdraw[j]
                         ### Log-Posterior of the proposed state
                         Mo1 <- try(Model(prop, Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal failed for",
                                        Data[["parm.names"]][j], ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop[j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                              Mo1[["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Proposal for",
                                        Data[["parm.names"]][j],
                                        "resulted in non-finite value(s).\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter,
                                        "Current:", round(Mo0[["parm"]][j]),
                                        "Proposed:", round(prop[j],5),
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0
                              }
                         else {
                              ### Accept/Reject
                              u <- log(runif(1)) < {Mo1[["LP"]] -
                                   Mo0[["LP"]]}
                              if(u == TRUE) {
                                   Mo0 <- Mo1
                                   Acceptance[j] <- Acceptance[j] + 1}}}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
.mcmcnuts <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     A <- Specs[["A"]]
     delta <- Specs[["delta"]]
     epsilon <- Specs[["epsilon"]]
     Lmax <- Specs[["Lmax"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     leapfrog <- function(theta, r, grad, epsilon, Model, Data, Mo0, Debug)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- try(Model(thetaprime, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed in leapfrog.\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(thetaprime, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal in leapfrog", 
                         "resulted in non-finite value(s).\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(thetaprime, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0}
          thetaprime <- Mo1[["parm"]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime,
               rprime=rprime,
               gradprime=gradprime,
               Mo1=Mo1)
          return(out)
          }
     stop.criterion <- function(thetaminus, thetaplus, rminus, rplus)
          {
          thetavec <- thetaplus - thetaminus
          criterion <- (thetavec %*% rminus >= 0) &&
               (thetavec %*% rplus >= 0)
          return(criterion)
          }
     build.tree <- function(theta, r, grad, logu, v, j, epsilon, joint0, Mo0)
          {
          if(j == 0) {
               ### Base case: Take a single leapfrog step in direction v
               leap <- leapfrog(theta=theta, r=r, grad=grad,
                    epsilon=v*epsilon, Model=Model, Data=Data, Mo0=Mo0,
                    Debug=Debug)
               rprime <- leap$rprime
               thetaprime <- leap$thetaprime
               Mo1 <- leap$Mo1
               gradprime <- leap$gradprime
               joint <- Mo1[["LP"]] - 0.5 * as.vector(rprime %*% rprime)
               ### Is the new point in the slice?
               nprime <- logu < joint
               ### Is the simulation wildly inaccurate?
               sprime <- logu - 1000 < joint
               # Set the return values---minus=plus for all things here,
               # since the "tree" is of depth 0.
               thetaminus <- thetaprime
               thetaplus <- thetaprime
               rminus <- rprime
               rplus <- rprime
               gradminus <- gradprime
               gradplus <- gradprime
               ### Compute the acceptance probability
               alphaprime <- min(1, exp(Mo1[["LP"]] - 0.5 *
                    as.vector(rprime %*% rprime) - joint0))
               nalphaprime <- 1}
          else {
               # Recursion: Implicitly build the height j-1 left and
               # right subtrees
               tree <- build.tree(theta=theta, r=r, grad=grad, logu=logu,
                    v=v, j=j-1, epsilon=epsilon, joint=joint0, Mo0=Mo0)
               thetaminus <- tree$thetaminus
               rminus <- tree$rminus
               gradminus <- tree$gradminus
               thetaplus <- tree$thetaplus
               rplus <- tree$rplus
               gradplus <- tree$gradplus
               thetaprime <- tree$thetaprime
               gradprime <- tree$gradprime
               Mo1 <- tree$Mo1
               nprime <- tree$nprime
               sprime <- tree$sprime
               alphaprime <- tree$alphaprime
               nalphaprime <- tree$nalphaprime
               ### If the first subtree stopping criterion is met, then stop
               if(sprime == 1) {
                    if(v == -1) {
                         tree <- build.tree(theta=thetaminus, r=rminus,
                              grad=gradminus, logu=logu, v=v, j=j-1,
                              epsilon=epsilon, joint0=joint0, Mo0=Mo0)
                         thetaminus <- tree$thetaminus
                         rminus <- tree$rminus
                         gradminus <- tree$gradminus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    else {
                         tree <- build.tree(theta=thetaplus, r=rplus,
                              grad=gradplus, logu=logu, v=v, j=j-1,
                              epsilon=epsilon, joint0=joint0, Mo0=Mo0)
                         thetaplus <- tree$thetaplus
                         rplus <- tree$rplus
                         gradplus <- tree$gradplus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    ### Choose a subtree to propagate a sample up from
                    temp <- nprime2 / (nprime + nprime2)
                    if(!is.finite(temp)) temp <- 0
                    if(runif(1) < temp) {
                         thetaprime <- thetaprime2
                         gradprime <- gradprime2
                         Mo1 <- Mo12}
                    ### Update the number of valid points
                    nprime <- nprime + nprime2
                    ### Update the stopping criterion
                    sprime <- sprime && sprime2 &&
                         stop.criterion(thetaminus, thetaplus, rminus,
                              rplus)
                    ### Update acceptance probability statistics
                    alphaprime <- alphaprime + alphaprime2
                    nalphaprime <- nalphaprime + nalphaprime2}}
          out <- list(thetaminus=thetaminus,
               rminus=rminus,
               gradminus=gradminus,
               thetaplus=thetaplus,
               rplus=rplus,
               gradplus=gradplus,
               thetaprime=thetaprime,
               gradprime=gradprime,
               Mo1=Mo1,
               nprime=nprime,
               sprime=sprime,
               alphaprime=alphaprime,
               nalphaprime=nalphaprime)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data,
          LogFile)
          {
          cat("\nFinding a reasonable initial value for epsilon...",
               file=LogFile, append=TRUE)
          epsilon <- 0.001
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta=theta0, r=r0, grad=grad0,
               epsilon=epsilon, Model=Model, Data=Data, Mo0=Mo0,
               Debug=Debug)
          if(!is.finite(leap$Mo1[["LP"]]))
               stop("LP is not finite in find.reasonable.epsilon().",
                    file=LogFile, append=TRUE)
          acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta=theta0, r=r0, grad=grad0,
                    epsilon=epsilon, Model=Model, Data=Data, Mo0=Mo0,
                    Debug=Debug)
               if(!is.finite(leap$Mo1[["LP"]]))
                    stop("LP is not finite in find.reasonable.epsilon().",
                         file=LogFile, append=TRUE)
               acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="",
               file=LogFile, append=TRUE)
          return(epsilon)
          }
     Count <- 0
     evals <- 0
     grad <- partial(Model, post[1,], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(theta0=post[1,], grad0=grad,
               Mo0=Mo0, Model=Model, Data=Data, LogFile=LogFile)
     DiagCovar[1,] <- epsilon
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Reset Dev, Mon, and thinned
     if(A < Iterations) {
          Dev <- matrix(Dev[1:(floor((Iterations-A)/Thinning)+1),])
          Mon <- matrix(Mo0[["Monitor"]], floor((Iterations-A)/Thinning)+1,
               length(Mo0[["Monitor"]]), byrow=TRUE)
          thinned <- matrix(0, floor((Iterations-A)/Thinning)+1, LIV)}
     ### Begin NUTS
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter > A) {
               if((iter-A) %% Thinning == 0) {
                    thinned[((iter-A)/Thinning+1),] <- post[iter,]
                    Dev[((iter-A)/Thinning+1)] <- Mo0[["Dev"]]
                    Mon[((iter-A)/Thinning+1),] <- Mo0[["Monitor"]]}}
          else if(A >= Iterations) {
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- post[iter,]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}}
          prop <- post[iter,]
          r0 <- runif(LIV) ### r0 is momenta
          ### Joint log-probability of theta and momenta r
          joint <- Mo0[["LP"]] - 0.5 * as.vector(r0 %*% r0)
          ### Resample u ~ U([0, exp(joint)])
          logu <- joint - rexp(1)
          ### Initialize Tree
          thetaminus <- prop
          thetaplus <- prop
          rminus <- r0
          rplus <- r0
          gradminus <- grad
          gradplus <- grad
          j <- 0 ### Initial height j=0
          n <- 1 ### Initially, the only valid point is the initial point
          s <- 1 ### Loop until s == 0
          while (s == 1) {
               ### Choose a direction: -1=backward, 1=forward.
               v <- 2*(runif(1) < 0.5) - 1
               ### Double the size of the tree.
               if(v == -1) {
                    tree <- build.tree(theta=thetaminus, r=rminus,
                         grad=gradminus, logu=logu, v=v, j=j,
                         epsilon=epsilon, joint0=joint, Mo0=Mo0)
                    thetaminus <- tree$thetaminus
                    rminus <- tree$rminus
                    gradminus <- tree$gradminus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               else {
                    tree <- build.tree(theta=thetaplus, r=rplus,
                         grad=gradplus, logu, v=v, j=j, epsilon=epsilon,
                         joint0=joint, Mo0=Mo0)
                    thetaplus <- tree$thetaplus
                    rplus <- tree$rplus
                    gradplus <- tree$gradplus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               ### Accept/Reject
               Count <- Count + 1
               if((sprime == 1) && (runif(1) < nprime/n)) {
                    post[iter,] <- thetaprime
                    Mo0 <- Mo1
                    grad <- gradprime
                    Acceptance <- Acceptance + 1
                    if(iter > A) {
                         if((iter-A) %% Thinning == 0) {
                              thinned[((iter-A)/Thinning+1),] <- Mo1[["parm"]]
                              Dev[((iter-A)/Thinning+1)] <- Mo1[["Dev"]]
                              Mon[((iter-A)/Thinning+1),] <- Mo1[["Monitor"]]}}
                    else if(A >= Iterations) {
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}}}
               ### Update number of observed valid points
               n <- n + nprime
               ### Decide if it is time to stop
               s <- sprime &&
                    stop.criterion(thetaminus, thetaplus, rminus, rplus)
               ### Increment depth
               j <- j + 1
               if(j*j >= Lmax) s <- 0}
          ### Adaptation of epsilon
          eta <- 1 / (iter - 1 + t0)
          Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
          if(iter <= A) {
               epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
               eta <- (iter-1)^-kappa
               epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                    eta * log(epsilon))
               DiagCovar <- rbind(DiagCovar, epsilon)}
          else epsilon <- epsilonbar
          }
     Acceptance <- round(Acceptance / Count * Iterations)
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcohss <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     VarCov, Debug, LogFile)
     {
     A <- Specs[["A"]]
     n <- Specs[["n"]]
     w <- 0.05 # as with Roberts & Rosenthal
     if(!is.symmetric.matrix(VarCov)) {
               cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
                    append=TRUE)
               VarCov <- as.symmetric.matrix(VarCov)}
          if(!is.positive.definite(VarCov)) {
               cat("\nNon-Positive-Definite Covar, correcting now...\n",
                    file=LogFile, append=TRUE)
               VarCov <- as.positive.definite(VarCov)}
     decomp.freq <- max(floor(Iterations / Thinning / 100), 10)
     cat("\nEigendecomposition will occur every", decomp.freq,
          "iterations.\n\n", file=LogFile, append=TRUE)
     S.eig <-try(eigen(VarCov), silent=!Debug[["DB.eigen"]])
     if(inherits(S.eig, "try-error")) {
          if(Debug[["DB.eigen"]] == TRUE)
               cat("\nWARNING: Eigendecomposition failed.\n",
                    file=LogFile, append=TRUE)
          S.eig <- NULL
          DiagCovar <- matrix(0, floor(Iterations/Thinning)+1, LIV)
          }
     else DiagCovar <- matrix(diag(S.eig$vectors),
               floor(Iterations/Thinning)+1, LIV, byrow=TRUE)
     tuning <- 1 #Tuning
     edge.scale <- 5 #Tuning
     if(A > Iterations)
          post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     else post <- matrix(Mo0[["parm"]], A, LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Eigenvectors of the Sample Covariance Matrix
          if({iter %% decomp.freq == 0} & {iter > 2} & {iter <= A}) {
               VarCov2 <- try({VarCov*n +
                    cov(post[1:(iter-1),,drop=FALSE])*(iter-1)}/{n+iter-1},
                    silent=TRUE)
               if(inherits(VarCov2, "try-error")) VarCov2 <- VarCov
               if(!is.symmetric.matrix(VarCov2))
                    VarCov2 <- as.symmetric.matrix(VarCov2)
               if(!is.positive.definite(VarCov2))
                    VarCov2 <- as.positive.definite(VarCov2)
               S.eig <- try(eigen(VarCov2), silent=!Debug[["DB.eigen"]])
               if(inherits(S.eig, "try-error")) {
                    if(Debug[["DB.eigen"]] == TRUE)
                         cat("\nWARNING: Eigendecomposition failed in",
                              "iteration", iter, ".\n",
                              file=LogFile, append=TRUE)
                    S.eig <- eigen(VarCov)}}
          ### Hypercube or Eigenvector
          if(runif(1) < w || is.null(S.eig)) {
               vals <- rep(tuning, LIV)
               vecs <- diag(1, nrow=LIV)
               }
          else {
               vals <- S.eig$values
               vecs <- S.eig$vectors}
          ### Slice Interval
          Mo0.1 <- try(Model(Mo0[["parm"]], Data),
               silent=!Debug[["DB.Model"]])
          if(inherits(Mo0.1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(Mo0[["parm"]], collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo0.1 <- Mo0}
          Mo0 <- Mo0.1
          y.slice <- Mo0[["LP"]] - rexp(1)
          L <- -1 * runif(LIV)
          U <- L + 1
          ### Rejection Sampling
          repeat {
               wt <- runif(LIV, min=L, max=U)
               v <- as.numeric(vecs %*% {edge.scale * wt * vals})
               prop <- Mo0[["parm"]] + v
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Rejection sampling failed.\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Rejection sampling resulted",
                              "in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0}
               if(Mo1[["LP"]] >= y.slice) break
               else if(all(abs(wt) < 1e-100)) {
                    Mo1 <- Mo0
                    break}
               L[wt < 0] <- wt[wt < 0]
               U[wt > 0] <- wt[wt > 0]}
          Mo0 <- Mo1
          if(iter <= A) post[iter,] <- Mo0[["parm"]]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]
               DiagCovar[t.iter,] <- diag(S.eig$vectors)}
          }
     if(A > 0) VarCov <- VarCov2
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcpcn <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     Debug, LogFile)
     {
     beta <- Specs[["beta"]]
     U <- chol(VarCov)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
               ",   Proposal: Multivariate,   LP: ",
               round(Mo0[["LP"]],1), "\n", sep="",
               file=LogFile, append=TRUE)
          ### Propose new parameter values
          prop <- as.vector(sqrt(1 - beta*beta)*Mo0[["parm"]] +
               beta*(rbind(rnorm(LIV)) %*% U))
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else {
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1}}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcram <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     Debug, LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     Block <- Specs[["B"]]
     Dist <- Specs[["Dist"]]
     gamma <- Specs[["gamma"]]
     n <- Specs[["n"]]
     B <- length(Block)
     if(B == 0) {
          if(!is.symmetric.matrix(VarCov)) {
               cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
                    append=TRUE)
               VarCov <- as.symmetric.matrix(VarCov)}
          if(!is.positive.definite(VarCov)) {
               cat("\nNon-Positive-Definite Covar, correcting now...\n",
                    file=LogFile, append=TRUE)
               VarCov <- as.positive.definite(VarCov)}
          Iden.Mat <- diag(LIV)
          S.z <- try(t(chol(VarCov)), silent=!Debug[["DB.chol"]])
          if(!inherits(S.z, "try-error")) {
               if(Debug[["DB.chol"]] == TRUE)
                    cat("\nWARNING: Cholesky decomposition failed for",
                         "proposal.\n", file=LogFile, append=TRUE)
               S <- S.z}
          else S <- Iden.Mat
          DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning)+1,
               LIV, byrow=TRUE)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Propose New Parameter Values
               if(Dist == "t") U <- rt(LIV, df=5)
               else U <- rnorm(LIV)
               prop <- Mo0[["parm"]] + rbind(U) %*% S
               ### Log-Posterior
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal resulted in non-finite",
                              "value(s).\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + 1}}
               ### Adaptation
               eta <- min(1, LIV*{n + iter}^(-gamma))
               VarCov.test <- S %*% {Iden.Mat +
                    eta*(min(1, exp(log.alpha)) - alpha.star) *
                    U %*% t(U) / sum(U*U)} %*% t(S)
               if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
                    !is.matrix(VarCov.test)) {VarCov.test <- VarCov}
               if(!is.symmetric.matrix(VarCov.test))
                    VarCov.test <- as.symmetric.matrix(VarCov.test)
               if(is.positive.definite(VarCov.test)) {
                    S.z <- try(t(chol(VarCov)), silent=!Debug[["DB.chol"]])
                    if(!inherits(S.z, "try-error")) {
                         VarCov <- VarCov.test
                         S <- S.z
                         }
                    else if(Debug[["DB.chol"]] == TRUE)
                         cat("\nWARNING: Cholesky decomposition failed for",
                              "proposal in iteration", iter, ".\n",
                              file=LogFile, append=TRUE)}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]
                    DiagCovar[t.iter,] <- diag(VarCov)}
               }
          }
     else {
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from ",
                    "number of blocks.", file=LogFile, append=TRUE)
          DiagCovar <- rep(0, LIV)
          Iden.Mat <- S <- S.z <- list()
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from ",
                         "block length.", file=LogFile, append=TRUE)
               if(!is.symmetric.matrix(VarCov[[b]])) {
                    cat("\nAsymmetric Covar block, correcting now...\n",
                         file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.symmetric.matrix(VarCov[[b]])}
               if(!is.positive.definite(VarCov[[b]])) {
                    cat("\nNon-Positive-Definite Covar block,",
                         "correcting now...\n", file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.positive.definite(VarCov[[b]])}
               Iden.Mat[[b]] <- diag(length(diag(VarCov[[b]])))
               S.z[[b]] <- try(t(chol(VarCov[[b]])),
                    silent=!Debug[["DB.chol"]])
               if(!inherits(S.z[[b]], "try-error")) S[[b]] <- S.z[[b]]
               else {
                    if(Debug[["DB.chol"]] == TRUE)
                         cat("\nWARNING: Cholesky decomposition failed for",
                              "proposal in block", b, ".\n", file=LogFile,
                              append=TRUE)
                    S[[b]] <- Iden.Mat[[b]]}
               DiagCovar[Block[[b]]] <- diag(VarCov[[b]])
          }
          DiagCovar <- matrix(DiagCovar, floor(Iterations/Thinning)+1,
               LIV, byrow=TRUE)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Blockwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose New Parameter Values
                    if(Dist == "t") U <- rt(length(Block[[b]]), df=5)
                    else U <- rnorm(length(Block[[b]]))
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- prop[Block[[b]]] +
                         rbind(U) %*% S[[b]]
                    ### Log-Posterior
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed in block",
                                   b, ".\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in block", b,
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         log.u <- log(runif(1))
                         log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                         if(!is.finite(log.alpha)) log.alpha <- 0
                         if(log.u < log.alpha) {
                              Mo0 <- Mo1
                              Acceptance <- Acceptance +
                                   length(Block[[b]]) / LIV}}
                    ### Adaptation
                    eta <- min(1, length(Block[[b]])*{n + iter}^(-gamma))
                    VarCov.test <- S[[b]] %*% {Iden.Mat[[b]] +
                         eta*(min(1, exp(log.alpha)) - alpha.star) *
                         U %*% t(U) / sum(U*U)} %*% t(S[[b]])
                    if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
                         !is.matrix(VarCov.test)) {VarCov.test <- VarCov[[b]]}
                    if(!is.symmetric.matrix(VarCov.test))
                         VarCov.test <- as.symmetric.matrix(VarCov.test)
                    if(is.positive.definite(VarCov.test)) {
                         S.z[[b]] <- try(t(chol(VarCov[[b]])),
                              silent=!Debug[["DB.chol"]])
                         if(!inherits(S.z[[b]], "try-error")) {
                              VarCov[[b]] <- VarCov.test
                              S[[b]] <- S.z[[b]]}
                         else if(Debug[["DB.chol"]] == TRUE)
                              cat("\nWARNING: Cholesky decomposition",
                                   "failed for proposal in block", b,
                                   "in iteration", iter, ".\n",
                                   file=LogFile, append=TRUE)}
                    DiagCovar[floor(iter / Thinning)+1,Block[[b]]] <- diag(VarCov[[b]])
                    }
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcrdmh <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     Debug, LogFile)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          s <- sample(c(-1,1), LIV, replace=TRUE)
          u1 <- runif(LIV, -1, 1)
          epsilon1 <- u1^s
          ### Random-Scan Componentwise Estimation
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- prop[j]*epsilon1[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed for",
                              Data[["parm.names"]][j], ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal for",
                              Data[["parm.names"]][j],
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    epsilon2 <- log(abs(Mo1[["parm"]][j] /
                         Mo0[["parm"]][j]))
                    if(!is.finite(epsilon2)) epsilon2 <- 0
                    ### Accept/Reject
                    u2 <- log(runif(1)) < (epsilon2 + Mo1[["LP"]] -
                         Mo0[["LP"]])
                    if(u2 == TRUE) {
                         Mo0 <- Mo1
                         Acceptance[j] <- Acceptance[j] + 1}}}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrefractive <- function(Model, Data, Iterations, Status, Thinning,
     Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, Debug,
     LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     m <- Specs[["m"]]
     w <- Specs[["w"]]
     r <- Specs[["r"]]
     alpha.star <- 0.65
     if(Adaptive < Iterations) DiagCovar <- matrix(w, nrow(thinned), LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          prop <- Mo0[["parm"]]
          p <- rnorm(LIV)
          a <- 1
          g <- partial(Model, prop, Data)
          for (i in 1:m) {
               if(t(p) %*% g > 0) {
                    u <- g / sqrt(sum(g*g))
                    r1 <- 1
                    r2 <- r
                    }
               else {
                    u <- -g / sqrt(sum(g*g))
                    r1 <- r
                    r2 <- 1}
               cos.theta.1 <- (t(p) %*% u) / sqrt(sum(p*p))
               cos.2.theta.1 <- cos.theta.1 * cos.theta.1
               cos.2.theta.2 <- 1 - (r1^2 / r2^2)*(1 - cos.2.theta.1)
               if(cos.2.theta.2 > 0) cos.theta.2 <- sqrt(cos.2.theta.2)
               else cos.theta.2 <- -sqrt(abs(cos.2.theta.2))
               if(cos.2.theta.2 < 0) p <- as.vector(p - 2*(t(p) %*% u) %*% u)
               else {
                    p <- (r1 / r2)*p -
                         sqrt(sum(p*p))*((r1 / r2)*cos.theta.1 -
                         cos.theta.2)*u
                    a <- (r1 / r2)^(LIV-1)*(cos.theta.1 / cos.theta.2)*a}
               prop <- prop + w*p
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal resulted in non-finite",
                              "value(s).\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               prop <- Mo1[["parm"]]}
          ### Accept/Reject
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]] + exp(a)
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log(runif(1)) < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(Adaptive < Iterations)
                    w <- w + (w / (alpha.star * (1 - alpha.star))) *
                         (1 - alpha.star) / iter
               }
          else if(Adaptive < Iterations)
               w <- abs(w - (w / (alpha.star * (1 - alpha.star))) *
                    alpha.star / iter)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]
               if(Adaptive < Iterations) DiagCovar[t.iter,] <- w}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrj <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     bin.n <- Specs[["bin.n"]]
     bin.p <- Specs[["bin.p"]]
     parm.p <- Specs[["parm.p"]]
     selectable <- Specs[["selectable"]]
     selected <- Specs[["selected"]]
     cur.parm <- cur.sel <- selected
     cur.parm[which(selectable == 0)] <- 1
     nonzero.post <- rep(0, LIV)
     p <- parm.p
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Propose a variable to include/exclude
          v.change <- sample(LIV, 1, prob=selectable)
          prop.sel <- cur.sel
          prop.parm <- cur.parm
          ### Change proposed size, but not above bin.n
          if(sum(cur.sel) < bin.n) {
               prop.sel[v.change] <- 1 - prop.sel[v.change]
               prop.parm[v.change] <- 1 - prop.parm[v.change]}
          else if(prop.sel[v.change] == 1) 
               prop.parm[v.change] <- prop.sel[v.change] <- 0
          ### Priors
          prior.cur <- sum(dbern(cur.sel, p[which(selectable == 1)], log=TRUE),
               dbinom(sum(cur.sel), bin.n, bin.p, log=TRUE))
          prior.prop <- sum(dbern(prop.sel, p[which(selectable == 1)], log=TRUE),
               dbinom(sum(prop.sel), bin.n, bin.p, log=TRUE))
          ### Hit-And-Run Proposal Parameters
          theta <- rnorm(LIV)
          theta <- theta / sqrt(sum(theta*theta))
          lambda <- runif(1)
          ### Random-Scan Componentwise Estimation (Within-Model)
          for (j in sample(which(cur.parm == 1))) {
               ### Propose new parameter values
               temp.post <- Mo0[["parm"]]
               temp.post[which(temp.post == 0)] <- nonzero.post[which(temp.post == 0)]
               temp.post[which(cur.parm == 0)] <- 0
               prop <- Mo0[["parm"]] <- temp.post
               prop[j] <- prop[j] + lambda*theta[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Within-model proposal failed for",
                              Data[["parm.names"]][j], ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Within-model proposal for",
                              Data[["parm.names"]][j],
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               ### Accept/Reject (Within-Model Move)
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               if(Mo0[["parm"]][j] != 0) nonzero.post[j] <- Mo0[["parm"]][j]
               Acceptance <- Acceptance + (u * (1 / sum(cur.parm)))}
          ### Random-Scan Componentwise Estimation (Between-Models)
          prop <- Mo0[["parm"]]
          prop[v.change] <- prop.sel[v.change]*(prop[v.change] +
               lambda*theta[v.change])
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Between-models proposal failed for",
                         Data[["parm.names"]][j], ".\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter,
                         "Current:", round(Mo0[["parm"]][j]),
                         "Proposed:", round(prop[j],5),
                         file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Between-models proposal for",
                         Data[["parm.names"]][j],
                         "resulted in non-finite value(s).\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter,
                         "Current:", round(Mo0[["parm"]][j]),
                         "Proposed:", round(prop[j],5),
                         file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          ### Accept/Reject (Between-Models Move)
          u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]] + prior.prop -
               prior.cur)
          if(u == TRUE) {
               Mo0 <- Mo1
               cur.sel <- prop.sel
               cur.parm <- prop.parm}
          if(Mo0[["parm"]][v.change] != 0)
               nonzero.post[v.change] <- Mo0[["parm"]][v.change]
          Acceptance <- Acceptance + (u * (1 / sum(prop.parm)))
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrss <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, Debug, LogFile)
     {
     m <- Specs[["m"]]
     w <- Specs[["w"]]
     reflections <- 0
     Norm <- function(x) return(sqrt(sum(x*x)))
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          prop <- Mo0[["parm"]]
          y.slice <- Mo0[["LP"]] - rexp(1)
          g <- partial(Model, prop, Data)
          p <- rnorm(LIV)
          reflections <- 0
          ### Take m Steps
          for (i in 1:m) {
               prop0 <- prop
               prop <- prop + w*p
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Stepping out proposal failed",
                              "in step", i, ".\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Stepping out proposal resulted",
                              "in non-finite value(s) in step", i, ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0}
               prop <- Mo1[["parm"]]
               ### Reflect at boundary
               if(y.slice > Mo1[["LP"]]) {
                    reflections <- reflections + 1
                    prop <- prop0
                    g <- partial(Model, prop, Data)
                    p <- p - 2*g*{(t(p) %*% g) / Norm(g)^2}}}
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Final proposal failed.\n",
                         file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Final proposal resulted",
                         "in non-finite value(s).\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          Mo0 <- Mo1
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrwm <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, Debug, LogFile)
     {
     Block <- Specs[["B"]]
     if(length(Block) == 0) {
          if(!is.symmetric.matrix(VarCov)) {
               cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
                    append=TRUE)
               VarCov <- as.symmetric.matrix(VarCov)}
          if(!is.positive.definite(VarCov)) {
               cat("\nNon-Positive-Definite Covar, correcting now...\n",
                    file=LogFile, append=TRUE)
               VarCov <- as.positive.definite(VarCov)}
          U <- chol(VarCov)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
               ### Propose new parameter values
               prop <- as.vector(Mo0[["parm"]] +
                    rbind(rnorm(LIV)) %*% U)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                   if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal resulted",
                              "in non-finite value(s).\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                        Mo0 <- Mo1
                        Acceptance <- Acceptance + 1}}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     else {
          B <- length(Block)
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from ",
                    "number of blocks.", file=LogFile, append=TRUE)
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from ",
                         "block length.", file=LogFile, append=TRUE)
               if(!is.symmetric.matrix(VarCov[[b]])) {
                    cat("\nAsymmetric Covar block, correcting now...\n",
                         file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.symmetric.matrix(VarCov[[b]])}
               if(!is.positive.definite(VarCov[[b]])) {
                    cat("\nNon-Positive-Definite Covar block,",
                         "correcting now...\n", file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.positive.definite(VarCov[[b]])}}
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile,
                         append=TRUE)
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]] +
                         rbind(rnorm(length(Block[[b]]))) %*%
                         chol(VarCov[[b]])
                    if({b == 1} & {iter %% Status == 0})
                         cat(",   Proposal: Blockwise,   LP: ",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in block", b,
                                   "failed.\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in block", b,
                                   "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else {
                         ### Accept/Reject
                         log.u <- log(runif(1))
                         log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                         if(!is.finite(log.alpha)) log.alpha <- 0
                         if(log.u < log.alpha) {
                              Mo0 <- Mo1
                              Acceptance <- Acceptance +
                                   length(Block[[b]]) / LIV}}
                    }
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcsamwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, Debug, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV,
          byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          else staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed for",
                              Data[["parm.names"]][j], ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal for",
                              Data[["parm.names"]][j],
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         Acceptance[j] <- Acceptance[j] + 1}}}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- tuning}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
.mcmcsgld <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     file <- Specs[["file"]]
     Nr <- Specs[["Nr"]]
     Nc <- Specs[["Nc"]]
     size <- Specs[["size"]]
     Acceptance <- matrix(Iterations, 1, LIV)
     con <- file(file, open="r")
     on.exit(close(con))
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Sample Data
          seek(con, 0)
          skip.rows <- sample.int(Nr - size, size=1)
          Data[["X"]] <- matrix(scan(file=con, sep=",", skip=skip.rows,
               nlines=size, quiet=TRUE), size, Nc, byrow=TRUE)
          ### Propose new parameter values
          g <- partial(Model, Mo0[["parm"]], Data)
          eta <- rnorm(LIV, 0, epsilon[iter])
          prop <- Mo0[["parm"]] + {epsilon[iter]/2}*g + eta
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal failed.\n", file=LogFile,
                         append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) {
               if(Debug[["DB.Model"]] == TRUE) {
                    cat("\nWARNING: Proposal resulted in non-finite",
                         "value(s).\n", file=LogFile, append=TRUE)
                    cat("  Iteration:", iter, "Proposal:\n",
                         paste("c(",paste(prop, collapse=","),")",
                         sep=""), "\n", file=LogFile, append=TRUE)}
               Mo1 <- Mo0
               }
          Mo0 <- Mo1
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcslice <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     Block <- Specs[["B"]]
     B <- length(Block)
     Bounds <- Specs[["Bounds"]]
     m <- Specs[["m"]]
     Type <- Specs[["Type"]]
     w <- Specs[["w"]]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Proceed by Block
          for (b in 1:B) {
               ### Random-Scan Componentwise Estimation
               if(Type[[b]] == "Continuous") {
                    for (j in sample(Block[[b]])) {
                         y.slice <- Mo0[["LP"]] - rexp(1)
                         u <- runif(1,0,w[[b]])
                         intL <- intR <- prop <- Mo0[["parm"]]
                         L <- intL[j] - u
                         R <- intR[j] + (w[[b]] - u)
                         ### Unlimited number of steps
                         if(is.infinite(m[[b]])) {
                              repeat {
                                   if(L <= Bounds[[b]][1]) break
                                   intL[j] <- L
                                   MoL <- try(Model(intL, Data),
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(MoL, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n",
                                                  file=LogFile, append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   else if(!is.finite(MoL[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n",
                                                  file=LogFile, append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   if(MoL[["LP"]] <= y.slice) break
                                   L <- L - w[[b]]}
                              repeat {
                                   if(R >= Bounds[[b]][2]) break
                                   intR[j] <- R
                                   MoR <- try(Model(intR, Data),
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(MoR, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the upper bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n",
                                                  file=LogFile, append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   else if(!is.finite(MoR[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the upper bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n",
                                                  file=LogFile, append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   if(MoR[["LP"]] <= y.slice) break
                                   R <- R + w[[b]]}
                              }
                         else if(m[[b]] > 1) {
                              ### Limited number of steps
                              J <- floor(runif(1,0,m[[b]]))
                              K <- (m[[b]] - 1) - J
                              while (J > 0) {
                                   if(L <= Bounds[[b]][1]) break
                                   intL[j] <- L
                                   MoL <- try(Model(intL, Data),
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(MoL, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n",
                                                  file=LogFile, append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   else if(!is.finite(MoL[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n",
                                                  file=LogFile, append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   if(MoL[["LP"]] <= y.slice) break
                                   L <- L - w[[b]]
                                   J <- J - 1}
                              while (K > 0) {
                                   if(R >= Bounds[[b]][2]) break
                                   intR[j] <- R
                                   MoR <- try(Model(intR, Data),
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(MoR, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the upper bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n",
                                                  file=LogFile, append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   else if(!is.finite(MoR[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n",
                                                  file=LogFile, append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   R <- R + w[[b]]
                                   K <- K - 1}
                              }
                         ### Shrink the interval to lower and upper bounds
                         if(L < Bounds[[b]][1]) L <- Bounds[[b]][1]
                         if(R > Bounds[[b]][2]) R <- Bounds[[b]][2]
                         ### Rejection Sampling
                         repeat {
                              L <- min(L,R)
                              R <- max(L,R)
                              prop[j] <- runif(1,L,R)
                              Mo1 <- try(Model(prop, Data),
                                   silent=!Debug[["DB.Model"]])
                              if(inherits(Mo1, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Rejection sampling",
                                             "failed for",
                                             Data[["parm.names"]][j], "\n",
                                             file=LogFile, append=TRUE)
                                        cat("  Iteration:", iter,
                                             "Current:", round(Mo0[["parm"]][j]),
                                             "Proposed:", round(prop[j],5),
                                             "\n", file=LogFile, append=TRUE)}
                                   Mo1 <- Mo0
                                   }
                              else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                   Mo1[["Monitor"]])))) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Rejection sampling for",
                                             Data[["parm.names"]][j],
                                             "resulted in non-finite value(s).\n",
                                             file=LogFile, append=TRUE)
                                        cat("  Iteration:", iter,
                                             "Current:", round(Mo0[["parm"]][j]),
                                             "Proposed:", round(prop[j],5),
                                             "\n", file=LogFile, append=TRUE)}
                                   Mo1 <- Mo0}
                              if(Mo1[["LP"]] >= y.slice) break
                              else if(abs(R-L) < 1e-100) break
                              if(Mo1[["parm"]][j] > Mo0[["parm"]][j])
                                   R <- Mo1[["parm"]][j]
                              else L <- Mo1[["parm"]][j]}
                         Mo0 <- Mo1
                         }
                    }
               else if(Type[[b]] == "Nominal") {
                    for (j in sample(Block[[b]])) {
                         y.slice <- Mo0[["LP"]] - rexp(1)
                         LP.grid <- theta <- Bounds[[b]][1]:Bounds[[b]][2]
                         for (i in 1:length(LP.grid)) {
                              prop[j] <- theta[i]
                              Mo1 <- try(Model(prop, Data),
                                   silent=!Debug[["DB.Model"]])
                              if(inherits(Mo1, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                        cat("\nWARNING: Evaluating",
                                             Data[["parm.names"]][j], "at",
                                             round(prop[j],5), "failed.\n",
                                             file=LogFile, append=TRUE)
                                   LP.grid[i] <- 0
                                   }
                              else if(!is.finite(Mo1[["LP"]])) {
                                   if(Debug[["DB.Model"]] == TRUE)
                                        cat("\nWARNING: Evaluating",
                                             Data[["parm.names"]][j], "at",
                                             round(prop[j],5), "resulted",
                                             "in non-finite value(s).\n",
                                             file=LogFile, append=TRUE)
                                   LP.grid[i] <- 0
                                   }
                              else if(Mo1[["LP"]] < y.slice)
                                   LP.grid[i] <- 0
                              else LP.grid[i] <- exp(Mo1[["LP"]])}
                         if(sum(LP.grid) > 0)
                              LP.grid <- LP.grid / sum(LP.grid)
                         else LP.grid <- rep(1/length(LP.grid),
                              length(LP.grid))
                         ### Rejection Sampling
                         repeat {
                              prop[j] <- theta[sample(1:length(LP.grid),1,
                                   prob=LP.grid)]
                              Mo1 <- try(Model(prop, Data),
                                   silent=!Debug[["DB.Model"]])
                              if(is.finite(Mo1[["LP"]])) {
                                   if(Mo1[["LP"]] >= y.slice) break}}
                         Mo0 <- Mo1
                         }
                    }
               else { ### Ordinal
                    for (j in sample(Block[[b]])) {
                         y.slice <- Mo0[["LP"]] - rexp(1)
                         intL <- intR <- prop <- Mo0[["parm"]]
                         L <- intL[j] - w[[b]]
                         R <- intR[j] + w[[b]]
                         ### Unlimited number of steps
                         if(is.infinite(m[[b]])) {
                              repeat {
                                   if(L <= Bounds[[b]][1]) break
                                   intL[j] <- L
                                   MoL <- try(Model(intL, Data),
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(MoL, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n", file=LogFile,
                                                  append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   else if(!is.finite(MoL[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n", file=LogFile,
                                                  append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   if(MoL[["LP"]] <= y.slice) break
                                   L <- L - w[[b]]}
                              repeat {
                                   if(R >= Bounds[[b]][2]) break
                                   intR[j] <- R
                                   MoR <- try(Model(intR, Data),
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(MoR, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the upper bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n", file=LogFile,
                                                  append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   else if(!is.finite(MoR[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the upper bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n", file=LogFile,
                                                  append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   if(MoR[["LP"]] <= y.slice) break
                                   R <- R + w[[b]]}
                              }
                         else if(m[[b]] > 1) {
                              ### Limited number of steps
                              J <- floor(runif(1,0,m[[b]]))
                              K <- (m[[b]] - 1) - J
                              while (J > 0) {
                                   if(L <= Bounds[[b]][1]) break
                                   intL[j] <- L
                                   MoL <- try(Model(intL, Data),
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(MoL, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n", file=LogFile,
                                                  append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   else if(!is.finite(MoL[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n", file=LogFile,
                                                  append=TRUE)
                                        L <- L + w[[b]]
                                        break}
                                   if(MoL[["LP"]] <= y.slice) break
                                   L <- L - w[[b]]
                                   J <- J - 1}
                              while (K > 0) {
                                   if(R >= Bounds[[b]][2]) break
                                   intR[j] <- R
                                   MoR <- try(Model(intR, Data),
                                       silent=!Debug[["DB.Model"]])
                                   if(inherits(MoR, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the upper bound failed for",
                                                  Data[["parm.names"]][j],
                                                  ".\n", file=LogFile,
                                                  append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   if(!is.finite(MoR[["LP"]])) {
                                        if(Debug[["DB.Model"]] == TRUE)
                                             cat("\nWARNING: Stepping out",
                                                  "the lower bound for",
                                                  Data[["parm.names"]][j],
                                                  "resulted in a non-finite",
                                                  "LP.\n", file=LogFile,
                                                  append=TRUE)
                                        R <- R - w[[b]]
                                        break}
                                   R <- R + w[[b]]
                                   K <- K - 1}
                              }
                         ### Shrink the interval to lower and upper bounds
                         if(L < Bounds[[b]][1]) L <- Bounds[[b]][1] 
                         if(R > Bounds[[b]][2]) R <- Bounds[[b]][2]
                         ### Rejection Sampling
                         repeat {
                              prop[j] <- sample(L:R,1)
                              Mo1 <- try(Model(prop, Data),
                                   silent=!Debug[["DB.Model"]])
                              if(inherits(Mo1, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Rejection sampling",
                                             "failed for",
                                             Data[["parm.names"]][j], "\n",
                                             file=LogFile, append=TRUE)
                                        cat("  Iteration:", iter,
                                             "Current:", round(Mo0[["parm"]][j]),
                                             "Proposed:", round(prop[j],5),
                                             "\n", file=LogFile, append=TRUE)}
                                   Mo1 <- Mo0
                                   }
                              else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                   Mo1[["Monitor"]])))) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Rejection sampling for",
                                             Data[["parm.names"]][j],
                                             "resulted in non-finite value(s).\n",
                                             file=LogFile, append=TRUE)
                                        cat("  Iteration:", iter,
                                             "Current:", round(Mo0[["parm"]][j]),
                                             "Proposed:", round(prop[j],5),
                                             "\n", file=LogFile, append=TRUE)}
                                   Mo1 <- Mo0}
                              if(Mo1[["LP"]] >= y.slice) break
                              else if(abs(R-L) < 1e-100) break
                              if(Mo1[["parm"]][j] > Mo0[["parm"]][j])
                                   R <- Mo1[["parm"]][j]
                              else L <- Mo1[["parm"]][j]}
                         Mo0 <- Mo1
                         }
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=.colVars(thinned))
     return(out)
     }
.mcmcsmwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, Debug, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          else staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed for",
                              Data[["parm.names"]][j], ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal for",
                              Data[["parm.names"]][j],
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         Acceptance[j] <- Acceptance[j] + 1}}}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
.mcmcthmc <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     m <- Specs[["m"]]
     invm <- as.inverse(m)
     U <- chol(m)
     Temperature <- Specs[["Temperature"]]
     gr <- partial(Model, Mo0[["parm"]], Data)
     sqrt.Temp <- sqrt(Temperature)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum1 <- momentum0 <- as.vector(rnorm(LIV) %*% U)
          kinetic0 <- t(momentum0) %*% invm %*% momentum0 / 2
          Mo0.1 <- Mo0
          for (l in 1:L) {
               if(2*(l-1) < L) momentum1 <- momentum1 * sqrt.Temp
               else momentum1 <- momentum1 / sqrt.Temp
               momentum1 <- momentum1 + (epsilon/2) * gr
               prop <- prop + as.vector(epsilon %*% invm) * momentum1
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed in leapfrog", l,
                              ".\n", file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal in leapfrog", l,
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(prop, collapse=","),")",
                              sep=""), "\n", file=LogFile, append=TRUE)}
                    Mo1 <- Mo0.1}
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal failed in leapfrog",
                                   l, ".\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0.1
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal in leapfrog",
                                   l, "resulted in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0.1}}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr <- partial(Model, prop, Data)
               momentum1 <- momentum1 + (epsilon/2) * gr
               if(2*l > L) momentum1 <- momentum1 / sqrt.Temp
               else momentum1 <- momentum1 * sqrt.Temp}
          momentum1 <- -momentum1
          kinetic1 <- t(momentum1) %*% invm %*% momentum1 / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               kinetic0 <- kinetic1
               Acceptance <- Acceptance + 1}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=matrix(epsilon, 1, LIV),
          Mon=Mon,
          thinned=thinned,
          VarCov=cov(thinned))
     return(out)
     }
.mcmctwalk <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
     LogFile)
     {
     xp0 <- SIV <- Specs[["SIV"]]
     n1 <- Specs[["n1"]]
     at <- Specs[["at"]]
     aw <- Specs[["aw"]]
     IntProd <- function(x) {return(sum(x*x))}
     DotProd <- function(x, y) {return(sum(x*y))}
     Simh1 <- function(dim, pphi, x, xp, beta)
          {
          phi <- runif(dim) < pphi
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + beta*(xp[i] - x[i]))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi)))
          }
     Simfbeta <- function(at)
          {
          if(runif(1) < (at-1)/(2*at))
               return(exp(1/(at + 1)*log(runif(1))))
          else
               return(exp(1/(1 - at)*log(runif(1))))
          }
     Simh2 <- function(dim, pphi, aw, x, xp)
          {
          u <- runif(dim)
          phi <- runif(dim) < pphi
          z <- (aw/(1+aw))*(aw*u^2 + 2*u -1)
          z <- z*phi
          return(list(rt=x + (x - xp)*z, nphi=sum(phi)))
          }
     Simh3 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))
          x + sigma*rnorm(dim)*phi
          return(list(rt=x + sigma*rnorm(dim)*phi, nphi=sum(phi),
               sigma=sigma))
          }
     G3U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd(h - xp)/(sigma^2))
          else
               return(0) 
          }
     Simh4 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))/3
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + sigma*rnorm(1))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi), sigma=sigma))
          }
     G4U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd((h - x))/(sigma^2))
          else
               return(0)
          }
     OneMove <- function(dim, Model, Data, x, U, xp, Up, at=at, aw=aw,
          pphi=pphi, F1=0.4918, F2=0.9836, F3=0.9918, Mo0.1, Mo0.2)
          {
          dir <- runif(1) ### Determine which set of points
          ker <- runif(1) ### Choose a kernel
          if(ker < F1) {
               ### Kernel h1: Traverse
               funh <- 1
               if(dir < 0.5) {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, xp, x, beta)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["LP"]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propUp <- Mo1.2[["LP"]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, x, xp, beta)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propU <- Mo1.1[["LP"]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F2) {
               ### Kernel h2: Walk
               funh <- 2
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh2(dim, pphi, aw, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh2(dim, pphi, aw, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propU <- Mo1.1[["LP"]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F3) {
               ### Kernel h3: Blow
               funh <- 3
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh3(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         W1 <- G3U(nphi, sigma,  yp, xp,  x)
                         W2 <- G3U(nphi, sigma,  xp, yp,  x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh3(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[["LP"]] * -1
                         W1 <- G3U(nphi, sigma, y, x, xp)
                         W2 <- G3U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else {
               ## Kernel h4: Hop
               funh <- 4
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh4(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         W1 <- G4U(nphi, sigma, yp, xp, x)
                         W2 <- G4U(nphi, sigma, xp, yp, x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh4(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data),
                         silent=!Debug[["DB.Model"]])
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[["LP"]] * -1
                         W1 <- G4U(nphi, sigma, y, x, xp)
                         W2 <- G4U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if(check1 & check2 & is.finite(A) & (dir < 0.5))
               Mo0.2 <- Mo1.2
          else if(check1 & check2 & is.finite(A) & (dir >= 0.5))
               Mo0.1 <- Mo1.1
          else if(!is.finite(A)) A <- 0
          return(list(y=y, propU=propU, yp=yp, propUp=propUp, A=A,
               funh=funh, nphi=nphi, Mo0.1=Mo0.1, Mo0.2=Mo0.2))
          }
     Runtwalk <- function(Iterations, dim, x0, xp0, pphi, at, aw,
          F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, Model, Data, Status,
          Thinning, Acceptance, Dev, Mon, Mo0, thinned, Debug, LogFile)
          {
          x <- x0 ### Primary vector of initial values
          xp <- xp0 ### Secondary vector of initial values
          Mo0.1 <- try(Model(x, Data), silent=!Debug[["DB.Model"]])
          Mo0.2 <- try(Model(xp, Data), silent=!Debug[["DB.Model"]])
          if(inherits(Mo0.1, "try-error") | inherits(Mo0.2, "try-error"))
               stop("Error in estimating the log-posterior.",
                    file=LogFile, append=TRUE)
          if(any(!is.finite(c(Mo0.1[["LP"]], Mo0.2[["LP"]]))))
               stop("The log-posterior is non-finite.", file=LogFile,
                    append=TRUE)
          if(identical(x, as.vector(Mo0.1[["parm"]])) &
               identical(xp, as.vector(Mo0.2[["parm"]]))) {
               U <- Mo0.1[["LP"]] * -1
               Up <- Mo0.2[["LP"]] * -1}
          else {
               cat("\nInitial values are out of support.", file=LogFile,
                    append=TRUE)
               cat("\n  Initial.Values=", x, file=LogFile, append=TRUE)
               cat("\n SIV=", xp, file=LogFile, append=TRUE)
               stop("Try re-specifying initial values.", file=LogFile,
                    append=TRUE)}
          if(any(abs(x - xp) <= 0))
               stop("\nBoth vectors of initial values are not unique.",
                    file=LogFile, append=TRUE)
          Acceptance <- 0
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate Subset,   LP: ",
                         round(Mo0.1[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Assign x and xp
               x <- as.vector(Mo0.1[["parm"]])
               xp <- as.vector(Mo0.2[["parm"]])
               ### Propose New Parameter Values
               move <- OneMove(dim=dim, Model, Data, x, U, xp, Up,
                    at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3,
                    Mo0.1=Mo0.1, Mo0.2=Mo0.2)
               ### Accept/Reject
               if(runif(1) < move$A) {
                    Mo0.1 <- move$Mo0.1
                    Mo0.2 <- move$Mo0.2
                    Acceptance <- Acceptance + 1 #move$nphi/dim
                    x <- move$y
                    U <- move$propU
                    xp <- move$yp
                    Up <- move$propUp
                    }
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0.1[["parm"]]
                    Dev[t.iter] <- Mo0.1[["Dev"]]
                    Mon[t.iter,] <- Mo0.1[["Monitor"]]}
               }
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=.colVars(thinned))
          return(out)
          }
     out <- Runtwalk(Iterations=Iterations, dim=LIV, x0=Mo0[["parm"]],
          xp0=xp0, pphi=min(LIV, n1)/LIV, at=6, aw=1.5, Model=Model,
          Data=Data, Status=Status, Thinning=Thinning,
          Acceptance=Acceptance, Dev=Dev, Mon=Mon, Mo0=Mo0,
          thinned=thinned, Debug=Debug, LogFile=LogFile)
     ### Output
     return(out)
     }
.mcmcuess <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     VarCov, Debug, LogFile)
     {
     A <- Specs[["A"]]
     Block <- Specs[["B"]]
     m <- Specs[["m"]]
     n <- Specs[["n"]]
     w <- 0.05
     B <- length(Block)
     if(B == 0) {
          if(!is.symmetric.matrix(VarCov)) {
               cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
                    append=TRUE)
               VarCov <- as.symmetric.matrix(VarCov)}
          if(!is.positive.definite(VarCov)) {
               cat("\nNon-Positive-Definite Covar, correcting now...\n",
                    file=LogFile, append=TRUE)
               VarCov <- as.positive.definite(VarCov)}
          decomp.freq <- max(LIV * floor(Iterations / Thinning / 100), 10)
          cat("\nEigendecomposition will occur every", decomp.freq,
               "iterations.\n\n", file=LogFile, append=TRUE)
          S.eig <-try(eigen(VarCov), silent=!Debug[["DB.eigen"]])
          if(inherits(S.eig, "try-error")) S.eig <- NULL
          obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
          obs.scatter <- tcrossprod(Mo0[["parm"]])*n
          DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning)+1,
               LIV, byrow=TRUE)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Eigenvectors of the Sample Covariance Matrix
               if({iter %% decomp.freq == 0} & {iter > 1} & {iter < A}) {
                    VarCov <- obs.scatter/{n + iter} -
                         tcrossprod(obs.sum/{n + iter})
                    S.eig <- eigen(VarCov)}
               ### Non-Adaptive or Adaptive
               if(runif(1) < w || is.null(S.eig)) {
                    v <- rnorm(LIV)
                    v <- v / sqrt(sum(v*v))
                    }
               else {
                    which.eig <- floor(1 + LIV * runif(1))
                    v <- S.eig$vectors[,which.eig] *
                         sqrt(abs(S.eig$values[which.eig]))}
               ### Slice Interval
               Mo0.1 <- try(Model(Mo0[["parm"]], Data),
                    silent=!Debug[["DB.Model"]])
               if(inherits(Mo0.1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed.\n", file=LogFile,
                              append=TRUE)
                         cat("  Iteration:", iter, "Proposal:\n",
                              paste("c(",paste(Mo0[["parm"]],
                              collapse=","),")",sep=""), "\n",
                              file=LogFile, append=TRUE)}
                    Mo0.1 <- Mo0}
               Mo0 <- Mo0.1
               y.slice <- Mo0[["LP"]] - rexp(1)
               L <- -runif(1)
               U <- L + 1
               if(m > 0) {
                    L.y <- try(Model(Mo0[["parm"]] + v*L, Data)[["LP"]],
                         silent=!Debug[["DB.Model"]])
                    if(inherits(L.y, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Stepping out the lower",
                                   "bound failed.\n", file=LogFile,
                                   append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(Mo0[["parm"]] + v*L,
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         L.y <- Mo0[["LP"]]
                         }
                    else if(!is.finite(L.y)) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Stepping out the lower",
                                   "bound resulted in non-finite LP.\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(Mo0[["parm"]] + v*L,
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         L.y <- Mo0[["LP"]]}
                    U.y <- try(Model(Mo0[["parm"]] + v*U, Data)[["LP"]],
                         silent=!Debug[["DB.Model"]])
                    if(inherits(U.y, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Stepping out the upper",
                                   "bound failed.\n", file=LogFile,
                                   append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(Mo0[["parm"]] + v*U,
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         U.y <- Mo0[["LP"]]
                         }
                    else if(!is.finite(U.y)) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Stepping out the upper",
                                   "bound resulted in non-finite LP.\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(Mo0[["parm"]] + v*U,
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         U.y <- Mo0[["LP"]]}
                    step <- 0
                    while({L.y > y.slice || U.y > y.slice} && step < m) {
                         step <- step + 1
                         if(runif(1) < 0.5) {
                              L <- L - 1
                              L.y <- try(Model(Mo0[["parm"]] + v*L,
                                   Data)[["LP"]], silent=!Debug[["DB.Model"]])
                              if(inherits(L.y, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Stepping out the lower",
                                             "bound failed.\n", file=LogFile,
                                             append=TRUE)
                                        cat("  Iteration:", iter, "Proposal:\n",
                                             paste("c(",paste(Mo0[["parm"]] +
                                             v*L, collapse=","),")",sep=""),
                                             "\n", file=LogFile, append=TRUE)}
                                   L.y <- Mo0[["LP"]]
                                   }
                              else if(!is.finite(L.y)) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Stepping out the lower",
                                             "bound resulted in non-finite LP.\n",
                                             file=LogFile, append=TRUE)
                                        cat("  Iteration:", iter, "Proposal:\n",
                                             paste("c(",paste(Mo0[["parm"]] +
                                             v*L, collapse=","),")",sep=""),
                                             "\n", file=LogFile, append=TRUE)}
                                   L.y <- Mo0[["LP"]]
                                   }
                              }
                         else {
                              U <- U + 1
                              U.y <- try(Model(Mo0[["parm"]] + v*U,
                                   Data)[["LP"]], silent=!Debug[["DB.Model"]])
                              if(inherits(U.y, "try-error")) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Stepping out the upper",
                                             "bound failed.\n", file=LogFile,
                                             append=TRUE)
                                        cat("  Iteration:", iter, "Proposal:\n",
                                             paste("c(",paste(Mo0[["parm"]] +
                                             v*U, collapse=","),")",sep=""),
                                             "\n", file=LogFile, append=TRUE)}
                                   U.y <- Mo0[["LP"]]
                                   }
                              else if(!is.finite(U.y)) {
                                   if(Debug[["DB.Model"]] == TRUE) {
                                        cat("\nWARNING: Stepping out the upper",
                                             "bound resulted in non-finite LP.\n",
                                             file=LogFile, append=TRUE)
                                        cat("  Iteration:", iter, "Proposal:\n",
                                             paste("c(",paste(Mo0[["parm"]] +
                                             v*U, collapse=","),")",sep=""),
                                             "\n", file=LogFile, append=TRUE)}
                                   U.y <- Mo0[["LP"]]}}}}
               ### Rejection Sampling
               repeat {
                    prop.offset <- runif(1, min=L, max=U)
                    prop <- Mo0[["parm"]] + prop.offset * v
                    Mo1 <- try(Model(prop, Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Rejection sampling failed.\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0
                         }
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]])))) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Rejection sampling resulted",
                                   "in non-finite value(s).\n",
                                   file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(prop, collapse=","),")",
                                   sep=""), "\n", file=LogFile, append=TRUE)}
                         Mo1 <- Mo0}
                    prop <- Mo1[["parm"]]
                    if(Mo1[["LP"]] >= y.slice) break
                    else if(abs(prop.offset < 1e-100)) {
                         Mo1 <- Mo0
                         break}
                    if(prop.offset < 0) L <- prop.offset
                    else U <- prop.offset}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]
                    DiagCovar[t.iter,] <- diag(S.eig$vectors)}
               obs.sum <- obs.sum + Mo1[["parm"]]
               obs.scatter <- obs.scatter + tcrossprod(Mo1[["parm"]])
               Mo0 <- Mo1}
          }
     else {
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from ",
                    "number of blocks.", file=LogFile, append=TRUE)
          S.eig <- obs.sum <- obs.scatter <- list()
          decomp.freq <- rep(0, length(B))
          DiagCovar <- rep(0, LIV)
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from block ",
                         "length.", file=LogFile, append=TRUE)
               if(!is.symmetric.matrix(VarCov[[b]])) {
                    cat("\nAsymmetric Covar block, correcting now...\n",
                         file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.symmetric.matrix(VarCov[[b]])}
               if(!is.positive.definite(VarCov[[b]])) {
                    cat("\nNon-Positive-Definite Covar block,",
                         "correcting now...\n", file=LogFile, append=TRUE)
                    VarCov[[b]] <- as.positive.definite(VarCov[[b]])}
               decomp.freq[b] <- max(length(Block[[b]]) *
                    floor(Iterations / Thinning / 100), 10)
               S.eig[[b]] <-try(eigen(VarCov[[b]]),
                    silent=!Debug[["DB.eigen"]])
               if(inherits(S.eig[[b]], "try-error")) S.eig[[b]] <- NULL
               obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n,
                    length(Block[[b]]), 1)
               obs.scatter[[b]] <- tcrossprod(Mo0[["parm"]][Block[[b]]])*n
               DiagCovar[Block[[b]]] <- diag(VarCov[[b]])}
          if(all(decomp.freq == decomp.freq[1]))
               cat("\nEigendecomposition will occur every", decomp.freq[1],
                    "iterations.\n\n", file=LogFile, append=TRUE)
          else cat("\nEigendecomposition frequency varies by block,",
                    "and will occur between\n",
                    min(decomp.freq), "and", max(decomp.freq),
                    "iterations.\n\n", file=LogFile, append=TRUE)
          DiagCovar <- matrix(DiagCovar, floor(Iterations / Thinning)+1,
               LIV, byrow=TRUE)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Blockwise,   LP: ",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Proceed by Block
               for (b in 1:B) {
                    ### Eigenvectors of the Sample Covariance Matrix
                    if({iter %% decomp.freq[b] == 0} & {iter > 1} &
                         {iter < A}) {
                         VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
                              tcrossprod(obs.sum[[b]]/{n + iter})
                         S.eig[[b]] <- eigen(VarCov[[b]])}
                    ### Non-Adaptive or Adaptive
                    if(runif(1) < w || is.null(S.eig[[b]])) {
                         v <- rnorm(length(Block[[b]]))
                         v <- v / sqrt(sum(v*v))
                         }
                    else {
                         which.eig <- floor(1 + length(Block[[b]]) * runif(1))
                         v <- S.eig[[b]]$vectors[,which.eig] *
                              sqrt(abs(S.eig[[b]]$values[which.eig]))}
                    ### Slice Interval
                    Mo0.1 <- try(Model(Mo0[["parm"]][Block[[b]]], Data),
                         silent=!Debug[["DB.Model"]])
                    if(inherits(Mo0.1, "try-error")) {
                         if(Debug[["DB.Model"]] == TRUE) {
                              cat("\nWARNING: Proposal for block", b,
                                   "failed.\n", file=LogFile, append=TRUE)
                              cat("  Iteration:", iter, "Proposal:\n",
                                   paste("c(",paste(Mo0[["parm"]][Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                                   file=LogFile, append=TRUE)}
                         Mo0.1 <- Mo0}
                    Mo0 <- Mo0.1
                    y.slice <- Mo0[["LP"]] - rexp(1)
                    L <- -runif(1)
                    U <- L + 1
                    if(m > 0) {
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- prop[Block[[b]]] + v*L
                         L.y <- try(Model(prop, Data)[["LP"]],
                              silent=!Debug[["DB.Model"]])
                         if(inherits(L.y, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Stepping out the lower",
                                        "bound failed for block", b, ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              L.y <- Mo0[["LP"]]
                              }
                         else if(!is.finite(L.y)) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Stepping out the lower",
                                        "bound resulted in non-finite LP",
                                        "for block", b, ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              L.y <- Mo0[["LP"]]}
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- prop[Block[[b]]] + v*U
                         U.y <- try(Model(prop, Data)[["LP"]],
                              silent=!Debug[["DB.Model"]])
                         if(inherits(U.y, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Stepping out the upper",
                                        "bound failed for block", b, ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              U.y <- Mo0[["LP"]]
                              }
                         else if(!is.finite(U.y)) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Stepping out the upper",
                                        "bound resulted in non-finite LP",
                                        "for block", b, ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              U.y <- Mo0[["LP"]]}
                         step <- 0
                         while({L.y > y.slice || U.y > y.slice} && step < m) {
                              step <- step + 1
                              if(runif(1) < 0.5) {
                                   L <- L - 1
                                   prop <- Mo0[["parm"]]
                                   prop[Block[[b]]] <- prop[Block[[b]]] + v*L
                                   L.y <- try(Model(prop, Data)[["LP"]],
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(L.y, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE) {
                                             cat("\nWARNING: Stepping out the",
                                                  "lower bound failed for",
                                                  "block", b, ".\n",
                                                  file=LogFile, append=TRUE)
                                             cat("  Iteration:", iter,
                                                  "Proposal:\n", paste("c(",
                                                  paste(prop[Block[[b]]],
                                                  collapse=","),")",sep=""),
                                                  "\n", file=LogFile,
                                                  append=TRUE)}
                                        L.y <- Mo0[["LP"]]
                                        }
                                   else if(!is.finite(L.y)) {
                                        if(Debug[["DB.Model"]] == TRUE) {
                                             cat("\nWARNING: Stepping out the",
                                                  "lower bound resulted in ",
                                                  "non-finite LP for block",
                                                  b, ".\n", file=LogFile,
                                                  append=TRUE)
                                             cat("  Iteration:", iter,
                                                  "Proposal:\n", paste("c(",
                                                  paste(prop[Block[[b]]],
                                                  collapse=","),")",sep=""),
                                                  "\n", file=LogFile,
                                                  append=TRUE)}
                                        L.y <- Mo0[["LP"]]
                                        }
                                   }
                              else {
                                   U <- U + 1
                                   prop <- Mo0[["parm"]]
                                   prop[Block[[b]]] <- prop[Block[[b]]] + v*U
                                   U.y <- try(Model(prop, Data)[["LP"]],
                                        silent=!Debug[["DB.Model"]])
                                   if(inherits(U.y, "try-error")) {
                                        if(Debug[["DB.Model"]] == TRUE) {
                                             cat("\nWARNING: Stepping out the",
                                                  "upper bound failed for",
                                                  "block", b, ".\n",
                                                  file=LogFile, append=TRUE)
                                             cat("  Iteration:", iter,
                                                  "Proposal:\n", paste("c(",
                                                  paste(prop[Block[[b]]],
                                                  collapse=","),")",sep=""),
                                                  "\n", file=LogFile,
                                                  append=TRUE)}
                                        U.y <- Mo0[["LP"]]
                                        }
                                   else if(!is.finite(U.y)) {
                                        if(Debug[["DB.Model"]] == TRUE) {
                                             cat("\nWARNING: Stepping out the",
                                                  "upper bound resulted in",
                                                  "non-finite LP for block",
                                                  b, ".\n", file=LogFile,
                                                  append=TRUE)
                                             cat("  Iteration:", iter,
                                                  "Proposal:\n", paste("c(",
                                                  paste(prop[Block[[b]]],
                                                  collapse=","),")",sep=""),
                                                  "\n", file=LogFile,
                                                  append=TRUE)}
                                        U.y <- Mo0[["LP"]]}}}}
                    ### Rejection Sampling
                    repeat {
                         prop.offset <- runif(1, min=L, max=U)
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- prop[Block[[b]]] + prop.offset*v
                         Mo1 <- try(Model(prop, Data),
                              silent=!Debug[["DB.Model"]])
                         if(inherits(Mo1, "try-error")) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Rejection sampling",
                                        "failed for block", b, ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0
                              }
                         else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                              Mo1[["Monitor"]])))) {
                              if(Debug[["DB.Model"]] == TRUE) {
                                   cat("\nWARNING: Rejection sampling",
                                        "resulted in non-finite",
                                        "value(s) for block", b, ".\n",
                                        file=LogFile, append=TRUE)
                                   cat("  Iteration:", iter, "Proposal:\n",
                                        paste("c(",paste(prop[Block[[b]]],
                                        collapse=","),")",sep=""), "\n",
                                        file=LogFile, append=TRUE)}
                              Mo1 <- Mo0}
                         prop <- Mo1[["parm"]]
                         if(Mo1[["LP"]] >= y.slice) break
                         else if(abs(prop.offset < 1e-100)) {
                              Mo1 <- Mo0
                              break}
                         if(prop.offset < 0) L <- prop.offset
                         else U <- prop.offset}
                    ### Save Thinned Samples
                    if(iter %% Thinning == 0) {
                         t.iter <- floor(iter / Thinning) + 1
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]
                         DiagCovar[t.iter,Block[[b]]] <- diag(S.eig[[b]]$vectors)}
                    obs.sum[[b]] <- obs.sum[[b]] + Mo1[["parm"]][Block[[b]]]
                    obs.scatter[[b]] <- obs.scatter[[b]] +
                         tcrossprod(Mo1[["parm"]][Block[[b]]])
                    Mo0 <- Mo1}
               }
          }
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
.mcmcusamwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, Debug, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Periodicity <- Specs[["Periodicity"]]
     Fit <- Specs[["Fit"]]
     Begin <- Specs[["Begin"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.",
               file=LogFile, append=TRUE)
     ivs <- Mo0[["parm"]]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV,
          byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed for",
                              Data[["parm.names"]][j], ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal for",
                              Data[["parm.names"]][j],
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         post[iter,] <- Mo0[["parm"]]
                         Acceptance[j] <- Acceptance[j] + 1}}}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- tuning}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
.mcmcusmwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, Debug, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Fit <- Specs[["Fit"]]
     Begin <- Specs[["Begin"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.",
               file=LogFile, append=TRUE)
     ivs <- Mo0[["parm"]]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP: ",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
               if(inherits(Mo1, "try-error")) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal failed for",
                              Data[["parm.names"]][j], ".\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]])))) {
                    if(Debug[["DB.Model"]] == TRUE) {
                         cat("\nWARNING: Proposal for",
                              Data[["parm.names"]][j],
                              "resulted in non-finite value(s).\n",
                              file=LogFile, append=TRUE)
                         cat("  Iteration:", iter,
                              "Current:", round(Mo0[["parm"]][j]),
                              "Proposed:", round(prop[j],5),
                              file=LogFile, append=TRUE)}
                    Mo1 <- Mo0
                    }
               else {
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         post[iter,] <- Mo0[["parm"]]
                         Acceptance[j] <- Acceptance[j] + 1}}}
           ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }

#End
