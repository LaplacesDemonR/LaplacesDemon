###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of the LaplacesDemon function is to use MCMC on the         #
# logarithm of the unnormalized joint posterior density of a Bayesian     #
# model.                                                                  #
###########################################################################

LaplacesDemon <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=10000, Status=100, Thinning=10, Algorithm="MWG",
     Specs=NULL, LogFile="", ...)
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
     if(is.null(Data$mon.names))
          stop("In Data, mon.names is NULL.", file=LogFile, append=TRUE)
     if(is.null(Data$parm.names))
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
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(!identical(length(Initial.Values), length(Data$parm.names))) {
          cat("WARNING: The length of Initial Values differed from",
               "Data$parm.names.\n", file=LogFile, append=TRUE)
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(any(!is.finite(Initial.Values))) {
          cat("WARNING: Initial Values contain non-finite values.\n",
               file=LogFile, append=TRUE)
          Initial.Values <- rep(0, length(Data$parm.names))}
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
     if(Algorithm %in% c("ADMG","AGG","AHMC","AIES","AM","AMM","AMWG",
          "CHARM","DEMC","DRAM","DRM","ESS","Experimental","GG","Gibbs",
          "HARM","HMC","HMCDA","IM","INCA","MALA","MCMCMC","MTM","MWG",
          "NUTS","OHSS","RAM","Refractive","RDMH","RJ","RSS","RWM","SAMWG",
          "SGLD","Slice","SMWG","THMC","twalk","UESS","USAMWG","USMWG")) {
          if(Algorithm == "ADMG") {
               Algorithm <- "Adaptive Directional Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), "Periodicity"))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "AGG") {
               Algorithm <- "Adaptive Griddy-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("Grid","dparm","smax","CPUs","Packages","Dyn.libs")))
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
                         length(Initial.Values)), L=2, Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("epsilon","L","Periodicity")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
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
               }
          else if(Algorithm == "AIES") {
               Algorithm <- "Affine-Invariant Ensemble Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                          append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("Nc","Z","beta","CPUs","Packages","Dyn.libs")))
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
               if(!identical(names(Specs), c("Adaptive","Periodicity")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "AMM") {
               Algorithm <- "Adaptive-Mixture Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2), B=NULL,
                         Periodicity=1, w=0.05)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("Adaptive","B","Periodicity","w")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.null(Specs[["B"]])) {
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Specs[["B"]])) {
                              Covar[[b]] <- diag(length(Specs[["B"]][[b]]))}}}
               Specs[["w"]] <- abs(Specs[["w"]])
               if(Specs[["w"]] <= 0 || Specs[["w"]] >= 1) {
                    Specs[["w"]] <- 0.05
                    cat("\nw was misspecified and changed to 0.05.\n",
                         file=LogFile, append=TRUE)}
               }
          else if(Algorithm == "AMWG") {
               Algorithm <- "Adaptive Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), "Periodicity"))
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
                    if(!identical(names(Specs), "alpha.star"))
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
               if(!identical(names(Specs), c("Nc","Z","gamma","w")))
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
               if(!identical(names(Specs), c("Adaptive","Periodicity")))
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
               if(!identical(names(Specs), c("B")))
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
               if(!identical(names(Specs),
                    c("Grid","dparm","CPUs","Packages","Dyn.libs")))
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
                    if(!identical(names(Specs), c("FC","MWG")))
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
                    #if(!is.null(MWG) & !is.vector(MWG) & !is.numeric(MWG))
                    #     stop("MWG must be a numeric vector.",
                    #          file=LogFile, append=TRUE)
                    }
               
               }
          else if(Algorithm == "HARM") {
               Algorithm <- "Hit-And-Run Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(alpha.star=NA, B=NULL)
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!identical(names(Specs), c("alpha.star","B")))
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
                         length(Initial.Values)), L=2)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), c("epsilon","L")))
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
               }
          else if(Algorithm == "HMCDA") {
               Algorithm <- "Hamiltonian Monte Carlo with Dual-Averaging"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("A","delta","epsilon","Lmax","lambda")))
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
               if(!identical(names(Specs), "mu"))
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
               if(!identical(names(Specs), c("Adaptive","Periodicity")))
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
               if(!identical(names(Specs),
                    c("A","alpha.star","gamma","delta","epsilon")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- abs(Specs[["A"]][1])
               Specs[["gamma"]] <- min(max(Specs[["gamma"]][1], 1),
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
               if(!identical(names(Specs),
                    c("lambda","CPUs","Packages","Dyn.libs")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["lambda"]] <- abs(Specs[["lambda"]])
               if(Specs[["CPUs"]] <= 1)
                    cat("\nCPUs must be at least 2. Attempting 2 CPUs...\n")
               Specs[["CPUs"]] <- max(2, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "MTM") {
               Algorithm <- "Multiple-Try Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(K=4, CPUs=1, Packages=NULL, Dyn.libs=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("K","CPUs","Packages","Dyn.libs")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["K"]] <- abs(round(Specs[["K"]]))
               if(Specs[["CPUs"]] < 1)
                    cat("\nCPUs must be at least 1.\n")
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "MWG") {
               Algorithm <- "Metropolis-within-Gibbs"
               Specs <- NULL
               }
          else if(Algorithm == "NUTS") {
               Algorithm <- "No-U-Turn Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), c("A","delta","epsilon")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- max(min(round(abs(Specs[["A"]])),
                    Iterations),1)
               Specs[["delta"]] <- max(min(abs(Specs[["delta"]]),
                    1), 1/Iterations)
               if(!is.null(Specs[["epsilon"]]))
                    Specs[["epsilon"]] <- abs(Specs[["epsilon"]][1])
               }
          else if(Algorithm == "OHSS") {
               Algorithm <- "Oblique Hyperrectangle Slice Sampler"
               Specs <- NULL
               }
          else if(Algorithm == "RAM") {
               Algorithm <- "Robust Adaptive Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                         Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("alpha.star","Dist","gamma","Periodicity")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["alpha.star"]] <- Specs[["alpha.star"]][1]
               if(Specs[["alpha.star"]] <= 0 ||
                    Specs[["alpha.star"]] >= 1) {
                    cat("\nalpha.star not in (0,1). Changed to 0.234.\n",
                         file=LogFile, append=TRUE)
                    Specs[["alpha.star"]] <- 0.234}
               if(Specs[["Dist"]] != "t" & Specs[["Dist"]] != "N") {
                    cat("\nDist was not t or N, and changed to N.\n",
                         file=LogFile, append=TRUE)
                    Specs[["Dist"]] <- "N"}
               Specs[["gamma"]] <- Specs[["gamma"]][1]
               if(Specs[["gamma"]] <= 0.5 || Specs[["gamma"]] > 1) {
                    cat("\ngamma not in (0.5,1]. Changed to 0.66.\n",
                         file=LogFile, append=TRUE)
                    Specs[["gamma"]] <- 0.66}
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
               if(!identical(names(Specs), c("Adaptive","m","w","r")))
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
               if(!identical(names(Specs),
                    c("bin.n","bin.p","parm.p","selectable","selected")))
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
               if(!identical(names(Specs), c("m","w")))
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
                    Specs[["w"]] <- ifelse(Specs[["w"]] <= 0, 1,
                         Specs[["w"]])}
               }
          else if(Algorithm == "RWM") {
               Algorithm <- "Random-Walk Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(B=list())
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), "B"))
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
               if(!identical(names(Specs), c("Dyn","Periodicity")))
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
               if(!identical(names(Specs),
                    c("epsilon","file","Nr","Nc","size")))
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
                    Specs <- list(m=Inf, w=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), c("m","w")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               if(length(Specs[["m"]]) == 1)
                    Specs[["m"]] <- rep(Specs[["m"]], length(Initial.Values))
               else if(length(Specs[["m"]]) != length(Initial.Values)) {
                    cat("\nm was misspecified, and is replaced with Inf.\n",
                         file=LogFile, append=TRUE)
                    m <- rep(Inf, length(Initial.Values))}
               if(any(Specs[["m"]] < 1)) {
                    cat("\nm was misspecified, and is replaced with 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["m"]] <- ifelse(Specs[["m"]] < 1, 1,
                         Specs[["m"]])}
               Specs[["w"]] <- abs(Specs[["w"]])
               if(length(Specs[["w"]]) == 1)
                    Specs[["w"]] <- rep(Specs[["w"]],
                         length(Initial.Values))
               else if(length(Specs[["w"]]) != length(Initial.Values)) {
                    cat("\nw was misspecified, and is replaced with 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["w"]] <- rep(1, length(Initial.Values))}
               if(any(Specs[["w"]] <= 0)) {
                    cat("\nw was misspecified, and is replaced with 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["w"]] <- ifelse(Specs[["w"]] <= 0, 1,
                         Specs[["w"]])}
               }
          else if(Algorithm == "SMWG") {
               Algorithm <- "Sequential Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), "Dyn"))
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
               if(!identical(names(Specs), c("epsilon","L","Temperature")))
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
               if(!identical(names(Specs), c("SIV","n1","at","aw")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["SIV"]])) {
                    cat("\nGenerating SIV...\n", file=LogFile, append=TRUE)
                    if(!is.null(Data$PGF))
                         Specs[["SIV"]] <- GIV(Model, Data, PGF=TRUE)
                    else Specs[["SIV"]] <- GIV(Model, Data)}
               if(!identical(length(Specs[["SIV"]]),
                    length(Initial.Values))) {
                    cat("\nGenerating SIV due to length mismatch.\n",
                         file=LogFile, append=TRUE)
                    if(!is.null(Data$PGF))
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
                    Specs=list(m=100)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs), "m"))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               }
          else if(Algorithm == "USAMWG") {
               Algorithm <- "Updating Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!identical(names(Specs),
                    c("Dyn","Periodicity","Fit","Begin")))
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
               if(!identical(names(Specs), c("Dyn","Fit","Begin")))
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
     if(!identical(length(Mo0[["Monitor"]]), length(Data$mon.names)))
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
          cat("     increase if apply functions are 'vectorized'.\n",
               file=LogFile, append=TRUE)}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of for loops\n",
               file=LogFile, append=TRUE)
          cat("     were found in the Model specification. Iteration speed will\n",
               file=LogFile, append=TRUE)
          cat("     increase if for loops are 'vectorized'.\n",
               file=LogFile, append=TRUE)}
     #########################  Initial Settings  #########################
     Acceptance <- 0
     if(!is.finite(Mo0[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n",
               file=LogFile, append=TRUE)
          if(!is.null(Data$PGF))
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
     if(!is.null(Data$n)) if(length(Data$n) == 1) N <- Data$n
     if(!is.null(Data$N)) if(length(Data$N) == 1) N <- Data$N
     if(!is.null(Data$y)) N <- nrow(matrix(Data$y))
     if(!is.null(Data$Y)) N <- nrow(matrix(Data$Y))
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
               tuning <- sqrt(diag(Covar)); VarCov <- Covar}
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != LIV)
                    tuning <- rep(ScaleF, LIV)
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning
               }
          else if(is.null(Covar)) {
               tuning <- rep(ScaleF, LIV)
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning}
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
          "Elliptical Slice Sampler",
          "Independence Metropolis",
          "Metropolis-Adjusted Langevin Algorithm",
          "Oblique Hyperrectangle Slice Sampler",
          "Robust Adaptive Metropolis",
          "Univariate Eigenvector Slice Sampler")) {
          ### Algorithms that require VarCov, but not tuning
          if(is.list(Covar)) VarCov <- Covar
          else if(is.matrix(Covar) & !is.list(Covar)) VarCov <- Covar
          else if(is.vector(Covar) & !is.list(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- abs(as.vector(Covar))
               }
          else if(is.null(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- rep(ScaleF, LIV)}
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
               tuning <- sqrt(diag(Covar))}
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != length(Initial.Values))
                    tuning <- rep(ScaleF, LIV)}
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
     if(Algorithm == "Adaptive Directional Metropolis-within-Gibbs") {
          mcmc.out <- ADMG(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, LogFile)}
     else if(Algorithm == "Adaptive Griddy-Gibbs") {
          mcmc.out <- AGG(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, LogFile)}
     else if(Algorithm == "Adaptive Hamiltonian Monte Carlo") {
          mcmc.out <- AHMC(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Affine-Invariant Ensemble Sampler") {
          mcmc.out <- AIES(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Adaptive Metropolis") {
          mcmc.out <- AM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & !is.list(VarCov)) {
          mcmc.out <- AMM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & is.list(VarCov)) {
          mcmc.out <- AMM.B(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- AMWG(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, LogFile)}
     else if(Algorithm == "Componentwise Hit-And-Run Metropolis") {
          mcmc.out <- CHARM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Delayed Rejection Adaptive Metropolis") {
          mcmc.out <- DRAM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Delayed Rejection Metropolis") {
          mcmc.out <- DRM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Differential Evolution Markov Chain") {
          mcmc.out <- DEMC(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Elliptical Slice Sampler") {
          mcmc.out <- Ess(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               VarCov, LogFile)}
     else if(Algorithm == "Experimental") {
#          mcmc.out <- Experimental(Model, Data, Iterations, Status,
#               Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0,
#               ScaleF, thinned, LogFile)}
          stop("Experimental function not found.", file=LogFile,
               append=TRUE)}
     else if(Algorithm == "Gibbs Sampler") {
          mcmc.out <- Gibbs(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, LogFile)}
     else if(Algorithm == "Griddy-Gibbs") {
          mcmc.out <- GG(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Hamiltonian Monte Carlo") {
          mcmc.out <- HMC(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               VarCov, LogFile)}
     else if(Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging") {
          mcmc.out <- HMCDA(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, LogFile)}
     else if(Algorithm == "Hit-And-Run Metropolis") {
          mcmc.out <- HARM(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, LogFile)}
     else if(Algorithm == "Independence Metropolis") {
          mcmc.out <- IM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               VarCov, LogFile)}
     else if(Algorithm == "Interchain Adaptation") {
          mcmc.out <- INCA(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Metropolis-Adjusted Langevin Algorithm") {
          mcmc.out <- MALA(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               VarCov, LogFile)}
     else if(Algorithm == "Metropolis-Coupled Markov Chain Monte Carlo") {
          mcmc.out <- MCMCMC(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Multiple-Try Metropolis") {
          mcmc.out <- MTM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, tuning,
               LogFile)}
     else if(Algorithm == "Metropolis-within-Gibbs") {
          mcmc.out <- MWG(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, LogFile)}
     else if(Algorithm == "No-U-Turn Sampler") {
          mcmc.out <- NUTS(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Oblique Hyperrectangle Slice Sampler") {
          mcmc.out <- OHSS(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               VarCov, LogFile)}
     else if(Algorithm == "Random Dive Metropolis-Hastings") {
          mcmc.out <- RDMH(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, LogFile)}
     else if(Algorithm == "Random-Walk Metropolis") {
          mcmc.out <- RWM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, VarCov, LogFile)}
     else if(Algorithm == "Refractive Sampler") {
          mcmc.out <- Refractive(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned,
               LogFile)}
     else if(Algorithm == "Reflective Slice Sampler") {
          mcmc.out <- RSS(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned,
               LogFile)}
     else if(Algorithm == "Reversible-Jump") {
          mcmc.out <- RJ(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Robust Adaptive Metropolis") {
          mcmc.out <- RAM(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               VarCov, LogFile)}
     else if(Algorithm == "Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- SAMWG(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data$parm.names, LogFile)}
     else if(Algorithm == "Sequential Metropolis-within-Gibbs") {
          mcmc.out <- SMWG(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               tuning, parm.names=Data$parm.names, LogFile)}
     else if(Algorithm == "Stochastic Gradient Langevin Dynamics") {
          mcmc.out <- SGLD(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "Slice Sampler") {
          mcmc.out <- Slice(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, LogFile)}
     else if(Algorithm == "Tempered Hamiltonian Monte Carlo") {
          mcmc.out <- THMC(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               LogFile)}
     else if(Algorithm == "t-walk") {
          mcmc.out <- twalk(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, LogFile)}
     else if(Algorithm == "Univariate Eigenvector Slice Sampler") {
          mcmc.out <- UESS(Model, Data, Iterations, Status, Thinning, Specs,
               Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
               VarCov, LogFile)}
     else if(Algorithm == "Updating Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- USAMWG(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data$parm.names, LogFile)}
     else if(Algorithm == "Updating Sequential Metropolis-within-Gibbs") {
          mcmc.out <- USMWG(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data$parm.names, LogFile)}
     else stop("The algorithm is unrecognized.", file=LogFile, append=TRUE)
     #########################  MCMC is Finished  #########################
     Acceptance <- mcmc.out$Acceptance
     Dev <- mcmc.out$Dev
     DiagCovar <- mcmc.out$DiagCovar
     Mon <- mcmc.out$Mon
     thinned <- mcmc.out$thinned
     VarCov <- mcmc.out$VarCov
     remove(mcmc.out)
     rownames(DiagCovar) <- NULL
     colnames(DiagCovar) <- Data$parm.names
     thinned <- matrix(thinned[-1,], nrow(thinned)-1, ncol(thinned))
     Dev <- matrix(Dev[-1,], nrow(Dev)-1, 1)
     Mon <- matrix(Mon[-1,], nrow(Mon)-1, ncol(Mon))
     if(is.matrix(VarCov) & !is.list(VarCov)) {
          colnames(VarCov) <- rownames(VarCov) <- Data$parm.names}
     else if(is.vector(VarCov) & !is.list(VarCov)) {
          names(VarCov) <- Data$parm.names}
     thinned.rows <- nrow(thinned)
     ### Warnings (After Updating)
     if(any(Acceptance == 0))
          cat("\nWARNING: All proposals were rejected.\n", file=LogFile,
               append=TRUE)
     ### Real Values
     thinned <- ifelse(!is.finite(thinned), 0, thinned)
     Dev <- ifelse(!is.finite(Dev), 0, Dev)
     Mon <- ifelse(!is.finite(Mon), 0, Mon)
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
     acf.temp <- matrix(1, trunc(10*log10(thinned.rows)), LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV) {
          temp0 <- acf(thinned[,j], lag.max=nrow(acf.temp), plot=FALSE)
          acf.temp[,j] <- abs(temp0$acf[2:{nrow(acf.temp)+1},,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning}
     Rec.Thin <- ifelse(is.na(Rec.Thin), nrow(acf.temp), Rec.Thin)
     ### Assess ESS for all deviance and monitor samples
     ESS2 <- ESS(Dev)
     ESS3 <- ESS(Mon)
     ### Assess ESS for stationary samples
     if(Stat.at < thinned.rows) {
          ESS4 <- ESS(thinned[Stat.at:thinned.rows,])
          ESS5 <- ESS(Dev[Stat.at:thinned.rows,])
          ESS6 <- ESS(Mon[Stat.at:thinned.rows,])}
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n", file=LogFile, append=TRUE)
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- apply(thinned, 2, sd)
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025), na.rm=TRUE)
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500), na.rm=TRUE)
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)), silent=TRUE)
     if(inherits(temp, "try-error"))
          temp <- MCSE(as.vector(Dev), method="sample.variance")
     Deviance[3] <- temp
     Deviance[4] <- ESS2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
          if(inherits(temp, "try-error")) 
               temp <- MCSE(Mon[,j], method="sample.variance")
          Monitor[3] <- temp
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data$mon.names[j]}
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(Stat.at < thinned.rows) {
          thinned2 <- matrix(thinned[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(thinned))
          Dev2 <- matrix(Dev[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Dev))
          Mon2 <- matrix(Mon[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Mon))
          Summ2[,1] <- colMeans(thinned2)
          Summ2[,2] <- apply(thinned2, 2, sd)
          Summ2[,3] <- 0
          Summ2[,4] <- ESS4
          Summ2[,5] <- apply(thinned2, 2, quantile, c(0.025), na.rm=TRUE)
          Summ2[,6] <- apply(thinned2, 2, quantile, c(0.500), na.rm=TRUE)
          Summ2[,7] <- apply(thinned2, 2, quantile, c(0.975), na.rm=TRUE)
          for (i in 1:ncol(thinned2)) {
               temp <- try(MCSE(thinned2[,i]), silent=TRUE)
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else Summ2[i,3] <- MCSE(thinned2[,i],
                    method="sample.variance")}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- try(MCSE(as.vector(Dev2)), silent=TRUE)
          if(inherits(temp, "try-error"))
               temp <- MCSE(as.vector(Dev2), method="sample.variance")
          Deviance[3] <- temp
          Deviance[4] <- ESS5
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
               temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
               if(inherits(temp, "try-error"))
                    temp <- MCSE(as.vector(Mon[,j]),
                    method="sample.variance")
               Monitor[3] <- temp
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon2[,j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon2[,j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon2[,j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data$mon.names[j]}
          }
     ### Column names to samples
     if(identical(ncol(Mon), length(Data$mon.names)))
          colnames(Mon) <- Data$mon.names
     if(identical(ncol(thinned), length(Data$parm.names))) {
          colnames(thinned) <- Data$parm.names}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     if(Algorithm %in% c("Adaptive Griddy-Gibbs",
          "Affine-Invariant Ensemble Sampler",
          "Componentwise Hit-And-Run Metropolis",
          "Delayed Rejection Metropolis",
          "Elliptical Slice Sampler",
          "Gibbs Sampler",
          "Griddy-Gibbs",
          "Hamiltonian Monte Carlo",
          "Hit-And-Run Metropolis",
          "Independence Metropolis",
          "Metropolis-Coupled Markov Chain Monte Carlo",
          "Metropolis-within-Gibbs",
          "Multiple-Try Metropolis",
          "No-U-Turn Sampler",
          "Random Dive Metropolis-Hastings",
          "Random-Walk Metropolis",
          "Reflective Slice Sampler",
          "Refractive Sampler",
          "Reversible-Jump",
          "Sequential Metropolis-within-Gibbs",
          "Slice Sampler",
          "Stochastic Gradient Langevin Dynamics",
          "Tempered Hamiltonian Monte Carlo",
          "t-walk") & {Stat.at < thinned.rows}) {
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
ADMG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     LogFile)
     {
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     obs.sum <- matrix(0, LIV, 1)
     obs.scatter <- matrix(0, LIV, LIV)
     s <- svd(VarCov)
     U <- colSums(s$u)
     tol <- LIV*max(s$d)*.Machine$double.eps
     problem <- any(s$d <= tol)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          AccRate <- as.vector(Acceptance) / iter
          if(problem == FALSE) lambda <- rnorm(LIV, U,
               0.01 + s$d * exp(2*s$d*(AccRate - 0.3)))
          else lambda <- rnorm(LIV, 0, ScaleF)
          lambda[which(AccRate < 0.05)] <- rnorm(length(which(AccRate < 0.05)),
               0, sqrt(0.0001 * ScaleF))
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- prop[j] + lambda[j]
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + Mo0[["parm"]]
          obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
          ### Adaptation
          if(iter %% Periodicity == 0) {
               VarCov <- obs.scatter/iter - tcrossprod(obs.sum/iter)
               diag(VarCov) <- diag(VarCov) + 1e-05
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               s <- svd(VarCov)
               U <- colSums(s$u)
               tol <- LIV*max(s$d)*.Machine$double.eps
               problem <- any(s$d <= tol)}
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
AGG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     LogFile)
     {
     Grid <- Specs[["Grid"]]
     dparm <- Specs[["dparm"]]
     smax <- Specs[["smax"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     AGGCP <- function(Model, Data, j, Mo0, Grid, tuning, smax)
          {
          G <- length(Grid[[j]])
          x <- Grid[[j]] * sqrt(2) * tuning[j]
          LP.grid <- rep(0, G)
          prop <- Mo0[["parm"]]
          theta <- prop[j] + x
          for (g in 1:G) {
               prop[j] <- theta[g]
               Mo1 <- Model(prop, Data)
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
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) Mo1 <- Mo0
          else tuning[j] <- min(max(sqrt(sum(LP.grid * x^2)),
                    1e-10), smax)
          Mo0 <- Mo1
          return(list(Mo0=Mo0, tuning=tuning))
          }
     AGGCPP <- function(Model, Data, j, Mo0, Grid, tuning, smax, cl)
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
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]])))) Mo1 <- Mo0
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
                         ",   Proposal: Componentwise,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample(LIV)) {
                    if(j %in% dparm) Mo0 <- GGDP(Model, Data, j, Mo0, Grid)
                    else {
                         agg <- AGGCP(Model, Data, j, Mo0, Grid, tuning,
                              smax)
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
                         ",   Proposal: Componentwise,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- GGDPP(Model, Data, j, Mo0, Grid, cl)
                    else {
                         agg <- AGGCPP(Model, Data, j, Mo0, Grid,
                              tuning, smax, cl)
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
          Mon=Mon, thinned=thinned, VarCov=apply(thinned, 2, var))
     return(out)
     }
AHMC <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     DiagCovar[1,] <- epsilon
     gr0 <- partial(Model, post[1,], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- post[iter,]
          momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
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
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Adaptation
          if({iter > 10} & {iter %% Periodicity == 0}) {
               acceptances <- apply(post[(iter-9):iter,], 2, function(x)
                    {length(unique(x))})
               eps.num <- which(acceptances <= 1)
               epsilon[eps.num] <- epsilon[eps.num] * 0.8
               eps.num <- which(acceptances > 7)
               epsilon[eps.num] <- epsilon[eps.num] * 1.2
               DiagCovar <- rbind(DiagCovar, epsilon)}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
AIES <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
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
               if(!is.null(Data$PGF)) {
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
                         cat(",   Proposal: Multivariate,   LP:",
                              round(Mo0[[1]][["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0[[i]]
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
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
                         Mo1[[i]][["Dev"]], Mo1[[i]][["Monitor"]]))))
                         Mo1[[i]] <- Mo0[[i]]
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[[i]][["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
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
          VarCov=apply(thinned, 2, var))
     return(out)
     }
AM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
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
               silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
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
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)}
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
AMM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Block <- Specs[["B"]]
     Periodicity <- Specs[["Periodicity"]]
     w <- Specs[["w"]]
     obs.sum <- matrix(0, LIV, 1)
     obs.scatter <- matrix(0, LIV, LIV)
     if(all(upper.triangle(VarCov) == 0)) prop.R <- NULL
     else prop.R <- ScaleF * chol(VarCov)
     tuning <- sqrt(0.0001 * ScaleF)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values from a mixture
          if(is.null(prop.R) || runif(1) < w) {
               prop <- rnorm(LIV, Mo0[["parm"]], tuning)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Non-Adaptive Component,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)}
          else {
               prop <- Mo0[["parm"]] +
                    as.vector(rbind(rnorm(LIV)) %*% prop.R)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Adaptive Component,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + Mo0[["parm"]]
          obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               VarCov <- obs.scatter/iter - tcrossprod(obs.sum/iter)
               diag(VarCov) <- diag(VarCov) + 1e-05
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               prop.R <- try(ScaleF * chol(VarCov), silent=TRUE)
               if(!is.matrix(prop.R)) prop.R <- NULL}
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
AMM.B <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Block <- Specs[["B"]]
     Periodicity <- Specs[["Periodicity"]]
     w <- Specs[["w"]]
     B <- length(Block)
     if(!identical(length(VarCov), B))
          stop("Number of components in Covar differs from number of blocks.")
     obs.scatter <- obs.sum <- prop.R <- list()
     for (b in 1:B) {
          if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
               stop("Diagonal of Covar[[",b,"]] differs from block length.")
          obs.sum[[b]] <- matrix(0, length(Block[[b]]), 1)
          obs.scatter[[b]] <- matrix(0, length(Block[[b]]), length(Block[[b]]))
          if(all(upper.triangle(VarCov[[b]]) == 0)) prop.R[[b]] <- NA
          else prop.R[[b]] <- ScaleF * chol(VarCov[[b]])}
     tuning <- sqrt(0.0001 * ScaleF)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Proceed by Block
          for (b in 1:B) {
               ### Propose new parameter values from a mixture
               prop <- Mo0[["parm"]]
               if(any(is.na(prop.R[[b]])) || runif(1) < w) {
                    prop[Block[[b]]] <- rnorm(length(Block[[b]]),
                         Mo0[["parm"]][Block[[b]]], tuning)
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Non-Adaptive Component,   LP:",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)}
               else {
                    prop[Block[[b]]] <- Mo0[["parm"]][[Block[[b]]]] +
                         as.vector(rbind(rnorm(length(Block[[b]]))) %*%
                              prop.R[[b]])
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Adaptive Component,   LP:",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)}
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + length(Block[[b]]) / LIV
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}}
               ### Update Sample and Scatter Sum
               obs.sum[[b]] <- obs.sum[[b]] + Mo0[["parm"]][Block[[b]]]
               obs.scatter[[b]] <- obs.scatter[[b]] +
                    tcrossprod(Mo0[["parm"]][Block[[b]]])
               ### Adapt the Proposal Variance
               if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
                    VarCov[[b]] <- obs.scatter[[b]]/iter -
                         tcrossprod(obs.sum[[b]]/iter)
                    diag(VarCov[[b]]) <- diag(VarCov[[b]]) + 1e-05
                    if(b == 1) DiagCovar <- rbind(DiagCovar, rep(0,LIV))
                    DiagCovar[nrow(DiagCovar),Block[[b]]] <- diag(VarCov[[b]])
                    prop.R[[b]] <- try(ScaleF * chol(VarCov[[b]]), silent=TRUE)
                    if(!is.matrix(prop.R[[b]])) prop.R[[b]] <- NA}
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
AMWG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     LogFile)
     {
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
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
CHARM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     if(is.na(alpha.star)) {
          Acceptance <- matrix(0, 1, LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + u}
               if(iter %% Thinning == 0) {
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
               VarCov=apply(thinned, 2, var))
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
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + tau[j]*lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         tau[j] <- tau[j] + (tau[j] / (alpha.star *
                         (1 - alpha.star))) * (1 - alpha.star) / iter
                         }
                    else {
                         tau[j] <- abs(tau[j] - (tau[j] / (alpha.star *
                         (1 - alpha.star))) * alpha.star / iter)
                         }
                    Acceptance[j] <- Acceptance[j] + u}
               if(iter %% Thinning == 0) {
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
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     }
DEMC <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
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
                         if(!is.null(Data$PGF)) {
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
                    cat(",   Proposal: Multivariate,   LP:",
                         round(Mo0[[1]][["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0[[i]]
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
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=thinned,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
DRAM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     DR <- 1
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
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
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov), silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
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
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- Model(prop, Data)
               if(any(!is.finite(c(Mo2[["LP"]], Mo2[["Dev"]],
                    Mo2[["Monitor"]]))))
                    Mo2 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[["LP"]] - Mo2[["LP"]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                    {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] - Mo0[["LP"]]))}
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
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)}
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
DRM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
     {
     DR <- 1
     U <- chol(VarCov)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% U,
               silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               prop <- Mo0[["parm"]]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVNz <- try(rbind(rnorm(LIV)) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
               else {
                    prop <- Mo0[["parm"]]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- Model(prop, Data)
               if(any(!is.finite(c(Mo2[["LP"]], Mo2[["Dev"]],
                    Mo2[["Monitor"]]))))
                    Mo2 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[["LP"]] - Mo2[["LP"]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                    {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] - Mo0[["LP"]]))}
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo2
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
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
Ess <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     LogFile)
     {
     Block <- Specs[["B"]]
     if(length(Block) == 0) {
          nu <- rnorm(LIV, 0, diag(VarCov))
          U <- chol(VarCov)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Propose new parameter values
               nu <- as.vector(rbind(rnorm(LIV)) %*% U)
               theta <- theta.max <- runif(1, 0, 2*pi)
               theta.min <- theta - 2*pi
               shrink <- TRUE
               log.u <- log(runif(1))
               ### Rejection Sampling
               while (shrink == TRUE) {
                    prop <- Mo0[["parm"]]*cos(theta) + nu*sin(theta)
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo0[["parm"]]
                              Dev[t.iter] <- Mo0[["Dev"]]
                              Mon[t.iter,] <- Mo0[["Monitor"]]}
                         shrink <- FALSE
                         }
                    else {
                         if(theta < 0) theta.min <- theta
                         else theta.max <- theta
                         theta <- runif(1, theta.min, theta.max)}}
               }
          }
     else {
          B <- length(Block)
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from number of blocks.")
          nu <- rep(NA, LIV)
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from block length.")
               nu[Block[[b]]] <- rnorm(length(Block[[b]]), 0,
                    diag(VarCov[[b]]))}
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
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
                    while (shrink == TRUE) {
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]]*cos(theta) +
                              nu[Block[[b]]]*sin(theta)
                         ### Log-Posterior of the proposed state
                         Mo1 <- Model(prop, Data)
                         if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                              Mo1[["Monitor"]]))))
                              Mo1 <- Mo0
                         ### Accept/Reject
                         log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                         if(!is.finite(log.alpha)) log.alpha <- 0
                         if(log.u < log.alpha) {
                              Mo0 <- Mo1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo0[["parm"]]
                                   Dev[t.iter] <- Mo0[["Dev"]]
                                   Mon[t.iter,] <- Mo0[["Monitor"]]}
                              shrink <- FALSE
                              }
                         else {
                              if(theta < 0) theta.min <- theta
                              else theta.max <- theta
                              theta <- runif(1, theta.min, theta.max)}}
                    }
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
Gibbs <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     LogFile)
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
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Gibbs Sampling of Full Conditionals
          prop <- FC(Mo0[["parm"]], Data)
          Mo0 <- Model(prop, Data)
          ### Metropolis-within-Gibbs
          if(MWGlen > 0) {
               ### Random-Scan Componentwise Estimation
               for (j in sample(MWG)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- rnorm(1, prop[j], tuning[j])
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + u}
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
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
GG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
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
                         ",   Proposal: Componentwise,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample(LIV)) {
                    if(j %in% dparm) Mo0 <- GGDP(Model, Data, j, Mo0, Grid)
                    else Mo0 <- GGCP(Model, Data, j, Mo0, Grid)
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
                         ",   Proposal: Componentwise,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               for (j in sample(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- GGDPP(Model, Data, j, Mo0, Grid, cl)
                    else Mo0 <- GGCPP(Model, Data, j, Mo0, Grid, cl)
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]}}}
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar,
          Mon=Mon, thinned=thinned, VarCov=apply(thinned, 2, var))
     return(out)
     }
### Griddy-Gibbs Continuous Parameter (Non-Parallelized)
GGCP <- function(Model, Data, j, Mo0, Grid)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     prop <- Mo0[["parm"]]
     theta <- prop[j] + Grid[[j]]
     for (g in 1:G) {
          prop[j] <- theta[g]
          Mo1 <- Model(prop, Data)
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
     Mo1 <- Model(prop, Data)
     if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Continuous Parameter (Parallelized)
GGCPP <- function(Model, Data, j, Mo0, Grid, cl)
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
     Mo1 <- Model(prop, Data)
     if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Discrete Parameter (Non-Parallelized)
#where j is which parameter, and Grid are discrete values
GGDP <- function(Model, Data, j, Mo0, Grid)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     prop <- Mo0[["parm"]]
     theta <- Grid[[j]]
     for (g in 1:G) {
          prop[j] <- theta[g]
          Mo1 <- Model(prop, Data)
          LP.grid[g] <- Mo1[["LP"]]
          theta[g] <- Mo1[["parm"]][j]}
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     prop[j] <- sample(theta, 1, prob=LP.grid)
     Mo1 <- Model(prop, Data)
     if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Discrete Parameter (Parallelized)
GGDPP <- function(Model, Data, j, Mo0, Grid, cl)
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
     Mo1 <- Model(prop, Data)
     if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
          Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
HARM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     Block <- Specs[["B"]]
     if(is.na(alpha.star) & {length(Block) == 0}) {
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
               ### Propose new parameter values
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1) * d
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                         ",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     else if(is.na(alpha.star) & {length(Block) > 0}) {
          B <- length(Block)
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
                    prop[Block[[b]]] <- prop[Block[[b]]] + runif(1) * d
                    if({b == 1} & {iter %% Status == 0}) 
                         cat(",   Proposal: Multivariate,   LP:",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + length(Block[[b]]) / LIV
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}
                         }
                    }
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
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
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Propose new parameter values
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1,0,tau) * d
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1
                    tau <- tau + (tau / (alpha.star *
                         (1 - alpha.star))) * (1 - alpha.star) / iter
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]
                         DiagCovar[t.iter,] <- tau}
                    }
               else {
                    tau <- abs(tau - (tau / (alpha.star *
                         (1 - alpha.star))) * alpha.star / iter)
                    if(iter %% Thinning == 0) DiagCovar[t.iter,] <- tau}
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
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
                         cat(",   Proposal: Multivariate,   LP:",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + length(Block[[b]]) / LIV
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
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     }
HMC <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     gr0 <- partial(Model, Mo0[["parm"]], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
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
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=matrix(epsilon, 1, LIV),
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
HMCDA <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     A <- Specs[["A"]]
     delta <- Specs[["delta"]]
     epsilon <- Specs[["epsilon"]]
     Lmax <- Specs[["Lmax"]]
     lambda <- Specs[["lambda"]]
     leapfrog <- function(theta, r, grad, epsilon, Model, Data)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- Model(thetaprime, Data)
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
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
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
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
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
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum1 <- momentum0 <- runif(LIV)
          joint <- Mo0[["LP"]] - 0.5 * as.vector(momentum0 %*% momentum0)
          L <- max(1, round(lambda / epsilon))
          L <- min(L, Lmax)
          gr1 <- gr0
          Mo0.1 <- Mo0
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Leapfrog Function
          for (l in 1:L) {
               momentum1 <- momentum1 + 0.5 * epsilon * gr1
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
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
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
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
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
IM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     LogFile)
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
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% U, silent=TRUE)
          if(!inherits(MVNz, "try-error")) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate,   LP:",
                        round(Mo0[["LP"]],1), "\n", sep="",
                        file=LogFile, append=TRUE)
               prop <- as.vector(mu) + as.vector(MVNz)}
          else {prop <- as.vector(Mo0[["parm"]])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
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
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
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
INCA <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
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
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov), silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               prop <- as.vector(post[iter,]) + as.vector(MVNz)}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:",
                         round(Mo0[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
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
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Save log.alpha
          if({iter %% Periodicity} == 0)
               tmpAlpha[Periodicity] <- min(1, exp(log.alpha))
          else tmpAlpha[(iter %% Periodicity)] <- min(1, exp(log.alpha))
          ### Shrinkage of Adaptive Proposal Variance
          if({iter < Adaptive} & {Acceptance > 5} & {Acceptance / iter < 0.05}) {
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
               DiagCovar[iter/Periodicity,] <- diag(VarCov)
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
MALA <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     VarCov, LogFile)
     {
     A <- Specs[["A"]]
     alpha.star <- Specs[["alpha.star"]]
     delta <- Specs[["delta"]]
     gamma.const <- Specs[["gamma"]]
     epsilon <- Specs[["epsilon"]]
     Gamm <- VarCov
     mu <- Mo0[["parm"]]
     sigma2 <- 1 / (LIV*LIV)
     DiagCovar <- matrix(diag(Gamm), nrow(thinned), LIV)
     Iden <- diag(LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          gr <- partial(Model, Mo0[["parm"]], Data)
          Dx <- {delta/max(delta, abs(gr))}*gr
          gamm <- min(gamma.const/iter, 1)
          Lambda <- Gamm + epsilon[2]*Iden
          prop <- as.vector(rmvn(1, Mo0[["parm"]] + {sigma2/2}*
               as.vector(tcrossprod(Lambda, t(Dx)))*Dx,
               sigma2*Lambda))
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]
                    DiagCovar[t.iter,] <- diag(Lambda)}}
          ### Adapt Gamma (first, since it uses mu[t] not [t+1])
          xmu <- Mo0[["parm"]] - mu
          Gamm.prop <- Gamm + gamm*{xmu %*% t(xmu) - Gamm}
          norm.Gamm <- norm(Gamm.prop, type="F")
          if(norm.Gamm <= A) Gamm <- Gamm.prop               
          else Gamm <- {A/Gamm.prop}*Gamm.prop
          ### Adapt mu
          mu.prop <- mu + gamm*(Mo0[["parm"]] - mu)
          norm.mu <- sqrt(sum(mu.prop*mu.prop))
          if(norm.mu <= A) mu <- mu.prop
          else mu <- {A/norm.mu}*mu.prop
          ### Adapt sigma
          sigma2 <- interval(sqrt(sigma2) +
               gamm*(min(exp(log.alpha),1) - alpha.star),
               epsilon[1], A, reflect=FALSE)^2
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
MCMCMC <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
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
                    ",   Proposal: Multivariate,   LP:",
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
          Mo1 <- parLapply(cl, 1:CPUs, function(x) Model(prop[x,], Data))
          for (i in 1:CPUs) {
               if(any(!is.finite(c(Mo1[[i]][["LP"]], Mo1[[i]][["Dev"]],
                    Mo1[[i]][["Monitor"]]))))
                    Mo1[[i]] <- Mo0[[i]]}
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
     cat("\nSwap Acceptance Rate:", round(Acceptance.swap / Iterations, 5), "\n")
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
MTM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, tuning, LogFile)
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
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop1 <- matrix(Mo0[["parm"]], K, LIV, byrow=TRUE)
               prop1[,j] <- rnorm(K, prop1[,j], tuning[j])
               ### Log-Posterior of the proposed states
               if(CPUs == 1) {
                    ### Non-parallel
                    for (k in 1:K) {
                         Mo1[[k]] <- Model(prop1[k,], Data)
                         if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
                         LP[k] <- LW[k] <- Mo1[[k]][["LP"]]
                         prop1[k,] <- Mo1[[k]][["parm"]]}
                    }
               else {
                    ### Parallel
                    Mo1 <- parLapply(cl, 1:K, function(x)
                         Model(prop1[x,], Data))
                    for (k in 1:K) {
                         if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
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
               Mo2 <- Model(prop5, Data)
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
                         Mo1[[k]] <- Model(prop4[k,], Data)
                         if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
                         denom[k] <- Mo1[[k]][["LP"]]}
                    }
               else {
                    ### Parallel
                    Mo1 <- parLapply(cl, 1:K, function(x)
                         Model(prop4[x,], Data))
                    for (k in 1:K) {
                         if(any(!is.finite(c(Mo1[[k]][["LP"]],
                              Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
                         denom[k] <- Mo1[[k]][["LP"]]}}
               denom <- logadd(denom)
               ### Accept/Reject
               u <- log(runif(1)) < (numerator - denom)
               if(u == TRUE) Mo0 <- Mo2
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
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
          VarCov=apply(thinned, 2, var))
     return(out)
     }
MWG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     LogFile)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
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
NUTS <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     A <- Specs[["A"]]
     delta <- Specs[["delta"]]
     epsilon <- Specs[["epsilon"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     leapfrog <- function(theta, r, grad, epsilon, Model, Data)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- Model(thetaprime, Data)
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
     build.tree <- function(theta, r, grad, logu, v, j, epsilon, joint0)
          {
          if(j == 0) {
               ### Base case: Take a single leapfrog step in direction v
               leap <- leapfrog(theta, r, grad, v*epsilon, Model, Data)
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
               tree <- build.tree(theta, r, grad, logu, v, j-1, epsilon,
                    joint0)
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
                         tree <- build.tree(thetaminus, rminus, gradminus,
                              logu, v, j-1, epsilon, joint0)
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
                         tree <- build.tree(thetaplus, rplus, gradplus,
                              logu, v, j-1, epsilon, joint0)
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
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
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
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
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
          epsilon <- find.reasonable.epsilon(post[1,], grad, Mo0, Model,
               Data, LogFile)
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
                    ",   Proposal: Multivariate,   LP:",
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
                    tree <- build.tree(thetaminus, rminus, gradminus,
                         logu, v, j, epsilon, joint)
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
                    tree <- build.tree(thetaplus, rplus, gradplus, logu,
                         v, j, epsilon, joint)
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
               j <- j + 1}
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
          VarCov=apply(thinned, 2, var))
     return(out)
     }
OHSS <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     VarCov, LogFile)
     {
     w <- 0.05 # as with Roberts & Rosenthal
     decomp.freq <- max(floor(Iterations / Thinning / 100), 10)
     S.eig <-try(eigen(VarCov), silent=TRUE)
     if(inherits(S.eig, "try-error")) S.eig <- NULL
     tuning <- 1 #Tuning
     edge.scale <- 5 #Tuning
     t.iter <- 1
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Eigenvectors of the Sample Covariance Matrix
          if({iter %% decomp.freq == 0} & {iter > 1}) {
               S.eig <- try(eigen(cov(thinned[1:(t.iter-1),,drop=FALSE])),
                    silent=TRUE)
               if(inherits(S.eig, "try-error")) S.eig <- eigen(diag(LIV))}
          ### Hypercube or Eigenvector
          if(runif(1) < w || is.null(S.eig)) {
               vals <- rep(tuning, LIV)
               vecs <- diag(1, nrow=LIV)
               }
          else {
               vals <- S.eig$values
               vecs <- S.eig$vectors}
          ### Slice Interval
          y.slice <- Model(Mo0[["parm"]], Data)[["LP"]] - rexp(1)
          L <- -1 * runif(LIV)
          U <- L + 1
          ### Rejection Sampling
          repeat {
               wt <- runif(LIV, min=L, max=U)
               v <- as.numeric(vecs %*% (edge.scale * wt * vals))
               prop <- Mo0[["parm"]] + v
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               y1 <- Mo1[["LP"]]
               if(y1 >= y.slice) break
               L[wt < 0] <- wt[wt < 0]
               U[wt > 0] <- wt[wt > 0]}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo1[["parm"]]
               Dev[t.iter] <- Mo1[["Dev"]]
               Mon[t.iter,] <- Mo1[["Monitor"]]
               DiagCovar <- rbind(DiagCovar, S.eig$vectors)}
          Mo0 <- Mo1
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
RAM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
     LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     Dist <- Specs[["Dist"]]
     gamma <- Specs[["gamma"]]
     Periodicity <- Specs[["Periodicity"]]
     if(!is.symmetric.matrix(VarCov)) {
          cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
               append=TRUE)
          VarCov <- as.symmetric.matrix(VarCov)}
     if(!is.positive.definite(VarCov)) {
          cat("\nNon-Positive-Definite Covar, correcting now...\n",
               file=LogFile, append=TRUE)
          VarCov <- as.positive.definite(VarCov)}
     Iden.Mat <- diag(LIV)
     S.z <- try(t(chol(VarCov)), silent=TRUE)
     if(!inherits(S.z, "try-error")) S <- S.z
     else S <- Iden.Mat
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose New Parameter Values
          if(Dist == "t") U <- rt(LIV, df=5)
          else U <- rnorm(LIV)
          prop <- Mo0[["parm"]] + rbind(U) %*% S
          ### Log-Posterior
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Adaptation
          if(iter %% Periodicity == 0) {
               eta <- min(1, LIV*iter^(-gamma))
               VarCov.test <- S %*% {Iden.Mat +
                    eta*(min(1, exp(log.alpha)) - alpha.star) *
                    U %*% t(U) / sum(U*U)} %*% t(S)
               if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
                    !is.matrix(VarCov.test)) {VarCov.test <- VarCov}
               if(!is.symmetric.matrix(VarCov.test))
                    VarCov.test <- as.symmetric.matrix(VarCov.test)
               if(is.positive.definite(VarCov.test)) {
                    S.z <- try(t(chol(VarCov)), silent=TRUE)
                    if(!inherits(S.z, "try-error")) {
                         VarCov <- VarCov.test
                         S <- S.z}}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))}
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
RDMH <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     LogFile)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               epsilon <- runif(1,-1,1)^sample(c(-1,1),1)
               prop[j] <- prop[j]*epsilon
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               epsilon <- log(abs(Mo1[["parm"]][j] / Mo0[["parm"]][j]))
               if(!is.finite(epsilon)) epsilon <- 0
               ### Accept/Reject
               u <- log(runif(1)) < (epsilon + Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
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
          VarCov=apply(thinned,2,var))
     return(out)
     }
Refractive <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, LogFile)
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
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
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
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
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
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]
                    if(Adaptive < Iterations) DiagCovar[t.iter,] <- w}
               }
          else if(Adaptive < Iterations) {
               w <- abs(w - (w / (alpha.star * (1 - alpha.star))) *
                    alpha.star / iter)
               if(iter %% Thinning == 0) DiagCovar[t.iter,] <- w}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
RJ <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
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
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
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
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
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
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
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
          if(iter %% Thinning == 0) {
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
          VarCov=apply(thinned, 2, var))
     return(out)
     }
RSS <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, LogFile)
     {
     m <- Specs[["m"]]
     w <- Specs[["w"]]
     reflections <- 0
     Norm <- function(x) return(sqrt(sum(x*x)))
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          prop <- Mo0[["parm"]]
          g <- partial(Model, prop, Data)
          p <- rnorm(LIV)
          reflections <- 0
          ### Take m Steps
          for (i in 1:m) {
               prop <- prop + w*p
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               prop <- Mo1[["parm"]]
               ### Reflect at boundary
               if(Mo0[["LP"]] > Mo1[["LP"]]) {
                    reflections <- reflections + 1
                    p <- p - 2*g*{(t(p) %*% g) / Norm(g)^2}}}
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          prop <- Mo1[["parm"]]
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo1[["parm"]]
               Dev[t.iter] <- Mo1[["Dev"]]
               Mon[t.iter,] <- Mo1[["Monitor"]]}
          Mo0 <- Mo1
          }
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
RWM <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     VarCov, LogFile)
     {
     Block <- Specs[["B"]]
     if(length(Block) == 0) {
          U <- chol(VarCov)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Propose new parameter values
               prop <- as.vector(Mo0[["parm"]] +
                    rbind(rnorm(LIV)) %*% U)
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
               }
          }
     else {
          B <- length(Block)
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from number of blocks.")
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from block length.")}
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
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- Mo0[["parm"]][[Block[[b]]]] +
                         rbind(rnorm(length(Block[[b]]))) %*%
                         chol(VarCov[[b]])
                    if({b == 1} & {iter %% Status == 0})
                         cat(",   Proposal: Multivariate,   LP:",
                              round(Mo0[["LP"]],1), "\n", sep="",
                              file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                         Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + 1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}
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
          VarCov=VarCov)
     return(out)
     }
SAMWG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          else staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn,1,sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
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
SGLD <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
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
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Sample Data
          seek(con, 0)
          skip.rows <- sample.int(Nr - size, size=1)
          Data$X <- matrix(scan(file=con, sep=",", skip=skip.rows,
               nlines=size, quiet=TRUE), size, Nc, byrow=TRUE)
          ### Propose new parameter values
          g <- partial(Model, Mo0[["parm"]], Data)
          eta <- rnorm(LIV, 0, epsilon[iter])
          prop <- Mo0[["parm"]] + {epsilon[iter]/2}*g + eta
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
               Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          Mo0 <- Mo1
          if(iter %% Thinning == 0) {
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
          VarCov=apply(thinned, 2, var))
     return(out)
     }
Slice <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     m <- Specs[["m"]]
     w <- Specs[["w"]]
     uni.slice <- function(x0, g, j, w=1, m=Inf, lower=-Inf, upper=Inf,
          gx0=NULL, Data)
          {
          logy <- gx0[["LP"]] - rexp(1)
          ### Find the initial sampling interval
          u <- runif(1,0,w)
          L <- x0[j] - u
          R <- x0[j] + (w-u)
          ### Expand the interval until its ends are outside the slice,
          ### or until the limit on steps is reached
          intL <- intR <- x0
          intL[j] <- L
          intR[j] <- R
          ### Unlimited number of steps
          if(is.infinite(m)) { 
               repeat {
                    if(L <= lower) break
                    intL[j] <- L
                    gL <- g(intL, Data)
                    if(!is.finite(gL[["LP"]])) {
                         L <- L + w
                         break}
                    if(gL[["LP"]] <= logy) break
                    L <- L - w}
               repeat {
                    if(R >= upper) break
                    intR[j] <- R
                    gR <- g(intR, Data)
                    if(!is.finite(gR[["LP"]])) {
                         R <- R - w
                         break}
                    if(gR[["LP"]] <= logy) break
                    R <- R + w}
               }
          else if(m > 1) {
               ### Limited number of steps
               J <- floor(runif(1,0,m))
               K <- (m-1) - J
               while (J > 0) {
                    if(L <= lower) break
                    intL[j] <- L
                    gL <- g(intL, Data)
                    if(!is.finite(gL[["LP"]])) {
                         L <- L + w
                         break}
                    if(gL[["LP"]] <= logy) break
                    L <- L - w
                    J <- J - 1}
               while (K > 0) {
                    if(R >= upper) break
                    intR[j] <- R
                    gR <- g(intR, Data)
                    if(!is.finite(gR[["LP"]])) {
                         R <- R - w
                         break}
                    R <- R + w
                    K <- K - 1}
               }
          ### Shrink the interval to lower and upper bounds
          if(L < lower) L <- lower
          if(R > upper) R <- upper
          #### Sample from the interval, shrinking it on each rejection
          prop <- x0
          repeat { 
               x1 <- runif(1,L,R)
               prop[j] <- x1
               gx1 <- g(prop, Data)
               if(is.finite(gx1[["LP"]])) {
                    x1 <- gx1[["parm"]][j]
                    if(gx1[["LP"]] >= logy) break
                    if(x1 > x0[j]) R <- x1
                    else L <- x1}
               }
          return(gx1)
          }
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Univariate Slice Sampling
               Mo0 <- Mo1 <- uni.slice(x0=Mo0[["parm"]], g=Model, j=j,
                    #w=1, m=Inf, lower=-Inf, upper=Inf, gx0=Mo0, Data)
                    w=w[j], m=m[j], lower=-Inf, upper=Inf, gx0=Mo0, Data)}
          if(iter %% Thinning == 0) {
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
          VarCov=apply(thinned, 2, var))
     return(out)
     }
SMWG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, LogFile)
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
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          else staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn,1,sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
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
THMC <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     Temperature <- Specs[["Temperature"]]
     gr <- partial(Model, Mo0[["parm"]], Data)
     sqrt.Temp <- sqrt(Temperature)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Multivariate,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum1 <- momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          Mo0.1 <- Mo0
          for (l in 1:L) {
               if(2*(l-1) < L) momentum1 <- momentum1 * sqrt.Temp
               else momentum1 <- momentum1 / sqrt.Temp
               momentum1 <- momentum1 + (epsilon/2) * gr
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr <- partial(Model, prop, Data)
               momentum1 <- momentum1 + (epsilon/2) * gr
               if(2*l > L) momentum1 <- momentum1 / sqrt.Temp
               else momentum1 <- momentum1 * sqrt.Temp}
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               kinetic0 <- kinetic1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=matrix(epsilon, 1, LIV),
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
twalk <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
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
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
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
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
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
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
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
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
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
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
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
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
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
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
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
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
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
          Thinning, Acceptance, Dev, Mon, Mo0, thinned, LogFile)
          {
          x <- x0 ### Primary vector of initial values
          xp <- xp0 ### Secondary vector of initial values
          Mo0.1 <- try(Model(x, Data), silent=TRUE)
          Mo0.2 <- try(Model(xp, Data), silent=TRUE)
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
                         ",   Proposal: Multivariate Subset,   LP:",
                         round(Mo0.1[["LP"]],1), "\n", sep="",
                         file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0.1[["parm"]]
                    Dev[t.iter] <- Mo0.1[["Dev"]]
                    Mon[t.iter,] <- Mo0.1[["Monitor"]]}
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
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- move$Mo0.1[["parm"]]
                         Dev[t.iter] <- move$Mo0.1[["Dev"]]
                         Mon[t.iter,] <- move$Mo0.1[["Monitor"]]}
                    x <- move$y
                    U <- move$propU
                    xp <- move$yp
                    Up <- move$propUp
                    }
               }
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     out <- Runtwalk(Iterations=Iterations, dim=LIV, x0=Mo0[["parm"]],
          xp0=xp0, pphi=min(LIV, n1)/LIV, at=6, aw=1.5, Model=Model,
          Data=Data, Status=Status, Thinning=Thinning,
          Acceptance=Acceptance, Dev=Dev, Mon=Mon, Mo0=Mo0,
          thinned=thinned, LogFile=LogFile)
     ### Output
     return(out)
     }
UESS <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
     VarCov, LogFile)
     {
     m <- Specs[["m"]]
     Norm <- function(x) return(sqrt(sum(x^2)))
     w <- 0.05
     decomp.freq <- max(LIV * floor(Iterations / Thinning / 100), 10)
     S.eig <-try(eigen(VarCov), silent=TRUE)
     if(inherits(S.eig, "try-error")) S.eig <- NULL
     obs.sum <- matrix(0, LIV, 1)
     obs.scatter <- matrix(0, LIV, LIV)
     scatter.N <- 0
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Eigenvectors of the Sample Covariance Matrix
          if({iter %% decomp.freq == 0} & {iter > 1})
               S.eig <- eigen(obs.scatter/scatter.N -
                    tcrossprod(obs.sum/scatter.N))
          ### Non-Adaptive or Adaptive
          if(runif(1) < w || is.null(S.eig)) {
               v <- rnorm(LIV)
               v <- v / Norm(v)
               }
          else {
               which.eig <- floor(1 + LIV * runif(1))
               v <- S.eig$vectors[,which.eig] *
                    sqrt(abs(S.eig$values[which.eig]))}
          ### Slice Interval
          y.slice <- Model(Mo0[["parm"]], Data)[["LP"]] - rexp(1)
          L <- -runif(1)
          U <- L + 1
          if(m > 0) {
               L.y <- Model(Mo0[["parm"]] + v*L, Data)[["LP"]]
               U.y <- Model(Mo0[["parm"]] + v*U, Data)[["LP"]]
               step <- 0
               while({L.y > y.slice || U.y > y.slice} && step < m) {
                    step <- step + 1
                    if(runif(1) < 0.5) {
                         L <- L - 1
                         L.y <- Model(Mo0[["parm"]] + v*L, Data)[["LP"]]
                         }
                    else {
                         U <- U + 1
                         U.y <- Model(Mo0[["parm"]] + v*U, Data)[["LP"]]}}}
          ### Rejection Sampling
          repeat {
               prop.offset <- runif(1, min=L, max=U)
               prop <- Mo0[["parm"]] + prop.offset * v
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               y1 <- Mo1[["LP"]]
               prop <- Mo1[["parm"]]
               if(y1 >= y.slice) break
               if(prop.offset < 0) L <- prop.offset
               else U <- prop.offset}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo1[["parm"]]
               Dev[t.iter] <- Mo1[["Dev"]]
               Mon[t.iter,] <- Mo1[["Monitor"]]
               DiagCovar <- rbind(DiagCovar, S.eig$vectors)
               obs.sum <- obs.sum + Mo1[["parm"]]
               obs.scatter <- obs.scatter + tcrossprod(Mo1[["parm"]])
               scatter.N <- scatter.N + 1}
          Mo0 <- Mo1}
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=cov(thinned))
     return(out)
     }
USAMWG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, LogFile)
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
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter,
                    ",   Proposal: Componentwise,   LP:",
                    round(Mo0[["LP"]],1), "\n", sep="",
                    file=LogFile, append=TRUE)
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn,1,sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo0[["parm"]]}
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
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
USMWG <- function(Model, Data, Iterations, Status, Thinning, Specs,
     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
     parm.names, LogFile)
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
                    ",   Proposal: Componentwise,   LP:",
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
          else dynsample <- as.vector(apply(Dyn,1,sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                    Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo0[["parm"]]}
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
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
