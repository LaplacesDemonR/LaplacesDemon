###########################################################################
# Consort                                                                 #
#                                                                         #
# The purpose of the Consort function is to consort with Laplace's Demon  #
# regarding an object of class demonoid.                                  #
###########################################################################

Consort <- function(object=NULL)
     {
     if(is.null(object)) stop("The object argument is empty.")
     if(!identical(class(object), "demonoid"))
          stop("Consort requires an object of class demonoid.")
     oname <- deparse(substitute(object))
     dname <- as.vector(strsplit(as.character(object$Call), "=")[3][[1]])
     cat("\n#############################################################\n")
     cat("# Consort with Laplace's Demon                              #\n")
     cat("#############################################################\n")
     print.demonoid(object)
     ### Check Acceptance.Rate
     Acc.Rate.Level <- 2
     if((object$Algorithm == "Adaptive Hamiltonian Monte Carlo") |
        (object$Algorithm == "Tempered Hamiltonian Monte Carlo")) {
          L <- object$Specs[["L"]]
          if(L == 1) {
               Acc.Rate.Low <- 0.5
               Acc.Rate.High <- 0.65}
          else {
               Acc.Rate.Low <- 0.6
               Acc.Rate.High <- 0.7}}
     else if(object$Algorithm %in% c("Adaptive Griddy-Gibbs",
          "Automated Factor Slice Sampler",
          "Elliptical Slice Sampler", "Griddy-Gibbs",
          "Oblique Hyperrectangle Slice Sampler",
          "Reflective Slice Sampler", "Slice Sampler",
          "Stochastic Gradient Langevin Dynamics",
          "Univariate Eigenvector Slice Sampler")) {
          Acc.Rate.Low <- 1
          Acc.Rate.High <- 1
          }
     else if(object$Algorithm == "Gibbs Sampler") {
          if(is.null(object$Specs[["MWG"]]))
               Acc.Rate.Low <- Acc.Rate.High <- 1
          else {
               Acc.Rate.Low <- 0.15
               Acc.Rate.High <- 0.5}
          }
     else if((object$Algorithm == "Metropolis-Adjusted Langevin Algorithm")) {
          Acc.Rate.Low <- 0.40
          Acc.Rate.High <- 0.80
          }
     else if(object$Algorithm == "Hamiltonian Monte Carlo") {
          L <- object$Specs[["L"]]
          if(L == 1) {
               Acc.Rate.Low <- 0.4
               Acc.Rate.High <- 0.8}
               
          else {
               Acc.Rate.Low <- 0.6
               Acc.Rate.High <- 0.7}
          }
     else if(object$Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging") {
          A <- object$Specs[["A"]]
          delta <- object$Specs[["delta"]]
          Lmax <- object$Specs[["Lmax"]]
          lambda <- object$Specs[["lambda"]]
          Acc.Rate.Low <- max(round(delta - 0.05, 2), 0.01)
          Acc.Rate.High <- min(round(delta + 0.05, 2), 1)
          }
     else if(object$Algorithm == "No-U-Turn Sampler") {
          A <- object$Specs[["A"]]
          delta <- object$Specs[["delta"]]
          Acc.Rate.Low <- max(round(delta - 0.05, 2), 0.01)
          Acc.Rate.High <- min(round(delta + 0.05, 2), 1)
          }
     else if(object$Algorithm == "Refractive Sampler") {
          m <- object$Specs[["m"]]
          Acc.Rate.Low <- 0.6
          Acc.Rate.High <- 0.7
          }
     else if((object$Algorithm %in% c("Reversible-Jump",
          "Multiple-Try Metropolis"))) {
          Acc.Rate.Low <- 0.10
          Acc.Rate.High <- 0.90
          }
     else {
          Acc.Rate.Low <- 0.15
          Acc.Rate.High <- 0.5}
     if(any(object$Acceptance.Rate == 0)) {
          cat("\nWARNING: Acceptance Rate = 0\n\n")
          cat(oname, " <- LaplacesDemon(Model, Data=", dname,
               ", Initial.Values,\n", sep="")
          cat("     Covar=NULL, Iterations=100000, Status=1000, ",
               "Thinning=100,\n", sep="")
          cat("     Algorithm=\"AMM\", Specs=list(Adaptive=500, ",
               "B=NULL, Periodicity=100, w=0.05))\n\n", sep="")
          stop("Try the above code before consorting again.")}
     if(any(object$Acceptance.Rate < Acc.Rate.Low)) {Acc.Rate.Level <- 1}
     else if(any(object$Acceptance.Rate > Acc.Rate.High)) {Acc.Rate.Level <- 3}
     LIV <- object$Parameters
     ### Check MCSE
     MCSE.crit <- 0.0627
     MCSE.temp <- object$Summary2[1:object$Parameters,"MCSE"] /
          object$Summary2[1:LIV,"SD"]
     MCSE.temp2 <- sum(!is.finite(MCSE.temp))
     MCSE.temp[which(!is.finite(MCSE.temp))] <- 0
     MCSE.tot <- 0
     if(MCSE.temp2 < LIV)
          MCSE.tot <- sum(MCSE.temp < MCSE.crit)
     ### Check ESS
     ESS.temp <- object$Summary2[1:LIV,"ESS"]
     ESS.temp[which(!is.finite(ESS.temp))] <- 0
     ESS.min <- min(ESS.temp)
     if(all(is.finite(object$Summary2[1:LIV,"ESS"])))
          ESS.worst <- which.min(object$Summary2[1:LIV,"ESS"])
     else ESS.worst <- which.min(object$Summary1[1:LIV,"ESS"])
     ESS.crit <- 100
     ### Check Stationarity
     Stationarity <- FALSE
     if(object$Rec.BurnIn.Thinned < object$Thinned.Samples)
          Stationarity <- TRUE
     ### Check Diminishing Adaptation (If Adaptive)
     if(nrow(object$CovarDHis) > 1)
          Dim.Adapt <- sum(diff(object$CovarDHis)) <= 0
     else Dim.Adapt <- TRUE
     ### Suggested Values
     Rec.Iterations <- trunc(object$Rec.Thinning / object$Thinning *
          object$Iterations)
     if(any(object$Acceptance.Rate == 0)) {
          Rec.Adaptive <- round(LIV * 1.0E7)
          Rec.Periodicity <- round(object$Rec.Thinning * 1.0E7)}
     else if(all(object$Acceptance.Rate > 0)) {
          Rec.Adaptive <- round(LIV / mean(object$Acceptance.Rate))
          Rec.Periodicity <- round(object$Rec.Thinning /
               mean(object$Acceptance.Rate))}
     if(Rec.Adaptive > Rec.Iterations) Rec.Adaptive <- Rec.Iterations - 1
     if(Rec.Periodicity > Rec.Iterations)
          Rec.Periodicity <- Rec.Iterations - 1
     Status.temp <- trunc(Rec.Iterations / {object$Minutes *
          Rec.Iterations / object$Iterations})
     if(Status.temp < Rec.Iterations) Rec.Status <- Status.temp
     else Rec.Status <- trunc(sqrt(Rec.Iterations))
     if(Rec.Status > Rec.Iterations) Rec.Status <- Rec.Iterations
     Rec.Thinning <- object$Rec.Thinning
     ### The Demonic Suggestion of Laplace's Demon
     cat("\nDemonic Suggestion\n\n")

     cat("Due to the combination of the following conditions,\n\n")

     cat("1. ", object$Algorithm, "\n", sep="")
     if(Acc.Rate.Level == 1) {
          cat("2. The acceptance rate (", min(object$Acceptance.Rate),
               ") is below ", Acc.Rate.Low, ".\n", sep="")}
     else if(Acc.Rate.Level == 2) {
          cat("2. The acceptance rate (", mean(object$Acceptance.Rate),
               ") is within the interval [", Acc.Rate.Low, ",",
               Acc.Rate.High, "].\n", sep="")}
     else {
          cat("2. The acceptance rate (", max(object$Acceptance.Rate),
               ") is above ", Acc.Rate.High, ".\n", sep="")}
     if(MCSE.tot < LIV) {
          cat("3. At least one target MCSE is >= ",
               MCSE.crit * 100,"% of its marginal posterior\n", sep="")
          cat("   standard deviation.\n")}
     else {
          cat("3. Each target MCSE is < ", MCSE.crit * 100,
               "% of its marginal posterior\n", sep="")
          cat("   standard deviation.\n")}
     if(ESS.min < ESS.crit) {
          cat("4. At least one target distribution has an ",
               "effective sample size\n", sep="")
          cat("   (ESS) less than ", ESS.crit,
               ". The worst mixing chain is: ", 
               rownames(object$Summary1)[ESS.worst], " (ESS=",
               object$Summary1[ESS.worst,"ESS"], ").\n", sep="")}
     else {
          cat("4. Each target distribution has an effective ",
               "sample size (ESS)\n", sep="")
          cat("   of at least ", ESS.crit, ".\n", sep="")}
     if(Stationarity == FALSE) {
          cat("5. At least one target distribution is not ",
               "stationary.\n\n", sep="")}
     else if(Stationarity == TRUE & object$Rec.BurnIn.Thinned > 0) {
          cat("5. Each target distribution became stationary by\n")
          cat("   ", object$Rec.BurnIn.Thinned + 1,
               " iterations.\n\n", sep="")}
     else {
          cat("5. Each target distribution became stationary by\n")
          cat("   ", object$Rec.BurnIn.Thinned + 1,
               " iteration.\n\n", sep="")}

     ### Determine if Laplace's Demon is appeased...
     Appeased <- FALSE
     if(!is.null(object$Specs[["Adaptive"]]))
          Adaptive <- object$Specs[["Adaptive"]]
     else Adaptive <- object$Iterations + 1
     if({Adaptive > object$Iterations} &
          {Acc.Rate.Level == 2} & {MCSE.tot >= LIV} &
          {ESS.min >= ESS.crit} & {Stationarity == TRUE}) {
          Appeased <- TRUE
          cat("Laplace's Demon has been appeased, and suggests\n")
          cat("the marginal posterior samples should be plotted\n")
          cat("and subjected to any other MCMC diagnostic deemed\n")
          cat("fit before using these samples for inference.\n\n")
          }
     else {
          if(object$Algorithm %in% c("Adaptive Directional Metropolis-within-Gibbs",
               "Adaptive Griddy-Gibbs",
               "Adaptive Metropolis-within-Gibbs",
               "Componentwise Hit-And-Run Metropolis",
               "Gibbs Sampler",
               "Griddy-Gibbs",
               "Metropolis-within-Gibbs",
               "Multiple-Try Metropolis",
               "Random Dive Metropolis-Hastings",
               "Reversible-Jump",
               "Sequential Adaptive Metropolis-within-Gibbs",
               "Sequential Metropolis-within-Gibbs",
               "Slice Sampler",
               "Updating Sequential Adaptive Metropolis-within-Gibbs",
               "Updating Sequential Metropolis-within-Gibbs")) {
               options(warn=-1)
               postcor <- try(cor(object$Posterior1), silent=TRUE)
               options(warn=0)
               if(!inherits(postcor, "try-error")) {
                    postcor <- try(quantile(abs(postcor)), silent=TRUE)
                    if(!inherits(postcor, "try-error")) {
                         cat("Quantiles of Absolute Posterior1 Correlation:\n")
                         print(postcor)
                         if(postcor["75%"] >= 0.5)
                              cat("\nPossibly excessive posterior correlation for a componentwise algorithm.")
                         cat("\n\n")}}}

          if(Dim.Adapt == FALSE) {
               cat("WARNING: Diminishing adaptation did not occur.\n")
               if(!object$Algorithm %in% c("Automated Factor Slice Sampler",
                    "Interchain Adaptation",
                    "Metropolis-Adjusted Langevin Algorithm",
                    "No-U-Turn Sampler", "Refractive Sampler",
                    "Sequential Adaptive Metropolis-within-Gibbs",
                    "Univariate Eigenvector Slice Sampler",
                    "Updating Sequential Adaptive Metropolis-within-Gibbs"))
                    cat("         A new algorithm will be suggested.\n\n")}
          
          cat("Laplace's Demon has not been appeased, and suggests\n")
          cat("copy/pasting the following R code into the R",
               " console,\n", sep="")
          cat("and running it.\n\n")

          if(object$Algorithm != "Interchain Adaptation")
               cat("Initial.Values <- as.initial.values(", oname, ")\n", sep="")
          if(object$Algorithm %in% c("Adaptive Metropolis-within-Gibbs",
               "Metropolis-within-Gibbs"))
               Time <- object$Iterations / object$Minutes
          else Time <- object$Iterations / object$Minutes / LIV
          if(Time >= 100) Fast <- TRUE
          else Fast <- FALSE

          if({Acc.Rate.Level == 2} & {MCSE.tot >= LIV} &
          {ESS.min >= ESS.crit} & {Stationarity == TRUE}) Ready <- TRUE
          else Ready <- FALSE

          Alg <- switch(object$Algorithm,
               "Adaptive Directional Metropolis-within-Gibbs"="ADMG",
               "Adaptive Griddy-Gibbs"="AGG",
               "Adaptive Hamiltonian Monte Carlo"="AHMC",
               "Adaptive Metropolis"="AM",
               "Adaptive Metropolis-within-Gibbs"="AMWG",
               "Adaptive-Mixture Metropolis"="AMM",
               "Affine-Invariant Ensemble Sampler"="AIES",
               "Automated Factor Slice Sampler"="AFSS",
               "Componentwise Hit-And-Run Metropolis"="CHARM",
               "Delayed Rejection Adaptive Metropolis"="DRAM",
               "Delayed Rejection Metropolis"="DRM",
               "Differential Evolution Markov Chain"="DEMC",
               "Elliptical Slice Sampler"="ESS",
               "Experimental"="Exper",
               "Gibbs Sampler"="Gibbs",
               "Griddy-Gibbs"="GG",
               "Hamiltonian Monte Carlo"="HMC",
               "Hamiltonian Monte Carlo with Dual-Averaging"="HMCDA",
               "Hit-And-Run Metropolis"="HARM",
               "Independence Metropolis"="IM",
               "Interchain Adaptation"="INCA",
               "Metropolis-Adjusted Langevin Algorithm"="MALA",
               "Metropolis-Coupled Markov Chain Monte Carlo"="MCMCMC",
               "Metropolis-within-Gibbs"="MWG",
               "Multiple-Try Metropolis"="MTM",
               "No-U-Turn Sampler"="NUTS",
               "Oblique Hyperrectangle Slice Sampler"="OHSS",
               "Preconditioned Crank-Nicolson"="pCN",
               "Random Dive Metropolis-Hastings"="RDMH",
               "Random-Walk Metropolis"="RWM",
               "Reflective Slice Sampler"="RSS",
               "Refractive Sampler"="Refractive",
               "Reversible-Jump"="RJ",
               "Robust Adaptive Metropolis"="RAM",
               "Sequential Adaptive Metropolis-within-Gibbs"="SAMWG",
               "Sequential Metropolis-within-Gibbs"="SMWG",
               "Slice Sampler"="Slice",
               "Stochastic Gradient Langevin Dynamics"="SGLD",
               "Tempered Hamiltonian Monte Carlo"="THMC",
               "t-walk"="t-walk",
               "Univariate Eigenvector Slice Sampler"="UESS",
               "Updating Sequential Adaptive Metropolis-within-Gibbs"="USAMWG",
               "Updating Sequential Metropolis-within-Gibbs"="USMWG")

          Componentwise <- 0
          if(Alg %in% c("ADMG","AFSS","AGG","AMWG","CHARM","GG","Gibbs",
               "MWG","RJ","SAMWG","SMWG","Slice","USAMWG","USMWG"))
               Componentwise <- 1
          if({(Alg == "ADMG") & !Dim.Adapt} |
               {(Alg == "ADMG") & !Ready}) {
               ### ADMG
               n <- object$Specs[["n"]] + object$Iterations
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"ADMG\", ",
                    "Specs=list(n=", n, ", Periodicity=", Rec.Periodicity,
                    "))\n\n", sep="")
               }
          else if({Alg == "AFSS"} |
               {(Alg == "AMWG") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "AMWG") & Dim.Adapt & !Fast & !Ready} |
               {(Alg == "AMWG") & !Dim.Adapt & Fast & !Ready} | 
               {(Alg == "AMWG") & !Dim.Adapt & !Fast & !Ready}) {
               ### AFSS
               if(Ready == TRUE) A <- 0
               else if(Alg == "AFSS") A <- object$Specs[["A"]]
               else A <- Inf
               block <- "NULL"
               if(Alg == "AFSS") m <- paste(oname, "$Specs$m", sep="")
               else m <- 100
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- "Block"
               if(Alg == "AFSS")
                    n <- object$Specs[["n"]] + object$Iterations
               else n <- 0
               if(Alg == "AFSS")
                   w <- paste(oname, "$CovarDHis[nrow(", oname,
                        "$CovarDHis),]", sep="")
               else w <- 1
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AFSS\", ",
                    "Specs=list(A=", A, ", B=", block, ", m=",
                    m, ",\n", sep="")
               cat("     n=", n, ", w=", w, "))\n\n", sep="")
               }
          else if(Alg == "AGG") {
               ### AGG
               Grid <- "GaussHermiteQuadRule(3)$nodes"
               CPUs <- substr(object$Call, 1, nchar(object$Call))[10]
               CPUs <- strsplit(CPUs, " ")
               pos <- grep("CPUs", CPUs[[1]])
               CPUs <- as.numeric(strsplit(CPUs[[1]][pos+2], ","))
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AGG\", Specs=", oname,
                    "$Specs))\n\n", sep="")
               }
          else if({(Alg == "AHMC") & Dim.Adapt & !Ready} |
               {(Alg == "HMC") & Dim.Adapt & !Ready}) {
               ### AHMC
               if(L > 1) {
                    L <- round(L*(Rec.Iterations/object$Iterations))
                    Rec.Iterations <- object$Iterations
                    Rec.Status <- object$Status
                    Rec.Thinning <- object$Thinning}
               m <- "NULL"
               if(!is.null(object$Specs[["m"]]) &
                    !identical(object$Specs[["m"]],list()))
                    m <- paste(oname, "$Specs$m", sep="")
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AHMC\",\n", sep="")
               cat("     Specs=list(epsilon=", oname, "$CovarDHis[nrow(",
                    oname, "$CovarDHis),],\n", sep="")
               cat("     L=", object$Specs[["L"]], ", m=", m,
                    ", Periodicity=", Rec.Periodicity, "))\n\n", sep="")
               }
          else if({Alg == "AIES"}) {
               Nc <- 2*LIV
               if(Acc.Rate.Level == 1)
                    Nc <- 2 + 20*(0.15 - object$Acceptance.Rate)
               if(is.null(object$Specs[["Z"]])) Z <- "NULL"
               else Z <- "object$Specs[[\"Z\"]]"
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AIES\", Specs=list(Nc=", Nc,
                    ", Z=", Z, ",\n", sep="")
               cat("     beta=", object$Specs[["beta"]],
                    ", CPUs=", object$Specs[["CPUs"]],
                    ", Packages=NULL, Dyn.libs=NULL))\n\n", sep="")
               }
          else if({(Alg == "AM") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "AM") & Dim.Adapt & !Fast & !Ready}) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AM\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive, ", Periodicity=",
                    Rec.Periodicity, "))\n\n", sep="")
               }
          else if({(Alg == "AHMC") & !Dim.Adapt} |
               {(Alg == "AM") & !Dim.Adapt & !Fast & Ready} |
               {(Alg == "AMM") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "AMM") & Dim.Adapt & !Fast & !Ready} |
               {(Alg == "DRAM") & !Dim.Adapt & !Fast & Ready} |
               {(Alg == "DRM") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "DRM") & Dim.Adapt & !Fast & !Ready} |
               {(Alg == "RAM") & !Dim.Adapt & !Fast & !Ready} |
               {(Alg == "RAM") & !Dim.Adapt & Fast & !Ready} |
               {(Alg == "RWM") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "RWM") & Dim.Adapt & !Fast & !Ready}) {
               ### AMM
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- paste(oname, "$Specs$B", sep="")
               if(Rec.Status > Rec.Iterations)
                    Rec.Status <- Rec.Iterations
               n <- object$Iterations
               if(!is.null(object$Specs[["n"]]))
                    n <- object$Specs[["n"]] + object$Iterations
               w <- 0.05
               if(!is.null(object$Specs[["w"]]))
                    w <- object$Specs[["w"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AMM\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive,
                    ", B=", block, ", n=", n, ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity, ", w=", w,
                    "))\n\n", sep="")
               }
          else if({(Alg == "AM") & !Dim.Adapt & Fast & Ready} |
               {(Alg == "AM") & !Dim.Adapt & Fast & !Ready} |
               {(Alg == "AMM") & !Dim.Adapt & Fast & Ready} |
               {(Alg == "AMM") & !Dim.Adapt & Fast & !Ready} |
               {(Alg == "DRAM") & !Dim.Adapt & Fast & Ready} |
               {(Alg == "DRAM") & !Dim.Adapt & Fast & !Ready} |
               {(Alg == "MWG") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "MWG") & Dim.Adapt & !Fast & !Ready}) {
               ### AMWG
               if(Componentwise == 0) {
                    Rec.Iterations <- max(nrow(object$Posterior1),
                         trunc(Rec.Iterations / LIV))
                    Rec.Status <- trunc(interval(trunc(Rec.Status / LIV),
                         1, Rec.Iterations))
                    Rec.Thinning <- trunc(Rec.Iterations /
                         nrow(object$Posterior1))
                    Rec.Iterations <- nrow(object$Posterior1) *
                         Rec.Thinning}
               if(Rec.Periodicity > Rec.Iterations)
                    Rec.Periodicity <- max(trunc(Rec.Iterations * 0.01),1)
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- paste(oname, "$Specs$B", sep="")
               n <- object$Iterations
               if(!is.null(object$Specs[["n"]]))
                    n <- object$Specs[["n"]] + object$Iterations
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AMWG\", ",
                    "Specs=list(B=", block, ", n=", n, ", Periodicity=",
                    Rec.Periodicity, "))\n\n", sep="")
               }
          else if(Alg == "CHARM" & Acc.Rate.Level == 2) {
               ### CHARM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"CHARM\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          if(Alg == "CHARM" & Acc.Rate.Level != 2) {
               ### CHARM (Adaptive)
               al <- 0.44
               if(!is.na(object$Specs[["alpha.star"]]))
                    al <- object$Specs[["alpha.star"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"CHARM\", ",
                    "Specs=list(alpha.star=", al, "))\n\n", sep="")
               }
          else if(Alg == "DEMC" & Dim.Adapt & !Ready) {
               ### DEMC
               Nc <- max(3, round(LIV * 0.03))
               gamma <- "NULL"
               if(!is.null(object$Specs[["gamma"]]))
                    gamma <- object$Specs[["gamma"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"DEMC\", ",
                    "Specs=list(Nc=", Nc, ", Z=", oname, "$Posterior1, ",
                    "gamma=", gamma,
                    ", w=", object$Specs[["w"]],"))\n\n", sep="")
               }
          else if({(Alg == "DRAM") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "DRAM") & Dim.Adapt & !Fast & !Ready}) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"DRAM\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive, ", Periodicity=",
                    Rec.Periodicity, "))\n\n", sep="")
               }
          else if({(Alg == "DRM") & Dim.Adapt & Fast & Ready} |
               {(Alg == "DRM") & Dim.Adapt & !Fast & Ready}) {
               ### DRM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"DRM\", ",
                    "Specs=NULL))\n\n", sep="")
               }
          else if(Alg == "ESS") {
               ### ESS
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- "Block"
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"ESS\", Specs=list(B=", block,
                    "))\n\n", sep="")
               }
          else if(Alg == "GG") {
               ### GG
               Grid <- object$Specs[["Grid"]]
               dparm <- object$Specs[["dparm"]]
               CPUs <- object$Specs[["CPUs"]]
               Packages <- object$Specs[["Packages"]]
               Dyn.libs <- object$Specs[["Dyn.libs"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"GG\", Specs=", oname,
                    "$Specs)\n\n", sep="")
               }
          else if(Alg == "Gibbs") {
               ### Gibbs
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"Gibbs\", ",
                    "Specs=", oname, "$Specs)\n\n", sep="")
               }
          else if(Alg == "HARM") {
               ### HARM
               al <- 0.234
               if(is.na(object$Specs[["alpha.star"]]) |
                    !is.null(object$Specs[["alpha.star"]]))
                    al <- object$Specs[["alpha.star"]]
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- "Block"
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"HARM\", ",
                    "Specs=list(alpha.star=", al, ", B=", block, ")\n\n",
                    sep="")
               }
          else if({(Alg == "AHMC") & Dim.Adapt & Ready} |
               {(Alg == "HMC") & Dim.Adapt & Ready}) {
               ### HMC
               if(L > 1) {
                    L <- round(L*(Rec.Iterations/object$Iterations))
                    Rec.Iterations <- object$Iterations
                    Rec.Status <- object$Status
                    Rec.Thinning <- object$Thinning}
               m <- "NULL"
               if(!is.null(object$Specs[["m"]]) &
                    !identical(object$Specs[["m"]],list()))
                    m <- paste(oname, "$Specs$m", sep="")
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"HMC\", ",
                    "Specs=list(epsilon=", oname, "$CovarDHis[1,], ",
                    sep="")
               cat("L=", L, ", m=", m, "))\n\n", sep="")
               }
          else if(Alg == "HMCDA" & Dim.Adapt) {
               ### HMCDA
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"HMCDA\", ",
                    "Specs=list(A=", round(Rec.Iterations/2),
                    ", delta=", object$Specs[["delta"]], ",\n", sep="")
               cat("     epsilon=",
                    round(object$CovarDHis[nrow(object$CovarDHis),1],3),
                    ", Lmax=", object$Specs[["Lmax"]],
                    ", lambda=", object$Specs[["lambda"]],
                    "))\n\n", sep="")
               }
          else if(Alg == "IM") {
               ### IM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"IM\", ",
                    "Specs=list(mu=apply(", oname,
                    "$Posterior1, 2, mean)))\n\n", sep="")
               }
          else if(Alg == "INCA") {
               ### INCA
               detectedcores <- detectCores()
               cat(oname, " <- LaplacesDemon.hpc(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"INCA\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive, ", Periodicity=",
                    Rec.Periodicity, "),\n", sep="")
               cat("     Chains=",detectedcores,
                    ", CPUs=", detectedcores,
                    ", Packages=NULL, Dyn.libs=NULL)\n\n", sep="")
               }
          else if(Alg == "MALA") {
               ### MALA
               A <- object$Specs[["A"]]
               alpha.star <- object$Specs[["alpha.star"]]
               if(Dim.Adapt & Ready) gamma <- 0
               else gamma <- object$Specs[["gamma"]]
               delta <- object$Specs[["delta"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"MALA\", ",
                    "Specs=list(A=", A, ", alpha.star=", alpha.star,
                    ",\n", sep="")
               cat("          gamma=", gamma,
                    ", delta=", delta, ", epsilon=c(1e-6, 1e-7)))\n\n",
                    sep="")
               }
          else if(Alg == "MCMCMC") {
               ### MCMCMC
               lambda <- object$Specs[["lambda"]]
               CPUs <- object$Specs[["CPUs"]]
               Packages <- object$Specs[["Packages"]]
               Dyn.libs <- object$Specs[["Dyn.libs"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname,"$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"MCMCMC\", Specs=list(lambda=",
                    lambda, ", CPUs=", CPUs,
                   ", Packages=NULL, Dyn.libs=NULL))\n\n", sep="")
               }
          else if(Alg == "MTM") {
               ### MTM
               K <- object$Specs[["K"]]
               CPUs <- object$Specs[["CPUs"]]
               Packages <- object$Specs[["Packages"]]
               Dyn.libs <- object$Specs[["Dyn.libs"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname,"$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"MTM\", Specs=list(K=", K,
                    ", CPUs=", CPUs,
                    ", Packages=NULL, Dyn.libs=NULL))\n\n", sep="")
               }
          else if({(Alg == "ADMG") & Dim.Adapt & Fast & Ready} |
               {(Alg == "ADMG") & Dim.Adapt & !Fast & Ready} |
               {(Alg == "AMWG") & Dim.Adapt & Fast & Ready} |
               {(Alg == "AMWG") & Dim.Adapt & !Fast & Ready} |
               {(Alg == "MWG") & Dim.Adapt & Fast & Ready} |
               {(Alg == "MWG") & Dim.Adapt & !Fast & Ready}) {
               ### MWG
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- paste(oname, "$Specs$B", sep="")
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"MWG\", ",
                    "Specs=list(B=", block, "))\n\n", sep="")
               }
          else if({Alg == "NUTS"} |
               {(Alg == "HMCDA") & !Dim.Adapt}) {
               ### NUTS
               if(Alg == "HMCDA") delta <- 0.6
               else delta <- object$Specs[["delta"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"NUTS\", ",
                    "Specs=list(A=", round(Rec.Iterations/2),
                    ", delta=", delta, ",\n", sep="")
               cat("     epsilon=",
                    max(round(object$CovarDHis[nrow(object$CovarDHis),1],3),
                    1e-10), ", Lmax=", object$Specs[["Lmax"]],
                    "))\n\n", sep="")
               }
          else if(Alg == "OHSS") {
               ### OHSS
               if(Ready == TRUE) A <- 0
               else A <- object$Specs[["A"]]
               n <- object$Specs[["n"]] + object$Iterations
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"OHSS\", ",
                    "Specs=list(A=", A, ", n=", n,"))\n\n", sep="")
               }
          else if(Alg == "pCN") {
               ### pCN
               beta <- object$Specs[["beta"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RWM\", ",
                    "Specs=list(beta=", beta, "))\n\n", sep="")
               }
          else if({(Alg == "AM") & !Dim.Adapt & !Fast & !Ready} | 
               {(Alg == "AMM") & !Dim.Adapt & !Fast & Ready} | 
               {(Alg == "AMM") & !Dim.Adapt & !Fast & !Ready} | 
               {(Alg == "DRAM") & !Dim.Adapt & !Fast & !Ready} |
               {(Alg == "RAM") & Dim.Adapt & Fast & !Ready} |
               {(Alg == "RAM") & Dim.Adapt & !Fast & !Ready}) {
               ### RAM
               al <- 0.234
               if(!is.null(object$Specs[["alpha.star"]]))
                    al <- object$Specs[["alpha.star"]]
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- paste(oname, "$Specs$B", sep="")
               Dist <- "N"
               if(!is.null(object$Specs[["Dist"]]))
                    Dist <- object$Specs[["Dist"]]
               gamma <- 0.66
               if(!is.null(object$Specs[["gamma"]]))
                    gamma <- object$Specs[["gamma"]]
               n <- object$Iterations
               if(!is.null(object$Specs[["n"]]))
                    n <- object$Specs[["n"]] + object$Iterations
               if(Rec.Status > Rec.Iterations) Rec.Status <- Rec.Iterations
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RAM\", ",
                    "Specs=list(alpha.star=", al, ", B=", block,
                    ", Dist=\"", Dist, "\",\n", sep="")
               cat("     gamma=", gamma, ", n=", n, "))\n\n", sep="")
               }
          else if(Alg == "RDMH") {
               ### RDMH
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RDMH\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          else if(Alg == "Refractive") {
               ### Refractive
               Adaptive <- object$Specs[["Adaptive"]]
               m <- object$Specs[["m"]]
               w <- object$Specs[["w"]]
               if(Adaptive < object$Iterations)
                    w <- object$CovarDHis[nrow(object$CovarDHis),1]
               if(!Dim.Adapt) Adaptive <- 1
               else Adaptive <- Rec.Iterations + 1
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"Refractive\", ",
                    "Specs=list(Adaptive=", Adaptive, ", m=", m,
                    ", w=", w, ", r=1.3))\n\n", sep="")
               }
          else if(Alg == "RJ" & (Acc.Rate.Level > 1)) {
               ### RJ
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RJ\", ",
                    "Specs=list(bin.n=", object$Specs[["bin.n"]],
                    ", bin.p=", object$Specs[["bin.p"]], ",\n", sep="")
               cat("     parm.p=", paste("c(",
                    paste(object$Specs[["parm.p"]], collapse=","),
                    ")", sep=""), ",\n", sep="")
               cat("     selectable=", paste("c(",
                    paste(object$Specs[["selectable"]], collapse=","),
                    ")", sep=""), ",\n", sep="")
               cat("     selected=1*(Initial.Values != 0)))\n\n", sep="")
               }
          else if(Alg == "RSS") {
               ### RSS
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RSS\", ",
                    "Specs=", oname, "$Specs)\n\n", sep="")
               }
          else if({(Alg == "AM") & Dim.Adapt & Fast & Ready} |
               {(Alg == "AM") & Dim.Adapt & !Fast & Ready} |
               {(Alg == "AMM") & Dim.Adapt & Fast & Ready} |
               {(Alg == "AMM") & Dim.Adapt & !Fast & Ready} |
               {(Alg == "DRAM") & Dim.Adapt & Fast & Ready} |
               {(Alg == "DRAM") & Dim.Adapt & !Fast & Ready} |
               {(Alg == "RAM") & !Dim.Adapt & !Fast & Ready} |
               {(Alg == "RAM") & !Dim.Adapt & Fast & Ready} |
               {(Alg == "RAM") & Dim.Adapt & Fast & Ready} |
               {(Alg == "RAM") & Dim.Adapt & !Fast & Ready} |
               {(Alg == "RWM") & Dim.Adapt & Fast & Ready} |
               {(Alg == "RWM") & Dim.Adapt & !Fast & Ready}) {
               ### RWM
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]], list()))
                    block <- "Block"
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RWM\", ",
                    "Specs=list(B=", block, "))\n\n", sep="")
               }
          else if({(Alg == "DEMC") & Dim.Adapt & Fast & Ready} |
               {(Alg == "DEMC") & Dim.Adapt & !Fast & Ready}) {
               ### RWM from DEMC
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=cov(", oname, "$Posterior2), Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RWM\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          else if({(Alg == "SAMWG") & !Ready} |
               {(Alg == "SMWG") & !Ready}) {
               ### SAMWG
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"SAMWG\", ",
                    "Specs=list(Dyn=Dyn, Periodicity=", Rec.Periodicity,
                    "))\n\n", sep="")
               }
          else if({(Alg == "SAMWG") & Ready} |
               {(Alg == "SMWG") & Ready}) {
               ### SMWG
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"SMWG\", ",
                    "Specs=list(Dyn=Dyn))\n\n", sep="")
               }
          else if(Alg == "SGLD") {
               ### SGLD
               eps <- object$Specs[["epsilon"]]
               Fi <- object$Specs[["file"]]
               nr <- object$Specs[["Nr"]]
               nc <- object$Specs[["Nc"]]
               size <- object$Specs[["size"]]
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"SGLD\", ",
                    "Specs=", oname, "$Specs)\n\n", sep="")
               }
          else if(Alg == "Slice") {
               ### Slice
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- paste(oname, "$Specs$B", sep="")
               Bounds <- "c(-Inf,Inf)"
               if(!is.null(object$Specs[["Bounds"]]) &
                    !identical(object$Specs[["Bounds"]],list()))
                    Bounds <- paste(oname, "$Specs$Bounds", sep="")
               m <- "Inf"
               if(!is.null(object$Specs[["m"]]) &
                    !identical(object$Specs[["m"]],list()))
                    m <- paste(oname, "$Specs$m", sep="")
               Type <- "\"Continuous\""
               if(!is.null(object$Specs[["Type"]]) &
                    !identical(object$Specs[["Type"]],list()))
                    Type <- paste(oname, "$Specs$Type", sep="")
               w <- 1
               if(!is.null(object$Specs[["w"]]) &
                    !identical(object$Specs[["w"]],list()))
                    w <- paste(oname, "$Specs$w", sep="")
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"Slice\", ",
                    "Specs=list(B=", block, ", Bounds=", Bounds, ",\n",
                    sep="")
               cat("     m=", m, ", Type=", Type, ", w=", w, "))\n\n",
                    sep="")
               }
          else if(Alg == "THMC") {
               ### THMC
               if(L > 1) {
                    L <- round(L*(Rec.Iterations/object$Iterations))
                    Rec.Iterations <- object$Iterations
                    Rec.Status <- object$Status
                    Rec.Thinning <- object$Thinning}
               m <- "NULL"
               if(!is.null(object$Specs[["m"]]) &
                    !identical(object$Specs[["m"]],list()))
                    m <- paste(oname, "$Specs$m", sep="")
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"THMC\", ",
                    "Specs=list(epsilon=", oname, "$CovarDHis[1,],\n",
                    sep="")
               cat("     L=", L, ", m=", m,
                    ", Temperature=", object$Specs[["Temperature"]],
                    "))\n\n", sep="")
               }
          else if(Alg == "t-walk" |
               {(Alg == "DEMC") & !Dim.Adapt}) {
               ### twalk
               if(Alg == "DEMC") {
                    n1 <- 4
                    at <- 6
                    aw <- 1.5
                    }
               else {
                    n1 <- object$Specs[["n1"]]
                    at <- object$Specs[["at"]]
                    aw <- object$Specs[["aw"]]}
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"twalk\", ",
                    "Specs=list(SIV=NULL, n1=", n1,
                    ", at=", at, ", aw=", aw, "))\n\n", sep="")
               }
          else if(Alg == "UESS") {
               ### UESS
               if(Ready == TRUE) A <- 0
               else A <- object$Specs[["A"]]
               block <- "NULL"
               if(!is.null(object$Specs[["B"]]) &
                    !identical(object$Specs[["B"]],list()))
                    block <- "Block"
               m <- object$Specs[["m"]]
               n <- object$Specs[["n"]] + object$Iterations
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"UESS\", ",
                    "Specs=list(A=", A, ", B=", block, ", m=", m,
                    ", n=", n, "))\n\n", sep="")
               }
          else if((Alg == "USAMWG") | (Alg == "USMWG")) {
               ### USAMWG or USMWG
               cat("A Demonic Suggestion will not be made.\n\n")
               }
          }
     cat("Laplace's Demon is finished consorting.\n")
     }

#End
