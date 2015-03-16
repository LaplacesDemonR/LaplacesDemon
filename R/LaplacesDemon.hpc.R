###########################################################################
# LaplacesDemon.hpc                                                       #
#                                                                         #
# The purpose of the LaplacesDemon.hpc function is to extend              #
# LaplacesDemon to parallel processing on multiple cores.                 #
###########################################################################

#I don't know how to 'shut up' the makeCluster function...
#no luck with capture.output, sink, or invisible

LaplacesDemon.hpc <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=10000, Status=100, Thinning=10, Algorithm="MWG",
     Specs=list(B=NULL), Debug=list(DB.chol=FALSE, DB.eigen=FALSE,
     DB.MCSE=FALSE, DB.Model=TRUE), LogFile="", Chains=2, CPUs=2,
     Type="PSOCK", Packages=NULL, Dyn.libs=NULL)
     {
     detectedCores <- max(detectCores(), as.integer(Sys.getenv("NSLOTS")),
          na.rm=TRUE)
     cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
          append=TRUE)
     if(CPUs > detectedCores) {
          cat("\nOnly", detectedCores, "will be used.\n", file=LogFile,
               append=TRUE)
          CPUs <- detectedCores}
     if(is.vector(Initial.Values)) {
          Initial.Values <- matrix(Initial.Values, Chains,
               length(Initial.Values), byrow=TRUE)
          cat("\nWarning: initial values were a vector, and are now a",
               file=LogFile, append=TRUE)
          cat("\n", Chains, "x", length(Initial.Values), "matrix.\n",
               file=LogFile, append=TRUE)}
     if(Algorithm == "INCA" && Chains != CPUs) {
          Chains <- CPUs
          cat("\nINCA:", Chains, "chains will be used\n")}
     cat("\nLaplace's Demon is preparing environments for CPUs...",
          file=LogFile, append=TRUE)
     cat("\n##################################################\n",
          file=LogFile, append=TRUE)
     cl <- makeCluster(CPUs, Type)
     cat("\n##################################################\n",
          file=LogFile, append=TRUE)
     on.exit({stopCluster(cl); cat("\n\nLaplace's Demon has finished.\n", file=LogFile, append=TRUE)})
     Packages <- c(Packages, "LaplacesDemon")
     varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
          ls(envir=parent.env(environment()))))
     clusterExport(cl, varlist=varlist, envir=environment())
     clusterSetRNGStream(cl)
     wd <- getwd()
     clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
          envir=environment())
     demon.wrapper <- function(x, ...)
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
          LaplacesDemon(Model, Data, Initial.Values[x,],
               Covar, Iterations, Status, Thinning,
               Algorithm, Specs, Debug,
               LogFile=paste(LogFile, ".", x, sep=""))
          }
     cat("\nStatus messages are not displayed for parallel processing.",
          file=LogFile, append=TRUE)
     cat("\nLaplace's Demon is beginning parallelization...\n",
          file=LogFile, append=TRUE)
     if(Algorithm == "INCA") {
          ### Start hpc server
          system(paste("Rscript -e 'library(parallel);library(LaplacesDemon);server_Listening(n=",CPUs,")'", sep=""), wait=FALSE)
          cat("Start hpc server...\n", file=LogFile, append=TRUE)
          ### Export chain number
          clusterExport(cl, varlist="Chains", envir=environment())
          ### Connect each process to server_Listening with 0.5s time delay
          clusterEvalQ(cl, con <- NULL)
          doCon <- function(i) {
               Sys.sleep(i/2)
               con <<- socketConnection("localhost", 19009, blocking=TRUE,
                    open="r+")}
          clusterExport(cl, varlist="doCon", envir=environment())
          expr <- NULL
          for (i in 1:CPUs) {
               tmp <- parse(text=paste("doCon(", i,")", sep=""))
               expr <- c(expr, tmp)}
          clusterApply(cl, expr, eval, env=.GlobalEnv)
          cat("\nOpen connections to hpc server...", file=LogFile,
               append=TRUE)}
     LaplacesDemon.out <- clusterApply(cl, 1:Chains, demon.wrapper,
          Model, Data, Initial.Values, Covar, Iterations, Status,
          Thinning, Algorithm, Specs, Debug)
     class(LaplacesDemon.out) <- "demonoid.hpc"
     if(Algorithm == "INCA") {
          ### Stop server_Listening
          clusterEvalQ(cl, {close(con)})
          cat("\nClose connections to hpc server...", file=LogFile,
               append=TRUE)
          }
     return(LaplacesDemon.out)
     }

#End
