###########################################################################
# PMC.RAM                                                                 #
#                                                                         #
# The purpose of the PMC.RAM function is to estimate the RAM required to  #
# update a given model and data in PMC.                                   #
###########################################################################

PMC.RAM <- function(Model, Data, Iterations, Thinning, M, N)
     {
     if(missing(Model))
          stop("The Model argument is required.")
     if(missing(Data))
          stop("The Data argument is required.")
     Const <- 1048600
     LIV <- length(Data[["parm.names"]])
     LM <- length(Data[["mon.names"]])
     alpha <- as.vector(object.size(matrix(rep(1/M, M), M, Iterations))) /
          Const
     Covar <- as.vector(object.size(array(0,
          dim=c(LIV,LIV,Iterations,M)))) / Const
     Data <- as.vector(object.size(Data)) / Const
     Deviance <- as.vector(object.size(rep(0,N))) / Const
     Initial.Values <- as.vector(object.size(matrix(0, M,
          length(Data[["parm.names"]])))) / Const
     LH <- as.vector(object.size(array(0, dim=c(N, Iterations, M)))) /
          Const
     LP <- as.vector(object.size(array(0, dim=c(N, Iterations, M)))) /
          Const
     Model <- as.vector(object.size(Model)) / Const
     Monitor <- as.vector(object.size(matrix(runif(N*LM), N, LM))) / Const
     Mu <- as.vector(object.size(array(0, dim=c(Iterations, LIV, M)))) /
          Const                  
     Posterior1 <- as.vector(object.size(array(0,
          dim=c(N, LIV, Iterations, M)))) / Const
     Posterior2 <- Posterior1[,,Iterations,1]
     Posterior2 <- as.vector(object.size(Posterior2)) / Const
     #Note: Posterior2 gets thinned, but at one point it's this large.
     Summary <- as.vector(object.size(matrix(0, LIV+1+LM, 7))) / Const
     W <- as.vector(object.size(matrix(0, N, Iterations))) / Const
     mem.list <- list(alpha=alpha,
          Covar=Covar,
          Data=Data,
          Deviance=Deviance,
          Initial.Values=Initial.Values,
          LH=LH,
          LP=LP,
          Model=Model,
          Monitor=Monitor, 
          Mu=Mu,
          Posterior1=Posterior1,
          Posterior2=Posterior2,
          Summary=Summary,
          W=W,
          Total=sum(alpha,Covar,Data,Deviance,Initial.Values,LH,LP,
               Model,Monitor,Mu,Posterior1,Posterior2,Summary,W))
     return(mem.list)
     }

#End
