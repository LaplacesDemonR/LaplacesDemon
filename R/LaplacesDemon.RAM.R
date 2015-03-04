###########################################################################
# LaplacesDemon.RAM                                                       #
#                                                                         #
# The purpose of the LaplacesDemon.RAM function is to estimate the RAM    #
# required to update a given model and data in LaplacesDemon.             #
###########################################################################

LaplacesDemon.RAM <- function(Model, Data, Iterations, Thinning,
     Algorithm="RWM")
     {
     if(missing(Model))
          stop("The Model argument is required.")
     if(missing(Data))
          stop("The Data argument is required.")
     Const <- 1048600
     LIV <- length(Data[["parm.names"]])
     LM <- length(Data[["mon.names"]])
     Covar <- 0
     if(Algorithm %in% c("ADMG","AFSS","AM","AMM","DRAM","DRM","ESS","IM",
          "INCA","MALA","OHSS","RWM","RAM","UESS")) {
          ### Covariance is required
          Covar <- Covar + as.vector(object.size(matrix(runif(LIV*LIV),
          LIV, LIV))) / Const
          }
     else if(Algorithm %in% c("AGG","AM","AMM","AMWG","DRAM","DRM","INCA",
          "MWG","RWM","SAMWG","SMWG","USAMWG","USMWG")) {
          ### Variance is required
          Covar <- Covar + as.vector(object.size(runif(LIV))) / Const}
     Data <- as.vector(object.size(Data)) / Const
     Deviance <- as.vector(object.size(runif(round(Iterations /
          Thinning)))) / Const
     Initial.Values <- as.vector(object.size(runif(LIV))) / Const
     Model <- as.vector(object.size(Model)) / Const
     Monitor <- as.vector(object.size(matrix(runif(Iterations*LM),
          round(Iterations / Thinning), LM))) / Const
     post <- 0
     if(Algorithm %in% c("AHMC","AM","DRAM","INCA","NUTS","OHSS"))
          post <- as.vector(object.size(matrix(runif(Iterations*LIV),
               Iterations, LIV))) / Const
     Posterior1 <- as.vector(object.size(matrix(runif(round(Iterations /
          Thinning)), round(Iterations / Thinning), LIV))) / Const
     Posterior2 <- as.vector(object.size(matrix(runif(round(Iterations /
          Thinning)), round(Iterations / Thinning), LIV))) / Const
     Summary1 <- as.vector(object.size(matrix(runif((LIV+1+LM)*7),
          LIV+1+LM, 7))) / Const
     Summary2 <- as.vector(object.size(matrix(runif((LIV+1+LM)*7),
          LIV+1+LM, 7))) / Const
     mem.list <- list(Covar=Covar,
          Data=Data,
          Deviance=Deviance,
          Initial.Values=Initial.Values,
          Model=Model,
          Monitor=Monitor,
          post=post,
          Posterior1=Posterior1,
          Posterior2=Posterior2,
          Summary1=Summary1,
          Summary2=Summary2,
          Total=sum(Covar,Data,Deviance,Initial.Values,Model,Monitor,
               post,Posterior1,Posterior2,Summary1,Summary2))
     return(mem.list)
     }

#End
