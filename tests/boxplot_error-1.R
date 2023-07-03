#### Test 3: square domain ####
#            locations in nodes
#            PDE
#            no covariates
#            no BC
#            order FE = 2
#            separate penalizations
library(fdaPDE)
library(microbenchmark)

rm(list=ls())
graphics.off()

# x = seq(0,1, length.out = 41)
# y = x
# locations = expand.grid(x,y)
# 
# mesh = create.mesh.2D(locations)
# plot(mesh)
# 
# nnodes = dim(mesh$nodes)[1]

#FEMbasis = create.FEM.basis(mesh)

# Test function
f = function(x, y, z)
{
  a1=1
  a2=4
  (a1*sin(2*pi*x)*cos(2*pi*y)+a2*sin(3*pi*x))*cos(pi*z)
}

seed_numbers = seq(from = 0, to = 10000, by = 1000)
flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
#flags = rbind(c(FALSE, FALSE),c(TRUE, FALSE),c(TRUE, TRUE))


nRef = seq(from = 1,to = 3, by = 1)

dim1 = dim(flags)[1]
dim2 = length(seed_numbers)
dim3 = length(nRef)

#RMSE_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))
RMSE_matrix <- array(0, dim = c(dim1, dim2, dim3))
#time_matrix <- array(0, nrow = dim(flags)[1], ncol = length(nRef))
time_matrix <- array(0, dim = c(dim1, dim2, dim3))
for (i in seq(from = 1,to = dim(flags)[1], by = 1)) { 
  h = 0
  for (n in nRef) {
    h = h + 1
    k = 0
    
    x = seq(0,1, length.out = (10*n)+1)
    y = x
    locations = expand.grid(x,y)
    
    mesh = create.mesh.2D(locations)
    plot(mesh)
    
    nnodes = dim(mesh$nodes)[1]
    
    FEMbasis = create.FEM.basis(mesh)
    
    for (seed_number in seed_numbers) {
      k = k + 1
      #if(i == 1)
      NumTimePoints=11
      #else
      #  NumTimePoints=13
      TimePoints=seq(0,2,length.out = NumTimePoints)
      
      SpacePoints=mesh$nodes
      SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))
      
      # Exact solution (pointwise at nodes)
      sol_exact=f(SpaceTimePoints[,1],SpaceTimePoints[,2],SpaceTimePoints[,3])
      
      set.seed(seed_number)
      ran = range(sol_exact)
      data = sol_exact + rnorm(length(sol_exact), mean=0, sd=0.05*abs(ran[2]-ran[1]))
      observations=matrix(data,nrow(SpacePoints),NumTimePoints)
    
      # Set PDE parameters
      PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)
      
      # Set smoothing parameter
      lambdaS = 1e-4
      lambdaT = 1e-4
      
      #### Test 3.1: Without GCV
      results <- microbenchmark(output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                                                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                                                            PDE_parameters=PDE_parameters,
                                                            #FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2],
                                FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2], FLAG_FINITEDIFFERENCES = flags[i,3]   )
                                ,times = 1)
      
    
      
      # RMSE evaluation on a fine grid
      xeval=runif(1000,0,1)  # quadrato con x tra 0 e 1
      yeval=runif(1000,0,1)  # quadrato con y tra 0 e 1
      teval=runif(1000,0,2)
      sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval))
      sol_exact=f(xeval,yeval,teval)
      RMSE<-function(f,g) sqrt(mean((f-g)^2))
      RMSE_matrix[i, k, h] = RMSE(sol_eval,sol_exact)
      time_matrix[i, k, h] = results$time/10^6 #ms
      
    }
  }

}

quartz()
labels <- c("B-SPLINES", "PARABOLIC", "PARABOLIC ITER", "FIN DIFF", "FIN DIFF ITER")
boxplot(data.frame(t(RMSE_matrix[,,2])), las=1, col='gold', names = labels)
title(main = "Square domain")
#boxplot(data.frame(t(time_matrix)), las=2, col='gold')

mediana_tempi <- apply(time_matrix, c(1, 3), median)
mediana_tempi_mat <- array(mediana_tempi, dim = c(5, 1, 3))

quartz()
par(mfrow = c(1, 1))
colors <- seq(from = 1, to = dim(flags)[1], by = 1)

# Ciclo for per il plottaggio dei dati
for (i in seq(from = 1, to = dim(flags)[1], by = 1)) {
  if (i == 1) {
    plot(nRef*10+1, mediana_tempi_mat[i,,], type = "l", ylim = c(10, 500000), col = colors[i],log = "y", main = "Computational costs", xlab = "#nodes", ylab = "time [ms]", lwd = 2)
  } else {
    lines(nRef*10+1, mediana_tempi_mat[i,,], col = colors[i], log = "y", lwd = 2)
  }
}
legend_vector = c("B-SPLINE", "PARABOLICO MONOLITICO", "PARABOLICO ITERATIVO", "DIFFERENZE FINITE MONOLITICO", "DIFFERENZE FINITE ITERATIVO")
legend("topleft", legend = legend_vector , col = colors, lty = 1, lwd = 2)

################################################################################################
################################################################################################
#### Test 3: square domain ####
#            locations in nodes
#            PDE
#            no covariates
#            no BC
#            order FE = 2
#            separate penalizations
library(fdaPDE)
library(microbenchmark)

rm(list=ls())
graphics.off()

# Test function
f = function(x, y, z)
{
  a1=1
  a2=4
  (a1*sin(2*pi*x)*cos(2*pi*y)+a2*sin(3*pi*x))*cos(pi*z)
}

#flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
flags = rbind(c(FALSE, FALSE),c(TRUE, FALSE),c(TRUE, TRUE))

nRef = seq(from = 1,to = 5, by = 1)

time_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(nRef))
seed_number = 1000000000


for (i in seq(from = 1,to = dim(flags)[1], by = 1)) { 
  k = 0
  for (n in nRef) {
    k = k + 1
    
    x = seq(0,1, length.out = (10*n)+1)
    y = x
    locations = expand.grid(x,y)
    
    mesh = create.mesh.2D(locations)
    plot(mesh)
    
    nnodes = dim(mesh$nodes)[1]
    
    FEMbasis = create.FEM.basis(mesh)
    
    
    NumTimePoints=11
    TimePoints=seq(0,2,length.out = NumTimePoints)
    
    SpacePoints=mesh$nodes
    SpaceTimePoints=cbind(rep(mesh$nodes[,1],NumTimePoints),rep(mesh$nodes[,2],NumTimePoints),rep(TimePoints,each=nrow(mesh$nodes)))
    
    # Exact solution (pointwise at nodes)
    sol_exact=f(SpaceTimePoints[,1],SpaceTimePoints[,2],SpaceTimePoints[,3])
    
    set.seed(seed_number)
    ran = range(sol_exact)
    data = sol_exact + rnorm(length(sol_exact), mean=0, sd=0.05*abs(ran[2]-ran[1]))
    observations=matrix(data,nrow(SpacePoints),NumTimePoints)
    
    # Set PDE parameters
    PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)
    
    # Set smoothing parameter
    lambdaS = 1e-4
    lambdaT = 1e-4
    

    
    results <- microbenchmark(output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                                                          FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                                                          PDE_parameters=PDE_parameters,
                                                          FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2])
                              #FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2], FLAG_FINITEDIFFERENCES = flags[i,3]   )
                              , times = 1)
    
    time_matrix[i, k] = results$time/10^6 #ms
  }
}
time_matrix

# Creazione della finestra grafica
quartz()
par(mfrow = c(1, 1))
colors <- seq(from = 1, to = dim(flags)[1], by = 1)

# Ciclo for per il plottaggio dei dati
for (i in seq(from = 1, to = dim(flags)[1], by = 1)) {
  if (i == 1) {
    plot(nRef, time_matrix[i,], type = "l",col = colors[i],log = "y", main = "Costi computazionali", xlab = "#nodes", ylab = "time [ms]", lwd = 2)
  } else {
    lines(nRef, time_matrix[i,], col = colors[i], log = "y", lwd = 2)
  }
}

legend_vector = c("B-SPLINE", "PARABOLICO MONOLITICO", "PARABOLICO ITERATIVO")
legend("topleft", legend = legend_vector , col = colors, lty = 1)


#### Test 2: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1
#            separate penalizations 
library(fdaPDE)
library(microbenchmark)

rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

# plot spatial locations
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis = create.FEM.basis(mesh)

f<-function(x,y,t)
{
  K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
  res=numeric(length =length(x))
  for(i in 1:length(x))
  {
    if(x[i]>=0 && y[i]>0)
      res[i]=cos(t[i])*(0.25*pi+x[i])+(y[i]-0.5)^2
    if(x[i]>=0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-0.25*pi-x[i])+(-y[i]-0.5)^2
    if(x[i]<0 && y[i]>0)
      res[i]=cos(t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    if(x[i]<0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
  }
  res
}


seed_numbers = seq(from = 0, to = 30000, by = 1000)

flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
#flags = rbind(c(FALSE, FALSE),c(TRUE, FALSE),c(TRUE, TRUE))
RMSE_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))
time_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))

for (i in seq(from = 1,to = dim(flags)[1], by = 1)) { 
  k = 0
  for (seed_number in seed_numbers) {
    k = k + 1
    NumTimeInstants=5
    TimePoints=seq(0,pi,length.out =NumTimeInstants)
    nnodes = nrow(mesh$nodes)
    
    space_time_locations = cbind(rep(TimePoints,each=nrow(locations)),rep(locations[,1],NumTimeInstants),rep(locations[,2],NumTimeInstants))
    sol_exact = f(space_time_locations[,2],space_time_locations[,3],space_time_locations[,1])
    
    ndata = length(sol_exact)
    
    # Create covariates
    set.seed(509875)
    cov1 = rnorm(ndata, mean = 1, sd = 2)
    
    # Add error to simulate data
    set.seed(seed_number)
    data = sol_exact + 2*cov1 
    data = data + rnorm(ndata, mean = 0, sd =  0.05*diff(range(sol_exact)))
    observations = matrix(data,nrow(locations),NumTimeInstants)
    
    # Set smoothing parameter
    lambdaS = 10^-2
    lambdaT = 10^-2.5
    
    results <- microbenchmark(    output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                                                              observations=observations, 
                                                              covariates = cov1,
                                                              FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                                  FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2], FLAG_FINITEDIFFERENCES = flags[i,3]   )
                              , times = 1)
    
    RMSE<-function(f,g) sqrt(mean((f-g)^2))
    sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = space_time_locations, lambdaS=lambdaS, lambdaT=lambdaT)
    
    RMSE_matrix[i,k] = RMSE(sol_eval,sol_exact)
    time_matrix[i, k] = results$time/10^6 #ms
  }
}


quartz()
labels <- c("B-SPLINES", "PARABOLIC", "PARABOLIC ITER", "FIN DIFF", "FIN DIFF ITER")
boxplot(data.frame(t(RMSE_matrix)), las=1, col='gold',  names = labels)
title(main = "C-shaped domain")

quartz()
labels <- c("B-SPLINES", "PARABOLIC", "PARABOLIC ITER", "FIN DIFF", "FIN DIFF ITER")
boxplot(data.frame(t(time_matrix)), las=1, col='gold',  names = labels)
title(main = "C-shaped domain")




################################ 2.5D ###############################################
# hub pointwise (examples with and without covariates) #
library(fdaPDE)
library(microbenchmark)
rm(list=ls())

data(hub2.5D)
mesh <- create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes,triangles = hub2.5D$hub2.5D.triangles)
FEMbasis <- create.FEM.basis(mesh)

plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

# Locations at nodes
nodesLocations=mesh$nodes

lambdaS=10^seq(from = -3,to = -1, by = 0.5) #-2.5 meglio
lambdaT=10^seq(from = -4,to = -2, by = 0.5) #2 meglio

## SPLINE lambdaS = -3 lambdaT = -1.5
## MON lambdaS = -2.5 lambdaT = -3
## ITER lambdaS = -2.5 lambdaT = -3

seed_numbers = seq(from = 0, to = 10000, by = 1000)
flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
#flags = rbind(c(FALSE, FALSE),c(TRUE, FALSE),c(TRUE, TRUE))

RMSE_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))
time_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))


for(i in seq(from = 1,to = dim(flags)[1], by = 1)){
  k = 0
  for (seed_number in seed_numbers) {
      k = k + 1
    
      # Exact data - Locations at nodes
      set.seed(seed_number)
      
      nnodes = nrow(mesh$nodes)
      a1 = rnorm(1,mean = 1, sd = 0.1)
      a2 = rnorm(1,mean = 1, sd = 0.1)
      a3 = rnorm(1,mean = 1, sd = 0.1)
      TimeNodes = 0:4
      
      locations = cbind(rep(TimeNodes,each=nnodes),rep(nodesLocations[,1],length(TimeNodes)),rep(nodesLocations[,2],length(TimeNodes)),rep(nodesLocations[,3],length(TimeNodes)))
      
      func = function(x)
      {
        (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
      }
      
      func_evaluation = func(locations)
      # Plot the exact solution
      #plot(FEM.time(coeff=array(func_evaluation,dim=c(length(func_evaluation),1,1)), FEMbasis = FEMbasis,time_mesh=TimeNodes,FLAG_PARABOLIC=T),t = 2)
      
      
      if(i==1){
        lambdaS = 10^(-3)
        lambdaT = 10^(-1.5)
      }
      else{
        lambdaS = 10^(-2.5)
        lambdaT = 10^(-3)
      }
      
      
      lambdaS_par=10^seq(-4, -3, 0.25)
      lambdaT_par=10^seq(1, 1.8, 0.2)
      
      cov1=4*sin(2*pi*locations[,2])*cos(2*pi*locations[,3])
      cov2=rnorm(nnodes*length(TimeNodes), mean=3, sd=0.1)*rep(exp(-TimeNodes/length(TimeNodes)),each=nnodes)
      W=cbind(cov1,cov2)
      
      
      # Fix betas
      beta_exact=c(0.45,0.3)
      
      ran=range(func_evaluation)
      
      data = func_evaluation +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))
      
      ran = range(func_evaluation+ W%*%beta_exact)
      datacov=func_evaluation+ W%*%beta_exact +rnorm(nnodes*length(TimeNodes),mean=0,sd=0.05*(ran[2]-ran[1]))
      
      data = matrix(data,nrow(mesh$nodes),length(TimeNodes))
      datacov = matrix(datacov,nrow(mesh$nodes),length(TimeNodes))
      
  
  
  
      results <- microbenchmark(output_CPP<-smooth.FEM.time(time_mesh = TimeNodes, observations= data,
                                  FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                                  FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2], FLAG_FINITEDIFFERENCES = flags[i,3] )
                                  ,times = 1)
      
      # Example of RMSE computation
      TimeNodesEval=seq(0,4,length.out = 9)
      eval_locations = cbind(rep(TimeNodesEval,each=nnodes),rep(nodesLocations[,1],length(TimeNodesEval)),rep(nodesLocations[,2],length(TimeNodesEval)),rep(nodesLocations[,3],length(TimeNodesEval)))
      sol_exact = func(eval_locations)
      RMSE<-function(f,g) sqrt(mean((f-g)^2))
      
      sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,locations = nodesLocations,time.instants = TimeNodesEval, lambdaS = output_CPP$optimization$lambda_position[1], lambdaT = output_CPP$optimization$lambda_position[2])
      RMSE_matrix[i,k] = RMSE(sol_eval,sol_exact)
      time_matrix[i,k] = results$time/10^6 #ms
      
  }
}



mediana_lambda <- apply(RMSE_matrix, c(1, 3), median)
mediana_lambda_mat <- array(mediana_lambda, dim = c(5, 5, 1))

mediana_lambdaS <- apply(mediana_lambda_mat, MARGIN=1, min)
mediana_lambdaS_mat <- array(mediana_lambda, dim = c(5, 1, 1))

  
quartz()
labels <- c("B-SPLINES", "PARABOLIC", "PARABOLIC ITER", "FIN DIFF", "FIN DIFF ITER")
#x11()
boxplot(data.frame(t(RMSE_matrix)), las=1, col='gold', names = labels)
title(main = "Hub 2.5D")
boxplot(data.frame(t(time_matrix)), las=2, col='gold')






######### 3D (These tests are slow!) #########
##### ATTUALMENTE NON FUNZIONA MISERIACCIA
#### sphere 3D pointwise (with or without covariates + locations at nodes or not + stochastic GCV) ####
library(fdaPDE)
rm(list=ls())

# Build mesh: Sphere
data(sphere3Ddata)
sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
plot(sphere3D)
FEMbasis <- create.FEM.basis(sphere3D)
nodesLocations=sphere3D$nodes
nnodes = nrow(sphere3D$nodes)
TimeLocations = seq(0,1,length.out = 5)
Locations = cbind(rep(TimeLocations,each=nnodes),rep(nodesLocations[,1],length(TimeLocations)),rep(nodesLocations[,2],length(TimeLocations)),rep(nodesLocations[,3],length(TimeLocations)))

seed_numbers = seq(from = 0, to = 10000, by = 1000)
flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
#flags = rbind(c(FALSE, FALSE),c(TRUE, FALSE),c(TRUE, TRUE))

RMSE_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))


for (i in seq(from = 1,to = dim(flags)[1], by = 1)) {
  k = 0
  for (seed_number in seed_numbers) {
    k = k + 1
    
    # Exact data - Locations at nodes
    set.seed(seed_number)
    # Exact test function
    a1 = rnorm(1,mean = 1, sd = 0.3)
    set.seed(seed_number+1)
    a2 = rnorm(1,mean = 1, sd = 0.3)
    set.seed(seed_number+2)
    a3 = rnorm(1,mean = 1, sd = 0.3)
    set.seed(seed_number+3)
    a4 = rnorm(1,mean = 1, sd = 0.3)
    
    func = function(x)
    {
      a1*sin(2*pi*(x[,1]*x[,2]))+a2*cos(2*pi*x[,2])+a3*cos(2*pi*x[,3])+a4*sin(2*pi*x[,4])
    }
    
    func_evaluation = func(Locations)
    ran=range(func_evaluation)
    
    # plot(FEM(func_evaluation[1:nnodes],FEMbasis))
    # Set smoothing parameter
    # Generate locations
    nloc = 1000
    loc=matrix(data=runif(3*nloc, min=-1,max=1),nrow=nloc,ncol=3,byrow=T)
    
    ind=NULL
    for(row in 1:nloc){
      normvec = (loc[row,1]^2+loc[row,2]^2+loc[row,3]^2)
      if(normvec>0.975)   # check points outside the sphere and remove them
        ind = c(ind,row)
    }
    
    loc=loc[-ind,]
    nloc=dim(loc)[1]
    timeloc = seq(0,1,length.out=5)
    loc = cbind(rep(timeloc,each=nloc),rep(loc[,1],length(timeloc)),rep(loc[,2],length(timeloc)),rep(loc[,3],length(timeloc)))
    
    
    # Exact test function - locations different from nodes
    func_evaluation2=func(loc)
    
    
    cov1=(4*sin(2*pi*Locations[,2])+6*sin((2*pi*Locations[,3])^2))*(1-exp(-Locations[,1]))/3
    cov2=cos(-2*pi*Locations[,4])+2*Locations[,1]*sin(2*pi*Locations[,2])/6
    
    cov1_nonod=(4*sin(2*pi*loc[,2])+6*sin((2*pi*loc[,3])^2))*(1-exp(-loc[,1]))/3
    cov2_nonod=cos(-2*pi*loc[,4])+2*loc[,1]*sin(2*pi*loc[,2])/6
    
    W=cbind(cov1,cov2)
    W2=cbind(cov1_nonod,cov2_nonod)
    # 
    #     lambdaS=10^seq(-5.0, -4.0, 0.25)
    #     lambdaT=10^seq(-1.5, -0.5, 0.25)
    #     
    #     lambdaS2=10^seq(-5.5, -4.5, 0.25)
    #     lambdaT2=10^seq(-1.5, -0.5, 0.25)
    #     
    #     lambdaS_par=10^seq(-4.8, -4.4, 0.1)
    #     lambdaT_par=10^seq(1.4, 1.8, 0.1)
    #     
    #     lambdaS_par2=10^seq(-4.4, -4.0, 0.1)
    #     lambdaT_par2=10^seq(1.4, 1.8, 0.1)
    lambdaS = 10^-2
    lambdaT = 10^-3
    
    
    if(i == 4 || i == 5){
    lambdaS=10^(-2.5)
    lambdaT=10^(-1)}
    
    beta_exact= c(0.7,2.0)
    ran = range(func_evaluation)
    data = func_evaluation +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))
    
    ran = range(func_evaluation2)
    data_noloc = func_evaluation2 +rnorm(length(func_evaluation2),mean=0,sd=0.05*(ran[2]-ran[1]))
    
    ran = range(func_evaluation+ W%*%beta_exact)
    datacov=func_evaluation+ W%%beta_exact +rnorm(length(func_evaluation),mean=0,sd=0.05*(ran[2]-ran[1]))
    
    data = matrix(data,nnodes,length(TimeLocations))
    data_noloc = matrix(data_noloc,nloc,length(timeloc))
    datacov = matrix(datacov,nnodes,length(TimeLocations))
    
    
    
    output_CPP<-smooth.FEM.time(observations= data, 
                                FEMbasis=FEMbasis, time_mesh = TimeLocations,
                                lambdaS=lambdaS, lambdaT=lambdaT,
    FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2], FLAG_FINITEDIFFERENCES = flags[i,3] )
    
    ####RMSE computation
    #sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,locations = nodesLocations, time.instants = TimeLocations, lambdaS = output_CPP$optimization$lambda_position[1], lambdaT = output_CPP$optimization$lambda_position[1])
    sol_exact = func_evaluation
    
    eval_locations = cbind(rep(TimeLocations,each=nnodes),rep(nodesLocations[,1],length(TimeLocations)),rep(nodesLocations[,2],length(TimeLocations)),rep(nodesLocations[,3],length(TimeLocations)))
    
    f_nodi = output_CPP$solution$f
    
    sol_eval = NaN
    if (length(sol_exact)!= length(f_nodi)){
      sol_eval = f_nodi[(nnodes+1):(length(f_nodi)-nnodes)]}
    else{
      sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,locations = nodesLocations, time.instants = TimeLocations)
    }
    
    RMSE<-function(f,g) sqrt(mean((f-g)^2))
    RMSE_matrix[i, k] = RMSE(sol_eval,sol_exact)
  }
}

quartz()
#x11()
boxplot(data.frame(t(RMSE_matrix)), las=2, col='gold')
labels <- c("B-SPLINES", "FD MON", "FD MON PREC", "FD ITER", "FD ITER PREC", names = labels)