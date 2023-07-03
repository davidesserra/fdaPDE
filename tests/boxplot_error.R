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

seed_numbers = seq(from = 1000,to = 10000, by = 2500)
flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
RMSE_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))
nRef = seq(from = 1,to = 4, by = 1)

for (i in seq(from = 1,to = dim(flags)[1], by = 1)) { 
  k = 0
  for (seed_number in seed_numbers) {
    for (n in nRef) {
      x = seq(0,1, length.out = (5*n)+1)
      y = x
      locations = expand.grid(x,y)
      
      mesh = create.mesh.2D(locations)
      plot(mesh)
      
      nnodes = dim(mesh$nodes)[1]
      
      FEMbasis = create.FEM.basis(mesh)
    
      k = k + 1
      if(k == 1)
        NumTimePoints=11
      else
        NumTimePoints=13
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
      output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                                  FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                                  PDE_parameters=PDE_parameters,
                                  FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2], FLAG_FINITEDIFFERENCES = flags[i,3]   )
      
      
      # RMSE evaluation on a fine grid
      xeval=runif(1000,0,1)  # quadrato con x tra 0 e 1
      yeval=runif(1000,0,1)  # quadrato con y tra 0 e 1
      teval=runif(1000,0,2)
      sol_eval=eval.FEM.time(output_CPP$fit.FEM.time,space.time.locations = cbind(teval,xeval,yeval))
      sol_exact=f(xeval,yeval,teval)
      RMSE<-function(f,g) sqrt(mean((f-g)^2))
      RMSE_matrix[i, k] = RMSE(sol_eval,sol_exact)
      
  }

}

quartz()
labels <- c("B-SPLINE", "PAR MONOLITICO", "PAR ITERATIVO", "FD MONOLITICO", "FD ITERATIVO")
boxplot(data.frame(t(RMSE_matrix)), las=1, col='gold', names = labels)

New_alg_mon <- c(0.1416270, 0.1386325, 0.1374639, 0.1439854)
New_alg_iter <- c(0.1416283, 0.1386323, 0.1374649, 0.1439876)

# Aggiunta della nuova riga alla matrice
RMSE_matrix <- rbind(RMSE_matrix, New_alg_mon)
RMSE_matrix <- rbind(RMSE_matrix, New_alg_iter)

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

flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
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
    
    #### Test 3.1: Without GCV
    results <- microbenchmark(output_CPP<-smooth.FEM.time(time_mesh = TimePoints, observations=observations, 
                                                          FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT,
                                                          PDE_parameters=PDE_parameters,
                                                          FLAG_PARABOLIC = flags[i,1], FLAG_ITERATIVE = flags[i,2], FLAG_FINITEDIFFERENCES = flags[i,3]   )
                              , times = 1)
    
    time_matrix[i, k] = results$time/10^6 #ms
  }
}
time_matrix

# Creazione della finestra grafica
quartz()
par(mfrow = c(1, 1))
colors <- seq(from = 1, to = dim(flags)[1]+2, by = 1)

# Ciclo for per il plottaggio dei dati
for (i in seq(from = 1, to = dim(flags)[1]+2, by = 1)) {
  if (i == 1) {
    plot(nRef*10+1, time_matrix[i,], type = "l", ylim = c(10, 500000), col = colors[i],log = "y", main = "Costi computazionali", xlab = "#nodes", ylab = "time [ms]", lwd = 2)
  } else {
    lines(nRef*10+1, time_matrix[i,], col = colors[i], log = "y", lwd = 2)
  }
}

legend_vector = c("B-SPLINE", "PARABOLICO MONOLITICO", "PARABOLICO ITERATIVO", "DIFFERENZE FINITE MONOLITICO", "DIFFERENZE FINITE ITERATIVO")
legend("topleft", legend = legend_vector , col = colors, lty = 1)

New_alg_mon <- c(56.59197, 519.39407, 1827.99189, 4779.248, 8617.9887)
New_alg_iter <- c(8.13768, 27.56307, 60.50222, 114.152, 201.5802)

# Aggiunta della nuova riga alla matrice
time_matrix <- rbind(time_matrix, New_alg_mon)
time_matrix <- rbind(time_matrix, New_alg_iter)

legend_vector = c("B-SP", "PAR MON", "PAR ITER", "FD MON", "FD ITER", "FD MON PREC", "FD ITER PREC")
legend("topleft", legend = legend_vector , col = colors, lty = 1)

################################################################################################
################################################################################################


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
beta_matrix <- matrix(0, nrow = dim(flags)[1], ncol = length(seed_numbers))

for (i in seq(from = 1,to = dim(flags)[1], by = 1)) { 
  k = 0
  for (seed_number in seed_numbers) {
    k = k + 1
    NumTimeInstants=5
    TimePoints=seq(0,pi,length.out =NumTimeInstants)
    
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
    lambdaT = 10^-2
    
    output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                                observations=observations, 
                                covariates = cov1,
                                FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT)
    
    beta_matrix[i, k] = output_CPP$beta
    
  }
  
}

quartz()
labels <- c("B-SPLINE", "PAR MONOLITICO", "PAR ITERATIVO", "FD MONOLITICO", "FD ITERATIVO")
boxplot(data.frame(t(beta_matrix)), las=1, col='gold', names = labels)



#### Test 4: semi circular domain ####
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

flags = rbind(c(FALSE, FALSE, FALSE),c(TRUE, FALSE, FALSE),c(TRUE, TRUE, FALSE),c(FALSE, FALSE, TRUE),c(FALSE, TRUE, TRUE))
nRef = seq(from = 1,to = 3, by = 1)

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
    
    #### Test 3.1: Without GCV
    results <- microbenchmark(output_CPP<-smooth.FEM.time(observations=observations,
                                                          incidence_matrix = incidence_matrix,
                                                          FEMbasis = FEMbasis, time_mesh = time_mesh,
                                                          lambdaS = lambdaS, lambdaT = lambdaT,
                                                          BC = BC,
                                                          PDE_parameters = PDE_parameters,
                                                          FLAG_PARABOLIC = flags[i,1],
                                                          FLAG_ITERATIVE = flags[i,2],
                                                          FLAG_FINITEDIFFERENCES = flags[i,3],
                                                          DOF.evaluation = NULL), times = 1)
    
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
    plot(nRef, time_matrix[i,], type = "l", ylim = c(10, 500000), col = colors[i],log = "y", main = "Costi computazionali", xlab = "#nodes", ylab = "time [ms]", lwd = 2)
  } else {
    lines(nRef, time_matrix[i,], col = colors[i], log = "y", lwd = 2)
  }
}

legend_vector = c("B-SPLINE", "PARABOLICO MONOLITICO", "PARABOLICO ITERATIVO", "DIFFERENZE FINITE MONOLITICO", "DIFFERENZE FINITE ITERATIVO")
legend("topleft", legend = legend_vector , col = colors, lty = 1)

################################################################################################
################################################################################################
