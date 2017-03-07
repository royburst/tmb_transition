
rm(list=ls())
options(scipen=999)
.libPaths('/home/j/temp/geospatial/packages')

library(INLA)
library(TMB)
library(data.table)
library(RandomFields)
library(raster)

# need to set to same directory as the template file, also pull from git
# Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/5_spde"))

# source some functions made for this bit
source('utils.R')


###############################################################
# SIMULATE AND SET UP THE DATA
# Simulate a surface, this returns a list of useful objects like samples and truth
simobj <- mortsim(  nu         = 2            ,  #  Matern smoothness parameter (alpha in INLA speak)
                    betas      = c(-3,-1,1.5,1) ,  #  Intercept coef and covariate coef For Linear predictors
                    scale      = 1            ,  #  Matern scale eparameter
                    Sigma2     = (.5) ^ 2      ,  #  Variance (Nugget)
                    rho        = 0.9          ,  #  AR1 term
                    l          = 51           ,  #  Matrix Length
                    n_clusters = 10           ,  #  number of clusters sampled ]
                    n_periods  = 4            ,  #  number of periods (1 = no spacetime)
                    mean.exposure.months = 100,  #  mean exposure months per cluster
                    extent = c(0,1,0,1)       ,  #  xmin,xmax,ymin,ymax
                    ncovariates = 3           ,  #  how many covariates to include?
                    seed   = NULL             ,
                    returnall=TRUE            )


# get samples from which to fit
dt <- simobj[["d"]]

# set some stuff up
dt[,id:=1:.N]
coords   <- cbind(dt$x,dt$y)
nperiod  <- length(unique(dt$period))

# MESH For now use same mesh per time point
# TODO CHANGE THIS
mesh_s <- inla.mesh.2d(
  loc=coords,
  max.edge=c(0.2,0.2),
  cutoff=0.05)

nodes <- mesh_s$n # get number of mesh nodes
spde <- inla.spde2.matern( mesh_s ) # Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
# ^ this gives us a linear reduction of \Sigma^{-1} as:
  # \Sigma^{-1} = \kappa^4 M_0 + 2\kappa^2M_1 + M_2
  # M_2 = M_1M_0^{-1}M_1
  # Where the Ms are all sparse matrices stored as "dgTMatrix"
  names(spde$param.inla)


# pull covariate(s) at mesh knots
covs <- raster::extract(simobj$cov.raster,cbind(mesh_s$loc[,1],mesh_s$loc[,2]))
names(spde$param.inla)

# Data to pass to TMB
X_xp = cbind( 1, covs)

Data = list(n_i=nrow(dt),                   # Total number of observations
            n_x=mesh_s$n,                   # Number of vertices in SPDE mesh
            n_t=nperiod,                    # Number of periods
            n_p=ncol(X_xp),                 # Number of columns in covariate matrix X
            x_s=mesh_s$idx$loc-1,           # Association of each cluster with a given vertex in SPDE mesh
            c_i=dt$deaths,                  # Number of observed deaths in the cluster (N+ in binomial likelihood)
            Exp_i=dt$exposures,             # Number of observed exposures in the cluster (N in binomial likelihood)
            s_i=dt$id-1,                    # no site specific effect in my model (different sampling locations in time)
            t_i=dt$period-1,                # Sample period ( starting at zero )
            X_xp=X_xp,                      # Covariate design matrix
            G0=spde$param.inla$M0,          # SPDE sparse matrix
            G1=spde$param.inla$M1,          # SPDE sparse matrix
            G2=spde$param.inla$M2)          # SPDE sparse matrix

# staring values for parameters
Parameters = list(alpha   =  rep(0,ncol(X_xp)),                     # FE parameters alphas
                  log_tau_E=1.0,                                    # log inverse of tau  (Epsilon)
                  #                  log_tau_O=1.0,                 # log inverse of tau (SP)
                  log_kappa=0.0,	                                  # Matern Range parameter
                  rho=0.5,                                          # Autocorrelation term
                  epsilon=matrix(1,ncol=nperiod,nrow=mesh_s$n),     # GP
                  sp=matrix(rnorm(mesh_s$n)))                       # RE for mesh points


# which parameters are random
Random = c("epsilon",'sp')

##########################################################
## FIT MODEL
# Make object
# Compile
TMB::compile("basic_spde.cpp")
dyn.load( dynlib('basic_spde') )

obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=TRUE, DLL='basic_spde')
#obj$env$beSilent()

# Run optimizer
message('running optimizer')
start_time = Sys.time()
opt0 = nlminb(start       =    obj$par,
              objective   =    obj$fn,
              gradient    =    obj$gr,
              lower       =    c(rep(-20,sum(names(obj$par)=='alpha')),rep(-10,2),-0.999),
              upper       =    c(rep(20 ,sum(names(obj$par)=='alpha')),rep( 10,2), 0.999),
              control     =    list(eval.max=1e4, iter.max=1e4, trace=0))
model.runtime = (Sys.time() - start_time)
# opt0[["final_gradient"]] = obj$gr( opt0$par )
# head(summary(SD0))

# Get standard errors
message('getting standard errors')
#Report0 = obj$report()
SD0 = sdreport(obj,getReportCovariance=TRUE)
#fe_var_covar <- SD0$cov.fixed

#### Prediction
message('making predictions')
mu    <- c(SD0$par.fixed[names(SD0$par.fixed)=='alpha'],SD0$par.random[names(SD0$par.random)=="epsilon"])
sigma <- SD0$cov

## simulate draws
require(MASS)
npar   <- length(mu)
ndraws <- 50
draws  <- t(mvrnorm(n=ndraws,mu=mu,Sigma=sigma))
# ^ look for a quicker way to do this..cholesky

# separate out the draws
epsilon_draws <- draws[rownames(draws)=='epsilon',]
alpha_draws   <- draws[rownames(draws)=='alpha',]

# get surface to project on to
pcoords = cbind(x=simobj$fullsamplespace$x, y=simobj$fullsamplespace$y)
groups_periods <- simobj$fullsamplespace$t

# use inla helper functions to project the spatial effect.
A.pred <- inla.spde.make.A(
  mesh = mesh_s,
  loc = pcoords,
  group = groups_periods)


# values of S at each cell (long by nperiods)
cell_s <- as.matrix(A.pred %*% epsilon_draws)

# extract cell values  from covariates
vals <- extract(simobj$cov.raster, pcoords[1:(nrow(simobj$fullsamplespace)/nperiod),])
vals <- (cbind(int = 1, vals))

cell_l <- vals %*% alpha_draws
message('assuming no time varying covariates')
cell_l <- do.call(rbind,list(cell_l,cell_l,cell_l,cell_l)) # since there are no time varying components

# add together linear and st components
pred <- cell_l + cell_s
pred <- plogis(pred)

# make them into time bins
len = nrow(pred)/nperiod

# make summary plots - median, 2.5% and 97.5%
summ <- t(apply(pred,1,quantile,probs=c(.5,.025,.975)))

median_ras <- rasterFromXYZT(data.table(pcoords,p=summ[,1],t=rep(1:nperiod,each=len)),"p","t")
lower_ras  <- rasterFromXYZT(data.table(pcoords,p=summ[,2],t=rep(1:nperiod,each=len)),"p","t")
upper_ras  <- rasterFromXYZT(data.table(pcoords,p=summ[,3],t=rep(1:nperiod,each=len)),"p","t")

draw1 <- rasterFromXYZT(data.table(pcoords,p=pred[,1],t=rep(1:nperiod,each=len)),"p","t")
draw2 <- rasterFromXYZT(data.table(pcoords,p=pred[,2],t=rep(1:nperiod,each=len)),"p","t")
draw3 <- rasterFromXYZT(data.table(pcoords,p=pred[,3],t=rep(1:nperiod,each=len)),"p","t")
draw4 <- rasterFromXYZT(data.table(pcoords,p=pred[,4],t=rep(1:nperiod,each=len)),"p","t")

# plot
pdf('mean_raster_tmb.pdf',width=12,height=6)
p<-2
zmax <- max(c(as.vector(simobj$r.true.mr[[p]]),as.vector(upper_ras[[1]])))
par(mfrow=c(2,4))
plot(simobj$r.true.mr[[p]],main='TRUTH w/ sample locs',zlim=c(0,zmax))
points(simobj$d$x[simobj$d$period==p],simobj$d$y[simobj$d$period==p])
plot(median_ras[[p]],      main='TMB MEDIAN',zlim=c(0,zmax))
plot(lower_ras[[p]],       main='TMB LOWER',zlim=c(0,zmax))
plot(upper_ras[[p]],       main='TMB UPPER',zlim=c(0,zmax))
plot(draw1[[p]],       main='TMB DRAW',zlim=c(0,zmax))
plot(draw2[[p]],       main='TMB DRAW',zlim=c(0,zmax))
plot(draw3[[p]],       main='TMB DRAW',zlim=c(0,zmax))
plot(draw4[[p]],       main='TMB DRAW',zlim=c(0,zmax))

dev.off()
