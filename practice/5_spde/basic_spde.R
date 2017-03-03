
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
                    betas      = c(-3,-1,1,1) ,  #  Intercept coef and covariate coef For Linear predictors
                    scale      = 2            ,  #  Matern scale eparameter
                    Sigma2     = (1) ^ 2      ,  #  Variance (Nugget)
                    rho        = 0.9          ,  #  AR1 term
                    l          = 51           ,  #  Matrix Length
                    n_clusters = 100          ,  #  number of clusters sampled ]
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
              upper       =    c(rep(20 ,sum(names(obj$par)=='alpha')),rep(10,2),0.999),
              control     =    list(eval.max=1e4, iter.max=1e4, trace=0))
model.runtime = (Sys.time() - start_time)
opt0[["final_gradient"]] = obj$gr( opt0$par )


# Get standard errors
message('getting standard errors')
Report0 = obj$report()
SD0 = try( sdreport(obj) )

#### Prediction
message('making predictions')
# get back a surface (use epsilon, then interpolate somehow? )
epsilon=SD0$par.random[names(SD0$par.random)=="epsilon"] # should be mesh_s$n*nperiods long (make sure its in the right order)
# later  figure out how to get draws of this, for now, must be mean.. see SD0$diag.cov.random (mesh_s$n*5 long.. )

# get surface to project on to
pcoords = cbind(x=simobj$fullsamplespace$x, y=simobj$fullsamplespace$y)
groups_periods <- simobj$fullsamplespace$t

# use inla helper functions to project the spatial effect.
A.pred <- inla.spde.make.A(
  mesh = mesh_s,
  loc = pcoords,
  group = groups_periods)

# values of S at each cell (long by nperiods)
cell_s <- as.matrix(A.pred %*% epsilon)

# fixed effects
npars <- sum(names(opt0$par)=='alpha')

# extract cell values  from covariates
vals <- extract(simobj$cov.raster, pcoords[1:(nrow(simobj$fullsamplespace)/nperiod),])
vals <- (cbind(int = 1, vals))

cell_l <- vals %*% opt0$par[1:npars]
cell_l = rep(cell_l,nperiod) # since there are no time varying components

# add together linear and st components
pred <- cell_l + cell_s
pred <- plogis(as.vector(pred))

# make them into time bins
len = length(pred)/nperiod
res = data.table(pcoords,
                 pred,
                 t=rep(1:nperiod,each=len))
mean_ras=rasterFromXYZT(res,"pred","t")

# make a latent field raster
res2 = data.table(pcoords,
                  epsilon=cell_s,
                  t=rep(1:nperiod,each=len))
names(res2)=c('x','y','epsilon','t')
epsilon_ras=rasterFromXYZT(res2,"epsilon","t")

# plot
pdf('mean_raster_tmb.pdf',width=12,height=6)
par(mfrow=c(1,2))
plot(simobj$r.true.mr,main='TRUTH')
plot(mean_ras,        main='TMB FIT')
dev.off()
