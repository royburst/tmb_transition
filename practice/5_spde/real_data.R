#### Try to run TMB on real u5m data (one age bin region)

############### SETUP
rm(list=ls())
gc()
options(scipen=999)
.libPaths('/home/j/temp/geospatial/packages')

library(INLA)
library(TMB)
library(data.table)
library(RandomFields)
library(raster)

dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
setwd(paste0(dir,"/practice/5_spde"))
system(paste0('cd ',dir,'\ngit pull origin develop'))
source('utils.R')

####################################################
## pull in data
## LETS DO WSSA AGE BIN 2
load('/share/geospatial/mbg/u5m/died/model_image_history/2017_04_04_22_07_02_bin2_wssa_0NA.RData')

# sort out the data a bit
modnames      <-  c('gam','gbm','ridge','enet','lasso')
child_fit_on  <-  '_cv_pred'
df            <-  rbind(df[,paste0(modnames) := lapply(modnames, function(x) get(paste0(x,child_fit_on)))])

#### Make objects needed for TMB to run
coords   <- cbind(df$longitude,df$latitude)

A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = coords,
                           group = df$period_id)

# make model frame, include rates as a covariate
X_xp = as.matrix(cbind(1, df[,c(modnames,'rates'),with=FALSE]))

# make SPDE matrices
spde <- inla.spde2.matern( mesh_s,alpha=2 ) ## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)

Data = list(n_i=nrow(df),                   ## Total number of observations
            n_x=mesh_s$n,                   ## Number of vertices in SPDE mesh
            n_t=length(unique(df$period_id)),                    ## Number of periods
            n_p=ncol(X_xp),                 ## Number of columns in covariate matrix X
  ##            x_s=mesh_s$idx$loc-1,           ## Association of each cluster with a given vertex in SPDE mesh
            c_i=df$died,                  ## Number of observed deaths in the cluster (N+ in binomial likelihood)
            Exp_i=df$N,             ## Number of observed exposures in the cluster (N in binomial likelihood)
            s_i=0:(nrow(X_xp)-1),                    ## no site specific effect in my model (different sampling locations in time)
            t_i=df$period_id-1,                ## Sample period ( starting at zero )
            X_xp=X_xp,                      ## Covariate design matrix
            G0=spde$param.inla$M0,          ## SPDE sparse matrix
            G1=spde$param.inla$M1,          ## SPDE sparse matrix
            G2=spde$param.inla$M2,          ## SPDE sparse matrix
            Aproj = A.proj,                 ## mesh to prediction point projection matrix
            options = 1)                 ## option1==1 use priors
              #spde=(spde$param.inla)[c('M1','M2','M3')])

Parameters = list(alpha   =  rep(0,ncol(X_xp)),                     ## FE parameters alphas
            log_tau_E=1.0,                                    ## log inverse of tau  (Epsilon)
            log_kappa=0.0,	                                  ## Matern Range parameter
            rho=0.5,
            epsilon=matrix(1,ncol=length(unique(df$period_id)),nrow=mesh_s$n))     ## GP locations

#### Try to fit using TMB
message('compiling')
templ <- "basic_spde" # _aoz" #spde2
TMB::compile(paste0(templ,".cpp"))
dyn.load( dynlib(templ) )

lower =  c(rep(-20,dim(X_xp)[2]),rep(-10,2),-0.999)
upper =  c(rep( 20,dim(X_xp)[2]),rep( 10,2), 0.999)

ptm <- proc.time()[3] # start timing

obj <- MakeADFun(data=Data, parameters=Parameters, random="epsilon", hessian=TRUE, DLL=templ)

## Run optimizer
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr,
                              lower       =    lower,
                              upper       =    upper,
                              control     =    list(eval.max=1e4, iter.max=1e4, trace=1))) # rel.tol=.01,step.max=10)))

fit_time <- proc.time()[3] - ptm


##### Make predictions
