#### Try to run TMB on real u5m data (one age bin region)

############### SETUP
rm(list=ls())
gc()
options(scipen=999)
.libPaths('/home/j/temp/geospatial/packages')

library(matrixcalc)
library(MASS)
library(INLA)
library(TMB)
library(data.table)
library(RandomFields)
library(raster)


# MBG FUNCTIONS
setwd('/share/code/geospatial/upstream/mbg/')
for(funk in list.files(recursive=TRUE,pattern='functions')) source(funk)

# TMB FUNCTIONS
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
setwd(paste0(dir,"/practice/5_spde"))
system(paste0('cd ',dir,'\ngit pull origin develop'))
source('utils.R')

# draws for prediction
ndraws <- 10

# make a chunky mesh or use the original?
use_orig_mesh <- TRUE

####################################################
## pull in data
## LETS DO WSSA AGE BIN 2
load('/share/geospatial/mbg/u5m/died/model_image_history/2017_04_04_22_07_02_bin2_wssa_0NA.RData')

# sort out the data a bit
modnames      <-  c('gam','gbm','ridge','enet','lasso')
child_fit_on  <-  '_cv_pred'
df            <-  rbind(df[,paste0(modnames) := lapply(modnames, function(x) get(paste0(x,child_fit_on)))])
df <- df[,c('N','died',modnames,'rates','period_id','latitude','longitude','weight'),with=FALSE]

# na omit stuff
nrow(df)
df <- na.omit(df)
nrow(df)

#### Make objects needed for TMB to run
coords   <- cbind(df$longitude,df$latitude)

# make model frame, include rates as a covariate
X_xp = as.matrix(cbind(1, df[,c(modnames,'rates'),with=FALSE]))

# try with a different mesh (too dense before?)
if(!use_orig_mesh){
  mesh_orig <- list(s=mesh_s,t=mesh_t)
  simple_polygon_list <-  load_simple_polygon(gaul_list        = get_gaul_codes('wssa'),
                                                buffer            = 0.4)
  simple_polygon   <- simple_polygon_list[[2]]
  mesh_s <- build_space_mesh(d = df,
                             simple = simple_polygon,
                             max_edge = 2,
                             mesh_offset = 2)
}

A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = coords,
                           group = df$period_id)


# n periods
nperiods <- length(unique(df$period_id))

# make SPDE matrices
spde <- inla.spde2.matern( mesh_s,alpha=2 ) ## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)

Data = list(n_i=nrow(df),                   ## Total number of observations
            n_x=mesh_s$n,                   ## Number of vertices in SPDE mesh
            n_t=nperiods,                    ## Number of periods
            n_p=ncol(X_xp),                 ## Number of columns in covariate matrix X
  ##            x_s=mesh_s$idx$loc-1,           ## Association of each cluster with a given vertex in SPDE mesh
            c_i=df$died,                  ## Number of observed deaths in the cluster (N+ in binomial likelihood)
            Exp_i=df$N,             ## Number of observed exposures in the cluster (N in binomial likelihood)
            s_i=0:(nrow(X_xp)-1),                    ## no site specific effect in my model (different sampling locations in time)
            t_i=df$period_id-1,                ## Sample period ( starting at zero )
            w_i=df$weight,
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


#########################################
#########################################
####  fit using TMB
message('compiling')
templ <- "basic_spde_weighted" # _aoz" #spde2
TMB::compile(paste0(templ,".cpp"))
dyn.load( dynlib(templ) )

lower =  c(rep(-20,dim(X_xp)[2]),rep(-10,2),-0.999)
upper =  c(rep( 20,dim(X_xp)[2]),rep( 10,2), 0.999)

######
ptm <- proc.time()[3] # start timing
#openmp(1)
obj <- MakeADFun(data=Data, parameters=Parameters, random="epsilon", hessian=TRUE, DLL=templ)

## Run optimizer
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr,
                              lower       =    lower,
                              upper       =    upper,
                              control     =    list(eval.max=1e4, iter.max=1e4, trace=1)))

fit_time <- proc.time()[3] - ptm
#########################################



#########################################
#########################################
##### Make predictions

ptm <- proc.time()[3]
SD0 = sdreport(obj,getReportCovariance=TRUE) #,bias.correct=TRUE)
tmb_sdreport_time <- proc.time()[3] - ptm

mu    <- c(SD0$value)
sigma <- SD0$cov

### simulate draws
npar   <- length(mu)

## make sigma symmetric
i <- 0
while(!is.symmetric.matrix(sigma)){
  sigma <- round(sigma, 10 - i)
  i <- i + 1
}
message(sprintf("rounded sigma to %i decimals to make it symmetric", 10 - i - 1))

# check we have a legit sigma matrix
is.positive.definite(sigma)

# predict out draws
draws  <- t(mvrnorm(n=ndraws,mu=mu,Sigma=sigma))

## separate out the draws
epsilon_draws <- draws[rownames(draws)=='Epsilon_xt',]
alpha_draws   <- draws[rownames(draws)=='alpha',]

# use one of the covariate layers as a 'fullsamplespace',
f_orig <- data.table(cbind(coordinates(cov_list[[1]][[1]]),t=1))
# add time periods
fullsamplespace <- copy(f_orig)
for(p in 2:nperiods){
  tmp <- f_orig
  tmp[,t := p]
  fullsamplespace <- rbind(fullsamplespace,tmp)
}

# pull out covariates in format we expect them
# a list of length periods with a brick of named covariates inside
new_cl <- list()
for(p in 1:nperiods){
  new_cl[[p]] <- list()
  for(n in names(cov_list)){
    new_cl[[p]][[n]] <- cov_list[[n]][[p]]
  }
  new_cl[[p]] <- brick(new_cl[[p]])
}

## get surface to project on to
pcoords <- cbind(x=fullsamplespace$x, y=fullsamplespace$y)
groups_periods <- fullsamplespace$t

## use inla helper functions to project the spatial effect.
A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = pcoords,
    group = groups_periods)

## values of S at each cell (long by nperiods)
cell_s <- as.matrix(A.pred %*% epsilon_draws)

## extract cell values  from covariates, deal with timevarying covariates here
vals <- list()
for(p in 1:nperiods){
  vals[[p]] <- extract(new_cl[[p]], pcoords[1:(nrow(fullsamplespace)/nperiods),])
  vals[[p]] <- (cbind(int = 1, vals[[p]]))
  vals[[p]] <- vals[[p]] %*% alpha_draws # same as getting cell_ll for each time period
}
cell_l <- do.call(rbind,vals)

## add together linear and st components
pred <- cell_l + cell_s

# save prediction timing
totalpredict_time <- proc.time()[3] - ptm

# plot it
summ  <- cbind(median=(apply(pred,1,median))) #,sd=(apply(pred,1,sd)))

# make a median raster
ras   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ[,1]), t=rep(1:nperiods,each=nrow(pred)/nperiods)),"p","t")

pdf('test2.pdf')
plot(ras)
dev.off()
#########################################


## Compare with INLA Predictions we have for this area as well
