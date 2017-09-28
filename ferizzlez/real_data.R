#### Try to run TMB on real u5m data (one age bin region)

############### SETUP
rm(list=ls())
gc()
options(scipen=999)
.libPaths('/home/j/temp/geospatial/geos_packages')

library(INLA)
require('TMB', lib.loc='/snfs2/HOME/azimmer/R/x86_64-pc-linux-gnu-library/3.3/') # COMPILED WITH METIS
library(data.table)
library(RandomFields)
library(raster)

## need to set to same directory as the template file, also pull from git
## Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))

# MBG FUNCTIONS
setwd('/share/code/geospatial/royburst/mbg/')
for(funk in list.files(recursive=TRUE,pattern='functions')) {
  message(funk)
  try(source(funk))
}

# TMB FUNCTIONS
setwd(paste0(dir,"/ferizzlez"))
source('../utils.R')

# draws for prediction
ndraws <- 10

# make a chunky mesh or use the original?
use_orig_mesh <- FALSE # TRUE

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
                             max_edge = 1,
                             mesh_offset = 2)
}

A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = coords,
                           group = df$period_id)


# n periods
nperiods <- length(unique(df$period_id))

# make SPDE matrices
spde <- inla.spde2.matern( mesh_s,alpha=2 ) ## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)


Data = list(num_i=nrow(df),                 ## Total number of observations
            num_s=mesh_s$n,                 ## Number of vertices in SPDE mesh
            num_t=nperiods,                  ## Number of periods
            num_z=1,
            y_i=df$died,                  ## Number of observed deaths in the cluster (N+ in binomial likelihood)
            n_i=df$N,               ## Number of observed exposures in the cluster (N in binomial likelihood)
            t_i=df$period_id-1,                ## Sample period ( starting at zero )
            w_i=df$weight,
            X_ij=X_xp,                      ## Covariate design matrix
            M0=spde$param.inla$M0,          ## SPDE sparse matrix
            M1=spde$param.inla$M1,          ## SPDE sparse matrix
            M2=spde$param.inla$M2,          ## SPDE sparse matrix
            Aproj = A.proj,                 ## mesh to prediction point projection matrix
            options = c(1))                 ## option1==1 use priors


  ## staring values for parameters
  Parameters = list(alpha_j   =  rep(0,ncol(X_xp)),                 ## FE parameters alphas
                    logtau=1.0,                                     ## log inverse of tau  (Epsilon)
                    logkappa=0.0,	                                  ## Matern Range parameter
                    trho=0.5,
                    zrho=0.5,
                    Epsilon_stz=matrix(1, nrow=mesh_s$n, ncol=nperiods))     ## GP locations


#########################################
#########################################
####  fit using TMB
system(paste0('cd ',dir,'\ngit pull origin develop'))

templ <- "model"
TMB::compile(paste0(templ,".cpp"))
dyn.load( dynlib(templ) )
openmp(10)
config(tape.parallel=0, DLL=templ)

obj <- MakeADFun(data=Data, parameters=Parameters,  map=list(zrho=factor(NA)), random="Epsilon_stz", hessian=TRUE, DLL=templ)


## Run optimizer
ptm <- proc.time()[3]
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr))
                        #      control     =    list(eval.max=1e4, iter.max=1e4, trace=1)))
tmb_fit_time <- proc.time()[3] - ptm

#########################################



#########################################
#########################################
##### Make predictions

# Get standard errors
## Report0 = obj$report()
ptm <- proc.time()[3]
SD0 = sdreport(obj,getReportCovariance=TRUE,bias.correct=TRUE)
## fe_var_covar <- SD0$cov.fixed
tmb_sdreport_time <- proc.time()[3] - ptm

##### Prediction
message('making predictions')
mu    <- c(SD0$value)
sigma <- SD0$cov

### simulate draws
require(MASS)
npar   <- length(mu)

## make sigma symmetric
require(matrixcalc)
i <- 0
while(!is.symmetric.matrix(sigma)){
  sigma <- round(sigma, 10 - i)
  i <- i + 1
}
if(i>9) stop('Too much rounding of sigma to make it symmetric, something is wrong.')

message(sprintf("rounded sigma to %i decimals to make it symmetric", 10 - i - 1))


## round more or add to diagonal to make sigma pos-def
if(!is.positive.definite(sigma)){
  if(fixsigma){
    message('Sigma was not positive definite, tryna fix it')
    i <- 0
    sigma2 <- sigma
    while(!is.positive.definite(sigma2) & i < 6){
      message(i)
      sigma2 <- round(sigma2, 10 - i)
      i <- i + 1
    }
    if(is.positive.definite(sigma2)){
      sigma <- sigma2
      message(sprintf("rounded sigma to %i decimals to make it pos-def", 10 - i - 1))
    }else{
      i <- 0
      message('adding to diag approach')
      while(!is.positive.definite(sigma)){
        message(i)
        sigma <- sigma + diag(1, nrow(sigma))
        i <- i + 1
      }
      message(sprintf("added %i to the diagonal to make sigma pos-def", i))
    }
  } else {
    stop('Sigma was not definite positive and you set fixsigma==FALSE. So this is an error.')
  }
} else{
  message('Sigma was positive definite hurray')
}

## now we can take draws
message('Predicting Draws')
draws  <- t(mvrnorm(n=ndraws,mu=mu,Sigma=sigma))
  ## ^ look for a quicker way to do this..cholesky

## separate out the draws
epsilon_draws <- draws[rownames(draws)=='Epsilon_stz',]
alpha_draws   <- draws[rownames(draws)=='alpha_j',]


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
  message(p)
  vals[[p]] <- raster::extract(new_cl[[p]], pcoords[1:(nrow(fullsamplespace)/nperiods),])
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

pdf('test3.pdf')
plot(ras)
dev.off()
#########################################


## Compare with INLA Predictions we have for this area as well
