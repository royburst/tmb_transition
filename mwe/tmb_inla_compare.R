## -------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------
## AUTHOR: ROY BURSTEIN
## DECEMBER 2017
## Minimal working example to compare geostatistical models in TMB and INLA
## Instructions:
##    Set the variables needed at the top, the rest of the script should run
##    and save a table comparing run times and parameter estimates from INLA and TMB
## This script does the following:
##    1. Input user settings
##    2. Load Packages
##    3. Load and setup data
##    4. Fit model in TMB
##    5. Make predictions from fit TMB model
##    6. Fit model in INLA
##    7. Make predictions from fit INLA model
##    8. Make and save in code directory a table and plot of comparisons
## -------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------
### 1. PLEASE SET THE FOLLOWING USER OPTIONS:
## -------------------------------------------------------------------------------------
## A. Relative paths for you machine
# Path the R packages, set to NULL if this doesnt apply to you
packloc <- '/home/j/temp/geospatial/geos_packages/'

# Path to code directory: where this script and model.cpp live
codedir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition/mwe")

# Path to data: where the shapefile, covariate rasters, and data file live.
#  Also where the outputs of this test will be saved. Can be same as codedir
datadir <- '/homes/royburst/tmb_mwe/'

## B. set some fitting options: scaling these is what affects run time the most
# Decimal Degrees of max edge for finite elements mesh
#  generally the smaller the better with 0.25 or smaller as max size for anything publishable
max_edge <- 4

# Number of draws to take for prediction. Test on 100, with 1000 as minimum for anything publishable
ndraws   <- 200

# Number of cores for parallelization. Make sure you've set it such that ndraws/ncores is an integer.
ncores   <- 4
## -------------------------------------------------------------------------------------



## -------------------------------------------------------------------------------------
### 2. QUICK ENVIRON SET UP STUFF
## -------------------------------------------------------------------------------------
# packages
if(!is.null(packloc)) .libPaths(packloc)
for(pack in c('INLA','TMB','data.table','parallel','raster'))
  library(pack, character.only=TRUE)

# set working directory
setwd(codedir)

# scipen
options(scipen=999)

# patch fix for inla mesher, if needed
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
## -------------------------------------------------------------------------------------





## -------------------------------------------------------------------------------------
### 3. LOAD DATA AND SET EVERYTHING UP FOR MODEL FITTING AND PREDICTION IN TMB AND INLA
## -------------------------------------------------------------------------------------
## Load the data
df <- fread(sprintf('%s/tmb_mwe_df.csv',datadir)) # main data for fitting
df[, died := round(died,1)] # clean it up for this test (no nonints in binomial k)
simple_polygon <- readRDS(sprintf('%s/shapefile.RDS',datadir)) # shapefile
cov_list <- readRDS(sprintf('%s/cov_rasters.RDS',datadir)) # list of raster covariates

# names of the fixed effects that will be used in the model
covs      <-  c('gam','gbm','lasso')

# matrix of coordinate locations for each row in the data
coords   <- cbind(df$longitude,df$latitude)

# model frame, with intercept and covariates
X_xp = as.matrix(cbind(1, df[,c(covs),with=FALSE]))

# construct the finite elements mesh
mesh_s <- inla.mesh.2d(
  boundary = inla.sp2segment(simple_polygon),
  loc      = cbind(df$longitude,df$latitude),
  max.edge = max_edge,
  offset   = 2,
  cutoff   = max_edge)
#  plot a picture of it
plot(mesh_s,asp=1);points(df$longitude,df$latitude,col=df$fyear)

# constuct a projection matrix from data to mesh nodes. used for prediction
A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = coords,
                           group = df$period_id)

# number of time periods in the data, in this MWE there are 4
nperiods <- length(unique(df$period_id))

# make SPDE matrices which we use in our likelihood functoin
spde <- inla.spde2.matern(mesh_s,alpha=2 )

# use one of the covariate layers as a 'fullsamplespace' dt, one row per pixel
#  basically used as the template for prediction
f_orig <- data.table(cbind(coordinates(cov_list[[1]][[1]]),t=1))
fullsamplespace <- copy(f_orig)
for(p in 2:nperiods){
  tmp <- f_orig
  tmp[,t := p]
  fullsamplespace <- rbind(fullsamplespace,tmp)
}

# pull out covariates in format we expect them:
#  a list of length periods with a raster-brick of named covariates inside
new_cl <- list()
for(p in 1:nperiods){
  new_cl[[p]] <- list()
  for(n in names(cov_list)){
    new_cl[[p]][[n]] <- cov_list[[n]][[p]]
  }
  new_cl[[p]] <- brick(new_cl[[p]])
}

# get space-time coordinates of the surface to project on to, in table form
pcoords <- cbind(x=fullsamplespace$x, y=fullsamplespace$y)
groups_periods <- fullsamplespace$t

# use inla helper functions to project the spatial effect, from mesh nodes to prediction template pixels
A.pred <- inla.spde.make.A(
    mesh  = mesh_s,
    loc   = pcoords,
    group = groups_periods)

## extract cell values from covariates into the predction template space
cov_vals <- list()
for(p in 1:nperiods){
  message(p)
  cov_vals[[p]] <- raster::extract(new_cl[[p]], pcoords[1:(nrow(fullsamplespace)/nperiods),])
  cov_vals[[p]] <- (cbind(int = 1, cov_vals[[p]]))
}

# Data list to be passed on to TMB, require in model.cpp
Data = list(num_i=nrow(df),               ## Total number of observations
            num_s=mesh_s$n,               ## Number of vertices in SPDE mesh
            num_t=nperiods,               ## Number of periods
            num_z=1,                      ## Z dimension, set to 1 for simplicity for now
            y_i=df$died,                  ## Number of observed deaths in the cluster (N+ in binomial likelihood)
            n_i=df$N,                     ## Number of observed exposures in the cluster (N in binomial likelihood)
            t_i=df$period_id-1,           ## Sample period ( starting at zero )
            w_i=df$weight,                ## Row wise data rate
            X_ij=X_xp,                    ## Covariate design matrix
            M0=spde$param.inla$M0,        ## SPDE sparse matrix
            M1=spde$param.inla$M1,        ## SPDE sparse matrix
            M2=spde$param.inla$M2,        ## SPDE sparse matrix
            Aproj = A.proj,               ## mesh to prediction point projection matrix
            options = c(1,1))             ## option1==1 use priors
                                          ## option2==1 turn ADREPORT off

# Set starting values for parameters
Parameters = list(alpha_j   = rep(0,ncol(X_xp)),  ## FE parameters alphas
                  logtau    = 1.0,                ## log inverse of tau, variance of GP
                  logkappa  = 0.0,	              ## Matern kappa
                  trho      = 0.5,                ## temporal AR1 rho
                  zrho      = 0.5,                ## Z AR1 rho, will be mapped out ignored for now since not using this dimension
                  Epsilon_stz=matrix(1, nrow=mesh_s$n, ncol=nperiods)) ## GP locations

## -------------------------------------------------------------------------------------



## -------------------------------------------------------------------------------------
### 4. Fit the model using TMB
## -------------------------------------------------------------------------------------
# load, and if necessary compile the CPP template
templ <- "model"
TMB::compile(paste0(templ,".cpp"))
dyn.load(dynlib(templ))

# set some options
openmp(ncores)
config(tape.parallel=0, DLL=templ)

# Run AD
obj <- MakeADFun(data=Data, parameters=Parameters, map=list(zrho=factor(NA)), random="Epsilon_stz", hessian=TRUE, DLL=templ)
runSymbolicAnalysis(obj) # METIS ordering

# Fit the model
ptm <- proc.time()[3]
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr))
fit_time_tmb <- proc.time()[3] - ptm

## -------------------------------------------------------------------------------------




## -------------------------------------------------------------------------------------
### 5. Make predictions with the fit TMB model
## -------------------------------------------------------------------------------------
# Get precision matrix using sdreport
ptm <- proc.time()[3]
SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE)
tmb_sdreport_time <- proc.time()[3] - ptm

# get the mean values for all fixed and random effects
mu    <- c(SD0$par.fixed,SD0$par.random)

# simulate mvn draws from the precision matrix
rmvnorm_prec <- function(mu, prec, n.sims) {
  z     <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv <- Cholesky(prec)
  mu + solve(as(L_inv, "pMatrix"), solve(t(as(L_inv, "Matrix")), z))
}
ptm2 <- proc.time()[3]
draws <- do.call('cbind',
          mclapply(rep(ndraws/ncores,ncores), rmvnorm_prec, mu = mu , prec = SD0$jointPrecision, mc.cores = ncores))
tmb_get_draws_time <- proc.time()[3] - ptm2

# separate out the draws by random and fixed effects
parnames      <- c(names(SD0$par.fixed), names(SD0$par.random))
epsilon_draws <- draws[parnames=='Epsilon_stz',]
alpha_draws   <- draws[parnames=='alpha_j',]

## Use the pred matrix to get RE values at each cell (long by nperiods)
cell_s <- as.matrix(A.pred %*% epsilon_draws)

# multiply covariate values by draws of coefficient values
tmb_vals <- list()
for(p in 1:nperiods)  tmb_vals[[p]] <- cov_vals[[p]] %*% alpha_draws

# rbind them to get fixed effects by cell (long by nperiods), and to match cell_s dimensions
cell_l <- do.call(rbind,tmb_vals)

# add together predictions from fixed and random effects
#  the resulting matrix should be length(fullsamplespace) by ndraws
pred_tmb <- cell_l + cell_s

# save prediction timing
totalpredict_time_tmb <- proc.time()[3] - ptm

## -------------------------------------------------------------------------------------




## -------------------------------------------------------------------------------------
### 6. FIT THE SAME MODEL USING R-INLA
## -------------------------------------------------------------------------------------
# same basic inla setup stuff first
space   <- inla.spde.make.index("space",
                                n.spde = spde$n.spde,
                                n.group = nperiods)
design_matrix <- data.frame(int = 1, df[, covs, with=F])
stack.obs <- inla.stack(tag     = 'est',
                        data    = list(died=df$died), ## response
                        A       = list(A.proj,1),
                        effects = list(space,  design_matrix))
formula <- formula(paste0('died ~ -1+int+',
(paste(covs, collapse = ' + ')),
' + f(space, model = spde, group = space.group, control.group = list(model = \'ar1\'))'))

# fit the model
ptm <- proc.time()[3]
inla.setOption("enable.inla.argument.weights", TRUE)
res_fit <- inla(formula,
                data              = inla.stack.data(stack.obs),
                control.predictor = list(A       = inla.stack.A(stack.obs),
                                         link    = 1,
                                         compute = FALSE),
                control.fixed   = list(expand.factor.strategy = 'inla'),
                control.inla    = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                control.compute = list(config = TRUE),
                family          = 'binomial',
                num.threads     = ncores,
                Ntrials         = df$N,
                weights         = df$weight,
                verbose         = TRUE,
                keep            = FALSE)
fit_time_inla <- proc.time()[3] - ptm
## -------------------------------------------------------------------------------------




## -------------------------------------------------------------------------------------
### 7. PREDICT FROM THE INLA FITTED MODEL
## -------------------------------------------------------------------------------------
# draw samples of parameter values
ptm <- proc.time()[3]
draws <- inla.posterior.sample(ndraws, res_fit)
inla_get_draws_time <- proc.time()[3] - ptm

# get parameter names
par_names <- rownames(draws[[1]]$latent)

# index to spatial field and fixed coefficient samples
s_idx <- grep('^space.*', par_names)
l_idx <- which(!c(1:length(par_names)) %in% grep('^space.*|Predictor', par_names))

# get samples as matrices
pred_s <- sapply(draws, function (x) x$latent[s_idx])
pred_l <- sapply(draws, function (x) x$latent[l_idx])
rownames(pred_l) <- res_fit$names.fixed

# get samples of random effects for all coo locations
s <- as.matrix(A.pred %*% pred_s)

# multiply covariate values by draws of coefficient values
inla_vals <- list()
for(p in 1:nperiods)  inla_vals[[p]] <- cov_vals[[p]] %*% pred_l

# rbind them to get fixed effects by cell (long by nperiods), and to match cell_s dimensions
l <- do.call(rbind,inla_vals)

# add together predictions from fixed and random effects
#  the resulting matrix should be length(fullsamplespace) by ndraws
pred_inla <- s+l

# save prediction timing
totalpredict_time_inla <- proc.time()[3] - ptm

## -------------------------------------------------------------------------------------






## -------------------------------------------------------------------------------------
### 8. MAKE AND SAVE COMPARISON METRICS AND PLOTS
## -------------------------------------------------------------------------------------
## A. Comparison Table
res <- data.table(st_mesh_nodes    =rep(nrow(epsilon_draws),2)) # number of mesh nodes
# some info on user options
res[,cores            := rep(ncores,2)]
res[,s_mesh_max_edge  := rep(max_edge,2)]
res[,periods          := c(4,4)]
res[,draws            := c(ndraws,ndraws)]

# time variables
res[,fit_time_seconds        := c(fit_time_inla,fit_time_tmb)]
res[,total_pred_time_seconds := c(totalpredict_time_inla,totalpredict_time_tmb)]
res[,get_draws_time_seconds  := c(inla_get_draws_time,tmb_sdreport_time+tmb_get_draws_time)]

# fe coefficients
for(i in 1:length(res_fit$names.fixed)){
  fn <- res_fit$names.fixed[i]
  res[[paste0('fixedeff_',fn,'_mean')]] <- c(res_fit$summary.fixed$mean[i], SD0$par.fixed[i])
  res[[paste0('fixedeff_',fn,'_sd')]]   <- c(res_fit$summary.fixed$sd[i], sqrt(SD0$cov.fixed[i,i]))
}

# hyperparameters
res[,hyperpar_logtau_mean   := c(res_fit$summary.hyperpar[1,1], SD0$par.fixed['logtau']) ]
res[,hyperpar_logtau_sd     := c(res_fit$summary.hyperpar[1,2], sqrt(SD0$cov.fixed['logtau','logtau'])) ]
res[,hyperpar_logkappa_mean := c(res_fit$summary.hyperpar[2,1], SD0$par.fixed['logkappa']) ]
res[,hyperpar_logkappa_sd   := c(res_fit$summary.hyperpar[2,2], sqrt(SD0$cov.fixed['logkappa','logkappa'])) ]
res[,hyperpar_rho_mean      := c(res_fit$summary.hyperpar[3,1], SD0$par.fixed['trho']) ]
res[,hyperpar_rho_sd        := c(res_fit$summary.hyperpar[3,2], sqrt(SD0$cov.fixed['trho','trho'])) ]

# finish it up
rr <- data.table(item=colnames(res))
rr <- cbind(rr,t(res))
names(rr) <- c('_','R-INLA','TMB')
rr$diff   <- rr[,2]-rr[,3]

# save it in the DATADIR directory
write.csv(rr,sprintf('%s/tmb_inla_comparison_metrics.csv',datadir))


## -----------------------------------------
## B. Plot predictions
## Helper function for turning an xyzt table into a raster
rasterFromXYZT <- function(table,z,t){
  n_periods = length(unique(table[,t]))
  table$t=table[,t]
  res=  stack(rasterFromXYZ(as.matrix(table[t==1,c('x','y',z),with=F])))
  if(n_periods>1)
    for(r in 2:n_periods)
      res=addLayer(res, rasterFromXYZ(as.matrix(table[t==r,c('x','y',z),with=F])))
  return(res)
}

# make some useful summaries of predictive draws first
summ_inla <- cbind(median=(apply(pred_inla,1,median)),sd=(apply(pred_inla,1,sd)))
summ_tmb  <- cbind(median=(apply(pred_tmb,1,median)) ,sd=(apply(pred_tmb,1,sd)))

# make rasters
ras_med_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,1]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
ras_med_tmb   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,1]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")
ras_sdv_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,2]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
ras_sdv_tmb   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,2]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")

# plot and save the rasters
cellIdx <- function (x) which(!is.na(getValues(x[[1]])))
library(grid); library(gridExtra); library(viridis)
pdf(sprintf('%s/tmb_inla_comparison_plots.pdf',datadir), height=10,width=14)
for(thing in c('median','stdev')){
    layout(matrix(1:20, 4, 5, byrow = TRUE))
    samp=sample(cellIdx(ras_med_inla[[1]]),1e4)
    for(i in 1:4){
      if(thing=='median'){
        rinla <-  ras_med_inla[[i]]
        rtmb  <-  ras_med_tmb[[i]]
      }
      if(thing=='stdev'){
        rinla <-  ras_sdv_inla[[i]]
        rtmb  <-  ras_sdv_tmb[[i]]
      }

      tmp <- subset(df,period_id==i)
      tmp$dat <- tmp$died/tmp$N

      # rasters
      par(mar = c(0, 0, 1.4, 1),bty='n')
      maxes <- max(c(as.vector(rtmb),as.vector(rinla)),na.rm=TRUE)
      plot(rinla-rtmb, maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend=T, main=paste0('DIFFERENCE ',thing))
      plot(rtmb,  maxpixel=1e7, col=rainbow(100), axes=FALSE, legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main='TMB',zlim=c(0,maxes))
      plot(rinla, maxpixel=1e7, col=rainbow(100), axes=FALSE, legend=FALSE, main='R-INLA',zlim=c(0,maxes))

      # scatter
      par(mar = c(4, 4, 2, 2),bty='n')
      plot(x=as.vector(rinla)[samp],y=as.vector(rtmb)[samp],xlab='R-INLA',ylab='TMB',cex=.01,pch=19,main='COMPARE')
      lines(x=c(0,maxes),y=c(0,maxes),col='red')

      # plot data points over the shapefile
      plot(simple_polygon, main='DATA LOCATIONS')
      points(x=tmp$longitude,y=tmp$latitude, pch=19, cex=tmp$dat)
    }
}
dev.off()
## -------------------------------------------------------------------------------------
