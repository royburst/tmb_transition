#### Try to run TMB on real u5m data (one age bin region)

# NOTES: RUN THESE LINES BEFORE LAUNCHING R
# source /homes/imdavis/intel/mkl/bin/mklvars.sh intel64
# export MKL_INTERFACE_LAYER=GNU,LP64
# export MKL_THREADING_LAYER=GNU
# export OMP_NUM_THREADS=25
# numactl --physcpubind=+0-9
# /homes/imdavis/R_mkl_geos/R-3.4.1-mkl_gcc484/R-3.4.1/bin/R
# OR JUST IN my basrc i have R_MKL $numcores

# OR JUST RUN THE SCRIPT:
 ### /homes/imdavis/R_mkl_geos/R-3.4.1-mkl_gcc484/R-3.4.1/bin/R < /homes/royburst/tmb_transition/ferizzlez/real_data_tmb_inla_compare.R --no-save --args .25 100 20


############### SETUP
rm(list=ls())
gc()
options(scipen=999)
.libPaths('/home/j/temp/geospatial/geos_packages')
##.libPaths('/home/j/temp/geospatial/packages')


library(INLA)
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
require('TMB', lib.loc='/snfs2/HOME/azimmer/R/x86_64-pc-linux-gnu-library/3.3/') # COMPILED WITH METIS
#require('mvnfast',lib.loc='/home/j/temp/geospatial/sandbox_packages/')
library(data.table)
#library(RandomFields)
library(raster)
library(viridis)
library(ggplot2)
library(MASS)
library(parallel)


## Print out information on session
sessionInfo()

# MBG FUNCTIONS
#setwd('/share/code/geospatial/royburst/mbg/')
for(funk in list.files(recursive=TRUE,pattern='functions')) {
  message(funk)
  try(source(funk))
}

# Set WD
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
setwd(paste0(dir,"/mwe"))
source('../utils.R')

# draws for prediction

# set some fitting options
max_edge <- 4     # mesh size
ndraws   <- 200   # number of draws for prediction
ncores   <- 4     # cores to fit on

####################################################
####################################################
### SET UP NECESSARY BEFORE RUNNING MODEL

## pull in data
## load in the data (TODO SAVE THIS IN REPO AS WELL)
df <- fread('/homes/royburst/tmb_mwe/tmb_mwe_df.csv') # main data for fitting
simple_polygon <- readRDS('/homes/royburst/tmb_mwe/shapefile.RDS') # shapefile
cov_list <- readRDS('/homes/royburst/tmb_mwe/cov_rasters.RDS') # shapefile
#load('/share/geospatial/mbg/u5m/died/model_image_history/2017_04_04_22_07_02_bin2_wssa_0NA.RData')

# fixed effects for the model
covs      <-  c('gam','gbm','lasso')

#### Make objects needed for TMB
coords   <- cbind(df$longitude,df$latitude)

# make model frame, include rates as a covariate
X_xp = as.matrix(cbind(1, df[,c(covs),with=FALSE]))

# construct the mesh and plot a picture of it
mesh_s <- inla.mesh.2d(
  boundary = inla.sp2segment(simple_polygon),
  loc      = cbind(df$longitude,df$latitude),
  max.edge = max_edge,
  offset   = 2,
  cutoff   = max_edge)
plot(mesh_s,asp=1);points(df$longitude,df$latitude,col=df$fyear)

# constuct a projection matrix from data to mesh nodes. used for prediction
A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = coords,
                           group = df$period_id)

# n periods
nperiods <- length(unique(df$period_id))

# make SPDE matrices which we use in our likelihood functoin
spde <- inla.spde2.matern(mesh_s,alpha=2 ) ## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)


# use one of the covariate layers as a 'fullsamplespace', basically used as the template for prediction
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

## use inla helper functions to project the spatial effect, from mesh nodes to prediction template
A.pred <- inla.spde.make.A(
    mesh  = mesh_s,
    loc   = pcoords,
    group = groups_periods)

## extract cell values  from covariates
cov_vals <- list()
for(p in 1:nperiods){
  message(p)
  cov_vals[[p]] <- raster::extract(new_cl[[p]], pcoords[1:(nrow(fullsamplespace)/nperiods),])
  cov_vals[[p]] <- (cbind(int = 1, cov_vals[[p]]))
}

# Data list passed on to TMB
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
                                          ## option2==1 ADREPORT off

## starting values for parameters
Parameters = list(alpha_j   = rep(0,ncol(X_xp)),  ## FE parameters alphas
                  logtau    = 1.0,                ## log inverse of tau, variance of GP
                  logkappa  = 0.0,	              ## Matern kappa
                  trho      = 0.5,                ## temporal AR1 rho
                  zrho      = 0.5,                ## ag eAR1 rho, ignored for now
                  Epsilon_stz=matrix(1, nrow=mesh_s$n, ncol=nperiods)) ## GP locations
#########################################



#########################################
#########################################
####  fit using TMB
templ <- "model"
TMB::compile(paste0(templ,".cpp"))
dyn.load( dynlib(templ) )
openmp(ncores)
config(tape.parallel=0, DLL=templ)

# Run AD
obj <- MakeADFun(data=Data, parameters=Parameters, map=list(zrho=factor(NA)), random="Epsilon_stz", hessian=TRUE, DLL=templ)

## Run optimizer
ptm <- proc.time()[3]
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr))
                        #      control     =    list(eval.max=1e4, iter.max=1e4, trace=1)))
fit_time_tmb <- proc.time()[3] - ptm

#########################################




#########################################
#########################################
##### Make predictions with TMB model output

# Get precision matrix
ptm <- proc.time()[3]
SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE)
tmb_sdreport_time <- proc.time()[3] - ptm

##### Prediction
message('making predictions')
mu    <- c(SD0$par.fixed,SD0$par.random)

### simulate draws from the precision matrix
rmvnorm_prec <- function(mu, prec, n.sims) {
  z     <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv <- Cholesky(prec)
  mu + solve(as(L_inv, "pMatrix"), solve(t(as(L_inv, "Matrix")), z))
}
ptm2 <- proc.time()[3]
draws <- do.call('cbind',
          mclapply(rep(ndraws/ncores,ncores), rmvnorm_prec, mu = mu , prec = SD0$jointPrecision, mc.cores = ncores))
tmb_get_draws_time <- proc.time()[3] - ptm2

## separate out the draws by random and fixed effects
parnames      <- c(names(SD0$par.fixed), names(SD0$par.random))
epsilon_draws <- draws[parnames=='Epsilon_stz',]
alpha_draws   <- draws[parnames=='alpha_j',]

## Use the pred matrix to get RE values at each cell ( long by nperiods)
cell_s <- as.matrix(A.pred %*% epsilon_draws)

# covariate values by alpha draws
tmb_vals <- list()
for(p in 1:nperiods)  tmb_vals[[p]] <- cov_vals[[p]] %*% alpha_draws

# fixed effects by cell
cell_l <- do.call(rbind,tmb_vals)

## add together FE and GP components
pred_tmb <- cell_l + cell_s

# save prediction timing
totalpredict_time_tmb <- proc.time()[3] - ptm




#########################################
#########################################
#### fit using INLA

# same basic inla setup stuff
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

# set fitting time
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




#########################################
#########################################
#### predict with INLA
ptm <- proc.time()[3]
draws <- inla.posterior.sample(ndraws, res_fit)
inla_get_draws_time <- proc.time()[3] - ptm

## get parameter names
par_names <- rownames(draws[[1]]$latent)

## index to spatial field and fixed coefficient samples
s_idx <- grep('^space.*', par_names)
l_idx <- which(!c(1:length(par_names)) %in% grep('^space.*|Predictor', par_names))

## get samples as matrices
pred_s <- sapply(draws, function (x) x$latent[s_idx])
pred_l <- sapply(draws, function (x) x$latent[l_idx])
rownames(pred_l) <- res_fit$names.fixed

## get samples of RE for all coo locations
s <- as.matrix(A.pred %*% pred_s)

## extract cell values  from covariates, deal with timevarying covariates here
inla_vals <- list()
for(p in 1:nperiods)  inla_vals[[p]] <- cov_vals[[p]] %*% pred_l

# fixed effects by cell
l <- do.call(rbind,inla_vals)

## add together FE and GP components
pred_inla <- s+l

##
totalpredict_time_inla <- proc.time()[3] - ptm




###########################################################
###########################################################
## OPTIONAL MAKE COMPARISONS PLOTS
len <- nrow(pred_inla)/nperiod

# make some useful files first
summ_inla <- cbind(median=(apply(pred_inla,1,median)),sd=(apply(pred_inla,1,sd)))
summ_tmb  <- cbind(median=(apply(pred_tmb,1,median)) ,sd=(apply(pred_tmb,1,sd)))

ras_med_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,1]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
ras_med_tmb   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,1]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")
ras_sdv_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,2]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
ras_sdv_tmb   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,2]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")

pdf(sprintf('/share/geospatial/royburst/sandbox/tmb/inla_compare_real_data/test_%sdegree_mesh_withspeedups2.pdf',as.character(max_edge)), height=10,width=14)


####
#### MAKE TABLE
res <- data.table(st_mesh_nodes    =rep(nrow(epsilon_draws),2))
res[,cores  := rep(ncores,2)]
res[,s_mesh_max_edge  := rep(max_edge,2)]
res[,periods        := c(4,4)]
res[,draws        := c(ndraws,ndraws)]

# time variables
res[,fit_time  := c(fit_time_inla,fit_time_tmb)]
res[,pred_time := c(totalpredict_time_inla,totalpredict_time_tmb)]
res[,pt_tmb_sdreport_time := c(NA,tmb_sdreport_time)]
res[,pt_get_draws_time := c(inla_get_draws_time,tmb_get_draws_time)]

# fe coefficients
for(i in 1:length(res_fit$names.fixed)){
  fn <- res_fit$names.fixed[i]
  res[[paste0('fixedeff_',fn,'_mean')]] <- c(res_fit$summary.fixed$mean[i], SD0$par.fixed[i])
  res[[paste0('fixedeff_',fn,'_sd')]]   <- c(res_fit$summary.fixed$sd[i], sqrt(SD0$cov.fixed[i,i]))
}

# hyperparameters
res[,hyperpar_logtau_mean := c(res_fit$summary.hyperpar[1,1], SD0$par.fixed['logtau']) ]
res[,hyperpar_logtau_sd := c(res_fit$summary.hyperpar[1,2], sqrt(SD0$cov.fixed['logtau','logtau'])) ]

res[,hyperpar_logkappa_mean := c(res_fit$summary.hyperpar[2,1], SD0$par.fixed['logkappa']) ]
res[,hyperpar_logkappa_sd := c(res_fit$summary.hyperpar[2,2], sqrt(SD0$cov.fixed['logkappa','logkappa'])) ]

res[,hyperpar_rho_mean := c(res_fit$summary.hyperpar[3,1], SD0$par.fixed['trho']) ]
res[,hyperpar_rho_sd := c(res_fit$summary.hyperpar[3,2], sqrt(SD0$cov.fixed['trho','trho'])) ]

rr <- data.table(item=colnames(res))
rr <- cbind(rr,t(res))
names(rr) <- c('_','R-INLA','TMB')
rr$diff <- rr[,2]-rr[,3]

library(gridExtra)
library(grid)
grid.table(rr)

####
#### MAKE PLOTS

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

      # residual
      #  tmp$inla<-extract(ras_med_inla[[i]],cbind(tmp$longitude,y=tmp$latitude))
      #  tmp$tmb<-extract(ras_med_tmb[[i]],cbind(tmp$longitude,y=tmp$latitude))
      #  tmp$dat <- tmp$died/tmp$N
      #  tmp$resid_inla <- tmp$dat-tmp$inla
      #  tmp$resid_tmb  <- tmp$dat-tmp$tmb
      #  tmp<-subset(tmp,dat<quantile(tmp$dat,.9))
      #  plot(x=tmp$dat,y=tmp$resid_inla, pch=19,col='red',cex=.1,main='RESIDUALS')
      #  points(x=tmp$dat,y=tmp$resid_tmb, pch=19,col='blue',cex=.1)
    }
}

layout(matrix(1, 1, 1, byrow = TRUE))
####
#### Compare mean and distribution of random effects
summ_gp_tmb  <- t(cbind((apply(epsilon_draws,1,quantile,probs=c(.1,.5,.9)))))
summ_gp_inla <- t(cbind((apply(pred_s,1,quantile,probs=c(.1,.5,.9)))))
  # all time-space random effects

plot_d <- data.table(tmb_median = summ_gp_tmb[,2],inla_median = summ_gp_inla[,2],
                     tmb_low    = summ_gp_tmb[,1],inla_low    = summ_gp_inla[,1],
                     tmb_up     = summ_gp_tmb[,3],inla_up     = summ_gp_inla[,3])

plot_d$period <- factor(rep(1:4,each=nrow(plot_d)/4))
plot_d$loc    <- rep(1:(nrow(plot_d)/4),rep=4)

if(nrow(plot_d)>2500)
  plot_d <- plot_d[sample(nrow(plot_d),2500,replace=F),]


ggplot(plot_d, aes(x=tmb_median,y=inla_median,col=period)) + theme_bw() +
  geom_point() + geom_line(aes(group=loc)) + geom_abline(intercept=0,slope=1,col='red') +
  ggtitle('Posterior Medians of Random Effects at Mesh Nodes, TMB v R-INLA. Connected dots same location different periods. ')

# plot locations where they are different, are they near or far from data?
plot_d[, absdiff := abs(tmb_median-inla_median)]
nodelocs <- do.call("rbind", replicate(4, mesh_s$loc, simplify = FALSE))
biggdiff <- unique(nodelocs[which(plot_d$absdiff>quantile(plot_d$absdiff,prob=0.80)),])

nodelocs <- cbind(nodelocs,plot_d)
if(nrow(nodelocs)>2500)
  nodelocs <- nodelocs[sample(nrow(nodelocs),2500,replace=FALSE),]

par(mfrow=c(2,2))
for(i in 1:4){
  plot(simple_polygon, main='Mesh nodes sized by abs difference TMB and R-INLA')
  points(x=df$longitude[df$period_id==i],y=df$latitude[df$period_id==i], pch=19, cex=0.1)
  points(x=nodelocs$V1[nodelocs$period==i],y=nodelocs$V2[nodelocs$period==i], pch=1, cex=nodelocs$absdiff[nodelocs$period==i]*5, col='red')
}

# catterpillar plot
plot_d <- plot_d[order(period,tmb_median)]
plot_d[,i := seq(1,.N), by = period]
ggplot(plot_d, aes(i, tmb_median, col=i)) + theme_bw() + # [seq(1, nrow(plot_d), 5)]
  geom_linerange(aes(ymin = tmb_low, ymax = tmb_up), col='blue', size=.8, alpha=.3) +
  geom_linerange(aes(x=i,ymin = inla_low, ymax = inla_up), col='red', size=.8, alpha=.3) +
  facet_wrap(~period) +
  ggtitle('Comparison of random effects (10% to 90% quantiles) ... RED == R-INLA ... BLUE == TMB')


dev.off()
