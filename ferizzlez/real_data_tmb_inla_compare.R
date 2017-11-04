#### Try to run TMB on real u5m data (one age bin region)

# NOTES: RUN THESE LINES BEFORE LAUNCHING R
# source /homes/imdavis/intel/mkl/bin/mklvars.sh intel64
# export MKL_INTERFACE_LAYER=GNU,LP64
# export MKL_THREADING_LAYER=GNU
# export OMP_NUM_THREADS=10
# numactl --physcpubind=+0-9
# /homes/imdavis/R_mkl_geos/R-3.4.1-mkl_gcc484/R-3.4.1/bin/R

# OR JUST RUN THE SCRIPT:
 ### /homes/imdavis/R_mkl_geos/R-3.4.1-mkl_gcc484/R-3.4.1/bin/R < /homes/royburst/tmb_transition/ferizzlez/real_data_tmb_inla_compare.R --no-save


############### SETUP
rm(list=ls())
gc()
options(scipen=999)
.libPaths('/home/j/temp/geospatial/geos_packages')
##.libPaths('/home/j/temp/geospatial/packages')


library(INLA)
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
require('TMB', lib.loc='/snfs2/HOME/azimmer/R/x86_64-pc-linux-gnu-library/3.3/') # COMPILED WITH METIS
require('mvnfast',lib.loc='/home/j/temp/geospatial/sandbox_packages/')
library(data.table)
library(RandomFields)
library(raster)
library(viridis)
library(ggplot2)

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
ndraws <- 250

# make a chunky mesh or use the original?
max_edge <- .25

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
  mesh_orig <- list(s=mesh_s,t=mesh_t)
  simple_polygon_list <-  load_simple_polygon(gaul_list        = get_gaul_codes('wssa'),
                                                buffer            = 0.4)
  simple_polygon   <- simple_polygon_list[[2]]
  mesh_s <- build_space_mesh(d = df,
                             simple = simple_polygon,
                             max_edge = max_edge,
                             mesh_offset = 2)


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
            options = c(1,1))                 ## option1==1 use priors
                                              ## option2==1 ADREPORT off

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
openmp(25)
config(tape.parallel=0, DLL=templ)

obj <- MakeADFun(data=Data, parameters=Parameters,  map=list(zrho=factor(NA)), random="Epsilon_stz", hessian=TRUE, DLL=templ)


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
##### Make predictions with TMB

# Get standard errors
## Report0 = obj$report()
ptm <- proc.time()[3]
SD0 = TMB::sdreport(obj,getJointPrecision=TRUE, getReportCovariance=TRUE,bias.correct=FALSE)
tmb_sdreport_time <- proc.time()[3] - ptm

# NOTE getReportCovariance returns for ADREPORTED variables which are used from transofrms, if none are reported, just do getJointPrecision
# So lets invert the precision matrix to get the joint covariance
# ALSO LOOK INTO THIS, SO WE DONT HAVE TO INVERT IT MAYBE https://cran.r-project.org/web/packages/sparseMVN/sparseMVN.pdf
ptm <- proc.time()[3]
sigma <- as.matrix(solve(SD0$jointPrecision))
tmb_invert_precision_matrix <- proc.time()[3] - ptm

##### Prediction
message('making predictions')
mu    <- c(SD0$par.fixed,SD0$par.random) #c(SD0$value)
#sigma <- SD0$cov

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
draws  <- t(mvnfast::rmvn(n=ndraws,mu=mu,sigma=sigma))
# sped up with mvnfast https://cran.r-project.org/web/packages/mvnfast/vignettes/mvnfast.html

## separate out the draws
parnames <- c(names(SD0$par.fixed), names(SD0$par.random))
epsilon_draws <- draws[parnames=='Epsilon_stz',]
alpha_draws   <- draws[parnames=='alpha_j',]


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
pred_tmb <- cell_l + cell_s

# save prediction timing
totalpredict_time_tmb <- proc.time()[3] - ptm



#########################################
#########################################
#### fit using INLA
A <- inla.spde.make.A(
  mesh = mesh_s,
  loc = coords,
  group = df$period_id
)
nperiod <- length(unique(df$period_id))
space   <- inla.spde.make.index("space",
                                n.spde = spde$n.spde,
                                n.group = nperiod)

inla.covs <- c(modnames, 'rates') ## add in malaria rates
design_matrix <- data.frame(int = 1, df[, inla.covs, with=F])
stack.obs <- inla.stack(tag='est',
                        data=list(died=df$died), ## response
                        A=list(A,1), ## proj matrix, not sure what the 1 is for
                        effects=list(
                          space,
                          design_matrix)
                        )

formula <- formula(paste0('died ~ -1+int+',
(paste(inla.covs, collapse = ' + ')),
' + f(space, model = spde, group = space.group, control.group = list(model = \'ar1\'))'))

ptm <- proc.time()[3]

inla.setOption("enable.inla.argument.weights", TRUE)
res_fit <- inla(formula,
                data = inla.stack.data(stack.obs),
                control.predictor = list(A = inla.stack.A(stack.obs),
                                         link = 1,
                                         compute = FALSE),
                control.fixed = list(expand.factor.strategy = 'inla'),
                control.inla = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                control.compute=list(config = TRUE),
                family = 'binomial',
                num.threads = 25,
                Ntrials = df$N,
                weights = df$weight,
                verbose = TRUE,
                keep = TRUE)
fit_time_inla <- proc.time()[3] - ptm


#########################################
#########################################
#### predict with INLA
ptm <- proc.time()[3]
draws <- inla.posterior.sample(ndraws, res_fit)

## get parameter names
par_names <- rownames(draws[[1]]$latent)

## index to spatial field and linear coefficient samples
s_idx <- grep('^space.*', par_names)
l_idx <- which(!c(1:length(par_names)) %in% grep('^space.*|Predictor', par_names))

## get samples as matrices
pred_s <- sapply(draws, function (x) x$latent[s_idx])
pred_l <- sapply(draws, function (x) x$latent[l_idx])
rownames(pred_l) <- res_fit$names.fixed

## replicate coordinates and years (from above in tmb predict)
pcoords <- cbind(x=fullsamplespace$x, y=fullsamplespace$y)
groups_periods <- fullsamplespace$t

## use inla helper functions to project the spatial effect.
A.pred <- inla.spde.make.A(mesh = mesh_s,
                           loc = pcoords,
                           group = groups_periods)

## get samples of s for all coo locations
s <- A.pred %*% pred_s
s <- as.matrix(s)


## extract cell values  from covariates, deal with timevarying covariates here
vals <- list()
for(p in 1:nperiods){
  message(p)
  vals[[p]] <- raster::extract(new_cl[[p]], pcoords[1:(nrow(fullsamplespace)/nperiods),])
  vals[[p]] <- (cbind(int = 1, vals[[p]]))
  vals[[p]] <- vals[[p]] %*% pred_l# same as getting cell_ll for each time period
}

l <- do.call(rbind,vals)

pred_inla <- s+l

## make them into time bins
len = nrow(pred_inla)/nperiod
totalpredict_time_inla <- proc.time()[3] - ptm




###########################################################
###########################################################
## MAKE COMPARISONS

# make some useful files first
summ_inla <- cbind(median=(apply(pred_inla,1,median)),sd=(apply(pred_inla,1,sd)))
summ_tmb  <- cbind(median=(apply(pred_tmb,1,median)) ,sd=(apply(pred_tmb,1,sd)))

ras_med_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,1]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
ras_med_tmb   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,1]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")
ras_sdv_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,2]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
ras_sdv_tmb   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,2]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")



pdf(paste0('/share/geospatial/royburst/sandbox/tmb/inla_compare_real_data/test_quarterdegree_mesh_withspeedups.pdf'), height=10,width=14)



####
#### MAKE TABLE
res <- data.table(st_mesh_nodes    =rep(nrow(epsilon_draws),2))
res[,s_mesh_max_edge  := rep(max_edge,2)]
res[,periods        := c(4,4)]

# time variables
res[,fit_time  := c(fit_time_inla,fit_time_tmb)]
res[,pred_time := c(totalpredict_time_inla,totalpredict_time_tmb)]

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

ggplot(plot_d, aes(x=tmb_median,y=inla_median,col=period)) + theme_bw() +
  geom_point() + geom_line(aes(group=loc)) + geom_abline(intercept=0,slope=1,col='red') +
  ggtitle('Posterior Medians of Random Effects at Mesh Nodes, TMB v R-INLA. Connected dots same location different periods. ')

# plot locations where they are different, are they near or far from data?
plot_d[, absdiff := abs(tmb_median-inla_median)]
nodelocs <- do.call("rbind", replicate(4, mesh_s$loc, simplify = FALSE))
biggdiff <- unique(nodelocs[which(plot_d$absdiff>quantile(plot_d$absdiff,prob=0.80)),])

nodelocs <- (cbind(nodelocs,plot_d))

par(mfrow=c(2,2))
for(i in 1:4){
  plot(simple_polygon, main='Mesh nodes sized by abs difference TMB and R-INLA')
  points(x=df$longitude[df$period_id==i],y=df$latitude[df$period_id==i], pch=19, cex=0.1)
  points(x=nodelocs$V1[nodelocs$period==i],y=nodelocs$V2[nodelocs$period==i], pch=1, cex=nodelocs$absdiff[nodelocs$period==i]*5, col='red')
}

# catterpillar plot
plot_d <- plot_d[order(period,tmb_median)]
plot_d$i <- rep(1:(nrow(plot_d)/4),4)
ggplot(plot_d, aes(i, tmb_median, col=i)) + theme_bw() + # [seq(1, nrow(plot_d), 5)]
  geom_linerange(aes(ymin = tmb_low, ymax = tmb_up), col='blue', size=.8, alpha=.3) +
  geom_linerange(aes(x=i,ymin = inla_low, ymax = inla_up), col='red', size=.8, alpha=.3) +
  facet_wrap(~period) +
  ggtitle('Comparison of random effects (10% to 90% quantiles) ... RED == R-INLA ... BLUE == TMB')

# COMPARE DISTRIBUTIONS
#par(mfrow=c(1,2))
#qqnorm(pred_s[777,]);abline(a=0,b=1, col='red')#
#qqplot(pred_s[777,],epsilon_draws[777,])
#lines(c(0,1),c(0,1),col='red')

dev.off()

#write.csv(res,'/share/geospatial/royburst/sandbox/tmb/inla_compare_real_data/tmb_inla_compare_real_data.csv')
