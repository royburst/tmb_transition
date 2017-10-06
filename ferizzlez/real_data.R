#### Try to run TMB on real u5m data (one age bin region)

# NOTES: MKL: export OMP_NUM_THREADS = 10; /homes/imdavis/R_mkl_geos/R-3.4.1-mkl_gcc484/R-3.4.1/bin/R
#


############### SETUP
rm(list=ls())
gc()
options(scipen=999)
.libPaths('/home/j/temp/geospatial/geos_packages')
##.libPaths('/home/j/temp/geospatial/packages')


library(INLA)
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
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
                             max_edge = 3,
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
##### Make predictions with TMB

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
pred_tmb <- cell_l + cell_s

# save prediction timing
totalpredict_time <- proc.time()[3] - ptm

## ##########
## plot it ##
## ##########
summ  <- cbind(median=(apply(pred_tmb,1,median)),sd=(apply(pred_tmb,1,sd)))

# make a median raster
ras   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ[,1]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")
pdf('test_real_data_tmb_median.pdf')
plot(ras,maxpixel=1e7)
dev.off()

# make a sd raster
ras   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ[,2]), t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")
pdf('test_real_data_tmb_sd.pdf')
plot(ras,maxpixel=1e5)
dev.off()

#########################################


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
                num.threads = 10,
                Ntrials = df$N,
                weights = df$weight,
                verbose = TRUE,
                keep = TRUE)
inla_fit_time <- proc.time()[3] - ptm


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
inla_totalpredict_time <- proc.time()[3] - ptm


## ##########
## plot it ##
## ##########
summ  <- cbind(median=(apply(pred_inla,1,median)),sd=(apply(pred_inla,1,sd)))

# make a median raster
ras   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ[,1]), t=rep(1:nperiods,each=nrow(pred)/nperiods)),"p","t")
pdf('test_real_data_inla_median.pdf')
plot(ras)
dev.off()

# make a median raster
ras   <- rasterFromXYZT(data.table(pcoords,p=plogis(summ[,2]), t=rep(1:nperiods,each=nrow(pred)/nperiods)),"p","t")
pdf('test_real_data_inla_sd.pdf')
plot(ras)
dev.off()

###########################################################


###########################################################
###########################################################
## TODO this is just copied from tmb_testing at the moment
### Summarize draws and compare
## make summary plots - median, 2.5% and 97.5%
summ_inla <- cbind(median=(apply(pred_inla,1,median)),sd=(apply(pred_inla,1,sd)))
## make summary plots - median, 2.5% and 97.5%
summ_tmb <- cbind(median=(apply(pred_tmb,1,median)),sd=(apply(pred_tmb,1,sd)))

## get error and SD
truth <- qlogis(as.vector(simobj$r.true.mr))


e_tmb  <- summ_tmb[,1]-truth
e_inla <- summ_inla[,1]-truth

m_diff  <- summ_tmb[,1] - summ_inla[,1]
sd_diff <- summ_tmb[,2] - summ_inla[,2]

m_tmb_r   <- rasterFromXYZT(data.table(pcoords,p=summ_tmb[,1], t=rep(1:nperiod,each=len)),"p","t")
m_inla_r  <- rasterFromXYZT(data.table(pcoords,p=summ_inla[,1],t=rep(1:nperiod,each=len)),"p","t")
e_tmb_r   <- rasterFromXYZT(data.table(pcoords,p=e_tmb, t=rep(1:nperiod,each=len)),"p","t")
e_inla_r  <- rasterFromXYZT(data.table(pcoords,p=e_inla,t=rep(1:nperiod,each=len)),"p","t")
sd_tmb_r  <- rasterFromXYZT(data.table(pcoords,p=summ_tmb[,2], t=rep(1:nperiod,each=len)),"p","t")
sd_inla_r <- rasterFromXYZT(data.table(pcoords,p=summ_inla[,2],t=rep(1:nperiod,each=len)),"p","t")
m_diff_r  <- rasterFromXYZT(data.table(pcoords,p=m_diff,t=rep(1:nperiod,each=len)),"p","t")
sd_diff_r <- rasterFromXYZT(data.table(pcoords,p=sd_diff,t=rep(1:nperiod,each=len)),"p","t")

emn <- min(c(e_inla,e_tmb))
emx <- max(c(e_inla,e_tmb))
smn <- min(c(summ_inla[,2],summ_tmb[,2]))
smx <- max(c(summ_inla[,2],summ_tmb[,2]))
mmn <- min(c(summ_inla[,1],summ_tmb[,1],truth))
mmx <- max(c(summ_inla[,1],summ_tmb[,1],truth))


## plot
print('making plots')
require(grDevices)
##pdf(sprintf('mean_error_tmb_inla_%i_clusts_%iexpMths_wo_priors.pdf', n.clust, n.expMths), height=20,width=16)
pdf("plot.pdf", height=20,width=16)

par(mfrow=c(4,3),
    mar = c(3, 3, 3, 9))

truthr<- simobj$r.true.mr
values(truthr) <- truth

## 1
plot(truthr[[1]],main='TRUTH',zlim=c(mmn,mmx))
points(simobj$d$x[simobj$d$period==1],simobj$d$y[simobj$d$period==1])

## 2
plot(as.vector(sd_tmb_r[[1]]),as.vector(sd_inla_r[[1]]), col = rainbow(11)[cut(pcoords[, 1], breaks = 10)], main = "Color by X")
legend("bottomright", legend = unique(cut(pcoords[, 1], breaks = 10)), col = rainbow(11), pch = 16)
abline(a = 0, b = 1)
## 3
plot(as.vector(sd_tmb_r[[1]]),as.vector(sd_inla_r[[1]]), col = rainbow(11)[cut(pcoords[, 2], breaks = 10)], main = "Color by Y")
legend("bottomright", legend = unique(cut(pcoords[, 2], breaks = 10)), col = rainbow(11), pch = 16)
abline(a = 0, b = 1)
## 4
plot(m_tmb_r[[1]],main='MEDIAN TMB',zlim=c(mmn,mmx))
## 5
plot(m_inla_r[[1]],main='MEDIAN INLA',zlim=c(mmn,mmx))
## 6
cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
brks <- c(seq(min(values(m_diff_r)), 0, length = 15), 0, seq(0, max(values(m_diff_r)),length = 15))[-c(15, 16)]
plot(m_diff_r[[1]],main='MEDIAN DIFFERENCE', col = cls, breaks = brks)
## 7
error.values <- range(c(values(e_tmb_r), values(e_inla_r)))
cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
brks <- c(seq(min(error.values), 0, length = 15), 0, seq(0, max(error.values),length = 15))[-c(15, 16)]
plot(e_tmb_r[[1]],main='TMB ERROR',zlim=c(emn,emx), col = cls, breaks = brks)
## 8
cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
brks <- c(seq(min(error.values), 0, length = 15), 0, seq(0, max(error.values),length = 15))[-c(15, 16)]
plot(e_inla_r[[1]],main='INLA ERROR',zlim=c(emn,emx), col = cls, breaks = brks)
## 9
plot(x = e_inla_r[[1]], y = e_tmb_r[[1]], xlab = "inla error", ylab="tmb error");abline(a = 0, b= 1)
## 10
plot(sd_tmb_r[[1]],main='TMB SD',zlim=c(smn,smx))
points(simobj$d$x[simobj$d$period==1],simobj$d$y[simobj$d$period==1])
## 11
plot(sd_inla_r[[1]],main='INLA SD',zlim=c(smn,smx))
points(simobj$d$x[simobj$d$period==1],simobj$d$y[simobj$d$period==1])
## 12
cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
brks <- c(seq(min(values(sd_diff_r)), 0, length = 15), 0, seq(0, max(values(sd_diff_r)),length = 15))[-c(15, 16)]
plot(sd_diff_r[[1]], main='SD DIFFERENCE', col = cls, breaks = brks)
points(simobj$d$x[simobj$d$period==1],simobj$d$y[simobj$d$period==1])

dev.off()
