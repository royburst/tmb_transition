rm(list=ls())
gc()
options(scipen=999)
.libPaths("/home/j/temp/geospatial/geos_packages")

library(INLA)
library(TMB)
library(data.table)
library(RandomFields)
library(raster)

## need to set to same directory as the template file, also pull from git
## Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/ferizzlez"))
## source("./inla_tmb_compare.R")

## source some functions made for this bit
source('../utils.R')

if( grepl('geos',Sys.info()['nodename'])) INLA:::inla.dynload.workaround()
###############################################################
## SIMULATE AND SET UP THE DATA
## Simulate a surface, this returns a list of useful objects like samples and truth
simobj <- mortsim(nu         = 2               ,  ##  Matern smoothness parameter (alpha in INLA speak)
                  betas      = c(-3,-1,1,1)    ,  ##  Intercept coef and covariate coef For Linear predictors
                  scale      = .1              ,  ##  Matern scale eparameter
                  Sigma2     = (.25) ^ 2       ,  ##  Variance (Nugget)
                  rho        = 0.9             ,  ##  AR1 term
                  l          = 50              ,  ##  Matrix Length
                  n_clusters = 500             ,  ##  number of clusters sampled ]
                  n_periods  = 4               ,  ##  number of periods (1 = no spacetime)
                  mean.exposure.months = 100   ,  ##  mean exposure months per cluster
                  extent = c(0,1,0,1)          ,  ##  xmin,xmax,ymin,ymax
                  ncovariates = 3              ,  ##  how many covariates to include?
                  seed   = NULL                ,
                  returnall=TRUE               ,
                  tvc     = TRUE) # time varying covariates. either way returns an nperiod length list of covariate rasters (if no tvc, they are dups)


## get samples from which to fit
dt <- simobj[["d"]]

## set some stuff up
dt[,id:=1:.N]
coords   <- cbind(dt$x,dt$y)
nperiod  <- length(unique(dt$period))

## MESH For now use same mesh per time point
## TODO CHANGE THIS
data.boundary <- cbind(c(0, 0, 1, 1), c(0, 1, 1, 0))
mesh_s <- inla.mesh.2d(,
                       data.boundary,
                       max.edge=c(0.2,0.2),
                       cutoff=0.05)

nodes <- mesh_s$n ## get number of mesh nodes
spde <- inla.spde2.matern( mesh_s,alpha=2 ) ## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
## ^ this gives us a linear reduction of \Sigma^{-1} as:
## \Sigma^{-1} = \kappa^4 M_0 + 2\kappa^2M_1 + M_2
## M_2 = M_1M_0^{-1}M_1
## Where the Ms are all sparse matrices stored as "dgTMatrix"
## names(spde$param.inla)

## setup prediction mesh needed to get fomr mesh to data locations within tmb function
## get surface to project on to
data.coords = cbind(x=simobj$d$x, x=simobj$d$y)
data.periods <- simobj$d$period

## use inla helper functions to project the spatial effect from mesh points to data points
A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = data.coords,
                           group = data.periods)

## pull covariate(s) at mesh knots
## covs <- raster::extract(simobj$cov.raster.list[[1]],cbind(mesh_s$loc[,1],mesh_s$loc[,2]))
#### NOTE ABOVE: TO CHANGE WHEN AARON SETS UP MESH IN LIKELIHOOD, RIGHT NOW ONLY TAKING FIRST YEARS VALUES, INCORRECT,
    ## LATER ON X_xp WILL BE STRAIGHT FROM dt SO THIS STEP WILL BE UNNECESSARY

## Data to pass to TMB
X_xp = as.matrix(cbind(1, dt[,names(simobj$cov.raster.list[[1]]),with=FALSE]))


Data = list(num_i=nrow(dt),                 ## Total number of observations
            num_s=mesh_s$n,                 ## Number of vertices in SPDE mesh
            num_t=nperiod,                  ## Number of periods
            num_z=1,
            y_i=dt$deaths,                  ## Number of observed deaths in the cluster (N+ in binomial likelihood)
            n_i=dt$exposures,               ## Number of observed exposures in the cluster (N in binomial likelihood)
            t_i=dt$period-1,                ## Sample period ( starting at zero )
            w_i=rep(1,nrow(dt)),
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
                  Epsilon_stz=array(1, c(mesh_s$n, nperiod)))     ## GP locations

##########################################################
### FIT MODEL
## Make object
## Compile
templ <- "model"
TMB::compile(paste0(templ,".cpp"))
dyn.load( dynlib(templ) )


obj <- MakeADFun(data=Data, parameters=Parameters, random="epsilon", hessian=TRUE, DLL=templ)

## Run optimizer
ptm <- proc.time()[3]
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr,
                              lower       =    lower,
                              upper       =    upper,
                              control     =    list(eval.max=1e4, iter.max=1e4, trace=1)))
tmb_fit_time <- proc.time()[3] - ptm


# Get standard errors
## Report0 = obj$report()
ptm <- proc.time()[3]
SD0 = sdreport(obj,getReportCovariance=TRUE)
## fe_var_covar <- SD0$cov.fixed
tmb_sdreport_time <- proc.time()[3] - ptm

##### Prediction
message('making predictions')
#mu    <- c(SD0$par.fixed[names(SD0$par.fixed)=='alpha'],SD0$par.random[names(SD0$par.random)=="epsilon"])
mu    <- c(SD0$value)

sigma <- SD0$cov

### simulate draws
require(MASS)
npar   <- length(mu)
ndraws <- 50

## make sigma symmetric
require(matrixcalc)
i <- 0
while(!is.symmetric.matrix(sigma)){
  sigma <- round(sigma, 10 - i)
  i <- i + 1
}
message(sprintf("rounded sigma to %i decimals to make it symmetric", 10 - i - 1))

## round more or add to diagonal to make sigma pos-def
i <- 0
sigma2 <- sigma
while(!is.positive.definite(sigma2) & i < 6){
  sigma2 <- round(sigma2, 10 - i)
  i <- i + 1
}
if(is.positive.definite(sigma2)){
  sigma <- sigma2
  message(sprintf("rounded sigma to %i decimals to make it pos-def", 10 - i - 1))
}else{
  i <- 0
  while(!is.positive.definite(sigma)){
    sigma <- sigma + diag(1, nrow(sigma))
    i <- i + 1
  }
  message(sprintf("added %i to the diagonal to make sigma pos-def", i))
}

## now we can take draws
draws  <- t(mvrnorm(n=ndraws,mu=mu,Sigma=sigma))
## ^ look for a quicker way to do this..cholesky

## separate out the draws
epsilon_draws <- draws[rownames(draws)=='Epsilon_xt',]
alpha_draws   <- draws[rownames(draws)=='alpha',]

## get surface to project on to
pcoords = cbind(x=simobj$fullsamplespace$x, y=simobj$fullsamplespace$y)
groups_periods <- simobj$fullsamplespace$t

## use inla helper functions to project the spatial effect.
A.pred <- inla.spde.make.A(
  mesh = mesh_s,
  loc = pcoords,
  group = groups_periods)


## values of S at each cell (long by nperiods)
cell_s <- as.matrix(A.pred %*% epsilon_draws)

## extract cell values  from covariates, deal with timevarying covariates here
vals <- list()
for(p in 1:nperiod){
  vals[[p]] <- extract(simobj$cov.raster.list[[p]], pcoords[1:(nrow(simobj$fullsamplespace)/nperiod),])
  vals[[p]] <- (cbind(int = 1, vals[[p]]))
  vals[[p]] <- vals[[p]] %*% alpha_draws # same as getting cell_ll for each time period
}
cell_l <- do.call(rbind,vals)

## add together linear and st components
pred_tmp <- cell_l + cell_s
## make them into time bins
len = nrow(pred_tmp)/nperiod

tmb_totalpredict_time <- proc.time()[3] - ptm

## eras_tmb <- rasterFromXYZT(data.table(pcoords,p=e,t=rep(1:nperiod,each=len)),"p","t")

####################################################################
## fit using inla
A <- inla.spde.make.A(
  mesh = mesh_s,
  loc = as.matrix(dt[, c('x', 'y'),with=F]),
  group = dt$period
)
space = inla.spde.make.index("space",
                             n.spde = spde$n.spde,
                             n.group = nperiod)

design_matrix <- data.frame(int = 1,dt[,simobj$fe.names,with=F])
stack.obs=inla.stack(
  tag='est',
  data=list(died=dt$deaths), ## response
  A=list(A,1), ## proj matrix, not sure what the 1 is for
  effects=list(
    space,
    design_matrix)
)
formula <-
formula(paste0('died ~ -1+int+',
(paste(simobj$fe.names,collapse = ' + ')),
' + f(space,
                   model = spde,
                   group = space.group,
                   control.group = list(model = \'ar1\'))'
))
ptm <- proc.time()[3]
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
                Ntrials = dt$exposures,
                verbose = TRUE,
                keep = TRUE)
inla_fit_time <- proc.time()[3] - ptm

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

## replicate coordinates and years
sspc=simobj$fullsamplespace
coords <- cbind(x=sspc$x,y=sspc$y)

## get samples of s for all coo locations
s <- A.pred %*% pred_s
s <- as.matrix(s)


## predict out linear effects
## extract cell values  from covariates, deal with timevarying covariates here
vals <- list()
for(p in 1:nperiod){
  vals[[p]] <- extract(simobj$cov.raster.list[[p]], pcoords[1:(nrow(simobj$fullsamplespace)/nperiod),])
  vals[[p]] <- (cbind(int = 1, vals[[p]]))
  vals[[p]] <- vals[[p]] %*% pred_l # same as getting cell_ll for each time period
}
l <- do.call(rbind,vals)


#vals=as.matrix(cbind(int=1,sspc[,simobj$fe.name,with=F]))
#l <- vals %*% pred_l

pred_inla <- s+l

## make them into time bins
len = nrow(pred_inla)/nperiod
inla_totalpredict_time <- proc.time()[3] - ptm



###########################################################
### Summarize draws and compare
## make summary plots - median, 2.5% and 97.5%
summ_inla <- cbind(median=(apply(pred_inla,1,median)),sd=(apply(pred_inla,1,sd)))
## make summary plots - median, 2.5% and 97.5%
summ_tmb <- cbind(median=(apply(pred_tmp,1,median)),sd=(apply(pred_tmp,1,sd)))

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
