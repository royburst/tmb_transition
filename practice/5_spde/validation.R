rm(list=ls())
gc()
options(scipen=999)
.libPaths('/home/j/temp/geospatial/packages')

library(INLA)
library(TMB)
library(data.table)
library(RandomFields)
library(raster)

## need to set to same directory as the template file, also pull from git
## Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
#system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/5_spde"))
## source("./inla_tmb_compare.R")

## source some functions made for this bit
source('utils.R')
iii  <- as.numeric(commandArgs()[3])


###############################################################
## SIMULATE AND SET UP THE DATA
## Simulate a surface, this returns a list of useful objects like samples and truth
simobj <- mortsim(nu         = 2               ,  ##  Matern smoothness parameter (alpha in INLA speak)
                  betas      = c(-3,-1,1,1)        ,  ##  Intercept coef and covariate coef For Linear predictors
                  scale      = .1              ,  ##  Matern scale eparameter
                  Sigma2     = (.25) ^ 2       ,  ##  Variance (Nugget)
                  rho        = 0.9             ,  ##  AR1 term
                  l          = 50             ,  ##  Matrix Length
                  n_clusters = 5000           ,  ##  number of clusters sampled ]
                  n_periods  = 4               ,  ##  number of periods (1 = no spacetime)
                  mean.exposure.months = 1000 ,  ##  mean exposure months per cluster
                  extent = c(0,1,0,1)          ,  ##  xmin,xmax,ymin,ymax
                  ncovariates = 3              ,  ##  how many covariates to include?
                  seed   = NULL                ,
                  returnall=TRUE                ,
                  tvc     = TRUE) # time varying covariates. either way returns an nperiod length list of covariate rasters (if no tvc, they are dups)

######
## pull out some useful data
simdat <- getsimdata(simobj,meshatdatalocs=FALSE,options=1)

#####
## TMB
tmb <- fit_n_pred_TMB(simdata = simdat, fixsigma = FALSE)

#####
## INLA
inla <- fit_n_pred_INLA(simdata = simdat)

######
## validate
v <- rbind(validate(tmb),validate(inla))

######
## write output
fwrite(v,sprintf('/home/j/temp/geospatial/tmb_testing/cluster_out/%i.csv',iii))
