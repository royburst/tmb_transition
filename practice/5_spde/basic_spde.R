
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
simobj <- mortsim()

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

# pull covariate(s) at mesh knots
covs <- raster::extract(simobj$cov.raster,cbind(mesh_s$loc[,1],mesh_s$loc[,2]))
