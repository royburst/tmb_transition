# Roy Burstein
# Trying to recreate lasso problem in Gelman et al, BDA3 Ex 14.10
# Bayesian lasso

rm(list=ls())

library(TMB)
library(data.table)
library(glmnet)

# need to set to same directory as the template file, also pull from git
# Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin master'))
setwd(paste0(dir,"/practice/4_lasso"))

######################
# load the data
data(state)
d=data.frame(state.x77)
