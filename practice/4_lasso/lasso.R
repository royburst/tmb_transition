# Roy Burstein
# Trying to recreate lasso problem in Gelman et al, BDA3 Ex 14.10
# Bayesian lasso

rm(list=ls())
options(scipen=999)
library(TMB)
library(data.table)
library(glmnet)
library(pomp)

# need to set to same directory as the template file, also pull from git
# Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/4_lasso"))

######################
# load the data
data(state)
d  <- data.table(state.x77)
y  <- d[,'Life Exp',with=FALSE][[1]]
X  <- as.matrix(d[,c('Population','Income','Illiteracy','Murder','HS Grad','Frost','Area'),with=FALSE])
# center and scale X and add an intercept
cs <- function(x){(x-mean(x))/sd(x-mean(x))}
X <- apply(X,2,cs)
y <- cs(y)
X <- cbind(intercept=rep(1,nrow(X)),X)

# Fit in TMB
TMB::compile("lasso.cpp")
dyn.unload(dynlib("lasso"))
dyn.load(dynlib("lasso"))

par   <- list("log_sd"=1,"log_L"=-1,"betas"=rep(1,ncol(X)))
dat   <- list("y"=y, "X"=X) # is lasso

obj  <- MakeADFun(data       = dat,
                  parameters = par,
                  DLL        = "lasso")

# why is the gr going so high?
opt <- nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr) #, control=list(iter.max=10000)) #,trace=1) )


# use a Sobol design to make a matrix
names(opt$par) <- c('log_sd','log_Lambda',colnames(X))
opt$par
opt$objective
