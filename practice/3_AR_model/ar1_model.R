# Roy Burstein
rm(list=ls())

library(TMB)


# need to set to same directory as the template file, also pull from git
# Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin master'))
setwd(paste0(dir,"/practice/3_AR_model"))

######################
# Simulate data

N       <- 1000 # number of observations
rho     <- 0.9 # temporal autocorrelation
s2      <- .1

# AR data
y    <- numeric(length = N)
y[1] <- rnorm(1,0,s2)
for(i in 2:N) y[i] <- rnorm(1,y[i-1]*rho,sqrt(s2))

# fit in R
arima(y,order=c(1,0,0))


# Fit in TMB
compile("ar1_model.cpp")

# build inputs and object
dyn.load(dynlib("ar1_model"))


obj  <- MakeADFun(data       = list("y"=y),
                  parameters = list("log_sd"=0, 'rho'=0),
                  DLL        = "ar1_model")
