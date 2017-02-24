# Roy Burstein
# Code modified from James Thorson
rm(list=ls())

library(TMB)
library(lme4)

# need to set to same directory as the template file, also pull from git
# Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
setwd(paste0(dir,"/practice/2_linear_mixed_model"))
system(paste0('cd ',dir,'\ngit pull origin master'))

######################
# Simulate data

N       <- 1000 # number of observations
J       <- 10   # number of groups
j       <- rep(1:J, each=N/J)  # group id
theta_j <- rnorm(length(unique(j)), mean=0, sd=1) # group effect
X0      <- 3

y_i = theta_j[j] + X0 + rnorm(N, mean=0, sd=1)

######################
# Run in R for example

summary(lmer(y_i ~ 1|factor(j)))

######################
# Run in TMB
######################

# Compile model
compile("linear_mixed_model.cpp")

# Build inputs
d <- list("n_data"    = N,
          "n_factors" = J,
          "j"         = j-1, # cpp indexes starting at zero
          "y"         = y_i)
p <- list("X0"      = -10,
          "log_SD0" =   2,
          "log_SDJ" =   2,
          "J"       = rep(0,J))
r <- c("J")
#if( Use_REML==TRUE ) Random = union( Random, "X0") # not sure what this is, look at JT notes

# Build object
dyn.load(dynlib("linear_mixed_model"))
obj <- MakeADFun(data=d, parameters=p, random=r, DLL="linear_mixed_model")

# Prove that function and gradient calls work
obj$fn(obj$par)
obj$gr(obj$par)

# Optimize
opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, control=list("trace"=1))

# Get reporting and SEs
obj$report()
sdreport(obj)
