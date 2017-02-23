# Roy Burstein
# Code modified from James Thorson
# Also of help: https://github.com/kaskr/adcomp/wiki/Tutorial
# Running through some intro model fitting using TMB
rm(list=ls())

library(devtools)
library(TMB)

# need to set to same directory as the template file
setwd("/homes/royburst/tmb_transition/practice/linear_model/")

############
# Example 1 -- normal model
############

# Data
x <- rnorm(n=100,mean=10,sd=2)

# Example Probability for example parameters
dnorm(x, mean=20, sd=2 )           # Likelihood for each datum
sum(dnorm(x, mean=14, sd=2, log=TRUE ))     # Log-likelihood for all data
sum(dnorm(x, mean=12, sd=2, log=TRUE ))     # Log-likelihood for all data
sum(dnorm(x, mean=10, sd=2, log=TRUE ))     # Log-likelihood for all data
sum(dnorm(x, mean=8, sd=2, log=TRUE ))     # Log-likelihood for all data

#######################################################
######### Method 1 -- Pre-made functions in R
lm <- lm(x ~ 1)
summary(lm)

#######################################################
######### Method 2 -- Optimize using R
##Step 1 -- define function - negll, so looking for most negative
nll_norm <- function(par, data){
  # Parameters
  mean_hat <- par[1]
  sd_hat   <- par[2]
  # Log-likelihood
  ll_i <- dnorm(data$y, mean=mean_hat, sd=sd_hat, log=TRUE)
  nll  <- -1 * sum(ll_i)
  return(nll)
}

# example of what this returns
nll_norm(par=c(12,2), data=list('y' = x))
nll_norm(par=c(10,2), data=list('y' = x))
nll_norm(par=c(8,2),  data=list('y' = x))


## step 2 -- minimize negative log-likelihood to estimate parameters
data  <- list('y' = x)
start <- c(10,1) # starting value
opt   <- optim(par=start, fn=nll_norm, data=data, lower=c(-Inf,0.01), upper=Inf, method="L-BFGS-B", hessian=TRUE)
print(opt$par) # mean and SD
print(sqrt(diag(solve(opt$hessian)))) # standard errors


##################################################
###### Method 3 -- Optimize using TMB
# compile template file
# file must be in the same directory as this script
compile("linear_model_v1.cpp")

# build inputs and object
dyn.load( dynlib("linear_model_v1") )

# make an R object representing the cpp function
# data and parameter names must correspond to those in the cpp script
obj  <- MakeADFun(data       = list("y_i"=x),
                  parameters = list("mean"=0, "log_sd"=0),
                  DLL        = "linear_model_v1")

# test, should correspond with the nll_norm() calls above ^^
obj$fn(list("mean"=12, "log_sd"=log(2)))  # Call TMB function for one value
obj$fn(list("mean"=10, "log_sd"=log(2)))  # Call TMB function for one value
obj$fn(list("mean"=8,  "log_sd"=log(2)))  # Call TMB function for one value
# test gradient
obj$gr(list("mean"=12, "log_sd"=log(2))) # First order derivatives (w.r.t mu and sigma)
obj$gr(list("mean"=10, "log_sd"=log(2))) # First order derivatives (w.r.t mu and sigma) # SHOULD BE CLOSEST TO ZERO
obj$gr(list("mean"=8, "log_sd"=log(2))) # First order derivatives (w.r.t mu and sigma)
# test hessian
obj$he(list("mean"=10, "log_sd"=log(2))) # Second order derivatives, i.e. the Hessian matrix

# optimize
opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr) # will marginal gradient close to 0
opt$diagnostics <- data.frame( "name"=names(obj$par), "Est"=c, "final_gradient"=as.vector(obj$gr(opt$par)))
opt$par # estimated parameters
print(sqrt(diag(solve(obj$he(opt$par))))) # standard errors using the hessian
sdreport(obj) # standard errors the easier way




##################################################
##################################################
# Example 2 -- normal model with fixed effects, using TMB of course
############

# Simulate some basic data
X <- cbind(rep(1,100),runif(n=100))
y <- 2*X[,1] + 10*X[,2] + rnorm(n=100,mean=0,sd=2)
summary(lm(y~X[,2]))# test

# compile template file
compile("linear_model_v2.cpp")

# build inputs and object
dyn.load(dynlib("linear_model_v2"))  # dynamically link the C++ code
Params <- list("b_j"=rep(0,ncol(X)), "log_sd"=0)
Data   <- list("y_i"=y, "X_ij"=X)
Obj    <- MakeADFun(data=Data, parameters=Params, DLL="linear_model_v2")

# test and optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics = data.frame( "name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport(Obj) # standard errors
