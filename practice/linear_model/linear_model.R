
# Code modified from James Thorson
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
nll_norm(par=c(13,2), data=list('y' = x))
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
# Step 1 -- make and compile template file
compile("linear_model_v1.cpp")

# Step 2 -- build inputs and object
dyn.load( dynlib("linear_model_v1") )
pars = list("mean"=0, "log_sd"=0)
data = list("y_i"=x)
obj = MakeADFun( data=data, parameters=pars, DLL="linear_model_v1")

# Step 3 -- test and optimize
obj$fn(obj$par)
obj$gr(obj$par)
opt = nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr )
opt$diagnostics = data.frame( "name"=names(obj$par), "Est"=obj$par, "final_gradient"=as.vector(obj$gr(opt$par)))
opt$par # estimated parameters
SD = sdreport(obj) # standard errors

















############
# Example 1 -- average CPUE for canary rockfish
############

# Define a new design matrix
X = cbind( "CA"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD<42,1,0), "OR"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD>=42&WCGBTS_Canary_example$BEST_LAT_DD<46,1,0), "WA"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD>=46,1,0), "Pass"=WCGBTS_Canary_example$PASS-1.5 )

# Step 1 -- make and compile template file
compile( "linear_model_v2.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("linear_model_v2") )
Params = list("b_j"=rep(0,ncol(X)), "log_sd"=0)
Data = list( "y_i"=CPUE, "X_ij"=X )
Obj = MakeADFun( data=Data, parameters=Params, DLL="linear_model_v2")

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics = data.frame( "name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors
