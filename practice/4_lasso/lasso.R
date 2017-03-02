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

# posterior mode using optimizer -- doesnt seem to converge
opt <- nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, trace=2)
names(opt$par) <- c('log_sd','log_Lambda',colnames(X))
opt$par
opt$objective

# use a Sobol design to make a matrix to sample from the posteriors
SimpleLM <- lm(y~-1+X)
WindowLower <- c(log_sd =log((1-0.4) * summary(SimpleLM)$sigma), log_L = log(0.01))
WindowUpper <- c(log_sd =log((1+0.4) * summary(SimpleLM)$sigma), log_L = log(4))
for(i in 1:ncol(X)){
  WindowLower <- append(WindowLower,unname(summary(SimpleLM)$coef[,1] - 4 * summary(SimpleLM)$coef[,2])[i])
  WindowUpper <- append(WindowUpper,unname(summary(SimpleLM)$coef[,1] + 4 * summary(SimpleLM)$coef[,2])[i])
  names(WindowLower)[2+i]= names(WindowUpper)[2+i]=colnames(X)[i]
}
draws = 100000
pd <- sobolDesign(lower = WindowLower,upper = WindowUpper, draws)

# evaluate posterior
lik<-numeric(draws)
for(d in 1:draws){
  if(d%%10000==0) message(d)
  lik[d] <- exp(-1*obj$fn(pd[d,])) # inverse because min for nll
}
pd[which(lik%in%max(lik)),] # max lik

# sample
samps <- 10000
ps <- pd[sample(1:draws,size=1000,replace=T,prob=lik/sum(lik)),]
ps$lambda <- exp(ps$log_L)
ps$sigma  <- exp(ps$log_sd)
ps <- ps[,-c(1,2)]

# posterior summaries
apply(ps,2,median)
apply(ps,2,mean)

pdf('lasso_posterior.pdf'); plot(ps,cex=.1); dev.off()
