# need to set to same directory as the template file, also pull from git
# Clone the git directory to your H drive and this should work for anyone
dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
setwd(paste0(dir,"/practice/1_linear_model"))
# git pull
system(paste0('cd ',dir,'\ngit pull origin develop'))

library(TMB)
set.seed(123)
x <- seq(0,10,length=50001)
data <- list(Y=rnorm(length(x))+x,x=x)
parameters <- list(a=0,b=0,logSigma=0)
compile( "lin_par.cpp" )
dyn.load(dynlib("lin_par"))
obj <- MakeADFun(data,parameters,DLL="lin_par")
obj$hessian <- TRUE
opt <- do.call("optim",obj)

#########################
# Experiment with parallelization
#########################

# Parallelized experimentations
#ben <- benchmark(obj, n=1e3, cores=1:4)
#plot(ben)

# Parallel
ben <- benchmark(obj, n=1, cores=1:30, expr=expression(do.call("optim",obj)))
png( file="Benchmark.png", width=6, height=6, res=200, units="in")
  plot(ben,ylim=c(0,30))
dev.off()

##########################
# Help file version
##########################

runExample("lin_par",thisR=TRUE) ## Create obj
ben <- benchmark(obj,n=100,cores=1:4)
plot(ben)
ben <- benchmark(obj,n=10,cores=1:4,expr=expression(do.call("optim",obj)))
plot(ben)
