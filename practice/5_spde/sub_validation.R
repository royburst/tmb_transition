




dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/5_spde"))


mean.exposure.months  <- 1000
n_periods             <- 4
intercept_coef        <- -3
rho                   <- 0.9
num_covariates        <- 3
# make_time_stamp
run_date <- gsub("-","_",Sys.time())
run_date <- gsub(":","_",run_date)
run_date <- gsub(" ","_",run_date)
dir.create(sprintf('/home/j/temp/geospatial/tmb_testing/cluster_out/%s',run_date))


n=1
for(n_clusters in c(50,100,200,1000,2500,5000,10000,20000,30000,40000,50000,60000,70000,80000,100000)){
  for(i in 1:5){
    qsub <- sprintf("qsub -e /share/temp/sgeoutput/royburst/errors -o /share/temp/sgeoutput/royburst/output -cwd -P proj_geospatial -l mem_free=8G -N job_%i r_shell.sh validation.R %i %i %i %i %s %i %i %g %g",n,n,n_clusters,mean.exposure.months,n_periods,run_date,num_covariates,intercept_coef,rho)
    system(qsub)
    n=n+1
  }
}
