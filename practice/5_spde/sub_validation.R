




dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/5_spde"))


mean.exposure.months  <- 1000
meshatdatalocs         <- TRUE
intercept_coef        <- -3
rho                   <- 0.9
num_covariates        <- 3
# make_time_stamp
run_date <- gsub("-","_",Sys.time())
run_date <- gsub(":","_",run_date)
run_date <- gsub(" ","_",run_date)
dir.create(sprintf('/home/j/temp/geospatial/tmb_testing/cluster_out/%s',run_date))

message(run_date)

n=1
for(n_clusters in c(50,100,200,1000,2500,5000,10000,20000,30000,40000,50000,60000,70000,80000,100000)){
  for(i in 1:5){
    qsub <- sprintf("qsub -e /share/temp/sgeoutput/royburst/errors -o /share/temp/sgeoutput/royburst/output -cwd -P proj_geospatial -l mem_free=8G -N job_%i r_shell.sh validation.R %i %i %i %i %s %i %g %g",n,n,n_clusters,mean.exposure.months,meshatdatalocs,run_date,num_covariates,intercept_coef,rho)
    system(qsub)
    n=n+1
  }
}


#############################
## Once that has run compile them all
d<-plotbenchmarks(rd=run_date,justreturncompileddata=TRUE) # just get the data

plotbenchmarks(rd=run_date,x="num_data_points",y="root_mean_squared_error")
plotbenchmarks(rd=run_date,x="num_data_points",y="coverage_95")
plotbenchmarks(rd=run_date,x="num_data_points",y="seconds_to_fit" )
plotbenchmarks(rd=run_date,x="num_data_points",y="rho_diff" )
plotbenchmarks(rd=run_date,x="num_data_points",y="corr_spearman" )
plotbenchmarks(rd=run_date,x="num_data_points",y="total_time" )
plotbenchmarks(rd=run_date,x="num_data_points",y="intercept_covered" )
plotbenchmarks(rd=run_date,x="num_data_points",y="rho_covered" )
