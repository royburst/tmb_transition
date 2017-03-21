




dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/5_spde"))



for(nc in c(10,50,100,500,1000,5000,10000,50000,100000)){
  for(i in 1:10){

    message(nc); message(i)

    qsub <- sprintf("qsub -e /share/temp/sgeoutput/royburst/errors -o /share/temp/sgeoutput/royburst/output -cwd -l mem_free=1G -P proj_geospatial -N job_%nc_%i r_shell.sh inla_tmb_benchmark.R %i %i",nc,i)

    system(qsub)

  }
}
