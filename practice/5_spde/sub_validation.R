




dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/5_spde"))



for(i in 1:10){
  qsub <- sprintf("qsub -e /share/temp/sgeoutput/royburst/errors -o /share/temp/sgeoutput/royburst/output -cwd -P proj_geospatial -N job_%i r_shell.sh validation.R %i",i,i)
  system(qsub)
}
