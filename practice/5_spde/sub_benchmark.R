




dir <- paste0("/homes/",Sys.info()['user'],"/tmb_transition")
system(paste0('cd ',dir,'\ngit pull origin develop'))
setwd(paste0(dir,"/practice/5_spde"))



for(nc in c(10,50,100,500,1000,5000,10000,50000,100000)){
  for(i in 1:10){

    message(nc); message(i)

    qsub <- sprintf("qsub -e /share/temp/sgeoutput/royburst/errors -o /share/temp/sgeoutput/royburst/output -cwd -l mem_free=8G -P proj_geospatial -N job_%i_%i r_shell.sh inla_tmb_benchmark.R %i %i",nc,i,nc,i)

    system(qsub)

  }
}


library(data.table)
### once it runs make a  plot
d<-data.table()
for(f in list.files(path = "./benchmarkoutput"))
  d <- rbind(d,fread(paste0('./benchmarkoutput/',f)))

d <- d[,.(  inla_fit_time=mean(inla_fit_time),
            inla_totalpredict_time=mean(inla_totalpredict_time),
            tmb_fit_time=mean(tmb_fit_time),
            tmb_sdreport_time=mean(tmb_sdreport_time),
            tmb_totalpredict_time=mean(tmb_totalpredict_time),
            sims = .N),
        by = n_clusters]
d <- d[order(n_clusters),]

pdf('fittimes.pdf')
plot(d$n_clusters,d$inla_fit_time+d$inla_totalpredict_time,ylab='Fit time in seconds',type='l',xlab='numclusters',col='red')
lines(d$n_clusters,d$tmb_fit_time+d$tmb_totalpredict_time,col='blue')
dev.off()
