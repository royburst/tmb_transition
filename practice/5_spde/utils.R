#######
## A function to make a random covariate
makeRandomCovariate <- function(extent = c(0,1,0,1),
                                mean   = 0,
                                sd=.1,
                                l=51,
                                scale=1,
                                ext=T){



  if(ext) {
    extension = abs(extent[2]-extent[1])/3
  } else { extension=0 }

  cmean = 500
  csd  = 500

  while(!(cmean>(mean-1)&cmean<(mean+1))&csd>sd){
    cov.raster<- raster(outer(seq(0,1, l = l),
                              seq(0,1, l = l),
                              FUN = function(x,y) rnorm(1,0,1)*x^(rnorm(1,2,5))-
                                1*extension*(abs(y))^(runif(1,1,4)/runif(1,1,2))),
                        xmn=extent[1]-extension, xmx=extent[2]+extension,
                        ymn=extent[3]-extension, ymx=extent[4]+extension)
    cmean=(mean(as.matrix(cov.raster)))
    csd  =(sd(as.matrix(cov.raster)))

    if(cmean==Inf|csd==Inf|is.na(cmean)|is.na(csd)) {
      cmean = 500
      csd = 500
    }
  }

  RMmodel = RMmatern(nu    = 1,
                     var   = sd,
                     scale = scale)

  z=RFsimulate(model=RMmodel,
               x=coordinates(cov.raster)[,1],
               y=coordinates(cov.raster)[,2])@data[,1]

  cov.raster=cov.raster+rasterFromXYZ(data.frame(x=coordinates(cov.raster)[,1],
                                                 y=coordinates(cov.raster)[,2],
                                                 z=z))
  return(cov.raster)

}





#######
## simulate mortality data

mortsim <- function(
  nu         = 2            ,  #  Matern smoothness parameter (alpha in INLA speak)
  betas      = c(-3,-1,1,1) ,  #  Intercept coef and covariate coef For Linear predictors
  scale      = 2            ,  #  Matern scale eparameter
  Sigma2     = (1) ^ 2      ,  #  Variance (Nugget)
  rho        = 0.9          ,  #  AR1 term
  l          = 51           ,  #  Matrix Length
  n_clusters = 100          ,  #  number of clusters sampled ]
  n_periods  = 4            ,  #  number of periods (1 = no spacetime)
  mean.exposure.months = 100,  #  mean exposure months per cluster
  extent = c(0,1,0,1)       ,  #  xmin,xmax,ymin,ymax
  ncovariates = 3           ,  #  how many covariates to include?
  seed   = NULL             ,
  returnall=TRUE            ) {


  require(fields)
  require(RandomFields)
  #require(ggplot2); theme_set(theme_bw())
  #require(geoR)
  require(data.table)
  require(raster)


  if(length(betas)!=ncovariates+1)
    stop('length(betas)!=ncovariates+1')

  # seed
  if(!is.null(seed)) set.seed(seed)


  # Generate a covariate field, something simple that changes linearly with location
  # cov. raster has a deeper extent..
  cov.raster=stack(makeRandomCovariate())
  if(ncovariates>1)
    for(i in 2:ncovariates)
      cov.raster <- addLayer(cov.raster,makeRandomCovariate())
  names(cov.raster) <- paste0('X',1:ncovariates)

  #plot(cov.raster,main="Theoretical Covariate")
  ## TODO: Incorporate more realistic covariates

  # make a l by l matrix (raster to sample from), exposures, and covariate values
  template<- raster(outer(seq(extent[1],extent[2], l = l),
                          seq(extent[3],extent[4], l = l)))
  samplespace = data.table('x'=rep(coordinates(template)[,1],n_periods),
                           'y'=rep(coordinates(template)[,2],n_periods),
                           't'=rep(seq(1:n_periods),each=l*l))

  # sample frame
  d=data.table( "x"=runif(n_clusters*n_periods, min=extent[1],max=extent[2]),
                "y"=runif(n_clusters*n_periods, min=extent[3],max=extent[4]),
                "exposures"=round(abs(rnorm(n=n_clusters*n_periods,mean=mean.exposure.months,sd=mean.exposure.months/5))), # exposure months
                "period"=rep(1:n_periods,each=n_clusters),
                'int'=1)

  ## TODO: use pop data to get clumpier (and more realistic) sampling locations
  datanames='int'
  for(cov in 1:ncovariates){
    d[,paste0('X',cov)]=raster::extract(cov.raster[[cov]],cbind(d$x,d$y)) # extract covariate value at sampling locations
    datanames=c(datanames,paste0('X',cov))
  }

  if(!is.null(seed)) set.seed(seed)
  # Define the matern object describing the underlying spatial model
  # can have varying underlying SD
  if(length(Sigma2)==1){
    RMmodel = RMmatern(nu    = nu,
                       var   = Sigma2,
                       scale = scale)
    # A different way of doing this (similar to kronecker?)
    # JT calls this the 'innovations' approach, which should be identical to making kronecker covariance matrix
    # Simulate S
    Epsilon1 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      Epsilon1[,t] = RFsimulate(model=RMmodel,
                                x=samplespace$x[samplespace$t==t],
                                y=samplespace$y[samplespace$t==t])@data[,1]
    }

    # Rho
    Epsilon2 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      if(t==1) Epsilon2[,t] = Epsilon1[,t]
      if(t>=2) Epsilon2[,t] = rho * Epsilon1[,t-1] + Epsilon1[,t]
    }

  } else {
    stopifnot(n_periods==length(Sigma2))

    RMmodel = list()
    for(i in 1:n_periods){
      RMmodel[[i]]= RMmatern(nu    = nu,
                             var   = Sigma2[[i]],
                             scale = scale)
    }
    # A different way of doing this (similar to kronecker?)
    # JT calls this the 'innovations' approach, which should be identical to making kronecker covariance matrix
    # Simulate S
    Epsilon1 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      Epsilon1[,t] = RFsimulate(model=RMmodel[[t]],
                                x=samplespace$x[samplespace$t==t],
                                y=samplespace$y[samplespace$t==t])@data[,1]
    }

    # Rho
    Epsilon2 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      if(t==1) Epsilon2[,t] = Epsilon1[,t]
      if(t>=2) Epsilon2[,t] = rho * Epsilon1[,t-1] + Epsilon1[,t]
    }

  }



  # as vector
  S=as.vector(Epsilon2)


  # Make S into an n_periods rasterStack to extract samples layer

  samplespace$S=S
  r.samplespace=  rasterFromXYZT(samplespace,'S','t')
  names(r.samplespace)=paste0('period',1:n_periods)

  # plot(r.samplespace)
  # sample 'clusters' or sample sites from the sample space
  d$S=NA
  for(t in 1:n_periods)
    d$S[d$period==t] =
    raster::extract(r.samplespace[[t]],cbind(d$x[d$period==t],d$y[d$period==t]))


  ## Calculate linear predictor
  d$linpred = as.matrix(d[,datanames,with=FALSE]) %*% betas
  d$p = plogis(d$linpred+d$S) # logit link

  # Simulate deaths
  d$deaths <- rbinom(n_clusters*n_periods,size=d$exposures, prob=d$p)
  d$mr = d$death/d$exposure

  # save also a raster of true P on the surface for model comparison
  true.mr=samplespace

  samplespace$int=1
  samplespace=cbind(samplespace,raster::extract(cov.raster,cbind(samplespace$x,samplespace$y)))
  samplespace$P  =  plogis(as.matrix(samplespace[,datanames,with=FALSE]) %*% betas+samplespace$S)


  r.true.mr=  stack(rasterFromXYZ(as.matrix(samplespace[t==1,c('x','y','P'),with=F])))
  if(n_periods>1)
    for(r in 2:n_periods)
      r.true.mr=addLayer(r.true.mr,
                         rasterFromXYZ(as.matrix(samplespace[t==r,c('x','y','P'),with=F])))


  names(r.true.mr)<-paste0('period',1:n_periods)

  if(returnall==T){
    return(list(d=d,
                r.true.mr=r.true.mr,
                RMmodel=RMmodel,
                cov.raster=cov.raster,
                S.raster=r.samplespace,
                fullsamplespace=samplespace,
                template=template,
                fe.names=datanames[-1],
                l.s.ratio=median(abs(as.vector(as.matrix(samplespace[,datanames,with=FALSE]) %*% betas))/abs(samplespace$S))))
  } else { return(d) }


}





#######
## Helper function for turning an xyzt table into a raster
rasterFromXYZT <- function(table,
                           z,t){
  require(data.table)
  n_periods = length(unique(table[,t]))
  table$t=table[,t]
  res=  stack(rasterFromXYZ(as.matrix(table[t==1,c('x','y',z),with=F])))
  if(n_periods>1)
    for(r in 2:n_periods)
      res=addLayer(res, rasterFromXYZ(as.matrix(table[t==r,c('x','y',z),with=F])))
  return(res)
}
