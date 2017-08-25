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
  returnall=TRUE            ,
  tvc      =FALSE ) { # are there time varying covariates?


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

  # Allow for the incorporation of time varying covariates
  covs <- list(cov.raster)
  if(n_periods > 1){
    for(t in 2:n_periods){
      # make ncovariates rasters with t*.1 surface, plus some noise
      if(tvc ){
        temp         <- cov.raster
        values(temp) <- as.vector(temp)+t*0.1+rnorm(length(temp),0,0.001)
        covs[[t]] <-temp
      } else {
        covs[[t]] <- cov.raster
      }
    }
  }

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

  ## Extract covariate values at data points
  datanames='int'
  dees <- list()
  for(per in 1:n_periods){
      dtmp <- d[d$period==per,]
      locs <- cbind(dtmp$x,dtmp$y)
      for(cov in 1:ncovariates){
        dtmp[,paste0('X',cov)]=raster::extract(covs[[per]][[cov]],locs) # extract covariate value at sampling locations
        datanames=c(datanames,paste0('X',cov))
      }
      dees[[per]] <- dtmp
  }
  datanames=unique(datanames)
  d <- do.call(rbind,dees)


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
  #samplespace=cbind(samplespace,raster::extract(cov.raster,cbind(samplespace$x,samplespace$y)))
  dees <- list()
  for(per in 1:n_periods){
      dtmp <- samplespace[samplespace$t==per,]
      locs <- cbind(dtmp$x,dtmp$y)
      for(cov in 1:ncovariates){
        dtmp[,paste0('X',cov)]=raster::extract(covs[[per]][[cov]],locs) # extract covariate value at sampling locations
      }
      dees[[per]] <- dtmp
  }
  samplespace <- do.call(rbind,dees)


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
                cov.raster.list=covs,
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

#######
## Take a simobj and return useful data inputs for modelling
getsimdata <- function(simobj,               # simobj is a list returned from mortsim
                       meshatdatalocs=FALSE, # make the mesh at data locations?
                       options=1){

  ## get samples from which to fit
  dt <- simobj[["d"]]

  ## set some stuff up
  dt[,id:=1:.N]
  coords   <- cbind(dt$x,dt$y)
  nperiod  <- length(unique(dt$period))

  ## MESH
  data.boundary <- cbind(c(0, 0, 1, 1), c(0, 1, 1, 0))
  if(meshatdatalocs){
       mesh_s <- inla.mesh.2d(loc=coords,
                          max.edge=c(0.2,0.2),
                          cutoff=0.05)
  } else {
       mesh_s <- inla.mesh.2d(,
                             data.boundary,
                             max.edge=c(0.2,0.2),
                             cutoff=0.05)
   }
  nodes <- mesh_s$n ## get number of mesh nodes
  spde <- inla.spde2.matern( mesh_s,alpha=2 ) ## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
  ## ^ this gives us a linear reduction of \Sigma^{-1} as:
  ## \Sigma^{-1} = \kappa^4 M_0 + 2\kappa^2M_1 + M_2
  ## M_2 = M_1M_0^{-1}M_1
  ## Where the Ms are all sparse matrices stored as "dgTMatrix"
  ## names(spde$param.inla)

  ## setup prediction mesh needed to get fomr mesh to data locations within tmb function
  ## get surface to project on to
  data.periods <- dt$period

  ## use inla helper functions to project the spatial effect from mesh points to data points
  A.proj <- inla.spde.make.A(mesh  = mesh_s,
                             loc   = coords,
                             group = data.periods)
  ##model frame
  X_xp = as.matrix(cbind(1, dt[,names(simobj$cov.raster.list[[1]]),with=FALSE]))


  Data = list(n_i=nrow(dt),                   ## Total number of observations
              n_x=mesh_s$n,                   ## Number of vertices in SPDE mesh
              n_t=nperiod,                    ## Number of periods
              n_p=ncol(X_xp),                 ## Number of columns in covariate matrix X
  ##            x_s=mesh_s$idx$loc-1,           ## Association of each cluster with a given vertex in SPDE mesh
              c_i=dt$deaths,                  ## Number of observed deaths in the cluster (N+ in binomial likelihood)
              Exp_i=dt$exposures,             ## Number of observed exposures in the cluster (N in binomial likelihood)
              s_i=dt$id-1,                    ## no site specific effect in my model (different sampling locations in time)
              t_i=dt$period-1,                ## Sample period ( starting at zero )
              X_xp=X_xp,                      ## Covariate design matrix
              G0=spde$param.inla$M0,          ## SPDE sparse matrix
              G1=spde$param.inla$M1,          ## SPDE sparse matrix
              G2=spde$param.inla$M2,          ## SPDE sparse matrix
              Aproj = A.proj,                 ## mesh to prediction point projection matrix
              options = options)                 ## option1==1 use priors
              #spde=(spde$param.inla)[c('M1','M2','M3')])

  ## staring values for parameters
  Parameters = list(alpha   =  rep(0,ncol(X_xp)),                     ## FE parameters alphas
                    log_tau_E=1.0,                                    ## log inverse of tau  (Epsilon)
                    log_kappa=0.0,	                                  ## Matern Range parameter
                    rho=0.5,
                    epsilon=matrix(1,ncol=nperiod,nrow=mesh_s$n))     ## GP locations


  return( list( dt      = dt,
                coords  = coords,
                nperiod = nperiod,
                mesh_s  = mesh_s,
                nodes   = nodes,
                spde    = spde,
                A.proj  = A.proj,
                X_xp    = X_xp,
                fullsamplespace = simobj$fullsamplespace,
                cov.rasters = simobj$cov.raster.list,
                tmbdata = Data,
                tmbpars = Parameters ) )
}



######
## Fit and Predict from TMB
fit_n_pred_TMB <- function( templ = "basic_spde", # string name of template
                            cores = 1, # need to debug parallel
                            ndraws = 50,
                            fixsigma = FALSE, # posdef thing, will change cov
                            bias.correct= TRUE,
                            simdata){ # returned from getsimdata

  ### FIT MODEL
  ## Make object
  ## Compile
  message('compiling')
#  templ <- "basic_spde" # _aoz" #spde2
  TMB::compile(paste0(templ,".cpp"))
  dyn.load( dynlib(templ) )

  # unload objects from the list to be easier to work with
  for(n in names(simdata))
    assign(n,simdata[[n]])

  #openmp(cores) # any nyumber other than 1 does not converge or speed up.
  # map to kill certain variables
  lower =  c(rep(-20,dim(X_xp)[2]),rep(-10,2),-0.999)
  upper =  c(rep( 20,dim(X_xp)[2]),rep( 10,2), 0.999)

  # cancel out rho if needed
  mapout <- list()
  if(nperiod == 1){
    lower  <- lower[-1]
    upper  <- upper[-1]
    mapout <- list(rho=factor(NA))
  }
  # make object
  ptm <- proc.time()[3]
  obj <- MakeADFun(data=tmbdata, parameters=tmbpars, map=mapout, random="epsilon", hessian=TRUE, DLL=templ)

  ## Run optimizer
  opt0 <- do.call("nlminb",list(start       =    obj$par,
                                objective   =    obj$fn,
                                gradient    =    obj$gr,
                                lower       =    lower,
                                upper       =    upper,
                                control     =    list(eval.max=1e4, iter.max=1e4, trace=1))) # rel.tol=.01,step.max=10)))

  fit_time <- proc.time()[3] - ptm

  # Get standard errors
  ## Report0 = obj$report()
  ptm <- proc.time()[3]
  SD0 = sdreport(obj,getReportCovariance=TRUE,bias.correct=bias.correct)
#
#  #### TEST PRECIS MATRIX#
#  SD1 = sdreport(obj,getJointPrecision=TRUE,bias.correct=FALSE)
#  sigma<-solve(SD1$jointPrecision)
#  idx <- which(!rownames(sigma)%in%c("log_tau_E", "log_kappa", "rho"))
#  sigma<-sigma[idx,idx]
#  sigma<-as.matrix(sigma)
#  ####################
#
  ## fe_var_covar <- SD0$cov.fixed
  tmb_sdreport_time <- proc.time()[3] - ptm

  ##### Prediction
  message('making predictions')
  mu    <- c(SD0$value)
  sigma <- SD0$cov

  ### simulate draws
  require(MASS)
  npar   <- length(mu)

  ## make sigma symmetric
  require(matrixcalc)
  i <- 0
  while(!is.symmetric.matrix(sigma)){
    sigma <- round(sigma, 10 - i)
    i <- i + 1
  }
  if(i>9) stop('Too much rounding of sigma to make it symmetric, something is wrong.')

  message(sprintf("rounded sigma to %i decimals to make it symmetric", 10 - i - 1))

  ## round more or add to diagonal to make sigma pos-def
  if(!is.positive.definite(sigma)){
    if(fixsigma){
      message('Sigma was not positive definite, tryna fix it')
      i <- 0
      sigma2 <- sigma
      while(!is.positive.definite(sigma2) & i < 6){
        message(i)
        sigma2 <- round(sigma2, 10 - i)
        i <- i + 1
      }
      if(is.positive.definite(sigma2)){
        sigma <- sigma2
        message(sprintf("rounded sigma to %i decimals to make it pos-def", 10 - i - 1))
      }else{
        i <- 0
        message('adding to diag approach')
        while(!is.positive.definite(sigma)){
          message(i)
          sigma <- sigma + diag(1, nrow(sigma))
          i <- i + 1
        }
        message(sprintf("added %i to the diagonal to make sigma pos-def", i))
      }
    } else {
      stop('Sigma was not definite positive and you set fixsigma==FALSE. So this is an error.')
    }
  } else{
    message('Sigma was positive definite hurray')
  }

  ## now we can take draws
  message('Predicting Draws')
  draws  <- t(mvrnorm(n=ndraws,mu=mu,Sigma=sigma))
  ## ^ look for a quicker way to do this..cholesky

  ## separate out the draws
  epsilon_draws <- draws[rownames(draws)=='Epsilon_xt',]
  alpha_draws   <- draws[rownames(draws)=='alpha',]

  ## get surface to project on to
  pcoords = cbind(x=fullsamplespace$x, y=fullsamplespace$y)
  groups_periods <- fullsamplespace$t

  ## use inla helper functions to project the spatial effect.
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = pcoords,
    group = groups_periods)


  ## values of S at each cell (long by nperiods)
  cell_s <- as.matrix(A.pred %*% epsilon_draws)

  ## extract cell values  from covariates, deal with timevarying covariates here
  vals <- list()
  for(p in 1:nperiod){
    vals[[p]] <- extract(cov.rasters[[p]], pcoords[1:(nrow(fullsamplespace)/nperiod),])
    vals[[p]] <- (cbind(int = 1, vals[[p]]))
    vals[[p]] <- vals[[p]] %*% alpha_draws # same as getting cell_ll for each time period
  }
  cell_l <- do.call(rbind,vals)

  ## add together linear and st components
  pred <- cell_l + cell_s

  totalpredict_time <- proc.time()[3] - ptm


  return( list(pred = pred,
                fit_time = fit_time,
                predict_time = totalpredict_time,
                modelobj = SD0,
                model = 'tmb'))

}





#######################
### INLA
## Fit and Predict from TMB
fit_n_pred_INLA <- function( cores = 1,
                            ndraws = 50,
                            simdata,
                            so=simobj){ # returned from getsimdata
  # unload objects from the list to be easier to work with
  for(n in names(simdata))
    assign(n,simdata[[n]])


  A <- inla.spde.make.A(
    mesh = mesh_s,
    loc = as.matrix(dt[, c('x', 'y'),with=F]),
    group = dt$period
  )
  space = inla.spde.make.index("space",
                               n.spde = spde$n.spde,
                               n.group = nperiod)

  design_matrix <- data.frame(int = 1,dt[,so$fe.names,with=F])
  stack.obs=inla.stack(
    tag='est',
    data=list(died=dt$deaths), ## response
    A=list(A,1), ## proj matrix, not sure what the 1 is for
    effects=list(
      space,
      design_matrix)
  )
  formula <-
  formula(paste0('died ~ -1+int+',
  (paste(so$fe.names,collapse = ' + ')),
  ' + f(space,
                     model = spde,
                     group = space.group,
                     control.group = list(model = \'ar1\'))'
  ))
  ptm <- proc.time()[3]
  res_fit <- inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           link = 1,
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla'),
                  control.inla = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                  control.compute=list(config = TRUE),
                  family = 'binomial',
                  num.threads = 10,
                  Ntrials = dt$exposures,
                  verbose = TRUE,
                  keep = TRUE)
  inla_fit_time <- proc.time()[3] - ptm

  ptm <- proc.time()[3]
  draws <- inla.posterior.sample(ndraws, res_fit)

  ## get parameter names
  par_names <- rownames(draws[[1]]$latent)

  ## index to spatial field and linear coefficient samples
  s_idx <- grep('^space.*', par_names)
  l_idx <- match(sprintf('%s.1', res_fit$names.fixed),
                 par_names)


  ## get samples as matrices
  pred_s <- sapply(draws, function (x) x$latent[s_idx])
  pred_l <- sapply(draws, function (x) x$latent[l_idx])
  rownames(pred_l) <- res_fit$names.fixed

  ## replicate coordinates and years
  pcoords = cbind(x=fullsamplespace$x, y=fullsamplespace$y)
  groups_periods <- fullsamplespace$t

  ## use inla helper functions to project the spatial effect.
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = pcoords,
    group = groups_periods)


  ## get samples of s for all coo locations
  s <- A.pred %*% pred_s
  s <- as.matrix(s)


  ## predict out linear effects
  ## extract cell values  from covariates, deal with timevarying covariates here
  vals <- list()
  for(p in 1:nperiod){
    vals[[p]] <- extract(cov.rasters[[p]], pcoords[1:(nrow(fullsamplespace)/nperiod),])
    vals[[p]] <- (cbind(int = 1, vals[[p]]))
    vals[[p]] <- vals[[p]] %*% pred_l # same as getting cell_ll for each time period
  }
  l <- do.call(rbind,vals)


  #vals=as.matrix(cbind(int=1,sspc[,simobj$fe.name,with=F]))
  #l <- vals %*% pred_l

  pred_inla <- s+l

  ## make them into time bins
  inla_totalpredict_time <- proc.time()[3] - ptm


  return( list( pred         = pred_inla,
                fit_time     = inla_fit_time,
                predict_time = inla_totalpredict_time,
                modelobj=res_fit,
                model        = 'inla') )

}


########################################
# compare outputs.
# get mse, absE, relE, morans I (error), coverage, 2 types of correlation, coverage
validate <- function(fit,       # fit object from tmb or inla
                     so=simobj    # the sim object
                     ){

## summarize predictions and truth and errors
summ  <- cbind(median=(apply(fit$pred,1,median)),sd=(apply(fit$pred,1,sd)))
truth <- qlogis(as.vector(so$r.true.mr))
err   <- summ[,1]-truth

# calculate some things of interest
vv_root_mean_squared_error <- mean(err^2)^(1/2)
vv_mean_absolute_error     <- mean(abs(err))
vv_mean_error              <- mean(err)
vv_relative_mean_error     <- mean(err)/mean(truth)
vv_corr_pearson            <- cor(cbind(summ[,1],truth),method='pearson')[1,2]
vv_corr_spearman           <- cor(cbind(summ[,1],truth),method='spearman')[1,2]

# coverage
# a pixel is considered covered if the truth lies with alpha percentage
vv_coverage_95 <- mean(truth > apply(fit$pred,1,quantile,0.025) & truth < apply(fit$pred,1,quantile,0.975))
vv_coverage_80 <- mean(truth > apply(fit$pred,1,quantile,0.100) & truth < apply(fit$pred,1,quantile,0.900))
vv_coverage_50 <- mean(truth > apply(fit$pred,1,quantile,0.250) & truth < apply(fit$pred,1,quantile,0.750))

# timing info
vv_seconds_to_fit     <- fit$fit_time
vv_seconds_to_predict <- fit$predict_time

# numbers of inputs (clusters, time bins, covariates, etc)
vv_time_periods       <- length(unique(so$d$period))
vv_num_data_points    <- nrow(so$d)
vv_mean_ss            <- round(mean(so$d$exposures),0)
vv_num_covariates     <- length(grep('X',names(so$d)))
vv_average_true_p     <- mean(truth)

# parameters # vars are from env (should change this.. )
vv_intercept_coef       <- intercept_coef
vv_rho                  <- rho

if(fit$model=='tmb'){
  vv_rho_diff           <- fit$modelobj$par.fixed[['rho']] - vv_rho
  vv_rho_covered        <- vv_rho > (fit$modelobj$par.fixed[['rho']] - 1.96*sqrt(fit$modelobj$cov.fixed[['rho','rho']])) &
                           vv_rho < (fit$modelobj$par.fixed[['rho']] + 1.96*sqrt(fit$modelobj$cov.fixed[['rho','rho']]))
  vv_intercept_diff     <- fit$modelobj$par.fixed[[1]] - vv_intercept_coef
  vv_intercept_covered  <- vv_intercept_coef > (fit$modelobj$par.fixed[[1]] - 1.96*sqrt(fit$modelobj$cov.fixed[[1,1]])) &
                           vv_intercept_coef < (fit$modelobj$par.fixed[[1]] + 1.96*sqrt(fit$modelobj$cov.fixed[[1,1]]))
}
if(fit$model=='inla'){
  vv_rho_diff           <- summary(fit$modelobj)$hyperpar[3,1] - vv_rho
  vv_rho_covered        <- vv_rho >  summary(fit$modelobj)$hyperpar[3,3] &
                           vv_rho <  summary(fit$modelobj)$hyperpar[3,5]
  vv_intercept_diff     <- summary(fit$modelobj)$fixed['int','mean'] - vv_intercept_coef
  vv_intercept_covered  <- vv_intercept_coef >  summary(fit$modelobj)$fixed['int',3] &
                           vv_intercept_coef <  summary(fit$modelobj)$fixed['int',5]
}


# return a one row datatable with sensible column names
res <- data.table(software=fit$model)
for(var in ls()[grep('vv_',ls())])
  res[[gsub('vv_','',var)]]=get(var)

return(res)
}


########################################
# comparison plots
comparison_plots <- function(filename='plot.pdf',
                             so=simobj,
                             tmb_obj=tmb,
                             inla_obj=inla){



  summ_inla <- cbind(median=(apply(inla_obj$pred,1,median)),sd=(apply(inla_obj$pred,1,sd)))
  ## make summary plots - median, 2.5% and 97.5%
  summ_tmb <- cbind(median=(apply(tmb_obj$pred,1,median)),sd=(apply(tmb_obj$pred,1,sd)))

  nperiod = length(unique(so$d$period))
  len = nrow(tmb_obj$pred)/nperiod

  ## get error and SD
  truth <- qlogis(as.vector(so$r.true.mr))


  e_tmb  <- summ_tmb[,1]-truth
  e_inla <- summ_inla[,1]-truth

  m_diff  <- summ_tmb[,1] - summ_inla[,1]
  sd_diff <- summ_tmb[,2] - summ_inla[,2]

  pcoords = cbind(x=so$fullsamplespace$x, y=so$fullsamplespace$y)

  m_tmb_r   <- rasterFromXYZT(data.table(pcoords,p=summ_tmb[,1], t=rep(1:nperiod,each=len)),"p","t")
  m_inla_r  <- rasterFromXYZT(data.table(pcoords,p=summ_inla[,1],t=rep(1:nperiod,each=len)),"p","t")
  e_tmb_r   <- rasterFromXYZT(data.table(pcoords,p=e_tmb, t=rep(1:nperiod,each=len)),"p","t")
  e_inla_r  <- rasterFromXYZT(data.table(pcoords,p=e_inla,t=rep(1:nperiod,each=len)),"p","t")
  sd_tmb_r  <- rasterFromXYZT(data.table(pcoords,p=summ_tmb[,2], t=rep(1:nperiod,each=len)),"p","t")
  sd_inla_r <- rasterFromXYZT(data.table(pcoords,p=summ_inla[,2],t=rep(1:nperiod,each=len)),"p","t")
  m_diff_r  <- rasterFromXYZT(data.table(pcoords,p=m_diff,t=rep(1:nperiod,each=len)),"p","t")
  sd_diff_r <- rasterFromXYZT(data.table(pcoords,p=sd_diff,t=rep(1:nperiod,each=len)),"p","t")

  emn <- min(c(e_inla,e_tmb))
  emx <- max(c(e_inla,e_tmb))
  smn <- min(c(summ_inla[,2],summ_tmb[,2]))
  smx <- max(c(summ_inla[,2],summ_tmb[,2]))
  mmn <- min(c(summ_inla[,1],summ_tmb[,1],truth))
  mmx <- max(c(summ_inla[,1],summ_tmb[,1],truth))


  ## plot
  print('making plots')
  require(grDevices)
  ##pdf(sprintf('mean_error_tmb_inla_%i_clusts_%iexpMths_wo_priors.pdf', n.clust, n.expMths), height=20,width=16)
  pdf(filename, height=20,width=16)

  par(mfrow=c(4,3),
      mar = c(3, 3, 3, 9))

  truthr<- so$r.true.mr
  values(truthr) <- truth

  ## 1
  plot(truthr[[1]],main='TRUTH',zlim=c(mmn,mmx))
  points(so$d$x[so$d$period==1],so$d$y[so$d$period==1])

  ## 2
  plot(as.vector(sd_tmb_r[[1]]),as.vector(sd_inla_r[[1]]), col = rainbow(11)[cut(pcoords[, 1], breaks = 10)], main = "Color by X")
  legend("bottomright", legend = unique(cut(pcoords[, 1], breaks = 10)), col = rainbow(11), pch = 16)
  abline(a = 0, b = 1)
  ## 3
  plot(as.vector(sd_tmb_r[[1]]),as.vector(sd_inla_r[[1]]), col = rainbow(11)[cut(pcoords[, 2], breaks = 10)], main = "Color by Y")
  legend("bottomright", legend = unique(cut(pcoords[, 2], breaks = 10)), col = rainbow(11), pch = 16)
  abline(a = 0, b = 1)
  ## 4
  plot(m_tmb_r[[1]],main='MEDIAN TMB',zlim=c(mmn,mmx))
  ## 5
  plot(m_inla_r[[1]],main='MEDIAN INLA',zlim=c(mmn,mmx))
  ## 6
  cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
  brks <- c(seq(min(values(m_diff_r)), 0, length = 15), 0, seq(0, max(values(m_diff_r)),length = 15))[-c(15, 16)]
  plot(m_diff_r[[1]],main='MEDIAN DIFFERENCE', col = cls, breaks = brks)
  ## 7
  error.values <- range(c(values(e_tmb_r), values(e_inla_r)))
  cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
  brks <- c(seq(min(error.values), 0, length = 15), 0, seq(0, max(error.values),length = 15))[-c(15, 16)]
  plot(e_tmb_r[[1]],main='TMB ERROR',zlim=c(emn,emx), col = cls, breaks = brks)
  ## 8
  cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
  brks <- c(seq(min(error.values), 0, length = 15), 0, seq(0, max(error.values),length = 15))[-c(15, 16)]
  plot(e_inla_r[[1]],main='INLA ERROR',zlim=c(emn,emx), col = cls, breaks = brks)
  ## 9
  plot(x = e_inla_r[[1]], y = e_tmb_r[[1]], xlab = "inla error", ylab="tmb error");abline(a = 0, b= 1)
  ## 10
  plot(sd_tmb_r[[1]],main='TMB SD',zlim=c(smn,smx))
  points(so$d$x[so$d$period==1],so$d$y[so$d$period==1])
  ## 11
  plot(sd_inla_r[[1]],main='INLA SD',zlim=c(smn,smx))
  points(so$d$x[so$d$period==1],so$d$y[so$d$period==1])
  ## 12
  cls <- c(colorRampPalette(c("blue", "white"))(15), colorRampPalette(c("white", "red"))(15))[-15]
  brks <- c(seq(min(values(sd_diff_r)), 0, length = 15), 0, seq(0, max(values(sd_diff_r)),length = 15))[-c(15, 16)]
  plot(sd_diff_r[[1]], main='SD DIFFERENCE', col = cls, breaks = brks)
  points(so$d$x[so$d$period==1],so$d$y[so$d$period==1])

  dev.off()

}


###########
# plotting funct
plotbenchmarks <- function( x=NULL, #x and y are the variables for axes
                            y=NULL,
                            rd=run_date,
                            justreturncompileddata=FALSE) {
  require(ggplot2)

  message(sprintf('compiling validation runs from  /home/j/temp/geospatial/tmb_testing/cluster_out/%s/',rd))
  d<-data.table()
  for(f in list.files(path = sprintf('/home/j/temp/geospatial/tmb_testing/cluster_out/%s/',rd)))
    if(length(grep('.csv',f))>0)
      d <- rbind(d,fread(sprintf('/home/j/temp/geospatial/tmb_testing/cluster_out/%s/%s',rd,f)))
  dd<-copy(d)

  if(!justreturncompileddata){
    #summarize
    d[,mean_ss:=round(mean_ss,-2)]
    d[,total_time:=seconds_to_predict+seconds_to_fit]
    d<-d[,.(yp=mean(get(y)),number_of_simulations=sum(.N)),by=.(software,get(x))]

    filepath=sprintf('/home/j/temp/geospatial/tmb_testing/cluster_out/%s/%s__%s_plot.pdf',rd,x,y)

    message(sprintf('SAVING PLOT AT /home/j/temp/geospatial/tmb_testing/cluster_out/%s/%s$s_plot.pdf',rd,x,y))
    pdf(filepath,height=8,width=8)
    g<-ggplot(d,aes(x=get,y=yp,col=software))+
      geom_line(size=1)+geom_point(aes(size=number_of_simulations))+
      ylab(y)+xlab(x)+
      theme_bw()
    print(g)
    dev.off()

  } else{
    return(dd)
  }
}
