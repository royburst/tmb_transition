
sdreport <- function(obj,par.fixed=NULL,hessian.fixed=NULL,getJointPrecision=FALSE,bias.correct=FALSE,
                     bias.correct.control=list(sd=FALSE, split=NULL, nsplit=NULL), ignore.parm.uncertainty = FALSE,
                     getReportCovariance=TRUE, skip.delta.method=FALSE){
  if(is.null(obj$env$ADGrad) & (!is.null(obj$env$random)))
    stop("Cannot calculate sd's without type ADGrad available in object for random effect models.")
  ## Make object to calculate ADREPORT vector
  obj2 <- MakeADFun(obj$env$data,
                    obj$env$parameters,
                    type = "ADFun",
                    ADreport = TRUE,
                    DLL = obj$env$DLL,
                    silent = obj$env$silent)
  r <- obj$env$random
  ## Get full parameter (par), Fixed effects parameter (par.fixed)
  ## and fixed effect gradient (gradient.fixed)
  if(is.null(par.fixed)){ ## Parameter estimate not specified - use best encountered parameter
    par <- obj$env$last.par.best
    if(!is.null(r))par.fixed <- par[-r] else par.fixed <- par
    gradient.fixed <- obj$gr(par.fixed)
  } else {
    gradient.fixed <- obj$gr(par.fixed) ## <-- updates last.par
    par <- obj$env$last.par
  }
  ## In case of empty parameter vector:
  if(length(par.fixed)==0) ignore.parm.uncertainty <- TRUE
  ## Get Hessian wrt. fixed effects (hessian.fixed) and check if positive definite (pdHess).
  if(ignore.parm.uncertainty){
      hessian.fixed <- NULL
      pdHess <- TRUE
      Vtheta <- matrix(0, length(par.fixed), length(par.fixed))
  } else {
      if(is.null(hessian.fixed)){
          hessian.fixed <- optimHess(par.fixed,obj$fn,obj$gr) ## Marginal precision of theta.
      }
      pdHess <- !is.character(try(chol(hessian.fixed),silent=TRUE))
      Vtheta <- try(solve(hessian.fixed),silent=TRUE)
      if(is(Vtheta, "try-error")) Vtheta <- hessian.fixed * NaN
  }
  ## Get random effect block of the full joint Hessian (hessian.random) and its
  ## Cholesky factor (L)
  if(!is.null(r)){
    hessian.random <- obj$env$spHess(par,random=TRUE)   ## Conditional prec. of u|theta
    L <- obj$env$L.created.by.newton
    if(!is.null(L)){ ## Re-use symbolic factorization if exists
      message("ABOUT TO UPDATE CHOLESKY")
      TMB:::updateCholesky(L,hessian.random)
      hessian.random@factors <- list(SPdCholesky=L)
    }
  }
  ## Get ADreport vector (phi)
  phi <- try(obj2$fn(par), silent=TRUE)    ## NOTE_1: obj2 forward sweep now initialized !
  if(is.character(phi) | length(phi)==0){
      phi <- numeric(0)
  }
  ADGradForward0Initialized <- FALSE
  ADGradForward0Initialize <- function() { ## NOTE_2: ADGrad forward sweep now initialized !
      obj$env$f(par, order = 0, type = "ADGrad")
      ADGradForward0Initialized <<- TRUE
  }
  doDeltaMethod <- function(chunk=NULL){
      ## ======== Determine case
      ## If no random effects use standard delta method
      simpleCase <- is.null(r)
      if(length(phi)==0){ ## Nothing to report
          simpleCase <- TRUE
      } else { ## Something to report - get derivatives
          if(is.null(chunk)){ ## Do all at once
              Dphi <- obj2$gr(par)   ##### NOTE THIS IS WHATS HANGING
          } else {
              ## Do *chunk* only
              ## Reduce to Dphi[chunk,] and phi[chunk]
              w <- rep(0, length(phi))
              phiDeriv <- function(i){
                  w[i] <- 1
                  obj2$env$f(par, order=1, rangeweight=w, doforward=0) ## See NOTE_1
              }
              Dphi <- t( sapply(chunk, phiDeriv) )
              phi <- phi[chunk]
          }
          if(!is.null(r)){
              Dphi.random <- Dphi[,r,drop=FALSE]
              Dphi.fixed <- Dphi[,-r,drop=FALSE]
              if(all(Dphi.random==0)){ ## Fall back to simple case
                  simpleCase <- TRUE
                  Dphi <- Dphi.fixed
              }
          }
      }
      ## ======== Do delta method
      ## Get covariance (cov)
      if(simpleCase){
          if(length(phi)>0){
              cov <- Dphi %*% Vtheta %*% t(Dphi)
          } else cov <- matrix(,0,0)
      } else {
          tmp <- solve(hessian.random,t(Dphi.random))
          tmp <- as.matrix(tmp)
          term1 <- Dphi.random%*%tmp ## first term.
          if(ignore.parm.uncertainty){
              term2 <- 0
          } else {
              ## Use columns of tmp as direction for reverse mode sweep
              f <- obj$env$f
              w <- rep(0, length(par))
              if(!ADGradForward0Initialized) ADGradForward0Initialize()
              reverse.sweep <- function(i){
                  w[r] <- tmp[,i]
                  -f(par, order = 1, type = "ADGrad", rangeweight = w, doforward=0)[-r]
              }
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
              A <- t(do.call("cbind",lapply(seq_along(phi), reverse.sweep))) + Dphi.fixed

              ptm <- proc.time()[3]
              A <- t(do.call("cbind",lapply(seq_along(phi), reverse.sweep))) + Dphi.fixed

              proc.time()[3] - ptm

              A    <- Dphi.fixed
              for(i in 1:length(temp))
                A[i,] <- A[i,] + temp[[i]]

              temp <- list()
              for(i in 1:70000){
                temp[[i]] <- rnorm(10)
              }
              term2 <- A %*% (Vtheta %*% t(A)) ## second term
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========## ========
          }
          cov <- term1 + term2
      }
      ##list(phi=phi, cov=cov)
      cov
  }
  if (!skip.delta.method) {
      if (getReportCovariance) { ## Get all
          cov <- doDeltaMethod()
          sd <- sqrt(diag(cov))
      } else {
          tmp <- lapply(seq_along(phi), doDeltaMethod)
          sd <- sqrt(unlist(tmp))
          cov <- NA
      }
  } else {
      sd <- rep(NA, length(phi))
      cov <- NA
  }
  ## Output
  ans <- list(value=phi,sd=sd,cov=cov,par.fixed=par.fixed,
              cov.fixed=Vtheta,pdHess=pdHess,
              gradient.fixed=gradient.fixed)
  ## ======== Calculate bias corrected random effects estimates if requested
  if(bias.correct){
      epsilon <- rep(0,length(phi))
      names(epsilon) <- names(phi)
      parameters <- obj$env$parameters
      parameters$TMB_epsilon_ <- epsilon ## Appends to list without changing attributes
      doEpsilonMethod <- function(chunk = NULL) {
          if(!is.null(chunk)) { ## Only do *chunk*
              mapfac <- rep(NA, length(phi))
              mapfac[chunk] <- chunk
              parameters$TMB_epsilon_ <- updateMap(parameters$TMB_epsilon_,
                                                   factor(mapfac) )
          }
          obj3 <- MakeADFun(obj$env$data,
                            parameters,
                            random = obj$env$random,
                            checkParameterOrder = FALSE,
                            DLL = obj$env$DLL,
                            silent = obj$env$silent)
          ## Get good initial parameters
          obj3$env$start <- c(par, epsilon)
          obj3$env$random.start <- expression(start[random])
          ## Test if Hessian pattern is un-changed
          h <- obj$env$spHess(random=TRUE)
          h3 <- obj3$env$spHess(random=TRUE)
          pattern.unchanged <- identical(h@i,h3@i) & identical(h@p,h3@p)
          ## If pattern un-changed we can re-use symbolic Cholesky:
          if(pattern.unchanged){
              if(!obj$env$silent)
                  cat("Re-using symbolic Cholesky\n")
              obj3$env$L.created.by.newton <- L
          } else {
              if( .Call("have_tmb_symbolic", PACKAGE = "TMB") )
                  runSymbolicAnalysis(obj3)
          }
          if(!is.null(chunk)) epsilon <- epsilon[chunk]
          par.full <- c(par.fixed, epsilon)
          i <- (1:length(par.full)) > length(par.fixed) ## epsilon indices
          grad <- obj3$gr(par.full)
          Vestimate <-
              if(bias.correct.control$sd) {
                  ## requireNamespace("numDeriv")
                  hess <- numDeriv::jacobian(obj3$gr, par.full)
                  -hess[i,i] + hess[i,!i] %*% Vtheta %*% hess[!i,i]
              } else
                  matrix(NA)
          estimate <- grad[i]
          names(estimate) <- names(epsilon)
          list(value=estimate, sd=sqrt(diag(Vestimate)), cov=Vestimate)
      }
      nsplit <- bias.correct.control$nsplit
      if(is.null(nsplit)) {
          split <- bias.correct.control$split
      } else {
          split <- split(seq_along(phi),
                         cut(seq_along(phi), nsplit))
      }
      if( is.null( split ) ){ ## Get all
          ans$unbiased <- doEpsilonMethod()
      } else {
          tmp <- lapply(split, doEpsilonMethod)
          m <- if (bias.correct.control$sd)
                   length(phi) else 1
          ans$unbiased <- list(value = rep(NA, length(phi)),
                               sd    = rep(NA, m),
                               cov   = matrix(NA, m, m))
          for(i in seq_along(split)) {
              ans$unbiased$value[ split[[i]] ] <- tmp[[i]]$value
              if (bias.correct.control$sd) {
                  ans$unbiased$sd   [ split[[i]] ] <- tmp[[i]]$sd
                  ans$unbiased$cov  [ split[[i]],
                                      split[[i]] ] <- tmp[[i]]$cov
              }
          }
      }
  }
  ## ======== Find marginal variances of all random effects i.e. phi(u,theta)=u
  if(!is.null(r)){
    if(is(L,"dCHMsuper")){ ## Required by inverse subset algorithm
      ihessian.random <- .Call("tmb_invQ", L, PACKAGE = "TMB")
      iperm <- invPerm(L@perm+1L)
      diag.term1 <- diag(ihessian.random)[iperm]
      if(ignore.parm.uncertainty){
          diag.term2 <- 0
      } else {
          f <- obj$env$f
          w <- rep(0, length(par))
          if(!ADGradForward0Initialized) ADGradForward0Initialize()
          reverse.sweep <- function(i){
              w[i] <- 1
              f(par, order = 1, type = "ADGrad", rangeweight = w, doforward=0)[r]
          }
          nonr <- setdiff(seq_along(par), r)
          tmp <- sapply(nonr,reverse.sweep)
          if(!is.matrix(tmp)) ## Happens if length(r)==1
              tmp <- matrix(tmp, ncol=length(nonr) )
          A <- solve(hessian.random, tmp)
          diag.term2 <- rowSums((A %*% Vtheta)*A)
      }
      ans$par.random <- par[r]
      ans$diag.cov.random <- diag.term1 + diag.term2
      if(getJointPrecision){ ## Get V(u,theta)^-1
          if(length(par.fixed) == 0) {
              ans$jointPrecision <- hessian.random
          }
          else if (!ignore.parm.uncertainty) {
              G <- hessian.random %*% A
              G <- as.matrix(G) ## Avoid Matrix::cbind2('dsCMatrix','dgeMatrix')
              M1 <- cbind2(hessian.random,G)
              M2 <- cbind2(t(G), as.matrix(t(A)%*%G)+hessian.fixed )
              M <- rbind2(M1,M2)
              M <- forceSymmetric(M,uplo="L")
              dn <- c(names(par)[r],names(par[-r]))
              dimnames(M) <- list(dn,dn)
              p <- invPerm(c(r,(1:length(par))[-r]))
              ans$jointPrecision <- M[p,p]
          }
          else {
              warning("ignore.parm.uncertainty ==> No joint precision available")
          }
      }
    } else {
      warning("Could not report sd's of full randomeffect vector.")
    }
  }
  ## Copy a few selected members of the environment 'env'. In
  ## particular we need the 'skeleton' objects that allow us to put
  ## results back in same shape as original parameter list.
  ans$env <- new.env(parent = emptyenv())
  ans$env$parameters <- obj$env$parameters
  ans$env$random <- obj$env$random
  class(ans) <- "sdreport"
  ans
}

##' Extract parameters, random effects and reported variables along
##' with uncertainties and optionally Chi-square statis
