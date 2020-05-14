##############################################################################
# custom functions
##############################################################################

#-----------------------------------------------------------------------------
# function to get parameter estimate from single-component exponential fit
#-----------------------------------------------------------------------------
getExpParam = function(x) {
  MASS::fitdistr(x, "exponential")$estimate[1]
}

#-----------------------------------------------------------------------------
# get log likelihood from single-component exponential fit
#-----------------------------------------------------------------------------
getExpLogLik = function(x) {
  MASS::fitdistr(x, "exponential")$loglik[1]
}

#-----------------------------------------------------------------------------
# get Nagelkerke's pseudo r-squared from null and alternative log-likelihoods
#-----------------------------------------------------------------------------
getNR2 = function(ll0, ll, n){
  
  # Cox & Snell / Maximum likelihood pseudo r-squared
  CoxSnell.R2 <- 1 - exp((2 * (ll0 - ll))/n) 
  
  # Nagelkerke's pseudo r-squared
  r2ML.max <- 1 - exp(ll0 * 2/n)
  Nagelkerke.R2 <- CoxSnell.R2/r2ML.max
  
  return(Nagelkerke.R2)
}

#-----------------------------------------------------------------------------
# run exponential mixture deconvolution with specified number of components
#-----------------------------------------------------------------------------
runExpMix = function(x, ncomp, scale){
  # initialize starting estimates
  lmin <- 10^-c(1:(ncomp-1))
  lmax <- 1-sum(lmin)
  lstart <- c(lmin, lmax)
  rstart <- (10^c(0:-(ncomp-1)))^2/scale
  # lstart <- c(.01,0.1,.89)
  # rstart <- c(1,0.01,.0001)/scale
  
  # get estimates
  out <- expRMM_EM(x/scale, lambda=lstart, rate=rstart, epsilon=1e-06, maxit=50000)
}

#-----------------------------------------------------------------------------
# function returns parameter estimates from exponential mixture deconvolution
# 'scale' option used for ensuring convergence
# expRMM_EM(idsites$D2[idsites$D2>0], lambda=c(1), rate=c(.0001), maxit=50000) # testing
#-----------------------------------------------------------------------------
fitExpMix = function(x, scale, mincomp=2, maxcomp=5, iterate=FALSE){
  
  n <- length(x)
  
  # fit single component exponential
  out.1 <- fitdistr(x, "exponential")
  
  if(iterate){
    prevout <- list("loglik"=-1e8)
    nproc <- mincomp:maxcomp
    
    for (i in nproc){
      
      # # initialize starting estimates
      # lmin <- 10^-c(1:(i-1))
      # lmax <- 1-sum(lmin)
      # lstart <- c(lmin, lmax)
      # rstart <- (10^c(0:-(i-1)))^2/scale
      # # lstart <- c(.01,0.1,.89)
      # # rstart <- c(1,0.01,.0001)/scale
      # 
      # # get estimates
      # out <- expRMM_EM(x/scale, lambda=lstart, rate=rstart, epsilon=1e-06, maxit=50000)
      
      out <- runExpMix(x, i, scale)
      
      # stop at i-1 components if loglik_i does not increase by >10
      if (2*out$loglik - 2*prevout$loglik < 9.21){
        out <- prevout
        break
      }
      
      prevout <- out
    }
  } else {
    # out <- expRMM_EM(x/scale, lambda=lstart, rate=rstart, epsilon=1e-06, maxit=50000)
    out <- runExpMix(x, maxcomp, scale)
  }
  
  
  # output as data frame
  lout <- out$lambda
  rout <- 1/out$rate*scale
  nparams <- length(lout)
  nr2 <- getNR2(out.1$loglik, out$loglik, n)
  
  param_labs <- paste0("p", 1:nparams)
  
  data.frame(param=param_labs, 
             lambda=lout, 
             rate=rout, 
             n=nparams, 
             loglik=out$loglik,
             loglik1=out.1$loglik,
             nr2=nr2)
}

#-----------------------------------------------------------------------------
# function simulates mixture of exponentials from exp_fits_wide data
#-----------------------------------------------------------------------------
simDists <- function(data, index){
  
  # get estimates for sample[index]
  testsim <- as.list(data[index,])
  
  # parse lambda and rate estimates
  lest <- unlist(testsim[c(2,4,6)])
  rest <- 1/unlist(testsim[c(3,5,7)])
  
  # simulate distances from mixture and return
  rexpmix(testsim$tot, lambda=lest, rate=rest)
}

#-----------------------------------------------------------------------------
# get chromosome number from input file
#-----------------------------------------------------------------------------
getChrNum <- function(x){
  chr <- gsub("\\..*", "", x)
  chr <- as.numeric(gsub("chr", "", chr))
}

#-----------------------------------------------------------------------------
# function takes distances to two nearest mutations and 4 exponential rates
# and returns rate index with greatest probability
#-----------------------------------------------------------------------------
assignCluster <- function(dist1, dist2, test_rates){
  dist <- min(dist1, dist2)
  match <- which.max(lapply(test_rates, function(x) pexp(dist+1, 1/x)-pexp(dist-1,1/x)))
  return(paste0("c", match))
}