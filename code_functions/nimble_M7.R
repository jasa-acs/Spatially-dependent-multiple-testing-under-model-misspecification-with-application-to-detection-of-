#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: nimble code
# Compile the nimble models for fitting the simulation study
#===============================================================

Hmat <- cbind(rep(1,M), matrix( rnorm(M*30), ncol = 30 )) # dummy EOF matrix

#===============================================================
# M7: EOF structure with skew-t discrepancy, 30 EOFs
#===============================================================

# User-defined nimbleFunction for the skewt distribution
dskewt <- nimbleFunction(
  run = function(x = double(0), mu = double(0, default = 0), 
                 sigma = double(0, default = 1),
                 delta = double(0, default = 0),
                 nu = double(0, default = 0.01),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    
    xi <- mu - (sigma/sqrt(1 - 2*delta^2/3.141593))*sqrt(2/3.141593)*delta
    omega <- sigma/sqrt(1 - 2*delta^2/3.141593)
    alpha <- delta/sqrt(1 - delta^2)
    df <- 1/nu
    
    z <- (x-xi)/omega
    part1 <- dt(x = z, df = df, log = TRUE )
    part2 <- pt(q = alpha*z*sqrt((df+1)/(df + z^2)), df = df + 1, log.p = TRUE )
    logdens <- log(2/omega) + part1 + part2
    
    if( log ){ return( logdens ) }
    else{ return( exp( logdens ) ) }
    
  })

rskewt <- nimbleFunction(
  run = function(n = integer(0, default = 1), mu = double(0, default = 0), 
                 sigma = double(0, default = 1),
                 delta = double(0, default = 0),
                 nu = double(0, default = 0.01) ) {
    
    returnType(double(0))
    
    xi <- mu - (sigma/sqrt(1 - 2*delta^2/3.141593))*sqrt(2/3.141593)*delta
    omega <- sigma/sqrt(1 - 2*delta^2/3.141593)
    alpha <- delta/sqrt(1 - delta^2)
    df <- 1/nu
    
    if(n != 1) print("rskewt only allows n = 1; using n = 1.")
    u1 <- rnorm( n = 1, mean = 0, sd = 1 )
    u2 <- rnorm( n = 1, mean = 0, sd = 1 )
    if ( u2 > alpha * u1 ){ 
      u1 <- -u1
    }
    z <- omega*u1
    v <- rchisq(n = 1, df = df)/df
    y <- z/sqrt(v) + xi
    
    return(y)
  })

registerDistributions( list(dskewt = list(
  BUGSdist = "dskewt(mu, sigma, delta, nu)"
)))

# Create the model
nim_code <- nimbleCode({
  # likelihood
  for(i in 1:M){
    z[i] ~ dbinom( size = N, prob = expit(logit_p_mn[i] + logit_p[i]) )
  }
  # process structure
  logit_p_mn[1:M] <- Hmat[1:M,1:Ne] %*% alpha[1:Ne]
  
  logit_p[1:M] <- sigma*logit_p_scl[1:M]
  for(i in 1:M){
    logit_p_scl[i] ~ dskewt( mu = 0, sigma = 1, delta = delta, nu = nu )
  }
  
  alpha[1] <- beta[1]
  alpha[2:Ne] <- tau*beta[2:Ne]
  # hyperparameters
  beta[1] ~ dnorm( 0, sd = 10 )
  for(j in 2:Ne){
    beta[j] ~ dnorm( 0, sd = 1 )
  }
  tau ~ dunif( 0, 100 )
  sigma ~ dunif( 0, 100 )
  delta ~ dunif( -1, 1 )
  nu ~ dunif( 0, 1 )
  
})

nim_model <- nimbleModel(
  code = nim_code, 
  constants = list( M = M, Ne = 30 + 1 ),
  inits = list( beta = rep(0, 30 + 1), sigma = 1, delta = 0, 
                nu = 1/100, logit_p_scl = rep(0,M), tau = 1 ),
  data = list( z = rbinom( M, 200, 0.1 ), N = 200, Hmat = Hmat )
)

compl_model <- compileNimble( nim_model )
conf_model <- configureMCMC( nim_model )
conf_model$addMonitors(c("logit_p", "logit_p_mn"))
nim_mcmc <- buildMCMC(conf_model)
nim_Cmcmc <- compileNimble(nim_mcmc, project = nim_model)
