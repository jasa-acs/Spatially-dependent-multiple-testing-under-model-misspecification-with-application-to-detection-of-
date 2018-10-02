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

#===============================================================
# M3: skew-t random effects
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

# Model, setup, compile
nim_code <- nimbleCode({
  # likelihood
  for(i in 1:M){
    z[i] ~ dbinom( size = N, prob = expit(logit_p[i]) )
    logit_p[i] ~ dskewt( mu = mu, sigma = sigma, delta = delta, nu = nu )
  }
  
  # hyperparameters
  mu ~ dnorm( 0, sd = 10 )
  sigma ~ dunif( 0, 100 )
  delta ~ dunif( -1, 1 )
  nu ~ dunif( 0, 1 )
  
})

nim_model <- nimbleModel(
  code = nim_code, 
  constants = list( M = M ),
  inits = list( mu = 0, sigma = 1, delta = 0, nu = 1/100, logit_p = rep(0,M) ),
  data = list( z = rbinom( M, 200, 0.1 ), N = 200 )
)

compl_model <- compileNimble( nim_model )
conf_model <- configureMCMC( nim_model )

# Look at configuration
conf_model$addMonitors('logit_p')

# Build/compile MCMC
nim_mcmc <- buildMCMC(conf_model)
nim_Cmcmc <- compileNimble(nim_mcmc, project = nim_model)

