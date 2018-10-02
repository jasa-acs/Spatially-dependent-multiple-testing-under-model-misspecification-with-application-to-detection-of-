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
# M4: CAR random effects
#===============================================================
Q <- diag(colSums(Wmat)) - Wmat
test <- eigen(Q)
Q.evals <- test$values
Q.evecs <- test$vectors

# MV Gaussian with non-full-rank precision
dmnorm_mdr <- nimbleFunction(
  run = function( x = double(1), mn = double(1), Q = double(2), 
                  Qevecs = double(2), Qevals = double(1),
                  tau = double(0, default = 1),
                  log = integer(0, default = 0)) {
    
    returnType(double(0))
    
    N <- length(Qevals)
    logDet <- -0.5*(N-1)*log(tau*tau) + 0.5*sum(log(Qevals[1:(N-1)]))
    logExp <- -0.5*(1/(tau*tau))*t(x - mn)%*%Q%*%(x - mn)
    if(log){ return( -0.5*(N-1)*log(2*3.141593) + logDet + logExp[1,1] ) }
    else{ return(exp( -0.5*(N-1)*log(2*3.141593) + logDet + logExp[1,1] )) }
  }
)

rmnorm_mdr <- nimbleFunction(
  run = function( n = integer(0, default = 1), mn = double(1), Q = double(2),
                  Qevecs = double(2), Qevals = double(1),
                  tau = double(0, default = 1) ) {
    
    returnType(double(1))
    if(n != 1) print("rmnorm_mdr only allows n = 1; using n = 1.")
    
    N <- length(Qevals)
    z <- rnorm( n = 1, mean = 0, sd = sqrt((tau*tau)/Qevals[1]) )
    y <- z*Qevecs[,1]
    for( i in 2:(N-1)){
      z <- rnorm( n = 1, mean = 0, sd = sqrt((tau*tau)/Qevals[i]) )
      y <- y + z*Qevecs[,i]
    }
    x <- y + mn
    return(x)
  }
)

registerDistributions( list(dmnorm_mdr = list(
  BUGSdist = "dmnorm_mdr(mn, Q, Qevecs, Qevals, tau)",
  types = c("value = double(1)", "mn = double(1)", "Q = double(2)", "Qevecs = double(2)", "Qevals = double(1)" )
)))

# Create the model 
nim_code <- nimbleCode({
  # likelihood
  for(i in 1:M) 
    z[i] ~ dbinom(size = N, prob = expit(logit_p[i]))
  # latent process      
  logit_p[1:M] ~ dmnorm_mdr( mn = mn[1:M], Q = Q[1:M, 1:M], Qevecs = Qevecs[1:M, 1:M], 
                             Qevals = Qevals[1:M], tau = tau )
  # process structure
  mn[1:M] <- ones[1:M]*mu
  # hyperparameters
  tau ~ dunif(0, 100)
})
nim_model <- nimbleModel( 
  code = nim_code, 
  constants = list(mu = 0, M = M, ones = rep(1,M) ),
  inits = list(tau = 1, logit_p = rep(0,M) ),
  data = list(z = rbinom(M,500,0.1), Q = Q, Qevals = Q.evals, Qevecs = Q.evecs, N = 500)
)

compl_model <- compileNimble( nim_model ) 
conf_model <- configureMCMC( nim_model ) 
conf_model$removeSamplers("logit_p[1:237]")
for(i in 1:floor(M/8)){
  conf_model$addSampler(target = paste("logit_p[",8*(i-1) + 1,":",8*i,"]",sep=""), type = "RW_block")
}
conf_model$addSampler(target = "logit_p[232:237]", type = "RW_block")
conf_model$addMonitors('logit_p')
nim_mcmc <- buildMCMC(conf_model)
nim_Cmcmc <- compileNimble(nim_mcmc, project = nim_model ) 
