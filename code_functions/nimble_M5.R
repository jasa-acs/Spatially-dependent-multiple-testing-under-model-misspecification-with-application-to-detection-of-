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
# M5: Hybrid CAR random effects
#===============================================================

# Precision for the hybrid CAR model
carhyb_prec <- nimbleFunction(     
  run = function( CAR_prec = double(2), Imat = double(2), lambda = double(0), tau = double(0)) {
    returnType(double(2))
    result <- (1/(tau*tau))*((1-lambda)*Imat + lambda*CAR_prec)
    return(result)
  })

# Create the model 
nim_code <- nimbleCode({
  # likelihood
  for(i in 1:M){
    z[i] ~ dbinom(size = N, prob = expit(logit_p[i]))
  }
  # latent process      
  logit_p[1:M] ~ dmnorm(mean = mn[1:M], prec = C[1:M, 1:M])
  # process structure
  mn[1:M] <- ones[1:M]*mu
  C[1:M,1:M] <- carhyb_prec(CAR_prec[1:M,1:M], Imat[1:M,1:M], lambda, tau)
  # hyperparameters
  mu ~ dnorm( mean = 0, sd = 100 )
  tau ~ dunif(0, 100)
  lambda ~ dunif(0, 1)
})
nim_model <- nimbleModel( 
  code = nim_code, 
  constants = list(M = M, Imat = diag(rep(1,M)),
                   CAR_prec = diag(rowSums(Wmat)) - Wmat, ones = rep(1,M) ),
  inits = list(mu = 0, lambda = 0.5, tau = 1, logit_p = rep(0,M),
               C = carhyb_prec(diag(rowSums(Wmat)) - Wmat, diag(rep(1,M)), 0.5, 1) ),
  data = list(z = rbinom(M, 500, 0.1), N = 500)
)

compl_model <- compileNimble( nim_model ) 
conf_model <- configureMCMC(nim_model ) 
# Look at configuration
conf_model$getSamplers()
conf_model$removeSamplers("logit_p[1:237]")
for(i in 1:floor(M/8)){
  conf_model$addSampler(target = paste("logit_p[",8*(i-1) + 1,":",8*i,"]",sep=""), type = "RW_block")
}
conf_model$addSampler(target = "logit_p[232:237]", type = "RW_block")
conf_model$addMonitors('logit_p')
nim_mcmc <- buildMCMC(conf_model)
nim_Cmcmc <- compileNimble(nim_mcmc, project = nim_model )


