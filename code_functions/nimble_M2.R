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
# M2: Gaussian random effects
#===============================================================
nim_code <- nimbleCode({
  # likelihood
  for(i in 1:M){
    z[i] ~ dbinom(size = N, prob = expit(tau*logit_p[i] + mu))
    logit_p[i] ~ dnorm( mean = 0, sd = 1 )
  }
  # hyperparameters
  mu ~ dnorm(0, sd = 100)
  tau ~ dunif(0, 100)
})

nim_model <- nimbleModel(
  code = nim_code, constants = list(M = M),
  inits = list(mu = 0, tau = 1, logit_p = rep(0,M) ), 
  data = list(z = rbinom(M, 100, 0.1), N = 100)
)
compl_model <- compileNimble(nim_model)
conf_model <- configureMCMC(nim_model)

conf_model$addMonitors("logit_p")
nim_mcmc <- buildMCMC(conf_model)
nim_Cmcmc <- compileNimble(nim_mcmc, project = nim_model )