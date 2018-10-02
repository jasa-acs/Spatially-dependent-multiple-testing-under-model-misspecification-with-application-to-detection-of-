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
# M6: Exponential GP random effects
#===============================================================

# Exponential covariance for spatial GP
spatial_expcov <- nimbleFunction(
  run = function(distances = double(2), phi = double(0), tau = double(0)) {
    returnType(double(2))
    result <- tau*tau*exp(-(distances/phi))
    return(result)
  })

# Create the model ========
nim_code <- nimbleCode({
  # likelihood
  for(i in 1:M)
    z[i] ~ dbinom(size = N, prob = expit(logit_p[i]))
  # latent process
  logit_p[1:M] ~ dmnorm(mn[1:M], cov = C[1:M, 1:M])
  # process structure
  mn[1:M] <- ones[1:M]*mu
  C[1:M,1:M] <- spatial_expcov(distances[1:M,1:M], phi, tau)
  # hyperparameters
  mu ~ dnorm(0, sd = 100)
  tau ~ dunif(0, 100)
  phi ~ dunif(0, phi.ub)
})
nim_model <- nimbleModel(
  code = nim_code,
  constants = list(M = M, phi.ub = max(xyz.dist)/2,
                   distances = xyz.dist, ones = rep(1,M) ),
  inits = list(mu = 0, phi = 0.5, tau = 1, logit_p = rep(0,M),
               C = spatial_expcov(xyz.dist, 1, 1) ),
  data = list(z = rbinom(M, 100, 0.1), N = 100)
)

compl_model <- compileNimble(nim_model)
conf_model <- configureMCMC(nim_model)
# Look at configuration
conf_model$removeSamplers("logit_p[1:237]")
for(t in 1:M){
  conf_model$addSampler(target = paste("logit_p[",t,"]",sep=""), type = "RW_block")
}

# conf_model$getSamplers()
conf_model$addMonitors('logit_p')
nim_mcmc <- buildMCMC(conf_model)
nim_Cmcmc <- compileNimble(nim_mcmc, project = nim_model )
