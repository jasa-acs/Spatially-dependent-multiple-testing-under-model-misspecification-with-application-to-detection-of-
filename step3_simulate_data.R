#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: simulate data
# Generate all simulated data sets for Section 3
#===============================================================

logit <- function(x){return(log(x/(1-x)))}
expit <- function(x){return(exp(x)/(1+exp(x)))}

# Constants
threshold <- 1
M <- nrow(wraf05_centroids)
N.ens <- c(50, 100, 400)
N.rep <- 100
# Bounds on the true proportion of false H0
schm.lb <- c(0.80, 0.45, 0.10)
schm.ub <- c(0.90, 0.55, 0.20)

#===============================================================
# TRUE STATE 1 (G-RE): Gaussian random effects
#===============================================================
# Storage
gaussian.pA.save <- gaussian.pN.save <- array(NA, dim=c(M,N.rep,3))
gaussian.zA.save <- gaussian.zN.save <- list()

# Hyperparameters
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.19) )
sd.A <- c( 0.72, 0.74, 0.775 )
sd.N <- c( 0.72, 0.74, 0.775 )

# Simulate pA and pN
set.seed(1)
for(schm in 1:3){
  for(t in 1:N.rep){
    
    propFalesH0 <- 0
    while( propFalesH0 > schm.ub[schm] | propFalesH0 < schm.lb[schm] ){
      beta.A <- rnorm( M, mean = 0, sd = sd.A[schm] )
      beta.N <- rnorm( M, mean = 0, sd = sd.N[schm] )
      
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
      
      propFalesH0 <- mean( pA/pN > 1 )
    }
    gaussian.pA.save[,t,schm] <- pA
    gaussian.pN.save[,t,schm] <- pN
  }
}
# Draw Binomial samples
for(scheme in 1:3){
  gaussian.zA.save[[scheme]] <- gaussian.zN.save[[scheme]] <- list()
  for(ens.size in 1:length(N.ens)){
    gaussian.zA.save[[scheme]][[ens.size]] <- gaussian.zN.save[[scheme]][[ens.size]] <- matrix(NA,M,N.rep)
    for(rep in 1:N.rep){
      gaussian.zA.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gaussian.pA.save[,rep,scheme])
      gaussian.zN.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gaussian.pN.save[,rep,scheme])
    }
  }
}

#===============================================================
# TRUE STATE 2 (NG-RE): Non-Gaussian random effects
#===============================================================
# Storage
gamma.pA.save <- gamma.pN.save <- array(NA, dim=c(M,N.rep,3))
gamma.zA.save <- gamma.zN.save <- list()

# Hyperparameters
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.18) )
shape.A <- c( 4, 3.75, 3.5 )
shape.N <- c( 4, 3.75, 3.5 )
scale.A <- c( 1.5/4, 0.4, 1.5/3.5 )
scale.N <- c( 1.5/4, 0.4, 1.5/3.5 )
shift.A <- c( 1.5, 1.5, 1.5 )
shift.N <- c( 1.5, 1.5, 1.5 )

# Simulate pA and pN
set.seed(2)
for(schm in 1:3){
  for(t in 1:N.rep){
    
    propFalesH0 <- 0
    while( propFalesH0 > schm.ub[schm] | propFalesH0 < schm.lb[schm] ){
      beta.A <- rgamma(M, shape = shape.A[schm], scale = scale.A[schm]) - shift.A[schm]
      beta.N <- rgamma(M, shape = shape.N[schm], scale = scale.N[schm]) - shift.N[schm]
      
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
      
      propFalesH0 <- mean( pA/pN > 1 )
    }
    gamma.pA.save[,t,schm] <- pA
    gamma.pN.save[,t,schm] <- pN
  }
}
# Draw Binomial samples
for(scheme in 1:3){
  gamma.zA.save[[scheme]] <- gamma.zN.save[[scheme]] <- list()
  for(ens.size in 1:length(N.ens)){
    gamma.zA.save[[scheme]][[ens.size]] <- gamma.zN.save[[scheme]][[ens.size]] <- matrix(NA,M,N.rep)
    for(rep in 1:N.rep){
      gamma.zA.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gamma.pA.save[,rep,scheme])
      gamma.zN.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gamma.pN.save[,rep,scheme])
    }
  }
}

#===============================================================
# TRUE STATE 3 (GP-S): GP effects with short range
#===============================================================
# Storage
gpshort.pA.save <- gpshort.pN.save <- array(NA, dim=c(M,N.rep,3))
gpshort.zA.save <- gpshort.zN.save <- list()

draw_gp <- function( distances, cov.model, params ){
  M <- nrow(distances)
  # Parameters
  tausq <- params[1] # Variance
  phi <- params[2] # Range
  kappa <- params[3] # Smoothness
  Cov.true <- tausq*cov.spatial( obj = distances, cov.model = cov.model,
                                 kappa = kappa, cov.pars = c(1, phi) ) 
  return(mvrnorm(1, mu = rep(0,M), Sigma = Cov.true))
}

# Rescale to lie in [0,1]
distance.use <- wraf05_R3_dist/max(wraf05_R3_dist)

# Hyperparameters
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.18) )
var.A <- c( 0.6, 0.6, 0.6 )
var.N <- c( 0.6, 0.6, 0.6 )
range.A <- c( 0.06, 0.06, 0.06 )
range.N <- c( 0.06, 0.06, 0.06 )
smooth.A <- c( 2, 2, 2 )
smooth.N <- c( 2, 2, 2 )

# Simulate pA and pN
set.seed(3)
for(schm in 1:3){
  for(t in 1:N.rep){
    
    propFalesH0 <- 0
    while( propFalesH0 > schm.ub[schm] | propFalesH0 < schm.lb[schm] ){
      beta.A <- draw_gp( distances = distance.use, cov.model = "matern", 
                         params = c(var.A[schm], range.A[schm], smooth.A[schm]) )
      beta.A <- beta.A - mean(beta.A)
      beta.N <- draw_gp( distances = distance.use, cov.model = "matern", 
                         params = c(var.N[schm], range.N[schm], smooth.N[schm]) )
      beta.N <- beta.N - mean(beta.N)
      
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
      
      propFalesH0 <- mean( pA/pN > 1 )
    }
    gpshort.pA.save[,t,schm] <- pA
    gpshort.pN.save[,t,schm] <- pN
  }
}
# Draw Binomial samples
for(scheme in 1:3){
  gpshort.zA.save[[scheme]] <- gpshort.zN.save[[scheme]] <- list()
  for(ens.size in 1:length(N.ens)){
    gpshort.zA.save[[scheme]][[ens.size]] <- gpshort.zN.save[[scheme]][[ens.size]] <- matrix(NA,M,N.rep)
    for(rep in 1:N.rep){
      gpshort.zA.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gpshort.pA.save[,rep,scheme])
      gpshort.zN.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gpshort.pN.save[,rep,scheme])
    }
  }
}

#===============================================================
# TRUE STATE 4 (GP-L): GP effects with long range
#===============================================================
# Storage
gplong.pA.save <- gplong.pN.save <- array(NA, dim=c(M,N.rep,3))
gplong.zA.save <- gplong.zN.save <- list()

draw_gp <- function( distances, cov.model, params ){
  M <- nrow(distances)
  # Parameters
  tausq <- params[1] # Variance
  phi <- params[2] # Range
  kappa <- params[3] # Smoothness
  Cov.true <- tausq*cov.spatial( obj = distances, cov.model = cov.model,
                                 kappa = kappa, cov.pars = c(1, phi) ) 
  return(mvrnorm(1, mu = rep(0,M), Sigma = Cov.true))
}

# Rescale to lie in [0,1]
distance.use <- wraf05_R3_dist/max(wraf05_R3_dist)

# Hyperparameters
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.18) )
var.A <- c( 0.6, 0.6, 0.6 )
var.N <- c( 0.6, 0.6, 0.6 )
range.A <- c( 0.1, 0.1, 0.1 )
range.N <- c( 0.1, 0.1, 0.1 )
smooth.A <- c( 2, 2, 2 )
smooth.N <- c( 2, 2, 2 )

# Simulate pA and pN
set.seed(4)
for(schm in 1:3){
  for(t in 1:N.rep){
    
    propFalesH0 <- 0
    while( propFalesH0 > schm.ub[schm] | propFalesH0 < schm.lb[schm] ){
      beta.A <- draw_gp( distances = distance.use, cov.model = "matern", 
                         params = c(var.A[schm], range.A[schm], smooth.A[schm]) )
      beta.A <- beta.A - mean(beta.A)
      beta.N <- draw_gp( distances = distance.use, cov.model = "matern", 
                         params = c(var.N[schm], range.N[schm], smooth.N[schm]) )
      beta.N <- beta.N - mean(beta.N)
      
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
      
      propFalesH0 <- mean( pA/pN > 1 )
    }
    gplong.pA.save[,t,schm] <- pA
    gplong.pN.save[,t,schm] <- pN
  }
}
# Draw Binomial samples
for(scheme in 1:3){
  gplong.zA.save[[scheme]] <- gplong.zN.save[[scheme]] <- list()
  for(ens.size in 1:length(N.ens)){
    gplong.zA.save[[scheme]][[ens.size]] <- gplong.zN.save[[scheme]][[ens.size]] <- matrix(NA,M,N.rep)
    for(rep in 1:N.rep){
      gplong.zA.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gplong.pA.save[,rep,scheme])
      gplong.zN.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], gplong.pN.save[,rep,scheme])
    }
  }
}

#===============================================================
# TRUE STATE 5 (EOF-G): EOF mean with Gaussian discrepancy
#===============================================================
# Storage
eofG.pA.save <- eofG.pN.save <- array(NA, dim=c(M,N.rep,3))
eofG.zA.save <- eofG.zN.save <- list()

eofALL <- ALL_hotJan_EOFs
eofNAT <- NAT_hotJan_EOFs

N.eof <- 30
Hmat.A <- eofALL[,1:N.eof]
Hmat.N <- eofNAT[,1:N.eof]

# Hyperparameters
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.19) )
sd.A <- c( 0.01, 0.01, 0.01 )
sd.N <- c( 0.01, 0.01, 0.01 )

coef.sd <- rep(c(3.5,1,0.05), c(5,5,20))

# Simulate pA and pN
set.seed(5)
for(schm in 1:3){
  for(t in 1:N.rep){
    
    propFalesH0 <- 0
    while( propFalesH0 > schm.ub[schm] | propFalesH0 < schm.lb[schm] ){
      alpha.A <- alpha.N <- rep(NA, N.eof)
      for(j in 1:N.eof){
        alpha.A[j] <- rnorm(1, mean = 0, sd = coef.sd[j])
        alpha.N[j] <- rnorm(1, mean = 0, sd = coef.sd[j])
      }
      
      beta.A <- rnorm( M, mean = 0, sd = sd.A[schm] )
      beta.N <- rnorm( M, mean = 0, sd = sd.N[schm] )
      
      pA <- expit( mu.A[schm] + Hmat.A %*% alpha.A + beta.A )
      pN <- expit( mu.N[schm] + Hmat.N %*% alpha.N + beta.N ) 
      
      propFalesH0 <- mean( pA/pN > 1 )
    }
    eofG.pA.save[,t,schm] <- pA
    eofG.pN.save[,t,schm] <- pN
  }
}

# Draw Binomial samples
for(scheme in 1:3){
  eofG.zA.save[[scheme]] <- eofG.zN.save[[scheme]] <- list()
  for(ens.size in 1:length(N.ens)){
    eofG.zA.save[[scheme]][[ens.size]] <- eofG.zN.save[[scheme]][[ens.size]] <- matrix(NA,M,N.rep)
    for(rep in 1:N.rep){
      eofG.zA.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], eofG.pA.save[,rep,scheme])
      eofG.zN.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], eofG.pN.save[,rep,scheme])
    }
  }
}

#===============================================================
# TRUE STATE 6 (EOF-NG): EOF mean with gamma discrepancy
#===============================================================
# Storage
eofNG.pA.save <- eofNG.pN.save <- array(NA, dim=c(M,N.rep,3))
eofNG.zA.save <- eofNG.zN.save <- list()

eofALL <- ALL_hotJan_EOFs
eofNAT <- NAT_hotJan_EOFs

N.eof <- 30
Hmat.A <- eofALL[,1:N.eof]
Hmat.N <- eofNAT[,1:N.eof]

# Hyperparameters
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.19) )

shape.A <- c( 5, 5, 5 )
shape.N <- c( 5, 5, 5 )
scale.A <- c( 2/5, 2/5, 2/5 )
scale.N <- c( 2/5, 2/5, 2/5 )
shift.A <- c( 2, 2, 2 )
shift.N <- c( 2, 2, 2 )

coef.sd <- rep(c(3.5,1,0.05), c(5,5,20))

# Simulate pA and pN
set.seed(6)
for(schm in 1:3){
  for(t in 1:N.rep){
    
    propFalesH0 <- 0
    while( propFalesH0 > schm.ub[schm] | propFalesH0 < schm.lb[schm] ){
      alpha.A <- alpha.N <- rep(NA, N.eof)
      for(j in 1:N.eof){
        alpha.A[j] <- rnorm(1, mean = 0, sd = coef.sd[j])
        alpha.N[j] <- rnorm(1, mean = 0, sd = coef.sd[j])
      }
      
      beta.A <- 0.02*(rgamma(M, shape = shape.A[schm], scale = scale.A[schm]) - shift.A[schm])
      beta.N <- 0.02*(rgamma(M, shape = shape.N[schm], scale = scale.N[schm]) - shift.N[schm])
      
      pA <- expit( mu.A[schm] + Hmat.A %*% alpha.A + beta.A )
      pN <- expit( mu.N[schm] + Hmat.N %*% alpha.N + beta.N ) 
      
      propFalesH0 <- mean( pA/pN > 1 )
    }
    eofNG.pA.save[,t,schm] <- pA
    eofNG.pN.save[,t,schm] <- pN
  }
}

# Draw Binomial samples
for(scheme in 1:3){
  eofNG.zA.save[[scheme]] <- eofNG.zN.save[[scheme]] <- list()
  for(ens.size in 1:length(N.ens)){
    eofNG.zA.save[[scheme]][[ens.size]] <- eofNG.zN.save[[scheme]][[ens.size]] <- matrix(NA,M,N.rep)
    for(rep in 1:N.rep){
      eofNG.zA.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], eofNG.pA.save[,rep,scheme])
      eofNG.zN.save[[scheme]][[ens.size]][,rep] <- rbinom(M, N.ens[ens.size], eofNG.pN.save[,rep,scheme])
    }
  }
}
