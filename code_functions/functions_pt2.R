#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: functions
# All behind-the-scenes functions.
#===============================================================

#===============================================================
# Part 2: Calculate the EOFs
#===============================================================

M <- 237
N.yr <- 56
xyz.dist <- wraf05_R3_dist/max(wraf05_R3_dist) # scaled R^3 distances

# Gaussian RE
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.19) )
sd.A <- c( 0.72, 0.74, 0.775 )
sd.N <- c( 0.72, 0.74, 0.775 )
schm.lb <- c(0.80, 0.45, 0.10)
schm.ub <- c(0.90, 0.55, 0.20)

gauss.HmatA.array <- gauss.HmatN.array <- array(NA, dim=c(M,M,3))
temp.pA <- temp.pN <- matrix(NA,M,N.yr)
set.seed(0)
for(schm in 1:3){
  for(i in 1:N.yr){
    pA <- 2
    pN <- 1
    while( mean(pA/pN > 1) < schm.lb[schm] || mean(pA/pN > 1) > schm.ub[schm] ){
      beta.A <- rnorm( M, mean = 0, sd = sd.A[schm] )
      beta.N <- rnorm( M, mean = 0, sd = sd.N[schm] )
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
    }
    temp.pA[,i] <- logit(pA)
    temp.pN[,i] <- logit(pN)
  }
  A.eig <- eigen( cov(t(temp.pA)) )
  N.eig <- eigen( cov(t(temp.pN)) )
  
  gauss.HmatA.array[,,schm] <- A.eig$vectors
  gauss.HmatN.array[,,schm] <- N.eig$vectors
}

# Gamma RE
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.18) )
shape.A <- c( 4, 3.75, 3.5 )
shape.N <- c( 4, 3.75, 3.5 )
scale.A <- c( 1.5/4, 0.4, 1.5/3.5 )
scale.N <- c( 1.5/4, 0.4, 1.5/3.5 )
shift.A <- c( 1.5, 1.5, 1.5 )
shift.N <- c( 1.5, 1.5, 1.5 )

gamma.HmatA.array <- gamma.HmatN.array <- array(NA, dim=c(M,M,3))
temp.pA <- temp.pN <- matrix(NA,M,N.yr)
set.seed(0)
for(schm in 1:3){
  for(i in 1:N.yr){
    pA <- 2
    pN <- 1
    while( mean(pA/pN > 1) < schm.lb[schm] || mean(pA/pN > 1) > schm.ub[schm] ){
      beta.A <- rgamma(M, shape = shape.A[schm], scale = scale.A[schm]) - shift.A[schm]
      beta.N <- rgamma(M, shape = shape.N[schm], scale = scale.N[schm]) - shift.N[schm]
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
    }
    temp.pA[,i] <- logit(pA)
    temp.pN[,i] <- logit(pN)
  }
  A.eig <- eigen( cov(t(temp.pA)) )
  N.eig <- eigen( cov(t(temp.pN)) )
  
  gamma.HmatA.array[,,schm] <- A.eig$vectors
  gamma.HmatN.array[,,schm] <- N.eig$vectors
}

# GP-S 
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

mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.18) )
var.A <- c( 0.6, 0.6, 0.6 )
var.N <- c( 0.6, 0.6, 0.6 )
range.A <- c( 0.06, 0.06, 0.06 )
range.N <- c( 0.06, 0.06, 0.06 )
smooth.A <- c( 2, 2, 2 )
smooth.N <- c( 2, 2, 2 )

gps.HmatA.array <- gps.HmatN.array <- array(NA, dim=c(M,M,3))
temp.pA <- temp.pN <- matrix(NA,M,N.yr)
set.seed(0)
for(schm in 1:3){
  for(i in 1:N.yr){
    pA <- 2
    pN <- 1
    while( mean(pA/pN > 1) < schm.lb[schm] || mean(pA/pN > 1) > schm.ub[schm] ){
      beta.A <- draw_gp( distances = xyz.dist, cov.model = "matern", 
                         params = c(var.A[schm], range.A[schm], smooth.A[schm]) )
      beta.A <- beta.A - mean(beta.A)
      beta.N <- draw_gp( distances = xyz.dist, cov.model = "matern", 
                         params = c(var.N[schm], range.N[schm], smooth.N[schm]) )
      beta.N <- beta.N - mean(beta.N)
      
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
    }
    temp.pA[,i] <- logit(pA)
    temp.pN[,i] <- logit(pN)
  }
  A.eig <- eigen( cov(t(temp.pA)) )
  N.eig <- eigen( cov(t(temp.pN)) )
  
  gps.HmatA.array[,,schm] <- A.eig$vectors
  gps.HmatN.array[,,schm] <- N.eig$vectors
}

# GP-L
mu.A <- c( logit(0.08), logit(0.08), logit(0.08) )
mu.N <- c( logit(0.03), logit(0.08), logit(0.18) )
var.A <- c( 0.6, 0.6, 0.6 )
var.N <- c( 0.6, 0.6, 0.6 )
range.A <- c( 0.1, 0.1, 0.1 )
range.N <- c( 0.1, 0.1, 0.1 )
smooth.A <- c( 2, 2, 2 )
smooth.N <- c( 2, 2, 2 )

gpl.HmatA.array <- gpl.HmatN.array <- array(NA, dim=c(M,M,3))
temp.pA <- temp.pN <- matrix(NA,M,N.yr)
set.seed(0)
for(schm in 1:3){
  for(i in 1:N.yr){
    pA <- 2
    pN <- 1
    while( mean(pA/pN > 1) < schm.lb[schm] || mean(pA/pN > 1) > schm.ub[schm] ){
      beta.A <- draw_gp( distances = xyz.dist, cov.model = "matern", 
                         params = c(var.A[schm], range.A[schm], smooth.A[schm]) )
      beta.A <- beta.A - mean(beta.A)
      beta.N <- draw_gp( distances = xyz.dist, cov.model = "matern", 
                         params = c(var.N[schm], range.N[schm], smooth.N[schm]) )
      beta.N <- beta.N - mean(beta.N)
      
      pA <- expit( mu.A[schm] + beta.A )
      pN <- expit( mu.N[schm] + beta.N ) 
    }
    temp.pA[,i] <- logit(pA)
    temp.pN[,i] <- logit(pN)
  }
  A.eig <- eigen( cov(t(temp.pA)) )
  N.eig <- eigen( cov(t(temp.pN)) )
  
  gpl.HmatA.array[,,schm] <- A.eig$vectors
  gpl.HmatN.array[,,schm] <- N.eig$vectors
}

# EOF-NG
eofng.HmatA.array <- eofng.HmatN.array <- array(NA, dim=c(M,N.yr,3))
for(t in 1:3){
  eofng.HmatA.array[,,t] <- ALL_hotJan_EOFs
  eofng.HmatN.array[,,t] <- NAT_hotJan_EOFs
}

