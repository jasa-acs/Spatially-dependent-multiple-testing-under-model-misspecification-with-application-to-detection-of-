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
# Part 1: code to conduct the multiple testing procedures
#===============================================================

# LRT for a traditional hypothesis test
logLik <- function(yA, yN, nA, nN, pA, pN) {
  dbinom(yA, nA, pA, log = TRUE) + dbinom(yN, nN, pN, log = TRUE)
}
restricted_mle <- function(pAhat, pNhat, rr) {
  # from Farrington and Manning, Stats in Med, 9:1447 (1990)
  # easily derivable from likelihood with simple plug-in constraint
  a <- 2
  b <- -(rr * (1 + pNhat) + 1 + pAhat)
  c <- rr*(pAhat + pNhat)
  pAtilde <- (-b - sqrt(b*b - 4*a*c)) / (2*a)
  return(pAtilde)
}
lrt.pval <- function( zA, nA, zN, nN, RR0 ){
  
  # Restricted log-likelihood
  if( zN/nN < zA/nA/RR0 ){
    pA.restr <- restricted_mle(zA/nA,zN/nN,RR0)
    restr.loglik <- logLik(zA,zN,nA,nN,pA.restr,pA.restr/RR0)
  }
  else{
    restr.loglik <- logLik(zA,zN,nA,nN,zA/nA,zN/nN)
  }
  
  # Unrestricted log-likelihood
  unrestr.loglik <- logLik(zA,zN,nA,nN,zA/nA,zN/nN)
  
  # Test
  lambda <- -2*(restr.loglik - unrestr.loglik)
  
  # P-val 
  pval <- pchisq(lambda, 1, lower.tail = FALSE)
  
  return(list(restr.loglik = restr.loglik,
              unrestr.loglik = unrestr.loglik,
              lambda = lambda, pval = pval ))
  
}

# Classify based on Benjamini-Hochberg FDR procedure
classify_BH <- function( Pvals, alpha = 0.05, region.names ){
  
  K <- length(Pvals)
  # Collect into a data frame
  test.df <- data.frame( position = 1:K, region = region.names, Pvals = Pvals )
  # Sort
  test.df <- test.df[order(test.df$Pvals),]
  
  # Identify the cutoff position for prespecified alpha
  if( min(test.df$Pvals) > alpha ){ ind.temp <- rep(0,K) }
  else{
    ind.temp <- rep(NA,K)
    for(j in 1:K){
      ifelse( test.df$Pvals[j] <= (j*alpha)/K, ind.temp[j] <- 1, ind.temp[j] <- 0)
    }
  }
  test.df$ind.reject <- ind.temp
  
  # Resort 
  classify.H0.df <- test.df[order(test.df$position),]
  return(classify.H0.df)
}

# Classify based on R1 criteria
classify_R1 <- function( RR.samp, direction = "leq", threshold = 1, threshold2 = NULL,
                         region.names = NULL, R1.alpha = 0.1 ){
  
  K <- ncol(RR.samp)  # number of clusters
  B <- nrow(RR.samp)  # number of MCMC samples
  
  if( is.null(region.names) ){ region.names <- 1:K }
  
  # Calculate the posterior probs of H0
  if( direction == "leq" ){
    post.prob.H0.df <- data.frame(
      position = 1:K,
      region = region.names,
      post.prob.H0 = colMeans( RR.samp <= threshold )
    )
  }
  if( direction == "geq" ){
    post.prob.H0.df <- data.frame(
      position = 1:K,
      region = region.names,
      post.prob.H0 = colMeans( RR.samp >= threshold )
    )
  }
  if( direction == "both" ){
    post.prob.H0.df <- data.frame(
      position = 1:K,
      region = region.names,
      post.prob.H0 = colMeans( RR.samp >= threshold & RR.samp <= threshold2 )
    )
  }
  if( direction == "out" ){
    post.prob.H0.df <- data.frame(
      position = 1:K,
      region = region.names,
      post.prob.H0 = 1 - colMeans( RR.samp >= threshold & RR.samp <= threshold2 )
    )
  }
  
  # Sort
  post.prob.H0.df.sorted <- post.prob.H0.df[order(post.prob.H0.df$post.prob.H0),]
  
  # R1 =================
  R1.threshold.stat <- rep(NA,K)
  for(j in 1:K){
    R1.threshold.stat[j] <- mean(post.prob.H0.df.sorted$post.prob.H0[1:j])
  }
  
  # Identify the cutoff position for prespecified alpha
  if( min(R1.threshold.stat) > R1.alpha ){
    R1.ind.reject <- rep(0,K)
  }
  else{
    R1.cutoff.ind <- which.max( R1.threshold.stat[R1.threshold.stat <= R1.alpha] )
    
    # Create a vector indicating which locations should be rejected
    R1.ind.reject <- rep(0,K)
    R1.ind.reject[1:R1.cutoff.ind] <- 1
  }
  post.prob.H0.df.sorted$threshold.stat <- R1.threshold.stat
  post.prob.H0.df.sorted$ind.reject <- R1.ind.reject
  
  # Unsort
  classify.H0.df <- post.prob.H0.df.sorted[order(post.prob.H0.df.sorted$position),]
  return( classify.H0.df )
}

# Classify based on R2 criteria
classify_R2 <- function( RR.samp, threshold = 1, region.names = NULL, R2.lambda = 1/9 ){
  
  K <- ncol(RR.samp)  # number of clusters
  B <- nrow(RR.samp)  # number of MCMC samples
  cutoff <- 1/(1 + R2.lambda)
  
  if( is.null(region.names) ){ region.names <- 1:K }
  
  # Calculate the posterior probs of H0
  post.prob.H0.df <- data.frame(
    position = 1:K,
    region = region.names,
    post.prob.H0 = colMeans( RR.samp < threshold )
  )
  
  post.prob.H0.df$ind.reject <- rep(0,K)
  post.prob.H0.df$ind.reject[post.prob.H0.df$post.prob.H0 < cutoff] <- 1
  
  # Unsort
  classify.H0.df <- post.prob.H0.df
  return( classify.H0.df )
}

# Classify based on R3 criteria
classify_R3 <- function( RR.samp, threshold = 1, region.names = NULL, R3.alpha = 0.1*ncol(RR.samp) ){
  
  K <- ncol(RR.samp)  # number of clusters
  B <- nrow(RR.samp)  # number of MCMC samples
  
  if( is.null(region.names) ){ region.names <- 1:K }
  
  # Calculate the posterior probs of H0
  post.prob.H0.df <- data.frame(
    position = 1:K,
    region = region.names,
    post.prob.H0 = colMeans( RR.samp < threshold )
  )
  
  # Sort
  post.prob.H0.df.sorted <- post.prob.H0.df[order(post.prob.H0.df$post.prob.H0),]
  
  # R3 =================
  R3.threshold.stat <- cumsum(post.prob.H0.df.sorted$post.prob.H0)
  
  # Identify the cutoff position for prespecified alpha
  if( min(R3.threshold.stat) > R3.alpha ){
    R3.ind.reject <- rep(0,K)
  }
  else{
    R3.cutoff.ind <- which.max( R3.threshold.stat[R3.threshold.stat <= R3.alpha] )
    
    # Create a vector indicating which locations should be rejected
    R3.ind.reject <- rep(0,K)
    R3.ind.reject[1:R3.cutoff.ind] <- 1
  }
  post.prob.H0.df.sorted$threshold.stat <- R3.threshold.stat
  post.prob.H0.df.sorted$ind.reject <- R3.ind.reject
  
  # Unsort
  classify.H0.df <- post.prob.H0.df.sorted[order(post.prob.H0.df.sorted$position),]
  return( classify.H0.df )
}

# Evaluate FDs/FDR/loss
evaluate_loss <- function( classify.obj, lambda, trueRR, threshold ){
  
  # No. of false discoveries: rejected H0 but RR <= threshold
  known.FDs <- sum( classify.obj$ind.reject == 1 & trueRR <= threshold )
  # No. of false negatives: failed to reject H0 but RR > threshold
  known.FNs <- sum( classify.obj$ind.reject == 0 & trueRR > threshold )
  # No. of true discoveries: rejected H0 and RR > threshold
  known.TDs <- sum( classify.obj$ind.reject == 1 & trueRR > threshold )
  # No. of true negatives: failed to reject H0 and RR <= threshold
  known.TNs <- sum( classify.obj$ind.reject == 0 & trueRR <= threshold )
  
  R1.loss <- lambda*known.FDs + known.FNs
  return( list( known.FDs = known.FDs, known.FNs = known.FNs,
                known.TDs = known.TDs, known.TNs = known.TNs,
                loss = R1.loss ))
}
evaluate_FDR <- function( classify.obj, trueRR, threshold ){
  
  # No. of false discoveries: rejected H0 but RR <= threshold
  known.FDs <- sum( classify.obj$ind.reject == 1 & trueRR <= threshold )
  # No. of false negatives: failed to reject H0 but RR > threshold
  known.FNs <- sum( classify.obj$ind.reject == 0 & trueRR > threshold )
  # No. of true discoveries: rejected H0 and RR > threshold
  known.TDs <- sum( classify.obj$ind.reject == 1 & trueRR > threshold )
  # No. of true negatives: failed to reject H0 and RR <= threshold
  known.TNs <- sum( classify.obj$ind.reject == 0 & trueRR <= threshold )
  
  ifelse( sum(classify.obj$ind.reject) == 0, FDR <- 0,
          FDR <- sum( classify.obj$ind.reject == 1 & trueRR <= threshold )/sum(classify.obj$ind.reject) )
  ifelse( sum(abs(classify.obj$ind.reject - 1)) == 0, FNR <- 0,
          FNR <- sum( classify.obj$ind.reject == 0 & trueRR > threshold )/(M - sum(classify.obj$ind.reject)) )
  ifelse( sum( trueRR > threshold ) == 0, power <- 1,
          power <- sum( classify.obj$ind.reject == 1 & trueRR > threshold )/sum( trueRR > threshold ) )
  return( list( known.FDs = known.FDs, known.FNs = known.FNs,
                known.TDs = known.TDs, known.TNs = known.TNs,
                FDR = FDR, FNR = FNR, power = power ))
}
evaluate_FD <- function( classify.obj, trueRR, threshold ){
  
  # No. of false discoveries: rejected H0 but RR <= threshold
  known.FDs <- sum( classify.obj$ind.reject == 1 & trueRR <= threshold )
  # No. of false negatives: failed to reject H0 but RR > threshold
  known.FNs <- sum( classify.obj$ind.reject == 0 & trueRR > threshold )
  # No. of true discoveries: rejected H0 and RR > threshold
  known.TDs <- sum( classify.obj$ind.reject == 1 & trueRR > threshold )
  # No. of true negatives: failed to reject H0 and RR <= threshold
  known.TNs <- sum( classify.obj$ind.reject == 0 & trueRR <= threshold )
  
  FD <- known.FDs
  FN <- known.FNs
  return( list( known.FDs = known.FDs, known.FNs = known.FNs,
                known.TDs = known.TDs, known.TNs = known.TNs,
                FD = FD, FN = FN ))
}

# WRAPPER FUNCTION
classify_tests_evaluate <- function(
  zA, zN, N,       # data
  pA, pN,          # true probabilities
  RR.samp,         # n.iter x M matrix of posterior samples of RR
  threshold = 1,   # null exceedance
  R1.alpha = 0.1,  # Control FDR at this level 
  R2.lambda = 9,   # Loss for FD, relative to FN
  R3.alpha = 5,    # Control FD at this level
  fit.LRT = TRUE   # Fit the Frequentist results?
){
  
  M <- ncol(RR.samp)
  
  # LRT fits ====================================
  if(fit.LRT){
    # Calculate P-values
    Pval.vec <- rep(NA, M)
    for(i in 1:M){
      Pval.vec[i] <- lrt.pval( zA = zA[i], nA = N, zN = zN[i], nN = N, RR0 = threshold )$pval
    }
    BH.results <- classify_BH( Pvals = Pval.vec, alpha = 1/10, region.names = 1:M )
    
    FWER.results <- BH.results
    FWER.results$ind.reject <- rep(0,M)
    FWER.results$ind.reject[Pval.vec <= (0.1/M)] <- 1
    
    unadj.results <- BH.results
    unadj.results$ind.reject <- rep(0,M)
    unadj.results$ind.reject[Pval.vec <= (0.1)] <- 1
    
    # Evaluate
    BH.loss <- evaluate_loss( classify.obj = BH.results, lambda = R2.lambda, trueRR = pA/pN, threshold = threshold )
    BH.FDR <- evaluate_FDR( classify.obj = BH.results, trueRR = pA/pN, threshold = threshold )
    BH.FD <- evaluate_FD( classify.obj = BH.results, trueRR = pA/pN, threshold = threshold )
    
    FWER.loss <- evaluate_loss( classify.obj = FWER.results, lambda = R2.lambda, trueRR = pA/pN, threshold = threshold )
    FWER.FDR <- evaluate_FDR( classify.obj = FWER.results, trueRR = pA/pN, threshold = threshold )
    FWER.FD <- evaluate_FD( classify.obj = FWER.results, trueRR = pA/pN, threshold = threshold )
    
    unadj.loss <- evaluate_loss( classify.obj = unadj.results, lambda = R2.lambda, trueRR = pA/pN, threshold = threshold )
    unadj.FDR <- evaluate_FDR( classify.obj = unadj.results, trueRR = pA/pN, threshold = threshold )
    unadj.FD <- evaluate_FD( classify.obj = unadj.results, trueRR = pA/pN, threshold = threshold )
    
  }
  
  # Classify using Bayesian fits ================
  R1.results <- classify_R1( RR.samp = RR.samp, threshold = threshold, R1.alpha = R1.alpha )
  R2.results <- classify_R2( RR.samp = RR.samp, threshold = threshold, R2.lambda = R2.lambda )
  R3.results <- classify_R3( RR.samp = RR.samp, threshold = threshold, R3.alpha = R3.alpha )
  
  # Evaluate
  R1.loss <- evaluate_loss( classify.obj = R1.results, lambda = R2.lambda, trueRR = pA/pN, threshold = threshold )
  R1.FDR <- evaluate_FDR( classify.obj = R1.results, trueRR = pA/pN, threshold = threshold )
  R1.FD <- evaluate_FD( classify.obj = R1.results, trueRR = pA/pN, threshold = threshold )
  
  R2.loss <- evaluate_loss( classify.obj = R2.results, lambda = R2.lambda, trueRR = pA/pN, threshold = threshold )
  R2.FDR <- evaluate_FDR( classify.obj = R2.results, trueRR = pA/pN, threshold = threshold )
  R2.FD <- evaluate_FD( classify.obj = R2.results, trueRR = pA/pN, threshold = threshold )
  
  R3.loss <- evaluate_loss( classify.obj = R3.results, lambda = R2.lambda, trueRR = pA/pN, threshold = threshold )
  R3.FDR <- evaluate_FDR( classify.obj = R3.results, trueRR = pA/pN, threshold = threshold )
  R3.FD <- evaluate_FD( classify.obj = R3.results, trueRR = pA/pN, threshold = threshold )
  
  # Summarize
  if(fit.LRT){
    results <- data.frame( Method = c("unadj", "FWER","BH","R1","R2","R3"),
                           FD = c(unadj.FD$FD, FWER.FD$FD, BH.FD$FD, R1.FD$FD, R2.FD$FD, R3.FD$FD),
                           FN = c(unadj.FD$FN, FWER.FD$FN, BH.FD$FN, R1.FD$FN, R2.FD$FN, R3.FD$FN),
                           TD = c(unadj.FD$known.TDs, FWER.FD$known.TDs, BH.FD$known.TDs, R1.FD$known.TDs, R2.FD$known.TDs, R3.FD$known.TDs),
                           TN = c(unadj.FD$known.TNs, FWER.FD$known.TNs, BH.FD$known.TNs, R1.FD$known.TNs, R2.FD$known.TNs, R3.FD$known.TNs),
                           Loss = c(unadj.loss$loss, FWER.loss$loss, BH.loss$loss, R1.loss$loss, R2.loss$loss, R3.loss$loss),
                           FDR = round(c(unadj.FDR$FDR, FWER.FDR$FDR, BH.FDR$FDR, R1.FDR$FDR, R2.FDR$FDR, R3.FDR$FDR),5),
                           FNR = round(c(unadj.FDR$FNR, FWER.FDR$FNR, BH.FDR$FNR, R1.FDR$FNR, R2.FDR$FNR, R3.FDR$FNR),5),
                           power = round(c(unadj.FDR$power, FWER.FDR$power, BH.FDR$power, R1.FDR$power, R2.FDR$power, R3.FDR$power),5) )
  }
  else{
    results <- data.frame( Method = c("unadj", "FWER","BH","R1","R2","R3"),
                           FD = c(NA, NA, NA, R1.FD$FD, R2.FD$FD, R3.FD$FD),
                           FN = c(NA, NA, NA, R1.FD$FN, R2.FD$FN, R3.FD$FN),
                           TD = c(NA, NA, NA, R1.FD$known.TDs, R2.FD$known.TDs, R3.FD$known.TDs),
                           TN = c(NA, NA, NA, R1.FD$known.TNs, R2.FD$known.TNs, R3.FD$known.TNs),
                           Loss = c(NA, NA, NA, R1.loss$loss, R2.loss$loss, R3.loss$loss),
                           FDR = round(c(NA, NA, NA, R1.FDR$FDR, R2.FDR$FDR, R3.FDR$FDR),5),
                           FNR = round(c(NA, NA, NA, R1.FDR$FNR, R2.FDR$FNR, R3.FDR$FNR),5),
                           power = round(c(NA, NA, NA, R1.FDR$power, R2.FDR$power, R3.FDR$power),5) )
    
  }
  
  return(results)
  
}
