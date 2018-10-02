#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: case studies
#===============================================================

# Constants
hot.threshold <- 5
wet.threshold <- 1
M <- nrow(casestudy_df)
R1.alpha <- 0.1
HmatA.Jan <- ALL_hotJan_EOFs
HmatN.Jan <- NAT_hotJan_EOFs
HmatA.Mar <- ALL_wetMarch_EOFs
HmatN.Mar <- NAT_wetMarch_EOFs

# Source in the nimble code
source("code_functions/nimble_RNB.R")

#===============================================================
# Fit to hot, wet
#===============================================================

# hot, ALL
compl_model$z <- casestudy_df$zA_hotJan2015
compl_model$N <- casestudy_df$nA_2015[1]
compl_model$Hmat <- cbind(rep(1,M),HmatA.Jan[,1:N.yr])
set.seed(0)
nim_Cmcmc$run(20000)
postsamp.hot.ALL <- as.matrix(nim_Cmcmc$mvSamples)[10001:20000,]

# hot, NAT
compl_model$z <- casestudy_df$zN_hotJan2015
compl_model$N <- casestudy_df$nN_2015[1]
compl_model$Hmat <- cbind(rep(1,M),HmatN.Jan[,1:N.yr])
set.seed(0)
nim_Cmcmc$run(20000)
postsamp.hot.NAT <- as.matrix(nim_Cmcmc$mvSamples)[10001:20000,]

# wet, ALL
compl_model$z <- casestudy_df$zA_wetMarch2015
compl_model$N <- casestudy_df$nA_2015[1]
compl_model$Hmat <- cbind(rep(1,M),HmatA.Mar[,1:N.yr])
set.seed(0)
nim_Cmcmc$run(20000)
postsamp.wet.ALL <- as.matrix(nim_Cmcmc$mvSamples)[10001:20000,]

# wet, NAT
compl_model$z <- casestudy_df$zN_wetMarch2015
compl_model$N <- casestudy_df$nN_2015[1]
compl_model$Hmat <- cbind(rep(1,M),HmatN.Mar[,1:N.yr])
set.seed(0)
nim_Cmcmc$run(20000)
postsamp.wet.NAT <- as.matrix(nim_Cmcmc$mvSamples)[10001:20000,]

# save(postsamp.hot.ALL, postsamp.hot.NAT, postsamp.wet.ALL, postsamp.wet.NAT,
#      file = "data/postsamp_RNB.RData")

#===============================================================
# Hot: Calculate RR and classify
#===============================================================
pA.samp.hot <- expit(postsamp.hot.ALL[,c(N.yr+1+3):c(N.yr+1+3+M-1)] + postsamp.hot.ALL[,c(N.yr+1+3+M):c(N.yr+1+3+2*M-1)])
pN.samp.hot <- expit(postsamp.hot.NAT[,c(N.yr+1+3):c(N.yr+1+3+M-1)] + postsamp.hot.NAT[,c(N.yr+1+3+M):c(N.yr+1+3+2*M-1)])
RR.samp.hot <- pA.samp.hot/pN.samp.hot

# H0: RR <= 5
classify.hot <- classify_R1( RR.samp = RR.samp.hot, threshold = 5,
                             region.names = casestudy_df$region, R1.alpha = 0.1 )

hot.results.df <- data.frame( region.names = casestudy_df$region, H0.leq5 = classify.hot$ind.reject )

#===============================================================
# Hot: Calculate RR and classify using BH/LRT
#===============================================================
hot.Pval.vec <- rep(NA, M)
for(i in 1:M){
  hot.Pval.vec[i] <- lrt.pval( zA = casestudy_df$zA_hotJan2015[i],
                               nA = casestudy_df$nA_2015[1],
                               zN = casestudy_df$zN_hotJan2015[i],
                               nN = casestudy_df$nN_2015[1], RR0 = 5 )$pval
}
hot.BH.results.leq1 <- classify_BH( Pvals = hot.Pval.vec, alpha = 1/10, region.names = casestudy_df$region )

#===============================================================
# Wet: Calculate RR and classify
#===============================================================
pA.samp.wet <- expit(postsamp.wet.ALL[,c(N.yr+1+3):c(N.yr+1+3+M-1)] + postsamp.wet.ALL[,c(N.yr+1+3+M):c(N.yr+1+3+2*M-1)])
pN.samp.wet <- expit(postsamp.wet.NAT[,c(N.yr+1+3):c(N.yr+1+3+M-1)] + postsamp.wet.NAT[,c(N.yr+1+3+M):c(N.yr+1+3+2*M-1)])
RR.samp.wet <- pA.samp.wet/pN.samp.wet

# H0: RR <= 1
classify.wet.leq1 <- classify_R1( RR.samp = RR.samp.wet, threshold = 1,
                                  region.names = casestudy_df$region, R1.alpha = 0.1 )

# H0: RR <= 2
classify.wet.leq2 <- classify_R1( RR.samp = RR.samp.wet, threshold = 2,
                                  region.names = casestudy_df$region, R1.alpha = 0.1 )

# H0: RR <= 2/3 or RR >= 3/2
classify.wet.out <- classify_R1( RR.samp = RR.samp.wet, threshold = 2/3,
                                 threshold2 = 3/2, direction = "out",
                                 region.names = casestudy_df$region, R1.alpha = 0.1 )

# H0: RR <= 1/2 or RR >= 2
classify.wet.out2 <- classify_R1( RR.samp = RR.samp.wet, threshold = 1/2,
                                  threshold2 = 2, direction = "out",
                                  region.names = casestudy_df$region, R1.alpha = 0.1 )

# H0: RR >= 1
classify.wet.geq1 <- classify_R1( RR.samp = RR.samp.wet, threshold = 1,
                                  region.names = casestudy_df$region,  direction = "geq", R1.alpha = 0.1 )

# H0: RR >= 1/2
classify.wet.geq05 <- classify_R1( RR.samp = RR.samp.wet, threshold = 1/2,
                                   region.names = casestudy_df$region, direction = "geq", R1.alpha = 0.1 )

#===============================================================
# Wet: Calculate RR and classify using BH/LRT
#===============================================================
wet.Pval.vec <- rep(NA, M)
for(i in 1:M){
  wet.Pval.vec[i] <- lrt.pval( zA = casestudy_df$zA_wetMarch2015[i],
                               nA = casestudy_df$nA_2015[1],
                               zN = casestudy_df$zN_wetMarch2015[i],
                               nN = casestudy_df$nN_2015[1], RR0 = 1 )$pval
}
wet.BH.results.leq1 <- classify_BH( Pvals = wet.Pval.vec, alpha = 1/10, region.names = casestudy_df$region )

#===============================================================
# Wet: combine results
#===============================================================

wet.results.df <- data.frame(
  region.names = casestudy_df$region,
  H0.leq1 = classify.wet.leq1$ind.reject,
  H0.leq2 = classify.wet.leq2$ind.reject,
  H0.out = classify.wet.out$ind.reject,
  H0.out2 = classify.wet.out2$ind.reject,
  H0.geq1 = classify.wet.geq1$ind.reject,
  H0.geq05 = classify.wet.geq05$ind.reject
)
