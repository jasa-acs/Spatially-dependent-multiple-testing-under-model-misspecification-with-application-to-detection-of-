#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: fit the simulation study
#===============================================================

M <- 237 # Sample size
N.yr <- 56 # Number of years to calculate EOFs with
Wmat <- Wmat05 # Adjacency matrix
xyz.dist <- wraf05_R3_dist/max(wraf05_R3_dist) # scaled R^3 distances
ITER <- c(12000, 12000, 12000, 7000, 7000, 7000, 12000, 12000, 12000) # Total MCMC iterations
stSAVE <- ITER - 5000 + 1 # Start post-burn in samples

# Constants
threshold <- 1
N.ens <- c(50, 100, 400)
N.rep <- 100
R1.alpha <- 0.1
R2.lambda <- 9
R3.alpha <- 0.1*M

# Combine
pA.list <- list(gaussian.pA.save, gamma.pA.save, gpshort.pA.save, gplong.pA.save, eofG.pA.save, eofNG.pA.save)
pN.list <- list(gaussian.pN.save, gamma.pN.save, gpshort.pN.save, gplong.pN.save, eofG.pN.save, eofNG.pN.save)
zA.list <- list(gaussian.zA.save, gamma.zA.save, gpshort.zA.save, gplong.zA.save, eofG.zA.save, eofNG.zA.save)
zN.list <- list(gaussian.zN.save, gamma.zN.save, gpshort.zN.save, gplong.zN.save, eofG.zN.save, eofNG.zN.save)

HmatA.use <- list(gauss.HmatA.array, gamma.HmatA.array, gps.HmatA.array, gpl.HmatA.array,
                  eofng.HmatA.array, eofng.HmatA.array)
HmatN.use <- list(gauss.HmatN.array, gamma.HmatN.array, gps.HmatN.array, gpl.HmatN.array,
                  eofng.HmatN.array, eofng.HmatN.array)

# Model names
mod.name <- c("RNB", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9")

# Number of EOFs
N.eofs <- c(N.yr, rep(NA,5), 30, 15, 45)

# Loop over: fitted models and replicates
for(mod_idx in 0:9){
  #==================================
  # Model 1 (and LRT)
  #==================================
  if(mod_idx == 0){

    for(rep_idx in 1:N.rep){
      # Storage
      sv_name <- paste("M1_replicate", rep_idx, "_results.RData", sep = "")
      results_all <- list()

      # Fits for the model/replicate combo (if not already done)
      if( !file.exists(sv_name) ){
        for( true.state in 1:6 ){
          results_all[[true.state]] <- list()

          pA.array <- pA.list[[true.state]]
          pN.array <- pN.list[[true.state]]
          zA.array <- zA.list[[true.state]]
          zN.array <- zN.list[[true.state]]

          # SCHEMES
          for( scheme in 1:3 ){
            results_all[[true.state]][[scheme]] <- list()

            # ENSEMBLE SIZES
            for( ens.size in 1:length(N.ens) ){

              pA <- pA.array[,rep_idx,scheme]
              pN <- pN.array[,rep_idx,scheme]
              zA <- zA.array[[scheme]][[ens.size]][,rep_idx]
              zN <- zN.array[[scheme]][[ens.size]][,rep_idx]

              # Fit to ALL simulations ================
              p.ALL <- matrix(NA, 5000, M)
              for(j in 1:M){
                p.ALL[,j] <- rbeta(n = 5000, shape1 = 1 + zA[j], shape2 = 1 + N.ens[ens.size] - zA[j])
              }

              # Fit to NAT simulations ================
              p.NAT <- matrix(NA, 5000, M)
              for(j in 1:M){
                p.NAT[,j] <- rbeta(n = 5000, shape1 = 1 + zN[j], shape2 = 1 + N.ens[ens.size] - zN[j])
              }

              # Calculate RR ==========================
              RR.samp <- p.ALL/p.NAT

              # Classify hypotheses and evaluate ======
              results_all[[true.state]][[scheme]][[ens.size]] <- classify_tests_evaluate(
                zA = zA, zN = zN, N = N.ens[ens.size], pA = pA, pN = pN, RR.samp = RR.samp, fit.LRT = TRUE,
                threshold = threshold, R1.alpha = R1.alpha, R2.lambda = R2.lambda, R3.alpha = R3.alpha )
            }
          }
        }
        save( results_all, file = sv_name )
      }
      cat(rep_idx, " ")
    }
  } # End Model 1 and LRT

  #==================================
  # RNB and models 2-9
  #==================================
  if(mod_idx > 0){

    for(rep_idx in 1:N.rep){

      # Source the correct nimble code
      source(paste("code_functions/nimble_", mod.name[mod_idx], ".R", sep = ""))

      # Storage
      sv_name <- paste(mod.name[mod_idx], "_replicate", rep_idx, "_results.RData", sep = "")
      results_all <- list()

      # Fits for the model/replicate combo (if not already done)
      if( !file.exists(sv_name) ){
        for( true.state in 1:6 ){
          results_all[[true.state]] <- list()

          pA.array <- pA.list[[true.state]]
          pN.array <- pN.list[[true.state]]
          zA.array <- zA.list[[true.state]]
          zN.array <- zN.list[[true.state]]

          # SCHEMES
          for( scheme in 1:3 ){
            results_all[[true.state]][[scheme]] <- list()

            # ENSEMBLE SIZES
            for( ens.size in 1:length(N.ens) ){

              pA <- pA.array[,rep_idx,scheme]
              pN <- pN.array[,rep_idx,scheme]
              zA <- zA.array[[scheme]][[ens.size]][,rep_idx]
              zN <- zN.array[[scheme]][[ens.size]][,rep_idx]

              # Fit to ALL simulations ================
              compl_model$z <- zA
              compl_model$N <- N.ens[ens.size]
              if( mod_idx %in% c(1,7:9) ){ # Additional assignments for models with EOFs
                compl_model$Hmat <- cbind(rep(1,M),HmatA.use[[true.state]][,1:N.eofs[mod_idx],scheme])
              }
              nim_Cmcmc$run(ITER[mod_idx])
              p.ALL <- as.matrix(nim_Cmcmc$mvSamples)[stSAVE[mod_idx]:ITER[mod_idx],]

              # Fit to NAT simulations ================
              compl_model$z <- zN
              compl_model$N <- N.ens[ens.size]
              if( mod_idx %in% c(1,7:9) ){ # Additional assignments for models with EOFs
                compl_model$Hmat <- cbind(rep(1,M),HmatN.use[[true.state]][,1:N.eofs[mod_idx],scheme])
              }
              nim_Cmcmc$run(ITER[mod_idx])
              p.NAT <- as.matrix(nim_Cmcmc$mvSamples)[stSAVE[mod_idx]:ITER[mod_idx],]

              # Calculate RR ==========================
              if( mod_idx == 1 ){
                pA.samp <- expit(p.ALL[,c(N.yr+1+3):c(N.yr+1+3+M-1)] + p.ALL[,c(N.yr+1+3+M):c(N.yr+1+3+2*M-1)])
                pN.samp <- expit(p.NAT[,c(N.yr+1+3):c(N.yr+1+3+M-1)] + p.NAT[,c(N.yr+1+3+M):c(N.yr+1+3+2*M-1)])
                RR.samp <- pA.samp/pN.samp
              }
              if( mod_idx == 2 ){
                RR.samp <- expit(p.ALL[,"mu"] + p.ALL[,1:M])/expit(p.NAT[,"mu"] + p.NAT[,1:M])
              }
              if( mod_idx %in% c(3,5) ){
                RR.samp <- expit(p.ALL[,2:(M+1)])/expit(p.NAT[,2:(M+1)])
              }
              if( mod_idx %in% c(4,6) ){
                RR.samp <- expit(p.ALL[,1:M])/expit(p.NAT[,1:M])
              }
              if( mod_idx %in% 7:9 ){
                pA.samp <- expit(p.ALL[,c(N.eofs[mod_idx]+1+2):c(N.eofs[mod_idx]+1+2+M-1)] + p.ALL[,c(N.eofs[mod_idx]+1+2+M):c(N.eofs[mod_idx]+1+2+2*M-1)])
                pN.samp <- expit(p.NAT[,c(N.eofs[mod_idx]+1+2):c(N.eofs[mod_idx]+1+2+M-1)] + p.NAT[,c(N.eofs[mod_idx]+1+2+M):c(N.eofs[mod_idx]+1+2+2*M-1)])
                RR.samp <- pA.samp/pN.samp
              }

              # Classify hypotheses and evaluate ======
              results_all[[true.state]][[scheme]][[ens.size]] <- classify_tests_evaluate(
                zA = zA, zN = zN, N = N.ens[ens.size], pA = pA, pN = pN, RR.samp = RR.samp, fit.LRT = FALSE,
                threshold = threshold, R1.alpha = R1.alpha, R2.lambda = R2.lambda, R3.alpha = R3.alpha )
            }
          }
        }
        save( results_all, file = sv_name )
      }
    } # End RNB and models 2-9

  }
}


#===============================================================
# Combine results
#===============================================================

# Aggregate over replicates (WRAF05)

# RNB fits
rnb.results.R1 <- combine_results2( mod.txt = "RNB", fit.txt = "RNB", row.idx = 4 )
rnb.results.R2 <- combine_results2( mod.txt = "RNB", fit.txt = "RNB", row.idx = 5 )
rnb.results.R3 <- combine_results2( mod.txt = "RNB", fit.txt = "RNB", row.idx = 6 )

# (M1) Independent fits
unadj.results <- combine_results( mod.txt = "M1", fit.txt = "Unadj (LRT)", row.idx = 1 )
fwer.results <- combine_results( mod.txt = "M1", fit.txt = "FWER (LRT)", row.idx = 2 )
BH.results <- combine_results( mod.txt = "M1", fit.txt = "BH (LRT)", row.idx = 3 )
indep.results.R1 <- combine_results( mod.txt = "M1", fit.txt = "(M1) Indep", row.idx = 4 )
indep.results.R2 <- combine_results( mod.txt = "M1", fit.txt = "(M1) Indep", row.idx = 5 )
indep.results.R3 <- combine_results( mod.txt = "M1", fit.txt = "(M1) Indep", row.idx = 6 )

# (M2) Gaussian fits
gauss.results.R1 <- combine_results( mod.txt = "M2", fit.txt = "(M2) Gauss", row.idx = 4 )
gauss.results.R2 <- combine_results( mod.txt = "M2", fit.txt = "(M2) Gauss", row.idx = 5 )
gauss.results.R3 <- combine_results( mod.txt = "M2", fit.txt = "(M2) Gauss", row.idx = 6 )

# (M3) skew-t fits
skewt.results.R1 <- combine_results( mod.txt = "M3", fit.txt = "(M3) Skew-t", row.idx = 4 )
skewt.results.R2 <- combine_results( mod.txt = "M3", fit.txt = "(M3) Skew-t", row.idx = 5 )
skewt.results.R3 <- combine_results( mod.txt = "M3", fit.txt = "(M3) Skew-t", row.idx = 6 )

# (M4) CAR fits
car.results.R1 <- combine_results2( mod.txt = "M4", fit.txt = "(M4) CAR", row.idx = 4 )
car.results.R2 <- combine_results2( mod.txt = "M4", fit.txt = "(M4) CAR", row.idx = 5 )
car.results.R3 <- combine_results2( mod.txt = "M4", fit.txt = "(M4) CAR", row.idx = 6 )

# (M5) hybrid CAR fits
carhyb.results.R1 <- combine_results2( mod.txt = "M5", fit.txt = "(M5) hybrid CAR", row.idx = 4 )
carhyb.results.R2 <- combine_results2( mod.txt = "M5", fit.txt = "(M5) hybrid CAR", row.idx = 5 )
carhyb.results.R3 <- combine_results2( mod.txt = "M5", fit.txt = "(M5) hybrid CAR", row.idx = 6 )

# (M6) GP fits
gp.results.R1 <- combine_results2( mod.txt = "M6", fit.txt = "(M6) GP", row.idx = 4 )
gp.results.R2 <- combine_results2( mod.txt = "M6", fit.txt = "(M6) GP", row.idx = 5 )
gp.results.R3 <- combine_results2( mod.txt = "M6", fit.txt = "(M6) GP", row.idx = 6 )

# (M7) EOF30 fits
eof30.results.R1 <- combine_results2( mod.txt = "M7", fit.txt = "(M7) EOF-30", row.idx = 4 )
eof30.results.R2 <- combine_results2( mod.txt = "M7", fit.txt = "(M7) EOF-30", row.idx = 5 )
eof30.results.R3 <- combine_results2( mod.txt = "M7", fit.txt = "(M7) EOF-30", row.idx = 6 )

# (M8) EOF10 fits
eof10.results.R1 <- combine_results2( mod.txt = "M8", fit.txt = "(M8) EOF-10", row.idx = 4 )
eof10.results.R2 <- combine_results2( mod.txt = "M8", fit.txt = "(M8) EOF-10", row.idx = 5 )
eof10.results.R3 <- combine_results2( mod.txt = "M8", fit.txt = "(M8) EOF-10", row.idx = 6 )

# (M9) EOF50 fits
eof50.results.R1 <- combine_results2( mod.txt = "M9", fit.txt = "(M9) EOF-50", row.idx = 4 )
eof50.results.R2 <- combine_results2( mod.txt = "M9", fit.txt = "(M9) EOF-50", row.idx = 5 )
eof50.results.R3 <- combine_results2( mod.txt = "M9", fit.txt = "(M9) EOF-50", row.idx = 6 )


# Combine
wraf05.R1 <- rbind( unadj.results, fwer.results, BH.results, indep.results.R1, gauss.results.R1,
                    skewt.results.R1, car.results.R1, carhyb.results.R1, gp.results.R1,
                    eof30.results.R1, eof10.results.R1, eof50.results.R1, rnb.results.R1 )
wraf05.R2 <- rbind( unadj.results, fwer.results, BH.results, indep.results.R2, gauss.results.R2,
                    skewt.results.R2, car.results.R2, carhyb.results.R2, gp.results.R2,
                    eof30.results.R2, eof10.results.R2, eof50.results.R2, rnb.results.R2 )
wraf05.R3 <- rbind( unadj.results, fwer.results, BH.results, indep.results.R3, gauss.results.R3,
                    skewt.results.R3, car.results.R3, carhyb.results.R3, gp.results.R3,
                    eof30.results.R3, eof10.results.R3, eof50.results.R3, rnb.results.R3 )

# save(wraf05.R1, wraf05.R2, wraf05.R3, file = "wraf05_simstud_results.RData")
