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
# Part 3: plotting functions; code to summarize the simulation
# study fits
#===============================================================

plotEOFs <- function( eofDF, eofplot, leg.text ){
  g <- grid.arrange(
    ggplot( eofDF ) + aes( long, lat, group = group ) + 
      geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),  
                    color = "black", fill = "white" ) +
      geom_polygon( aes( fill = eofplot ) ) + 
      geom_path( color = "black" ) + coord_equal() + ylab(NULL) + xlab(NULL) +
      scale_fill_gradientn( colours = brewer.pal(11, "RdBu")[11:1], name = leg.text,
                            limits = c(-max(abs(eofplot)), max(abs(eofplot))) ) + 
      scale_y_continuous( limits = c( -60, 85 ) ) +
      theme( panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_rect( fill = "#deebf7" ),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             legend.position = "left") + 
      guides(fill = guide_colorbar(barheight = 10)) ,
    ggplot( eofDF ) + aes( long, lat, group = group ) + 
      geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),  
                    color = "black", fill = "white" ) +
      geom_polygon( aes( fill = eofplot ) ) + 
      scale_fill_gradientn( colours = brewer.pal(11, "RdBu")[11:1], name = leg.text,
                            limits = c(-max(abs(eofplot)), max(abs(eofplot))) ) + 
      geom_path( color = "black" ) + coord_equal() + ylab(NULL) + xlab(NULL) +
      coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90,-65) ) +
      theme( panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_rect( fill = "#deebf7" ),
             legend.position = "none", axis.text = element_blank(),
             axis.ticks = element_blank() ),
    ncol = 2, widths = c(1.6, 0.25)
    
  )
  return(g)
}

combine_results <- function( mod.txt, fit.txt, 
                          truestate.txt = c("G-RE","NG-RE","GP-S","GP-L","EOF-G","EOF-NG"),
                          row.idx = 5  ){
  
  n.ts <- 6 # length(results.list)
  n.ens <- 3 # length(results.list[[1]][[1]])
  n.schm <- 3 # length(results.list[[1]])
  n.reps <- 100 # length(results.list[[1]][[1]][[1]]) 
  
  FD.mat <- FN.mat <- TD.mat <- TN.mat <- array(NA, dim=c(n.reps, n.ens, n.schm, n.ts))
  Loss.mat <- FDR.mat <- FNR.mat <- power.mat <- array(NA, dim=c(n.reps, n.ens, n.schm, n.ts))
  for(l in 1:n.reps){ # Replicates
    load(paste(mod.txt, "_replicate", l, "_results.RData", sep = ""))
    for(i in 1:n.ts){ # True state
     for(j in 1:n.ens){ # Ensemble size
      for(k in 1:n.schm){ # Scheme
          FD.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$FD[row.idx]
          FN.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$FN[row.idx]
          TD.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$TD[row.idx]
          TN.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$TN[row.idx]
          Loss.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$Loss[row.idx]
          FDR.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$FDR[row.idx]
          FNR.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$FNR[row.idx]
          power.mat[l,j,k,i] <- results_all[[i]][[k]][[j]]$power[row.idx]
        }
      }
    }
  }
  
  FD.vec <- FN.vec <- TD.vec <- TN.vec <- NULL
  Loss.vec <- FDR.vec <- FNR.vec <- power.vec <- NULL
  for( j in 1:n.ens ){
    for( k in 1:n.schm ){
      FD.vec <- c( FD.vec, colMeans( FD.mat[ , j, k, ] ) )
      FN.vec <- c( FN.vec, colMeans( FN.mat[ , j, k, ] ) )
      TD.vec <- c( TD.vec, colMeans( TD.mat[ , j, k, ] ) )
      TN.vec <- c( TN.vec, colMeans( TN.mat[ , j, k, ] ) )
      FDR.vec <- c( FDR.vec, colMeans( FDR.mat[ , j, k, ] ) )
      FNR.vec <- c( FNR.vec, colMeans( FNR.mat[ , j, k, ] ) )
      Loss.vec <- c( Loss.vec, colMeans( Loss.mat[ , j, k, ] ) )
      power.vec <- c( power.vec, colMeans( power.mat[ , j, k, ] ) )
    }
  }
  
  output <- data.frame(
    Fit = rep( fit.txt, n.schm*n.ens*n.ts ),
    TrueState = ordered( rep( truestate.txt, n.schm*n.ens ), levels = truestate.txt ),
    Scheme = rep( c( rep( "Scheme 1", n.ts ), 
                     rep( "Scheme 2", n.ts ), 
                     rep( "Scheme 3", n.ts) ), n.ens ),
    EnsSize = factor( c( rep( 50, n.schm*n.ts ), rep( 100, n.schm*n.ts), 
                         rep( 400, n.schm*n.ts ) ) ),
    FD = FD.vec, FN = FN.vec, TD = TD.vec, TN = TN.vec,
    Loss = Loss.vec, FDR = FDR.vec, FNR = FNR.vec, power = power.vec
  )
  
  return(output)
  
}

combine_results2 <- function( mod.txt, fit.txt, 
                             truestate.txt = c("G-RE","NG-RE","GP-S","GP-L","EOF-G","EOF-NG"),
                             row.idx = 5  ){
  
  n.ts <- 6 # length(results.list)
  n.ens <- 3 # length(results.list[[1]][[1]])
  n.schm <- 3 # length(results.list[[1]])
  n.reps <- 100 # length(results.list[[1]][[1]][[1]]) 
  
  FD.mat <- FN.mat <- TD.mat <- TN.mat <- array(NA, dim=c(n.reps, n.ens, n.schm, n.ts))
  Loss.mat <- FDR.mat <- FNR.mat <- power.mat <- array(NA, dim=c(n.reps, n.ens, n.schm, n.ts))
  for(l in 1:n.reps){ # Replicates
    for(i in 1:n.ts){ # True state
      for(j in 1:n.ens){ # Ensemble size
        for(k in 1:n.schm){ # Scheme
          
          load(paste(mod.txt, "_replicate", l, 
                     "_trueState", i,
                     "_scheme", k,
                     "_ensSize", j,
                     "_results.RData", sep = ""))
          
          FD.mat[l,j,k,i] <- results_all$FD[row.idx]
          FN.mat[l,j,k,i] <- results_all$FN[row.idx]
          TD.mat[l,j,k,i] <- results_all$TD[row.idx]
          TN.mat[l,j,k,i] <- results_all$TN[row.idx]
          Loss.mat[l,j,k,i] <- results_all$Loss[row.idx]
          FDR.mat[l,j,k,i] <- results_all$FDR[row.idx]
          FNR.mat[l,j,k,i] <- results_all$FNR[row.idx]
          power.mat[l,j,k,i] <- results_all$power[row.idx]
        }
      }
    }
  }
  
  FD.vec <- FN.vec <- TD.vec <- TN.vec <- NULL
  Loss.vec <- FDR.vec <- FNR.vec <- power.vec <- NULL
  for( j in 1:n.ens ){
    for( k in 1:n.schm ){
      FD.vec <- c( FD.vec, colMeans( FD.mat[ , j, k, ] ) )
      FN.vec <- c( FN.vec, colMeans( FN.mat[ , j, k, ] ) )
      TD.vec <- c( TD.vec, colMeans( TD.mat[ , j, k, ] ) )
      TN.vec <- c( TN.vec, colMeans( TN.mat[ , j, k, ] ) )
      FDR.vec <- c( FDR.vec, colMeans( FDR.mat[ , j, k, ] ) )
      FNR.vec <- c( FNR.vec, colMeans( FNR.mat[ , j, k, ] ) )
      Loss.vec <- c( Loss.vec, colMeans( Loss.mat[ , j, k, ] ) )
      power.vec <- c( power.vec, colMeans( power.mat[ , j, k, ] ) )
    }
  }
  
  output <- data.frame(
    Fit = rep( fit.txt, n.schm*n.ens*n.ts ),
    TrueState = ordered( rep( truestate.txt, n.schm*n.ens ), levels = truestate.txt ),
    Scheme = rep( c( rep( "Scheme 1", n.ts ), 
                     rep( "Scheme 2", n.ts ), 
                     rep( "Scheme 3", n.ts) ), n.ens ),
    EnsSize = factor( c( rep( 50, n.schm*n.ts ), rep( 100, n.schm*n.ts), 
                         rep( 400, n.schm*n.ts ) ) ),
    FD = FD.vec, FN = FN.vec, TD = TD.vec, TN = TN.vec,
    Loss = Loss.vec, FDR = FDR.vec, FNR = FNR.vec, power = power.vec
  )
  
  return(output)
  
}


plot_R1_R3 <- function( plotDF, target = 0.1, met.type = c("Mean FDR", "Mean power") ){
  grid.arrange(
    ggplot( plotDF[plotDF$Scheme == "Scheme 1",], aes(x = Fit, y = Metric, group = EnsSize)) + 
      facet_grid( MetType ~ TrueState, scales = "free_y" ) +
      geom_hline( data = data.frame(MetType = met.type, Metric = c(target, NA) ), 
                  aes( yintercept = Metric ), color = "gray40") +
      geom_line(aes(color = EnsSize)) + geom_point(aes(color = EnsSize, shape = EnsSize), size = 2.5) + 
      scale_color_manual(values = brewer.pal(5,"Set1"), name = "Ensemble size") + xlab(NULL) + ylab(NULL) +
      theme( legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
             strip.text = element_text(size = 14), axis.text.y = element_text(size = 12),
             title = element_text(size = 14), legend.text = element_text(size = 16),
             legend.title = element_text(size = 16), panel.spacing.y = unit(1, "lines")) +
      ggtitle("Scheme 1: approximately 85% true rejections"),
    ggplot( plotDF[plotDF$Scheme == "Scheme 2",], aes(x = Fit, y = Metric, group = EnsSize)) + 
      facet_grid( MetType ~ TrueState, scales = "free_y" ) +
      geom_hline( data = data.frame(MetType = met.type, Metric = c(target, NA) ), 
                  aes( yintercept = Metric ), color = "gray40") +
      geom_line(aes(color = EnsSize)) + geom_point(aes(color = EnsSize, shape = EnsSize), size = 2.5) + 
      scale_color_manual(values = brewer.pal(5,"Set1"), name = "Ensemble size") + xlab(NULL) + ylab(NULL) +
      theme( legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
             strip.text = element_text(size = 14), axis.text.y = element_text(size = 12),
             title = element_text(size = 14), legend.text = element_text(size = 16),
             legend.title = element_text(size = 16), panel.spacing.y = unit(1, "lines")) +
      ggtitle("Scheme 2: approximately 50% true rejections"),
    ggplot( plotDF[plotDF$Scheme == "Scheme 3",], aes(x = Fit, y = Metric, group = EnsSize)) + 
      facet_grid( MetType ~ TrueState, scales = "free_y" ) +
      geom_hline( data = data.frame(MetType = met.type, Metric = c(target, NA) ), 
                  aes( yintercept = Metric ), color = "gray40") +
      geom_line(aes(color = EnsSize)) + geom_point(aes(color = EnsSize, shape = EnsSize), size = 2.5) + 
      scale_color_manual(values = brewer.pal(5,"Set1"), name = "Ensemble size") + xlab(NULL) + ylab(NULL) +
      theme( legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
             strip.text = element_text(size = 14), axis.text.y = element_text(size = 12),
             title = element_text(size = 14), legend.text = element_text(size = 16),
             legend.title = element_text(size = 16), panel.spacing.y = unit(1, "lines")) +
      guides(colour = guide_legend( override.aes = list( shape = c(16,17,15) ) )) + scale_shape( guide = FALSE) +
      ggtitle("Scheme 3: approximately 15% true rejections"),
    ncol = 1, heights = c(1,1,1.1)
  )
}

plot_R2 <- function( plotDF ){
  grid.arrange(
    ggplot( plotDF[plotDF$Scheme == "Scheme 1",], aes(x = Fit, y = Loss, group = EnsSize)) + 
      facet_grid( . ~ TrueState ) + geom_line(aes(color = EnsSize)) + geom_point(aes(color = EnsSize, shape = EnsSize), size = 2.5) + 
      ggtitle("Scheme 1: approximately 85% true rejections") +
      scale_color_manual(values = brewer.pal(5,"Set1"), name = "Ensemble size") + xlab(NULL) + ylab("Mean Loss\n") +
      theme( legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
             strip.text = element_text(size = 14), axis.text.y = element_text(size = 12),
             title = element_text(size = 14), legend.text = element_text(size = 16),
             legend.title = element_text(size = 16), panel.spacing.y = unit(1, "lines")) +
      ggtitle("Scheme 1: (approximately 85% true rejections)"),
    ggplot( plotDF[plotDF$Scheme == "Scheme 2",], aes(x = Fit, y = Loss, group = EnsSize)) + 
      facet_grid( . ~ TrueState ) + geom_line(aes(color = EnsSize)) + 
      geom_point(aes(color = EnsSize, shape = EnsSize), size = 2.5) + ggtitle("Scheme 2: approximately 50% true rejections") +
      scale_color_manual(values = brewer.pal(5,"Set1"), name = "Ensemble size") + xlab(NULL) + ylab("Mean Loss\n") +
      theme( legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
             strip.text = element_text(size = 14), axis.text.y = element_text(size = 12),
             title = element_text(size = 14), legend.text = element_text(size = 16),
             legend.title = element_text(size = 16), panel.spacing.y = unit(1, "lines")) +
      ggtitle("Scheme 2: (approximately 50% true rejections)"),
    ggplot( plotDF[plotDF$Scheme == "Scheme 3",], aes(x = Fit, y = Loss, group = EnsSize)) + 
      facet_grid( . ~ TrueState ) + geom_line(aes(color = EnsSize)) + geom_point(aes(color = EnsSize, shape = EnsSize), size = 2.5) + 
      ggtitle("Scheme 3: approximately 15% true rejections") +
      scale_color_manual(values = brewer.pal(5,"Set1"), name = "Ensemble size") + xlab(NULL) + ylab("Mean Loss\n") +
      theme( legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
             strip.text = element_text(size = 14), axis.text.y = element_text(size = 12),
             title = element_text(size = 14), legend.text = element_text(size = 16),
             legend.title = element_text(size = 16), panel.spacing.y = unit(1, "lines")) +
      guides(colour = guide_legend( override.aes = list( shape = c(16,17,15) ) )) + scale_shape( guide = FALSE) +
      ggtitle("Scheme 3: (approximately 15% true rejections)"),
    ncol = 1, heights = c(1,1,1.1)
  )
}
