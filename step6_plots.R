#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: plots for the final paper
#===============================================================

#===============================================================
# Figure 1 (and A.1, A.2): comparing loss functions
#===============================================================

# Use artifical "P-values", with M = 100 tests
M <- 100
rank <- (1:M)/100

set.seed(32)
test2 <- sort(rbeta(M, 2, 0.7)) # Clustered around 1
set.seed(22356)
test3 <- sort(rbeta(M, 0.7, 2)) # Clustered around 0
# set.seed(533704)
set.seed(18061)
test4 <- sort(rbeta(M, 0.35, 0.35)) # Clustered around 0 and 1

testAvg2 <- testSum2 <- testAvg3 <- testSum3 <- testAvg4 <- testSum4 <- rep(NA,M)
for(i in 1:M){
  testAvg2[i] <- mean(test2[1:i])
  testSum2[i] <- sum(test2[1:i])
  testAvg3[i] <- mean(test3[1:i])
  testSum3[i] <- sum(test3[1:i])
  testAvg4[i] <- mean(test4[1:i])
  testSum4[i] <- sum(test4[1:i])
}

for(i in 2:4){
  text2.x <- 0.12 
  text2.y <- 0.96 
  if( i == 2 ){
    testUse <- test2
    testSumUse <- testSum2
    testAvgUse <- testAvg2
    pdf(file.path( outputDir, "FigureA2_compare_mostly1.pdf"), height = 6.5, width = 9)
  }
  if( i == 3 ){
    testUse <- test3
    testSumUse <- testSum3
    testAvgUse <- testAvg3
    pdf(file.path( outputDir, "FigureA1_compare_mostly0.pdf"), height = 6.5, width = 9)
  }
  if( i == 4 ){
    testUse <- test4
    testSumUse <- testSum4
    testAvgUse <- testAvg4
    pdf(file.path( outputDir, "Figure1_compare_bimodal.pdf"), height = 6.5, width = 9)
  }

  par(mar = c(5,7,5,9), mfrow=c(1,1))
  hist( testUse, axes=FALSE, breaks = 10, col = "gray93", border = NA,
        xlab=NA, ylab=NA, main=NA, ylim=c(0,35))
  #-------------
  par(new = T)
  plot( testUse, testUse, pch=15, cex = 0.75, xlab="P(H = 0|Z)", ylab=NA, col=4, axes = FALSE, ylim=c(0,1), xlim=c(0,1) )
  grid( lwd = 1.5 )
  mtext(side = 2, line = 3, "lambda (relative loss of a FD vs FN)", col = 4)
  axis(side = 1, at = seq(-0.2,1.2,0.1) )
  axis(side = 3, at = c(-0.2,1.2), tck = 0)
  axis(side = 2, col = 4, col.ticks = 4, col.axis = 4, at = seq(-0.2,1.2,0.1),
       labels = c(NA, NA, "0","1/9","1/4","3/7","2/3","1","3/2","7/3", "4", "9", "Inf", NA, NA)[15:1] )
  points(testUse[testUse <= 0.2], testUse[testUse <= 0.2], pch="O", col=4, cex = 1)
  text(x = text2.x, y = text2.y, label = "Horizontal threshold for", 
       cex = 1, col="gray35", adj = 0)
  text(x = text2.x, y = text2.y, label = "\n\nrelative loss = 4", cex = 1, col=4, adj = 0)
  text(x = text2.x, y = text2.y, label = "\n\n\n\nand", cex = 1, col="gray35", adj = 0)
  text(x = text2.x, y = text2.y, label = "\n\n\n\n\n\ncumulative avg = 0.2", cex = 1, col=3, adj = 0)
  text(x = text2.x, y = text2.y, label = "\n\n\n\n\n\n\n\nand", cex = 1, col="gray35", adj = 0)
  text(x = text2.x, y = text2.y, label = "\n\n\n\n\n\n\n\n\n\ncumulative sum = 0.2M", cex = 1, col=2, adj = 0)
  #-------------
  par(new = T)
  plot( testUse, testAvgUse, cex = 0.75, xlab=NA, ylab=NA, type="p", 
        pch=17, col = 3, axes = FALSE, ylim = c(0,1) )
  points(testUse[testAvgUse <= 0.2], testAvgUse[testAvgUse <= 0.2], pch="O", col=3, cex = 1)
  axis(side = 4, col = 3, col.ticks = 3, col.axis = 3, at = seq(-0.1,1.1,0.1))
  mtext(side = 4, line = 5, "Cumulative average", col = 3, adj = 0)
  abline(h = 0.2, lwd=1.5 )
  #-------------
  par(new = T)
  plot( testUse, testSumUse, cex = 0.75, axes=F, xlab=NA, ylab=NA, type="p", 
        pch=16, col=2, ylim=c(0,M))
  points(testUse[testSumUse <= 0.2*M], testSumUse[testSumUse <= 0.2*M], pch="O", col=2, cex = 1)
  axis(side = 4, line = 2.5, col = 2, col.ticks = 2, col.axis = 2, at = seq(-10,200,10))
  mtext(side = 4, line = 5, "Cumulative sum", col = 2, adj = 1)
  dev.off()


}

#===============================================================
# Figure 2: correlations
#===============================================================
dist_cut <-  cut( c(lowerTriangle(wraf05_R3_dist)), 
                  breaks = seq(0, 12742, by = 400), dig.lab = 5)

corrDF_plot <- data.frame( dist = c(lowerTriangle(wraf05_R3_dist)),
                           dist_cut = dist_cut,
                           ALLcor = c(lowerTriangle(ALL_hotJan_cor)),
                           NATcor = c(lowerTriangle(NAT_hotJan_cor)) )

g1 <- ggplot( corrDF_plot, aes( x = dist_cut, y = ALLcor) ) +
  geom_boxplot( outlier.colour = NA ) + ylim(-0.6, 1) + ggtitle("Factual scenario") +
  geom_hline( yintercept = 0, color = "red", alpha = 0.6 ) +
  ylab("Correlation") + xlab(expression("Distance in " ~ R^3 ~ " (km)")) +
  scale_x_discrete( drop = FALSE,
    breaks = levels(dist_cut)[seq(1, length(levels(dist_cut)), 3)],
    labels = seq(0, 12742, by = 400)[seq(1, length(levels(dist_cut)), 3)] ) +
  theme(axis.text.x = element_text(angle=45, hjust=1) )
g2 <- ggplot( corrDF_plot, aes( x = dist_cut, y = NATcor) ) +
  geom_boxplot( outlier.colour = NA ) + ylim(-0.6, 1) + ggtitle("Counterfactual scenario") +
  geom_hline( yintercept = 0, color = "red", alpha = 0.6 ) +
  ylab("Correlation") + xlab(expression("Distance in " ~ R^3 ~ " (km)")) +
  scale_x_discrete( drop = FALSE,
    breaks = levels(dist_cut)[seq(1, length(levels(dist_cut)), 3)],
    labels = seq(0, 12742, by = 400)[seq(1, length(levels(dist_cut)), 3)] ) +
  theme(axis.text.x = element_text(angle=45, hjust=1) )

pdf(file.path( outputDir, "Figure2_corr_hotJan.pdf"), height = 4, width = 9)
grid.arrange(g1, g2, ncol = 2)
dev.off()

#===============================================================
# Figures A.3, A.4: EOFs
#===============================================================
tempDF.A <- merge(wraf05_plotDF,
                  data.frame(ShortName = wraf05_centroids$region,
                             eof1 = ALL_hotJan_EOFs[,1],
                             eof2 = ALL_hotJan_EOFs[,2],
                             eof3 = ALL_hotJan_EOFs[,3],
                             eof4 = ALL_hotJan_EOFs[,4] ), by = "ShortName")
tempDF.N <- merge(wraf05_plotDF,
                  data.frame(ShortName = wraf05_centroids$region,
                             eof1 = NAT_hotJan_EOFs[,1],
                             eof2 = NAT_hotJan_EOFs[,2],
                             eof3 = NAT_hotJan_EOFs[,3],
                             eof4 = NAT_hotJan_EOFs[,4] ), by = "ShortName")

gA1 <- plotEOFs( tempDF.A, tempDF.A$eof1, "EOF1" )
gA2 <- plotEOFs( tempDF.A, tempDF.A$eof2, "EOF2" )
gA3 <- plotEOFs( tempDF.A, tempDF.A$eof3, "EOF3" )
gA4 <- plotEOFs( tempDF.A, tempDF.A$eof4, "EOF4" )

pdf(file.path( outputDir, "FigureA3_ALL_hotJan_EOFs.pdf" ), height = 14, width = 11)
grid.arrange( gA1, gA2, gA3, gA4, ncol = 1 )
dev.off()

gN1 <- plotEOFs( tempDF.N, tempDF.N$eof1, "EOF1" )
gN2 <- plotEOFs( tempDF.N, tempDF.N$eof2, "EOF2" )
gN3 <- plotEOFs( tempDF.N, tempDF.N$eof3, "EOF3" )
gN4 <- plotEOFs( tempDF.N, tempDF.N$eof4, "EOF4" )

pdf(file.path( outputDir, "FigureA4_NAT_hotJan_EOFs.pdf" ), height = 14, width = 11)
grid.arrange( gN1, gN2, gN3, gN4, ncol = 1 )
dev.off()


#===============================================================
# Figure 3: simstudy results
#===============================================================

# Plots: WRAF05
FDRpower_only <- wraf05.R1[,c(1:4,10,12)] # pick out categorical variables and FDR, power
FDRpower_only <- cbind(rbind(FDRpower_only[,1:4], FDRpower_only[,1:4]), c(FDRpower_only$FDR, FDRpower_only$power))
names(FDRpower_only)[5] <- "Metric"
FDRpower_only$MetType <- c(rep("Mean FDR", nrow(FDRpower_only)/2), rep("Mean power", nrow(FDRpower_only)/2))

pdf(file.path( outputDir, "Figure3_wraf05_R1_results.pdf" ), height = 19, width = 15)
plot_R1_R3( plotDF = FDRpower_only[FDRpower_only$Fit != "Unadj (LRT)",], target = 0.1, met.type = c("Mean FDR", "Mean power") )
dev.off()

pdf(file.path( outputDir, "FigureA5_wraf05_R2_results.pdf" ), height = 19, width = 15)
plot_R2( plotDF = wraf05.R2[wraf05.R2$Fit != "Unadj (LRT)",] )
dev.off()

FDFN_only <- wraf05.R3[,c(1:4,5,6)] # pick out categorical variables and FD, FN
FDFN_only <- cbind(rbind(FDFN_only[,1:4], FDFN_only[,1:4]), c(FDFN_only$FD, FDFN_only$FN))
names(FDFN_only)[5] <- "Metric"
FDFN_only$MetType <- c(rep("Mean FD", nrow(FDFN_only)/2), rep("Mean FN", nrow(FDFN_only)/2))

pdf(file.path( outputDir, "FigureA6_wraf05_R3_results.pdf" ), height = 19, width = 15)
plot_R1_R3( plotDF = FDFN_only[FDFN_only$Fit != "Unadj (LRT)",], target = 23.7, met.type = c("Mean FD", "Mean FN") )
dev.off()

# Plots: WRAF2
load("data/wraf2_simstud_results.RData")

FDRpower_only2 <- wraf2.R1[,c(1:4,10,12)] # pick out categorical variables and FDR, power
FDRpower_only2 <- cbind(rbind(FDRpower_only2[,1:4], FDRpower_only2[,1:4]), c(FDRpower_only2$FDR, FDRpower_only2$power))
names(FDRpower_only2)[5] <- "Metric"
FDRpower_only2$MetType <- c(rep("Mean FDR", nrow(FDRpower_only2)/2), rep("Mean power", nrow(FDRpower_only2)/2))

pdf(file.path( outputDir, "FigureA7_wraf2_R1_results.pdf" ), height = 19, width = 15)
plot_R1_R3( plotDF = FDRpower_only2, target = 0.1, met.type = c("Mean FDR", "Mean power") )
dev.off()

pdf(file.path( outputDir, "FigureA8_wraf2_R2_results.pdf" ), height = 19, width = 15)
plot_R2( plotDF = wraf2.R2 )
dev.off()

FDFN_only2 <- wraf2.R3[,c(1:4,5,6)] # pick out categorical variables and FD, FN
FDFN_only2 <- cbind(rbind(FDFN_only2[,1:4], FDFN_only2[,1:4]), c(FDFN_only2$FD, FDFN_only2$FN))
names(FDFN_only2)[5] <- "Metric"
FDFN_only2$MetType <- c(rep("Mean FD", nrow(FDFN_only2)/2), rep("Mean FN", nrow(FDFN_only2)/2))

pdf(file.path( outputDir, "FigureA9_wraf2_R3_results.pdf" ), height = 19, width = 15)
plot_R1_R3( plotDF = FDFN_only2, target = 6.8, met.type = c("Mean FD", "Mean FN") )
dev.off()

#===============================================================
# Figure 4: Hot, wet results
#===============================================================

# HOT: Match regions with reject/fail to reject
hot.lab <- rep(NA, nrow(wraf05_plotDF))
for(i in 1:length(levels(wraf05_plotDF$ShortName))){
  if(sum(hot.results.df$region.names == levels(wraf05_plotDF$ShortName)[i])>0){
    hot.lab[wraf05_plotDF$ShortName == levels(wraf05_plotDF$ShortName)[i]] <- hot.results.df$H0.leq5[hot.results.df$region.names == levels(wraf05_plotDF$ShortName)[i]]
  }
}

# WET: Match regions with reject/fail to reject
wet.lab <- rep(NA, nrow(wraf05_plotDF))
for(i in 1:length(levels(wraf05_plotDF$ShortName))){
  if(sum(wet.results.df$region.names == levels(wraf05_plotDF$ShortName)[i])>0){
    wet.lab[wraf05_plotDF$ShortName == levels(wraf05_plotDF$ShortName)[i]] <- wet.results.df$H0.leq1[wet.results.df$region.names == levels(wraf05_plotDF$ShortName)[i]]

  }
}

pdf(file.path( outputDir, "Figure4_HotJan_WetMarch.pdf" ), height = 10, width = 11)
grid.arrange(
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( hot.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#fc4e2a") ) +
    scale_y_continuous( limits = c( -60, 85 ) ) +
    ggtitle("Conclusive evidence for 5X increase in probability of a hot January in 2015") +
    theme(panel.grid.major = element_blank(),
          text = element_text(size = 16),
          panel.grid.minor = element_blank(),
          panel.background = element_rect( fill = "#deebf7"),
          legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_blank(),
          axis.ticks = element_blank()),
  # Antarctica only
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( hot.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#fc4e2a") ) +
    coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90,-65) ) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_rect( fill = "#deebf7" ),
           legend.position = "none", axis.text = element_blank(),
           axis.ticks = element_blank() ),
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( wet.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#4292c6") ) +
    scale_y_continuous( limits = c( -60, 85 ) ) +
    ggtitle("Conclusive evidence for an increase in probability of a wet March in 2015") +
    theme(panel.grid.major = element_blank(),
          text = element_text(size = 16),
          panel.grid.minor = element_blank(),
          panel.background = element_rect( fill = "#deebf7"),
          legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_blank(),
          axis.ticks = element_blank()),
  # Antarctica only
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( wet.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#4292c6") ) +
    coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90,-65) ) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_rect( fill = "#deebf7" ),
           legend.position = "none", axis.text = element_blank(),
           axis.ticks = element_blank() ),  ncol = 2, widths = c(1.6, 0.25)
)
dev.off()

#===============================================================
# Figure 4a: Hot, wet results with BH
#===============================================================

# HOT: Match regions with reject/fail to reject
hot.lab <- rep(NA, nrow(wraf05_plotDF))
for(i in 1:length(levels(wraf05_plotDF$ShortName))){
  if(sum(hot.BH.results.leq1$region == levels(wraf05_plotDF$ShortName)[i])>0){
    hot.lab[wraf05_plotDF$ShortName == levels(wraf05_plotDF$ShortName)[i]] <- hot.BH.results.leq1$ind.reject[hot.BH.results.leq1$region == levels(wraf05_plotDF$ShortName)[i]]
  }
}

# WET: Match regions with reject/fail to reject
wet.lab <- rep(NA, nrow(wraf05_plotDF))
for(i in 1:length(levels(wraf05_plotDF$ShortName))){
  if(sum(wet.BH.results.leq1$region == levels(wraf05_plotDF$ShortName)[i])>0){
    wet.lab[wraf05_plotDF$ShortName == levels(wraf05_plotDF$ShortName)[i]] <- wet.BH.results.leq1$ind.reject[wet.BH.results.leq1$region == levels(wraf05_plotDF$ShortName)[i]]
    
  }
}

pdf(file.path( outputDir, "Figure4a_BH_HotJan_WetMarch.pdf" ), height = 10, width = 11)
grid.arrange(
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( hot.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#fc4e2a") ) +
    scale_y_continuous( limits = c( -60, 85 ) ) +
    ggtitle("Conclusive evidence for 5X increase in probability of a hot January in 2015") +
    theme(panel.grid.major = element_blank(),
          text = element_text(size = 16),
          panel.grid.minor = element_blank(),
          panel.background = element_rect( fill = "#deebf7"),
          legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_blank(),
          axis.ticks = element_blank()),
  # Antarctica only
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( hot.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#fc4e2a") ) +
    coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90,-65) ) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_rect( fill = "#deebf7" ),
           legend.position = "none", axis.text = element_blank(),
           axis.ticks = element_blank() ),
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( wet.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#4292c6") ) +
    scale_y_continuous( limits = c( -60, 85 ) ) +
    ggtitle("Conclusive evidence for an increase in probability of a wet March in 2015") +
    theme(panel.grid.major = element_blank(),
          text = element_text(size = 16),
          panel.grid.minor = element_blank(),
          panel.background = element_rect( fill = "#deebf7"),
          legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_blank(),
          axis.ticks = element_blank()),
  # Antarctica only
  ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
    geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                  color = "black", fill = "white" ) +
    geom_polygon( aes( fill = as.factor( wet.lab ) ) ) +
    geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
    scale_fill_manual( labels = c("No", "Yes"),
                       name = "", values = c("#969696", "#4292c6") ) +
    coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90,-65) ) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_rect( fill = "#deebf7" ),
           legend.position = "none", axis.text = element_blank(),
           axis.ticks = element_blank() ),  ncol = 2, widths = c(1.6, 0.25)
)
dev.off()

#===============================================================
# Figure 5: Wet results, multiple H0
#===============================================================

wet.results.df$mult.cat <- rep(NA, nrow(wet.results.df))
for(i in 1:nrow(wet.results.df)){
  if( wet.results.df$H0.out[i] == 1 ){ # Reject: RR > 3/2 or RR < 2/3
    if( wet.results.df$H0.leq1[i] == 1 ){
      wet.results.df$mult.cat[i] <- "Conclusive that both 2/3 < RR < 3/2 and RR > 1"
    }
    if( wet.results.df$H0.geq1[i] == 1 ){
      wet.results.df$mult.cat[i] <- "Conclusive that both 2/3 < RR < 3/2 and RR < 1"
    }
    if( wet.results.df$H0.leq1[i] == 0 & wet.results.df$H0.geq1[i] == 0 ){
      wet.results.df$mult.cat[i] <- "No difference: conclusive that 2/3 < RR < 3/2"
    }
  }
  else{
    if( wet.results.df$H0.leq1[i] == 1 ){
      ifelse( wet.results.df$H0.leq2[i] == 1,
              wet.results.df$mult.cat[i] <- "Conclusive that RR > 2",
              wet.results.df$mult.cat[i] <- "Conclusive that RR > 1")
    }
    if( wet.results.df$H0.geq1[i] == 1 ){
      ifelse( wet.results.df$H0.geq05[i] == 1,
              wet.results.df$mult.cat[i] <- "Conclusive that RR < 1/2",
              wet.results.df$mult.cat[i] <- "Conclusive that RR < 1")
    }
  }
  if( is.na( wet.results.df$mult.cat[i]) ){
    wet.results.df$mult.cat[i] <- "Inconclusive"
  }
}

wet.results.df$mult.cat <- ordered(wet.results.df$mult.cat,
                                   levels = c("Conclusive that RR < 1/2", "Conclusive that RR < 1",
                                              "Conclusive that both 2/3 < RR < 3/2 and RR < 1",
                                              "No difference: conclusive that 2/3 < RR < 3/2",
                                              "Conclusive that both 2/3 < RR < 3/2 and RR > 1",
                                              "Conclusive that RR > 1", "Conclusive that RR > 2", "Inconclusive") )

# Relabel
wet.results.df$mult.cat2 <- as.character(wet.results.df$mult.cat)
wet.results.df$mult.cat2[wet.results.df$mult.cat == "Conclusive that RR < 1/2"] <- "Conclusive for 2X decrease"
wet.results.df$mult.cat2[wet.results.df$mult.cat == "Conclusive that RR < 1"] <- "Conclusive for decrease"
wet.results.df$mult.cat2[wet.results.df$mult.cat == "Conclusive that both 2/3 < RR < 3/2 and RR < 1"] <- "Most likely no change; some evidence of decrease  "
wet.results.df$mult.cat2[wet.results.df$mult.cat == "No difference: conclusive that 2/3 < RR < 3/2"] <- "Conclusive for no change"
wet.results.df$mult.cat2[wet.results.df$mult.cat == "Conclusive that both 2/3 < RR < 3/2 and RR > 1"] <- "Most likely no change; some evidence of increase"
wet.results.df$mult.cat2[wet.results.df$mult.cat == "Conclusive that RR > 1"] <- "Conclusive for increase"
wet.results.df$mult.cat2[wet.results.df$mult.cat == "Conclusive that RR > 2"] <- "Conclusive for 2X increase"
wet.results.df$mult.cat2 <- ordered(wet.results.df$mult.cat2,
                                    levels = c("Conclusive for 2X decrease", "Conclusive for decrease",
                                               "Most likely no change; some evidence of decrease  ",
                                               "Conclusive for no change",
                                               "Most likely no change; some evidence of increase",
                                               "Conclusive for increase", "Conclusive for 2X increase", "Inconclusive") )

# Match regions with reject/fail to reject
wet.lab2 <- rep(NA, nrow(wraf05_plotDF))
for(i in 1:length(levels(wraf05_plotDF$ShortName))){
  if(sum(wet.results.df$region.names == levels(wraf05_plotDF$ShortName)[i])>0){
    wet.lab2[wraf05_plotDF$ShortName == levels(wraf05_plotDF$ShortName)[i]] <- as.character(wet.results.df$mult.cat2[wet.results.df$region.names == levels(wraf05_plotDF$ShortName)[i]])

  }
}
wet.lab2 <- ordered(wet.lab2,
                    levels = c("Conclusive for 2X decrease",
                               "Conclusive for decrease",
                               "Most likely no change; some evidence of decrease  ",
                               "Conclusive for no change",
                               "Most likely no change; some evidence of increase",
                               "Conclusive for increase",
                               "Conclusive for 2X increase",
                               "Inconclusive") )

# Ordered as with wet.lab2
col.use <- c( "#08519c", "#4292c6", "#74c476", "#ffff33", "#fd8d3c", "#ef3b2c", "#a50f15",
              "#969696" )

wetplot <- ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
  geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                color = "black", fill = "white" ) +
  geom_polygon( aes( fill = as.factor( wet.lab2 ) ) ) +
  geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
  scale_fill_manual( name = "", drop = FALSE,
                     values = col.use[c(7:1,8)] ) +
  scale_y_continuous( limits = c( -60, 85 ) ) +
  ggtitle("Conclusive evidence for changes in probability of a wet March in 2015") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 16),
        panel.grid.minor = element_blank(),
        panel.background = element_rect( fill = "#deebf7"),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + guides( fill = guide_legend( ncol = 2 ) )

wetplot2 <- ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
  geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                color = "black", fill = "white" ) +
  geom_polygon( aes( fill = as.factor( wet.lab2 ) ) ) +
  geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
  scale_fill_manual( name = "", drop = FALSE,
                     values = col.use[c(7:1,8)] ) +
  coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90,-65) ) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect( fill = "#deebf7" ), axis.text = element_blank(),
         axis.ticks = element_blank() ) + guides( fill = guide_legend( ncol = 2 ) )

# extract legend
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

wetlegend <- g_legend(wetplot)

# Version 2: wider range for "no change"
wet.results.df$mult.cat.v2 <- rep(NA, nrow(wet.results.df))
for(i in 1:nrow(wet.results.df)){
  if( wet.results.df$H0.out2[i] == 1 ){ # Reject: RR > 2 or RR < 1/2
    if( wet.results.df$H0.leq1[i] == 1 ){
      wet.results.df$mult.cat.v2[i] <- "Conclusive that both 1/2 < RR < 2 and RR > 1"
    }
    if( wet.results.df$H0.geq1[i] == 1 ){
      wet.results.df$mult.cat.v2[i] <- "Conclusive that both 1/2 < RR < 2 and RR < 1"
    }
    if( wet.results.df$H0.leq1[i] == 0 & wet.results.df$H0.geq1[i] == 0 ){
      wet.results.df$mult.cat.v2[i] <- "No difference: conclusive that 1/2 < RR < 2"
    }
  }
  else{
    if( wet.results.df$H0.leq1[i] == 1 ){
      ifelse( wet.results.df$H0.leq2[i] == 1,
              wet.results.df$mult.cat.v2[i] <- "Conclusive that RR > 2",
              wet.results.df$mult.cat.v2[i] <- "Conclusive that RR > 1")
    }
    if( wet.results.df$H0.geq1[i] == 1 ){
      ifelse( wet.results.df$H0.geq05[i] == 1,
              wet.results.df$mult.cat.v2[i] <- "Conclusive that RR < 1/2",
              wet.results.df$mult.cat.v2[i] <- "Conclusive that RR < 1")
    }
  }
  if( is.na( wet.results.df$mult.cat.v2[i]) ){
    wet.results.df$mult.cat.v2[i] <- "Inconclusive"
  }
}

wet.results.df$mult.cat.v2 <- ordered(wet.results.df$mult.cat.v2,
                                      levels = c("Conclusive that RR < 1/2", "Conclusive that RR < 1",
                                                 "Conclusive that both 1/2 < RR < 2 and RR < 1",
                                                 "No difference: conclusive that 1/2 < RR < 2",
                                                 "Conclusive that both 1/2 < RR < 2 and RR > 1",
                                                 "Conclusive that RR > 1", "Conclusive that RR > 2", "Inconclusive") )

# Relabel
wet.results.df$mult.cat.v22 <- as.character(wet.results.df$mult.cat.v2)
wet.results.df$mult.cat.v22[wet.results.df$mult.cat.v2 == "Conclusive that RR < 1/2"] <- "Conclusive for 2X decrease"
wet.results.df$mult.cat.v22[wet.results.df$mult.cat.v2 == "Conclusive that RR < 1"] <- "Conclusive for decrease"
wet.results.df$mult.cat.v22[wet.results.df$mult.cat.v2 == "Conclusive that both 1/2 < RR < 2 and RR < 1"] <- "Most likely no change; some evidence of decrease  "
wet.results.df$mult.cat.v22[wet.results.df$mult.cat.v2 == "No difference: conclusive that 1/2 < RR < 2"] <- "Conclusive for no change"
wet.results.df$mult.cat.v22[wet.results.df$mult.cat.v2 == "Conclusive that both 1/2 < RR < 2 and RR > 1"] <- "Most likely no change; some evidence of increase"
wet.results.df$mult.cat.v22[wet.results.df$mult.cat.v2 == "Conclusive that RR > 1"] <- "Conclusive for increase"
wet.results.df$mult.cat.v22[wet.results.df$mult.cat.v2 == "Conclusive that RR > 2"] <- "Conclusive for 2X increase"
wet.results.df$mult.cat.v22 <- ordered(wet.results.df$mult.cat.v22,
                                       levels = c("Conclusive for 2X decrease", "Conclusive for decrease",
                                                  "Most likely no change; some evidence of decrease  ",
                                                  "Conclusive for no change",
                                                  "Most likely no change; some evidence of increase",
                                                  "Conclusive for increase", "Conclusive for 2X increase", "Inconclusive") )

# Match regions with reject/fail to reject
wet.lab22 <- rep(NA, nrow(wraf05_plotDF))
for(i in 1:length(levels(wraf05_plotDF$ShortName))){
  if(sum(wet.results.df$region.names == levels(wraf05_plotDF$ShortName)[i])>0){
    wet.lab22[wraf05_plotDF$ShortName == levels(wraf05_plotDF$ShortName)[i]] <- as.character(wet.results.df$mult.cat.v22[wet.results.df$region.names == levels(wraf05_plotDF$ShortName)[i]])

  }
}
wet.lab22 <- ordered(wet.lab22,
                     levels = c("Conclusive for 2X decrease",
                                "Conclusive for decrease",
                                "Most likely no change; some evidence of decrease  ",
                                "Conclusive for no change",
                                "Most likely no change; some evidence of increase",
                                "Conclusive for increase",
                                "Conclusive for 2X increase",
                                "Inconclusive") )


wetplotv2 <- ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
  geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                color = "black", fill = "white" ) +
  geom_polygon( aes( fill = as.factor( wet.lab22 ) ) ) +
  geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
  scale_fill_manual( name = "", drop = FALSE,
                     values = col.use[c(7:1,8)] ) +
  scale_y_continuous( limits = c( -60, 85 ) ) +
  ggtitle("Conclusive evidence for changes in probability of a wet March in 2015") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 16),
        panel.grid.minor = element_blank(),
        panel.background = element_rect( fill = "#deebf7"),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + guides( fill = guide_legend( ncol = 2 ) )

wetplot2v2 <- ggplot( wraf05_plotDF ) + aes( long, lat, group = group ) +
  geom_polygon( data = map_data( "world", interior = FALSE ), aes( x = long, y = lat, group = group ),
                color = "black", fill = "white" ) +
  geom_polygon( aes( fill = as.factor( wet.lab22 ) ) ) +
  geom_path( color = "black", size = 0.35 ) + coord_equal() + ylab(NULL) + xlab(NULL) +
  scale_fill_manual( name = "", drop = FALSE,
                     values = col.use[c(7:1,8)] ) +
  coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90,-65) ) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect( fill = "#deebf7" ), axis.text = element_blank(),
         axis.ticks = element_blank() ) + guides( fill = guide_legend( ncol = 2 ) )

# extract legend
wetlegend2 <- g_legend(wetplotv2)

pdf(file.path( outputDir, "Figure5_wetMarch_v1_v2.pdf" ), height = 10.2, width = 11)
grid.arrange(arrangeGrob(wetplotv2 + theme(legend.position="none"),
                         wetplot2v2 + theme(legend.position="none"),
                         wetplot + theme(legend.position="none"),
                         wetplot2 + theme(legend.position="none"),
                         ncol = 2, widths = c(1.6, 0.25) ),
             wetlegend2, nrow=2, heights=c(2, 0.25))
dev.off()



