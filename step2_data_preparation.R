#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: data perparation
# Data preparation for the simulation study and case studies.
#===============================================================

#===============================================================
# Step (A): create a data frame for plotting the WRAF05 regions
#===============================================================

# Load the shape file
wraf.shp <- readShapeSpatial("data/shape_files/region_fx-WRAF05-v4-0_WRAF_All-Hist_est1_v4-0_run000_000000-000000.shp")
reg.names <- wraf.shp$ShortName
wraf.shp@data$id <- rownames(wraf.shp@data)
wraf.shp.points <- fortify(wraf.shp, region="id")
wraf.shp.df <- join(wraf.shp.points, wraf.shp@data, by="id")
wraf.shp@data$id <- rownames(wraf.shp@data)
wraf.shp.points <- fortify(wraf.shp, region="id")

# Create a data frame with points and an ordering
wraf.shp.df <- join(wraf.shp.points, wraf.shp@data, by="id")

# Thin the points to reduce file size of plots
wraf05_plotDF <- wraf.shp.df[wraf.shp.df$order %% 20 == 0,1:12]

# Adjust longitude of the Chukotka region
wraf05_plotDF$long[wraf05_plotDF$ShortName == "Chukotka" & wraf05_plotDF$long < 0] <- (
  wraf05_plotDF$long[wraf05_plotDF$ShortName == "Chukotka" & wraf05_plotDF$long < 0] + 360
)

#===============================================================
# Step (B): calculate the lat/lon centroid of each WRAF region,
# as well as sets up the adjacency matrix needed for the CAR 
# model M3.
#===============================================================

# Download the region weight info  (if not already done)
if( !exists("region_wts_nc") ){
  download.file( url = paste("http://portal.nersc.gov/c20c/data/LBNL/CAM5-1-1degree/All-Hist/est1/v2-0/fx/atmos-WRAF05-v4-0/region/run000/region_fx-WRAF05-v4-0_CAM5-1-1degree_All-Hist_est1_v2-0_run000_000000-000000.nc", sep="" ),
                 destfile = "region_weights.nc" )
  region_wts_nc <- nc_open("region_weights.nc")
}
wraf_region <- ncvar_get(region_wts_nc, "region") # 237 total regions
wraf_lon <- ncvar_get(region_wts_nc, "lon")
wraf_lon[wraf_lon > 180] <- wraf_lon[wraf_lon > 180] - 360
wraf_lat <- ncvar_get(region_wts_nc, "lat")
wraf_nameshort <- ncvar_get(region_wts_nc, "nameshort")
M <- length(wraf_nameshort)

# Region centroids
# Pick out grid cells in each region that have more than 0.33 weight
wraf_long <- NULL
for(k in 1:length(wraf_nameshort)){
  for( i in 1:dim(wraf_region)[1]){
    for( j in 1:dim(wraf_region)[2]){
      if( is.na(wraf_region[i,j,k]) == FALSE & wraf_region[i,j,k] > 0.33 ){
        wraf_long <- rbind(wraf_long, c(wraf_lon[i],wraf_lat[j],wraf_nameshort[k]))
      }
    }
  }
}

wraf05_map <- data.frame( longitude = as.numeric(wraf_long[,1]),
  latitude = as.numeric(wraf_long[,2]), region = as.factor(wraf_long[,3]) )
# Far eastern Russian region that straddles the -180 line
wraf05_map$longitude[wraf05_map$region == wraf_nameshort[195] & wraf05_map$longitude < 0 ] <- wraf05_map$longitude[wraf05_map$region == wraf_nameshort[195] & wraf05_map$longitude < 0 ] + 360

# Calculate the central (median) latitude and longitude
cent.lat <- cent.lon <- rep(NA,M)
for(i in 1:M){
  cent.lat[i] <- median(wraf05_map$latitude[wraf05_map$region == wraf_nameshort[i]])
  cent.lon[i] <- median(wraf05_map$longitude[wraf05_map$region == wraf_nameshort[i]])
}

# Combine into a data frame
wraf05_centroids <- data.frame( region = wraf_nameshort, central.lon = cent.lon, central.lat = cent.lat )

# R^3 distance between centroids
xyz.crds <- matrix(NA,nrow(wraf05_centroids),3)
# Transform degrees to radians
mean.lat.radians <- wraf05_centroids$central.lat*(pi/180)
mean.lon.radians <- wraf05_centroids$central.lon*(pi/180)

for(i in 1:nrow(wraf05_centroids)){
  xyz.crds[i,1] <- 6371*cos(mean.lat.radians[i])*cos(mean.lon.radians[i])
  xyz.crds[i,2] <- 6371*cos(mean.lat.radians[i])*sin(mean.lon.radians[i])
  xyz.crds[i,3] <- 6371*sin(mean.lat.radians[i])
}
wraf05_R3_dist <- as.matrix(dist(xyz.crds, diag=TRUE, upper=TRUE))
rownames(wraf05_R3_dist) <- colnames(wraf05_R3_dist) <- wraf05_centroids$region

# Neighbors: calculate the adjacency matrix for CAR model
Wmat05 <- matrix(NA,M,M)
for(i in 1:M){
  Wmat05[i,i] <- 0
  if(i < M){
    for(j in (i+1):M){
      # Only search over regions that are "close"
      if(wraf05_R3_dist[i,j] < 6000){
        # Count the number of shared pixels (weights for both are > 0.2)
        shrd.px <- sum( wraf_region[,,i] > 0.2 & is.na(wraf_region[,,i]) == FALSE
                        & wraf_region[,,j] > 0.2 & is.na(wraf_region[,,j]) == FALSE )
        ifelse(shrd.px > 0, Wmat05[i,j] <- 1, Wmat05[i,j] <- 0)
      }
      else{
        Wmat05[i,j] <- 0
      }
      Wmat05[j,i] <- Wmat05[i,j]
    }
  }
}
rownames(Wmat05) <- colnames(Wmat05) <- wraf_nameshort

# Adjustments: link together continents
Wmat05_final <- Wmat05

# 5 Papua New Guinea: add 232 Far North and North Queensland
Wmat05_final[5,232] <- Wmat05_final[232,5] <- 1

# 10 east GCC: remove 193 North Arab League
Wmat05_final[10,193] <- Wmat05_final[193,10] <- 0

# 51 Iberia: add 189 Morocco
Wmat05_final[51,189] <- Wmat05_final[189,51] <- 1

# 48 Madagascar: add 188 Mozambique
Wmat05_final[48,188] <- Wmat05_final[188,48] <- 1

# 94 Kalimantan: add 153 west Indonesia, 154 east Indonesia
Wmat05_final[94,153] <- Wmat05_final[153,94] <- 1
Wmat05_final[94,154] <- Wmat05_final[154,94] <- 1

# 195 Chukotka: add 209 north Alaska
Wmat05_final[195,209] <- Wmat05_final[209,195] <- 1

# 216 Namibia: remove 217 north central SADC
Wmat05_final[216,217] <- Wmat05_final[217,216] <- 0

# 219 Northeast Greenland NP: add 133 northwest EEA
Wmat05_final[219,133] <- Wmat05_final[133,219] <- 1

# 235 East ACD: add 90 northern east China and 163 Heilongjiang and Jilin
Wmat05_final[235,90] <- Wmat05_final[90,235] <- 1
Wmat05_final[235,163] <- Wmat05_final[163,235] <- 1

# 17 South Patagonia: add 181 Antarctic peninsula
Wmat05_final[17,181] <- Wmat05_final[181,17] <- 1

# Check
mean(rowSums(Wmat05_final))
sum(rowSums(Wmat05_final) - colSums(Wmat05_final))
Wmat05 <- Wmat05_final
rownames(Wmat05) <- colnames(Wmat05) <- wraf_nameshort

# # Check that corresponding precision is rank M-1
# prec.eig <- eigen(diag(rowSums(Wmat05)) - Wmat05)
# prec.eig$values

#===============================================================
# Step (C): using the aggregated CAM5.1 runs, calculate
# empirical correlations and EOFs.
# Also: calculate the binomial variates used in Section 4 for
# the appication
#===============================================================

# The following code calculates the threshold, count variables,
# exceedance probability estimates, empirical correlation, and
# EOFs for wet Marches and hot Januarys, for both ALL and NAT 
# scenarios. 
#
# Note: these quantities are calculated only for 1959-2014;
# the 2015 simulations are used for the application.
#
# The same quantities for other events/months can be calculated
# similarly.

# Constants
M <- 237 # Total number of regions
N.yr.rep <- length(1959:2014)

# Storage
ALL_hotcountJan <- NAT_hotcountJan <- matrix( NA, M, N.yr.rep )
ALL_wetcountMarch <- NAT_wetcountMarch <- matrix( NA, M, N.yr.rep )
hot_thresh <- wet_thresh <- rep( NA, M ) # Separate threshold for each region

for(k in 1:M){
  
  # Pick out data from month/region of interest
  ALL_tas_sub <- ALL_tas_raw[1 + 12*(0:55),k,] # Pick out all the Januarys for region k
  NAT_tas_sub <- NAT_tas_raw[1 + 12*(0:55),k,]
  
  ALL_pr_sub <- ALL_pr_raw[3 + 12*(0:55),k,] # Pick out all the Marches for region k
  NAT_pr_sub <- NAT_pr_raw[3 + 12*(0:55),k,]
  
  # Calculate anomalies (make each year mean-zero)
  ALL_tas_anom <- ALL_tas_sub - matrix(rowMeans(ALL_tas_sub, na.rm=TRUE), ncol = 1) %x% matrix(rep(1, 400), nrow = 1)
  NAT_tas_anom <- NAT_tas_sub - matrix(rowMeans(NAT_tas_sub, na.rm=TRUE), ncol = 1) %x% matrix(rep(1, 400), nrow = 1)
  
  ALL_pr_anom <- ALL_pr_sub - matrix(rowMeans(ALL_pr_sub, na.rm=TRUE), ncol = 1) %x% matrix(rep(1, 400), nrow = 1)
  NAT_pr_anom <- NAT_pr_sub - matrix(rowMeans(NAT_pr_sub, na.rm=TRUE), ncol = 1) %x% matrix(rep(1, 400), nrow = 1)
  
  # Calculate the threshold: based on the ALL simulations from 1985 to
  # 2014 (i.e., years 27-56), using only the 50 ensemble members that 
  # cover this entire period.
  hot_thresh[k] <- quantile( ALL_tas_anom[27:56,c(1:10,36:50,61:70,86:100)], 0.9, na.rm = TRUE )
  wet_thresh[k] <- quantile( ALL_pr_anom[27:56,c(1:10,36:50,61:70,86:100)], 0.9, na.rm = TRUE )
  
  # Calculate count variables 
  for(j in 1:N.yr.rep){
    
    ALL_hotcountJan[k,j] <- sum( ALL_tas_anom[j,] > hot_thresh[k], na.rm = TRUE )
    ALL_wetcountMarch[k,j] <- sum( ALL_pr_anom[j,] > wet_thresh[k], na.rm = TRUE )
    
    NAT_hotcountJan[k,j] <- sum( NAT_tas_anom[j,] > hot_thresh[k], na.rm = TRUE )
    NAT_wetcountMarch[k,j] <- sum( NAT_pr_anom[j,] > wet_thresh[k], na.rm = TRUE )
    
  }
}

# Calculate the ensemble sizes for each year
n.ens <- rowSums(is.na(ALL_pr_anom) == FALSE)

# Calculate Bayesian posterior means using M1
# Z ~ binomial(n.ens, p), p ~ beta(1,1)
# p|Z ~ beta(1 + Z, n.ens - Z + 1)
# Want: E(p|Z) = (Z + 1)/(n.ens + 2)
ALL_hotprobJan <- NAT_hotprobJan <- matrix( NA, M, N.yr.rep )
ALL_wetprobMarch <- NAT_wetprobMarch <- matrix( NA, M, N.yr.rep )
for(k in 1:M){
  for(j in 1:N.yr.rep){
    ALL_hotprobJan[k,j] <- (ALL_hotcountJan[k,j] + 1)/(n.ens[j] + 2)
    NAT_hotprobJan[k,j] <- (NAT_hotcountJan[k,j] + 1)/(n.ens[j] + 2)
    
    ALL_wetprobMarch[k,j] <- (ALL_wetcountMarch[k,j] + 1)/(n.ens[j] + 2)
    NAT_wetprobMarch[k,j] <- (NAT_wetcountMarch[k,j] + 1)/(n.ens[j] + 2)
  }
}

# Calculate empirical correlation and corresponing EOFs
ALL_hotJan_cor <- cor( t(ALL_hotprobJan) ) 
NAT_hotJan_cor <- cor( t(NAT_hotprobJan) ) 
ALL_wetMarch_cor <- cor( t(ALL_wetprobMarch) ) 
NAT_wetMarch_cor <- cor( t(NAT_wetprobMarch) ) 

ALL_hotJan_EOFs <- eigen( cov( t(ALL_hotprobJan) ) )$vectors[,1:56] # Only 56 unique eigenvectors
NAT_hotJan_EOFs <- eigen( cov( t(NAT_hotprobJan) ) )$vectors[,1:56] # Only 56 unique eigenvectors
ALL_wetMarch_EOFs <- eigen( cov( t(ALL_wetprobMarch) ) )$vectors[,1:56] # Only 56 unique eigenvectors
NAT_wetMarch_EOFs <- eigen( cov( t(NAT_wetprobMarch) ) )$vectors[,1:56] # Only 56 unique eigenvectors

# Finally, calculate the count variables for 2015
ALL_hotJan2015_counts <- NAT_hotJan2015_counts <- rep(NA,M)
ALL_wetMarch2015_counts <- NAT_wetMarch2015_counts <- rep(NA,M)

# 2015/01 = month 673; 2015/03 = month 675
for(k in 1:M){
  hotThresh <- quantile( ALL_tas_raw[1 + 12*(26:55),k,c(1:10,36:50,61:70,86:100)], 0.9, na.rm = TRUE )
  wetThresh <- quantile( ALL_pr_raw[3 + 12*(26:55),k,c(1:10,36:50,61:70,86:100)], 0.9, na.rm = TRUE )
  
  ALL_hotJan2015_counts[k] <- sum( ALL_tas_raw[673,k,] > hotThresh, na.rm = TRUE )
  ALL_wetMarch2015_counts[k] <- sum( ALL_pr_raw[675,k,] > wetThresh, na.rm = TRUE )
  
  NAT_hotJan2015_counts[k] <- sum( NAT_tas_raw[673,k,] > hotThresh, na.rm = TRUE )
  NAT_wetMarch2015_counts[k] <- sum( NAT_pr_raw[675,k,] > wetThresh, na.rm = TRUE )
}
n.ALL2015 <- sum(is.na(ALL_tas_raw[673,k,]) == FALSE) # Same as month 675
n.NAT2015 <- sum(is.na(NAT_tas_raw[673,k,]) == FALSE) # Same as month 675

# Name and save
casestudy_df <- data.frame(
  region = dimnames(ALL_tas_raw)[[2]],
  nA_2015 = rep(n.ALL2015, M), nN_2015 = rep(n.NAT2015, M),
  zA_hotJan2015 = ALL_hotJan2015_counts,
  zA_wetMarch2015 = ALL_wetMarch2015_counts,
  zN_hotJan2015 = NAT_hotJan2015_counts,
  zN_wetMarch2015 = NAT_wetMarch2015_counts
)

# Add region names 
rownames(ALL_hotJan_cor) <- dimnames(ALL_tas_raw)[[2]]
rownames(NAT_hotJan_cor) <- dimnames(ALL_tas_raw)[[2]]
rownames(ALL_wetMarch_cor) <- dimnames(ALL_tas_raw)[[2]]
rownames(NAT_wetMarch_cor) <- dimnames(ALL_tas_raw)[[2]]
colnames(ALL_hotJan_cor) <- dimnames(ALL_tas_raw)[[2]]
colnames(NAT_hotJan_cor) <- dimnames(ALL_tas_raw)[[2]]
colnames(ALL_wetMarch_cor) <- dimnames(ALL_tas_raw)[[2]]
colnames(NAT_wetMarch_cor) <- dimnames(ALL_tas_raw)[[2]]

rownames(ALL_hotJan_EOFs) <- dimnames(ALL_tas_raw)[[2]]
rownames(NAT_hotJan_EOFs) <- dimnames(ALL_tas_raw)[[2]]
rownames(ALL_wetMarch_EOFs) <- dimnames(ALL_tas_raw)[[2]]
rownames(NAT_wetMarch_EOFs) <- dimnames(ALL_tas_raw)[[2]]