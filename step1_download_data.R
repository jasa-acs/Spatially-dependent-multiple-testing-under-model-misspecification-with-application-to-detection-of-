#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: download data
# Accessing CAM5.1 runs and aggregating to the WRAF regions
#===============================================================

# The data used for the analyses in Section 4 of the paper
# consist of simulations of the Community Atmospheric 
# Model, version 5.1, run in its ~1 degree configuration.
# The simulations are provided as part of the C20C+ Detection
# and Attribution project, a subproject of the World Climate 
# Research Programme's (WCRP) Climate Variability Programme's 
# (CLIVAR) Climate of the 20th Century Plus Project (C20C+).
# Full details on the simulations are available at
# http://portal.nersc.gov/c20c/
#
# The files are available through the online data portal. 
# Simulations are aggregated monthly for each of the ALL and
# NAT scenarios, with one file for each of 400 model runs,
# separately for surface temperature (tas) and precipitation
# (pr).
# Details on the runs are as follows:
# 001:010, 036:050, 061:070, 086:100 cover 1959/01 - 2015/06
#   Exceptions: run001 for ALL is 1959/01 - 2014/09
#               run002 for NAT is 1959/01 - 2014/10
# 011:035, 051:060, 071:85 cover 1996/01 - 2015/06
#   Exceptions: run025 for NAT is 1959/01 - 2014/10
# 101:400 cover 2010/01 - 2013/12
#
# The data from each run is on a longitude/latitude grid of
# size 288 x 192. A separate file of region weights identifies
# which grid cells comprise each of the WRAF05 regions (see 
# below). Finally, since we are calculating area averages, a 
# slight adjustment is made (using "cell_area") because the grid
# cells have varying area as one moves from the equator to the
# poles.
#
# NOTE: the following code is extremely time-consuming. A pre-
# processed file consisting of the final output 
# ("mon_atmos_WRAF05.RData") is available as supplementary 
# material at [url blinded].
#===============================================================

# Constants
M <- 237 # Total number of regions
N.mth <- 678 # Total number of months (1959/01-2015/06)
N.runs <- 400 # Total number of runs

# Storage: array dimensions are (month) x (region) x (run)
ALL_tas_raw <- NAT_tas_raw <- ALL_pr_raw <- NAT_pr_raw <- array(NA, dim=c(N.mth, M, N.runs))

# Region weights
download.file( url = paste("http://portal.nersc.gov/c20c/data/LBNL/CAM5-1-1degree/All-Hist/est1/v2-0/fx/atmos-WRAF05-v4-0/region/run000/region_fx-WRAF05-v4-0_CAM5-1-1degree_All-Hist_est1_v2-0_run000_000000-000000.nc", sep="" ),
               destfile = "region_weights.nc", mode = "wb" )
region_wts_nc <- nc_open("region_weights.nc")
region_wts <- ncvar_get(region_wts_nc, "region")
N.lon <- dim(region_wts)[1] # 288 lon grid points (approx 1 degree grid)
N.lat <- dim(region_wts)[2] # 192 lat grid points (approx 1 degree grid)
lat.vals <- ncvar_get(region_wts_nc, "lat")

# Area of grid cells changes with latitude: calculate an adjustment term
cell_area <- matrix(NA, N.lon, N.lat)
for(j in 1:N.lat){
  cell_area[,j] <- cos(lat.vals[j]/180*pi)
}
cell_area <- cell_area/sum(cell_area) # Normalize

# Function to aggregate over regions
region_mean <- function(x){ sum(x * reg_wt * cell_area, na.rm = TRUE)/sum(reg_wt * cell_area, na.rm = TRUE) }

# Download and process ALL tas and pr (loop over runs)
nersc.ALL.dir <- "http://portal.nersc.gov/c20c/data/LBNL/CAM5-1-1degree/All-Hist/est1/v2-0/mon/atmos/"
for(rn in 1:N.runs){
  
  # Setup
  if( rn %in% 1:9 ){ rn.txt <- paste("00", rn, sep="") }
  if( rn %in% 10:99 ){ rn.txt <- paste("0", rn, sep="") }
  if( rn >= 100 ){ rn.txt <- paste(rn) }
  
  # Pick out months included
  if( rn %in% c(2:10,36:50,61:70,86:100) ){
    ALL.strt <- c(1, 649, 661, 670) # corresponding to 1959/01, 2013/01, 2014/01, 2014/10
    ALL.end <- c(648, 660, 669, 678) # corresponding to 2012/12, 2013/12, 2014/09, 2015/06
    ALL.file.txt <- c("195901-201212.nc", "201301-201312.nc", "201401-201409.nc", "201410-201506.nc")
  }
  if( rn %in% c(1) ){
    ALL.strt <- c(1, 649, 661) # corresponding to 1959/01, 2013/01, 2014/01
    ALL.end <- c(648, 660, 669) # corresponding to 2012/12, 2013/12, 2014/09
    ALL.file.txt <- c("195901-201212.nc", "201301-201312.nc", "201401-201409.nc")
  }
  if( rn %in% c(11:35,51:60,71:85) ){
    ALL.strt <- c(445, 649, 661, 670) # corresponding to 1996/01, 2013/01, 2014/01, 2014/10
    ALL.end <- c(648, 660, 669, 678) # corresponding to 2012/12, 2013/12, 2014/09, 2015/06
    ALL.file.txt <- c("199601-201212.nc", "201301-201312.nc", "201401-201409.nc", "201410-201506.nc")
  }
  if( rn %in% 101:400 ){
    ALL.strt <- 613 # corresponding to 2010/01
    ALL.end <- 660 # corresponding to 2013/12
    ALL.file.txt <- "201001-201312.nc"
  }
  
  for(k in 1:length(ALL.file.txt)){
    download.file( url = paste(nersc.ALL.dir, "pr/run", rn.txt, 
                               "/pr_Amon_CAM5-1-1degree_All-Hist_est1_v2-0_run",
                               rn.txt, "_", ALL.file.txt[k], sep="" ), mode = "wb", 
                   destfile = paste( getwd(), "/pr_temp_ALL.nc", sep="" ) )
    download.file( url = paste(nersc.ALL.dir, "tas/run", rn.txt, 
                               "/tas_Amon_CAM5-1-1degree_All-Hist_est1_v2-0_run",
                               rn.txt, "_", ALL.file.txt[k], sep="" ), mode = "wb",
                   destfile = paste( getwd(), "/tas_temp_ALL.nc", sep="" ) )
    
    # TEMPERATURE
    ALLnc_tas <- nc_open("tas_temp_ALL.nc")
    ALL_tas <- ncvar_get(ALLnc_tas, "tas")

    # PRECIPITATION
    ALLnc_pr <- nc_open("pr_temp_ALL.nc")
    ALL_pr <- ncvar_get(ALLnc_pr, "pr")

    # Aggregate by region
    for( m in 1:M ){ # Loop over regions
      reg_wt <- region_wts[,,m]
      tas_temp <- apply( ALL_tas, 3, region_mean ) # Aggregate all months
      pr_temp <- apply( ALL_pr, 3, region_mean ) # Aggregate all months

      # Store
      ALL_tas_raw[ALL.strt[k]:ALL.end[k],m,rn] <- tas_temp
      ALL_pr_raw[ALL.strt[k]:ALL.end[k],m,rn] <- pr_temp
    }
    
    # Clean-up
    nc_close(ALLnc_tas)
    nc_close(ALLnc_pr)
  }
  cat(rn, " ")

}

# Download and process NAT tas and pr (loop over runs)
nersc.NAT.dir <- "http://portal.nersc.gov/c20c/data/LBNL/CAM5-1-1degree/Nat-Hist/CMIP5-est1/v2-0/mon/atmos/"
for(rn in 1:N.runs){
  
  # Setup
  if( rn %in% 1:9 ){ rn.txt <- paste("00", rn, sep="") }
  if( rn %in% 10:99 ){ rn.txt <- paste("0", rn, sep="") }
  if( rn >= 100 ){ rn.txt <- paste(rn) }
  
  # Pick out months included
  if( rn %in% c(1,3:10,36:50,61:70,86:90,96:98) ){
    NAT.strt <- c(1, 649, 661, 671) # corresponding to 1959/01, 2013/01, 2014/01, 2014/11
    NAT.end <- c(648, 660, 670, 678) # corresponding to 2012/12, 2013/12, 2014/10, 2015/06
    NAT.file.txt <- c( "195901-201212.nc", "201301-201312.nc", "201401-201410.nc", "201411-201506.nc")
  }
  if( rn %in% c(2) ){
    NAT.strt <- c(1, 649, 661) # corresponding to 1959/01, 2013/01, 2014/01, 2014/11
    NAT.end <- c(648, 660, 670) # corresponding to 2012/12, 2013/12, 2014/10, 2015/06
    NAT.file.txt <- c( "195901-201212.nc", "201301-201312.nc", "201401-201410.nc" )
  }
  if( rn %in% c(91:95,99:100) ){
    NAT.strt <- c(1, 649, 661, 671, 672) # corresponding to 1959/01, 2013/01, 2014/01, 2014/11, 2014/12
    NAT.end <- c(648, 660, 670, 671, 678) # corresponding to 2012/12, 2013/12, 2014/10, 2014/11, 2015/06
    NAT.file.txt <- c( "195901-201212.nc", "201301-201312.nc", "201401-201410.nc", "201411-201411.nc", "201412-201506.nc")
  }
  if( rn %in% c(11:24,26:35,51:60,71:85) ){
    NAT.strt <- c(445, 649, 661, 671) # corresponding to 1996/01, 2013/01, 2014/01, 2014/11
    NAT.end <- c(648, 660, 670, 678) # corresponding to 2012/12, 2013/12, 2014/10, 2015/06
    NAT.file.txt <- c( "199601-201212.nc", "201301-201312.nc", "201401-201410.nc", "201411-201506.nc")
  }
  if( rn %in% c(25) ){
    NAT.strt <- c(445, 649, 661) # corresponding to 1996/01, 2013/01, 2014/01
    NAT.end <- c(648, 660, 670) # corresponding to 2012/12, 2013/12, 2014/10
    NAT.file.txt <- c( "199601-201212.nc", "201301-201312.nc", "201401-201410.nc" )
  }
  if( rn %in% 101:400 ){
    NAT.strt <- 613 # corresponding to 2010/01
    NAT.end <- 660 # corresponding to 2013/12
    NAT.file.txt <- "201001-201312.nc"
  }
  
  for(k in 1:length(NAT.file.txt)){
    download.file( url = paste(nersc.NAT.dir, "pr/run", rn.txt, 
                               "/pr_Amon_CAM5-1-1degree_Nat-Hist_CMIP5-est1_v2-0_run",
                               rn.txt, "_", NAT.file.txt[k], sep="" ), mode = "wb",
                   destfile = paste( getwd(), "/pr_temp_NAT.nc", sep="" ) )
    download.file( url = paste(nersc.NAT.dir, "tas/run", rn.txt, 
                               "/tas_Amon_CAM5-1-1degree_Nat-Hist_CMIP5-est1_v2-0_run",
                               rn.txt, "_", NAT.file.txt[k], sep="" ), mode = "wb",
                   destfile = paste( getwd(), "/tas_temp_NAT.nc", sep="" ) )
    
    # TEMPERATURE
    NATnc_tas <- nc_open("tas_temp_NAT.nc")
    NAT_tas <- ncvar_get(NATnc_tas, "tas")
    
    # PRECIPITATION
    NATnc_pr <- nc_open("pr_temp_NAT.nc")
    NAT_pr <- ncvar_get(NATnc_pr, "pr")
    
    # Aggregate by region
    for( m in 1:M ){ # Loop over regions
      reg_wt <- region_wts[,,m]
      
      if(rn %in% c(91:95,99:100) & k == 4 ){
        tas_temp <- sum(NAT_tas * reg_wt * cell_area, na.rm = TRUE)/sum(reg_wt * cell_area, na.rm = TRUE)
        pr_temp <- sum(NAT_pr * reg_wt * cell_area, na.rm = TRUE)/sum(reg_wt * cell_area, na.rm = TRUE)
      }
      else{
        tas_temp <- apply( NAT_tas, 3, region_mean ) # Aggregate all months
        pr_temp <- apply( NAT_pr, 3, region_mean ) # Aggregate all months
      }
      
      # Store
      NAT_tas_raw[NAT.strt[k]:NAT.end[k],m,rn] <- tas_temp
      NAT_pr_raw[NAT.strt[k]:NAT.end[k],m,rn] <- pr_temp
    }
    
    # Clean-up
    nc_close(NATnc_tas)
    nc_close(NATnc_pr)
  }
  cat(rn, " ")
}

# Add variable names
region_names <- ncvar_get(region_wts_nc, "nameshort")

mnth.names <- paste(rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct",
                          "Nov","Dec"),56), rep(1959:2014, each=12),sep="")
mnth.names <- c(mnth.names, "Jan2015", "Feb2015", "Mar2015","Apr2015","May2015","Jun2015")

dimnames(ALL_tas_raw) <- list( mnth.names, region_names, paste("run", 1:N.runs, sep = "") )
dimnames(ALL_pr_raw) <- list( mnth.names, region_names, paste("run", 1:N.runs, sep = "") )
dimnames(NAT_tas_raw) <- list( mnth.names, region_names, paste("run", 1:N.runs, sep = "") )
dimnames(NAT_pr_raw) <- list( mnth.names, region_names, paste("run", 1:N.runs, sep = "") )