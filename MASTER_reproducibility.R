#===============================================================
# Spatially-Dependent Multiple Testing Under Model 
# Misspecification, with Application to Detection of 
# Anthropogenic Influence on Extreme Climate Events
# [Authors and affiliation blinded]
# January, 2018
# RESUBMISSION, version 2
#===============================================================

#===============================================================
# Reproducibility: MASTER script
#===============================================================

# This file contains instructions for reproducing the data, all
# analyses, and plots contained in the paper. 

# Download the repository from Bitbucket [url blinded].
# The following script assumes the working directory has
# been set to this folder. All of the necessary data and code
# are now available.

# As denoted below, Steps 1 and 4 are computationally 
# intensive and take a *very* long time to run. As such, pre-
# processed data files are available for these steps. Otherwise,
# code within an individual step (e.g., plotting in step 6) 
# assumes code for the previous steps has been run.

#===============================================================
# Step 0: Install necessary packages; setup
# Mirror CRAN on January 23, 2018 using the checkpoint package,
# as well as the R version used for all analyses
#===============================================================
library(checkpoint)
checkpoint( snapshotDate = "2018-01-23", R.version = "3.3.3" )

# Necessary packages
library(ncdf4)
library(maptools)
library(rgeos)
library(plyr)
library(gridExtra)
library(ggplot2)
library(ggmap)
library(geosphere)
library(colorspace)
library(RColorBrewer)
library(nimble)
library(geoR)
library(MASS)
library(ggalt)
library(gdata)

# Output directory (for plots)
outputDir <- "output"

#===============================================================
# Step 1: step1_download_data.R
# Master script to download and aggregate the CAM5.1 runs to
# the WRAF05 regions.
#===============================================================

# !!! NOTE: RUNNING THIS SCRIPT IS EXTREMELY TIME-CONSUMING !!!
# As an alternative, a pre-processed .RData file is available in 
# the /data folder ("mon_atmos_WRAF05.RData"). Data in this file
# are:
# ALL_tas_raw, ALL_pr_raw, NAT_tas_raw, NAT_pr_raw
# Each of these are a 3-dimensional array, containing monthly
# surface temperature (tas) and precipitation rate (pr) in the
# factual (ALL) and counterfactual (NAT) scenarios, for
# each WRAF05 region and each climate ensemble member. The 
# arrays are of size 678 x 237 x 400 (month x region x run).

# Raw data
source("step1_download_data.R")

# Alternatively: load pre-processed data
# load("data/mon_atmos_WRAF05.RData")

#===============================================================
# Step 2: step2_data_preparation.R
# Master script to conduct data processing for the analyses.
#===============================================================

source("step2_data_preparation.R")

# Note that there are several warning messages that do not adversely
# affect the code:
### Warning messages:
###   1: use rgdal::readOGR or sf::st_read 
###   2: use rgdal::readOGR or sf::st_read 

#===============================================================
# Now, all behind-the-scenes functions can be loaded.
#===============================================================
source("code_functions/functions_pt1.R") # Multiple testing functions
source("code_functions/functions_pt2.R") # Generate EOFs
source("code_functions/functions_pt3.R") # Plotting and summary functions

#===============================================================
# Step 3: step3_simulate_data.R
# Master script to generate simulated data.
#===============================================================

source("step3_simulate_data.R")

#===============================================================
# Step 4: source nimble_code.R and run step4_fit_simstud.R
# Compile nimble models and run the simulation study.
#===============================================================

# !!! NOTE: Running step4_fit_simstud.R will take a LONG time !!!
#
# This step involves running MCMC for 
# 10 models x 100 replicates x 6 true states x 3 schemes  
#   x 3 ensemble sizes x 2 scenarios = 108000 Bayesian models,
# so running the code in serial would likely take many days, or
# even possibly weeks.
# In practice, these 108000 tasks were broken up by model and
# replicate (each comprising 6 x 3 x 3 x 2 = 108 MCMC runs). The 
# MCMC for each group of tasks was completed using parallel 
# computing on the laboratory supercomputer (the necessary code
# for doing so is omitted for brevity).

# Note: these files are only for fitting the WRAF05 regions (M=237).

source("step4_fit_simstud.R")

# Alternatively: load pre-processed data
# load("data/wraf05_simstud_results.RData")

#===============================================================
# Step 5: step5_case_studies.R
# Master script to fit models and summarize the case studies.
#===============================================================

source("step5_case_studies.R")

#===============================================================
# Step 6: step6_plots.R
# Script to generate all plots in the paper draft.
#===============================================================

source("step6_plots.R")
