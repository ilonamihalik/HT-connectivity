
#############################################
# SCRIPT 3 - Modelling functional connectivity for bears among spawning salmon waterways in Haíɫzaqv (Heiltsuk) Territory, coastal British Columbia

# This code reproduces the steps for determining the Least Cost Paths from our current map output from Circuitscape. These pathways represent the pathways of highest conductance 
# (based on the cumulative current map) between each focal node.
#############################################

library(raster)
library(sf)
library(leastcostpath)
library(gdistance)
library(dplyr)

setwd("")

### 1) Bring in focal nodes (salmon points) and cumulative current output (from circuitscape) layers  --------------------------------

# Bring in points
reachpoints.tif <- raster("salmonnodes_BWG.asc")
reachpoints <- rasterToPoints(reachpoints.tif, spatial = TRUE)

current_map <- raster("BWGsalmon_cum_curmap.asc")

### 2) Create a transition layer from cumulative current layer
curr_as_transition <- transition(current_map, transitionFunction=mean, 8) # transition layer to account for diagonal connections from 8 neighborhood rule
curr_as_transition <- geoCorrection(curr_as_transition, type="c",multpl=F) #type=c for least cost paths


### 3) Calculate Least Cost Paths from cumulative current layer
current_LCP<-create_FETE_lcps(curr_as_transition,reachpoints,cost_distance=TRUE,parallel=TRUE)

##Save LCP as shapefile. Added corridor buffers in ArcMap (but could do in R). 
#writeOGR(current_LCP, dsn = ".", layer = "", driver = "ESRI Shapefile", overwrite = TRUE)
