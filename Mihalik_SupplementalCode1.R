#############################################
# SCRIPT 1 - Modelling functional connectivity for bears among spawning salmon waterways in Haíɫzaqv (Heiltsuk) Territory, coastal British Columbia

# This code subsets from the BC Vegetation Resource Inventory Database, and 
# prepares the 'Landcover class' and 'Water' spatial layers for our cumulative resistance surface.
# Note - the creation of our cumulative resistance surface was done in both R and ArcMap.
#############################################

setwd("")

library(sf)
library(tidyverse)
library(bcdata)

### 1) Subset from Vegetation Resource Inventory spatial data, collect subset ------------------------------------------
#Bring in VRI
VRI_all <- bcdc_query_geodata("2ebb35d8-c82f-4a17-9c96-612ac3532d55") %>%
   filter(INTERSECTS(midcoastTSA_sf)) %>% # Subset to central coast region to reduce space
   dplyr::select(FEATURE_ID, 
                 BCLCS_LEVEL_3, 
                 BCLCS_LEVEL_2, 
                 BCLCS_LEVEL_4,
                 BCLCS_LEVEL_1,
                 HARVEST_DATE,
                 VRI_LIVE_S,
                 PROJ_AGE_1, 
                 ) %>%
   collect()

VRI_all <- sf::st_read("./VRI 2023/VRI_midcoast_2024.shp")
VRI_all <- sf::st_transform(VRI_all, crs = 3005)

#Create 'Years since harvest' column to determine age of clearcut
VRI_all$HARVEST_YR <- format(as.Date(VRI_all$HARVEST_DATE, format = "%d/%m/%Y"), "%Y")
VRI_all$YR_SINCE_HARVST <- (2021 - as.numeric(VRI_all$HARVEST_YR))

#Subset only Treed-coniferous, -broadleaf, or --mixed forests
VRI_all_treesub <- VRI_all %>%
  dplyr::select(BCLCS_LEVEL_3, YR_SINCE_HARVST, HARVEST_YR, PROJ_AGE_1) %>%
  dplyr::filter(BCLCS_LEVEL_3 %in% c("TB", "TM", "TC"))

# 2) Create 'Mature Forest' layer ---------------------------------------------------
#Subset based on the following criteria:
#Years since harvest >= 75, OR = NA, OR Projected Age >= 75
VRI_forests_mature <- VRI_all_treesub %>%
  dplyr::filter(is.na(YR_SINCE_HARVST)|YR_SINCE_HARVST >= 75|PROJ_AGE_1 >= 75)

write_sf(VRI_forests_mature, "./mature_forests_VRI.shp")


# 3) Create 'Regenerating Clearcut' layer ---------------------------------------------------
# Subset based on the following criteria:
# Years since harvest >= 10 AND <= 75

VRI_forests_regenerating <- VRI_all_treesub %>%
  dplyr::filter(YR_SINCE_HARVST >= 10 & YR_SINCE_HARVST < 75)

write_sf(VRI_forests_regenerating, "./regenerating_forests_VRI.shp")


# 4) Create 'Bryoids, Shrublands, and Herbs' layer ---------------------------------------------------
# Subset based on following criteria:
# BCLCS Level 3 != TC, TB, TM, SI, RO, EL

VRI_shrubsherbs_all <- VRI_all %>%
  dplyr::filter(!BCLCS_LEVEL_3 %in% c("TB", "TM", "TC", "SI", "RO", "EL")) %>%
  dplyr::filter(!is.na(BCLCS_LEVEL_3))

write_sf(VRI_shrubsherbs_all, "./shrubland_VRI.shp", overwrite = TRUE)


# 5) Create 'Barren, Exposed Rock' layer ---------------------------------------------------
#Subset based on following criteria:
#BCLCS Level 3 = RO, EL
VRI_exposed_rock <- VRI_all %>%
  dplyr::filter(BCLCS_LEVEL_3 %in% c("RO", "EL")) %>%
  dplyr::filter(!is.na(BCLCS_LEVEL_3))

write_sf(VRI_exposed_rock, "./exposed_rock_VRI.shp")


# 6) Create 'Snow, Ice' layer ---------------------------------------------------
#Subset based on following criteria:
#BCLCS Level 3 = SI
VRI_snow_ice <- VRI_all %>%
  dplyr::filter(BCLCS_LEVEL_3 == "SI") %>%
  dplyr::filter(!is.na(BCLCS_LEVEL_3))

write_sf(VRI_snow_ice, "./Snow_Ice_VRI.shp")


# 7) Create Water layer ----------------------------------------------------------------------
#To be used in ArcMap to calculate 'distance to shore' using multiple ring buffer analysis 

# Subset only water
VRI_water <- VRI_all %>%
  dplyr::filter(BCLCS_LEVEL_1 == "W") %>%
  dplyr::filter(!is.na(BCLCS_LEVEL_1))

write_sf(VRI_water, "./Snow_Ice_VRI.shp")

# 8) [THIS STEP WAS DONE IN ARCMAP] Dissolved each Landcover layer created above in ArcMap ---------------------------------------

# 9) Bring in dissolved Landcover layers to combine

VRI_forests <- read_sf("./mature_forests_VRI_dissolved.shp", stringsAsFactors = TRUE)
VRI_clearcut <- read_sf("./regenerating_forests_VRI_dissolved.shp", stringsAsFactors = TRUE)
VRI_rock <- read_sf("./exposed_rock_VRI_dissolved.shp", stringsAsFactors = TRUE)
VRI_shrubs <- read_sf("./shrubland_VRI_dissolve.shp", stringsAsFactors = TRUE)
VRI_snow <- read_sf("./Snow_Ice_VRI_dissolve.shp", stringsAsFactors = TRUE)
VRI_water <- read_sf("./wateronly_VRI_buffers_noland.shp", stringsAsFactors = TRUE)

# Add new column to each to specify (and differentiate) type
VRI_forests$type <- 100
VRI_clearcut$type <- 200
VRI_rock$type <- 300
VRI_shrubs$type <- 400
VRI_snow$type <- 500
VRI_water$type <- 600

# Combine Mature Forests, Regenerating Forest, Snow_Ice, Rock, Bryoids layers into one "Landcover class" layer
VRI_combined <- rbind(VRI_clearcut, VRI_forests, VRI_rock, VRI_shrubs, VRI_snow)
write_sf(VRI_combined, "./all_land_VRI.shp")
