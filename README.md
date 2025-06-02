---
## Title of paper associated with this repository: Modelling functional connectivity for bears among spawning salmon waterways in Haíɫzaqv (Heiltsuk) Territory, coastal British Columbia
---

The data and R scripts below were used in the creation of our Landcover Class resistance surface, Validation model, as well as Figures 5 and 6 in 'Modelling functional connectivity for bears among spawning salmon waterways in Haíɫzaqv (Heiltsuk) Territory, coastal British Columbia' 

Please contact the corresponding author, Ilona Mihalik, at ilonammillie@gmail.com if you have any questions.

## Description of the data and code files

Mihalik_SupplementalCode1.R: subsets from the BC Vegetation Resource Inventory Database, and prepares the 'Landcover class' and 'Water' spatial layers for our cumulative resistance surface. Please note - the creation of our cumulative resistance surface included steps within both R (v4.2.2) and ArcMap 10.7.1.

Mihalik_SupplementalCode2.R: reproduces the steps for validating our cumulative resistance layer using an independent genetic dataset of individual bears from the same study area. These steps are outlined in Section 2.6 of our Methods and include: preparing and subsetting from the genetic recapture dataset, the creation of the network graph (nodes and edges) of transits between individual bears, the logistic regression analysis, and the creation of Figures 5 and 6. 

genetic_individuals.csv: these data are called in to the above script (Mihalik_SupplementalCode2.R). They include individual bear (black and grizzly) detections from 2015-2019, and are a subset of the larger independent project dataset. Note: unique numbers replace the real UTME and UTMN coordinates within this subset, as the sampling locations cannot be shared (due to agreements with the Haíɫzaqv Nation).

validation_resistances_3columns.csv: these data are called in to the above script (Mihalik_SupplementalCode2.R). They include the output effective resistances from Circuitscape, using our cumulative resistance surface and hair snag sites as focal nodes. The resistances are converted to effective conductance values to include in the Firth's penalized-likelihood logistic regression model. 

Mihalik_SupplementalCode3.R: reproduces the Least Cost Paths analysis to identify the most important pathways from the Circuitscape cumulative current map.

## Data used from other sources
All original spatial data are publicly available and were downloaded from the BC Data Catalogue (https://catalogue.data.gov.bc.ca) in September 2022. We obtained landcover, water, and forestry harvest year data from the provincial 2021 Vegetation Resources Inventory (VRI; https://catalogue.data.gov.bc.ca/dataset/2ebb35d8-c82f-4a17-9c96-612ac3532d55). We derived the terrain ruggedness index from a digital elevation model (DEM; https://catalogue.data.gov.bc.ca/dataset/digital-elevation-model-for-british-columbia-cded-1-250-000) using the Raster Terrain Analysis plugin in QGIS (following methods from Riley et al. 1999). Finally, spatial data for Pacific salmon spawning zones were obtained from the BC Historical Fish Distribution (https://catalogue.data.gov.bc.ca/dataset/bc-historical-fish-distribution-zones-50-000) spatial layer.

## Code/Software
The analyses here were conducted in R (v4.2.2; R Core Team 2022) and ArcMap 10.7.1 (ESRI 2019). 
Esri (2019) ArcGIS Desktop: Release 10.7.1. 
R Core Team (2022) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
