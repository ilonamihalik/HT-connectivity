#############################################
# SCRIPT 2 - Modelling functional connectivity for bears among spawning salmon waterways in Haíɫzaqv (Heiltsuk) Territory, coastal British Columbia

# This code reproduces the steps for validating our cumulative resistance layer using an 
# independent genetic dataset of individual bears from the same study area. These steps include: preparing and subsetting from the genetic recapture dataset, the creation of the network graph (nodes and edges)
# of transits between individual bears, the logistic regression analysis, and the creation of Figures 5 and 6 within the manuscript.
# Note: the genetic recapture dataset provided includes unique numbers to replace the real UTME and UTMN coordinates, as the sampling locations cannot be provided (due to agreements with the Haíɫzaqv Nation). 
#############################################

library(dplyr)
library(sf)
library(igraph)
library(ggeffects)
library(marginaleffects)
library(emmeans)
library(envalysis)
library(ggplot2)
library(brglm2)
library(gtsummary)
library(finalfit)
library(MASS)

setwd("")


### 1) Subset from genetic recapture dataset ------------------------------------------
# year = sampling year
# site_id = sampling year and site
# revisitid = sampling year and site and revisit number... sites were ‘revisited’ twice, approximately 10-14 days apart, within a sampling season

allsamples <- read.csv("./genetic_individuals.csv",  fileEncoding='latin1', stringsAsFactors=FALSE)

allsamples$site_id <- gsub("2015", "", allsamples$site_id)
allsamples$site_id <- gsub("2016", "", allsamples$site_id)
allsamples$site_id <- gsub("2017", "", allsamples$site_id)
allsamples$site_id <- gsub("2018", "", allsamples$site_id)
allsamples$site_id <- gsub("2019", "", allsamples$site_id)

##Subsets for each year
HT2015 <- allsamples %>%
  filter(year == "2015")
HT2015 <- transform(HT2015[grep("^\\d+$", HT2015$site_id),,drop=F], site_id= as.numeric(as.character(site_id)))

HT2016 <- allsamples %>%
  filter(year == "2016") 
HT2016 <- transform(HT2016[grep("^\\d+$", HT2016$site_id),,drop=F], site_id= as.numeric(as.character(site_id)))

HT2017 <- allsamples %>%
  filter(year == "2017") 
HT2017 <- transform(HT2017[grep("^\\d+$", HT2017$site_id),,drop=F], site_id= as.numeric(as.character(site_id)))

HT2018 <- allsamples %>%
  filter(year == "2018") 
HT2018 <- transform(HT2018[grep("^\\d+$", HT2018$site_id),,drop=F], site_id= as.numeric(as.character(site_id)))

HT2019 <- allsamples %>%
  filter(year == "2019") 
HT2019 <- transform(HT2019[grep("^\\d+$", HT2019$site_id),,drop=F], site_id= as.numeric(as.character(site_id)))

allsamples <- transform(allsamples[grep("^\\d+$", allsamples$site_id),,drop=F], site_id= as.numeric(as.character(site_id)))


##Determine unique revisits within each year
# 2015
unique2015 <- HT2015 %>%
  distinct(individual, revisitid, .keep_all = TRUE)

unique2015IND <- HT2015 %>%
  distinct(individual)

# 2016
unique2016 <- HT2016 %>%
  distinct(individual, revisitid, .keep_all = TRUE)

unique2016IND <- HT2016 %>%
  distinct(individual)

# 2017
unique2017 <- HT2017 %>%
  distinct(individual, revisitid, .keep_all = TRUE)

unique2017IND <- HT2017 %>%
  distinct(individual)


# 2018
unique2018 <- HT2018 %>%
  distinct(individual, revisitid, .keep_all = TRUE)

unique2018IND <- HT2018 %>%
  distinct(individual)

# 2019
unique2019 <- HT2019 %>%
  distinct(individual, revisitid, .keep_all = TRUE)

unique2019IND <- HT2019 %>%
  distinct(individual)

#By revisits
unique.ind.revisits2015 <- unique(unique2015[c("individual","site_id")])
unique.ind.revisits2016 <- unique(unique2016[c("individual","site_id")])
unique.ind.revisits2017 <- unique(unique2017[c("individual","site_id")])
unique.ind.revisits2018 <- unique(unique2018[c("individual","site_id")])
unique.ind.revisits2019 <- unique(unique2019[c("individual","site_id")])

### 2) Create edge lists for unique transits based on revisit info ------------------------------------------
#2015
edgelist2015 <- unique.ind.revisits2015 %>% group_by(individual) %>%
  filter(n()>=2) %>% group_by(individual) %>%
  do(data.frame(t(combn(.$site_id, 2)), stringsAsFactors=FALSE)) %>%
  rename(to = X1, from = X2, type = individual)

edgelist2015$from <- gsub("2015", "", edgelist2015$from)
edgelist2015$to <- gsub("2015", "", edgelist2015$to)
edgelist2015 <- edgelist2015[c("to", "from", "type")]

#2016
edgelist2016 <- unique.ind.revisits2016 %>% group_by(individual) %>%
  filter(n()>=2) %>% group_by(individual) %>%
  do(data.frame(t(combn(.$site_id, 2)), stringsAsFactors=FALSE)) %>%
  rename(to = X1, from = X2, type = individual)

edgelist2016$from <- gsub("2016", "", edgelist2016$from)
edgelist2016$to <- gsub("2016", "", edgelist2016$to)
edgelist2016 <- edgelist2016[c("to", "from", "type")]

#2017
edgelist2017 <- unique.ind.revisits2017 %>% group_by(individual) %>%
  filter(n()>=2) %>% group_by(individual) %>%
  do(data.frame(t(combn(.$site_id, 2)), stringsAsFactors=FALSE)) %>%
  rename(to = X1, from = X2, type = individual)


edgelist2017$from <- gsub("2017", "", edgelist2017$from)
edgelist2017$to <- gsub("2017", "", edgelist2017$to)
edgelist2017 <- edgelist2017[c("to", "from", "type")]

#2018
edgelist2018 <- unique.ind.revisits2018 %>% group_by(individual) %>%
  filter(n()>=2) %>% group_by(individual) %>%
  do(data.frame(t(combn(.$site_id, 2)), stringsAsFactors=FALSE)) %>%
  rename(to = X1, from = X2, type = individual)


edgelist2018$from <- gsub("2018", "", edgelist2018$from)
edgelist2018$to <- gsub("2018", "", edgelist2018$to)
edgelist2018 <- edgelist2018[c("to", "from", "type")]

#2019
edgelist2019 <- unique.ind.revisits2019 %>% group_by(individual) %>%
  filter(n()>=2) %>% group_by(individual) %>%
  do(data.frame(t(combn(.$site_id, 2)), stringsAsFactors=FALSE)) %>%
  rename(to = X1, from = X2, type = individual)

edgelist2019$from <- gsub("2019", "", edgelist2019$from)
edgelist2019$to <- gsub("2019", "", edgelist2019$to)
edgelist2019 <- edgelist2019[c("to", "from", "type")]

edgelist20152019 <- rbind(edgelist2015, edgelist2016, edgelist2017, edgelist2018, edgelist2019)
edgelist20152019$to <- edgelist20152019$to <- sub("^0+", "", edgelist20152019$to)
edgelist20152019$from <- edgelist20152019$from <- sub("^0+", "", edgelist20152019$from)


### 3) Create nodes list of unique sites ------------------------------------------
uniquesitelist20152019 <- dplyr::select(allsamples, -individual, -year) %>%
  dplyr::distinct(site_id, utme, utmn, utm_zone) 

uniquesitelist20152019 <- uniquesitelist20152019[!duplicated(uniquesitelist20152019$site_id),]


### 4) Create graph of nodes and edges --------------------------------------------
G <- graph_from_data_frame(edgelist20152019, directed = TRUE, vertices = uniquesitelist20152019)
E(G)$weight <- 1
G_2 <- igraph::simplify(G, remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type = "ignore") )

##All possible paths (edges) between nodes (transited and untransited)
allnodepairs <- data.frame(t(combn(uniquesitelist20152019$site_id, m = 2)))

## "Present" or "absent" transit from hair data
allnodepairs$transited <- c(0, 1)[mapply(function(from, to)
  any(G_2df$from == from & G_2df$to == to), allnodepairs$X1, allnodepairs$X2) + 1]

allnodepairs <- allnodepairs %>%
  dplyr::rename(From = X1) %>%
  dplyr::rename(To = X2) 

## Adding columns where To/From are switched as well
allnodepairs$From2 <- allnodepairs$To
allnodepairs$To2 <- allnodepairs$From

### Connect to Resistance layer output from Circuitscape (ran using hair snag sites as nodes)
##Bring in output resistance file (3 columns)
validresistances <- read.csv("./validation_resistances_3columns.csv")

##Create conductance column (inverse of resistance)
validresistances$conductance <- lapply(validresistances$Resistance, function(x) 1/x)
validresistances$conductance <- as.numeric(validresistances$conductance)

##Join to graph present/absent data
allnodepairs_conductance <- merge(validresistances, allnodepairs[, c("From", "To", "transited")], by.x=c("From", "To"), by.y=c("From", "To"), all.x = T )
allnodepairs_conductance <- subset(allnodepairs_conductance, select = -Resistance)

##Add in opposite connections as well
allnodepairs_conductance2 <- merge(allnodepairs, validresistances[, c("From", "To", "conductance")], by.x=c("From2", "To2"), by.y=c("From", "To"), all.x = T )
allnodepairs_conductance2 <- subset(allnodepairs_conductance2, select = -c(From2, To2))
allnodepairs_conductance2 <- allnodepairs_conductance2[c("From", "To", "conductance", "transited")]

allnodepairs_conductance_final <- rbind(allnodepairs_conductance, allnodepairs_conductance2) #Merge frontwards and backwards files
allnodepairs_conductance_final <- allnodepairs_conductance_final[complete.cases(allnodepairs_conductance_final), ] #remove NAs

allnodepairs_conductance <- allnodepairs_conductance_final #Rename to fit code below

##Join UTM data to spatialize
allnodepairs_conductance <- allnodepairs_conductance %>%
  left_join(dplyr::select(uniquesitelist20152019, utme, utmn, utm_zone, site_id), by = c("From" = "site_id"))

allnodepairs_conductance <- allnodepairs_conductance %>%
  dplyr::rename(From_utme = utme) %>%
  dplyr::rename(From_utmn = utmn) %>%
  dplyr::rename(From_utm_zone = utm_zone)

allnodepairs_conductance <- allnodepairs_conductance %>%
  left_join(dplyr::select(uniquesitelist20152019, utme, utmn, utm_zone, site_id), by = c("To" = "site_id"))

allnodepairs_conductance <- allnodepairs_conductance %>%
  dplyr::rename(To_utme = utme) %>%
  dplyr::rename(To_utmn = utmn) %>%
  dplyr::rename(To_utm_zone = utm_zone)


### 5) Calculate distances between hair snag sites to include as variable --------------------------------------------
node_to_sf <- st_as_sf(allnodepairs_conductance, coords = c("To_utme", "To_utmn"), crs = st_crs(3005))
node_from_sf <- st_as_sf(allnodepairs_conductance, coords = c("From_utme", "From_utmn"), crs = st_crs(3005))
node_distances <- st_distance(node_to_sf, node_from_sf, by_element = TRUE)
node_distances <- units::set_units(node_distances, "km")

##Join to allnodepairs_conductance
allnodepairs_conductance$distance_km <- node_distances
allnodepairs_conductance$distance_m <- units::set_units(node_distances, "m")

allnodepairs_conductance$distance_m <- units::drop_units(allnodepairs_conductance$distance_m)

# Transform "transited" into a factor and setting labels with factor()
allnodepairs_conductance$transited <- factor(allnodepairs_conductance$transited,
                                 levels = c(0,1),
                                 labels = c("no transit", "transit")
)


### 6) MODEL Firth's biased-reduced logistic regression --------------------------------------------
##Center and scale conductance and distance
allnodepairs_conductance$conductance_s <- (allnodepairs_conductance$conductance - mean(allnodepairs_conductance$conductance))/sd(allnodepairs_conductance$conductance)
allnodepairs_conductance$distance_m_s <- (allnodepairs_conductance$distance_m - mean(allnodepairs_conductance$distance_m))/sd(allnodepairs_conductance$distance_m)

##Write model - probability of transit by effective conductance and distance
logmodel3 <- glm(transited ~ conductance_s + distance_m_s,
                 data = allnodepairs_conductance,
                 family = binomial(link = "logit"),
                 method = "brglmFit"
)

##Summary of output
summary(logmodel3)

##Odds Ratio (OR) with P values
logmodel3_OR <- (tbl_regression(logmodel3, exponentiate = TRUE,
                                estimate_fun = purrr::partial(style_ratio, digits = 5),
                                pvalue_fun = purrr::partial(style_sigfig, digits = 3)))


### Figure 5: Marginal effects plot -------------------------------------------------------------
# Marginal Effects plot 
# Holding at coastal male black bear home range size (from Hatler et al. 2008, and also used by Service et al. 2019 paper) -- 84km2 radius
#84*1000000
#sqrt(8.4e+07/pi) = 5170.883m (5.17km)
#scaled distance value at 5169m = -1.619577

##Predicted response
logmodel3_pred_avg <- predict_response(logmodel3,
                                       terms = c("conductance_s [all]", "distance_m_s [-1.619577]"))
##Unscale conductance
sd(allnodepairs_conductance$conductance) #0.05415069
mean(allnodepairs_conductance$conductance) #0.3830321
logmodel3_pred_avg$x <- (logmodel3_pred_avg$x*0.05415069) + 0.3830321

##Data for plotting
logmodel3_data <- logmodel3$model
logmodel3_data$conductance_s <- (logmodel3_data$conductance_s*0.05415069) + 0.3830321

##Plot model using ggplot
logmodel3_plot_avg <- ggplot(logmodel3_pred_avg,
                             aes(x = x, y = predicted)) +
  geom_line() +
  xlab("Effective conductance (Ĝ)") + ylab("Probability of transit") +
  ylim(0, 1) +
  geom_line(size=1, colour = "gray45") +
  geom_hline(yintercept=1.00, color = "gray10") +
  geom_hline(yintercept=0.50, linetype="dashed", color = "gray45") +
  geom_hline(yintercept=0.00, color = "gray10") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, color = NULL), alpha = .15) +
  geom_point(data = logmodel3_data, aes(x = conductance_s, y = transited),
             size = 3, shape = 21, alpha = 0.5) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.position = "none")


### Figure 6: Odds ratio and confidence intervals -------------------------------------------------------------
or_CI <- round(exp(cbind(coef(logmodel3), confint(logmodel3))), digits=3) %>% 
  as.data.frame()
or_CI <- or_CI %>% 
  mutate(variable=rownames(or_CI)) # extract the variables from rownames

or_CI <- rename(or_CI, c("AOR" = "V1",
                         "lower_bound" = `2.5 %`,
                         "upper_bound" = `97.5 %`))
##Reorder variables
col_order <- c("variable", "AOR", "lower_bound", "upper_bound")
or_CI <- or_CI[, col_order] 

##Plot
plot_logit_model <- or_CI[-1,] %>% #remove row number 1 (The intercept) 
  ggplot(aes(x = reorder(variable, AOR), y = AOR)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_bound, 
                    ymax  = upper_bound),
                width = 0.1, 
                color = "gray50") +
  geom_point(size = 3.5, color = "blue") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  xlab("") + ylab("Odds ratio (95% CI, log scale)") +
  coord_flip(ylim = c(0, 3)) + 
  theme(axis.title = element_text(size = 14)) + 
  theme(axis.text = element_text(size = 14)) +
  scale_x_discrete(breaks=c("conductance_s", "distance_m_s"), labels = c("Conductance", "Distance"))
plot_logit_model

