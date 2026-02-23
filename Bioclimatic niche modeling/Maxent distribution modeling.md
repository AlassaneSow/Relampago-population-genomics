Loading and cropping WorldClim variables to just contain the SE US
```r
library(terra)
tiffiles <- list.files("S:/Smith Lab/Alassane/Main/Relampago Blight/WorldClim/wc2.1_2.5m_bio/", pattern = ".tif$", full.names = TRUE)
bioclim <-rast(tiffiles)
se_extent <- ext(-91,-75,25,37)
bioclim_se <- crop(bioclim, se_extent)
```
Removing correlated bioclim variables
```r
library(terra)
library(corrplot)
library(caret)
set.seed(12)
sample_bioclim <- spatSample(bioclim_se, size=10000,
                             method= "random", 
                             na.rm = TRUE, 
                             as.df = TRUE)
cor_matrix <- cor(sample_bioclim, method = "pearson")
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.cex = 0.7, tl.col = "black")
highCorVars <- findCorrelation(cor_matrix, cutoff = 0.7)
bioclim_reduced <- bioclim_se[[-highCorVars]]
removed <- as.data.frame(names(bioclim_se)[highCorVars]) 
names(bioclim_reduced)
writeRaster(bioclim_reduced, 
            filename = "bioclim_reduced.tif", 
            overwrite = TRUE)
```
These variables had a spearmens correlation coefficiet > 0.7 and were removed: 
wc2.1_2.5m_bio_1 - Annual Mean Temperature
wc2.1_2.5m_bio_10 - Mean Temperature of Warmest Quarter
wc2.1_2.5m_bio_11 - Mean Temperature of the Coldest Quarter
wc2.1_2.5m_bio_14 - Percipitation of Driest Month
wc2.1_2.5m_bio_15 - Percipitation Seasonality
wc2.1_2.5m_bio_16 - Percipitation of Wettest Quarter
wc2.1_2.5m_bio_17 - Percipitation of Driest Quarter
wc2.1_2.5m_bio_18 - Precipitation of Warmest Quarter
wc2.1_2.5m_bio_4 - Temperature Seasonality
wc2.1_2.5m_bio_6 - Min Temperature of Coldest Month
wc2.1_2.5m_bio_7 - Temperature Annual Range (BIO5-BIO6)

We are left with eight variables:
wc2.1_2.5m_bio_12 - Annual Precipitation
wc2.1_2.5m_bio_13 - Precipitation of Wettest Month
wc2.1_2.5m_bio_19 -  Precipitation of Coldest Quarter
wc2.1_2.5m_bio_2 - Mean Diurnal Range (Mean of monthly (max temp - min temp))
wc2.1_2.5m_bio_3 - Isothermality 
wc2.1_2.5m_bio_5 -  Max Temperature of Warmest Month
wc2.1_2.5m_bio_8 - Mean Temperature of Wettest Quarter
wc2.1_2.5m_bio_9 - Mean Temperature of Driest Quarter
Now we can use maxent
```r
setwd("S:/Smith Lab/Alassane/Main/Relampago Blight/WorldClim/wc2.1_2.5m_bio/")
library(maxnet)
library(terra)
library(raster)
occ <- read.csv("S:/Smith Lab/Alassane/Main/Relampago Blight/WorldClim/Parvodontia_occurences.csv")
occ <- occ[, colSums(is.na(occ)) == 0]

# Extract presence values
pres_vals <- terra::extract(bioclim_reduced, occ[, c("longitude","latitude")])

# Sample background values
bg_points <- terra::spatSample(bioclim_reduced, size=10000, method="random", xy=TRUE)
bg_vals <- terra::extract(bioclim_reduced, bg_points[, c("x","y")])

# Create response vector
y <- c(rep(1, nrow(pres_vals)), rep(0, nrow(bg_vals)))

# Combine predictor values
X <- rbind(pres_vals, bg_vals)
rows_complete <- complete.cases(X)
X_clean <- X[rows_complete, ]
X_clean <- X_clean[, !colnames(X_clean) %in% "ID"]
y_clean <- y[rows_complete]
# Fit MaxEnt
model <- maxnet(y_clean, X_clean)
# Predict Habitat Suitability
bioclim_for_predict <- bioclim_reduced[[colnames(X_clean)]]
prediction <- predict(bioclim_for_predict, model, type="cloglog", na.rm=TRUE)
plot(prediction, main="MaxEnt Habitat Suitability (SE US)")

pred_df <- as.data.frame(prediction, xy = TRUE)
colnames(pred_df) <- c("x", "y", "suitability")
head(pred_df)

ggplot(pred_df, aes(x = x, y = y, fill = suitability)) +
  geom_raster() +
  geom_sf(data = us_states, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  coord_sf(xlim = c(-90, -78), ylim = c(25, 35), expand = FALSE) +
  scale_fill_viridis_c(option = "E") +
  scale_size_continuous(limits = c(0,1))+
  labs(title = "MaxEnt Habitat Suitability",
  fill = "Suitability") +
  theme_minimal() -> fullsample


##########All the above is great but it is clearly skewed by the collections near Gainesville
ggplot(occ, aes(y = latitude, x = longitude)) +
  geom_point(data = occ) +
  geom_sf(data = us_states, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  coord_sf(xlim = c(-90, -78), ylim = c(25, 35), expand = FALSE) +
  labs(title = "Parvodontia relampaga collections")+
  theme_minimal()
##Below I randomly sample 100 of these points and then re-do the maxent analysis
subsample <- occ[sample(nrow(occ), 100), ]
#plot subsample
ggplot(subsample, aes(y = latitude, x = longitude)) +
  geom_point(subsample = occ) +
  geom_sf(data = us_states, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  coord_sf(xlim = c(-90, -78), ylim = c(25, 35), expand = FALSE) +
  labs(title = "Parvodontia relampaga collections")+
  theme_minimal()
#Re-do analysis and data prep
sub_pres_vals <- terra::extract(bioclim_reduced, subsample[, c("longitude","latitude")])

# Create response vector
y <- c(rep(1, nrow(sub_pres_vals)), rep(0, nrow(bg_vals)))

# Combine predictor values
X <- rbind(sub_pres_vals, bg_vals)
rows_complete <- complete.cases(X)
X_clean <- X[rows_complete, ]
X_clean <- X_clean[, !colnames(X_clean) %in% "ID"]
y_clean <- y[rows_complete]
# Fit MaxEnt
model <- maxnet(y_clean, X_clean)
# Predict Habitat Suitability
bioclim_for_predict <- bioclim_reduced[[colnames(X_clean)]]
prediction <- predict(bioclim_for_predict, model, type="cloglog", na.rm=TRUE)
plot(prediction, main="MaxEnt Habitat Suitability (SE US)")

pred_df <- as.data.frame(prediction, xy = TRUE)
colnames(pred_df) <- c("x", "y", "suitability")

ggplot(pred_df, aes(x = x, y = y, fill = suitability)) +
  geom_raster() +
  geom_sf(data = us_states, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  coord_sf(xlim = c(-90, -78), ylim = c(25, 35), expand = FALSE) +
  scale_fill_viridis_c(option = "E") +
  scale_size_continuous(limits = c(0,1))+
  labs(title = "MaxEnt Habitat Suitability (100 Random Samples) ",
  fill = "Suitability") +
  theme_minimal() -> subsample

# looks about the same!
```
