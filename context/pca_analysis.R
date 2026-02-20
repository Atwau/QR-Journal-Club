# Analaysis of ponds data 
# principle component analysis
# MONDAY 17 11 2025

#-------------------------------------------------------------------------------
#loading necessary packages
library(devtools)
library(ggbiplot)
library(tidyverse)
library(visdat)


# ---Principal Component Analysis---
# 1. Use ponds data
# loading ponds data

pd <- read.csv('Ponds.csv')
view(pd)
vis_miss(pd)+
  theme(plot.margin = margin(t=60, r=30, b=10, l=10))

# 2. Run PCA on numeric columns only
# scale. = TRUE standardizes variables (important since they're in cm)
pd_clean <- pd %>% drop_na(3:9) # removing missing values



ponds_pca <- prcomp(pd_clean[, 3:9],
                   center = TRUE,
                   scale. = TRUE) # scaling is standardizing the data

# 3. View PCA results
print(ponds_pca)          # Basic information about PCA, the correlations/loadings
summary(ponds_pca)        # Proportion of variance explained, key in interpretation
ponds_pca$x               # Obtaining Principal component actual scores

# 4. Plot the PCA results (base R)
biplot(ponds_pca)         # Creates a biplot showing variables and observations


# plotting the PCA biplot using ggbiplot
ggbiplot(ponds_pca, obs.scale = 1, var.scale = 1,
         groups = pd_clean$Station, point.size = 2,
         varname.size = 5,
         varname.color = 'blue',
         varname.adjust = 1.2,
         ellipse = F,
         circle = T, # shows the central point of the plot, may not be necessary
         ellipse.prob = 0.80,
         ellipse.linewidth = 1,
         ellipse.fill = TRUE,
         ellipse.alpha = 0.25) +
  labs(fill = "Station", color = "Station") +
  theme_bw() +
  theme(legend.direction = 'vertical', legend.position = 'right')
