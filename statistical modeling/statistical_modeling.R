# Statistical modeling
# Dr. Atwau Pius

library(tidyverse)
library(dplyr)

# importing diabetes dataset for modeling

dm <- read.csv("statistical modeling/datasets/diabetic_data.csv")
dm_ids_map <- read.csv("statistical modeling/datasets/IDS_mapping.csv")

# exploration 

glimpse(dm)
view(dm)
colnames(dm)
unique(dm$discharge_disposition_id)

# exploration IDS_mapping
glimpse(dm_ids_map)
view(dm_ids_map)
