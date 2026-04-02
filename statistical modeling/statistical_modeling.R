# Statistical modeling
# Dr. Atwau Pius

library(tidyverse)
library(dplyr)

# importing diabetes dataset for modeling

dm <- read.csv("statistical modeling/datasets/diabetic_data.csv")

# exploration 

glimpse(dm)
view(dm)
colnames(dm)
unique(dm$discharge_disposition_id)
