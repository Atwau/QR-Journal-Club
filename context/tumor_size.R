
# Multivariate analysis
# Week12
# Atwau Pius
# Biostatisitcs class, Biology department, by Dr. Raphael Wangalwa

#-------------------------------------------------------------------------------

# loading necessary packages
library(tidyverse)
library(dplyr)
library(visdat)


#-------------------------------------------------------------------------------
# setting working directory

ts <- read.csv("D:/Medicine/MSc Biochemistry/Sem 1 Yr 1/Biostatistics and Research Methods/Biology department/biostatistics/R_Biostatistical_Analysis_Dr_Wnagalwa/Tumor_Size.csv")

#-------------------------------------------------------------------------------
# Exploring the dataset

vis_dat(ts)
vis_miss(ts)

glimpse(ts)

# transforming the gender and treatment to factor variables

ts_mut <- ts %>% mutate(
  Treatment = factor(Treatment, levels = c("Control", "Drug A", "Drug B")),
  Gender = factor(Gender, levels = c("Men", "Women")),
  Tumor_size = as.numeric(Tumor_size)
)

ts_mut

vis_dat(ts_mut)


#-------------------------------------------------------------------------------
# Performing 2 way anova

# Two-way ANOVA (factorial with interaction)

# Use aov() or lm(); both are fine. We'll keep aov() for simplicity.
ts_aov <- aov(Tumor_size ~ Treatment * Gender, data = ts_mut)
summary(ts_aov)  # Classical ANOVA table (Type I SS)

# We use * to show that these two factors are interacting

# Interaction plot using ggplot
summary_ts <- ts_mut %>%
  group_by(Treatment, Gender) %>%
  summarise(
    mean = mean(Tumor_size, na.rm = TRUE),
  )

# using colored plots...
p <- ggplot(summary_ts, aes(x = Treatment, y = mean, group = Gender, color = Gender)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    title = "Interaction Plot (Mean): Tumor Size by Treatment and Gender",
    x = "Treatment",
    y = "Tumor Size"
  )
print(p)
