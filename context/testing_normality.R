library(tidyverse) # 
library(janitor)
library(nortest) # contains additional nrmality tests
library(dlookr)
library(visdat)

# bias against scientific notation
# use options(scipen = 999) to restrict display to decimals
# I don't like scientific notation of display
options(scipen = 999)

# loading kidney disease practice dataset withread.csv
# then nesting it in clean_names function from janitor package to automate
# formating of column names to R acceptable format
kd <- clean_names(read.csv('datasets/kidney_disease_dataset.csv'))

#-------------------------------------------------------------------------------

# explorative data analysis
View(kd)
glimpse(kd)

unique(kd$target)

#-------------------------------------------------------------------------------
# Research question:/ hypothsis
# Does age have any influence on risk of getting kidney disease?
# Old age increases the risk of kidney disease
# Old age has no influence on risk of getting kidney disease

# test for normality of age
age_normality <- kd$'age_of_the_patient' %>% 
  normality()

# checking for sanity/complete
vis_dat(kd)
vis_miss(kd)
getwd()

#===============================================================================

#------STATISTICAL TESTS FOR NORMALITY------------------------------------------

# Shapiro-Wilk test
# suitable fr data between 3- 5000 observations
shapiro.test(kd$blood_pressure_mm_hg)


# Kolmogorov–Smirnov test (KS test)
ks.test(kd$blood_pressure_mm_hg, 'pnorm', mean=mean(kd$blood_pressure_mm_hg), 
        sd=sd(kd$blood_pressure_mm_hg))

# Lilliefors test - variant for KS
lillie.test(kd$blood_pressure_mm_hg)

# Anderson Darling test
ad.test(kd$blood_pressure_mm_hg)

#===============================================================================

#---------------Using plots to test for normality-------------------------------

qqnorm(kd$blood_pressure_mm_hg)
qqline(kd$blood_pressure_mm_hg, col='red')

ggplot(kd, aes(x=blood_pressure_mm_hg))+
  labs(title='Density plot of blood pressure to visualize normality')+
  geom_density()
ggplot(kd, aes(x=blood_pressure_mm_hg))+
  geom_histogram()+
  labs(title='Histogram/Frequency plot for blood pressure')

#-------------------------------------------------------------------------------

# I want to analyze for normality using categories of either hypertension or
# diabetes status, to see if it will change the result of my normality tests

# First I have to select variable of interest into a dataframe
# Let me explore the dataset to see what variables I can pick out

str(kd)

# List of variables that I want to focus my analysis

# I will perhaps think about DM, and HTN, then perhaps create data frames which
# will contain variables which seem to be risk factors for DM, HTN respectively

# age, blood pressure, albumin in urine, sugar in urine, random blood sugar levels,
# blood urea, serum creatinine, sodium levels, haemoglobin levels, potassium levels,
# packed cell volumes, hypertension yes_no, diabetes yes_no, coronoary artery
# disease yes_no, eGFR, urine protein to urea ratio, urine output per day, 
# serum albumin levels, cholesterol levels, 


#===============================================================================

#----Section 5.6 r4ds

by_age <- group_by(kd, hypertension_yes_no)%>% print()
view(by_age)

by_htn <- kd %>% select(age_of_the_patient, diabetes_mellitus_yes_no, hypertension_yes_no, 
                       estimated_glomerular_filtration_rate_e_gfr, urine_output_ml_day, 
                       urine_protein_to_creatinine_ratio, blood_pressure_mm_hg) %>% 
  group_by(hypertension_yes_no) %>% summarise(mean_BP=mean(blood_pressure_mm_hg),
                                                   min_BP=min(blood_pressure_mm_hg),
                                                   max_BP=max(blood_pressure_mm_hg)) %>% 
  print()
view(by_htn)
