# Diabetes data set

# Overview of the dataset, from dataset description

# Pregnancies: Number of times pregnant
# Glucose: Plasma glucose concentration a 2 hours in an oral glucose tolerance test
# BloodPressure: Diastolic blood pressure (mm Hg)
# SkinThickness: Triceps skin fold thickness (mm)
# Insulin: 2-Hour serum insulin (mu U/ml)
# BMI: Body mass index (weight in kg/(height in m)^2)
# DiabetesPedigreeFunction: Diabetes pedigree function
# Age: Age (years)
# Outcome: Class variable (0 or 1) 
# (class value 1 is interpreted as "tested positive for
# diabetes")

#-------------------------------------------------------------------------------

# loading necessary libraries
library(tidyverse)
library(dlookr)
library(ggstatsplot)
library(visdat)
library(caret) # for generating confusion matrix to access the accuracy of model
library(pROC) # for determining model performance using area under the curvev (AUC)
# Receiver Operating Characteristic Curve
library(gridExtra) # for arranging multiple plots in grids for concise display
library(gtsummary) # for creating publication ready tables
library(performance) # for determining model performance, can also give plots
library(flextable) # for exporting gtsummary tables in word
library(MASS) # I want to use the stepAIC() to pick the best model
#-------------------------------------------------------------------------------
# importbroom# importing the diabetes data set
dm <- read.csv('diabetes.csv')

#-------------------------------------------------------------------------------

# exploring the diabetes data
view(dm)
glimpse(dm)
names(dm)
plot(dm)
spec(dm)


#-------------------------------------------------------------------------------
#m fitting model
?glm
# Fit multiple logistic regression
dm_model <- glm(
  Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness +
    Insulin + BMI + DiabetesPedigreeFunction + Age,
  data = dm,
  family = binomial
)
summary(dm_model)

# Fit multiple logistic regression
model_diabetes <- glm(
  Outcome ~ .,
  data = dm,
  family = binomial
)
summary(model_diabetes)

#-------------------------------------------------------------------------------
# best model using stepAIC()
# both sides
best_dm_model <- stepAIC(dm_model, direction='both')
summary(best_dm_model)

# forward model
forward_dm_model <- stepAIC(dm_model, direction='forward')
summary(forward_dm_model)

# backward model
backward_dm_model <- stepAIC(dm_model, direction='backward')
summary(backward_dm_model)

#-------------------------------------------------------------------------------
# odds ratios with CI
OR_table <- exp(cbind(OR = coef(dm_model), confint(dm_model)))

print(OR_table)

#-------------------------------------------------------------------------------
# TUESDAY CLASS

# Coefficients, SE, and CI
coef_table <- summary(dm_model)$coefficients
ci_table  <- confint(dm_model)

# Combine everything in one data frame
complete_table <- cbind(
  Estimate = coef(dm_model),
  Std_Error = coef_table[, "Std. Error"],
  OR = exp(coef(dm_model)),
  Lower_95_CI = exp(ci_table[, 1]),
  Upper_95_CI = exp(ci_table[, 2]),
  p_value = coef_table[, "Pr(>|z|)"]
)

print(complete_table)

#rounding off to 4 decimals
round(complete_table, 4)

# converting table into data fra,e

complete_table = as.data.frame(complete_table) %>% 
  round(4)

getwd()

write.csv(complete_table, file = 'dm_model_output.csv')

#-------------------------------------------------------------------------------

# model performance
# Use the function performance to view pseudo R2 results performance(dm_model)
performance(dm_model)

# model performance using AUC of ROC
roc_result <- roc(Outcome ~ fitted.values(dm_model), data = dm)

# Basic ROC plot
plot(roc_result, main = "ROC Curve")
auc(roc_result)

check_model(dm_model)

# Generating OR results using gtsummary package
OR_tab <- tbl_regression(dm_model, exponentiate = TRUE)
OR_tab

dm_linear <- glm(
  BMI ~ Glucose + Age,
  data = dm,
  family = lin
)

# dm model as linear regression
dm_linear_one <- lm(
  BMI ~ Glucose + Age,
  data = dm
)

summary(dm_linear_one)
#-------------------------------------------------------------------------------
# docs exploration
?stepAIC
