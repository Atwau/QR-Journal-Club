#===============================================================================
# MASTER DATA ANALYSIS SCRIPT
# Integrating: Cleaning, EDA, Diagnostics, Regression, PCA, and ANOVA
#===============================================================================

# 1. ENVIRONMENT SETUP
library(tidyverse)  # Data manipulation
library(janitor)    # Name cleaning
library(visdat)     # Missing data visualization
library(dlookr)     # Normality and EDA
library(performance)# Model diagnostics
library(gtsummary)  # Publication-ready tables
library(pROC)       # ROC/AUC analysis
library(MASS)       # Stepwise model selection
library(ggbiplot)   # PCA visualization

# Disable scientific notation for cleaner output
options(scipen = 999)

#-------------------------------------------------------------------------------

# 2. DATA IMPORT & CLEANING
# Using kidney disease logic for name formatting
raw_data <- read.csv('your_dataset.csv')
df <- clean_names(raw_data) 

# Handle missing values using the ponds analysis logic
df_clean <- df %>% drop_na()

# Factor transformation logic from tumor size analysis
# Replace 'GroupVar' and 'CategoryVar' with your actual column names
# df_clean <- df_clean %>% mutate(
#   GroupVar = factor(GroupVar),
#   CategoryVar = factor(CategoryVar)
# )

#-------------------------------------------------------------------------------

# 3. EXPLORATORY DATA ANALYSIS (EDA)
glimpse(df_clean)      # Structural overview
vis_dat(df_clean)      # Visualize data types
vis_miss(df_clean)     # Check for remaining gaps

# Descriptive Summary (Gapminder style)
summary_stats <- df_clean %>%
  group_by(across(where(is.factor))) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

#-------------------------------------------------------------------------------

# 4. DIAGNOSTIC TESTING (NORMALITY)
# Quantitative tests
# shapiro.test(df_clean$numeric_column) 

# Visual normality checks
# qqnorm(df_clean$numeric_column)
# qqline(df_clean$numeric_column, col = 'red')

# dlookr normality summary
normality_table <- df_clean %>% normality()

#-------------------------------------------------------------------------------

# 5. STATISTICAL MODELING

# A. Logistic Regression & Stepwise Selection
full_model <- glm(Outcome ~ ., data = df_clean, family = binomial)
best_model <- stepAIC(full_model, direction = 'both', trace = FALSE)
summary(best_model)

# B. Two-Way ANOVA (Factorial with Interaction)
# res_aov <- aov(Numeric_Var ~ Factor1 * Factor2, data = df_clean)
# summary(res_aov)

# C. Principal Component Analysis (PCA)
# Selecting only numeric columns for PCA
numeric_cols <- df_clean %>% select(where(is.numeric))
pca_results <- prcomp(numeric_cols, center = TRUE, scale. = TRUE)
summary(pca_results)

#-------------------------------------------------------------------------------

# 6. PERFORMANCE & VISUALIZATION

# Performance metrics and ROC
model_performance <- performance(best_model)
roc_obj <- roc(df_clean$Outcome ~ fitted(best_model))
plot(roc_obj, main = "ROC Curve")

# PCA Biplot
# ggbiplot(pca_results, ellipse = TRUE, groups = df_clean$GroupVar) + theme_bw()

# Interaction Plot (Tumor Size style)
# ggplot(df_clean, aes(x = Factor1, y = Numeric_Var, color = Factor2, group = Factor2)) +
#   stat_summary(fun = mean, geom = "line") +
#   stat_summary(fun = mean, geom = "point")

#-------------------------------------------------------------------------------

# 7. EXPORTING RESULTS
# Regression Table (OR and CI)
final_table <- tbl_regression(best_model, exponentiate = TRUE) %>% as_gt()

# Save model coefficients
# write.csv(as.data.frame(summary(best_model)$coefficients), "model_results.csv")