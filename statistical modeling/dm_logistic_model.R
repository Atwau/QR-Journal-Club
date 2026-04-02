# Readmission logistic regression model
# Dr. Atwau Pius

# --- 1. Load Libraries and Data ---
library(dplyr)
library(stats)

# Load the dataset, treating '?' as NA
dm <- read.csv("statistical modeling/datasets/diabetic_data.csv", na.strings = "?", stringsAsFactors = FALSE)

# --- 2. Data Preprocessing ---
# Create a binary target variable: 1 if readmitted (<30 or >30), 0 if NO
dm$readmitted_binary <- ifelse(dm$readmitted == "NO", 0, 1)
view(dm)

# Select relevant features and remove rows with missing values for simplicity
modeling_df <- dm %>%
  select(readmitted_binary, time_in_hospital, num_lab_procedures, 
         num_medications, number_diagnoses, gender, diabetesMed) %>%
  filter(!is.na(gender)) %>%
  mutate(
    gender = as.factor(gender),
    diabetesMed = as.factor(diabetesMed)
  )

# --- 3. Fit Logistic Regression Model ---
# We model the probability of readmission based on hospital stay and medication usage
model <- glm(readmitted_binary ~ time_in_hospital + num_lab_procedures + 
               num_medications + number_diagnoses + gender + diabetesMed, 
             data = modeling_df, family = binomial)

# --- 4. Model Results ---
# Display the summary of coefficients and significance levels
summary(model)

# --- 5. Model Evaluation (Optional: Accuracy) ---
# Predict probabilities on the same dataset
modeling_df$predicted_prob <- predict(model, type = "response")
modeling_df$predicted_class <- ifelse(modeling_df$predicted_prob > 0.5, 1, 0)

# Calculate Accuracy
accuracy <- mean(modeling_df$predicted_class == modeling_df$readmitted_binary)
cat("Model Accuracy:", round(accuracy, 4), "\n")

# --- 6. Calculate Odds Ratios ---
# Convert log-odds coefficients to odds ratios for direct interpretation
odds_ratios <- exp(coef(model))
print("Odds Ratios:")
print(odds_ratios)
