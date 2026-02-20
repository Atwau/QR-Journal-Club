#===============================================================================
# PHASE 1: ENVIRONMENT SETUP & DEPENDENCIES
#===============================================================================

# 1. Load reticulate and configure the environment
if (!require("reticulate")) install.packages("reticulate")
library(reticulate)

# It is best practice to create/use a specific environment for your project
# This prevents version conflicts with other Python projects
virtualenv_create("python_in_R")
use_virtualenv("python_in_R", required = TRUE)

# 2. Bulk install Python dependencies
# meddatasets relies on the data science stack (pandas, numpy, etc.)
py_packages <- c("meddatasets", "pandas", "numpy", "scikit-learn", 
                 "matplotlib", "seaborn", "xgboost")

py_install(py_packages, envname = "python_in_R", pip = TRUE)

# Verify configuration
py_config()

#-------------------------------------------------------------------------------
# All you need to run every time you start R
library(reticulate)
use_virtualenv("python_in_R", required = TRUE)

#===============================================================================
# PHASE 2: PYTHON INTEGRATION (Using meddatasets)
#===============================================================================

# 3. Import the Python module
mds <- import("meddatasets")

# 4. List available Python datasets
all_py_datasets <- mds$get_available_datasets()
print(all_py_datasets)

# 5. Load and convert a Python dataset to an R DataFrame
# Note: Use '$' to access Python methods
cvd <- mds$load_dataset("cardiovascular_risk")

#===============================================================================
# PHASE 3: R ANALYSIS (Tidyverse & Visualization)
#===============================================================================

library(tidyverse)
library(visdat)

# Visualizing data types and missing values
# 
vis_dat(cvd)

# Quick look at structure and summary
glimpse(cvd)
summary(cvd)

#===============================================================================
# PHASE 4: NATIVE R ALTERNATIVE (MedDataSets CRAN)
#===============================================================================

# If you prefer native R datasets without the Python overhead:
if (!require("MedDataSets")) install.packages("MedDataSets")
library(MedDataSets)

# View list of datasets available in the R-native package
data(package = "MedDataSets")

# Example: Loading a specific R dataset from the package
# data("name_of_dataset")