
# importing modules for this analysis

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
import meddatasets as md
import os # helps check working directory

#-------------------------------------------------------------------------------

# 1. ENVIRONMENT SETUP
# Setting visual style
sns.set_theme(style="whitegrid")
# Equivalent to options(scipen = 999)
pd.set_option('display.float_format', lambda x: '%.4f' % x)

os.getcwd()  # Check current working directory

#-------------------------------------------------------------------------------

# 2. DATA IMPORT & CLEANING
# Equivalent to read.csv and janitor::clean_names
df = pd.read_csv('datasets/diabetes.csv')
df.columns = [col.lower().replace(' ', '_') for col in df.columns]

# Handling missing values (ponds_analysis.R logic)
df_clean = df.dropna().copy()
df_clean
df_clean.describe()

# 3. EXPLORATORY DATA ANALYSIS (EDA)
print(df_clean.info())      # glimpse() equivalent
print(df_clean.describe())  # summary() equivalent

# Visual audit (visdat equivalent)
plt.figure(figsize=(10, 6))
sns.heatmap(df_clean.isnull(), cbar=False, cmap='viridis')
plt.title("Missing Data Map")
plt.show()

# 4. DIAGNOSTIC TESTING (NORMALITY)
# Shapiro-Wilk test
stat, p = stats.shapiro(df_clean['target_column'])
print(f'Shapiro-Wilk p-value: {p:.4f}')

# QQ-Plot
sm.qqplot(df_clean['target_column'], line='s')
plt.show()

# 5. STATISTICAL MODELING

# A. Logistic Regression (glycosylated_haem_logistic_registration.R logic)
# 'Outcome ~ Var1 + Var2'
formula = 'outcome ~ ' + ' + '.join(df_clean.columns.drop('outcome'))
model = smf.logit(formula=formula, data=df_clean).fit()
print(model.summary())

# Odds Ratios and 95% CI
params = model.params
conf = model.conf_int()
conf['OR'] = params
conf.columns = ['Lower CI', 'Upper CI', 'OR']
print(np.exp(conf))

# B. Two-Way ANOVA with Interaction (tumor_size.R logic)
# model_anova = smf.ols('tumor_size ~ C(treatment) * C(gender)', data=df_clean).fit()
# sm.stats.anova_lm(model_anova, typ=1)

# C. Principal Component Analysis (pca_analysis.R logic)
features = df_clean.select_dtypes(include=[np.number])
x = StandardScaler().fit_transform(features)
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
print(f"Explained variance: {pca.explained_variance_ratio_}")

# 6. PERFORMANCE & VISUALIZATION

# ROC Curve
fpr, tpr, thresholds = roc_curve(df_clean['outcome'], model.predict())
roc_auc = auc(fpr, tpr)

plt.plot(fpr, tpr, label=f'ROC curve (area = {roc_auc:.2d})')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.show()

# 7. EXPORTING RESULTS
# df_clean.to_csv('cleaned_data_output.csv', index=False)
