# 📌 Research Overview
This research portfolio presents a statistical analysis of 541 polyp cases to differentiate between neoplastic and non-neoplastic lesions. By evaluating clinical predictors such as White Light Endoscopy (WLE) findings, polyp size, and location, this project aims to provide a reproducible workflow for clinical decision-making in colorectal cancer prevention.

The project utilizes an R-based workflow (workflowr) to ensure transparency and reproducibility in medical statistics.

## Key Predictors Analyzed:
- WLE (White Light Endoscopy): Categorized as neoplastic or non-neoplastic.
- Size: Divided into $\le 5\text{mm}$ and $> 5\text{mm}$.
- Morphology (Paris Classification): Focuses on types IIa and Is-Ip.
- Location (Right/Left Colon): To identify anatomical correlations with neoplastic growth.
- Demographics: Age groups categorized by a 45-year threshold.

## Statistical Framework:
- Multivariate Logistic Regression: Used to calculate the probability of neoplastic outcomes.
- Data Partitioning: The dataset was split into Training ($n=378$) and Validation ($n=163$) sets, maintaining a consistent neoplastic ratio (approx. 57%) across groups.
- Decision Curve Analysis (DCA): Implemented to evaluate the clinical net benefit of the model against default strategies.

This analysis was performed using R. Key libraries used include:
- Data Manipulation: tidyverse, dplyr, readxl
- Medical Statistics: gtsummary, rms, pROC
- Visualization: ggplot2, dcurves (for dca)
- Reproducibility: workflowr
