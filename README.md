Ambient Fine Particulate Matter and Asthma Emergency Visits in Lima, Peru
Does Residential Proximity to Major Roads Matter?

Overview
This repository contains the R code used for the statistical analyses reported in the article "Ambient fine particulate matter and asthma emergency visits in Lima, Peru: Does Residential Proximity to Major Roads Matter?"

The study investigates the association between short-term ambient fine particulate matter (PM2.5) exposure and asthma emergency department (ED) visits in Lima, Peru, and evaluates whether residential proximity to major roads modifies this association.

The analytical workflow uses a time-stratified case-crossover design and conditional logistic regression with distributed lag non-linear models (DLNMs). The code covers:

Exposure-response functional form selection.
Lag-response specification and sensitivity analyses.
Assessment of non-linearity and spline flexibility.
Primary model fitting.
Effect modification by residential proximity to major roads.
Stratified and double-stratified analyses.
Study Population and Exposure Assessment
The study analyzed asthma ED visits from 10 public hospitals in Lima, Peru, between 2017 and 2023. Daily mean PM2.5, temperature, and relative humidity were assigned to geocoded residences using a high-resolution (1 × 1 km) multi-source exposure model.

Data Availability
The individual-level health data used in this study are not publicly available due to confidentiality and data protection requirements. Access to these data may be requested from the participating hospitals, subject to their authorization, applicable ethical requirements, and institutional data access procedures.

The R code used for the statistical analyses is publicly available in this repository.

Requirements
The analysis was developed in R and uses the following packages:

dlnm
survival
data.table
lmtest
splines (included with R)
The required packages can be installed using:


r
install.packages(c("dlnm", "survival", "data.table", "lmtest"))
Running the Analysis
To run the analysis:

Install R and the required packages.
Obtain authorized access to the study data.
Prepare the dataset according to the variable names and structure expected by the script.
Load the prepared dataset into R as bd_integral.
Run analysis.R.
The analysis requires a preprocessed dataset named bd_integral, containing the outcome, exposure variables, lagged exposure measurements, meteorological covariates, matching strata, and other variables used in the statistical models. The dataset is not included in this repository.

Reproducibility
This repository provides the R code used for the statistical analyses reported in the accompanying article. Reproducing the analyses requires access to the appropriate input data, the necessary computational environment, and compatible versions of R and the required packages.

Citation
Please cite the accompanying scientific article when using this code.
