
# Ambient Fine Particulate Matter and Asthma Emergency Visits in Lima, Peru

**Does Residential Proximity to Major Roads Matter?**

## Overview

This repository contains the R code used for the statistical analyses reported in the article *"Ambient fine particulate matter and asthma emergency visits in Lima, Peru: Does Residential Proximity to Major Roads Matter?"*

The study investigates the association between short-term ambient fine particulate matter (PM2.5) exposure and asthma emergency department (ED) visits in Lima, Peru, and evaluates whether residential proximity to major roads modifies this association.

The analytical workflow uses a time-stratified case-crossover design and conditional logistic regression with distributed lag non-linear models (DLNMs). The code covers:

- Exposure-response functional form selection.
- Lag-response specification and sensitivity analyses.
- Assessment of non-linearity and spline flexibility.
- Primary model fitting.
- Effect modification by residential proximity to major roads.
- Stratified and double-stratified analyses.

## Study Population and Exposure Assessment

The study analyzed asthma ED visits from 10 public hospitals in Lima, Peru, between 2017 and 2023. Daily mean PM2.5, temperature, and relative humidity were assigned to geocoded residences using a high-resolution (1 × 1 km) multi-source exposure model.

## Data Availability

The individual-level health data used in this study are not publicly available due to confidentiality and data protection requirements. Access may be requested from the participating hospitals, subject to their authorization, applicable ethical requirements, and institutional data access procedures.

The R code used for the statistical analyses is publicly available in this repository.

## Requirements

The analyses were conducted in R version 4.3.2 using the following packages:

- `dlnm`
- `survival`
- `data.table`
- `lmtest`
- `splines` (included with R)

Install the required packages using:

```r
install.packages(c("dlnm", "survival", "data.table", "lmtest"))
```

## Running the Analysis

The analysis requires a preprocessed dataset named `bd_integral`, containing the outcome indicator, matching-stratum identifier, PM2.5 exposure and lagged values, meteorological variables and their lagged values, holiday indicator, residential proximity indicators, and variables used for stratified analyses.

To run the analysis:

1. Install R and the required packages.
2. Obtain authorized access to the study data.
3. Prepare the analytical dataset with the variable names and structure required by `analysis.R`.
4. Load the dataset into R as `bd_integral`.
5. Run `analysis.R`.

The dataset is not included in this repository.

## Reproducibility

Reproducing the analyses requires access to the appropriate input data and a compatible R environment with the required packages.

## Citation

Please cite the accompanying scientific article when using this code.
