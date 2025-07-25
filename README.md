# Anopheles stephensi Climate Attribution

This repository contains an R-based pipeline to model the historical and recent distribution of *Anopheles stephensi* using boosted regression trees (BRTs), and assess the contribution of climate change using counterfactual climate simulations.

## Overview

The `main.R` script includes:

- Preprocessing of *An. stephensi* presence records and background points
- Extraction of environmental covariates (climate, land use, vegetation productivity)
- Spatial cross-validation and BRT model training
- Temporal projections of ecological suitability (1901–2019)
- Climate attribution using ISIMIP3a observed vs counterfactual data
- Variable importance analysis
- Evaluation with prevalence-adjusted Sørensen index (SI<sub>ppc</sub>)

## Requirements

Install the required R packages:

```r
install.packages(c(
  "raster", "sp", "sf", "dismo", "gbm", "blockCV",
  "stringr", "dplyr", "RColorBrewer", "readxl"
))
```
## Data

Occurrence: data/data_Ansteph2.csv

Background: data/background_not_stephensi.csv

Environmental variables: ISIMIP3a datasets, LUH2 land-use rasters in env_data folder

