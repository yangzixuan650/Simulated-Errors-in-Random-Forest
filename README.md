SIMERF: Simulated Errors Random Forest

## Overview

SIMERF (Simulated Errors Random Forest) is a novel machine learning framework developed as part of my PhD research in Statistics at the University of Western Ontario. The project addresses a common but often overlooked challenge in environmental prediction problems: measurement error in predictor variables.

Traditional machine learning algorithms generally assume that predictor variables are measured without error. However, environmental and geospatial datasets frequently contain uncertainty arising from sensor inaccuracies, interpolation methods, remote sensing products, and data aggregation processes. Ignoring these errors can lead to biased predictions and reduced model reliability.

SIMERF introduces a simulation-based framework that explicitly incorporates predictor uncertainty into the Random Forest learning process, leading to more robust predictive performance when measurement errors are present.

The methodology was developed and evaluated for wildfire occurrence prediction in northwestern Ontario, Canada.

## Research Motivation

Measurement error is pervasive in environmental and ecological datasets. Examples include:

- Remote sensing products with retrieval uncertainty
- Weather observations with instrumentation errors
- Interpolated climate surfaces
- Spatially aggregated environmental variables
- Missing data imputation procedures

While measurement error has been extensively studied in classical statistical models, relatively little work has focused on its impact within modern machine learning frameworks.

This project aims to bridge that gap by developing machine learning methods that explicitly account for uncertainty in predictor variables while maintaining the flexibility and predictive power of Random Forests.


## Methodology

The repository contains implementations of:

1. SIMERF (Simulated Errors Random Forest)

A Random Forest framework that incorporates predictor uncertainty through repeated simulation of measurement errors in the input variables during model training and prediction.

2. Monte Carlo Aggregated Random Forest (MC-RF)

A benchmark approach that uses Monte Carlo simulations of uncertain predictors and aggregates predictions across multiple Random Forest realizations.


## Study Area

The methodology was developed for wildfire occurrence prediction in two fire management regions of northwestern Ontario, Canada:

- Dryden Fire Management Area
Located in northwestern Ontario
Characterized by boreal forest ecosystems
Historically affected by wildfire activity
- Fort Frances Fire Management Area
Adjacent to the Dryden region
Contains diverse forest and climatic conditions
Important operational area for wildfire management



## Interactive Visualization

An interactive R Shiny application has been developed to visualize wildfire occurrence forecasts and model outputs.

R Shiny App: https://zixuanyang.shinyapps.io/DisplayMaps/<img width="802" height="80" alt="image" src="https://github.com/user-attachments/assets/16625da6-b4e0-43f7-90ff-2c68feef876a" />






