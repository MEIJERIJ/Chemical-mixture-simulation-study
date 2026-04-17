# title publication

> reference to publication


## About this repository
The repository contains R scripts used to implement the various methods discussed, set up the simulation study and generate the results for the case study. The data used is part of the [Teenager HBM Study - 3M site](https://archief.onderzoek.omgeving.vlaanderen.be/Onderzoek-3520355). Due to data protection agreements, this dataset is not publicly available. Nevertheless, the provided code and workflow offer a clear illustration of how the different multi‑pollutant methods were applied and can serve as a guide for reproducing the steps.

## Contents Overview

### [Simulation study](/Simulation%20study/)
This folder contains all information to generate the results of the simulation study. The R scripts for the generate the simulated data are identical for both scenarios.
- [Methods](/Simulation%20study/Method) contains all method discussed and applied on the simulated data.
- [Dense linear data-generating mechanism](/Simulation%20study/Dense%20linear%20data-generating%20mechanism) contains all results and R code to the generate/analyse the dense setting.
- [Sparse linear data-generating mechanism](/Simulation%20study/Sparse%20linear%20data-generating%20mechanism) contains all results and R code to the generate/analyse the sparse setting.

### [Case study](/Case%20study/)
This folder contains the information related to the case study such as the generated visualsions, R code used to do the data cleaning and the R code for the applied [methods](/Case%20study/Method).
