## R Code for the simulation study

This directory contains all the R Markdown files and/or R scripts necessary to generate your simulations and save the results. You should also store any functions you write for analysis here. Be sure to provide an overview of each file in the table below.

* * *

| File | Description |
|:-------|:-------------------|
| Prep work | 
|`Error_Term_Simulation.Rmd` | Muyang's code to simulate from LME models from December 2020. |
| Main Functions | 
| `stb.R`   | Script with simultaneous CIs for Q-Q plots|
| `mahalanobis_resid.R`   | Script with Mahalanobis residuals for LMEs|
| `lesaffre_verbeke.R` | Script with `lesaffre_verbeke()` command, diagnostic from Singer et al.|
| `Helper_Functions.R` | Store all helper functions used|
| `Simulation_Functions.R` | Store all simulation functions used|
| Our simualtion codes | 
| `Error_normality_Yucheng.Rmd` | Analysis of residuals from models with nonnormality|
| `Error_linearity_Yicheng.rmd` | Analysis of residuals from models with nonlinearity|
| `Error_heteroscedasticity_Yicheng.rmd` | Analysis of residuals from models with heteroscedasticity|
| `Error_Omitting_Fixed_Yicheng.rmd` | Analysis of residuals from models with missing fixed effect|
| `RE_normality_Yucheng.Rmd` | Analysis of random effects from models with nonnormality|
| `Error_hetero_and_non-normal_sim_Yucheng.Rmd` | Combined Scenario #1 Hetero and non-normal|
| `Error_hetero_nonlinear_Yicheng.rmd` | Combined Scenario #2 Hetero and non-linear|
| `datasets_different_ICC&confounding.Rmd` | Store datasets with different ICC and confounding specified in the comments|
| Simulations Overall | 
| `simulation big-picture sketch.md` | Simulation Settings|
| `foreach_template.Rmd` | A template for using foreach and design_matrix|
| `main.R` | Execuate simualtions as a whole|

* * *

#### List of Essential Packages 

library(lmtest)

library(dplyr)

library(tidyr)

library(lme4)

library(HLMdiag)

library(nlme)

library(PearsonDS)

library(MASS)

library(doParallel)

library(foreach)

library(purrr)

library(readr)

Not essential now but for future plotting: library(STB); library(ggplot2); library(gridExtra); library(qqplotr); 

