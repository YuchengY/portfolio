**Assumptions:**

|                    | Error term        | Random Effects                | 
|--------------------|-------------------|-------------------------------|
|  Normality         | Shapiro wilk test | Mahalanobis distance qq-plot  |
| Heteroscedasticity | Bp test  | NA                            |
|    Linearity       | ANOVA             | NA                            | 

Potential things to explore: heteroscedasticity+non-normal error term, omitting fixed effect terms etc.

**Scenarios:**

| Error_Normality                                   | Error_Heteroscedasticity                   | Error_Linearity    | RE_normality                               | Hetero non_normal | Hetero non_linear |Missing Fixed Effect |
|---------------------------------------------------|---------------------------------------------|--------------------|--------------------------------------------|-------------------| -------------------|-------------------|
| Extremely skewed: skewness=3, kurtosis = 11       | Heteroscedasticity factor = 2 | Squared            | Intercept moderately skewed, slope normal  | Hetero2+skewed  | Squared Z_ij; hfactor=2             | Reduced (Missing Z_ij)          |
| Moderately skewed: skewness=1.5, kurtosis = 6     | Heteroscedasticity factor = 4  | NA                 | Intercept normal, slope moderately skewed  | Hetero4+skewed                |Squared Z_ij; hfactor=4             |NA                |
| Slightly skewed: skewness=0.8, kurtosis = 4       | Heteroscedasticity factor = 8 | NA               | Intercept and slope both moderately skewed | Hetero8+skewed             |Squared Z_ij; hfactor=8             |NA                |
| Bimodal: means shift up/down by 1, proportion 40% | NA  | NA                 | NA                                   | NA                |NA                |NA                |
| Reference: normal distribution                    | Reference: constant variance (factor = 0)                | Reference: linear  | Reference: bivariate normal                | Reference: Hetero0, Hetero2, Hetero4, Hetero8               | Reference: Linear + hfactor = 0| Reference: Full model|


**Datasest:**

|                                   | Social science datasets                                                                                               | Longitudinal datasets                                                           |
|-----------------------------------|-----------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
| Number of groups                  | Up to 50, try datasets with small number of groups                                                                    | 20-30                                                                           |
| Number of observations per group  | Same sized: 25/group; Balanced: 20-30/group; Unbalanced: 2-50/group                                                                  | 5-8, Bottom 4 obs per group, have groups with omitted obs (2 or 3)              |
| ICC                               | High: 0.9; Middle: 0.5; Low: 0.1                                                                                      | NA                                                                              |
| confounding                       | High: error variance 4, RE variance 1; Middle: error variance 4, RE variance 4; Low: error variance 1, RE variance 4; | NA                                                                              |
| Covariance structure              | NA                                                                                                                    | Conventional covariance structure, have time as one fixed effect; AR error term |
| Random effects                    | Just intercept; intercept + slope                                                                                     | NA                                                                              |
