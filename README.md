# Statistical post-processing of visibility ensemble forecasts of the EUPPBench benchmark dataset[^1]

This package contains the implementation of MLP and POLR models for the post-processing of WMO-reported visibility ensemble forecasts of the EUPPBench dataset. 
The default architecture for both of the applied methods are the same that are used by Baran and Lakatos[^2].

## Data description: 

**Observations** are reported according to the WMO suggestions in values: 

Y = {0, 100, 200, . . . , 4900, 5000, 6000, 7000, . . . , 29000, 30000, 35000, 40000, . . . , 65000, 70000}.

**Forecasts**: 
- members: 51-member operational ECMWF ensemble and high-resolution (HRES) 
* the matching of forecasts (given in 1 m steps) and observations is performed by rounding down to the closest
reported value.
+ Time period:  2017 – 2018 initialized at 0000 UTC
- forecast horizon of 120 h and temporal resolution of 6 h
* 42 SYNOP stations in Germany and France 


## Applied models

Proportional odds logistic regression (POLR), Multilayer perceptron (MLP)

Variants: 

- Local estimation (L)
* Regional estimaton (R)
+ Semi-local estimation (C)

|      | Regional training | Local training | Semi-local training |
|:----:|:-----------------:|:--------------:|:-------------------:|
| POLR |       POLR-R: `polr.R`      |     POLR-L: `polrLoc.R ()`   |        POLR-C:    `polrGroup.R`      |
|  MLP |       MLP-R: `mlp.R`        |      MLP-L: `mlpLoc.R`      |        MLP-C:  `mlpGroup.R`       |


## Other scripts

`verif_scores.R`: Verification scores for the evaluation of forecast skill 

`utils.R`: Utility functions

## Additional information

This package is still a work in progress the list of fully implemented scripts are countinously updated: 

- [x] `polr.R`
- [ ] `polrLoc.R`
- [x] `polrGroup.R`
- [x] `mlp.R` 
- [x] `mlpLoc.R`
- [x] `mlpGroup.R`


[^1]: Demaeyer, J., Bhend, J., Lerch,
  S., Primo, C., Van Schaeybroeck, B., Atencia, A., Ben Bouallègue, Z.,
  Chen, J., Dabernig, M., Evans, G., Faganeli Pucer, J., Hooper, B., Horat, N.,
  Jobst, D., Merše, J., Mlakar, P., Möller, A., Mestre, O., Taillardat,
  M. and Vannitsem, S. (2023) The EUPPBench postprocessing benchmark dataset
  v1.0. *Earth Syst. Sci. Data*, doi:10.5194/essd-2022-465.

[^2]: Baran, S. and Lakatos, M. (2023) Statistical post-processing of visibility ensemble forecasts. *arXiv:\/}2305.15325* (submitted)
