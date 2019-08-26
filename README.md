This repository includes data and analysis scripts to accompany:

# Nitrogen pollution interacts with thermal stress to increase coral bleaching across the seascape

### Author of analysis and code: Mary K Donovan
### Currently submitted for peer review

-----

### Description:
This work analyzes patterns of coral bleaching around the island of Moorea, French Polynesia in relationship to patterns of heat stress and nitrogen pollution. Analyses include calculating spatial patterns of heat stress, mapping patterns of bleaching, and modeling the relationship between bleaching prevalence and severity as a function of heat stress, nitrogen, and their interaction with a Bayesian hierarchical model.

### Contents:
#### Scripts:
* **1_temperature_patterns.R:** R script that combines contemporaneous temperature patterns from in situ loggers at 6 Long Term Ecological Research sites from the Moorea Coral Reef LTER, and calculates metrics of heat stress during a mild bleaching event in 2016.
* **2_site_map_nutirent_patterns.R:** R script that creates a map of nutrient patterns around Moorea
* **3_bleaching_patterns.R:** R script that creates data inputs for maps of bleaching patterns around Moorea
* **4a_bayes_model_binom.R:** R script that creates a JAGS model object and runs a Bayesian Heirarchical model with JAGS in R for bleaching prevalence assuming a binomial response
* **4b_bayes_model_beta.R:** R script that creates a JAGS model object and runs a Bayesian Heirarchical model with JAGS in R for bleaching severity assuming a beta response
* **5_model_outputs.R:** R script that summarizes JAGS model outputs and creates figures
