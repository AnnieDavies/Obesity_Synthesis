# Obesity_Synthesis
Data and codes for the paper "Interventions to prevent obesity in children aged 5-18 years: a complex synthesis" Davies et al (2024).

The folder Data contains the data to be analysed.

The folder functions contains R functions used in the three R scripts to perform the analysis.

The folder charts contains the age- and sex-specific LMS parameters for different reference charts.

To perform the analysis, R scripts should be run in the following order:
- CreateMappedData.R - maps the BMI and percentile data onto the zBMI scale
- CreateAnalysisData.R - calculates the contrast level observations and formats the data to be read into JAGS
- RunModel.R - runs the analysis using JAGS

The text files contain JAGS models for three versions of the analysis:
- model_primary_RE: random effects (primary) analysis
- model_secondary_FE: fixed effects (secondary) analysis
- model_SSVS: model to perform stochastic search variable selection
