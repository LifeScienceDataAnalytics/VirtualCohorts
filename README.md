# VirtualCohorts
Code for generation of virtual cohorts from Bayesian Network generated on longitudinal data for Azheimer's Disease (AD) and Parkinson's Disease (PD)
Data used for AD is ADNI and PD is PPMI
For ADNI data
1_preProcessing.R: preprocesses ADNI data, deals with missing values, creation of auxilliary variables and imputation
2_AutoencodingAndBNCreation.Rmd: autoencoded models generated for group of variables at each time point for e.g. one autoencoded variable generated for 7 brain volumes. Bayesian network is created.
3_classification.R: classifier to classify virtual and real patients using coonservative and non-conservative methods
4_densityAndBarPlots.R: density and histogram generated for real and virtual patients 
5_counterfactualCognition.R: counterfactual scenario generated for virtual patients by changing the baseline score of demented subject to a normal
6_hybridBNCreation.R: generation of hybrid bayesian network and comparison with the one generated using dicrete data
