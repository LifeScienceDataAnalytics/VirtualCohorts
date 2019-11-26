# VirtualCohorts
Code for generation of virtual cohorts from Bayesian Network generated on longitudinal data for Azheimer's Disease (AD) and Parkinson's Disease (PD)
Data used for AD is ADNI and PD is PPMI

ADNI
* 1_preProcessing.R: preprocesses ADNI data, deals with missing values, creation of auxilliary variables and imputation
* 2_AutoencodingAndBNCreation.Rmd: autoencoded models generated for group of variables at each time point for e.g. one autoencoded variable generated for 7 brain volumes. Bayesian network is created.
* 3_classification.R: classifier to classify virtual and real patients using coonservative and non-conservative methods
* 4_densityAndBarPlots.R: density and histogram generated for real and virtual patients 
* 5_counterfactualCognition.R: counterfactual scenario generated for virtual patients by changing the baseline score of demented subject to a normal
* 6_hybridBNCreation.R: generation of hybrid bayesian network and comparison with the one generated using discrete data


PPMI 

* utility.R: functions used in the data analysis pipeline
* 1_getData.R: Read, filter and subset longitudinal and non- longitudinal data of PD group from PPMI data
* 2_BN.R: Data pre-processing for Bayesian Model. This includes data discretization, stipulating blacklist and         whitelist, building Bayesian network
* 3_VirtualPatient.R: Simulation of virtual patients. Classifier to classify virtual and real patients using conservative and non-conservative methods
* 4_onePatient.R: Sequential prediction of longitudinal variables in bayesian network
* 5_Counterfactual.R: counterfactual scenario generated for virtual patients by changing the baseline score of PD subject to a normal
* 6_hybridBNCreation.R:generation of hybrid bayesian network and comparison with the one generated using discrete data

