# This is the code library for the manuscript "Protein profiles predict treatment responses to the PI3K inhibitor umbralisib in patients with chronic lymphocytic leukemia".

Due to the ethical policy, the data is not shared here, but it can be required from the corresponding author Sigrid S Skåland: sigrid.skanland@ous-research.no

Here lists the description of each file in the repository

1_1_pflow_lasso_LOOCV.R: This is for the predictive modelling of treatment response of the drug umbralisib, using Lasso algorithm and (phospho) proteomics data. A nested leave-one-out cross-validation was constructed to evaluate the performance.

1_2_pflow_validation.R: This is for testing the selected proteins in 1_1 on an independent test set. Given that the test set only contains one class, Wilconxon's test were used to examine the predictions.

2_1_clustering.R: This is an unsupervised clustering within the non-responders using the k-means algorithm and DSS data, in order to explore potential treatment options for those who are likely not to respond to umbralisib.

2_2_correlation_plot.R: This is for exploring the correlation between (phospho) proteins and the selected drug combinations in 2_1.

my_function.r: Personal function for build the coding environment, such as loading required packages.

original_yogi_functions.r: This is the original codes from the package "yogiroc": https://github.com/jweile/yogiroc. Given that I need to edit it for better illustrations, the original codes were copied and sourced here.