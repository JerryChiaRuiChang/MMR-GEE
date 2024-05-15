# MMR-GEE
Code for the simulation study and the analysis of the Pro-CCM data corresponding to draft manuscript "Estimating Marginal Treatment Effect in Cluster Randomized Trials with Multi-level Missing Outcomes".
The data that support the findings in this paper are openly available at https://doi.org/10.7910/DVN/IIDE2B.

### `data-application` folder
Contains R code for data cleaning and running the analyses for the Pro-CCM data:
- `data-application.R`: main code file that (1) cleans and reorganizes the data from `Proactive community case management decreased malaria prevalence in rural Madagascar: results from a cluster randomized trial` (Ratovoson et. al., 2022) and (2) fits four estimators (CC-GEE, IPW-GEE, MIPW-GEE-EM, and MMR-GEE) to the cleaned data.

