# MMR-GEE
Code for the simulation study and the analysis of the Pro-CCM data corresponding to draft manuscript "Estimating Marginal Treatment Effect in Cluster Randomized Trials with Multi-level Missing Outcomes".
The data that support the findings in this paper are openly available at https://doi.org/10.7910/DVN/IIDE2B.

### `data-application` folder
Contains R code for data cleaning and running the analyses for the Pro-CCM data:
- `data-application.R`: main code file that (1) cleans and reorganizes the data from Proactive community case management decreased malaria prevalence in rural Madagascar: results from a cluster randomized trial (Ratovoson et. al., 2022) and (2) fits four estimators (CC-GEE, IPW-GEE, MIPW-GEE-EM, and MMR-GEE) to the cleaned data.
- `data-application-bootstrap.R`: code for obtaining standard error estimates through the cluster bootstrap approach.
- `EM.R`: is the helper function for implementing the EM algorithm. 
- `MR.R`: is the helper function for estimating the multiply robust weights.
  
### `simulation` folder
Contains R code for running the simulation based on the Pro-CCM study:
- `run-methods.R`: main code file for (1) estimating beta parameters and (2) obtaining bootstrap standard error using the cluster bootstrap approach for the simulation.
- `data-generation`: subfolder includes `data-generation1.R`, `data-generation2.R`, `data-generation3.R` for simulating the data under the Org-Pro-CCM, Alt-1, and Alt-2 design, respectively.
- `methods`: subfolder includes `EM.R` for implementing the EM algorithm, `MR.R` for estimating the multiply robust weights, `run-GEE-no-EM.R` for implementing CC-GEE, IPW-GEE, MIPW-GEE, and MMR-GEE without the EM algorithm, `run-GEE-EM.R` for implementing CC-GEE, IPW-GEE, MIPW-GEE, and MMR-GEE with the EM algorithm.

