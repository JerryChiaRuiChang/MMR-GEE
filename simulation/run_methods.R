######################################################################################################
# run_methods.R is used for (1) estimating beta parameters and (2) obtaining bootstrap standard 
# error using the cluster bootstrap approach for the simulation. The codes
# provided here runs the simulation under the Alt-2 design, which can be 
# easily modified for both the Org-Pro-CCM and Alt-1 design.
# ----------------------------------------------------------------------------------------------------
# The code requires the following files:
# data-generation folder: data_generation1.R/data_generation2.R/data_generation3.R  
# methods folder: run_GEE_no_EM.R/run_GEE_EM.R/MR.R/EM.R
######################################################################################################

# load required packages/functions
library(geepack)
library(broom.mixed)
library(dplyr)
library(CRTgeeDR)
library(SimCorMultRes)

source("data-generation/data-generation2.R")
source("methods/run-GEE-EM.R")
source("methods/EM.R")
source("methods/MR.R")
# use run-GEE-no-EM.R for the Org-Pro-CCM and Alt-1 designs
#source("methods/run-GEE-no-EM.R") 

# function to generate bootstrapped data
bootstrap <- function(data, data_subcluster) {
  sampling_id = sample(unique(data$idfkt2), size = length(unique(data$idfkt2)), replace = TRUE)
  boot_data_subcluster = data.frame()
  boot_data = data.frame()
  for (i in 1:length(unique(data$idfkt2))) {
    temp <- data %>% dplyr::filter(idfkt2 == sampling_id[[i]])
    temp_subcluster = data_subcluster %>% filter(idfkt2 == sampling_id[i])
    temp_subcluster$idfkt3 = i
    temp_subcluster$idmenage3 = unlist(lapply(1:30, function(x) paste(i, x, sep = "-")))
    
    temp = merge(data, temp_subcluster[,c("idfkt2", "idmenage2", "idfkt3", "idmenage3")], by = c("idfkt2", "idmenage2"))
    
    boot_data_subcluster = rbind(boot_data_subcluster, temp_subcluster)
    boot_data = rbind(boot_data, temp)
  }
  return(list(boot_data, boot_data_subcluster))
}

# main function to run the simulation
run_methods <- function(index) {
  # Generate Data 
  data <- data_generation(index)
  data_subcluster <- data[[2]]
  data <- data[[1]]
  
  # Rename subcluster-level missingness indicators
  # C_truth is the true subcluster-level missingness indicator
  # C is the observed cluster-level missingness indicator
  data <- data %>% rename(C_truth = C)
  data_subcluster <- data_subcluster %>% rename(C_truth = C)
  
  C = data %>% group_by(idmenage2) %>% dplyr::summarize(C = as.integer(sum(R) != 0))
  data = merge(data, C, by = "idmenage2")
  data_subcluster = merge(data_subcluster, C, by = "idmenage2")
  
  # Original design
  #PS1 = PS2_2 =  "R ~ A + sexe + age + dorm_moust + subclusters_size + sub_educ + IRS + A_subclusters_size + A_sub_educ"
  #PS2_1 = "C ~ subclusters_size + sub_educ"
  
  # Alternative design 1
  PS1 = PS2_2 =  "R ~ A * age + log_subclusters_size"
  PS2_1 = "C ~ A*(avg_sub_sex + IRS)"
  
  # Alternative design 2
  #PS1 = PS2_2 =  "R ~ A * age + log_subclusters_size"
  #PS2_1 = "C ~ A*(avg_sub_sex + IRS) + log_subclusters_size"
  
  results = run_GEE(data, data_subcluster, PS1, PS2_1, PS2_2, 
                    cluster_id = "idfkt2", subcluster_id = "idmenage2")
  

  boot_results = matrix(, nrow = 0, ncol = 6)
  current = 1
  while (nrow(boot_results) < 200) {
     set.seed(current)
     bootstrap_data <- bootstrap(data, data_subcluster)
     boot_data <- bootstrap_data[[1]]
     boot_data_subcluster <- bootstrap_data[[2]]
     
    tryCatch({
      temp <- run_GEE(boot_data, boot_data_subcluster, PS1, PS2_1, PS2_2, 
                      cluster_id = "idfkt3", subcluster_id = "idmenage3")
      boot_results <- rbind(boot_results, temp[[1]])
      boot_results = unique(boot_results)
    }, error = function(e) {
       # Code to handle the error
      print(paste("Error occurred at current iteration:", nrow(boot_results)))
      print(conditionMessage(e))
    })
    current = current + 1
  }
  
  return(list(results, boot_results))
}

# Example of running one iteration
index = 1
set.seed(index)
print(paste0("running iteration ", index))
results = run_methods(index)
# beta estimates
beta = results[[1]][[1]]
# bootstrap standard error
boot_results = results[[2]]
SD = apply(boot_results, 2, sd, na.rm = TRUE)
# constructing 95% CI
lower = exp(beta + qnorm(0.025) * SD)
upper = exp(beta + qnorm(0.975) * SD)

  