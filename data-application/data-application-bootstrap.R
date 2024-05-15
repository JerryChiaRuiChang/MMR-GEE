######################################################################################################
# data-application-bootstrap.R is used for obtaining standard error estimates for the analysis of the
# Pro-CCM study through the cluster bootstrap approach. 
# ----------------------------------------------------------------------------------------------------
# Required files: 
# (1) clean_data.Rdata: cleaned version of the data from Ratovson et. al. (2022), which was obtained
# from running data-application.R 
# (2) EM.R for running the EM algorithm
# (3) MR.R for estimating the multiply robust weights
# ----------------------------------------------------------------------------------------------------
# Output: 
# Standard error estimates and 95% CI
######################################################################################################

# Load required packages
library(geepack)
library(broom.mixed)
library(dplyr)
library(CRTgeeDR)

# load the cleaned version of the dataset
load("clean_data.Rdata")

# Source the helper functions EM.R and MR.R for running the EM algorithm and the proposed 
# MMR-GEE estimator
source("data-application/EM.R")
source("data-application/MR.R")

######################################################################################################
# run_bootstrap: main function to run the cluster bootstrap approach
######################################################################################################
run_bootstrap <- function(index) {
  set.seed(index)
  # Sample clusters with replacement 
  sampling_id1 = sample(unique(new_data$idfkt), size = length(unique(new_data$idfkt)), replace = TRUE)
  full_data_subclusters = data.frame()
  full_data = data.frame()
  
  # Rename the id of each cluster
  for (i in 1:22) {
    temp_subclusters = new_data_subclusters %>% filter(idfkt == sampling_id1[i])
    temp_subclusters$idfkt2 = i
    temp_subclusters$idmenage2 = unlist(lapply(temp_subclusters$idmenage, function(x) paste(i, x, sep = "-")))
    temp = merge(new_data, temp_subclusters[,c("idfkt", "idmenage", "idfkt2", "idmenage2")], by = c("idfkt", "idmenage"))
    full_data_subclusters = rbind(full_data_subclusters, temp_subclusters)
    full_data = rbind(full_data, temp)
  }
  
  # CC-GEE  
  complete_data <- full_data %>% filter(R == 1 & C == 1)
  CC_GEE <- geeglm(tdr1 ~ Grp, id = idfkt2, data = complete_data, family = binomial(link = "logit"), 
                     corstr="exchangeable")
  CC_GEE_results <- CC_GEE$coefficients[2]
  
  # IPW-GEE  
  complete_data <- full_data 
  PS1_form = "R ~ Grp + sexe + age + niv_scol + dorm_moust + subclusters_size + 
  avg_sub_sex + sub_educ + IRS + Grp:sexe + Grp:age + Grp:niv_scol + 
  Grp:dorm_moust + Grp:subclusters_size + Grp:avg_sub_sex + 
  Grp:IRS"
  PS1 <- glm(PS1_form, data = complete_data, family = binomial(link = "logit"))
  complete_data$PS <- fitted(PS1)
    
  # normalize weights
  complete_data <- complete_data %>% filter(R == 1 & C == 1) %>% mutate(w = 1/PS) %>% mutate(w = w/sum(w))
  IPW_GEE <- geeglm(tdr1 ~ Grp, id = idfkt2, data = complete_data, weights = w, 
                    family = binomial("logit"), corstr="exchangeable")
  IPW_GEE_results <- IPW_GEE$coefficients[2]
    
  # MIPW-GEE-EM: 
  complete_data <- full_data
  complete_data <- full_data 
  PS2_1_form = "C ~ Grp + subclusters_size + sub_educ + IRS + Grp:IRS"
  PS2_2_form = "R ~ Grp + sexe + age + niv_scol + dorm_moust + subclusters_size + 
  avg_sub_sex + sub_educ + IRS + Grp:subclusters_size + Grp:avg_sub_sex + Grp:sub_educ"
  EM_results <- EM(complete_data, Z = PS2_1_form, X = PS2_2_form, id = "idmenage2")
  predict_PS3_1 <- EM_results[[3]]
  predict_PS3_2 <- EM_results[[4]]
  complete_data$PS3_1 <- predict_PS3_1
  complete_data$PS3_2 <- predict_PS3_2
  
  # normalize weights
  complete_data <- complete_data %>% filter(R == 1 & C == 1) %>% mutate(weights = 1/(PS3_1*PS3_2)) %>% mutate(weights = weights/sum(weights))
  MIPW_GEE_EM <- geeglm(tdr1 ~ Grp, id = idfkt2, data = complete_data, weights = weights, 
                        family = binomial("logit"), corstr="exchangeable")
  MIPW_GEE_EM_results <- MIPW_GEE_EM$coefficients[2]
    
  # MMR-GEE results:
  complete_data <- full_data
  P1 = c(PS2_1_form, c("C ~ Grp + Grp:subclusters_size + Grp:sub_educ"))
  P2 = c(PS2_2_form, c("R ~ Grp + poly(age, 2) + niv_scol + dorm_nat + IRS + Grp:IRS"))
  pred_prob_cluster = matrix(, nrow = nrow(complete_data), ncol = 4)
  pred_prob_subject = matrix(, nrow = nrow(complete_data), ncol = 4)
  count = 0
  for (i in 1:2) {
    for (j in 1:2) {
      count = count + 1
      if (i == 1 & j == 1) {
        pred_prob_cluster[,count] <- predict_PS3_1
        pred_prob_subject[,count] <- predict_PS3_2
      } else {
        EM_results <- EM(complete_data, P1[[i]], P2[[j]], id = "idmenage2")
        pred_prob_cluster[,count] <- EM_results[[3]]
        pred_prob_subject[,count] <- EM_results[[4]] 
      }
    }
  }
  
  ghat <- compute_ghat(pred_prob_cluster, pred_prob_subject)
  r <- complete_data$R
  out <- try(modNR(ghat, r))
  weights <- as.vector(r * (1 + t(out$rhohat) %*% ghat)^(-1) * (sum(r))^(-1))
  complete_data$weights <- weights
  complete_data <- complete_data %>% filter(R == 1 & C == 1) 
  
  MMR_GEE <- geeglm(tdr1 ~ Grp, id = idfkt2, data = complete_data, weights = weights, 
                      family = binomial("logit"), corstr="exchangeable")
  MMR_GEE_results <- MMR_GEE$coefficients[2]
    
  results <- c(CC_GEE_results, IPW_GEE_results, MIPW_GEE_EM_results, MMR_GEE_results)
  return(results)
}

# Run 1,000 cluster bootstrap
boot_results = matrix(NA, nrow = 1000, ncol = 4)
for (index in 1:1000) {
  tryCatch({
    boot_results[index,] = run_bootstrap(index)
    print(paste0("Bootstrap ", index, " completed"))
  }, error = function(e) {
    # Code to handle the error
    print(paste("Error occurred at current iteration: ", index))
    print(conditionMessage(e))
  })
}

# Obtain standard error estimates and 95% CI
beta_est = c(-0.2856294, -0.2050506, -0.1660465, -0.1959094)
# some of the bootstrap estimates did not converge so we excluded them
SD = apply(boot_results[which(boot_results[,3] < Inf),], 2, sd, na.rm = TRUE)
lower_CI = exp(beta_est - 1.96 * SD)
upper_CI = exp(beta_est + 1.96 * SD)