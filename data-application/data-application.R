######################################################################################################
# data-application.R is used for the analysis of the Pro-CCM study. The file does the following jobs:
# (1) cleans and reorganizes the data from Ratovson et. al. (2022)
# (2) fits four estimators (CC-GEE, IPW-GEE, MIPW-GEE-EM, and MMR-GEE) to the cleaned data 
# ----------------------------------------------------------------------------------------------------
# Description of datasets:
# proccm.intervention.data: longitudinal data about different outcomes
# survey.data: baseline survey and trial data (the only one used for analysis)
# health.center.data: longitudinal data on the population level 
# ----------------------------------------------------------------------------------------------------
# Required files: 
# (1) EM.R for running the EM algorithm
# (2) MR.R for estimating the multiply robust weights
# ----------------------------------------------------------------------------------------------------
# Output: 
# clean_data.Rdata: cleaned version of the dataset 
# CC_GEE, IPW_GEE, MIPW_GEE_EM, MMR_GEE: results from fitting the four estimators
######################################################################################################
# load the data
#load("Ratovson et al 2021_Databases.RData")
baseline_data <- survey.data[survey.data$Period == 0,]
endline_data <- survey.data[survey.data$Period == 1,]

# Load required packages
library(lme4)
library(geepack)
library(broom.mixed)
library(dplyr)
library(CRTgeeDR)

# Source the helper functions EM.R and MR.R for running the EM algorithm and the proposed MMR-GEE estimator
source("data-application/EM.R")
source("data-application/MR.R")

######################################################################################################
# Reorganize and clean data, yielding three datasets
# new_data: individaul-level data
# new_data_subclusters: subcluster-level data
# new_data_clusters: cluster-level data
# All three datasets were saved as clean_data.R 
######################################################################################################
baseline_id <- unique(baseline_data$id_unique)
endline_id <- unique(endline_data$id_unique)
new_data <- baseline_data

# create missingness indicators R_{ijk} and add newcomers indicators
new_data$tdr1 <- NA
for (id in endline_id) {
  if (id %in% baseline_id) {
    new_data[new_data$id_unique == id,]$tdr1 = endline_data[endline_data$id_unique == id,]$tdr 
  }
}
new_data <- new_data %>% mutate(R = ifelse(is.na(tdr1), 0, 1))

# Add newcomers indicators
newcomers_id <- endline_id[!(endline_id %in% baseline_id)]
newcomers_data <- endline_data %>% filter(id_unique %in% newcomers_id)
newcomers_data <- rename(newcomers_data, tdr1 = tdr)
newcomers_data$tdr <- NA
newcomers_data$R <- 1
newcomers_data$new <- 1
new_data$new <- 0
new_data <- rbind(new_data, newcomers_data)

# Convert age into four groups - 0-4, 5-14, 15+
new_data <- fastDummies::dummy_cols(new_data, select_columns = "Classage")
new_data <- new_data %>% rename("age_0_5" = "Classage_[0,5)", "age_5_15" = "Classage_[5,15)", "age_15_plus" = "Classage_15+")
new_data <- new_data %>% group_by(idmenage) %>% mutate(Class_age_0_5 = mean(age_0_5), Class_age_5_15 = mean(age_5_15), 
                                                       Class_age_15_plus = mean(age_15_plus)) %>% ungroup()

new_data <- fastDummies::dummy_cols(new_data, select_columns = "niv_scol")
new_data <- new_data %>% rename("educ0" = "niv_scol_jamais scolarise", "educ1" = "niv_scol_Ecole primaire", 
                                "educ2" = "niv_scol_Ecole secondaire", "educ3" = "niv_scol_Niveau superieur")

# Create subcluster-level covariates and missingness indicators C_ij
new_data_subclusters <- 
  new_data %>% 
  group_by(idmenage) %>% 
  summarise(idfkt = first(idfkt), subclusters_size = n(), avg_sub_sex = mean(as.integer(sexe)-1), 
            avg_age_0_5 = mean(age_0_5), avg_age_5_15 = mean(age_5_15), avg_age_15_plus = mean(age_15_plus), 
            sub_educ = max(as.integer(niv_scol)-1), avg_sub_dorm_moust = mean(dorm_moust), 
            avg_sub_dorm_nat = mean(dorm_nat), Grp = first(Grp), IRS = first(IRS), C = ifelse(sum(R) == 0, 0, 1))

# Create cluster-level covariates
new_data_clusters <- 
  new_data %>% 
  group_by(idfkt) %>% 
  summarise(clusters_size = n(), avg_sex = mean(as.integer(sexe)-1), avg_age = mean(age), avg_educ = mean(as.integer(niv_scol)-1),
            avg_dorm_moust = mean(dorm_moust), avg_dorm_nat = mean(dorm_nat))

new_data_subclusters <- inner_join(new_data_subclusters, new_data_clusters, by = "idfkt")
new_data <- inner_join(new_data, new_data_subclusters, by = "idmenage") %>%
  dplyr::select(-ends_with(".x"))
new_data <- new_data %>% rename(Grp = Grp.y, IRS = IRS.y, idfkt = idfkt.y)
save(new_data, new_data_subclusters, new_data_clusters, file = "clean_data.Rdata")

# Get some summary statistics of the clean data
# 31% of overall individuals missing, 22.3% of missing subclusters 
mis_pct <- c(1 - mean(new_data$R), 1 - mean(new_data_subclusters$C))


######################################################################################################
# Fit four estimators (CC-GEE, IPW-GEE, MIPW-GEE-EM, and MMR-GEE) to the cleaned data
######################################################################################################

# CC-GEE:
complete_data <- new_data %>% filter(R == 1 & C == 1)
CC_GEE <- geeglm(tdr1 ~ Grp, id = idfkt, data = complete_data, 
                 family = binomial(link = "logit"), corstr="exchangeable")


# IPW-GEE: the covariates for the PS model were selected based on backward step-wise procedure using the AIC
complete_data <- new_data 
PS1 <- glm(R ~ Grp*(sexe + age + niv_scol + dorm_moust + dorm_nat + subclusters_size + 
                      avg_sub_sex + sub_educ + IRS), 
           data = complete_data, family = binomial(link = "logit"))
PS1 <- stepAIC(PS1, direction = 'backward')
complete_data$PS <- fitted(PS1)

# normalize weights
complete_data <- complete_data %>% filter(R == 1 & C == 1) %>% mutate(w = 1/PS) %>% mutate(w = w/sum(w))
IPW_GEE <- geeglm(tdr1 ~ Grp, id = nomfkt, data = complete_data, weights = w, 
                  family = binomial("logit"), corstr="exchangeable")


# MIPW-GEE-EM: The covariates for the PS model were selected based on backward step-wise procedure using the AIC
# and the EM algorithm is implemented by the EM.R helper functions
est_nuisance_EM1 <- list()
est_nuisance_EM2 <- list()
complete_data <- new_data
PS2_1 <- glm(C ~ Grp*(subclusters_size + avg_sub_sex + sub_educ + IRS),
             data = new_data_subclusters, family = binomial(link = "logit"))
PS2_1 <- stepAIC(PS2_1, direction = 'backward')
PS2_2 <- glm(R ~ Grp*(sexe + age + niv_scol + dorm_moust + dorm_nat + 
                        subclusters_size + avg_sub_sex + sub_educ +  IRS),
             data = complete_data[complete_data$C == 1,], family = binomial(link = "logit"))
PS2_2 <- stepAIC(PS2_2, direction = 'backward')
Z = formula(PS2_1)
X = formula(PS2_2)
EM_results <- EM(complete_data, Z, X, id = "idmenage")
est_nuisance_EM1[[1]] <- EM_results[[1]]
est_nuisance_EM2[[1]] <- EM_results[[2]]
predict_PS3_1 <- EM_results[[3]]
predict_PS3_2 <- EM_results[[4]]
complete_data$PS3_1 <- predict_PS3_1
complete_data$PS3_2 <- predict_PS3_2

# normalize weights
complete_data <- complete_data %>% filter(R == 1 & C == 1) %>% mutate(weights = 1/(PS3_1*PS3_2)) %>% mutate(weights = weights/sum(weights))
MIPW_GEE_EM <- geeglm(tdr1 ~ Grp, id = nomfkt, data = complete_data, weights = weights, 
                      family = binomial("logit"),  corstr="exchangeable")

# MMR-GEE: the multiply robust is estimated by the MR.R helper functions
complete_data <- new_data
P1 = c(formula(PS2_1), c("C ~ Grp + Grp:subclusters_size + Grp:sub_educ"))
P2 = c(formula(PS2_2), c("R ~ Grp + poly(age, 2) + niv_scol + dorm_nat + IRS + Grp:IRS"))
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
      EM_results <- EM(complete_data, P1[[i]], P2[[j]], id = "idmenage")
      est_nuisance_EM1[[count]] <- EM_results[[1]]
      est_nuisance_EM2[[count]] <- EM_results[[2]]
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
MMR_GEE <- geeglm(tdr1 ~ Grp, id = nomfkt, data = complete_data, weights = weights, 
                  family = binomial("logit"),  corstr="exchangeable")

# Summary of the four estimates
beta = c(CC_GEE$coefficients[2], IPW_GEE$coefficients[2], MIPW_GEE_EM$coefficients[2], MMR_GEE$coefficients[2])
odds_ratio = exp(beta)

# beta = c(-0.2856294, -0.2050506, -0.1660465, -0.1959094)
# OR = c(0.7515411, 0.8146061, 0.8470069, 0.8220867)