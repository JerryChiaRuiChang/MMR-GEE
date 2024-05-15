################################################################################
# data-generation3.R is used to generate the data for the Alt-2 design.
# (1) resamples baseline data from the pro-CCM study
# (1a) resamples 50 cluster with replacement
# (1b) resamples 30 subclusters with replacement
# (2) random assignment of treatment
# (3) outcome generation
# (4) Missing data generation under the alternative design 2
################################################################################

# load the required pacakges
library(geepack)
library(broom.mixed)
library(dplyr)
library(CRTgeeDR)
library(SimCorMultRes)

# main function for generating the data
data_generation <- function(i) {
  load("data-generation/clean_data.Rdata")
  ######################################################################################################
  # Step 1: sample cluster and subcluster id with replacement 
  ######################################################################################################
  # (sampling_id1, sampling_id2): cluster- and subcluster-level bootstrap id
  sampling_id1 = sample(unique(new_data$idfkt), size = 50, replace = TRUE)
  sampling_id2 = c()
  for (i in 1:50) {
    temp <- new_data %>% dplyr::filter(idfkt == sampling_id1[[i]])
    sampling_id2[[i]] = sample(unique(temp[temp$subclusters_size < 4,]$idmenage), size = 30, replace = TRUE)
  }
  
  full_data_subclusters = data.frame()
  full_data = data.frame()
  for (i in 1:50) {
    # Subset data by bootstrap id 
    temp_subclusters = new_data_subclusters %>% filter(idfkt == sampling_id1[i])
    temp_subclusters = temp_subclusters[match(sampling_id2[[i]], temp_subclusters$idmenage),]
    temp_subclusters$idfkt2 = i
    temp_subclusters$idmenage2 = unlist(lapply(1:30, function(x) paste(i, x, sep = "-")))
    
    temp = merge(new_data, temp_subclusters[,c("idfkt", "idmenage", "idfkt2", "idmenage2")], by = c("idfkt", "idmenage"))
    
    full_data_subclusters = rbind(full_data_subclusters, temp_subclusters)
    full_data = rbind(full_data, temp)
  }
  
  ######################################################################################################
  # Step 2: Create random assignment of treatment 
  ######################################################################################################
  trt_assignment <- sample(unique(full_data_subclusters$idfkt2), 11, replace = FALSE)
  full_data_subclusters$A <- 0
  full_data_subclusters[full_data_subclusters$idfkt2 %in% trt_assignment,]$A = 1
  full_data$A <- 0
  full_data[full_data$idfkt2 %in% trt_assignment,]$A = 1
  
  ######################################################################################################
  # Step 3: Generate binary outcome based on logistic model with cluster- and subcluster-level random intercept
  # Beta coefficients are based on GLM using AIC with backward selection
  # Cluster- and subcluster-level random intercepts are based on GLMM with Y ~ A
  # Covariates:  
  # Main effect: A, sexe, age, educ1, educ2, dorm_moust, subclusters_size, sub_educ, IRS, 
  # Interaction: A:educ1, A:educ2, A:dorm_moust, A:sub_educ, A:IRS
  ######################################################################################################
  
  boot_data <- full_data %>% dplyr::select(idfkt2, idmenage2, id_unique, A, sexe, age, niv_scol,
                                           subclusters_size, avg_sub_sex, IRS)
  # Converte age from months to year
  boot_data$sexe <- as.integer(boot_data$sexe) - 1
  boot_data$age <- boot_data$age/12
  boot_data$IRS <- as.integer(boot_data$IRS) - 1
  boot_data$niv_scol <- as.integer(boot_data$niv_scol) - 1
  boot_data$log_subclusters_size = log(boot_data$subclusters_size)
  boot_data <- boot_data %>% dplyr::mutate(Intercept = 1, A_sexe = A * sexe, A_age = A * age,
                                           A_niv_scol = A * niv_scol, A_subclusters_size = A * subclusters_size,
                                           A_avg_sub_sex = A * avg_sub_sex, A_IRS = A * IRS)
  boot_data_subcluster <- full_data_subclusters
  boot_data_subcluster$IRS <- as.integer(boot_data_subcluster$IRS) - 1
  boot_data_subcluster$log_subclusters_size = log(boot_data_subcluster$subclusters_size)
  boot_data_subcluster <- boot_data_subcluster %>% dplyr::mutate(Intercept = 1, A_subclusters_size = A * subclusters_size, 
                                                                 A_avg_sub_sex = A * avg_sub_sex, A_IRS = A * IRS)
  
  
  # Function to create a block exchangeable matrix with varying cluster sizes
  block_exchangeable_matrix <- function(cluster_sizes) {
    n_clusters <- length(cluster_sizes)
    matrix_size <- sum(cluster_sizes)
    
    result_matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
    
    start_row <- 1
    
    for (cluster_size in cluster_sizes) {
      block_matrix <- matrix(rep(1, cluster_size^2), nrow = cluster_size)
      result_matrix[start_row:(start_row + cluster_size - 1),
                    start_row:(start_row + cluster_size - 1)] <- block_matrix
      
      start_row <- start_row + cluster_size
    }
    
    return(result_matrix)
  }
  
  # Function to simulate correlated binary outcomes using the SimCorMultRes package 
  simulate_binary <- function(i, data, subclusters_data, alpha0, alpha1) {
    data <- data[data$idfkt2 == i, ]
    subclusters_data <- subclusters_data[subclusters_data$idfkt2 == i, ]
    # cluster size
    cluster_size <- nrow(data)
    # intercept
    beta_intercepts <- -0.5
    # regression parameter associated with the covariate
    beta_coefficients <- c(0.5, c(0, -0.5, -0.1), c(0, 1, 0.5), c(0, -0.5, -0.1), c(0, 1, 0.5))
    # create nested exchangeable correlation matrix for the NORTA method
    latent_correlation_matrix <- (1 - alpha0) * diag(cluster_size) + 
      (alpha0 - alpha1) * block_exchangeable_matrix(subclusters_data$subclusters_size) + 
      alpha1 * matrix(rep(1, cluster_size^2), nrow = cluster_size)
    # Covariate
    x <- as.matrix(data[,c("A", "sexe", "age", "niv_scol", "subclusters_size", "avg_sub_sex", "IRS",
                           "A_sexe", "A_age", "A_niv_scol", "A_subclusters_size", "A_avg_sub_sex", "A_IRS")])
    # simulation of clustered binary responses
    simulated_binary_dataset <- rbin(clsize = cluster_size, intercepts = beta_intercepts,
                                     betas = beta_coefficients, xformula = ~x, cor.matrix = latent_correlation_matrix,
                                     link = "logit")
    return(simulated_binary_dataset$Ysim)
  }
  
  # Simulate outcomes
  alpha0 = (0.8181^2 + 0.8638^2)/(0.8181^2 + 0.8638^2 + pi^2/3)
  alpha1 = (0.8181^2)/(0.8181^2 + 0.8638^2 + pi^2/3)
  boot_data$Y = unlist(lapply(1:50, function(x) simulate_binary(x, boot_data, boot_data_subcluster, alpha0, alpha1)))
  
  
  ######################################################################################################
  # Step 4: Generate missing data
  ######################################################################################################
  
  # Generate cluster-level missing outcomes
  gamma = c(4, -0.5, c(-1, -0.5, -0.5), c(0, -3.5, -2))
  boot_data_subcluster$lambda_prob = plogis(as.vector(as.matrix(boot_data_subcluster[, c("Intercept", "A", "log_subclusters_size", "avg_sub_sex",
                                                                                         "IRS", "A_subclusters_size", "A_avg_sub_sex",
                                                                                         "A_IRS")]) %*% gamma))
  
  boot_data_subcluster$C = unlist(lapply(boot_data_subcluster$lambda_prob, function(x) rbinom(1, 1, x)))
  boot_data <- merge(boot_data, boot_data_subcluster[, c("idmenage2", "C", "lambda_prob")], by = "idmenage2")
  
  # Generate ind-level missing outcomes
  eta = c(0, -0.5, c(0.2, 0, 1), c(0.2, 0, 0))
  boot_data$phi_prob = plogis(as.vector(as.matrix(boot_data[, c("Intercept", "A", "age", "niv_scol", "log_subclusters_size",
                                                                "A_age", "A_niv_scol", "A_subclusters_size")]) %*% eta))
  boot_data$R = unlist(lapply(boot_data$phi_prob, function(x) rbinom(1, 1, x)))
  boot_data[boot_data$C == 0,]$R = 0
  
  return(list(boot_data, boot_data_subcluster))
}