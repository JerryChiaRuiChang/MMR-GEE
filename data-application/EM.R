######################################################################################################
# EM.R is the helper function for implementing the EM algorithm. The file contains three functions:
# (1) E_step: function to compute expectation 
# (2) Q: function for maximization 
# (3) EM: main function for implementing the EM algorithm
# ----------------------------------------------------------------------------------------------------
# Input for EM(data, Z, X, id): 
# (1) data
# (2) Z: the formula for (sub)cluster-level PS models
# (3) X: the formula for individual-level PS models
# (4) id: (sub)cluster-level ID
# ----------------------------------------------------------------------------------------------------
# Output: 
# (1) (gamma_current, eta_current): estimated gamma and eta
# (2) (PS1, PS2): the fitted probability based on gamma_current and eta_current
# (3) converged: indicator of whether the EM converges or not (1 = converged, 0 = not converged)
######################################################################################################
# load required package
library(sjmisc)
library(dplyr)

# E_step: function for the the E-step
E_step <- function(subject_data, subject_covariates, cluster_data, cluster_covariates, gamma_iter, eta_iter, rounding = FALSE) {
  # compute expit(eta_iter Z*) and expit(gamma _iter X*) 
  cluster_data$expit_cluster <- subject_data$expit_cluster <- as.vector(plogis(gamma_iter %*%  t(cluster_covariates)))
  subject_data$expit_subject <- as.vector(plogis(eta_iter %*%  t(subject_covariates)))
  
  cluster_data <-  cluster_data %>% distinct(cluster, .keep_all = TRUE)
  
  # Compute C_iter
  cluster_data$product_one_minus_expit_subject <- aggregate(expit_subject ~ cluster, data = subject_data, FUN = function(x) prod(1 - x))[,2]
  cluster_data = dplyr::mutate(cluster_data, C_iter = C + 
                                 ((1 - C) * (product_one_minus_expit_subject * expit_cluster)/(1 - (expit_cluster * (1 - product_one_minus_expit_subject)))))
  
  return(cluster_data$C_iter)
}

# Q: function for the M-step
Q <- function(subject_data, subject_covariates, cluster_data, cluster_covariates, gamma_iter, eta_iter, C_iter) {
  C_iter <- rep(C_iter, dplyr::distinct(cluster_data)$SIZE)
  cluster_data$C_iter <- subject_data$C_iter <- C_iter
  
  # compute expit(etaX*) and expit(gammaZ*) 
  cluster_data$expit_cluster <- subject_data$expit_cluster <- as.vector(plogis(gamma_iter %*%  t(cluster_covariates)))
  subject_data$expit_subject <- as.vector(plogis(eta_iter %*%  t(subject_covariates)))
  
  # compute sum log(1 - expit(etaX*)) for each cluster
  cluster_data <-  cluster_data %>% distinct(cluster, .keep_all = TRUE)
  
  cluster_data$sum_log_one_minus_expit_subject <- aggregate(expit_subject ~ cluster, data = subject_data, FUN = function(x) sum(log(1 - x)))[,2]
  subject_data$sum_log_one_minus_expit_subject <- rep(cluster_data$sum_log_one_minus_expit_subject, dplyr::distinct(cluster_data)$SIZE)
  
  # Compute complete data log likelihood
  subject_data <- dplyr::mutate(subject_data, first_term = R * (log(expit_subject)) + (1-R) * (log(1 - expit_subject)))
  
  ## R = 1 and expit_subject = 1 produces NaN
  if (sum(subject_data$R == 1 & subject_data$expit_subject == 1) > 0) {
    subject_data[subject_data$R == 1 & subject_data$expit_subject == 1,]$first_term <- 0 
    print("first scenario")
  }
  
  ## R = 0 and expit_subject = 0 produces NaN
  if (sum(subject_data$R == 0 & subject_data$expit_subject == 0) > 0) {
    subject_data[subject_data$R == 0 & subject_data$expit_subject == 0,]$first_term <- 0 
    print("second scenario")
  }
  
  ## R = 0 and expit_subject = 1 produces -Inf
  if (sum(subject_data$R == 0 & subject_data$expit_subject == 1) > 0) {
    subject_data[subject_data$R == 0 & subject_data$expit_subject == 1,]$first_term <- -1000
    print("third scenario")
  }
  
  ## R = 1 and expit_subject = 0 produces -Inf
  if (sum(subject_data$R == 1 & subject_data$expit_subject == 0) > 0) {
    subject_data[subject_data$R == 1 & subject_data$expit_subject == 0,]$first_term <- -1000
    print("fourth scenario")
  }
  
  first_term <- aggregate(first_term ~ cluster, data = subject_data, FUN = sum)[,2]
  cluster_data$first_term <- first_term
  cluster_data <- dplyr::mutate(cluster_data,
                                first_term = C * (log(expit_cluster) + first_term),
                                second_term = (1 - C_iter) * (1 - C)  * log(1 - expit_cluster),
                                third_term = C_iter * (1 - C)  * (log(expit_cluster) + sum_log_one_minus_expit_subject))
  
  return(-sum(cluster_data$first_term + cluster_data$second_term + cluster_data$third_term ))
}


# EM: main function running the EM algorithm
EM <- function(data, Z, X, id, gamma_current = NULL, eta_current = NULL) {
  full_data <- data
  full_data$intercept <- 1
  full_data <- rename(full_data, cluster = id)
  subject_data <- full_data
  cluster_data <- full_data
  
  gamma_wrong <- glm(formula = Z, data = cluster_data, family = binomial)
  eta_wrong <- glm(formula = X, data = subject_data, family = binomial)
  
  # Get cluster-level data
  cluster_covariates <-  as.data.frame(model.matrix(gamma_wrong))
  cluster_data <- as.data.frame(cbind(cluster = cluster_data$cluster, SIZE = cluster_data$subclusters_size, 
                                      cluster_covariates, C = cluster_data$C))
  Z = colnames(cluster_covariates)
  
  # Get subject-level data
  subject_covariates <- as.data.frame(model.matrix(eta_wrong))
  subject_data <- as.data.frame(cbind(cluster = subject_data$cluster, SIZE = subject_data$subclusters_size, 
                                      subject_covariates, R = subject_data$R, C = subject_data$C))
  X = colnames(subject_covariates)
  
  iter = 0
  converged = FALSE
  # initalize eta and gamma
  if (is.null(gamma_current) & is.null(eta_current)) {
    gamma_current <- gamma_wrong$coefficients
    eta_current <-  eta_wrong$coefficients 
  }
  param_EM = matrix(c(gamma_current, eta_current), 
                    ncol =  ncol(cluster_covariates) + ncol(subject_covariates))
  
  while(iter < 1501) {
    # Compute C_iter: conditional expectation of C_truth at the current iteration
    C_iter = E_step(subject_data, subject_covariates, cluster_data, cluster_covariates, gamma_current, eta_current)
    gamma_current = optim(gamma_current, Q, subject_data = subject_data, subject_covariates = subject_covariates,
                          cluster_data = cluster_data, cluster_covariates = cluster_covariates, 
                          C_iter = C_iter, eta_iter = eta_current, method =  "Nelder-Mead")
    gamma_current = gamma_current$par
    eta_current = optim(eta_current, Q, subject_data = subject_data, subject_covariates = subject_covariates,
                        cluster_data = cluster_data, cluster_covariates = cluster_covariates, 
                        C_iter = C_iter, gamma_iter = gamma_current, method =  "Nelder-Mead")
    eta_current = eta_current$par
    param_EM = rbind(param_EM, c(gamma_current, eta_current))
    
    if (dist(param_EM[c(nrow(param_EM)-1,nrow(param_EM)),]) < 0.000001) {
      converged = TRUE
      break
    }
    iter = iter + 1
    print(iter)
  }
  
  # Get estimated gamma and eta
  PS1 <- as.vector(plogis(gamma_current %*%  t(cluster_covariates)))
  PS2 <- as.vector(plogis(eta_current %*%  t(subject_covariates)))
  return(list(gamma_current, eta_current, PS1, PS2, converged))
}
