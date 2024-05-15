######################################################################################################
# run_GEE_EM.R runs GEE with the following methods:
# (1) GEE using the full data (ideal case)
# (2) CC-GEE
# (3) IPW-GEE
# (4) MIPW-GEE with EM
# ==> relies on functions from EM.R
# (5) MMR-GEE with EM 
# ==> relies on functions from EM.R
# ==> relies on the following functions from MR.R
# (5a) compute_ghat_no_EM
# (5b) modNR
######################################################################################################

run_GEE <- function(data, data_subcluster, PS1, PS2_1, PS2_2, cluster_id, subcluster_id,
                        gamma_current = NULL, eta_current = NULL) {
  
  data$id = data[, cluster_id]
  ######################################################################################################
  # Truth
  ######################################################################################################
  Full_GEE <- geeglm(Y ~ A, id = id, data = data, family = binomial(link = "logit"))
  
  ######################################################################################################
  # CC-GEE
  ######################################################################################################
  complete_data <- data %>% filter(R == 1)
  CC_GEE <- geeglm(Y ~ A, id = id, data = complete_data, family = binomial(link = "logit"))
  
  ######################################################################################################
  # IPW-GEE 
  ######################################################################################################
  complete_data <- data
  PS1 <- glm(PS1, data = complete_data, family = binomial(link = "logit"))
  predict_PS1 <- predict(PS1, complete_data, type = "response")
  complete_data$PS<- predict_PS1
  
  # normalize weights
  complete_data <- complete_data %>% filter(R == 1) %>% mutate(weights = 1/PS) %>% mutate(weights = weights/sum(weights))
  IPW_GEE <- geeglm(Y ~ A, id = id, data = complete_data, weights = weights, family = binomial(link = "logit"))
  
  ######################################################################################################
  # MIPW-GEE-no-EM
  ######################################################################################################
  complete_data <- data
  PS2_1 <- glm(PS2_1, data = data_subcluster, family = binomial(link = "logit"))
  predict_PS2_1 <- predict(PS2_1, complete_data, type = "response")
  complete_data$PS2_1 <- predict_PS2_1 
  
  PS2_2 <- glm(PS2_2, data = complete_data[complete_data$C == 1,], family = binomial(link = "logit"))
  predict_PS2_2 <- predict(PS2_2, complete_data, type = "response")
  complete_data$PS2_2<- predict_PS2_2
  
  # normalize weights
  complete_data <- complete_data %>% filter(R == 1) %>% mutate(weights = 1/(PS2_1 * PS2_2)) %>% mutate(weights = weights/sum(weights))
  MIPW_GEE <- geeglm(Y ~ A, id = id, data = complete_data, weights = weights, family = binomial(link = "logit"))
  
  
  ######################################################################################################
  # MIPW-GEE-EM results: 
  ######################################################################################################
  est_nuisance_EM1 <- list()
  est_nuisance_EM2 <- list()
  complete_data <- data
  EM_results <- EM(complete_data, formula(PS2_1), formula(PS2_2), id = subcluster_id)
  est_nuisance_EM1[[1]] <- EM_results[[1]]
  est_nuisance_EM2[[1]] <- EM_results[[2]]
  predict_PS3_1 <- EM_results[[3]]
  predict_PS3_2 <- EM_results[[4]]
  complete_data$PS3_1 <- predict_PS3_1
  complete_data$PS3_2 <- predict_PS3_2
  # normalize weights
  complete_data <- complete_data %>% filter(R == 1 & C == 1) %>% mutate(weights = 1/(PS3_1*PS3_2)) %>% mutate(weights = weights/sum(weights))
  MIPW_GEE_EM <- geeglm(Y ~ A, id = id, data = complete_data, weights = weights, family = binomial("logit"))
  
  ######################################################################################################
  # MMR-GEE
  ######################################################################################################
  # Fit PS model without specifying the treatment-covairates interaction terms 
  complete_data <- data
  P1 = c(formula(PS2_1), c("C ~ A + A:subclusters_size"))
  P2 = c(formula(PS2_2), c("R ~ A + poly(age, 2) + A:age"))
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
        EM_results <- EM(complete_data, P1[[i]], P2[[j]], id = subcluster_id)
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
  complete_data <- complete_data %>% filter(R == 1 & C == 1)
  MMR_GEE <- geeglm(Y ~ A, id = idfkt2, data = complete_data, weights = complete_data$weights, family = binomial("logit"))
  
  TE = c(Full_GEE$coefficients[2], CC_GEE$coefficients[2], IPW_GEE$coefficients[2], 
         MIPW_GEE$coefficients[2], MIPW_GEE_EM$coefficients[2], MMR_GEE$coefficients[2])
  
  results = list(TE, est_nuisance_EM1[[1]], est_nuisance_EM2[[1]], PS2_1$coefficients, PS2_2$coefficients)
  return(results)
}