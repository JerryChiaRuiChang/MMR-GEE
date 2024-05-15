######################################################################################################
# run_GEE_no_EM.R runs GEE with the following methods:
# (1) GEE using the full data (ideal case)
# (2) CC-GEE
# (3) IPW-GEE
# (4) MIPW-GEE without EM
# (5) MMR-GEE without EM 
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
  # MMR-GEE-no-EM
  # PS3_1: C ~ A + A:subclusters_size
  # PS3_2: R ~ A + poly(age, 2) + A:age
  ######################################################################################################
  complete_data <- data
  PS3_1 <- glm(C ~ A + A:subclusters_size, data = data_subcluster, family = binomial(link = "logit"))
  predict_PS3_1 <- predict(PS3_1, complete_data, type = "response")
  
  PS3_2 <- glm(R ~ A + poly(age, 2) + A:age, data = complete_data[complete_data$C == 1,], family = binomial(link = "logit"))
  predict_PS3_2 <- predict(PS3_2, complete_data, type = "response")
  pred_prob_cluster <- cbind(predict_PS2_1, predict_PS3_1)
  pred_prob_subject <- cbind(predict_PS2_2, predict_PS3_2)
  ghat <- compute_ghat_no_EM(pred_prob_cluster, pred_prob_subject)
  
  r <- complete_data$R
  out <- try(modNR(ghat, r))
  weights <- as.vector(r * (1 + t(out$rhohat) %*% ghat)^(-1) * (sum(r))^(-1))
  complete_data$weights <- weights
  complete_data <- complete_data %>% filter(R == 1 & C == 1)
  MMR_GEE <- geeglm(Y ~ A, id = id, data = complete_data, weights = weights, family = binomial(link = "logit"))
  
  
  TE = c(Full_GEE$coefficients[2], CC_GEE$coefficients[2], IPW_GEE$coefficients[2], 
         MIPW_GEE$coefficients[2], MMR_GEE$coefficients[2])
  
  return(TE)
}