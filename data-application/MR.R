######################################################################################################
# MR.R is the helper function for estimating the multiply robust weights. The file contains two functions:
# (1) compute_ghat: function to compute ghat matrix 
# (2) modNR: main function for running the modified newton raphson algorithm to obtain rho hat
# ----------------------------------------------------------------------------------------------------
# Input for modNR <- function(ghat, r, maxit = 10000, eps = 10e-6): 
# (1) ghat
# (2) r: the missingness indicator
# (3) maxit: maximum number of iterations
# (4) eps: the tolerance level for convergence
# ----------------------------------------------------------------------------------------------------
# Output: 
# rhohat: estimated rho
######################################################################################################

# compute_ghat: function to compute ghat matrix
compute_ghat <- function(pred_prob_cluster, pred_prob_subject) {
  N <- nrow(pred_prob_cluster)
  K <- ncol(pred_prob_cluster)
  L <- ncol(pred_prob_subject)
  
  ghat <- matrix(0, nrow = K, ncol = N)
  index = 0
  for (i in 1:K) {
    index = index + 1
    pihats <- pred_prob_cluster[,i]
    lhats <- pred_prob_subject[,i]
    ghat[index, ] <- pihats* lhats - mean(pihats*lhats)
  }
  return(ghat)
}

# modNR: main function for running the modified newton raphson algorithm to obtain rho hat
modNR <- function(ghat, r, maxit = 10000, eps = 10e-6) {
  N <- length(r)
  M <- sum(r)
  observed <- which(as.logical(r))
  
  # function to minimize
  f <- function(x) {
    -(1 / N) * sum(as.vector(log(1 + t(x) %*% ghat[, observed])))
  }
  
  # step 0
  rhohat <- matrix(rep(0, nrow(ghat)), ncol = 1)
  l <- 0
  tau <- 1
  
  converged <- F
  for (j in 1:maxit) {
    # step 1
    delta1 <- rep(0, nrow(ghat))
    for (i in observed) {
      tmp <- ghat[, i] / as.numeric(1 + t(rhohat) %*% ghat[, i])
      delta1 <- delta1 + tmp
    }
    delta1 <- as.matrix(delta1)
    
    sumouter <- matrix(0, nrow = nrow(ghat), ncol = nrow(ghat))
    for (i in observed) {
      tmp <- (ghat[, i] %*% t(ghat[, i])) /
        as.numeric(1 + t(rhohat) %*% ghat[, i])^2
      sumouter <- sumouter + tmp
    }
    delta2 <- solve(-sumouter,tol=1e-20) %*% delta1
    
    if (base::norm(delta2, type = 'F') < eps) {
      converged <- T
      break
    }
    
    # step 2
    condA <- condB <- F
    maxit2 <- 0 # in case of a runaway while loop
    while (maxit2 <= 10000) {
      delta <- tau * delta2
      rhotmp <- rhohat - delta
      condA <- sum((1 + t(rhotmp) %*% ghat[, observed]) <= 0) == 0
      if(condA) {condB <- f(rhotmp) > f(rhohat)}
      ifelse(condA && !condB, break, tau <- tau/2)
      maxit2 <- maxit2 + 1
      if (maxit2 == 10000)
        stop("step 2 conditions not met after 10,000 restarts")
    }
    
    # step 3
    rhohat <- rhohat - delta
    l <- l + 1
    tau <- 1 # "numerical behavior not affected with this instead" - Han (2014)
    # also seems to converge faster with tau <- 1
  }
  
  if (!converged)
    warning(paste0('Search for rhohat did not converge with eps = ', eps,
                   ', maxit = ', maxit))
  
  out <- list(
    call = call,
    rhohat = rhohat,
    converged = converged,
    args = list(
      maxit = maxit,
      eps = eps
    )
  )
  return(out)
}
