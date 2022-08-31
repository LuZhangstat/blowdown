Block_COV <- function(coords, ind_ls, phi){
  nb = length(ind_ls)-1
  C_B <- matrix(NA, nrow = nb, ncol = nb)
  for(i in 1:nb){
    for(j in i:nb){
      #cat(i, " ", j, "\n")
      C_B[i, j] <- 0.0
      D_temp <- rdist(coords[ind_ls[i]:(ind_ls[i+1]-1), ],
                      coords[ind_ls[j]:(ind_ls[j+1]-1), ])
      C_B[i, j] <- mean(exp(-phi * D_temp))
      # for(k in ind_ls[i]:(ind_ls[i+1]-1)){
      #   for(l in ind_ls[j]:(ind_ls[j+1]-1)){
      #     C_B[i, j] <- C_B[i, j] + 
      #       exp(-phi * sqrt(sum((coords[k, ]- coords[l, ])^2)))
      #   }
      # }
      # n_pixels <- (ind_ls[i+1] - ind_ls[i]) * (ind_ls[j+1] - ind_ls[j])
      # C_B[i, j] <- C_B[i, j] / n_pixels
      if(i != j){
        C_B[j, i] <- C_B[i, j]
      }
    }
  }
  return(C_B)
}

log_lik <- function(sigmasq, tausq, y, HX, C_B, mu_beta, V_beta, Dh){
  COV = sigmasq *C_B + HX %*% tcrossprod(V_beta, HX) + tausq * diag(Dh)
  Chol_C = chol(COV)
  std_u <- forwardsolve(Chol_C, (y - HX%*%mu_beta), transpose = TRUE, 
                        upper.tri = TRUE)
  ll <- -determinant(Chol_C)$modulus - 0.5 * sum(std_u^2)
}

sample_beta_omega_no_pred <- 
  function(y, HX, C_B, sigmasq, tausq, mu_beta, V_beta, Dh){
    n = length(y) + length(mu_beta)
    invV_beta = solve(V_beta)
    invM = crossprod(Diagonal(x = sqrt(1/ (tausq * Dh))) %*%
                       cbind(HX, Diagonal(nrow(HX)))) +
      bdiag(invV_beta, chol2inv(chol(C_B)) / sigmasq)
    m = c(crossprod(HX, y/(tausq * Dh)), y/(tausq * Dh)) +
      c(invV_beta %*% mu_beta, rep(0, length(y)))
    
    cholinvM <- chol(invM)
    Mm = backsolve(cholinvM, 
                   forwardsolve(cholinvM, m, 
                                transpose = TRUE, upper.tri = TRUE))
    u <- backsolve(cholinvM, rnorm(n))
    return(Mm + u)
  }

rmvn <- function(N, mu = 0, V = matrix(1)){
  P <- length(mu)
  if(any(is.na(match(dim(V), P))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(N * P), ncol = P) %*% D + rep(mu, rep(N, P)))
}

sample_beta_omega <- function(phi_ls, sigmasq_ls, tausq_ls, coords_A, 
                              coords_AU, ind_ls_B, ind_ls_BU, HX, 
                              mu_beta, V_beta, flat_prior = FALSE){
  
  # phi_ls: vector of posterior samples of phi
  # sigmasq_ls: vector of posterior samples of sigmasq
  # tausq_ls: vector of posterior samples of tausq
  
  n_sam <- length(phi_ls)   # number of posterior samples
  nb = length(ind_ls_B)-1   # number of observed responses
  nbu = length(ind_ls_BU)-1 # number of polygons to predict
  p = length(mu_beta)       # number of predictors
  if(flat_prior){
    invV_beta = matrix(0.0, nrow = p, ncol = p)
    invVmu_beta = rep(0.0, p)
  }else{
    invV_beta = solve(V_beta)
    invVmu_beta = invV_beta %*% mu_beta
  }
  
  
  # preallocation #
  C_B_ls <- array(NA, dim = c(nb, nb, n_sam))
  C_BU_ls <- array(NA, dim = c(nb, nbu, n_sam))
  C_UU_ls <- array(NA, dim = c(nbu, nbu, n_sam))
  omega_BU_ls <- matrix(NA, nrow = nbu, ncol = n_sam) # preallocate the posterior samples of omega^u_B
  omega_B_ls <- matrix(NA, nrow = nb, ncol = n_sam)   # preallocate the posterior samples of omega_B
  beta_ls <- matrix(NA, nrow = p, ncol = n_sam)      # preallocate the posterior samples of beta
  
  
  # generate C_B_ls #
  for(i in 1:nb){
    for(j in i:nb){
      #cat(i, "\t", j, "\n")
      D_temp <- rdist(coords_A[ind_ls_B[i]:(ind_ls_B[i+1]-1), ],
                      coords_A[ind_ls_B[j]:(ind_ls_B[j+1]-1), ])
      for(k in 1:n_sam){
        C_B_ls[i, j, k] <- sigmasq_ls[[k]] * mean(exp(-phi_ls[[k]] * D_temp))
        if(i != j){
          C_B_ls[j, i, k] <- C_B_ls[i, j, k]
        }
      }
    }
  }
  
  # generate C_BU_ls #
  for(i in 1:nb){
    for(j in 1:nbu){
      #cat(i, "\t", j, "\n")
      D_temp <- rdist(coords_A[ind_ls_B[i]:(ind_ls_B[i+1]-1), ],
                      coords_AU[ind_ls_BU[j]:(ind_ls_BU[j+1]-1), ])
      for(k in 1:n_sam){
        C_BU_ls[i, j, k] <- sigmasq_ls[[k]] * mean(exp(-phi_ls[[k]] * D_temp))
      }
    }
  }
  
  # generate C_UU_ls #
  for(i in 1:nbu){
    for(j in i:nbu){
      #cat(i, "\t", j, "\n")
      D_temp <- rdist(coords_AU[ind_ls_BU[i]:(ind_ls_BU[i+1]-1), ],
                      coords_AU[ind_ls_BU[j]:(ind_ls_BU[j+1]-1), ])
      for(k in 1:n_sam){
        C_UU_ls[i, j, k] <- sigmasq_ls[[k]] * mean(exp(-phi_ls[[k]] * D_temp))
        if(i != j){
          C_UU_ls[j, i, k] <- C_UU_ls[i, j, k]
        }
      }
    }
  }
  
  # generate posterior samples of beta omega and omega^u_B
  for(k in 1:n_sam){
    Chol_CB <- chol(C_B_ls[, , k])
    
    # generate posterior samples of beta and omega
    invM = crossprod(Diagonal(x = sqrt(1/ (tausq_ls[k] * Dh))) %*% 
                       cbind(HX, Diagonal(nb))) + 
      bdiag(invV_beta, chol2inv(Chol_CB))
    m = c(crossprod(HX, y / (tausq_ls[k] * Dh)), y / (tausq_ls[k] * Dh)) + 
      c(invVmu_beta, rep(0, nb))
    
    cholinvM <- chol(invM)
    Mm = backsolve(cholinvM, 
                   forwardsolve(cholinvM, m, 
                                transpose = TRUE, upper.tri = TRUE))
    u <- backsolve(cholinvM, rnorm(nb + p))
    beta_ls[, k] = Mm[1:p] + u[1:p]
    omega_B_ls[, k] = Mm[-c(1:p)] + u[-c(1:p)]
    
    
    # generate posterior samples of ometa^u_B
    RTCBU <- forwardsolve(Chol_CB, C_BU_ls[, , k], 
                          transpose = TRUE, upper.tri = TRUE)
    RTW <- forwardsolve(Chol_CB, omega_B_ls[, k], 
                        transpose = TRUE, upper.tri = TRUE)
    
    # covariance of the conditional posterior predictive distribution
    C <- C_UU_ls[, , k] - crossprod(RTCBU) 
    # mean of the conditional posterior predictive distribution
    mpred <- c(crossprod(RTCBU, RTW))
    
    # generate posterior predictive samples of omega_BU_ls
    omega_BU_ls[, k] <- c(rnorm(nbu) %*% chol(C) + mpred)
  }
  return(list(beta_ls = beta_ls,
              omega_B_ls = omega_B_ls,
              omega_BU_ls = omega_BU_ls))
}



sample_beta_omega_h <- function(phi_ls, sigmasq_ls, tausq_ls, coords_A, 
                              coords_AU, hA, hAU, ind_ls_B, ind_ls_BU, HX, 
                              mu_beta, V_beta, flat_prior = FALSE){
  
  # phi_ls: vector of posterior samples of phi
  # sigmasq_ls: vector of posterior samples of sigmasq
  # tausq_ls: vector of posterior samples of tausq
  
  n_sam <- length(phi_ls)   # number of posterior samples
  nb = length(ind_ls_B)-1   # number of observed responses
  nbu = length(ind_ls_BU)-1 # number of polygons to predict
  p = length(mu_beta)       # number of predictors
  if(flat_prior){
    invV_beta = matrix(0.0, nrow = p, ncol = p)
    invVmu_beta = rep(0.0, p)
  }else{
    invV_beta = solve(V_beta)
    invVmu_beta = invV_beta %*% mu_beta
  }
  
  
  # preallocation #
  C_B_ls <- array(NA, dim = c(nb, nb, n_sam))
  C_BU_ls <- array(NA, dim = c(nb, nbu, n_sam))
  C_UU_ls <- array(NA, dim = c(nbu, nbu, n_sam))
  omega_BU_ls <- matrix(NA, nrow = nbu, ncol = n_sam) # preallocate the posterior samples of omega^u_B
  omega_B_ls <- matrix(NA, nrow = nb, ncol = n_sam)   # preallocate the posterior samples of omega_B
  beta_ls <- matrix(NA, nrow = p, ncol = n_sam)      # preallocate the posterior samples of beta
  
  
  # generate C_B_ls #
  for(i in 1:nb){
    for(j in i:nb){
      #cat(i, "\t", j, "\n")
      D_temp <- rdist(coords_A[ind_ls_B[i]:(ind_ls_B[i+1]-1), ],
                      coords_A[ind_ls_B[j]:(ind_ls_B[j+1]-1), ])
      for(k in 1:n_sam){
        C_B_ls[i, j, k] <- sigmasq_ls[[k]] * 
          sum((hA[ind_ls_B[i]:(ind_ls_B[i+1]-1)]) * (exp(-phi_ls[[k]] * D_temp)
           %*% hA[ind_ls_B[j]:(ind_ls_B[j+1]-1)]))
        if(i != j){
          C_B_ls[j, i, k] <- C_B_ls[i, j, k]
        }
      }
    }
  }

  
  # generate C_BU_ls #
  for(i in 1:nb){
    for(j in 1:nbu){
      #cat(i, "\t", j, "\n")
      D_temp <- rdist(coords_A[ind_ls_B[i]:(ind_ls_B[i+1]-1), ],
                      coords_AU[ind_ls_BU[j]:(ind_ls_BU[j+1]-1), ])
      for(k in 1:n_sam){
        C_BU_ls[i, j, k] <- sigmasq_ls[[k]] * 
          sum((hA[ind_ls_B[i]:(ind_ls_B[i+1]-1)]) * 
             (exp(-phi_ls[[k]] * D_temp) %*% 
             hAU[ind_ls_BU[j]:(ind_ls_BU[j+1]-1)]))
      }
    }
  }
  
  #t <- proc.time()
  # generate C_UU_ls #
  for(i in 1:nbu){
    for(j in i:nbu){
      #cat(i, "\t", j, "\n")
      D_temp <- rdist(coords_AU[ind_ls_BU[i]:(ind_ls_BU[i+1]-1), ],
                      coords_AU[ind_ls_BU[j]:(ind_ls_BU[j+1]-1), ])
      for(k in 1:n_sam){
        C_UU_ls[i, j, k] <- sigmasq_ls[[k]] * (
          sum((hAU[ind_ls_BU[i]:(ind_ls_BU[i+1]-1)]) * 
            (exp(-phi_ls[[k]] * D_temp) %*% 
              hAU[ind_ls_BU[j]:(ind_ls_BU[j+1]-1)])))
        if(i != j){
          C_UU_ls[j, i, k] <- C_UU_ls[i, j, k]
        }
      }
    }
  }
  #proc.time() - t
  
  # generate posterior samples of beta omega and omega^u_B
  for(k in 1:n_sam){
    Chol_CB <- chol(C_B_ls[, , k])
    
    # generate posterior samples of beta and omega
    invM = crossprod(Diagonal(x = sqrt(1/ (tausq_ls[k] * Dh))) %*% 
                       cbind(HX, Diagonal(nb))) + 
      bdiag(invV_beta, chol2inv(Chol_CB))
    m = c(crossprod(HX, y / (tausq_ls[k] * Dh)), y / (tausq_ls[k] * Dh)) + 
      c(invVmu_beta, rep(0, nb))
    
    cholinvM <- chol(invM)
    Mm = backsolve(cholinvM, 
                   forwardsolve(cholinvM, m, 
                                transpose = TRUE, upper.tri = TRUE))
    u <- backsolve(cholinvM, rnorm(nb + p))
    beta_ls[, k] = Mm[1:p] + u[1:p]
    omega_B_ls[, k] = Mm[-c(1:p)] + u[-c(1:p)]
    
    
    # generate posterior samples of ometa^u_B
    RTCBU <- forwardsolve(Chol_CB, C_BU_ls[, , k], 
                          transpose = TRUE, upper.tri = TRUE)
    RTW <- forwardsolve(Chol_CB, omega_B_ls[, k], 
                        transpose = TRUE, upper.tri = TRUE)
    
    # covariance of the conditional posterior predictive distribution
    C <- C_UU_ls[, , k] - crossprod(RTCBU) 
    # mean of the conditional posterior predictive distribution
    mpred <- c(crossprod(RTCBU, RTW))
    
    # generate posterior predictive samples of omega_BU_ls
    omega_BU_ls[, k] <- c(rnorm(nbu) %*% chol(C) + mpred)
  }
  return(list(beta_ls = beta_ls,
              omega_B_ls = omega_B_ls,
              omega_BU_ls = omega_BU_ls))
}



pred_sample_y <- function(beta_ls, omega_BU_ls, tausq_ls, HXU, DhU){
  # beta_ls: posterior samples of beta 
  # omega_BU_ls: posterior samples of omega on polygons for prediction
  # tausq_ls: posterior samples of tausq
  # HXU: design matrix on polygons for prediction
  # DhU: scales of noises on ploygons for predition
  
  n_sam <- length(tausq_ls)
  nbu <- nrow(omega_BU_ls)
  yU_ls <- matrix(NA, nrow = nbu, ncol = n_sam)
  
  for(i in 1:n_sam){
    noise <- DhU * rnorm(nbu) * sqrt(tausq_ls[i])
    yU_ls[, i] <- HXU %*% beta_ls[, i] + omega_BU_ls[, i] + noise
  }
  return(yU_ls)
}