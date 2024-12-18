fit_modbart <- function(X_reg, X_bart, Y, hypers = NULL, opts = NULL, 
                        length_scale = NULL,
                        x_pred, x_pred_bart, y_pred,
                        num_burn, num_thin, num_save) {
  
  if(is.null(opts)) opts <- Opts(update_sigma = FALSE)
  
  

  
  ## Initialize params
  
  grid_size <- length(y_pred)
  fitted_lm <- lm(Y ~ X_reg)
  beta         <- coef(fitted_lm)
  alpha        <- beta[1]
  beta         <- beta[-1]
  sigma        <- sigma(fitted_lm)
  alpha_probit <- 0
  if(is.null(hypers)) {
    hypers     <- Hypers(X = X_bart, Y = rnorm(nrow(X_bart)),
                         num_tree = 50, k = 1)
    hypers$sigma_hat          <- 1
    hypers$length_scale <- length_scale
    if(is.null(length_scale)) hypers$length_scale <- 2 * sd(Y) / pi
    hypers$shape_length_scale <- 1
    hypers$rate_length_scale  <- 1 / hypers$length_scale^2
  }
  
  ## Output record
  dens_array <- array(NA, c(nrow(x_pred), length(y_pred), num_save))
  reject_size <- numeric(num_save)
  beta_out <- matrix(NA, nrow = num_save, ncol = ncol(X_reg))
  alpha_out <- numeric(num_save)
  alpha_probit_out <- numeric(num_save)
  sigma_out <- numeric(num_save)
  sigma_mu <- numeric(num_save)
  length_scale <- numeric(num_save)
  s <- matrix(NA, nrow = num_save, ncol = length(hypers$group))
  

  ## Make forest
  forest <- MakeForest(hypers = hypers, opts = opts)
  
  state <- list(alpha = alpha, beta = beta, sigma = sigma, 
                alpha_probit = alpha_probit, forest = forest, 
                X_reg = X_reg, X_bart = X_bart, Y = Y)
  

  
  ## Run Loop
  for(i in 1:num_burn) {
    state <- iterate_gibbs(state)
    cat(paste("\rFinishing warmup", i, "\t\t\t"))
  }
  
  for(i in 1:num_save) {
    for(j in 1:num_thin) {
      state <- iterate_gibbs(state)
    }
    
    end_params = state$forest$get_params()
    
    alpha_probit_out[i] <- state$alpha_probit
    alpha_out[i] <- state$alpha
    beta_out[i,] <- state$beta
    sigma_out[i] <- state$sigma
    sigma_mu[i] <- end_params$sigma_mu
    length_scale[i] <- end_params$length_scale
    s[i,] <- end_params$s
    
    
    
    for(j in 1:nrow(x_pred)) {
      X_test <- matrix(rep(x_pred[j,], grid_size), 
                       nrow = grid_size, byrow = TRUE)
      X_test_bart <- matrix(rep(x_pred_bart[j,], grid_size), 
                            nrow = grid_size, byrow = TRUE)
      preds    <- state$alpha_probit + state$forest$predict(X_test_bart, y_pred)
      tilts  <- pnorm(preds)
      dens_est <- dnorm(x = y_pred, 
                        mean = state$alpha + as.numeric(X_test %*% state$beta),
                        sd = state$sigma) * tilts
      dens_est <- dens_est / sum(dens_est * (y_pred[2] - y_pred[1]))
      dens_array[j,,i] <- dens_est
    }
    
    reject_size[i] <- nrow(state$X_tot_reg)
    
    cat(paste("\rFinishing save", i, "\t\t\t"))
    
  }
  

  
  return(list(f_pred = dens_array, 
              reject_size = reject_size, 
              end_params = state$forest$get_params(), 
              alpha = alpha_out, 
              beta = beta_out, 
              sigma = sigma_out, 
              alpha_probit = alpha_probit_out, 
              sigma_mu = sigma_mu, 
              length_scale = length_scale, 
              s = s))
  
}

augment_modbart <- function(forest, Y, X_reg, X_bart,
                            alpha_lm, beta_lm, sigma_lm, 
                            alpha_probit) {
  
  Y_reject <- numeric(0)
  N <- length(Y)
  idx_remain <- 1:N
  X_tot_reg <- X_reg[-idx_remain, ,drop = FALSE]
  X_tot_bart <- X_bart[-idx_remain, , drop = FALSE]
  
  
  while(length(idx_remain) > 0) {
    num_remain <- length(idx_remain)
    mu <- alpha_lm + as.numeric(X_reg[idx_remain,,drop = FALSE] %*% beta_lm)
    omega_prop <- rnorm(n = num_remain, mean = mu, sd = sigma_lm)
    psi_prop <- forest$predict(X_bart[idx_remain,,drop=FALSE], omega_prop)
    psi_prop <- alpha_probit + as.numeric(psi_prop)
    accepted <- which(runif(num_remain) < pnorm(psi_prop))
    idx_remain <- idx_remain[-accepted]
    if(length(idx_remain) > 0) {
      X_tot_reg <- rbind(X_tot_reg, X_reg[idx_remain,,drop=FALSE])
      X_tot_bart <- rbind(X_tot_bart, X_bart[idx_remain,,drop=FALSE])
      Y_reject <- c(Y_reject, omega_prop[-accepted])
    }
  }
  
  return(list(X_tot_reg = X_tot_reg, X_tot_bart = X_tot_bart, Y_reject = Y_reject))
}

update_lm_params <- function(Y, X) {
  updating_lm <- lm(Y ~ X)
  dof <- nrow(X) - ncol(X) - 1
  SSE <- sigma(updating_lm)^2 * dof
  s_squared <- SSE / dof
  sigma <- sqrt(1 / rgamma(1, dof / 2, SSE / 2))
  betahat <- coef(updating_lm)
  varmat <- vcov(updating_lm) * sigma^2 / s_squared
  beta <- MASS::mvrnorm(1, betahat, varmat)
  alpha <- beta[1]
  beta <- beta[-1]
  
  return(list(alpha = alpha, beta = beta, sigma = sigma))
}


augment_probit <- function(forest, A, Y, X, alpha_probit) {
  psi <- forest$predict(X,Y)
  psi <- alpha_probit + as.numeric(psi)
  omega <- rtruncnorm(n = length(psi), 
                      a = ifelse(A == 0, -Inf, 0), 
                      b = ifelse(A == 0, 0, Inf),
                      mean = psi)
  
}

iterate_gibbs <- function(state) {
  
  ## Augment the data
  c(state$X_tot_reg, state$X_tot_bart, state$Y_reject) %<-% 
    augment_modbart(forest = state$forest, 
                    Y = state$Y, 
                    X_reg = state$X_reg,
                    X_bart = state$X_bart,
                    alpha_lm = state$alpha, 
                    beta_lm = state$beta, 
                    sigma_lm = state$sigma, 
                    alpha_probit = state$alpha_probit)
  
  
  ## Update regression parameters
  Y_aug <- c(state$Y, state$Y_reject)
  X_aug_reg <- rbind(state$X_reg, state$X_tot_reg)
  c(state$alpha, state$beta, state$sigma) %<-% 
    update_lm_params(Y_aug, X_aug_reg)
  
  ## Update Probit
  A <- c(rep(1, length(state$Y)), rep(0, length(state$Y_reject)))
  X_aug_bart   <- rbind(state$X_bart, state$X_tot_bart)
  omega <- augment_probit(state$forest, 
                          A, 
                          Y_aug, 
                          X_aug_bart, 
                          state$alpha_probit)
  
  ## Update for trees and alpha_probit with Normal(1,1) prior
  tmp          <- state$forest$do_gibbs(X_aug_bart, omega - state$alpha_probit, 
                                        Y_aug, X_aug_bart, Y_aug, 1)
  R            <- as.numeric(omega - tmp)
  state$alpha_probit <- rnorm(1, (1 + sum(R)) / (length(R) + 1),
                              1 / sqrt(length(R) + 1))
  
  return(state)
  
}