library(grf)
library(doParallel)
library(logger)

compute_trimming_df = function(X, Y, W, id, crossfit = F, K = 5){
  if (!crossfit){
    df_trim = fit_regressions(X, Y, W, X, Y, W, id)
    return(df_trim)
  }
  else {
    fold_id = sample(1:K, nrow(X), replace = T)
    
    registerDoParallel(cores = 4)
    df_trim = foreach (idx = 1:K, .combine = "rbind") %dopar% {
      # log_info('Fitting fold {idx}')

      fit_regressions(X[fold_id != idx,], Y[fold_id != idx], W[fold_id != idx], 
                      X[fold_id == idx,], Y[fold_id == idx], W[fold_id == idx], id[fold_id == idx])    
      } %>% data.frame()
    
    return(df_trim)
  }
}

fit_regressions = function(X_train, Y_train, W_train, X_test, Y_test, W_test, id){

  #Estimate conditional means and variances
  print("control mean and var")
  mu0 = predict(regression_forest(X_train, Y_train, sample.weights = 1-W_train), X_test)
  sigma0 = predict(regression_forest(X_train, Y_train^2, sample.weights = 1-W_train), X_test)-mu0^2

  print("Treat mean and var")
  mu1 = predict(regression_forest(X_train, Y_train, sample.weights = W_train), X_test)
  sigma1 = predict(regression_forest(X_train, Y_train^2, sample.weights = W_train), X_test)-mu1^2
  
  print("propensity")
  e_hat = predict(regression_forest(X_train, W_train), X_test)
  
  #Assemble data frame for trimming 
  df_trim = data.frame(id, W_test, Y_test, mu0, mu1,
                       sigma0, sigma1, e_hat)
  
  colnames(df_trim) = c("id", "W", "Y", "mu0", "mu1", "sigma0", "sigma1", "e_hat")
  df_trim = df_trim %>% mutate(k = sigma1/e_hat + sigma0/(1-e_hat))
  
  return(df_trim)
}


#Compute estimator when trimming based on mask
trimmed_estimator = function(df_trim, mask){

  aipw_scores = with(df_trim, mu1 + W*(Y - mu1)/e_hat - mu0 - (1-W)*(Y - mu0)/(1-e_hat))
  
  tau_hat = sum(aipw_scores*mask)/sum(mask)
  # se_hat = sqrt(mean(df_trim$k*mask)/mean(mask)^2)
  se_hat = sqrt(var(aipw_scores*mask)/mean(mask)^2)
  se_hat = se_hat/sqrt(nrow(df_trim))
  return(c(tau_hat, se_hat))
}

#Trim delta fraction of data
trimmed_estimator_delta = function(df_trim, delta){
  gamma = quantile(df_trim$k, delta)
  mask = df_trim$k <= gamma
  return(trimmed_estimator(df_trim, mask))
}


#Run previous function over a grid of delta values
compute_trimmed_estimates = function(df_trim, delta_grid){
  result_df = foreach (d = delta_grid, .combine = "rbind") %do% {
    results = trimmed_estimator_delta(df_trim, d)
    c(d, results[1], results[2], results[1] - 1.96*results[2], results[1] + 1.96*results[2])
  } %>% data.frame()
  
  colnames(result_df) = c("sample_frac", "tau_hat", "se_hat", "CI_lower", "CI_upper")
  rownames(result_df) = NULL
  return(result_df)
}

#Functions that combine cross-fitting and not cross-fitting
mixed_trimmed_estimator_delta = function(df_trim_mask, df_trim_estim, delta){
  gamma = quantile(df_trim_mask$k, delta)
  mask = df_trim_mask$k <= gamma
  return(trimmed_estimator(df_trim_estim, mask))
}

compute_mixed_trimmed_estimates = function(df_trim_mask, df_trim_estim, delta_grid){
  result_df = foreach (d = delta_grid, .combine = "rbind") %do% {
    results = mixed_trimmed_estimator_delta(df_trim_mask, df_trim_estim, d)
    c(d, results[1], results[2], results[1] - 1.96*results[2], results[1] + 1.96*results[2])
  } %>% data.frame()
  
  colnames(result_df) = c("sample_frac", "tau_hat", "se_hat", "CI_lower", "CI_upper")
  rownames(result_df) = NULL
  return(result_df)
}


smoothing = function(t, epsilon = 0.01){
  if (epsilon == 0){
    return(t > 0)
  } else {
    return(1/(1+exp(-t/epsilon)))
  }
}

bootstrap_quantile = function(df_trim, delta_grid, alpha = 0.95, B = 1000){
  
  result_df = compute_trimmed_estimates(df_trim, delta_grid)
  
  bootstrap_statistics = foreach (b = 1:B, .combine = "rbind") %do% {
    
    boot_idx = sample(1:nrow(df_trim), size = nrow(df_trim), replace = T)
    df_trim_boot = df_trim[boot_idx, ]
    
    tau_hats = numeric(length(delta_grid))
    se_hats = numeric(length(delta_grid))
    
    for (i in seq_along(delta_grid)) {
      estimates = trimmed_estimator_delta(df_trim_boot, delta_grid[i])
      tau_hats[i] = estimates[1]
      se_hats[i] = estimates[2]
    }
    (tau_hats - result_df$tau_hat) / se_hats
  } %>% data.frame()
  
  q = quantile(apply(abs(bootstrap_statistics), 1, max), alpha)
  return(q)
}


mixed_bootstrap_quantile = function(df_trim_mask, df_trim_estim, delta_grid, alpha = 0.95, B = 1000){
  
  result_df = compute_mixed_trimmed_estimates(df_trim_mask, df_trim_estim, delta_grid)
  
  bootstrap_statistics = foreach (b = 1:B, .combine = "rbind") %do% {
    
    boot_idx = sample(1:nrow(df_trim_mask), size = nrow(df_trim_mask), replace = T)
    df_trim_mask_boot = df_trim_mask[boot_idx, ]
    df_trim_estim_boot = df_trim_estim[boot_idx, ]
    
    tau_hats = numeric(length(delta_grid))
    se_hats = numeric(length(delta_grid))
    
    for (i in seq_along(delta_grid)) {
      estimates = mixed_trimmed_estimator_delta(df_trim_mask_boot, df_trim_estim_boot, delta_grid[i])
      tau_hats[i] = estimates[1]
      se_hats[i] = estimates[2]
    }
    (tau_hats - result_df$tau_hat) / se_hats
  } %>% data.frame()
  
  q = quantile(apply(abs(bootstrap_statistics), 1, max), alpha)
  return(q)
}

  