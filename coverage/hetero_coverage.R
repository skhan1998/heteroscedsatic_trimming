library(tidyverse)
library(grf)
library(latex2exp)
library(MASS)
library(logger)

source("../trimming_functions.R")
source("../plot_theme.R")

set.seed(1)

d = 3
tau = function(x) 2*(x[3] < 0) + (x[3] > 0) + 5*x[2]
e = function (x) 1 / (1 + 2*exp(x[2] - x[1]))
t_df = 5
t_var = t_df/(t_df - 2)
  
dgp = function(n){
  mu = rep(0, d)
  S = diag(rep(1, d))
  X = mvrnorm(n, mu, S)
  
  Y0 = 2*X[,1] - X[,2] + rnorm(n, mean = 0, sd = 1+max(0, X[,2])) + rt(n, t_df)

  Y1 = Y0 + apply(X, 1, tau)
  
  propensity = apply(X, 1, e)
  propensity = pmin(propensity, 0.95)
  propensity = pmax(propensity, 0.05)
  
  Z = rbinom(n, 1, propensity)
  Y = Z*Y1 + (1-Z)*Y0
  
  return(list("X" = X, "Y" = Y, "e" = propensity, "Z" = Z))
}

#Nuisances for above DGP
sigma_0 = function (x) 1 + max(x[2], 0) + t_var
sigma_1 = function (x) 1 + max(x[2], 0) + t_var 
k = function (x) sigma_0(x)/(1-e(x)) + sigma_1(x)/e(x)

estimate_gamma_bar = function(q, N = 1000000) { 
    mu = rep(0, d)
    S = diag(rep(1, d))
    X = mvrnorm(N, mu, S)
    k_values = apply(X, 1, k)
    return(quantile(k_values, q))
}




# n_grid = c(500, 1000, 2000, 4000, 8000, 16000, 32000)

n_grid = c(500, 1000, 2000, 4000, 8000)
N = 1000
delta_grid = c(1, 0.95, 0.9)
gamma_bar_grid = sapply(delta_grid, estimate_gamma_bar)

registerDoParallel(cores = 4)

coverage_results = foreach (n = n_grid, .combine = "rbind") %:%
  foreach (i = 1:N, .combine = "rbind") %dopar% {
    log_info('Running trial {i}')
    
    data = dgp(n)
    df_trim = compute_trimming_df(data$X, data$Y, data$Z, seq(1, n), crossfit = F) %>% arrange(id)
    df_trim_cross = compute_trimming_df(data$X, data$Y, data$Z, seq(1, n), crossfit = T) %>% arrange(id)

    result_df = compute_mixed_trimmed_estimates(df_trim, df_trim_cross, delta_grid)
    
    estimands = foreach (j = 1:length(delta_grid), .combine = "rbind") %do% {
      ites = apply(X, 1, tau)
      k_values = apply(X, 1, k)
      
      tau_A_hat = mean(ites[df_trim$k < quantile(df_trim$k, delta_grid[j])])
      tau_A_bar = mean(ites[k_values < gamma_bar_grid[j]])
      c(delta_grid[j], tau_A_hat, tau_A_bar)
    } %>% data.frame()
    
    colnames(estimands) = c("sample_frac", "tau_A_hat", "tau_A_bar")
    
    result_df = result_df %>% left_join(estimands, by = "sample_frac")
    
    result_df = result_df %>% mutate(lower_cover_bar = CI_lower < tau_A_bar,
                                     upper_cover_bar = CI_upper > tau_A_bar,
                                     lower_cover_hat = CI_lower < tau_A_hat,
                                     upper_cover_hat = CI_upper > tau_A_hat,                                     
                                     cover_bar = lower_cover_bar & upper_cover_bar,
                                     cover_hat = lower_cover_hat & upper_cover_hat,
                                     width = CI_upper - CI_lower)
    
    #Check simultaneous coverage
    q = mixed_bootstrap_quantile(df_trim, df_trim_cross, delta_grid)

    result_df = result_df %>% mutate(simul_lower_cover_bar = tau_hat - q*se_hat < tau_A_bar,
                                     simul_upper_cover_bar = tau_hat + q*se_hat > tau_A_bar,
                                     simul_cover_bar = simul_lower_cover_bar & simul_upper_cover_bar,
                                     simul_lower_cover_hat = tau_hat - q*se_hat < tau_A_hat,
                                     simul_upper_cover_hat = tau_hat + q*se_hat > tau_A_hat,
                                     simul_cover_hat = simul_lower_cover_hat & simul_upper_cover_hat)
    
    
    #Compare to Crump propensity trim 
    crump_mask = 1/(df_trim$e_hat*(1-df_trim$e_hat)) < 1/(0.9*0.1)
    crump_result = trimmed_estimator(df_trim_cross, crump_mask)
    
    result_df$crump_width = 2*1.96*crump_result[2]

    result_df$n = n
    result_df$trial = i
    
    result_df %>% dplyr::select(sample_frac, n, trial, 
                                cover_bar, cover_hat, 
                                simul_cover_bar, simul_cover_hat,
                                tau_A_bar, tau_A_hat, width,
                                crump_width)
    
  } %>% data.frame()

saveRDS(coverage_results, "mixed_coverage_results.RDS")

coverage_results %>%
  pivot_longer(cols = c(cover_bar,cover_hat)) %>%
  group_by(sample_frac, n, name) %>%
  summarise(coverage = mean(value)) %>%
  pivot_wider(names_from = n, values_from = coverage) %>%
  arrange(name, sample_frac) %>%
  write.table("coverage.csv",
            row.names = FALSE,
            sep = "&",
            quote = FALSE)

coverage_results %>%
  group_by(sample_frac, n) %>%
  summarise(avg = mean(width)) %>%
  pivot_wider(names_from = n, values_from = c(avg)) %>% 
  mutate(across(1:5, round, 3)) %>% 
  write.table("widths.csv",
            row.names = FALSE,
            sep = "&",
            quote = FALSE)

coverage_results %>%
  group_by(sample_frac, n) %>%
  summarise(avg = mean(crump_width)) %>%
  pivot_wider(names_from = n, values_from = c(avg)) %>% 
  mutate(across(1:5, round, 3)) %>% 
  write.table("crump_widths.csv",
              row.names = FALSE,
              sep = "&",
              quote = FALSE)



coverage_results %>%
  group_by(n, trial) %>%
  summarise(total_cover_bar = sum(simul_cover_bar),
            total_cover_hat = sum(simul_cover_hat)) %>%
  mutate(total_cover_bar = total_cover_bar == length(delta_grid),
         total_cover_hat = total_cover_hat == length(delta_grid)) %>%
  ungroup(trial) %>%
  summarise(simul_coverage_bar = mean(total_cover_bar),
            simul_coverage_hat = mean(total_cover_hat)) %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>% 
  write.table("simul_coverage.csv",
              row.names = FALSE,
              sep = "&",
              quote = FALSE)


