library(tidyverse)
library(grf)
library(doParallel)
library(ggallin)
library(latex2exp)

source("../plot_theme.R")
source("../trimming_functions.R")
set.seed(1)

args = commandArgs(trailingOnly=TRUE)

SAMPLE_FRACTION = 0.5
NUM_TREES = 2000
DATASET = args[1]

#load data
practice_year = read.csv(sprintf("practice_year_%s.csv", DATASET),
                         sep = ",") %>% filter(year == 3)
patient_year = read.csv(sprintf("patient_year_%s.csv", DATASET),
                        sep = ",") %>% filter(year == 3)
practice = read.csv(sprintf("practice_%s.csv", DATASET),
                    sep = ",")
patient = read.csv(sprintf("patient_%s.csv", DATASET),
                   sep = ",") %>% filter(id.patient %in% patient_year$id.patient)

subset_frac = 0.2
subset_ids = sample(patient$id.patient, size = round(subset_frac * nrow(patient)))
patient = patient %>% filter(id.patient %in% subset_ids)
patient_year = patient_year %>% filter(id.patient %in% subset_ids)

#assemble everything in one data frame
df = left_join(patient, patient_year, by = "id.patient") %>%
  left_join(practice_year %>% dplyr::select(id.practice, Z), by = "id.practice")

practice_df = left_join(practice_year, practice, by = "id.practice") %>%
  dplyr::select(-c(year, post, n.patients))

#Fit practice level propensity
X_practice = data.matrix(practice_df %>% dplyr::select(X1, X2, X3, X4, X5, X6, X7, X8, X9))
Z_practice = practice_df$Z
practice_df$e_hat_practice = predict(regression_forest(X_practice, Z_practice))$predictions

df = df %>%
  left_join(practice_df %>% dplyr::select(id.practice, e_hat_practice), by = "id.practice")

X = data.matrix(df %>% dplyr::select(V1, V2, V3, V4, V5))
Y = df$Y
Z = df$Z

start = Sys.time()
df_trim = compute_trimming_df(X, Y, Z, df$id.patient, crossfit = F)
df_trim_cross = compute_trimming_df(X, Y, Z, df$id.patient, crossfit = T)
finish = Sys.time()
print(finish - start)

colnames(df_trim) = c("id.patient", "W", "Y", "mu0", "mu1", "sigma0", "sigma1", "e_hat", "k")
colnames(df_trim_cross) = c("id.patient", "W", "Y", "mu0", "mu1", "sigma0", "sigma1", "e_hat", "k")

df_trim = df_trim %>% left_join(df %>% dplyr::select(id.patient, e_hat_practice), by = "id.patient")
df_trim$e_hat = df_trim$e_hat_practice

df_trim_cross = df_trim_cross %>% left_join(df %>% dplyr::select(id.patient, e_hat_practice), by = "id.patient")
df_trim_cross$e_hat = df_trim_cross$e_hat_practice

df_trim = df_trim %>% mutate(k = sigma1/e_hat + sigma0/(1-e_hat)) %>% arrange(id.patient)
df_trim_cross = df_trim_cross %>% mutate(k = sigma1/e_hat + sigma0/(1-e_hat)) %>% arrange(id.patient)

delta_grid = c(1, 0.9, 0.8, 0.7)

result_df = foreach (d = delta_grid, .combine = "rbind") %do% {
  results = mixed_trimmed_estimator_delta(df_trim, df_trim_cross, d)
  c(d, results[1], results[2], results[1] - 1.96*results[2], results[1] + 1.96*results[2])
} %>% data.frame()

colnames(result_df) = c("sample_frac", "tau_hat", "se_hat", "CI_lower", "CI_upper")
rownames(result_df) = NULL
print(result_df)

tau_hat = trimmed_estimator_delta(df_trim, 1)[1]

registerDoSEQ()

#Run bootstrap

B = 1000
bootstrap_statistics = foreach (b = 1:B, .combine = "rbind") %do% {
  print(b)

  boot_idx = sample(1:nrow(df), size = nrow(df), replace = T)
  df_trim_boot = df_trim[boot_idx, ]
  df_trim_boot_cross = df_trim_cross[boot_idx, ]

  tau_hats = numeric(length(delta_grid))
  se_hats = numeric(length(delta_grid))

  for (i in seq_along(delta_grid)) {
    estimates = mixed_trimmed_estimator_delta(df_trim_boot, df_trim_boot_cross, delta_grid[i])

    tau_hats[i] = estimates[1]
    se_hats[i] = estimates[2]
  }
  (tau_hats - result_df$tau_hat) / se_hats
} %>% data.frame()


#Compute simultaneous confidence intervals
q = quantile(apply(abs(bootstrap_statistics), 1, max), 0.95)
result_df$CIlower_corrected = result_df$tau_hat - result_df$se_hat*q
result_df$CIupper_corrected = result_df$tau_hat + result_df$se_hat*q

colnames(result_df) = c("sample_frac", "tau_hat", "se_hat",
                        "CIlower_naive", "CIupper_naive",
                        "CIlower_corrected", "CIupper_corrected")

#Update pointwise intervals to use bootstrap quantiles
bootstrap_quantiles = sapply(1:4,
                             function (x) quantile(abs(bootstrap_statistics[,x]), 0.95))

result_df$CIlower_naive = result_df$tau_hat - bootstrap_quantiles*result_df$se_hat
result_df$CIupper_naive = result_df$tau_hat + bootstrap_quantiles*result_df$se_hat

#Plot confidence intervals
result_df %>%
  pivot_longer(cols = c("CIlower_naive", "CIupper_naive",
                        "CIlower_corrected", "CIupper_corrected"),
               names_to = c("side", "type"),
               names_sep = "_") %>%
  pivot_wider(names_from = c("side"), values_from = c("value")) %>%
  ggplot(aes(x = as.factor(sample_frac))) +
  geom_errorbar(aes(ymin = CIlower, ymax = CIupper, color = type),
                width = 0.25) +
  geom_point(aes(y = tau_hat)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  labs(x = TeX("Fraction of sample, $1-\\delta$"),
       y = "Medicare spending per person, $",
       color = "Method") +
  ylim(-80, 75) + 
  scale_color_manual(values = color_map, labels = c("Simultaneous CI", "Pointwise CI")) +
  theme_custom()

ggsave(sprintf("acic_corrected_cis_%s.pdf", DATASET),
       width = 4, height = 3)

#Plot AIPW scores vs. k
data.frame(df_trim) %>%
  sample_frac(0.5) %>%
  ggplot(aes(x = k,
             y = mu1 + W*(Y - mu1)/e_hat - mu0 - (1-W)*(Y - mu0)/(1-e_hat))) +
  geom_point(size = 0.4) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = pseudolog10_trans,
                     breaks = c(-100000, -10000, -1000, -100, -10, 0,
                                10, 100, 1000, 10000, 10000, 100000)) +
  labs(x = TeX("\\hat{k}"), y = "AIPW score") +
  theme_custom()
ggsave(sprintf("acic_heterogeneity_%s.pdf", DATASET), 
       width = 4, height = 3, useDingbats = T)

