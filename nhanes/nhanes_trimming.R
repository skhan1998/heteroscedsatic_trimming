library(tidyverse)
library(grf)
library(latex2exp)

source("../trimming_functions.R")
source("../plot_theme.R")

set.seed(1)

#load data
df = read.csv("nhnes_data.txt", sep = " ") %>% 
  mutate(age_quantile = ntile(age, 10))

X_cols = c(
  "male",
  "edu.lt9",
  "edu.9to11",
  "edu.hischl",
  "edu.somecol",
  "edu.college",
  "edu.unknown",
  "income",
  "white",
  "black",
  "mexicanam",
  "otherhispan",
  "otherrace",
  "age_quantile"
)

X = data.matrix(df %>% select(X_cols))
Y = df$lead
W = df$smoking

df_trim = compute_trimming_df(X, Y, W, df$id, crossfit = F) %>% arrange(id)
df_trim_cross = compute_trimming_df(X, Y, W, df$id, crossfit = T) %>% arrange(id)

df_trim = df_trim %>%
  mutate(overlap = 1 / (e_hat * (1 - e_hat)))

hetero_objective = function(gamma) {
  mask = df_trim$k < gamma
  var_hat = mean(df_trim$k * mask) / mean(mask) ^ 2
  return(var_hat)
}

homo_objective = function(gamma) {
  mask = df_trim$overlap < gamma
  var_hat = mean(df_trim$overlap * mask) / mean(mask) ^ 2
  return(var_hat)
}

#Optimize cutoffs
gamma_hat_hetero = optimize(hetero_objective, interval = c(min(df_trim$k), max(df_trim$k)))$minimum
gamma_hat_homo = optimize(homo_objective, interval = c(min(df_trim$overlap), max(df_trim$overlap)))$minimum

df_trim = df_trim %>% mutate(
  hetero = gamma_hat_hetero - k,
  homo = gamma_hat_homo - overlap,
  naive = 1 / (0.1 * 0.9) - overlap
)

methods = c("naive", "homo", "hetero", "none")

result_df = foreach (i = 1:length(methods), .combine = "rbind") %do% {
  method = methods[i]
  
  if (method == "none") {
    mask = rep(1, nrow(df_trim))
  } else {
    mask = with(df_trim, get(method)) > 0 
  }
  
  results = trimmed_estimator(df_trim_cross, mask)
  
  c(results[1],
    results[2],
    results[1] - 1.96 * results[2],
    results[1] + 1.96 * results[2],
    mean(mask))
} %>% data.frame() %>% 
  mutate(method = methods)

colnames(result_df) = c("tau_hat",
                        "se_hat",
                        "CI_lower",
                        "CI_upper",
                        "sample_frac",
                        "method")

rownames(result_df) = NULL

#Plot confidence intervals and write table to file

method_labels = c("Optimal hetero. trim",
                  "Optimal homo. trim",
                  "0.1/0.9 rule",
                  "No trimming")

result_df %>%
  mutate(method_label = case_when(method == "hetero" ~ method_labels[1],
                                  method == "homo" ~ method_labels[2],
                                  method == "naive" ~ method_labels[3],
                                  method == "none" ~ method_labels[4]),
         method_label = factor(method_label, levels = c(method_labels[1],
                                                        method_labels[3],
                                                        method_labels[2],
                                                        method_labels[4]))) %>%
  ggplot(aes(x = sample_frac, color = as.factor(method_label))) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.05) +
  geom_point(aes(y = tau_hat)) +
  # geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
  #               position = position_dodge(width = 0.5)) +
  # geom_point(aes(y = tau_hat),
  #            position = position_dodge(width = 0.5)) +
  # geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  labs(x = TeX("Fraction of sample, $1-\\delta$"),
       y = "Blood lead level",
       color = "Method"
  ) +
  xlim(0.4, 1.1) + 
  theme_custom()

ggsave(
  "mixed_nhanes_trimming.pdf",
  width = 6,
  height = 3,
  units = "in")


write_df = result_df %>%
  mutate_at(vars(-method), function (x) round(x, 3)) %>%
  arrange(method) %>%
  select(method, tau_hat, se_hat, CI_lower, CI_upper, sample_frac) %>% 
  mutate(width = CI_upper - CI_lower)

write.table(write_df,
            "mixed_nhanes_trimming.csv",
            row.names = FALSE,
            sep = "&",
            quote = FALSE)

#Plot curve of point estimates and standard errors 

delta_grid = seq(0.1, 1.0, 0.01)

grid_result_df = foreach (d = delta_grid, .combine = "rbind") %do% {
  hetero_trim = trimmed_estimator_delta(df_trim, d)
  homo_trim = trimmed_estimator(df_trim, df_trim$overlap < quantile(df_trim$overlap, d))
  c(d, hetero_trim, homo_trim)
} %>% data.frame()

colnames(grid_result_df) = c("delta", "tau_hetero", "se_hetero", "tau_homo", "se_homo")

grid_result_df %>%
  pivot_longer(cols = c("tau_hetero", "se_hetero", "tau_homo", "se_homo"),
                                names_to = c("quantity", "type"),
                                names_sep = "_") %>%
  mutate(quantity = case_when(quantity == "se" ~ "Standard error",
                              quantity == "tau" ~ "Point estimate"),
         type = case_when(type == "hetero" ~ "Heteroscedastic",
                          type == "homo" ~ "Homoscedastic")) %>%
  ggplot(aes(x = delta, y = value, color = type)) +
  geom_line() +
  facet_wrap(quantity~., scales = "free", ncol = 1) +
  labs(x = TeX("Fraction of sample, $1-\\delta$"), y = "Value", color = "Method") +
  theme_custom()

ggsave(
  "mixed_nhanes_sequence.pdf",
  width = 6,
  height = 4,
  units = "in")



covariate_df = df %>% mutate(
  education = case_when(
    edu.lt9 == 1 ~ "lt9",
    edu.9to11 == 1 ~ "9to11",
    edu.hischl == 1 ~ "hischl",
    edu.somecol == 1 ~ "somecol",
    edu.college == 1 ~ "college",
    edu.unknown == 1 ~ "unknown"
  )
) %>%
  select(id, education) %>%
  right_join(df_trim, by = "id") %>%
  mutate(hetero_trim = if_else(k < gamma_hat_hetero, 1, 0),
         homo_trim = if_else(overlap < gamma_hat_homo, 1, 0)) %>%
  select(education, hetero_trim, homo_trim, e_hat, sigma0, sigma1)

summary_df = covariate_df %>%
  group_by(education) %>%
  summarise(count = n(),
            hetero_frac = sum(hetero_trim),
            homo_frac = sum(homo_trim),
            e_hat_avg = mean(e_hat),
            sigma_0_avg = mean(sigma0),
            sigma_1_avg = mean(sigma1)) %>% 
  mutate(hetero_frac = hetero_frac/sum(hetero_frac),
         homo_frac = homo_frac/sum(homo_frac))

colnames(summary_df) = c("Education", "Count", "Count in hetero sub.pop",
                         "Count in homo sub pop", "Avg. e_hat", "Avg. sigma0", "Avg. sigma1")

summary_df %>%
  mutate(Count = Count/sum(summary_df$Count)) %>% 
  mutate_at(vars(-"Education"), function (x) round(x, 3)) %>%
  write.table("normalized_mixed_covariate_summary.csv",
            row.names = FALSE,
            sep = "&",
            quote = FALSE)



colnames(summary_df) = c("Education", "Count", "hetero_count",
                         "homo_count", "Avg. e_hat", "Avg. sigma0", "Avg. sigma1")

summary_df %>% 
  select(hetero_count, homo_count, Education) %>% 
  filter(Education != "unknown") %>% 
  mutate(Education = factor(Education, 
                               levels = c("lt9", "9to11", "hischl", "somecol", "college"),
                               ordered = TRUE),
         hetero_count = hetero_count/sum(summary_df$hetero_count),
         homo_count = homo_count/sum(summary_df$homo_count)) %>% 
  pivot_longer(cols = c(hetero_count, homo_count)) %>% 
  ggplot(aes(x= Education, y = value, fill = name)) + 
  geom_bar(stat = "identity", position = "dodge")
  

#Plot curve of objectives 
gamma_grid = seq(min(df_trim$k), max(df_trim$k), 1)

obj_result_df = foreach (g = gamma_grid, .combine = "rbind") %do%{
  obj_het = hetero_objective(g)
  se_hat_het = trimmed_estimator(df_trim, df_trim$k < g)[2]

  obj_hom = homo_objective(g)
  se_hat_hom = trimmed_estimator(df_trim, df_trim$overlap <= g)[2]

  c(g, sqrt(obj_het/nrow(df_trim)), se_hat_het, sqrt(obj_hom/nrow(df_trim)), se_hat_hom)
} %>% data.frame()

colnames(obj_result_df) = c("gamma", "objective_het", "se_het", "objective_hom", "se_hom")

obj_result_df %>%
  pivot_longer(cols = c("objective_het", "se_het", "objective_hom", "se_hom"),
               names_to = c("quantity", "method"),
               names_sep = "_") %>%
  mutate(quantity = case_when(quantity == "objective" ~ "Heteroscedastic objective",
                              quantity == "se" ~ "Estimated std. error")) %>%
  filter(method == "het") %>%
  ggplot(aes(x = gamma, y = value, color = quantity)) +
  geom_line() +
  labs(x = TeX("Gamma, $\\gamma$"),
       y = "Objective value/Estimated std. error",
       color = "Quantity") +
  theme_custom()

ggsave(
  "nhanes_obj.pdf",
  width = 6,
  height = 3,
  units = "in")

  
  
  
  