library(tidyverse)
library(MASS)
library(svglite)
library(scales)

source("scripts/data_setup.R")
source("scripts/colors.R")
source("scripts/contrasts.R") 

scenarios <- data.frame(
  PopSpecies = factor(c("B", "M", "T", "B", "M", "T"), levels = c("B", "M", "T")), 
  PhenotypeCageSpecies = factor(c("B", "B", "B", "T", "T", "T"), levels = c("B", "M", "T")),
  trays = 1 
)


set.seed(123)
n_boot <- 1000
results_list <- list()
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

for(i in 1:n_boot) {
  boot_data <- tab[sample(nrow(tab), replace = TRUE), ]
  
  boot_model <- tryCatch({
    glm.nb(numFlies ~ PopSpecies + PhenotypeCageSpecies + PopSpecies:PhenotypeCageSpecies + offset(log(trays)),
           contrasts = list(PopSpecies = contrast_matrix, PhenotypeCageSpecies = contrast_matrix), 
           data = boot_data)
  }, error = function(e) return(NULL))
  
  if(is.null(boot_model)) next 
  
  preds <- predict(boot_model, newdata = scenarios, type = "response")
  w_B_AA <- preds[1]; w_B_Aa <- preds[2]; w_B_aa <- preds[3]
  w_T_AA <- preds[4]; w_T_Aa <- preds[5]; w_T_aa <- preds[6]
  
  # niche size (c)
  total_prod_B <- sum(preds[1:3])
  total_prod_T <- sum(preds[4:6])
  val_c <- total_prod_B / (total_prod_B + total_prod_T)
  
  # asymmetry (k)
  s_B <- (w_B_AA - w_B_aa) / w_B_AA
  s_T <- (w_T_aa - w_T_AA) / w_T_aa
  
  k_val <- NA
  if(s_B > 0.001 & s_T > 0.001) k_val <- s_T / s_B
  
  invade_A <- (val_c * (w_B_Aa/w_B_aa)) + ((1-val_c) * (w_T_Aa/w_T_aa)) > 1
  invade_a <- (val_c * (w_B_Aa/w_B_AA)) + ((1-val_c) * (w_T_Aa/w_T_AA)) > 1
  
  results_list[[i]] <- data.frame(Iteration = i, c = val_c, k = k_val, Stable = (invade_A & invade_a),
                                  Invade_A = invade_A, Invade_a = invade_a)
  setTxtProgressBar(pb, i)
}
close(pb)

sim_data <- do.call(rbind, results_list) %>% filter(!is.na(k))

probs <- sim_data %>% summarise(P_Stable = mean(Stable), P_A = mean(Invade_A), P_a = mean(Invade_a))

plot_data <- sim_data %>% filter(k < 15)

cloud_stats <- plot_data %>%
  summarise(
    c_med = median(c), 
    k_med = median(k),
    c_min = quantile(c, 0.25), 
    c_max = quantile(c, 0.75),
    k_min = quantile(k, 0.25), 
    k_max = quantile(k, 0.75)
  )

ymax_plot <- max(4.0, ceiling(cloud_stats$k_max * 1.1)) 

calc_s <- 0.24 
calculate_levene_q <- function(c, s_base, k_asymmetry) {
  s1 <- s_base; s2 <- min(s_base * k_asymmetry, 10) 
  w1_AA <- 1 + 2*s1; w1_Aa <- 1 + s1; w1_aa <- 1
  w2_AA <- 1;        w2_Aa <- 1 + s2; w2_aa <- 1 + 2*s2
  invade_A <- (c * (w1_Aa / w1_aa)) + ((1 - c) * (w2_Aa / w2_aa))
  invade_a <- (c * (w1_Aa / w1_AA)) + ((1 - c) * (w2_Aa / w2_AA))
  if (invade_A > 1 && invade_a > 1) {
    delta_func <- function(q) {((c * s1) / (1 + 2 * s1 * q)) - (((1 - c) * s2) / (1 + 2 * s2 * (1 - q)))}
    return(tryCatch(uniroot(delta_func, c(0.0001, 0.9999))$root, error=function(e) NA))
  } else { if (invade_A <= 1) return(0.0); if (invade_a <= 1) return(1.0) }
  return(NA)
}

grid_data <- expand.grid(Availability_c = seq(0.01, 0.99, length.out = 150), Selection_k = seq(0.1, ymax_plot, length.out = 150))
grid_data$q_hat <- mapply(calculate_levene_q, grid_data$Availability_c, calc_s, grid_data$Selection_k)

subtitle_text <- paste0("Prob. of Stability: ", round(probs$P_Stable * 100, 1), "%")

p <- ggplot() +
  geom_tile(data = grid_data, aes(x = Availability_c, y = Selection_k, fill = q_hat)) +
  scale_fill_gradient2(low = "#EDB72D", mid = "#9BAB96", high = "#499FFF", midpoint = 0.5, limits = c(0, 1)) +
  geom_contour(data = grid_data, aes(x = Availability_c, y = Selection_k, z = q_hat), breaks = c(0.01, 0.99), color = "black", linetype = "dashed", size = 0.5) +
  
  geom_point(data = plot_data, aes(x = c, y = k), color = "white", alpha = 0.2, size = 1) +
  
  geom_errorbar(data = cloud_stats, aes(x = c_med, ymin = k_min, ymax = k_max), color = "black", width = 0.02, size = 0.5) +
  geom_errorbar(data = cloud_stats, aes(y = k_med, xmin = c_min, xmax = c_max), color = "black", width = 0.05, size = 0.5) +
  geom_point(data = cloud_stats, aes(x = c_med, y = k_med), color = "black", size = 2, shape = 18) +
  
  coord_cartesian(ylim = c(0, ymax_plot), xlim = c(0.0, 1.0), expand = FALSE) +
  labs(title = "Stability Analysis (Median + IQR)",
       subtitle = subtitle_text,
       caption = "Bars represent the Interquartile Range (25th-75th percentile)",
       x = "Effective Niche Size (c)", y = "Selection Asymmetry (k)") +
  theme_classic(base_size = 14) + theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), legend.position = "right")

print(p)
ggsave("levene_stability_bootstrap_IQR.svg", p, width = 9, height = 8, scale = 0.5)