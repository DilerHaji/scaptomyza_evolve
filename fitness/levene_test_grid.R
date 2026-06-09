library(tidyverse)
library(MASS)      
library(svglite)
library(scales)
library(patchwork)

tab$PopSpecies <- factor(tab$PopSpecies, levels = c("B", "M", "T"))
tab$PhenotypeCageSpecies <- factor(tab$PhenotypeCageSpecies, levels = c("B", "M", "T"))

main_model <- glm.nb(numFlies ~ PopSpecies + PhenotypeCageSpecies + PopSpecies:PhenotypeCageSpecies + offset(log(trays)),
                     contrasts = list(PopSpecies = contrast_matrix, PhenotypeCageSpecies = contrast_matrix), 
                     data = tab)

anov <- anova(main_model, test = "Chisq")
p_val_interaction <- anov$`Pr(>Chi)`[which(rownames(anov) == "PopSpecies:PhenotypeCageSpecies")]
p_label <- ifelse(p_val_interaction < 0.001, "p < 0.001", paste0("p = ", round(p_val_interaction, 3)))

set.seed(123)
n_boot <- 1000 

scenarios_boot <- data.frame(
  PopSpecies = factor(c("B", "M", "T", "B", "M", "T"), levels = c("B", "M", "T")), 
  PhenotypeCageSpecies = factor(c("B", "B", "B", "T", "T", "T"), levels = c("B", "M", "T")),
  trays = 1 
)

results_list <- list()
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

for(i in 1:n_boot) {
  # resample data
  boot_data <- tab[sample(nrow(tab), replace = TRUE), ]
  
  # fit model (handle convergence errors)
  boot_model <- tryCatch({
    glm.nb(numFlies ~ PopSpecies + PhenotypeCageSpecies + PopSpecies:PhenotypeCageSpecies + offset(log(trays)),
           contrasts = list(PopSpecies = contrast_matrix, PhenotypeCageSpecies = contrast_matrix), 
           data = boot_data)
  }, error = function(e) return(NULL))
  
  if(is.null(boot_model)) next 
  
  # predict
  preds <- predict(boot_model, newdata = scenarios_boot, type = "response")
  
  # calculate c and k
  total_prod_B <- sum(preds[1:3])
  total_prod_T <- sum(preds[4:6])
  val_c <- total_prod_B / (total_prod_B + total_prod_T)
  
  w_B_AA <- preds[1]; w_B_aa <- preds[3]
  w_T_AA <- preds[4]; w_T_aa <- preds[6]
  
  s_B <- (w_B_AA - w_B_aa) / w_B_AA
  s_T <- (w_T_aa - w_T_AA) / w_T_aa
  
  k_val <- NA
  if(s_B > 0.001 & s_T > 0.001) k_val <- s_T / s_B
  
  results_list[[i]] <- data.frame(Iteration = i, c = val_c, k = k_val)
  setTxtProgressBar(pb, i)
}
close(pb)

sim_data <- do.call(rbind, results_list) %>% filter(!is.na(k))

plot_data_base <- sim_data %>% filter(k < 8)

cloud_stats_base <- plot_data_base %>%
  summarise(
    c_med = median(c), k_med = median(k),
    c_min = quantile(c, 0.25), c_max = quantile(c, 0.75),
    k_min = quantile(k, 0.25), k_max = quantile(k, 0.75)
  )


# helper
calculate_levene_q <- function(c, s_base, k_asymmetry) {
  s1 <- s_base; s2 <- min(s_base * k_asymmetry, 10) 
  w1_AA <- 1 + 2*s1; w1_Aa <- 1 + s1; w1_aa <- 1
  w2_AA <- 1;        w2_Aa <- 1 + s2; w2_aa <- 1 + 2*s2
  
  invade_A <- (c * (w1_Aa / w1_aa)) + ((1 - c) * (w2_Aa / w2_aa))
  invade_a <- (c * (w1_Aa / w1_AA)) + ((1 - c) * (w2_Aa / w2_AA))
  
  if (invade_A > 1 && invade_a > 1) {
    delta_func <- function(q) {((c * s1) / (1 + 2 * s1 * q)) - (((1 - c) * s2) / (1 + 2 * s2 * (1 - q)))}
    return(tryCatch(uniroot(delta_func, c(0.0001, 0.9999))$root, error=function(e) NA))
  } else { 
    if (invade_A <= 1) return(0.0)
    if (invade_a <= 1) return(1.0) 
  }
  return(NA)
}

# loop over s levels
s_levels <- seq(0,1,0.05)
grid_list <- list()
ann_list <- list()

for(s_val in s_levels) {
  temp_grid <- expand.grid(Availability_c = seq(0.01, 0.99, length.out = 100), 
                           Selection_k = seq(0.1, 8.0, length.out = 100))
  temp_grid$q_hat <- mapply(calculate_levene_q, temp_grid$Availability_c, s_val, temp_grid$Selection_k)
  temp_grid$s_facet <- paste0("s = ", s_val)
  grid_list[[as.character(s_val)]] <- temp_grid
  

  s1_vec <- rep(s_val, nrow(sim_data))
  s2_vec <- pmin(s_val * sim_data$k, 10)
  
  w1_Aa_w1_aa <- (1 + s1_vec) / 1
  w2_Aa_w2_aa <- (1 + s2_vec) / (1 + 2*s2_vec)
  w1_Aa_w1_AA <- (1 + s1_vec) / (1 + 2*s1_vec)
  w2_Aa_w2_AA <- (1 + s2_vec) / 1
  
  invade_A <- (sim_data$c * w1_Aa_w1_aa) + ((1 - sim_data$c) * w2_Aa_w2_aa) > 1
  invade_a <- (sim_data$c * w1_Aa_w1_AA) + ((1 - sim_data$c) * w2_Aa_w2_AA) > 1
  
  pct_stable <- mean(invade_A & invade_a)
  
  ann_list[[as.character(s_val)]] <- data.frame(
    s_facet = paste0("s = ", s_val),
    c = 0.02, k = 7.5, # Top-Left Position
    label = paste0("p(stable): ", round(pct_stable * 100, 1), "%")
  )
}

all_grids <- do.call(rbind, grid_list)
ann_df <- do.call(rbind, ann_list)

all_grids$s_facet <- factor(all_grids$s_facet, levels = paste0("s = ", s_levels))
ann_df$s_facet <- factor(ann_df$s_facet, levels = paste0("s = ", s_levels))


exp_data_faceted <- list()
stats_faceted <- list()

for(s_val in s_levels) {
  d <- plot_data_base; d$s_facet <- paste0("s = ", s_val)
  exp_data_faceted[[as.character(s_val)]] <- d
  
  s <- cloud_stats_base; s$s_facet <- paste0("s = ", s_val)
  stats_faceted[[as.character(s_val)]] <- s
}

final_cloud_data <- do.call(rbind, exp_data_faceted)
final_stats_data <- do.call(rbind, stats_faceted)
final_cloud_data$s_facet <- factor(final_cloud_data$s_facet, levels = paste0("s = ", s_levels))
final_stats_data$s_facet <- factor(final_stats_data$s_facet, levels = paste0("s = ", s_levels))


p_grid <- ggplot() +
  geom_tile(data = all_grids, aes(x = Availability_c, y = Selection_k, fill = q_hat)) +
  scale_fill_gradient2(low = "#EDB72D", mid = "#9BAB96", high = "#499FFF",  
                       midpoint = 0.5, limits = c(0, 1), name = expression(hat(q))) +
  geom_contour(data = all_grids, aes(x = Availability_c, y = Selection_k, z = q_hat), 
               breaks = c(0.01, 0.99), color = "black", linetype = "dashed", size = 0.3) +
  
  geom_point(data = final_cloud_data, aes(x = c, y = k), 
             alpha = 0.15, size = 0.8, color = "white") +
  
  geom_errorbar(data = final_stats_data, aes(x = c_med, ymin = k_min, ymax = k_max), 
                width = 0.02, size = 1, color = "black") +
  geom_errorbar(data = final_stats_data, aes(y = k_med, xmin = c_min, xmax = c_max), 
                width = 0.1, size = 1, color = "black") +
  geom_point(data = final_stats_data, aes(x = c_med, y = k_med), 
             size = 3, shape = 23, fill = "white", stroke = 1) +
  
  geom_text(data = ann_df, aes(x = c, y = k, label = label), 
            hjust = 0, vjust = 1, fontface = "bold", size = 3.5) +
  
  facet_wrap(~s_facet, ncol = 7) + 
  
  coord_cartesian(ylim = c(0, 8), xlim = c(0.0, 1.0), expand = FALSE) +
  labs(title = "",
       subtitle = "",
       x = "Effective Niche Size (c)", 
       y = "Selection Asymmetry (k)") +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey90"),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))

print(p_grid)
ggsave("Levene_Grid_Annotated.png", p_grid, width = 8, height = 4, scale = 1.5)
