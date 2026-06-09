library(tidyverse)
library(MASS)
library(svglite)
library(scales)

source("scripts/data_setup.R")
source("scripts/colors.R")
source("scripts/contrasts.R")

# setup 
scenarios <- data.frame(
  PopSpecies = factor(c("B", "M", "T", "B", "M", "T"), levels = c("B", "M", "T")),
  PhenotypeCageSpecies = factor(c("B", "B", "B", "T", "T", "T"), levels = c("B", "M", "T")),
  trays = 1
)

# bootstrap

# Levene (soft) invasion of allele A when rare:
# sum_i c_i * (w_Aa,i / w_aa,i)  >  1
# (Prout 1968, weighted arithmetic mean of relative heterozygote fitnesses)
#
# Dempster (hard) invasion of allele A when rare:
# sum_i c_i * w_Aa,i  >  sum_i c_i * w_aa,i
# (equivalent: arithmetic-mean fitness of Aa exceeds that of aa)

set.seed(123)
n_boot <- 1000
results_list <- list()
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

c_design <- 0.5  #  enforced niche weight in B+T cages

for (i in 1:n_boot) {
  boot_data <- tab[sample(nrow(tab), replace = TRUE), ]

  boot_model <- tryCatch({
    glm.nb(numFlies ~ PopSpecies + PhenotypeCageSpecies +
             PopSpecies:PhenotypeCageSpecies + offset(log(trays)),
           contrasts = list(PopSpecies = contrast_matrix,
                            PhenotypeCageSpecies = contrast_matrix),
           data = boot_data)
  }, error = function(e) return(NULL))

  if (is.null(boot_model)) next

  preds <- predict(boot_model, newdata = scenarios, type = "response")
  w_B_AA <- preds[1]; w_B_Aa <- preds[2]; w_B_aa <- preds[3]
  w_T_AA <- preds[4]; w_T_Aa <- preds[5]; w_T_aa <- preds[6]

  # Productivity-weighted c (as used in the original levene_test.R)
  total_prod_B <- sum(preds[1:3])
  total_prod_T <- sum(preds[4:6])
  c_prod <- total_prod_B / (total_prod_B + total_prod_T)

  # Selection coefficients
  s_B <- (w_B_AA - w_B_aa) / w_B_AA
  s_T <- (w_T_aa - w_T_AA) / w_T_aa
  k_val <- if (s_B > 0.001 & s_T > 0.001) s_T / s_B else NA

  inv_A_lev <- (c_design * (w_B_Aa / w_B_aa)) +
               ((1 - c_design) * (w_T_Aa / w_T_aa)) > 1
  inv_a_lev <- (c_design * (w_B_Aa / w_B_AA)) +
               ((1 - c_design) * (w_T_Aa / w_T_AA)) > 1
  stable_lev <- inv_A_lev & inv_a_lev

  W_Aa_arith <- c_design * w_B_Aa + (1 - c_design) * w_T_Aa
  W_AA_arith <- c_design * w_B_AA + (1 - c_design) * w_T_AA
  W_aa_arith <- c_design * w_B_aa + (1 - c_design) * w_T_aa

  inv_A_dem <- W_Aa_arith > W_aa_arith
  inv_a_dem <- W_Aa_arith > W_AA_arith
  stable_dem <- inv_A_dem & inv_a_dem

  inv_A_lev_prod <- (c_prod * (w_B_Aa / w_B_aa)) +
                    ((1 - c_prod) * (w_T_Aa / w_T_aa)) > 1
  inv_a_lev_prod <- (c_prod * (w_B_Aa / w_B_AA)) +
                    ((1 - c_prod) * (w_T_Aa / w_T_AA)) > 1
  stable_lev_prod <- inv_A_lev_prod & inv_a_lev_prod

  results_list[[i]] <- data.frame(
    Iteration   = i,
    c_prod      = c_prod,
    k           = k_val,
    # Levene soft (design c = 0.5)
    Stable_Levene_design         = stable_lev,
    Invade_A_Levene_design       = inv_A_lev,
    Invade_a_Levene_design       = inv_a_lev,
    # Dempster hard (design c = 0.5)
    Stable_Dempster_design       = stable_dem,
    Invade_A_Dempster_design     = inv_A_dem,
    Invade_a_Dempster_design     = inv_a_dem,
    # Levene soft (productivity-weighted c, matches original levene_test.R)
    Stable_Levene_productivity   = stable_lev_prod,
    # Arithmetic-mean fitnesses (for downstream equilibrium calculation)
    W_AA_arith = W_AA_arith,
    W_Aa_arith = W_Aa_arith,
    W_aa_arith = W_aa_arith
  )

  setTxtProgressBar(pb, i)
}
close(pb)

sim_data <- do.call(rbind, results_list)

probs <- sim_data %>%
  summarise(
    P_Levene_design       = mean(Stable_Levene_design),
    P_Dempster_design     = mean(Stable_Dempster_design),
    P_Levene_productivity = mean(Stable_Levene_productivity),
    # Bootstrap replicates where Levene predicts stability but Dempster does not
    P_Levene_only         = mean(Stable_Levene_design & !Stable_Dempster_design),
    # Bootstrap replicates where Dempster predicts stability but Levene does not
    P_Dempster_only       = mean(Stable_Dempster_design & !Stable_Levene_design),
    P_Both                = mean(Stable_Levene_design & Stable_Dempster_design),
    P_Neither             = mean(!Stable_Levene_design & !Stable_Dempster_design)
  )


sim_data <- sim_data %>%
  mutate(
    q_dempster = ifelse(
      Stable_Dempster_design,
      (W_Aa_arith - W_aa_arith) / (2 * W_Aa_arith - W_AA_arith - W_aa_arith),
      NA_real_
    )
  )

plot_probs <- data.frame(
  Model = factor(c("Levene (soft)", "Dempster (hard)"),
                 levels = c("Levene (soft)", "Dempster (hard)")),
  P_stable = c(probs$P_Levene_design, probs$P_Dempster_design)
)

p1 <- ggplot(plot_probs, aes(x = Model, y = P_stable, fill = Model)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", 100 * P_stable)),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Levene (soft)" = "#499FFF",
                               "Dempster (hard)" = "#EDB72D")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(y = "Bootstrap probability of protected polymorphism",
       x = NULL,
       title = "Levene vs Dempster protection at c = 0.5",
       subtitle = "Same bootstrap replicates, both invasion criteria evaluated") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

print(p1)
ggsave("levene_vs_dempster_protection.svg", p1, width = 6, height = 5, scale = 0.8)


eq_df <- sim_data %>%
  filter(!is.na(q_dempster)) %>%
  transmute(q_hat = q_dempster, Model = "Dempster (hard)")

if (nrow(eq_df) > 0) {
  p2 <- ggplot(eq_df, aes(x = q_hat, fill = Model)) +
    geom_histogram(bins = 40, color = "black") +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    scale_fill_manual(values = c("Dempster (hard)" = "#EDB72D")) +
    labs(x = "Equilibrium allele frequency q*  (Dempster, when stable)",
         y = "Bootstrap count",
         title = "Dempster-predicted equilibrium allele frequencies",
         subtitle = paste0("Only replicates where Dempster predicts stability (n = ",
                           nrow(eq_df), ")")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  print(p2)
  ggsave("dempster_equilibrium_q.svg", p2, width = 7, height = 5, scale = 0.8)
}