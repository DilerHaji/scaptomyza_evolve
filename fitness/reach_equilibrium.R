library(tidyverse)
library(MASS)
library(svglite)
library(scales)

source("scripts/data_setup.R")
source("scripts/colors.R")
source("scripts/contrasts.R")

tab$PopSpecies <- factor(tab$PopSpecies, levels = c("B", "M", "T"))
tab$PhenotypeCageSpecies <- factor(tab$PhenotypeCageSpecies, levels = c("B", "M", "T"))

main_model <- glm.nb(numFlies ~ PopSpecies + PhenotypeCageSpecies +
                       PopSpecies:PhenotypeCageSpecies + offset(log(trays)),
                     contrasts = list(PopSpecies = contrast_matrix,
                                      PhenotypeCageSpecies = contrast_matrix),
                     data = tab)

scenarios <- data.frame(
  PopSpecies = factor(c("B", "M", "T", "B", "M", "T"), levels = c("B", "M", "T")),
  PhenotypeCageSpecies = factor(c("B", "B", "B", "T", "T", "T"), levels = c("B", "M", "T")),
  trays = 1
)
preds <- predict(main_model, newdata = scenarios, type = "response")

preds_v <- as.numeric(preds)  
w <- list(
  B = c(AA = preds_v[1], Aa = preds_v[2], aa = preds_v[3]),
  T = c(AA = preds_v[4], Aa = preds_v[5], aa = preds_v[6])
)

max_w <- max(unlist(w))
w$B <- w$B / max_w
w$T <- w$T / max_w

c_design <- 0.5
W_AA <- c_design * w$B["AA"] + (1 - c_design) * w$T["AA"]
W_Aa <- c_design * w$B["Aa"] + (1 - c_design) * w$T["Aa"]
W_aa <- c_design * w$B["aa"] + (1 - c_design) * w$T["aa"]

if ((W_Aa > W_AA) & (W_Aa > W_aa)) {
  denom <- 2 * W_Aa - W_AA - W_aa
  q_dem_eq <- (W_Aa - W_aa) / denom
  cat("Dempster predicted equilibrium q* =", round(q_dem_eq, 4), "\n")
} else {
  q_dem_eq <- NA
  cat("Dempster predicts no interior equilibrium\n")
}

step_dempster <- function(p, w_list, c_val = 0.5) {
  wbar_B <- p^2 * w_list$B["AA"] + 2 * p * (1 - p) * w_list$B["Aa"] + (1 - p)^2 * w_list$B["aa"]
  wbar_T <- p^2 * w_list$T["AA"] + 2 * p * (1 - p) * w_list$T["Aa"] + (1 - p)^2 * w_list$T["aa"]
  W0     <- c_val * wbar_B + (1 - c_val) * wbar_T

  num <- c_val * p * (p * w_list$B["AA"] + (1 - p) * w_list$B["Aa"]) +
         (1 - c_val) * p * (p * w_list$T["AA"] + (1 - p) * w_list$T["Aa"])
  as.numeric(num / W0)
}

step_levene <- function(p, w_list, c_val = 0.5) {
  wbar_B <- p^2 * w_list$B["AA"] + 2 * p * (1 - p) * w_list$B["Aa"] + (1 - p)^2 * w_list$B["aa"]
  wbar_T <- p^2 * w_list$T["AA"] + 2 * p * (1 - p) * w_list$T["Aa"] + (1 - p)^2 * w_list$T["aa"]

  pB <- p * (p * w_list$B["AA"] + (1 - p) * w_list$B["Aa"]) / wbar_B
  pT <- p * (p * w_list$T["AA"] + (1 - p) * w_list$T["Aa"]) / wbar_T

  as.numeric(c_val * pB + (1 - c_val) * pT)
}

simulate <- function(step_fn, p0, n_gen, w_list, c_val = 0.5) {
  p <- numeric(n_gen + 1)
  p[1] <- p0
  for (g in 1:n_gen) p[g + 1] <- step_fn(p[g], w_list, c_val)
  data.frame(generation = 0:n_gen, p = p)
}

n_gen_long <- 200
p0_values  <- c(0.2, 0.5, 0.8)

traj_list <- list()
for (p0 in p0_values) {
  tr_d <- simulate(step_dempster, p0, n_gen_long, w)
  tr_d$model <- "Dempster (hard)"
  tr_d$p0    <- p0
  tr_l <- simulate(step_levene, p0, n_gen_long, w)
  tr_l$model <- "Levene (soft)"
  tr_l$p0    <- p0
  traj_list[[paste0("d_", p0)]] <- tr_d
  traj_list[[paste0("l_", p0)]] <- tr_l
}
traj_df <- do.call(rbind, traj_list)
traj_df$p0_label <- paste0("p[0] == ", traj_df$p0)

eq_vals <- traj_df %>%
  group_by(model, p0) %>%
  summarise(p_eq = p[generation == n_gen_long], .groups = "drop")


gen_to_90 <- traj_df %>%
  left_join(eq_vals, by = c("model", "p0")) %>%
  group_by(model, p0) %>%
  mutate(dev_0   = abs(p[1] - p_eq),
         dev_t   = abs(p - p_eq),
         closed  = ifelse(dev_0 > 0, 1 - dev_t / dev_0, 1)) %>%
  summarise(t_90 = {
    hits <- which(closed >= 0.9)
    if (length(hits) == 0) NA_integer_ else generation[min(hits)]
  }, .groups = "drop")



set.seed(123)
n_boot   <- 100
n_gen_bs <- 200
p0_bs    <- 0.5

boot_traj <- list()
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
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

  preds_b  <- predict(boot_model, newdata = scenarios, type = "response")
  preds_bv <- as.numeric(preds_b)
  w_b <- list(
    B = c(AA = preds_bv[1], Aa = preds_bv[2], aa = preds_bv[3]) / max(preds_bv),
    T = c(AA = preds_bv[4], Aa = preds_bv[5], aa = preds_bv[6]) / max(preds_bv)
  )
  tr_d <- simulate(step_dempster, p0_bs, n_gen_bs, w_b)
  tr_d$iter <- i; tr_d$model <- "Dempster (hard)"
  tr_l <- simulate(step_levene,   p0_bs, n_gen_bs, w_b)
  tr_l$iter <- i; tr_l$model <- "Levene (soft)"
  boot_traj[[paste0("d_", i)]] <- tr_d
  boot_traj[[paste0("l_", i)]] <- tr_l
  setTxtProgressBar(pb, i)
}
close(pb)

boot_df <- do.call(rbind, boot_traj)

boot_summary <- boot_df %>%
  group_by(model, generation) %>%
  summarise(p_median = median(p),
            p_lo     = quantile(p, 0.05),
            p_hi     = quantile(p, 0.95),
            .groups = "drop")



p_traj <- ggplot(traj_df, aes(x = generation, y = p, color = factor(p0), linetype = model)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
           fill = "grey75", alpha = 0.4) +
  annotate("text", x = 5, y = 0.97, label = "10-gen\nexperimental\nwindow",
           size = 3, fontface = "italic") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50") +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Dark2", name = expression(p[0])) +
  scale_linetype_manual(values = c("Dempster (hard)" = "solid",
                                   "Levene (soft)"   = "dashed"),
                        name = NULL) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, n_gen_long), expand = c(0, 0)) +
  labs(x = "Generation",
       y = "Allele frequency p",
       title = "Approach to equilibrium under Dempster (hard) and Levene (soft) selection",
       subtitle = paste0("Empirical fitness architecture (reciprocal assay, s ~ 0.24). ",
                         "Our experiment samples the shaded 10-generation window.")) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "right")

print(p_traj)
ggsave("reach_equilibrium_trajectories.svg", p_traj, width = 9, height = 5, scale = 0.9)
ggsave("reach_equilibrium_trajectories.png", p_traj, width = 9, height = 5, scale = 0.9, dpi = 300)



p_env <- ggplot(boot_summary, aes(x = generation, y = p_median,
                                  color = model, fill = model)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
           fill = "grey75", alpha = 0.4) +
  annotate("text", x = 5, y = 0.97, label = "10-gen\nwindow",
           size = 3, fontface = "italic") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50") +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.25, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Dempster (hard)" = "#EDB72D",
                                "Levene (soft)"   = "#499FFF")) +
  scale_fill_manual(values  = c("Dempster (hard)" = "#EDB72D",
                                "Levene (soft)"   = "#499FFF")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, n_gen_bs), expand = c(0, 0)) +
  labs(x = "Generation",
       y = "Allele frequency p (bootstrap median, 90% envelope)",
       title = "Bootstrap uncertainty in approach trajectory",
       subtitle = paste0("Starting from p0 = 0.5. ",
                         n_boot, " non-parametric bootstrap replicates of empirical fitness.")) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "right",
        legend.title = element_blank())

print(p_env)
ggsave("reach_equilibrium_envelope.svg", p_env, width = 9, height = 5, scale = 0.9)
ggsave("reach_equilibrium_envelope.png", p_env, width = 9, height = 5, scale = 0.9,
       dpi = 300)



closure_df <- traj_df %>%
  left_join(eq_vals, by = c("model", "p0")) %>%
  mutate(dev_0  = abs(p[generation == 0][1] - p_eq),
         dev_t  = abs(p - p_eq)) %>%
  group_by(model, p0) %>%
  mutate(dev_0  = abs(first(p[generation == 0]) - first(p_eq)),
         closed = ifelse(dev_0 > 1e-9, 1 - dev_t / dev_0, 1)) %>%
  ungroup() %>%
  filter(p0 == 0.5, generation > 0)  # pick one initial condition for clarity

p_clos <- ggplot(closure_df, aes(x = generation, y = closed, color = model)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
           fill = "grey75", alpha = 0.4) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  annotate("text", x = n_gen_long, y = 0.92, label = "90% closure",
           hjust = 1, size = 3, color = "grey30") +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Dempster (hard)" = "#EDB72D",
                                "Levene (soft)"   = "#499FFF")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 200),
                     limits = c(1, n_gen_long)) +
  labs(x = "Generation (log scale)",
       y = "Fraction of initial deviation from equilibrium closed",
       title = "Transient vs equilibrium regime",
       subtitle = paste0("Our 10-generation experiment (shaded) closes only a fraction ",
                         "of the deviation to the predicted equilibrium.")) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "right",
        legend.title = element_blank())

print(p_clos)
ggsave("reach_equilibrium_closure.svg", p_clos, width = 9, height = 5, scale = 0.9)
ggsave("reach_equilibrium_closure.png", p_clos, width = 9, height = 5, scale = 0.9, dpi = 300)

closure_at_10 <- closure_df %>%
  filter(generation == 10) %>%
  dplyr::select(model, closed)

