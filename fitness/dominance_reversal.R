# Computes the population-level "dominance coefficients" h_B and h_T from the
# reciprocal performance fitness assay, treating B-evolved = AA,
# B+T-evolved = Aa, T-evolved = aa as a single-locus fiction.

library(tidyverse)
library(MASS)

source("scripts/data_setup.R")
source("scripts/contrasts.R")

scenarios <- data.frame(
  PopSpecies = factor(c("B", "M", "T", "B", "M", "T"), levels = c("B", "M", "T")),
  PhenotypeCageSpecies = factor(c("B", "B", "B", "T", "T", "T"), levels = c("B", "M", "T")),
  trays = 1
)

compute_h <- function(preds) {
  w_B_AA <- preds[1]; w_B_Aa <- preds[2]; w_B_aa <- preds[3]
  w_T_AA <- preds[4]; w_T_Aa <- preds[5]; w_T_aa <- preds[6]
  h_B <- (w_B_Aa - w_B_aa) / (w_B_AA - w_B_aa)       # toward AA (favored on B)
  h_T <- (w_T_Aa - w_T_AA) / (w_T_aa - w_T_AA)       # toward aa (favored on T)
  c(h_B = h_B, h_T = h_T,
    reversal_strict = (h_B > 0.5) & (h_T > 0.5),
    convex_fitness_hull = ((w_B_Aa + w_T_Aa) / 2) >
                          ((w_B_AA + w_T_AA + w_B_aa + w_T_aa) / 4))
}

# Point estimate from full data
main_model <- glm.nb(numFlies ~ PopSpecies + PhenotypeCageSpecies +
                       PopSpecies:PhenotypeCageSpecies + offset(log(trays)),
                     contrasts = list(PopSpecies = contrast_matrix,
                                      PhenotypeCageSpecies = contrast_matrix),
                     data = tab)
preds_main <- as.numeric(predict(main_model, newdata = scenarios, type = "response"))
h_point <- compute_h(preds_main)

# Bootstrap
set.seed(123)
n_boot <- 1000
boot_h <- matrix(NA_real_, n_boot, 4,
                 dimnames = list(NULL, c("h_B", "h_T", "reversal_strict", "convex_fitness_hull")))

pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
for (i in 1:n_boot) {
  boot_data <- tab[sample(nrow(tab), replace = TRUE), ]
  boot_model <- tryCatch({
    glm.nb(numFlies ~ PopSpecies + PhenotypeCageSpecies +
             PopSpecies:PhenotypeCageSpecies + offset(log(trays)),
           contrasts = list(PopSpecies = contrast_matrix,
                            PhenotypeCageSpecies = contrast_matrix),
           data = boot_data)
  }, error = function(e) NULL)
  if (is.null(boot_model)) next
  preds_b <- predict(boot_model, newdata = scenarios, type = "response")
  boot_h[i, ] <- compute_h(as.numeric(preds_b))
  setTxtProgressBar(pb, i)
}
close(pb)

boot_df <- as.data.frame(boot_h) %>% drop_na()

# Plot
library(ggplot2)
library(svglite)

p1 <- ggplot(boot_df, aes(x = h_B, y = h_T)) +
  annotate("rect", xmin = 0.5, xmax = Inf, ymin = 0.5, ymax = Inf,
           fill = "#EDB72D", alpha = 0.25) +
  annotate("text", x = 1.05, y = 1.05, hjust = 1, vjust = 1,
           label = "dominance reversal\n(both h > 0.5)", size = 3.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_point(aes(x = h_point["h_B"], y = h_point["h_T"]),
             data = data.frame(x = 1), color = "red", size = 3, shape = 18) +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 1.2)) +
  labs(x = expression(h[B] ~~ "(dominance toward AA on host B)"),
       y = expression(h[T] ~~ "(dominance toward aa on host T)"),
       title = "Dominance-reversal diagnostic",
       subtitle = paste0("Bootstrap n = ", nrow(boot_df),
                         ". Red diamond = point estimate. ",
                         "P(both h > 0.5) = ",
                         round(mean(boot_df$reversal_strict), 3))) +
  theme_classic(base_size = 13) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

print(p1)
ggsave("dominance_reversal_bootstrap.svg", p1, width = 7, height = 7, scale = 0.8)
ggsave("dominance_reversal_bootstrap.png", p1, width = 7, height = 7, scale = 0.8,
       dpi = 300)