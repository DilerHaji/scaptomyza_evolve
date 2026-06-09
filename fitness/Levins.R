library(tidyverse)
library(MASS)
library(ggplot2)

centroids <- cloud_plot_data %>%
  group_by(PopSpecies) %>%
  summarise(
    Mean_B = mean(Fitness_B),
    Mean_T = mean(Fitness_T)
  )

B_spec_point <- centroids %>% filter(PopSpecies == "Barbarea")
T_spec_point <- centroids %>% filter(PopSpecies == "Turritus")
M_gen_point  <- centroids %>% filter(PopSpecies == "Mixture")

W_B_max <- B_spec_point$Mean_B
W_T_max <- T_spec_point$Mean_T
W_Mix_B <- M_gen_point$Mean_B
W_Mix_T <- M_gen_point$Mean_T

levins_eq <- function(alpha) {
  (W_Mix_B / W_B_max)^alpha + (W_Mix_T / W_T_max)^alpha - 1
}

alpha_sol <- tryCatch({
  uniroot(levins_eq, c(0.01, 20))$root
}, error = function(e) NA)

curve_data <- tibble(x = seq(0, W_B_max, length.out = 200)) %>%
  mutate(
    y = W_T_max * (1 - (x / W_B_max)^alpha_sol)^(1/alpha_sol)
  )

vectors <- bind_rows(B_spec_point, T_spec_point)

plot_limit <- max(c(W_B_max, W_T_max)) * 1.05

g_levins_fit <- ggplot(cloud_plot_data, aes(x = Fitness_B, y = Fitness_T)) +
  
  
  geom_segment(data = vectors,
               aes(x = 1.5, y = 1.5, xend = Mean_B, yend = Mean_T),
               color = "gray70", linetype = "dashed", size = 0.6) +
  
  geom_point(aes(color = PopSpecies), alpha = 0.05, size = 1, shape = 16, stroke = 0) +
  geom_density_2d(aes(color = PopSpecies), size = 0.3, alpha = 0.6, bins = 3) +
  
  geom_line(data = curve_data, aes(x = x, y = y),
            color = "black", size = 1.2) +
  
  geom_point(data = centroids, aes(x = Mean_B, y = Mean_T, fill = PopSpecies), 
             size = 4.5, shape = 21, color = "black", stroke = 1.5) +
  
  annotate("text", x = W_B_max * 0.5, y = W_T_max * 0.5, 
           label = paste0("Levins Fit: α = ", round(alpha_sol, 2)), 
           size = 4, fontface = "italic", color = "gray20") +
  
  scale_color_manual(values = c("Barbarea" = "#499FFF", "Mixture" = "#B01754", "Turritus" = "#EDB72D")) +
  scale_fill_manual(values = c("Barbarea" = "#499FFF", "Mixture" = "#B01754", "Turritus" = "#EDB72D")) +
  
  labs(
    x = "Fitness on Barbarea (Specialist Host B)",
    y = "Fitness on Turritus (Specialist Host T)",
    title = "Levins' Fitness Set",
    subtitle = "Curve represents the theoretical trade-off boundary (Power Law Fit)",
    color = "Evolution Treatment",
    fill = "Evolution Treatment"
  ) +
  
  coord_fixed(xlim = c(1.8, plot_limit), 
              ylim = c(1.8, plot_limit)) + 
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold"),
    legend.position = "right"
  )

print(g_levins_fit)
ggsave("fig_levins_exact_fit.png", g_levins_fit, width = 7, height = 6, dpi = 300)