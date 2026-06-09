library(tidyverse)
library(MASS)
library(lme4)
library(multcomp)
library(statmod)
library(mgcv)
library(merTools)
library(ggpubr)
library(ggrepel)
library(emmeans)
library(ggh4x)
library(ggtext)
library(ggsignif)   
library(svglite)

source("scripts/data_setup.R")
source("scripts/colors.R")
source("scripts/viz_data.R")
source("scripts/contrasts.R")
source("scripts/helper.R")

glm_model <- glm.nb(numFlies ~ 
                      PopSpecies + 
                      PhenotypeCageSpecies + 
                      PopSpecies:PhenotypeCageSpecies + 
                      offset(log(trays)),
                    contrasts = list(PopSpecies = contrast_matrix, 
                                     PhenotypeCageSpecies = contrast_matrix), 
                    data = tab)

gdat <- edit_plot_data(glm_sim_nb_slim(glm_model, n_sims = 100, dframe = full))
pred_se <- predict(glm_model, newdata = full, type = "response", se.fit = TRUE)
gdat$pred <- predict(glm_model, newdata = full, type = "response")
gdat$pred_upper <- gdat$pred + (1.96 * pred_se$se.fit)
gdat$pred_lower <- gdat$pred - (1.96 * pred_se$se.fit)

gdat <- gdat %>%
  mutate(highlight = case_when(
    PopSpecies == "Barbarea" & PhenotypeCageSpecies == "Tested\n on Barbarea" ~ "Home Treatment",
    PopSpecies == "Turritus" & PhenotypeCageSpecies == "Tested\n on Turritus" ~ "Home Treatment",
    PopSpecies == "Mixture"  & PhenotypeCageSpecies == "Tested\n on Mixture"  ~ "Home Treatment",
    TRUE ~ "Away Treatment"
  ))

tab2 <- tab %>%
  mutate(
    PopSpecies = recode(PopSpecies, "M"="Mixture", "T"="Turritus", "B"="Barbarea"),
    PhenotypeCageSpecies = recode(PhenotypeCageSpecies, "M"="Tested\n on Mixture", "T"="Tested\n on Turritus", "B"="Tested\n on Barbarea")
  )

x_levels <- c("Tested\n on Barbarea", "Tested\n on Mixture", "Tested\n on Turritus")
pop_levels <- c("Barbarea", "Mixture", "Turritus")

gdat$PhenotypeCageSpecies <- factor(gdat$PhenotypeCageSpecies, levels = x_levels)
gdat$PopSpecies           <- factor(gdat$PopSpecies, levels = pop_levels)
tab2$PhenotypeCageSpecies <- factor(tab2$PhenotypeCageSpecies, levels = x_levels)
tab2$PopSpecies           <- factor(tab2$PopSpecies, levels = pop_levels)

avg_fitness <- gdat %>%
  group_by(PopSpecies) %>%
  summarise(avg_y = mean(pred), .groups = "drop")


g_final <- ggplot(gdat, aes(x = PhenotypeCageSpecies, y = pred, col = PopSpecies, 
                            shape = highlight, size = highlight, group = PopSpecies)) +
  
  geom_point(data = tab2, 
             mapping = aes(x = PhenotypeCageSpecies, y = normFlies, col = PopSpecies, group = PopSpecies), 
             alpha = 0.2, shape = 16, size = 1, inherit.aes = FALSE, 
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.05)) + 
  
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pred_upper, ymax = pred_lower),
                position = position_dodge(0.5), width = 0.2, size = 0.5) +
  geom_line(linetype = "dashed", position = position_dodge(0.5), size = 0.5) + 
  
  geom_signif(
    y_position = 4.5, xmin = 1.9, xmax = 2.1, 
    annotation = "NS", tip_length = 0.02, 
    col = "black", size = 0.6, textsize = 3, fontface = "bold"
  ) +
  
  annotate("text", x = 2, y = 8, label = "Local Adaptation (p = 0.004)", 
           size = 3, color = "black", fontface = "bold") +
  annotate("text", x = 1.3, y = 0.3, label = "Error bars: 95% CI (1.96*SE)", 
           size = 2.5, color = "black", fontface = "italic") +
  
  geom_segment(data = avg_fitness,
               aes(x = 3.55, xend = 3.75, y = avg_y, yend = avg_y, col = PopSpecies),
               size = 2, inherit.aes = FALSE) +
  
  annotate("text", x = 3.65, y = 8, label = "Global\nAvg.", 
           size = 2.8, fontface = "italic", lineheight = 0.9) +
  
  labs(x = "", y = "F1 adults after viability selection\n (per female per seedling)", 
       col = "Evolution Treatment", shape = "Test Context") + 
  
  scale_color_manual(values = custom_colors2) +  
  scale_size_manual(values = c("Home Treatment" = 4, "Away Treatment" = 2)) +
  scale_shape_manual(values = c("Home Treatment" = 16, "Away Treatment" = 17)) +
  
  guides(
    size = "none",
    shape = guide_legend(override.aes = list(size = 4)),
    col = guide_legend(override.aes = list(size = 4, linetype = 0))
  ) + 
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(color = c("#499FFF", "#B01754", "#EDB72D"), size = 8, face = "bold"), 
    legend.position = "right",
    
    coord = coord_cartesian(clip = "off"), 
    plot.margin = margin(10, 50, 10, 10)
  )

ggsave("fig2b_annotated.png", g_final, scale = 0.65, units = "in", width = 8.5, height = 5, dpi = 300)

ggsave("fig2b_annotated.svg", g_final, scale = 0.65, units = "in", width = 8.5, height = 5)


marginal_emm <- emmeans(glm_model, ~ PopSpecies, type = "response")
bar_data <- as.data.frame(marginal_emm)

library(dplyr)
bar_data <- bar_data %>%
  mutate(PopSpecies = recode(PopSpecies, 
                             "M" = "Mixture", 
                             "T" = "Turritus", 
                             "B" = "Barbarea"))

pop_levels <- c("Barbarea", "Mixture", "Turritus") 
bar_data$PopSpecies <- factor(bar_data$PopSpecies, levels = pop_levels)

g_bar <- ggplot(bar_data, aes(x = PopSpecies, y = response, fill = PopSpecies)) +
  
  geom_col(width = 0.7, alpha = 0.9, color = "black", size = 0.3) +
  
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, size = 0.6) +
  
  labs(x = "Evolution Treatment", 
       y = "Global Average Fitness\n(Marginal Mean F1 Adults)",
       title = "Overall Performance") +
  
  scale_fill_manual(values = custom_colors2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(face = "bold", size = 11, color = c("#499FFF", "#B01754", "#EDB72D")),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  )

ggsave("fig2b_marginal_barplot.png", g_bar, width = 4, height = 5, dpi = 300)

ggsave("fig2b_marginal_barplot.svg", g_bar, width = 4, height = 5)

empirical_stats <- tab2 %>%
  group_by(PopSpecies, PhenotypeCageSpecies) %>%
  summarise(
    N = n(),
    Emp_Mean = mean(normFlies, na.rm = TRUE),
    Emp_SD   = sd(normFlies, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    Emp_SE = Emp_SD / sqrt(N),
    Emp_Lower = Emp_Mean - (qt(0.975, N - 1) * Emp_SE),
    Emp_Upper = Emp_Mean + (qt(0.975, N - 1) * Emp_SE)
  )

emm_obj <- emmeans(glm_model, ~ PopSpecies + PhenotypeCageSpecies, 
                   type = "response", 
                   offset = 0)

model_stats <- as.data.frame(emm_obj) %>%
  dplyr::rename(
    Mod_Pred  = response,
    Mod_Lower = asymp.LCL,
    Mod_Upper = asymp.UCL
  ) %>%
  dplyr::select(PopSpecies, PhenotypeCageSpecies, Mod_Pred, Mod_Lower, Mod_Upper)

model_stats <- model_stats %>%
  mutate(
    PopSpecies = recode(PopSpecies, 
                        "M"="Mixture", 
                        "T"="Turritus", 
                        "B"="Barbarea"),
    PhenotypeCageSpecies = recode(PhenotypeCageSpecies, 
                                  "M"="Tested\n on Mixture", 
                                  "T"="Tested\n on Turritus", 
                                  "B"="Tested\n on Barbarea")
  )

summary_table <- left_join(empirical_stats, model_stats, 
                           by = c("PopSpecies", "PhenotypeCageSpecies")) %>%
  dplyr::select(
    Evolution_Treatment = PopSpecies,
    Test_Environment    = PhenotypeCageSpecies,
    
    Empirical_Mean      = Emp_Mean,
    Empirical_95_Lower  = Emp_Lower,
    Empirical_95_Upper  = Emp_Upper,
    
    Model_Prediction    = Mod_Pred,
    Model_95_Lower      = Mod_Lower,
    Model_95_Upper      = Mod_Upper
  ) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3)))


write.csv(summary_table, "evolution_treatment_summary_table.csv", row.names = FALSE)