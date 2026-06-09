library(tidyr)
library(dplyr)
library(ggplot2)
library(MASS)
library(lme4)
library(multcomp)
library(statmod)
library(MASS)
library(mgcv)
library(merTools)
library(ggpubr)
library(ggrepel)
library(emmeans)
library(ggh4x)
library(ggtext)
library(ggsignif)

source("scripts/data_setup.R")


source("scripts/colors.R")
source("scripts/viz_data.R")
gbar
gcum
gtime + theme(panel.grid = element_blank())


source("scripts/contrasts.R")


source("scripts/helper.R")

glm <- glm.nb(numFlies ~
                PopSpecies + 
                PhenotypeCageSpecies + 
                PopSpecies:PhenotypeCageSpecies + 
                offset(log(trays)),
              contrasts = list(PopSpecies = contrast_matrix, PhenotypeCageSpecies = contrast_matrix), data = tab)


summary(glm)
emm <- emmeans(glm, ~ PopSpecies * PhenotypeCageSpecies, type = "response")
pairs(emm, by = "PhenotypeCageSpecies", adjust = "none")
test_for_overdispersion(glm)


gdat <- edit_plot_data(glm_sim_nb_slim(glm, n_sims = 100, dframe = full))
gdat$pred <- predict(glm, newdata = full, type = "response")
pred_se <- predict(glm, newdata = full, type = "response", se.fit = TRUE)
pred_se$se.fit
gdat$pred_upper <- gdat$pred + (1.96*pred_se$se.fit)
gdat$pred_lower <- gdat$pred - (1.96*pred_se$se.fit)

tab2 <- tab
tab2$PopSpecies <- as.character(tab2$PopSpecies)
tab2$PhenotypeCageSpecies <- as.character(tab2$PhenotypeCageSpecies)
tab2[tab2$PopSpecies %in% "M", "PopSpecies"] <- rep("Mixture", length(tab2[tab2$PopSpecies %in% "M", "PopSpecies"]))
tab2[tab2$PopSpecies %in% "T", "PopSpecies"] <- rep("Turritus", length(tab2[tab2$PopSpecies %in% "T", "PopSpecies"]))
tab2[tab2$PopSpecies %in% "B", "PopSpecies"] <- rep("Barbarea", length(tab2[tab2$PopSpecies %in% "B", "PopSpecies"]))
tab2[tab2$PhenotypeCageSpecies %in% "M", "PhenotypeCageSpecies"] <- rep("Tested\n on Mixture", length(tab2[tab2$PhenotypeCageSpecies %in% "M", "PhenotypeCageSpecies"]))
tab2[tab2$PhenotypeCageSpecies %in% "T", "PhenotypeCageSpecies"] <- rep("Tested\n on Turritus", length(tab2[tab2$PhenotypeCageSpecies %in% "T", "PhenotypeCageSpecies"]))
tab2[tab2$PhenotypeCageSpecies %in% "B", "PhenotypeCageSpecies"] <- rep("Tested\n on Barbarea", length(tab2[tab2$PhenotypeCageSpecies %in% "B", "PhenotypeCageSpecies"]))


ggplot(gdat, aes(x = PhenotypeCageSpecies, y = pred, col = PopSpecies)) +
  geom_point(tab2, mapping = aes(x = PhenotypeCageSpecies, y = normFlies,  col = PopSpecies), shape = 16, size = 1, alpha = 0.6, inherit.aes = F, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1)) + 
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_line(gdat[gdat$PopSpecies %in% "Mixture", ], mapping = aes(x = PhenotypeCageSpecies, y = pred, col = PopSpecies, group = PopSpecies), linetype = "dashed", inherit.aes = F, position = position_dodge(0.5), size = 0.5) + 
  geom_line(gdat[!gdat$PhenotypeCageSpecies %in% "Tested\n on Mixture" & !gdat$PopSpecies %in% "Mixture", ], mapping = aes(x = PhenotypeCageSpecies, y = pred, col = PopSpecies, group = PopSpecies), linetype = "dashed", inherit.aes = F, position = position_dodge(0.5), size = 0.5) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                position = position_dodge(0.5), width = 0.2, size = 0.6) +
  labs(x = "", y = "F1 adults after viability selection\n (per female per seedling)", col = "Source Treatment") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(color = c("#499FFF", "#B01754", "#EDB72D"), size = 10, face = "bold"), 
        legend.position = "none") + 
  scale_color_manual(values = custom_colors2)



y = 6.3
yend1 = 6.3
yend2 = 6
yann = 6.4
space = 0.6
space2 = 1.1
space3 = 

gdat$highlight <- c()
gdat[
  (gdat$PopSpecies %in% "Barbarea" & gdat$PhenotypeCageSpecies %in% "Tested\n on Barbarea") | 
  (gdat$PopSpecies %in% "Turritus" & gdat$PhenotypeCageSpecies %in% "Tested\n on Turritus") |
  (gdat$PopSpecies %in% "Mixture" & gdat$PhenotypeCageSpecies %in% "Tested\n on Mixture"), "highlight"] <- "Home Treatment"
gdat[is.na(gdat$highlight), "highlight"] <- "Away Treatment"

g <- ggplot(gdat, aes(x = PhenotypeCageSpecies, y = pred, col = PopSpecies, shape = highlight, size = highlight, group = PopSpecies)) +
  geom_point(tab2, mapping = aes(x = PhenotypeCageSpecies, y = normFlies,  col = PopSpecies), alpha = 0.2, shape = 16, size = 1, inherit.aes = F, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.05)) + 
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pred_upper, ymax = pred_lower),
                position = position_dodge(0.5), width = 0.2, size = 0.5) +
  geom_line(linetype = "dashed", position = position_dodge(0.5), size = 0.5) + 
  annotate("text", x = 2, y = 8, label = "Local Adaptation (p = 0.004)", size = 3, color = "black", fontface = "bold") +  # Top left
  annotate("text", x = 1.3, y = 0.3, label = "Error bars: 95% CI (1.96*SE)", size = 2.5, color = "black", fontface = "italic") +  # Bottom left
  labs(x = "", y = "F1 adults after viability selection\n (per female per seedling)", col = "Evolution Treatment", shape = "Test Context") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(color = c("#499FFF", "#B01754", "#EDB72D"), size = 8, face = "bold"), 
        legend.position = "right") + 
  scale_color_manual(values = custom_colors2) +  
  scale_size_manual(values = c("Home Treatment" = 4, "Away Treatment" = 2)) +
  scale_shape_manual(values = c("Home Treatment" = 16, "Away Treatment" = 17)) +
  guides(
    size = "none",
    shape = guide_legend(override.aes = list(size = 4)),
    col = guide_legend(override.aes = list(size = 4, linetype = 0))
  )


ggsave("fig2b.png", g, scale = 0.65, units = "in", width = 7.5, height = 5)




tab$Matching <- ifelse(tab$PopSpecies == tab$PhenotypeCageSpecies, "Match", "Mismatch")
glm_match <- glm.nb(numFlies ~ PopSpecies + Matching + offset(log(trays)), data = tab)
summary(glm_match)


weights_vec <- c(-0.5, 1, -0.5) 

step1_emm <- contrast(emm, 
                      method = list("M_vs_Others" = weights_vec), 
                      by = "PhenotypeCageSpecies")

final_interaction_test <- contrast(step1_emm, 
                                   method = list("Home_Advantage_Interaction" = weights_vec),
                                   by = NULL)

summary(final_interaction_test, side = ">")
