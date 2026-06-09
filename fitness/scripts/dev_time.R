gtdat$date3 <- as.numeric(sub(".Feb", "", gtdat$date2))

gtdat %>% 
  group_by(PhenotypeCageSpecies) %>% 
  summarise(val = mean( normFlies ))
  
  
  summarise(p_value = kruskal.test(val ~ date)$p.value)


model <- glm(date3 ~ PhenotypeCageSpecies, data = gtdat, weights = numFlies)
anova(model, test = "F")
summary(model)


timing_data <- gtdat %>%
  filter(!is.na(val) & val > 0) %>%
  mutate(date_numeric = as.numeric(sub("Feb", "", date2))) %>%
  group_by(Rep, PopSpecies, PhenotypeCageSpecies) %>%
  summarise(
    weighted_mean_day = weighted.mean(date_numeric, val, na.rm = TRUE),
    total_flies = sum(val, na.rm = TRUE),
    .groups = 'drop'
  )

timing_data %>%
  group_by(PhenotypeCageSpecies, PopSpecies) %>%
  summarise(
    mean_emergence = mean(weighted_mean_day),
    se_emergence = sd(weighted_mean_day)/sqrt(n()),
    .groups = 'drop'
  )

###### Plot ######

timing_data <- gtdat %>%
  filter(!is.na(val) & val > 0) %>%
  mutate(date_numeric = as.numeric(sub("Feb", "", date2))) %>%
  group_by(PopNumber, PopSpecies, PhenotypeCageSpecies) %>%
  summarise(
    weighted_mean_day = weighted.mean(date_numeric, val, na.rm = TRUE),
    total_flies = sum(val, na.rm = TRUE),
    .groups = 'drop'
  )


summary_data <- timing_data %>%
  group_by(PhenotypeCageSpecies, PopSpecies) %>%
  summarise(
    mean_day = mean(weighted_mean_day, na.rm = TRUE),
    se_day = sd(weighted_mean_day, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

summary_data$PhenotypeCageSpecies <- gsub("Tested on ", "Tested\n on ", summary_data$PhenotypeCageSpecies)

summary_data[
  (summary_data$PopSpecies %in% "Barbarea" & summary_data$PhenotypeCageSpecies %in% "Tested\n on Barbarea") | 
    (summary_data$PopSpecies %in% "Turritus" & summary_data$PhenotypeCageSpecies %in% "Tested\n on Turritus") |
    (summary_data$PopSpecies %in% "Mixture" & summary_data$PhenotypeCageSpecies %in% "Tested\n on Mixture"), "highlight"] <- "Home Treatment"
summary_data[is.na(summary_data$highlight), "highlight"] <- "Away Treatment"

g <- ggplot(summary_data, aes(x = PhenotypeCageSpecies, y = mean_day, col = PopSpecies, shape = highlight, size = highlight)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_day - se_day, ymax = mean_day + se_day),
                position = position_dodge(width = 0.8), 
                width = 0.25, size = 0.8) +
  scale_color_manual(values = custom_colors2) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                     labels = function(x) paste0("Feb\n ", x)) +
  scale_size_manual(values = c("Home Treatment" = 4, "Away Treatment" = 2)) +
  scale_shape_manual(values = c("Home Treatment" = 16, "Away Treatment" = 17)) +
  labs(y = "Weighted Mean\n Emergence Date ",
       x = "",
       col = "Evolution Treatment") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(color = c("#499FFF", "#B01754", "#EDB72D"), size = 6, face = "bold"), 
        legend.position = "none") + 
  guides(
    size = "none",
    shape = guide_legend(override.aes = list(size = 4)),
    col = guide_legend(override.aes = list(size = 4, linetype = 0))
  )



####### Stats ######

# Weighted linear model (no random effect, since singular fit)
m_fixed <- lm(weighted_mean_day ~ PhenotypeCageSpecies * PopSpecies,
              data = timing_data,
              weights = total_flies)

# ANOVA table
anova(m_fixed)
summary(m_fixed)

library(car)
anova_tab <- Anova(m_fixed, type = 2)

p_x   <- signif(anova_tab["PhenotypeCageSpecies", "Pr(>F)"], 3)
p_col <- signif(anova_tab["PopSpecies", "Pr(>F)"], 3)
p_int <- signif(anova_tab["PhenotypeCageSpecies:PopSpecies", "Pr(>F)"], 3)
coefs <- summary(m_fixed)$coefficients
coefs[grep(":", rownames(coefs)), ] 


annot_text <- paste0(
  "Test Host Plant: p = ", p_x, " \n",
  "Evolution Host Plant: p = ", p_col, "\n",
  "Interaction: p = ", p_int, "\n"
)

g2 <- g +
  annotate("text",
           x = 1.85, y = min(summary_data$mean_day) + 0.4,
           label = annot_text,
           hjust = 0, vjust = 1,
           size = 2.8,
           fontface = "italic")

ggsave("fig2C.png", g2, scale = 0.8, units = "in", width = 5, height = 4)



