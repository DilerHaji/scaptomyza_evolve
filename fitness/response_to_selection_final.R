library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)

rs <- read.csv("SelectionExperiment_Gen23andFinal_numFlies.csv")
rs$comb <- paste(rs$PlantSpecies, rs$Rep, sep = "")

rs[rs$Generation %in% "final", "numFlies"] <- rs[rs$Generation %in% "final", "numFlies"] + 200

################ Raw counts for all treatments ####################

rss <- rs[rs$Generation %in% c("2", "3"),]
rs2 <- setNames(
  aggregate(rss[, "numFlies"], by = list(rss[, "PlantSpecies"], rss[, "Generation"], rss[, "comb"]), FUN = sum),
  c("PlantSpecies", "Generation", "comb", "numFlies"))

rfinal <- rs[!rs$Generation %in% c("2", "3"),]
rfinal <- rfinal[, c("PlantSpecies", "Generation", "comb", "numFlies")]

rs3 <- rbind(rs2, rfinal)


################ Raw counts relative to mixture ####################

#rsr <- rs3
rsr <- rs

rsr$RelativeFlies <- c()

rsrmean <- setNames(
  aggregate(rsr$numFlies, by = list(rsr$PlantSpecies, rsr$Generation, rsr$comb), FUN = mean),
  c("PlantSpecies", "Generation", "comb", "numFlies")) 
rsrs <- spread(rsrmean, key = PlantSpecies, value = numFlies)

rsr[rsr$Generation %in% "2", "RelativeFlies"] <- rsr[rsr$Generation %in% "2", "numFlies"]/mean(rsrs[rsrs$Generation %in% "2", "M"], na.rm = TRUE)
rsr[rsr$Generation %in% "3", "RelativeFlies"] <- rsr[rsr$Generation %in% "3", "numFlies"]/mean(rsrs[rsrs$Generation %in% "3", "M"], na.rm = TRUE)
rsr[rsr$Generation %in% "final", "RelativeFlies"] <- rsr[rsr$Generation %in% "final", "numFlies"]/mean(rsrs[rsrs$Generation %in% "final", "M"], na.rm = TRUE)


rsr[rsr$PlantSpecies == "B", "PlantSpecies"] <- "Barbarea flies"
rsr[rsr$PlantSpecies == "T", "PlantSpecies"] <- "Turritus flies"
rsr[rsr$PlantSpecies == "M", "PlantSpecies"] <- "Mixture flies"
rsr[rsr$Generation == 2, "Generation"] <- 2
rsr[rsr$Generation == 3, "Generation"] <- 3
rsr[rsr$Generation == "final", "Generation"] <- 10

rsr$Generation <- as.numeric(rsr$Generation)

rsrmean2 <- setNames(
  aggregate(rsr$RelativeFlies, by = list(rsr$PlantSpecies, rsr$Generation), FUN = mean),
  c("PlantSpecies", "Generation", "RelativeFlies")) 


custom_colors <- c("Turritus flies" = "#EDB72D", "Barbarea flies" = "#499FFF", "Mixture flies" = "#B01754")
d = 1
ggplot(rsr, aes(x = Generation, y = RelativeFlies, col = PlantSpecies)) + 
  geom_point(position = position_jitterdodge(dodge.width = d, jitter.width = 0.01), shape = 16, size = 1) + 
  geom_point(rsrmean2, mapping = aes(x = Generation, y = RelativeFlies, col = PlantSpecies, group = PlantSpecies), inherit.aes = F, size = 3) + 
  geom_line(rsrmean2, mapping = aes(x = Generation, y = RelativeFlies, col = PlantSpecies, group = PlantSpecies), inherit.aes = F, linetype = "dashed") + 
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  labs(x = "Generation", y = "Population viability\n relative to the average\n mixed-host environment", col = "") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(), 
        legend.position = "none")


###### Updated version ######

rsr4 <- rsr[rsr$PlantSpecies %in% c("Barbarea flies", "Turritus flies"), ]
rsr4[rsr4$Generation %in% c(2,3), "Generation"] <- "2/3"
rsr4_se <- aggregate(rsr4$RelativeFlies, by = list(rsr4$PlantSpecies, rsr4$Generation), FUN = function(x){   sd(x[!is.na(x)]) / sqrt(length(x[!is.na(x)])) } )
rsr4_m <- aggregate(rsr4$RelativeFlies, by = list(rsr4$PlantSpecies,  rsr4$Generation), FUN = function(x){mean(x[!is.na(x)]) } )

rsr5 <- data.frame(PlantSpecies = rsr4_se$Group.1, Generation = rsr4_se$Group.2, se = rsr4_se$x, m = rsr4_m$x)

# d = 1
d = 0.2
g <- ggplot(rsr5, aes(x = factor(Generation, levels = c("2/3", "10")), y = m, col = PlantSpecies, group = PlantSpecies)) + 
  geom_point(position = position_jitterdodge(dodge.width = d, jitter.width = 0), shape = 16, size = 5) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_jitterdodge(dodge.width = d, jitter.width = 0), width = 0.1) +  
  geom_line(linetype = "dashed", position = position_jitterdodge(dodge.width = d, jitter.width = 0)) + 
  scale_color_manual(values = custom_colors) +
  scale_x_discrete(expand = c(0.2, 0.01, 0.2, 0)) +  # Reduce space on sides
  labs(x = "Generation", y = "Population viability\n relative to mean\n mixed-host plant environment", col = "") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(), 
        legend.position = "none")

ggsave("fig1a.png", g, scale = 1, width = 3, height = 2.5)



########## Statistical test ############

df <- rsr4

t.test(rsr4[rsr4$Generation == "2/3" &  rsr4$PlantSpecies %in% "Barbarea flies", "RelativeFlies"],  
       rsr4[rsr4$Generation == "10" & rsr4$PlantSpecies %in% "Barbarea flies", "RelativeFlies"])

t.test(rsr4[rsr4$Generation == "2/3" &  rsr4$PlantSpecies %in% "Turritus flies", "RelativeFlies"],  
       rsr4[rsr4$Generation == "10" & rsr4$PlantSpecies %in% "Turritus flies", "RelativeFlies"])

t.test(sample(rsr4[rsr4$Generation == "2/3", "RelativeFlies"], 8),  
       rsr4[rsr4$Generation == "10", "RelativeFlies"], paired = T)



barbarea_23 <- df[df$PlantSpecies == "Barbarea flies" & df$Generation == "2/3", ]
barbarea_10 <- df[df$PlantSpecies == "Barbarea flies" & df$Generation == "10", ]

turritus_23 <- df[df$PlantSpecies == "Turritus flies" & df$Generation == "2/3", ]
turritus_10 <- df[df$PlantSpecies == "Turritus flies" & df$Generation == "10", ]

barbarea_23 <- barbarea_23[order(barbarea_23$comb), ]
barbarea_10 <- barbarea_10[order(barbarea_10$comb), ]

turritus_23 <- turritus_23[order(turritus_23$comb), ]
turritus_10 <- turritus_10[order(turritus_10$comb), ]

barbarea_23_means <- aggregate(RelativeFlies ~ comb, barbarea_23, mean)
barbarea_10_means <- aggregate(RelativeFlies ~ comb, barbarea_10, mean)

turritus_23_means <- aggregate(RelativeFlies ~ comb, turritus_23, mean)
turritus_10_means <- aggregate(RelativeFlies ~ comb, turritus_10, mean)

barbarea_test <- t.test(barbarea_23_means$RelativeFlies, 
                        barbarea_10_means$RelativeFlies, 
                        paired = TRUE)

turritus_test <- t.test(turritus_23_means$RelativeFlies, 
                        turritus_10_means$RelativeFlies, 
                        paired = TRUE)

####### Sign test #######

barbarea_23 <- df[df$PlantSpecies == "Barbarea flies" & df$Generation == "2/3", ]
barbarea_10 <- df[df$PlantSpecies == "Barbarea flies" & df$Generation == "10", ]

turritus_23 <- df[df$PlantSpecies == "Turritus flies" & df$Generation == "2/3", ]
turritus_10 <- df[df$PlantSpecies == "Turritus flies" & df$Generation == "10", ]

barbarea_23 <- barbarea_23[order(barbarea_23$comb), ]
barbarea_10 <- barbarea_10[order(barbarea_10$comb), ]

turritus_23 <- turritus_23[order(turritus_23$comb), ]
turritus_10 <- turritus_10[order(turritus_10$comb), ]

barbarea_23_means <- aggregate(RelativeFlies ~ comb, barbarea_23, mean)
barbarea_10_means <- aggregate(RelativeFlies ~ comb, barbarea_10, mean)

turritus_23_means <- aggregate(RelativeFlies ~ comb, turritus_23, mean)
turritus_10_means <- aggregate(RelativeFlies ~ comb, turritus_10, mean)

# Barbarea
barbarea_diff <- barbarea_10_means$RelativeFlies - barbarea_23_means$RelativeFlies
barbarea_positive <- sum(barbarea_diff > 0)
barbarea_total <- length(barbarea_diff[barbarea_diff != 0])  

barbarea_sign_test <- binom.test(barbarea_positive, barbarea_total, p = 0.5)

# Turritus  
turritus_diff <- turritus_10_means$RelativeFlies - turritus_23_means$RelativeFlies
turritus_positive <- sum(turritus_diff > 0)
turritus_total <- length(turritus_diff[turritus_diff != 0]) 

turritus_sign_test <- binom.test(turritus_positive, turritus_total, p = 0.5)



########## Chi-Square test using relative prop ####### 

rs <- read.csv("SelectionExperiment_Gen23andFinal_numFlies.csv")
library(dplyr)

gen <- c(2,3)

data_summary <- rs %>%
  group_by(Generation, PlantSpecies, Rep) %>%
  summarise(mean_numFlies = mean(numFlies), .groups = "keep")

mixture_averages <- data_summary %>%
  filter(PlantSpecies == "M") %>%
  group_by(Generation) %>%
  summarise(mixture_baseline = mean(mean_numFlies), .groups = "drop")


calculate_relative_to_mixture <- function() {
  early_data <- data_summary %>%
    filter(Generation %in% gen, PlantSpecies %in% c("B", "T")) %>%
    group_by(PlantSpecies, Rep) %>%
    summarise(early_flies = mean(mean_numFlies), .groups = "drop")
  
  final_data <- data_summary %>%
    filter(Generation == "final", PlantSpecies %in% c("B", "T")) %>%
    select(PlantSpecies, Rep, final_flies = mean_numFlies)
  
  early_mixture <- mean(mixture_averages$mixture_baseline[mixture_averages$Generation %in% c(2, 3)])
  final_mixture <- mixture_averages$mixture_baseline[mixture_averages$Generation == "final"]
  
  early_data$early_relative <- early_data$early_flies / early_mixture
  final_data$final_relative <- final_data$final_flies / final_mixture
  
  relative_data <- merge(early_data, final_data, by = c("PlantSpecies", "Rep"))
  
  return(list(data = relative_data, early_mix = early_mixture, final_mix = final_mixture))
}

relative_results <- calculate_relative_to_mixture()
relative_data <- relative_results$data

relative_data$early_above_M <- relative_data$early_relative > 1.0
relative_data$final_above_M <- relative_data$final_relative > 1.0

early_above_count <- sum(relative_data$early_above_M)
early_below_count <- sum(!relative_data$early_above_M)
final_above_count <- sum(relative_data$final_above_M)
final_below_count <- sum(!relative_data$final_above_M)

time_contingency <- matrix(c(early_above_count, early_below_count,
                             final_above_count, final_below_count),
                           nrow = 2, byrow = TRUE,
                           dimnames = list(c("Early", "Final"), 
                                           c("Above_M", "Below_M")))

chi_test <- chisq.test(time_contingency)


early_prop_above <- early_above_count / (early_above_count + early_below_count)
final_prop_above <- final_above_count / (final_above_count + final_below_count)


fisher_test <- fisher.test(time_contingency)



########## McNemar's test ####### 

mcnemar_table <- matrix(c(
  sum(relative_data$early_above_M & relative_data$final_above_M),   # Above→Above
  sum(relative_data$early_above_M & !relative_data$final_above_M),  # Above→Below
  sum(!relative_data$early_above_M & relative_data$final_above_M),  # Below→Above  
  sum(!relative_data$early_above_M & !relative_data$final_above_M)  # Below→Below
), nrow = 2, byrow = TRUE,
dimnames = list(c("Early_Above", "Early_Below"), c("Final_Above", "Final_Below")))

discordant_pairs <- mcnemar_table[1,2] + mcnemar_table[2,1]

if (discordant_pairs > 0) {
  mcnemar_result <- mcnemar.test(mcnemar_table)
  cat("\nMcNemar's Test Results:\n")
  cat("McNemar's chi-square:", round(mcnemar_result$statistic, 4), "\n")
  cat("P-value:", round(mcnemar_result$p.value, 4), "\n")
  cat("Significant (p < 0.05):", ifelse(mcnemar_result$p.value < 0.05, "YES", "NO"), "\n")
  
  # Direction of change
  improved <- mcnemar_table[2,1]  # Below→Above
  declined <- mcnemar_table[1,2]  # Above→Below
  
  cat("\nDirection of change:\n")
  cat("Replicates that improved (Below M → Above M):", improved, "\n")
  cat("Replicates that declined (Above M → Below M):", declined, "\n")
  
  if (improved > declined) {
    cat("Overall direction: Improvement relative to M baseline\n")
  } else if (declined > improved) {
    cat("Overall direction: Decline relative to M baseline\n")
  } else {
    cat("Overall direction: No net change\n")
  }
} else {
  cat("No discordant pairs found - no individual replicates changed status\n")
}


for (species in c("B", "T")) {
  species_name <- ifelse(species == "B", "Barbarea", "Turritus")
  species_data <- relative_data[relative_data$PlantSpecies == species, ]
  
  species_mcnemar <- matrix(c(
    sum(species_data$early_above_M & species_data$final_above_M),
    sum(species_data$early_above_M & !species_data$final_above_M),
    sum(!species_data$early_above_M & species_data$final_above_M),
    sum(!species_data$early_above_M & !species_data$final_above_M)
  ), nrow = 2, byrow = TRUE)
  
  cat("\n", species_name, " McNemar's table:\n")
  print(species_mcnemar)
  
  species_discordant <- species_mcnemar[1,2] + species_mcnemar[2,1]
  if (species_discordant > 0) {
    species_test <- mcnemar.test(species_mcnemar)
    cat("P-value:", round(species_test$p.value, 4), "\n")
    cat("Significant:", ifelse(species_test$p.value < 0.05, "YES", "NO"), "\n")
  } else {
    cat("No changes in this species\n")
  }
}




########## Hemcnemar_data <- contingency_data %>%
  mutate(
    early_dominant = ifelse(early_flies > median(c(early_flies)), "High", "Low"),
    final_dominant = ifelse(final_flies > median(c(final_flies)), "High", "Low")
  )ritability ############

df <- rs3

df <- df %>%
  mutate(gen_group = case_when(
    Generation %in% c("2", "3") ~ "G0",
    Generation == "final" ~ "G10",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(gen_group))

colnames(df) <- c("treatment", "Generation", "comb", "viability", "gen_group")

df <- df %>%
  group_by(gen_group) %>%
  mutate(log_ratio = log(viability / exp(mean(log(viability[treatment == "M"]))))) %>%
  ungroup()

avg_logratio <- df %>%
  group_by(gen_group, treatment) %>%
  summarise(mean_log_ratio = mean(log_ratio), .groups = "drop")

mean_G0_all <- mean(avg_logratio$mean_log_ratio[avg_logratio$gen_group == "G0"])

S <- avg_logratio %>%
  filter(gen_group == "G0" & treatment %in% c("B", "T")) %>%
  mutate(S = mean_log_ratio - mean_G0_all) %>%
  dplyr::select(treatment, S)

R <- avg_logratio %>%
  filter(treatment %in% c("B", "T")) %>%
  pivot_wider(names_from = gen_group, values_from = mean_log_ratio) %>%
  mutate(R = G10 - G0) %>%
  dplyr::select(treatment, R)

G_values <- c(7, 8)

h2_results <- R %>%
  left_join(S, by = "treatment") %>%
  rowwise() %>%
  mutate(
    h2_G7 = R / (S * G_values[1]),
    h2_G8 = R / (S * G_values[2])
  ) %>%
  ungroup()

h2_results
