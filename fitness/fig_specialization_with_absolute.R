library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)

rs <- read.csv("SelectionExperiment_Gen23andFinal_numFlies.csv")
rs$comb <- paste(rs$PlantSpecies, rs$Rep, sep = "")

m_baseline <- rs %>%
  group_by(PlantSpecies, Generation, comb) %>%
  summarise(numFlies = mean(numFlies), .groups = "drop") %>%
  filter(PlantSpecies == "M") %>%
  group_by(Generation) %>%
  summarise(M_mean = mean(numFlies), .groups = "drop")

rsr <- rs %>%
  left_join(m_baseline, by = "Generation") %>%
  mutate(
    RelativeFlies = numFlies / M_mean,
    GenBin = ifelse(Generation %in% c("2", "3"), "G2/G3", "G10"),
    GenBin = factor(GenBin, levels = c("G2/G3", "G10")),
    Treatment = recode(PlantSpecies, "B" = "Barbarea", "M" = "Mixture", "T" = "Turritus"),
    Treatment = factor(Treatment, levels = c("Barbarea", "Mixture", "Turritus"))
  )

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
fg <- rsr %>%
  filter(PlantSpecies %in% c("B", "T")) %>%
  group_by(Treatment, GenBin) %>%
  summarise(mean = mean(RelativeFlies), sem = se(RelativeFlies), .groups = "drop") %>%
  mutate(ci_lo = mean - 1.96 * sem, ci_hi = mean + 1.96 * sem)

bt_g10 <- rsr %>% filter(PlantSpecies %in% c("B", "T"), GenBin == "G10") %>% pull(RelativeFlies)
tt <- t.test(bt_g10, mu = 1)
pval <- tt$p.value
cohen_d <- (mean(bt_g10) - 1) / sd(bt_g10)

rel_lo <- 0.50
rel_hi <- 1.70
abs_lo <- 0
abs_hi <- 1000
to_rel <- function(x) rel_lo + (x - abs_lo) / (abs_hi - abs_lo) * (rel_hi - rel_lo)
to_abs <- function(x) abs_lo + (x - rel_lo) / (rel_hi - rel_lo) * (abs_hi - abs_lo)

trt_order <- c("Barbarea", "Mixture", "Turritus")
bar_width <- 0.22
dx <- seq(-bar_width, bar_width, length.out = 3)

bg <- rsr %>%
  group_by(GenBin, Treatment) %>%
  summarise(abs_mean = mean(numFlies), .groups = "drop") %>%
  mutate(
    xcenter = as.numeric(GenBin) + dx[match(Treatment, trt_order)],
    scaled_mean = to_rel(abs_mean)
  )

colors <- c("Barbarea" = "#499FFF", "Turritus" = "#EDB72D", "Mixture" = "#7FA878")

points_rel <- rsr %>% filter(PlantSpecies %in% c("B", "T"))

dw <- 0.45

g <- ggplot() +
  geom_rect(
    data = bg,
    aes(xmin = xcenter - bar_width/2, xmax = xcenter + bar_width/2,
        ymin = rel_lo, ymax = scaled_mean, fill = Treatment),
    alpha = 0.22, color = NA, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = colors["Mixture"], linewidth = 1.1) +
  geom_point(
    data = points_rel,
    aes(x = GenBin, y = RelativeFlies, color = Treatment),
    position = position_jitterdodge(dodge.width = dw, jitter.width = 0.08),
    alpha = 0.35, size = 1.3
  ) +
  geom_line(
    data = fg,
    aes(x = GenBin, y = mean, color = Treatment, group = Treatment),
    position = position_dodge(width = dw),
    linetype = "dashed", linewidth = 0.6
  ) +
  geom_errorbar(
    data = fg,
    aes(x = GenBin, ymin = ci_lo, ymax = ci_hi, color = Treatment),
    position = position_dodge(width = dw),
    width = 0.12, linewidth = 0.7
  ) +
  geom_point(
    data = fg,
    aes(x = GenBin, y = mean, color = Treatment),
    position = position_dodge(width = dw), size = 5
  ) +
  annotate("text", x = 1.68, y = 1.32, label = "T ÷ B+T",
           color = colors["Turritus"], fontface = "bold", size = 3.4) +
  annotate("text", x = 1.68, y = 1.16, label = "B ÷ B+T",
           color = colors["Barbarea"], fontface = "bold", size = 3.4) +
  annotate("text", x = 1.52, y = 0.96, label = "B+T ÷ B+T",
           color = colors["Mixture"], fontface = "bold", size = 3.4) +
  annotate("segment", x = 2.22, xend = 2.22, y = 1.19, yend = 1.49,
           linewidth = 0.7, color = "grey30") +
  annotate("segment", x = 2.22, xend = 2.18, y = 1.19, yend = 1.19,
           linewidth = 0.7, color = "grey30") +
  annotate("segment", x = 2.22, xend = 2.18, y = 1.49, yend = 1.49,
           linewidth = 0.7, color = "grey30") +
  annotate("text", x = 2.24, y = 1.34, label = "Uniform\nSelection",
           hjust = 0, size = 3.1, lineheight = 0.9) +
  annotate("point", x = 2, y = 1.0, size = 1.3, color = "black") +
  annotate("segment", x = 2.03, xend = 2.24, y = 1.00, yend = 0.87,
           linewidth = 0.35, linetype = "dotted", color = "grey30") +
  annotate("text", x = 2.24, y = 0.83, label = "Subdivided\nSelection",
           hjust = 0, size = 3.1, lineheight = 0.9) +
  annotate("text", x = 0.55, y = 1.66,
           label = "Selection for Specialization",
           hjust = 0, size = 3.8, fontface = "bold") +
  annotate("text", x = 0.55, y = 1.60,
           label = sprintf("(one-sample t vs B+T: p = %.3f, d = %.2f)", pval, cohen_d),
           hjust = 0, size = 3.2) +
  annotate("text", x = 0.55, y = 1.53,
           label = "Error bars: 95% CI (1.96 × SEM)",
           hjust = 0, size = 2.7, fontface = "italic") +
  annotate("text", x = 0.55, y = 0.55,
           label = "Background bars: mean absolute F1 counts per replicate (right axis)",
           hjust = 0, size = 2.6, fontface = "italic", color = "grey40") +
  annotate("text", x = 0.35, y = 1.70, label = "C",
           hjust = 0, vjust = 1, size = 8, fontface = "bold") +
  scale_y_continuous(
    name = "Population Viability (Relative to B+T)",
    limits = c(rel_lo, rel_hi),
    breaks = c(0.75, 1.00, 1.25, 1.50),
    sec.axis = sec_axis(~ to_abs(.), name = "Absolute F1 adults (per replicate)",
                        breaks = c(0, 250, 500, 750, 1000))
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.35))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey60", linewidth = 0.8),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.y.left = element_text(margin = margin(r = 8)),
    axis.title.y.right = element_text(margin = margin(l = 8), color = "grey40"),
    axis.text.y.right = element_text(color = "grey40"),
    legend.position = "none",
    plot.margin = margin(10, 10, 20, 10)
  ) +
  coord_cartesian(clip = "off")

g_final <- g +
  annotation_custom(
    textGrob("Evolution Start", gp = gpar(fontsize = 9, col = "grey30")),
    xmin = 1, xmax = 1, ymin = 0.38, ymax = 0.38
  ) +
  annotation_custom(
    textGrob("Evolution End", gp = gpar(fontsize = 9, col = "grey30")),
    xmin = 2, xmax = 2, ymin = 0.38, ymax = 0.38
  )

cat(sprintf("One-sample t vs 1: t=%.2f df=%d p=%.4f d=%.2f\n",
            tt$statistic, tt$parameter, pval, cohen_d))

ggsave("fig_specialization_with_absolute.png", g_final,
       width = 7.2, height = 5.6, dpi = 300)
ggsave("fig_specialization_with_absolute.svg", g_final,
       width = 7.2, height = 5.6)

strip_svg_masks <- function(path) {
  library(xml2)
  doc <- read_xml(path)
  ns <- c(svg = "http://www.w3.org/2000/svg")
  for (node in xml_find_all(doc, ".//svg:clipPath | .//svg:mask", ns)) xml_remove(node)
  for (node in xml_find_all(doc, "//*[@clip-path]")) xml_attr(node, "clip-path") <- NULL
  for (node in xml_find_all(doc, "//*[@mask]")) xml_attr(node, "mask") <- NULL
  write_xml(doc, path)
}
strip_svg_masks("fig_specialization_with_absolute.svg")

print(g_final)
