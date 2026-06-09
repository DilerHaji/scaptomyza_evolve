pmiss2 <- aggregate(tab[, c("X11.Feb", 
                            "X12.Feb",
                            "X13.Feb", 
                            "X14.Feb",
                            "X15.Feb", 
                            "X16.Feb", 
                            "X17.Feb",
                            "X18.Feb", 
                            "X19.Feb",
                            "X20.Feb")], by = list(tab$PopSpecies, tab$PhenotypeCageSpecies), FUN = function(x){ round((1-(length(x[is.na(x) | x == 0]) / length(x)))*100) })

colnames(pmiss2) <- c("PopSpecies", "PhenotypeCageSpecies", "X11.Feb", "X12.Feb", "X13.Feb", "X14.Feb", "X15.Feb", "X16.Feb", "X17.Feb", "X18.Feb", "X19.Feb", "X20.Feb")

pmiss2$tray_count <- aggregate(tab[, "Group"], by = list(tab$PopSpecies, tab$PhenotypeCageSpecies), FUN = function(x){ length(x[!is.na(x)]) })[,3]

pmiss2g <- gather(pmiss2, key = date, value = counts, 3:12)
pmiss2g$date2 <- sub(x = pmiss2g$date, "X", "")
pmiss2g$PopSpecies <- sub(x = pmiss2g$PopSpecies, "^B$", "Barbarea")
pmiss2g$PopSpecies <- sub(x = pmiss2g$PopSpecies, "^T$", "Turritus")
pmiss2g$PopSpecies <- sub(x = pmiss2g$PopSpecies, "^M$", "Mixture")
pmiss2g$PhenotypeCageSpecies <- sub(x = pmiss2g$PhenotypeCageSpecies, "^B$", "Tested on Barbarea")
pmiss2g$PhenotypeCageSpecies <- sub(x = pmiss2g$PhenotypeCageSpecies, "^T$", "Tested on Turritus")
pmiss2g$PhenotypeCageSpecies <- sub(x = pmiss2g$PhenotypeCageSpecies, "^M$", "Tested on Mixture")

pmiss2g <- pmiss2g %>%
  group_by(PopSpecies, PhenotypeCageSpecies) %>%
  arrange(date2) %>%
  mutate(cumulative_counts = cumsum(counts)) %>%
  mutate(total_counts = sum(counts)) %>%
  mutate(cumulative_percentage = cumulative_counts / total_counts * 100)


pmiss2g$normCount <- pmiss2g$counts/pmiss2g$tray_count
pmiss2g_bar <- gather(pmiss2g, key = var, value = val, which(colnames(pmiss2g) %in% c("tray_count", "counts", "normCount")))
pmiss2g_bar$var <- sub(x = pmiss2g_bar$var, "^counts$", "Number of Flies Emerging")
pmiss2g_bar$var <- sub(x = pmiss2g_bar$var, "^tray_count$", "Number of Plants")
pmiss2g_bar$var <- sub(x = pmiss2g_bar$var, "^normCount$", "Number of Flies Emerging Per Plant")


gbar <- ggplot(pmiss2g_bar, aes(x = PopSpecies, y = val, fill = var)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_grid(rows = vars(var), cols = vars(PhenotypeCageSpecies), scales = "free_y") + 
  theme_bw() + 
  labs(y = "Total Count", 
       x = "Evolution\n Treatment") + 
  scale_fill_brewer(palette = "Paired") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold"),
    strip.text.y = element_blank(),
    legend.title = element_blank()
  )

max_dates <- pmiss2g %>%
  group_by(PhenotypeCageSpecies, date2) %>%
  summarise(avg_counts = mean(counts), .groups = "drop") %>%
  group_by(PhenotypeCageSpecies) %>%
  slice_max(avg_counts, n = 1) %>%
  dplyr::select(PhenotypeCageSpecies, date2, avg_counts)


gtime <- ggplot(pmiss2g, aes(x = date2, y = counts, col = PopSpecies, group = PopSpecies)) + 
  geom_point(size = 2) + 
  geom_line() + 
  geom_vline(data = max_dates, aes(xintercept = date2), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap2(~PhenotypeCageSpecies, ncol = 1, strip = strip_themed(
    text_x = list(
      element_text(color = "#499FFF", face = "bold", size = 10),
      element_text(color = "#B01754", face = "bold", size = 10), 
      element_text(color = "#EDB72D", face = "bold", size = 10)
    )
  )) + 
  theme_bw() + 
  labs(col = "Treatment of\n evolved parents", 
       y = "Total number of F1 flies emerging", 
       x = "Emergence Date") + 
  scale_color_manual(values = custom_colors2) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.background = element_blank()
  )

pmiss2g$name <- paste(pmiss2g$PopSpecies,pmiss2g$PhenotypeCageSpecies)
g <- ggplot(pmiss2g, aes(x = date2, y = counts, col = PopSpecies, group = name)) + 
  geom_point(size = 2) + 
  geom_line() + 
  facet_wrap(~PhenotypeCageSpecies, ncol = 1) + 
  theme_bw() + 
  labs(col = "Evolved \nSource \nPopulation", 
       y = "Number of Flies Emerging", 
       x = "Collection Date", 
       shape = "Number of \nplants") + 
  scale_color_manual(values = custom_colors2) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold"), panel.grid = element_blank()
  )

ggsave(filename = "timep.png", plot = g, units = "in", width = 6, height = 6, scale = 0.8)


gcum <- ggplot(pmiss2g, aes(x = date2, y = cumulative_percentage, col = PhenotypeCageSpecies, group = interaction(PopSpecies, PhenotypeCageSpecies))) + 
  geom_point() + 
  geom_line() + 
  theme_bw() + 
  labs(col = "Test \nTreatment", 
       y = "Cumulative Percentage of Flies Emerging", 
       x = "Collection Date", 
       shape = "Number of \nplants") + 
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 75, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "black") + 
  scale_color_manual(values = custom_colors3) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold")
  )



gtdat <- gather(tab, key = date, value = val, 7:16)
gtdat <- aggregate(gtdat$numFlies, by = list(gtdat$Rep, gtdat$Group,gtdat$PopSpecies, gtdat$PhenotypeCageSpecies, gtdat$date), sum)
colnames(gtdat) <- c("Rep", "Tray", "PopSpecies", "PhenotypeCageSpecies", "date", "count")
gtdat$date2 <- sub(x = gtdat$date, "X", "")

gtdat <- gather(tab, key = date, value = val, 7:16)
gtdat$date2 <- sub(x = gtdat$date, "X", "")
gtdat$line <- paste(gtdat$Rep, gtdat$Group)
is.na(gtdat$val) <- 0
gtdat$PopSpecies <- sub(x = gtdat$PopSpecies, "^B$", "Barbarea")
gtdat$PopSpecies <- sub(x = gtdat$PopSpecies, "^T$", "Turritus")
gtdat$PopSpecies <- sub(x = gtdat$PopSpecies, "^M$", "Mixture")
gtdat$PhenotypeCageSpecies <- sub(x = gtdat$PhenotypeCageSpecies, "^B$", "Tested on Barbarea")
gtdat$PhenotypeCageSpecies <- sub(x = gtdat$PhenotypeCageSpecies, "^T$", "Tested on Turritus")
gtdat$PhenotypeCageSpecies <- sub(x = gtdat$PhenotypeCageSpecies, "^M$", "Tested on Mixture")

gtime2 <- ggplot(gtdat, aes(x = date2, y = val, col = PopSpecies, group = line)) + 
  geom_line(linewidth = 1, alpha = 0.5) + 
  facet_grid(rows = vars(PhenotypeCageSpecies), cols = vars(PopSpecies) ) + 
  theme_bw() + 
  labs(col = "Evolved \nSource \nPopulation", 
       y = "Number of Flies Emerging", 
       x = "Collection Date", 
       shape = "Number of \nplants") + 
  scale_color_manual(values = custom_colors2) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.background = element_blank(), 
    strip.text = element_text()
  )

write.csv(x = gtdat, file = "gtdat.csv", quote = F, row.names = F)


