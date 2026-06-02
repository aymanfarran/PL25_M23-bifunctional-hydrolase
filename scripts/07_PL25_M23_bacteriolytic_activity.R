library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

df <- tribble(
  ~Strain, ~Treatment, ~Replicate, ~CFU_ml,
  "B. subtilis", "Negative control", 1, 5e8,
  "B. subtilis", "Negative control", 2, 3e8,
  "B. subtilis", "Negative control", 3, 7e8,
  "B. subtilis", "PL25_M23",         1, 4e5,
  "B. subtilis", "PL25_M23",         2, 7e5,
  "B. subtilis", "PL25_M23",         3, 2e5,

  "S. warneri", "Negative control", 1, 5e8,
  "S. warneri", "Negative control", 2, 8e8,
  "S. warneri", "Negative control", 3, 4e8,
  "S. warneri", "PL25_M23",         1, 2e5,
  "S. warneri", "PL25_M23",         2, 8e5,
  "S. warneri", "PL25_M23",         3, 4e5,

  "MRSA", "Negative control", 1, 5e8,
  "MRSA", "Negative control", 2, 8e8,
  "MRSA", "Negative control", 3, 2e8,
  "MRSA", "PL25_M23",         1, 4e8,
  "MRSA", "PL25_M23",         2, 8e8,
  "MRSA", "PL25_M23",         3, 2e8,

  "PL25", "Negative control", 1, 5e8,
  "PL25", "Negative control", 2, 3e8,
  "PL25", "Negative control", 3, 7e8,
  "PL25", "PL25_M23",         1, 7e4,
  "PL25", "PL25_M23",         2, 3e4,
  "PL25", "PL25_M23",         3, 7e4
)

strain_order <- c("B. subtilis", "S. warneri", "MRSA", "PL25")
df$Strain <- factor(df$Strain, levels = strain_order)

summary_df <- df %>%
  group_by(Strain, Treatment) %>%
  summarise(
    mean_cfu = mean(CFU_ml),
    sd_cfu = sd(CFU_ml),
    .groups = "drop"
  )

reduction_df <- summary_df %>%
  select(Strain, Treatment, mean_cfu) %>%
  pivot_wider(names_from = Treatment, values_from = mean_cfu) %>%
  mutate(
    fold_reduction = `Negative control` / PL25_M23,
    log10_reduction = log10(fold_reduction)
  )

dir.create("results/19_bacteriolytic_activity", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

write.csv(df, "results/19_bacteriolytic_activity/PL25_M23_CFU_replicates.csv",
          row.names = FALSE)
write.csv(summary_df, "results/19_bacteriolytic_activity/PL25_M23_CFU_summary.csv",
          row.names = FALSE)
write.csv(reduction_df, "results/19_bacteriolytic_activity/PL25_M23_log10_reduction.csv",
          row.names = FALSE)

p1 <- ggplot(df, aes(Strain, CFU_ml, color = Treatment, shape = Treatment)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.58),
    size = 3,
    stroke = 0.8
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    aes(group = Treatment),
    position = position_dodge(width = 0.58),
    width = 0.36,
    linewidth = 0.55
  ) +
  scale_y_log10(
    limits = c(1e4, 1e9),
    breaks = 10^(4:9),
    labels = function(x) parse(text = paste0("10^", log10(x)))
  ) +
  scale_color_manual(values = c(
    "Negative control" = "#4D4D4D",
    "PL25_M23" = "#D62728"
  )) +
  labs(
    tag = "a",
    x = NULL,
    y = expression(CFU~mL^{-1}),
    color = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.tag = element_text(face = "bold", size = 15)
  )

p2 <- ggplot(reduction_df, aes(Strain, log10_reduction)) +
  geom_col(fill = "#2F6DB0", width = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f", log10_reduction)),
    vjust = -0.35,
    size = 3.7
  ) +
  scale_y_continuous(
    limits = c(0, 4.5),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    tag = "b",
    x = NULL,
    y = expression(Log[10]~reduction~"(CFU"~mL^{-1}*")")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.tag = element_text(face = "bold", size = 15))

figure <- p1 / p2 + plot_layout(heights = c(1.5, 1))

ggsave(
  "figures/09_PL25_M23_bacteriolytic_activity.png",
  figure,
  width = 7,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  "figures/09_PL25_M23_bacteriolytic_activity.pdf",
  figure,
  width = 7,
  height = 7,
  bg = "white"
)

print(reduction_df)
