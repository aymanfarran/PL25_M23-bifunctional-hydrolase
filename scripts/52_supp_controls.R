#!/usr/bin/env Rscript
# Supplementary Figure S1 — assay controls and the salt-assay exclusion.
#
#   a  additive-only controls: 1 mM EDTA and 1 mM ZnCl2 without enzyme, to
#      show that neither additive is bactericidal on its own and so cannot
#      account for the killing in Figure 2a
#   b  untreated viability at 0 vs 1.0 M NaCl, all seven strains, justifying
#      the exclusion of ATCC 43079 and ATCC 49619 from the salt assay
#   c  mean log10 CFU reduction across every strain and condition tested
#
# Also writes Supplementary Table S1: every biological-replicate value behind
# Figures 1, 2 and S1, in long form.
#
# Usage:  Rscript scripts/52_supp_controls.R
#
# Inputs:  results/21_activity_7strain/*.tsv   (from 50_activity_read.R)
# Outputs: figures/S1_activity_controls.{pdf,png}
#          results/21_activity_7strain/TableS1_all_values.tsv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(patchwork); library(stringr)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
root <- normalizePath(file.path(here, ".."))
dat  <- file.path(root, "results", "21_activity_7strain")
figs <- file.path(root, "figures")
dir.create(figs, showWarnings = FALSE)

INK <- "black"; GREY_GRID <- "#E8E8E8"

strain_levels <- c(
  "Streptococcus equi subsp. zooepidemicus", "Streptococcus pneumoniae",
  "Staphylococcus aureus (MRSA)", "Staphylococcus aureus (MSSA)",
  "Mammaliicoccus sciuri", "Staphylococcus aureus", "Virgibacillus salarius")
strain_labels <- c(
  "Streptococcus equi subsp. zooepidemicus" = "italic('S. equi')~'43079'",
  "Streptococcus pneumoniae"                = "italic('S. pneumoniae')~'49619'",
  "Staphylococcus aureus (MRSA)"            = "'MRSA BAA-1026'",
  "Staphylococcus aureus (MSSA)"            = "'MSSA 29213'",
  "Mammaliicoccus sciuri"                   = "italic('M. sciuri')~'29061'",
  "Staphylococcus aureus"                   = "italic('S. aureus')~'BAA-976'",
  "Virgibacillus salarius"                  = "italic('V. salarius')~'PL25'")
lab_parse <- function(x) parse(text = strain_labels[x])
fct <- function(x) factor(x, levels = strain_levels)

theme_pub <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(
      axis.text  = element_text(colour = INK, size = base),
      axis.title = element_text(colour = INK, size = base + 1),
      axis.line  = element_line(colour = INK, linewidth = 0.35),
      axis.ticks = element_line(colour = INK, linewidth = 0.35),
      panel.grid.major.y = element_line(colour = GREY_GRID, linewidth = 0.3),
      legend.key.size = unit(11, "pt"),
      legend.text  = element_text(size = base - 0.5, colour = INK),
      legend.title = element_blank(),
      legend.background = element_blank(),
      plot.tag = element_text(face = "bold", size = base + 5, colour = INK))
}

primary  <- read_tsv(file.path(dat, "primary_raw.tsv"),  show_col_types = FALSE)
dose     <- read_tsv(file.path(dat, "dose_raw.tsv"),     show_col_types = FALSE)
timekill <- read_tsv(file.path(dat, "timekill_raw.tsv"), show_col_types = FALSE)
edta     <- read_tsv(file.path(dat, "edta_raw.tsv"),     show_col_types = FALSE)
salt     <- read_tsv(file.path(dat, "salt_raw.tsv"),     show_col_types = FALSE)

# ── S1a  additive-only controls ──────────────────────────────────────
# Reduction relative to untreated, within each biological replicate. Values
# near zero mean the additive alone did not kill.
ctrl <- edta |>
  group_by(strain, bio_rep) |>
  mutate(red = log10_cfu[condition == "Untreated"] - log10_cfu) |>
  ungroup() |>
  filter(condition %in% c("EDTA only", "Zn only")) |>
  mutate(strain = fct(strain),
         condition = factor(condition, c("EDTA only", "Zn only")))

# The default PDF font has no micro sign or subscript two, so every label
# carrying them is built as a plotmath expression rather than literal text.
ctrl_legend <- c(expression("1 mM EDTA only"), expression("1 mM ZnCl"[2]*" only"))

pS1a <- ggplot(ctrl, aes(strain, red, shape = condition)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#888888") +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymin = after_stat(y), ymax = after_stat(y)), width = 0.5,
               linewidth = 0.4, colour = INK,
               position = position_dodge(0.6)) +
  geom_point(size = 1.7, fill = "white", colour = INK, stroke = 0.4,
             position = position_dodge(0.6)) +
  scale_shape_manual(values = c(21, 24), labels = ctrl_legend) +
  scale_x_discrete(labels = lab_parse) +
  scale_y_continuous(expression("Log"[10]~"CFU/mL reduction vs untreated"),
                     limits = c(-0.8, 1.2)) +
  theme_pub() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = c(0.22, 0.92), legend.direction = "horizontal")

# ── S1b  salt-assay exclusion ────────────────────────────────────────
excl <- salt |> filter(!treated) |>
  mutate(strain = fct(strain), NaCl = factor(NaCl_M, c(0, 1), c("0 M", "1.0 M")))
excl_drop <- excl |> group_by(strain) |>
  summarise(d = mean(log10_cfu[NaCl == "0 M"]) - mean(log10_cfu[NaCl == "1.0 M"]),
            .groups = "drop") |>
  mutate(lab = ifelse(d > 1, sprintf("-%.1f log", d), ""))

pS1b <- ggplot(excl, aes(strain, log10_cfu, fill = NaCl)) +
  stat_summary(fun = mean, geom = "col", position = position_dodge(0.72),
               width = 0.64, colour = INK, linewidth = 0.3) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
               position = position_dodge(0.72), width = 0.2,
               linewidth = 0.3, colour = INK) +
  geom_point(position = position_dodge(0.72), shape = 21, size = 1.1,
             fill = "white", colour = INK, stroke = 0.35) +
  geom_text(data = excl_drop, aes(strain, 1.1, label = lab), inherit.aes = FALSE,
            size = 2.5, colour = INK, fontface = "bold") +
  scale_fill_manual(values = c("0 M" = "#D9D9D9", "1.0 M" = "#4D4D4D")) +
  scale_x_discrete(labels = lab_parse) +
  scale_y_continuous(expression("Untreated viability, log"[10]~"CFU/mL"),
                     limits = c(0, 9.4), expand = expansion(0)) +
  theme_pub() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = c(0.62, 0.95), legend.direction = "horizontal")

# ── S1c  full condition x strain summary ─────────────────────────────
grid <- bind_rows(
  dlog_named <- primary |>
    transmute(strain, red = log10_untreated - log10_treated,
              cond = "100 ug/mL"),
  dose |> filter(dose_ug_mL > 0) |>
    transmute(strain, red = log10_untreated - log10_treated,
              cond = paste0(dose_ug_mL, " ug/mL")),
  timekill |> filter(time_min > 0) |>
    transmute(strain, red = log10_untreated - log10_treated,
              cond = paste0(time_min, " min")),
  edta |> group_by(strain, bio_rep) |>
    mutate(red = log10_cfu[condition == "Untreated"] - log10_cfu) |>
    ungroup() |> filter(condition != "Untreated") |>
    transmute(strain, red, cond = condition)) |>
  group_by(strain, cond) |>
  summarise(red = mean(red), .groups = "drop") |>
  mutate(strain = fct(strain),
         cond = factor(cond, c("50 ug/mL", "100 ug/mL", "300 ug/mL",
                               "30 min", "60 min", "120 min", "180 min",
                               "EDTA only", "Zn only", "PL25_M23",
                               "PL25_M23 + EDTA", "PL25_M23 + Zn",
                               "PL25_M23 + EDTA + Zn")))

# Dose labels need the micro sign, so they too are plotmath expressions.
cond_labels <- c(
  "50 ug/mL"  = expression(50~mu*"g/mL"),
  "100 ug/mL" = expression(100~mu*"g/mL"),
  "300 ug/mL" = expression(300~mu*"g/mL"),
  "30 min" = "30 min", "60 min" = "60 min",
  "120 min" = "120 min", "180 min" = "180 min",
  "EDTA only" = "EDTA only", "Zn only" = "Zn only",
  "PL25_M23" = "PL25_M23",
  "PL25_M23 + EDTA" = "PL25_M23 + EDTA",
  "PL25_M23 + Zn" = "PL25_M23 + Zn",
  "PL25_M23 + EDTA + Zn" = "PL25_M23 + EDTA + Zn")

pS1c <- ggplot(grid, aes(cond, strain, fill = red)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", red),
                colour = red > 2.4), size = 2.2, show.legend = FALSE) +
  scale_fill_gradient(low = "#FFFFFF", high = "#1A1A1A",
                      name = expression("Log"[10]~"reduction"),
                      limits = c(0, 5)) +
  scale_colour_manual(values = c("FALSE" = INK, "TRUE" = "white")) +
  scale_x_discrete(labels = cond_labels) +
  scale_y_discrete(labels = lab_parse, limits = rev(strain_levels)) +
  theme_pub() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        legend.position = "right", legend.title = element_text(size = 8))

figS1 <- (pS1a | pS1b) / pS1c +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(figs, "S1_activity_controls.pdf"), figS1, width = 10, height = 8)
ggsave(file.path(figs, "S1_activity_controls.png"), figS1, width = 10, height = 8,
       dpi = 300, bg = "white")

# ── Supplementary Table S1 ───────────────────────────────────────────
tbl <- bind_rows(
  primary  |> transmute(assay = "primary", strain, bio_rep,
                        condition = "PL25_M23 100 ug/mL",
                        log10_untreated, log10_treated),
  dose     |> transmute(assay = "dose", strain, bio_rep,
                        condition = paste0("PL25_M23 ", dose_ug_mL, " ug/mL"),
                        log10_untreated, log10_treated),
  timekill |> transmute(assay = "time-kill", strain, bio_rep,
                        condition = paste0(time_min, " min at 100 ug/mL"),
                        log10_untreated, log10_treated),
  edta     |> transmute(assay = "metal dependence", strain, bio_rep,
                        condition, log10_untreated = NA_real_,
                        log10_treated = log10_cfu),
  salt     |> transmute(assay = "salt tolerance", strain, bio_rep,
                        condition = paste0(ifelse(treated, "PL25_M23 ", "untreated "),
                                           NaCl_M, " M NaCl"),
                        log10_untreated = NA_real_, log10_treated = log10_cfu)) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

write.table(tbl, file.path(dat, "TableS1_all_values.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

cat("wrote figures/S1_activity_controls.{pdf,png}\n")
cat(sprintf("wrote results/21_activity_7strain/TableS1_all_values.tsv (%d rows)\n", nrow(tbl)))
