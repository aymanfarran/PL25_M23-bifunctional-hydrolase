#!/usr/bin/env Rscript
# Main wet-lab figures for PL25_M23.
#
#   Figure 1  Antibacterial activity
#     a  spectrum at 100 ug/mL, all seven strains      dot plot, delta log10
#     b  dose-response, 0-300 ug/mL                    line plot, delta log10
#     c  time-kill, 0-180 min at 100 ug/mL             line plot, delta log10
#
#   Figure 2  Metal dependence and salt tolerance
#     a  PL25_M23 +/- EDTA +/- ZnCl2                   paired dot plot
#     b  salt tolerance, 0-1.0 M NaCl, five strains    line plot, log10 CFU/mL
#
# Panel 2b plots treated log10 CFU/mL rather than a reduction: untreated
# controls exist only at 0 and 1.0 M NaCl, so a matched reduction cannot be
# computed at 0.25 and 0.5 M.
#
# The experimental unit is the biological replicate. Technical duplicates were
# averaged before log transformation upstream (50_activity_read.R), giving
# n = 3 per condition.
#
# Statistics
#   1a  paired t-test, treated vs untreated
#   1b  linear mixed model on delta log10 with strain and dose as fixed effects
#       and biological replicate as a random intercept, Dunnett contrasts
#       against 0 ug/mL within each strain
#   1c  as 1b, against 0 min
#   2a  two planned contrasts only: PL25_M23 vs +EDTA, and +EDTA vs +EDTA+Zn
#   2b  paired t-test of each salt level against 0 M
#
# Usage:  Rscript scripts/51_activity_figures.R
#
# Inputs:  results/21_activity_7strain/*.tsv   (from 50_activity_read.R)
# Outputs: figures/F1_antibacterial_activity.{pdf,png}
#          figures/F2_metal_salt.{pdf,png}
#          results/21_activity_7strain/figure_stats.tsv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(patchwork); library(stringr)
  library(emmeans)
})
emm_options(msg.interaction = FALSE, msg.nesting = FALSE)

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
salt_strains <- strain_levels[3:7]   # 43079 and 49619 do not grow at 1.0 M NaCl

# Axis labels: ATCC designations, with the genus italicised.
strain_labels <- c(
  "Streptococcus equi subsp. zooepidemicus" = "italic('S. equi')~'43079'",
  "Streptococcus pneumoniae"                = "italic('S. pneumoniae')~'49619'",
  "Staphylococcus aureus (MRSA)"            = "'MRSA BAA-1026'",
  "Staphylococcus aureus (MSSA)"            = "'MSSA 29213'",
  "Mammaliicoccus sciuri"                   = "italic('M. sciuri')~'29061'",
  "Staphylococcus aureus"                   = "italic('S. aureus')~'BAA-976'",
  "Virgibacillus salarius"                  = "italic('V. salarius')~'PL25'")
lab_parse <- function(x) parse(text = strain_labels[x])

SHAPES <- c(21, 22, 23, 24, 25, 15, 17)
FILLS  <- c("#FFFFFF", "#D9D9D9", "#B0B0B0", "#8A8A8A", "#636363", "#3B3B3B", "#000000")
names(SHAPES) <- names(FILLS) <- strain_levels

theme_pub <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(
      axis.text  = element_text(colour = INK, size = base),
      axis.title = element_text(colour = INK, size = base + 1),
      axis.line  = element_line(colour = INK, linewidth = 0.35),
      axis.ticks = element_line(colour = INK, linewidth = 0.35),
      panel.grid.major.y = element_line(colour = GREY_GRID, linewidth = 0.3),
      legend.key.size = unit(11, "pt"),
      legend.text  = element_text(size = base - 0.5, colour = INK, hjust = 0),
      legend.title = element_blank(),
      legend.background = element_blank(),
      plot.tag = element_text(face = "bold", size = base + 5, colour = INK))
}

stars <- function(p) ifelse(is.na(p), "",
  ifelse(p <= 1e-4, "****", ifelse(p <= 1e-3, "***",
  ifelse(p <= 1e-2, "**", ifelse(p <= 0.05, "*", "ns")))))
paired_p <- function(a, b) tryCatch(t.test(a, b, paired = TRUE)$p.value,
                                    error = function(e) NA_real_)

# Every panel reports several tests at once, so p-values are Holm-adjusted
# across all comparisons shown in that panel. Symbols on the figures and the
# values in figure_stats.tsv are the adjusted ones; the raw p-values are kept
# alongside them for transparency.
holm <- function(p) p.adjust(p, method = "holm")

stat_log <- list()
note <- function(...) stat_log[[length(stat_log) + 1]] <<- tibble(...)

fct  <- function(x, lv = strain_levels) factor(x, levels = lv)
dlog <- function(d) mutate(d, red = log10_untreated - log10_treated)

primary  <- read_tsv(file.path(dat, "primary_raw.tsv"),  show_col_types = FALSE)
dose     <- read_tsv(file.path(dat, "dose_raw.tsv"),     show_col_types = FALSE)
timekill <- read_tsv(file.path(dat, "timekill_raw.tsv"), show_col_types = FALSE)
edta     <- read_tsv(file.path(dat, "edta_raw.tsv"),     show_col_types = FALSE)
salt     <- read_tsv(file.path(dat, "salt_raw.tsv"),     show_col_types = FALSE)

# ═════════════════════════════════════════════════════════════════════
# Figure 1a — antibacterial spectrum at 100 ug/mL
# ═════════════════════════════════════════════════════════════════════
spec <- dlog(primary) |> mutate(strain = fct(strain))
spec_stats <- spec |> group_by(strain) |>
  summarise(m = mean(red), sd = sd(red), top = max(red),
            p_raw = paired_p(log10_untreated, log10_treated), .groups = "drop") |>
  mutate(p = holm(p_raw), sym = stars(p))          # 7 tests in this panel
pwalk(spec_stats, \(strain, m, sd, top, p_raw, p, sym)
      note(figure = "1", panel = "a", strain = strain,
           comparison = "treated vs untreated at 100 ug/mL",
           estimate = m, p_raw = p_raw, p_holm = p, symbol = sym))

# The dashed reference line marks a 3-log reduction, the conventional
# bactericidal threshold, and is labelled in every panel that carries it.
three_log <- function(x, y = 3.13, size = 2.4) {
  list(geom_hline(yintercept = 3, linetype = "22", linewidth = 0.3, colour = "#888888"),
       annotate("text", x = x, y = y, label = "3-log reduction (99.9%)",
                size = size, colour = "#666666", hjust = 1, vjust = 0))
}

p1a <- ggplot(spec, aes(strain, red)) +
  three_log(x = 7.45) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymin = after_stat(y), ymax = after_stat(y)), width = 0.45,
               linewidth = 0.4, colour = INK) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
               width = 0.16, linewidth = 0.3, colour = INK) +
  geom_point(aes(shape = strain, fill = strain), size = 2, colour = INK,
             stroke = 0.4, position = position_jitter(width = 0.09, height = 0, seed = 1)) +
  geom_text(data = spec_stats, aes(strain, top + 0.42, label = sym),
            inherit.aes = FALSE, size = 3, fontface = "bold", colour = INK) +
  scale_shape_manual(values = SHAPES, guide = "none") +
  scale_fill_manual(values = FILLS, guide = "none") +
  scale_x_discrete(labels = lab_parse) +
  scale_y_continuous(expression("Log"[10]~"CFU/mL reduction vs untreated"),
                     limits = c(0, 5.0), breaks = 0:5) +
  theme_pub() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

# ═════════════════════════════════════════════════════════════════════
# Figures 1b and 1c — dose-response and time-kill
# ═════════════════════════════════════════════════════════════════════
# Both are fitted the same way: a repeated-measures model on delta log10 with
# strain and level as crossed fixed effects and biological replicate as a
# fixed block within strain, then Dunnett contrasts against the zero level
# within each strain.
#
# Replicate is a fixed block rather than a random intercept because at n = 3
# the replicate variance is estimated as zero and the random-effects fit is
# singular. The response is already a within-replicate difference against the
# matched untreated control, so most between-replicate variation has been
# removed before modelling.
level_panel <- function(df, xvar, xlab, breaks, fig, panel) {
  d <- dlog(df) |> mutate(strain = fct(strain), lvl = factor(.data[[xvar]]),
                          bio = factor(bio_rep))
  fit <- lm(red ~ strain * lvl + strain:bio, data = d)
  # Contrasts are taken unadjusted, then Holm-corrected across every
  # comparison in the panel rather than only within each strain.
  cmp <- emmeans(fit, ~ lvl | strain) |>
    contrast("trt.vs.ctrl", ref = 1, adjust = "none") |>
    summary(infer = TRUE) |> as_tibble() |>
    mutate(p_holm = holm(p.value))
  pwalk(cmp, \(contrast, strain, estimate, p.value, p_holm, ...)
        note(figure = fig, panel = panel, strain = as.character(strain),
             comparison = paste(xvar, contrast), estimate = estimate,
             p_raw = p.value, p_holm = p_holm, symbol = stars(p_holm)))

  ggplot(d, aes(.data[[xvar]], red, shape = strain, fill = strain)) +
    three_log(x = max(breaks), size = 2.2) +
    geom_point(size = 0.85, colour = "#9A9A9A", stroke = 0.25,
               position = position_jitter(width = 0, height = 0.04, seed = 2)) +
    stat_summary(fun = mean, geom = "line", aes(group = strain),
                 linewidth = 0.45, colour = INK) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
                 width = 0, linewidth = 0.3, colour = INK) +
    stat_summary(fun = mean, geom = "point", size = 1.9, colour = INK, stroke = 0.35) +
    scale_shape_manual(values = SHAPES, labels = lab_parse) +
    scale_fill_manual(values = FILLS, labels = lab_parse) +
    scale_x_continuous(xlab, breaks = breaks) +
    scale_y_continuous(expression("Log"[10]~"CFU/mL reduction vs untreated"),
                       limits = c(-0.6, 5.4)) +
    theme_pub()
}

p1b <- level_panel(dose, "dose_ug_mL", expression("PL25_M23 ("*mu*"g/mL)"),
                   c(0, 50, 100, 300), "1", "b")
p1c <- level_panel(timekill, "time_min", "Time (min)",
                   c(0, 30, 60, 120, 180), "1", "c")

fig1 <- p1a / (p1b | p1c) +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(figs, "F1_antibacterial_activity.pdf"), fig1, width = 8.2, height = 7.2)
ggsave(file.path(figs, "F1_antibacterial_activity.png"), fig1, width = 8.2, height = 7.2,
       dpi = 300, bg = "white")

# ═════════════════════════════════════════════════════════════════════
# Figure 2a — metal dependence
# ═════════════════════════════════════════════════════════════════════
# Main panel shows only the four enzyme-containing conditions. The untreated,
# EDTA-only and Zn-only controls establish that the additives are not
# themselves bactericidal and are reported in Supplementary Figure S1.
ed_levels <- c("PL25_M23", "PL25_M23 + EDTA", "PL25_M23 + EDTA + Zn", "PL25_M23 + Zn")
ed_wrap   <- c("PL25_M23", "+ EDTA", "+ EDTA\n+ Zn", "+ Zn")

ed <- edta |>
  group_by(strain, bio_rep) |>
  mutate(red = log10_cfu[condition == "Untreated"] - log10_cfu) |>
  ungroup() |> filter(condition %in% ed_levels) |>
  mutate(condition = factor(condition, ed_levels), strain = fct(strain))

# Two planned contrasts: does chelation remove activity, and does zinc
# restore it. Nothing else is tested.
wide <- ed |> select(strain, bio_rep, condition, red) |>
  pivot_wider(names_from = condition, values_from = red)
ed_stats <- wide |> group_by(strain) |>
  summarise(raw_edta = paired_p(`PL25_M23`, `PL25_M23 + EDTA`),
            raw_resc = paired_p(`PL25_M23 + EDTA`, `PL25_M23 + EDTA + Zn`),
            .groups = "drop")
# Holm across all 14 comparisons drawn in this panel, both contrasts together.
adj <- holm(c(ed_stats$raw_edta, ed_stats$raw_resc))
ed_stats <- ed_stats |>
  mutate(p_edta = adj[seq_len(n())], p_resc = adj[n() + seq_len(n())])
pwalk(ed_stats, \(strain, raw_edta, raw_resc, p_edta, p_resc) {
  note(figure = "2", panel = "a", strain = strain,
       comparison = "PL25_M23 vs +EDTA", estimate = NA_real_,
       p_raw = raw_edta, p_holm = p_edta, symbol = stars(p_edta))
  note(figure = "2", panel = "a", strain = strain,
       comparison = "+EDTA vs +EDTA+Zn", estimate = NA_real_,
       p_raw = raw_resc, p_holm = p_resc, symbol = stars(p_resc))
})

# Brackets carrying the two contrasts, drawn above each strain facet:
# PL25_M23 -> +EDTA at the lower level, +EDTA -> +EDTA+Zn above it.
brk <- ed_stats |>
  left_join(summarise(group_by(ed, strain), t = max(red), .groups = "drop"), "strain") |>
  transmute(strain,
            sym_edta = stars(p_edta), sym_resc = stars(p_resc),
            y1 = t + 0.35, y2 = t + 0.95)

p2a <- ggplot(ed, aes(condition, red)) +
  geom_hline(yintercept = 3, linetype = "22", linewidth = 0.3, colour = "#888888") +
  geom_line(aes(group = bio_rep), colour = "#AAAAAA", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymin = after_stat(y), ymax = after_stat(y)), width = 0.5,
               linewidth = 0.4, colour = INK) +
  geom_point(aes(shape = strain, fill = strain), size = 1.5, colour = INK, stroke = 0.35) +
  geom_segment(data = brk, aes(x = 1, xend = 2, y = y1, yend = y1),
               inherit.aes = FALSE, linewidth = 0.25, colour = INK) +
  geom_text(data = brk, aes(1.5, y1 + 0.06, label = sym_edta),
            inherit.aes = FALSE, size = 2.2, colour = INK, vjust = 0) +
  geom_segment(data = brk, aes(x = 2, xend = 3, y = y2, yend = y2),
               inherit.aes = FALSE, linewidth = 0.25, colour = INK) +
  geom_text(data = brk, aes(2.5, y2 + 0.06, label = sym_resc),
            inherit.aes = FALSE, size = 2.2, colour = INK, vjust = 0) +
  facet_wrap(~ strain, nrow = 2, labeller = as_labeller(strain_labels, label_parsed)) +
  scale_shape_manual(values = SHAPES, guide = "none") +
  scale_fill_manual(values = FILLS, guide = "none") +
  scale_x_discrete(labels = ed_wrap) +
  scale_y_continuous(expression("Log"[10]~"CFU/mL reduction vs untreated"), limits = c(-0.4, 5.6)) +
  theme_pub() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6.4, angle = 30, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.8, colour = INK))

# ═════════════════════════════════════════════════════════════════════
# Figure 2b — salt tolerance
# ═════════════════════════════════════════════════════════════════════
# Treated log10 CFU/mL, not a reduction: untreated controls were run only at
# 0 and 1.0 M NaCl, so a salt-matched reduction cannot be computed at 0.25 or
# 0.5 M. Comparing treated viability across salt levels needs no untreated
# control, so each level is still tested against 0 M NaCl.
salt_tr <- salt |> filter(treated, strain %in% salt_strains) |>
  mutate(strain = fct(strain, salt_strains))
salt_ref <- salt_tr |> filter(NaCl_M == 0) |> select(strain, bio_rep, ref = log10_cfu)
salt_stats <- salt_tr |> filter(NaCl_M > 0) |>
  left_join(salt_ref, c("strain", "bio_rep")) |>
  group_by(strain, NaCl_M) |>
  summarise(p_raw = paired_p(log10_cfu, ref), m = mean(log10_cfu), .groups = "drop") |>
  mutate(p = holm(p_raw))                          # 15 tests in this panel
pwalk(salt_stats, \(strain, NaCl_M, p_raw, m, p)
      note(figure = "2", panel = "b", strain = strain,
           comparison = paste(NaCl_M, "M NaCl vs 0 M"), estimate = m,
           p_raw = p_raw, p_holm = p, symbol = stars(p)))

# Only significant comparisons are marked. Fourteen of the fifteen are ns
# after correction, and labelling them all would bury the one that is not;
# the caption states that every unmarked comparison was ns.
salt_marks <- salt_stats |> filter(stars(p) != "ns") |> mutate(sym = stars(p))

p2b <- ggplot(salt_tr, aes(NaCl_M, log10_cfu, shape = strain, fill = strain)) +
  geom_point(size = 0.85, colour = "#9A9A9A", stroke = 0.25) +
  stat_summary(fun = mean, geom = "line", aes(group = strain),
               linewidth = 0.45, colour = INK) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
               width = 0, linewidth = 0.3, colour = INK) +
  stat_summary(fun = mean, geom = "point", size = 2, colour = INK, stroke = 0.35) +
  geom_text(data = salt_marks, aes(NaCl_M, m + 0.12, label = sym),
            inherit.aes = FALSE, size = 3, fontface = "bold", colour = INK, vjust = 0) +
  scale_shape_manual(values = SHAPES[salt_strains], labels = lab_parse) +
  scale_fill_manual(values = FILLS[salt_strains], labels = lab_parse) +
  scale_x_continuous("NaCl (M)", breaks = c(0, 0.25, 0.5, 1.0)) +
  scale_y_continuous(expression("Residual viable bacteria, log"[10]~"CFU/mL")) +
  theme_pub() +
  theme(legend.position = "right")

fig2 <- (p2a | p2b) + plot_layout(widths = c(1.55, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(figs, "F2_metal_salt.pdf"), fig2, width = 10.4, height = 4.6)
ggsave(file.path(figs, "F2_metal_salt.png"), fig2, width = 10.4, height = 4.6,
       dpi = 300, bg = "white")

all_stats <- bind_rows(stat_log) |> arrange(figure, panel, strain)
write.table(all_stats, file.path(dat, "figure_stats.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("wrote figures/F1_antibacterial_activity.{pdf,png}\n")
cat("wrote figures/F2_metal_salt.{pdf,png}\n")
cat(sprintf("wrote results/21_activity_7strain/figure_stats.tsv (%d tests)\n", nrow(all_stats)))
