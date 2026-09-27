#!/usr/bin/env Rscript
# Read the wet-lab workbook and write one tidy TSV per assay.
#
# The workbook stores each assay as three biological-replicate blocks of the
# same seven strains, with two technical replicates per condition. Blocks are
# stacked vertically in the three killing sheets and placed side by side in the
# EDTA and salt sheets, so each assay needs its own geometry.
#
# Conditions: PL25_M23 against OD600-normalised cells in PBS pH 7.4 at 37 degC,
# 60 min unless the assay varies time, Miles-Misra drop-plate enumeration.
#
# Usage:  Rscript scripts/50_activity_read.R [path/to/M23_New_Data.xlsx]
#
# Outputs (results/21_activity_7strain/):
#   primary_raw.tsv  dose_raw.tsv  timekill_raw.tsv  edta_raw.tsv  salt_raw.tsv

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(purrr); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
root <- normalizePath(file.path(here, ".."))
src  <- if (length(args)) args[1] else file.path(root, "M23_New_Data.xlsx")
out  <- file.path(root, "results", "21_activity_7strain")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Worksheet strain label -> the label used downstream. Order matters: it is the
# order the strains appear in every sheet.
strain_levels <- c(
  "Streptococcus equi subsp. zooepidemicus",
  "Streptococcus pneumoniae",
  "Staphylococcus aureus (MRSA)",
  "Staphylococcus aureus (MSSA)",
  "Mammaliicoccus sciuri",
  "Staphylococcus aureus",
  "Virgibacillus salarius"
)

# The workbook uses non-breaking spaces inside some labels.
norm <- function(x) str_squish(str_replace_all(as.character(x), " ", " "))

read_raw <- function(sheet) {
  suppressMessages(read_excel(src, sheet = sheet, col_names = FALSE,
                              .name_repair = "minimal")) |> as.data.frame()
}

# Collapse the technical replicates in `cols` of one row to a single value for
# that biological replicate: average the raw CFU counts first, then take log10.
# The experimental unit is the biological replicate, so the technical duplicates
# are averaged rather than treated as independent observations.
tech_log <- function(df, row, cols) {
  v <- suppressWarnings(as.numeric(unlist(df[row, cols])))
  v <- v[!is.na(v)]
  if (!length(v)) NA_real_ else log10(mean(v))
}

# ---- Vertically stacked sheets: primary, dose, time-kill -------------
# Each block holds the seven strains in order, two rows each
# (untreated then treated).
read_stacked <- function(sheet, blocks, level_cols, level_name) {
  df <- read_raw(sheet)
  map_dfr(seq_along(blocks), function(b) {
    map_dfr(seq_along(strain_levels), function(s) {
      r_untr <- blocks[b] + 2 * (s - 1)
      map_dfr(names(level_cols), function(lv) {
        cols <- level_cols[[lv]]
        tibble(
          strain    = strain_levels[s],
          bio_rep   = b,
          !!level_name := as.numeric(lv),
          log10_untreated = tech_log(df, r_untr,     cols),
          log10_treated   = tech_log(df, r_untr + 1, cols)
        )
      })
    })
  })
}

# Block starts are the first untreated row of each biological replicate, and
# the column pairs are the two technical replicates for each level. Both are
# readxl's own indices, which differ from the raw worksheet because readxl
# trims leading empty rows and columns.
primary <- read_stacked("Antibacterial activity ", c(3, 21, 39),
                        list("100" = 4:5), "dose_ug_mL")

dose <- read_stacked("Dose killing", c(3, 20, 37),
                     list("0" = 4:5, "50" = 7:8, "100" = 11:12, "300" = 15:16),
                     "dose_ug_mL")

timekill <- read_stacked("Time killing", c(3, 20, 37),
                         list("0" = 4:5, "30" = 7:8, "60" = 11:12,
                              "120" = 15:16, "180" = 19:20),
                         "time_min")

# ---- Side-by-side sheets: EDTA and salt ------------------------------
# The three biological replicates sit in three column groups rather than three
# row blocks. Condition labels carry micro-sign and subscript-two characters
# whose encoding does not survive a literal string comparison, so rows are
# classified from ASCII features of the label instead.
read_sidebyside <- function(sheet, classify) {
  df <- read_raw(sheet)
  groups <- list(list(strain = 1, cond = 3, vals = 4:5),
                 list(strain = 7, cond = 9, vals = 10:11),
                 list(strain = 13, cond = 15, vals = 16:17))
  map_dfr(seq_along(groups), function(b) {
    g <- groups[[b]]
    sc <- norm(df[[g$strain]]); cc <- norm(df[[g$cond]])
    keep <- which(sc %in% strain_levels & !is.na(cc) & nzchar(cc) & cc != "Treatment")
    map_dfr(keep, function(r) {
      tibble(strain = sc[r], bio_rep = b,
             log10_cfu = tech_log(df, r, g$vals)) |>
        bind_cols(classify(cc[r]))
    })
  })
}

# EDTA labels are distinguished by which of PL25_M23, EDTA and ZnCl2 appear.
edta <- read_sidebyside("EDTA + ZnCl2", function(lab) {
  enzyme <- grepl("PL25", lab, fixed = TRUE)
  chel   <- grepl("EDTA", lab, fixed = TRUE)
  zinc   <- grepl("ZnCl", lab, fixed = TRUE)
  tibble(condition = if (!enzyme && !chel && !zinc) "Untreated"
                     else if (!enzyme &&  chel) "EDTA only"
                     else if (!enzyme &&  zinc) "Zn only"
                     else if ( chel &&  zinc)   "PL25_M23 + EDTA + Zn"
                     else if ( chel)            "PL25_M23 + EDTA"
                     else if ( zinc)            "PL25_M23 + Zn"
                     else                       "PL25_M23")
})

# Salt labels carry the molarity as a decimal number immediately before "M NaCl".
salt <- read_sidebyside("Salt-tolerance", function(lab) {
  tibble(treated = grepl("PL25", lab, fixed = TRUE),
         NaCl_M  = as.numeric(str_match(lab, "([0-9.]+)\\s*M NaCl")[, 2]))
})

# S. equi and S. pneumoniae were excluded from the salt assay: neither grows
# at 1.0 M NaCl, so their wells were never filled.
salt <- filter(salt, !is.na(log10_cfu))

write_tsv_ <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  %-20s %3d rows\n", name, nrow(x)))
}

cat(sprintf("reading %s\n\n", src))
write_tsv_(primary,  "primary_raw.tsv")
write_tsv_(dose,     "dose_raw.tsv")
write_tsv_(timekill, "timekill_raw.tsv")
write_tsv_(edta,     "edta_raw.tsv")
write_tsv_(salt,     "salt_raw.tsv")
cat(sprintf("\nwrote 5 files to %s\n", out))
