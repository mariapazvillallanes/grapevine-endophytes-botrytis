# ============================================================
# 07 – Disease severity
# Clean Zenodo reproducibility script
# ============================================================

local({
# Figure 4 — Disease severity in grapevine leaves
# Recreates the supplied Figure_4.pdf from the cleaned raw data.
#
# Input:
#   07_disease_severity.csv
#
# Output:
#   Figure_4.pdf
#   Figure_4.png
#
# Notes:
# - The source workbook supplied for this reconstruction contains NC, NC+BC,
#   1013+BC and 1017+BC, but not Switch®+BC.
# - To reproduce the supplied Figure 4 exactly, 07_disease_severity.csv includes the
#   Switch®+BC observations shown in that figure (0% necrotic area at 7, 10,
#   and 14 dpi; n = 3 at each time).
# - Statistical grouping letters are carried over from the supplied final
#   figure rather than recalculated here.

# ----------------------------- Packages ---------------------------------------
packages <- c("ggplot2", "dplyr", "readr")
missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) stop("Missing required R packages: ", paste(missing, collapse = ", "))
library(ggplot2)
library(dplyr)
library(readr)

# ------------------------- Input and output paths ------------------------------

# The script can be run either from the project root (ZENODO)
# or from the /scripts directory.

if (dir.exists("data")) {
  project_dir <- "."
} else if (dir.exists("../data")) {
  project_dir <- ".."
} else {
  stop(
    "Project directory not found. ",
    "Run the script from the ZENODO project root or from the /scripts directory."
  )
}

input_file <- file.path(
  project_dir,
  "data",
  "07_disease_severity.csv"
)

output_dir <- file.path(project_dir, "output", "07_disease_severity")

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file
  )
}

# ------------------------------- Data -----------------------------------------

dat <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(decimal_mark = ","),
  show_col_types = FALSE
) %>%
  dplyr::filter(Treatment %in% c("NC_BC", "Switch_BC", "Pp1013_BC", "Pa1017_BC")) %>%
  dplyr::mutate(
    Treatment = dplyr::recode(
      Treatment,
      "NC_BC" = "NC+BC",
      "Switch_BC" = "Switch®+BC",
      "Pp1013_BC" = "1013+BC",
      "Pa1017_BC" = "1017+BC"
    )
  )

treatment_levels <- c("NC+BC", "1017+BC", "1013+BC", "Switch®+BC")

dat <- dat %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    Time_dpi = as.numeric(Time_dpi)
  )

# Deterministic horizontal offsets for the three biological replicates.
# This avoids random jitter and makes the exported figure reproducible.
rep_offsets <- c(`1` = -0.085, `2` = 0, `3` = 0.085)

dat <- dat %>%
  mutate(x_rep = Time_dpi + unname(rep_offsets[as.character(Replicate)]))

# ------------------------- Summary statistics ---------------------------------
sum_dat <- dat %>%
  group_by(Treatment, Time_dpi) %>%
  summarise(
    n = n(),
    mean = mean(Necrotic_pct, na.rm = TRUE),
    sd = sd(Necrotic_pct, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    sd = ifelse(is.na(sd), 0, sd),
    se = ifelse(is.na(se), 0, se)
  )

# ---------------------- Statistical grouping letters --------------------------
# Grouping displayed in the supplied Figure 4:
# NC+BC = a; 1017+BC = ab; 1013+BC = b; Switch®+BC = b
letters_dat <- expand.grid(
  Time_dpi = c(7, 10, 14),
  Treatment = treatment_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    letter = case_when(
      Treatment == "NC+BC"      ~ "a",
      Treatment == "1017+BC"    ~ "ab",
      Treatment == "1013+BC"    ~ "b",
      Treatment == "Switch®+BC" ~ "b"
    )
  ) %>%
  left_join(sum_dat, by = c("Treatment", "Time_dpi")) %>%
  mutate(
    # Vertical positions tuned to match the supplied figure
    y_lab = case_when(
      Treatment == "NC+BC"   ~ mean + se + 0.95,
      Treatment == "1017+BC" ~ 1.95,
      Treatment == "1013+BC" ~ 1.55,
      Treatment == "Switch®+BC" ~ 1.55
    ),
    # Horizontal offsets only for the crowded low-value labels
    x_lab = Time_dpi + case_when(
      Treatment == "1017+BC"    ~ -0.31,
      Treatment == "1013+BC"    ~  0.00,
      Treatment == "Switch®+BC" ~  0.31,
      TRUE                       ~  0.00
    )
  )

# ----------------------------- Appearance -------------------------------------
# Colours recovered from the supplied final PDF.
cols <- c(
  "NC+BC"      = "#9E1B1B",
  "1017+BC"    = "#2B66AD",
  "1013+BC"    = "#1B7F39",
  "Switch®+BC" = "#4D4D4D"
)

# Times New Roman is used in the original figure.
# On Windows this family name normally works directly.
font_family <- "Times New Roman"

# ------------------------------- Plot -----------------------------------------
p <- ggplot() +

  # Mean trajectories
  geom_line(
    data = sum_dat,
    aes(x = Time_dpi, y = mean, colour = Treatment, group = Treatment),
    linewidth = 0.80,
    lineend = "round"
  ) +

  # Mean ± SE
  geom_errorbar(
    data = sum_dat,
    aes(
      x = Time_dpi,
      ymin = pmax(mean - se, 0),
      ymax = mean + se,
      colour = Treatment
    ),
    width = 0.16,
    linewidth = 0.70
  ) +

  # Individual biological replicates: open circles
  geom_point(
    data = dat,
    aes(x = x_rep, y = Necrotic_pct, colour = Treatment),
    shape = 21,
    fill = "white",
    size = 2.55,
    stroke = 0.75,
    show.legend = FALSE
  ) +

  # Treatment means: filled circles
  geom_point(
    data = sum_dat,
    aes(x = Time_dpi, y = mean, colour = Treatment),
    size = 3.35
  ) +

  # Grouping letters
  geom_text(
    data = letters_dat,
    aes(x = x_lab, y = y_lab, label = letter, colour = Treatment),
    family = font_family,
    size = 4.1,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = cols,
    breaks = treatment_levels,
    labels = treatment_levels,
    name = "Treatment"
  ) +

  scale_x_continuous(
    breaks = c(7, 10, 14),
    limits = c(6.30, 14.70),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    breaks = seq(0, 50, 10),
    limits = c(-2.0, 51.5),
    expand = c(0, 0)
  ) +

  labs(
    x = "Sampling time (dpi)",
    y = "Necrotic leaf area (%)"
  ) +

  theme_classic(base_family = font_family, base_size = 12) +

  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),
    axis.title.y = element_text(size = 14, margin = margin(r = 8)),
    axis.text = element_text(size = 12, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.65),
    axis.ticks = element_line(colour = "black", linewidth = 0.55),
    axis.ticks.length = unit(3.5, "pt"),

    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.key.width = unit(1.35, "cm"),
    legend.key.height = unit(0.62, "cm"),

    plot.margin = margin(t = 8, r = 18, b = 8, l = 8)
  ) +

  guides(
    colour = guide_legend(
      override.aes = list(
        linewidth = 0.80,
        size = 3.35,
        shape = 16
      )
    )
  ) +

  coord_cartesian(clip = "off")

# ------------------------------- Export ---------------------------------------
# Dimensions match the supplied PDF aspect ratio (5:3).
ggsave(
  filename = file.path(
    output_dir,
    "Figure_4.pdf"
  ),
  plot = p,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_4.png"
  ),
  plot = p,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

print(p)

cat("\nFigure 4 analysis completed successfully.\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")

})

# ---- Supplementary Tables S14-S16 ----
local({
# ============================================================
# Disease severity analysis – Supplementary Tables S14–S16
# ============================================================

# Input:
#   data/07_disease_severity.csv
#
# Experimental structure:
#   7 treatments
#   3 sampling times: 7, 10 and 14 dpi
#   3 independent biological replicates per treatment × time
#
# Table S14:
#   Individual replicate values plus mean, SD and SE.
#
# Tables S15–S16:
#   Statistical analysis restricted to the four pathogen-inoculated
#   treatments:
#     NC+BC, Switch®+BC, 1013+BC, 1017+BC
#
#   At each sampling time:
#     Kruskal–Wallis test
#     Dunn post hoc comparisons with BH adjustment

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "rstatix",
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

invisible(
  lapply(
    required_packages,
    library,
    character.only = TRUE
  )
)

# -------------------------
# Input and output paths
# -------------------------

if (dir.exists("data")) {
  project_dir <- "."
} else if (dir.exists("../data")) {
  project_dir <- ".."
} else {
  stop(
    "Project directory not found. ",
    "Run the script from the ZENODO project root or from the /scripts directory."
  )
}

input_file <- file.path(
  project_dir,
  "data",
  "07_disease_severity.csv"
)

output_dir <- file.path(project_dir, "output", "07_disease_severity")

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file
  )
}

# -------------------------
# Read and validate data
# -------------------------

severity <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(
    decimal_mark = ","
  ),
  show_col_types = FALSE,
  na = c(
    "",
    "NA",
    "NaN"
  )
) %>%
  dplyr::mutate(
    Treatment =
      as.character(Treatment),
    Time_dpi =
      as.integer(Time_dpi),
    Replicate =
      as.integer(Replicate),
    Necrotic_pct =
      as.numeric(Necrotic_pct)
  )

required_columns <- c(
  "Treatment",
  "Time_dpi",
  "Replicate",
  "Necrotic_pct"
)

missing_columns <- setdiff(
  required_columns,
  names(severity)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

treatment_levels <- c(
  "NC",
  "NC_BC",
  "Switch_BC",
  "Pp1013",
  "Pa1017",
  "Pp1013_BC",
  "Pa1017_BC"
)

treatment_labels <- c(
  "NC" =
    "Negative control",
  "NC_BC" =
    "NC+BC",
  "Switch_BC" =
    "Switch®+BC",
  "Pp1013" =
    "1013",
  "Pa1017" =
    "1017",
  "Pp1013_BC" =
    "1013+BC",
  "Pa1017_BC" =
    "1017+BC"
)

severity <- severity %>%
  dplyr::mutate(
    Treatment = factor(
      Treatment,
      levels = treatment_levels
    )
  )

expected_n <- 7 * 3 * 3

if (nrow(severity) != expected_n) {
  warning(
    "Expected ",
    expected_n,
    " rows but found ",
    nrow(severity),
    "."
  )
}

replicate_check <- severity %>%
  dplyr::count(
    Treatment,
    Time_dpi
  )

if (any(replicate_check$n != 3L)) {
  stop(
    "Every treatment × time combination must contain exactly 3 biological replicates."
  )
}

# -------------------------
# Table S14
# -------------------------

table_s14 <- severity %>%
  dplyr::mutate(
    Treatment_label =
      unname(
        treatment_labels[
          as.character(Treatment)
        ]
      )
  ) %>%
  dplyr::select(
    Treatment_label,
    Time_dpi,
    Replicate,
    Necrotic_pct
  ) %>%
  tidyr::pivot_wider(
    names_from = Replicate,
    values_from = Necrotic_pct,
    names_prefix = "Replicate_"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Mean =
      mean(
        c_across(
          starts_with("Replicate_")
        )
      ),
    SD =
      stats::sd(
        c_across(
          starts_with("Replicate_")
        )
      ),
    SE =
      SD / sqrt(3),
    `Mean ± SD` =
      sprintf(
        "%.2f ± %.2f",
        Mean,
        SD
      ),
    `Mean ± SE` =
      sprintf(
        "%.2f ± %.2f",
        Mean,
        SE
      )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Treatment =
      factor(
        Treatment_label,
        levels =
          unname(
            treatment_labels[
              treatment_levels
            ]
          )
      )
  ) %>%
  dplyr::arrange(
    Treatment,
    Time_dpi
  ) %>%
  dplyr::transmute(
    Treatment =
      as.character(Treatment),
    `Time (dpi)` =
      Time_dpi,
    `Replicate 1 (%)` =
      Replicate_1,
    `Replicate 2 (%)` =
      Replicate_2,
    `Replicate 3 (%)` =
      Replicate_3,
    `Mean (%)` =
      Mean,
    SD =
      SD,
    SE =
      SE,
    `Mean ± SD`,
    `Mean ± SE`
  )

# -------------------------
# Pathogen-inoculated subset
# -------------------------

pathogen_levels <- c(
  "NC_BC",
  "Switch_BC",
  "Pp1013_BC",
  "Pa1017_BC"
)

stats_data <- severity %>%
  dplyr::filter(
    Treatment %in% pathogen_levels
  ) %>%
  dplyr::mutate(
    Treatment = factor(
      as.character(Treatment),
      levels = pathogen_levels
    )
  )

# -------------------------
# Table S15
# Kruskal-Wallis at each time
# -------------------------

kw_detailed <- stats_data %>%
  dplyr::group_by(
    Time_dpi
  ) %>%
  rstatix::kruskal_test(
    Necrotic_pct ~ Treatment
  ) %>%
  dplyr::ungroup()

table_s15 <- kw_detailed %>%
  dplyr::transmute(
    `Sampling time` =
      paste0(
        Time_dpi,
        " dpi"
      ),
    Variable =
      "Necrotic leaf area (%)",
    n =
      12L,
    `Kruskal–Wallis χ²` =
      round(
        statistic,
        3
      ),
    df =
      df,
    `p-value` =
      round(
        p,
        4
      ),
    Significance =
      dplyr::if_else(
        p < 0.05,
        "*",
        "ns"
      )
  )

# -------------------------
# Table S16
# Dunn + BH at each time
# -------------------------

dunn_detailed <- stats_data %>%
  dplyr::group_by(
    Time_dpi
  ) %>%
  rstatix::dunn_test(
    Necrotic_pct ~ Treatment,
    p.adjust.method = "BH"
  ) %>%
  dplyr::ungroup()

group_labels <- c(
  "NC_BC" =
    "NC+BC",
  "Switch_BC" =
    "Switch®+BC",
  "Pp1013_BC" =
    "1013+BC",
  "Pa1017_BC" =
    "1017+BC"
)

# rstatix reports the Dunn statistic using group1 - group2.
# The final Supplementary Table S16 presents the same listed comparisons
# using the opposite sign convention, so the displayed Z statistic is
# multiplied by -1. Raw and BH-adjusted p-values are unchanged.
table_s16 <- dunn_detailed %>%
  dplyr::mutate(
    group1_label =
      unname(
        group_labels[
          group1
        ]
      ),
    group2_label =
      unname(
        group_labels[
          group2
        ]
      ),
    Comparison =
      paste0(
        group1_label,
        " vs ",
        group2_label
      ),
    Significance =
      dplyr::if_else(
        p.adj < 0.05,
        "*",
        "ns"
      )
  ) %>%
  dplyr::transmute(
    `Sampling time` =
      paste0(
        Time_dpi,
        " dpi"
      ),
    Comparison,
    n1 =
      3L,
    n2 =
      3L,
    `Dunn Z statistic` =
      round(
        -statistic,
        3
      ),
    `Raw p-value` =
      round(
        p,
        3
      ),
    `BH-adjusted p-value` =
      round(
        p.adj,
        3
      ),
    Significance
  )

# -------------------------
# Export
# -------------------------

readr::write_csv(
  table_s14,
  file.path(
    output_dir,
    "Table_S14_Necrotic_Leaf_Area.csv"
  )
)

readr::write_csv(
  table_s15,
  file.path(
    output_dir,
    "Table_S15_Necrotic_KW.csv"
  )
)

readr::write_csv(
  table_s16,
  file.path(
    output_dir,
    "Table_S16_Necrotic_Dunn.csv"
  )
)

readr::write_csv(
  kw_detailed,
  file.path(
    output_dir,
    "Table_S15_Kruskal_Wallis_detailed.csv"
  )
)

readr::write_csv(
  dunn_detailed,
  file.path(
    output_dir,
    "Table_S16_Dunn_BH_detailed.csv"
  )
)


cat(
  "\nDisease severity analysis completed successfully.\n"
)
cat(
  "Input file:",
  input_file,
  "\n"
)
cat(
  "Generated Tables S14-S16 in:",
  output_dir,
  "\n"
)

})
