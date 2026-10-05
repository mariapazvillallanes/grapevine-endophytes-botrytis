# ============================================================
# 09 – qPCR analysis and expression figures
# Clean Zenodo reproducibility script
# ============================================================

local({
# ============================================================
# qPCR analysis – Supplementary Tables S17–S22
# Plant random intercept + Satterthwaite inference
# ============================================================

# Input:
#   data/09_qpcr_deltaCt.csv
#
# Primary model, fitted separately for each gene:
#   DeltaCt ~ Treatment * Time + (1 | Plant_ID)
#
# - REML estimation
# - Type III F tests with Satterthwaite denominator df
# - Directed contrasts: Satterthwaite t tests and BH adjustment
#   across 9 predefined contrasts within each gene and time
# - Interaction contrasts: Satterthwaite t tests and BH adjustment
#   across 10 tests within each gene
# - Missing DeltaCt values remain missing; no imputation

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "lme4",
  "lmerTest",
  "emmeans",
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

options(
  contrasts = c(
    "contr.sum",
    "contr.poly"
  )
)

emmeans::emm_options(
  lmer.df = "satterthwaite"
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
  "09_qpcr_deltaCt.csv"
)

output_dir <- file.path(project_dir, "output", "09_qpcr_analysis")

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
# Read primary DeltaCt data
# -------------------------

raw <- readr::read_delim(
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
  dplyr::rename(
    TreatmentCode = `Treatment code`,
    Time = `Time (dpi)`,
    DeltaCt = `ΔCt`
  ) %>%
  dplyr::mutate(
    TreatmentCode = factor(
      TreatmentCode,
      levels = c(
        "NC",
        "NC_BC",
        "1013",
        "1013_BC",
        "1017",
        "1017_BC"
      )
    ),
    Time = factor(
      Time,
      levels = c(
        3,
        5,
        7,
        10,
        14
      )
    ),
    Gene = factor(Gene),
    Plant_ID = interaction(
      TreatmentCode,
      Plant,
      drop = TRUE
    )
  )

analysis_data <- raw %>%
  dplyr::filter(
    !is.na(DeltaCt)
  )

sig_code <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# -------------------------
# Table S17
# Conventional relative expression
# -------------------------

# Table S17 in the final Supplementary Material was derived from
# the treatment-level mean DeltaCt values rounded to three decimals,
# matching the values reported in Table S18.
#
# Therefore:
#   log2FC = rounded mean DeltaCt of NC - rounded mean DeltaCt of treatment
#
# The SD reported in Table S17 is the SD of the original DeltaCt values
# within each treatment x gene x time cell (the sign change does not alter SD).

delta_summary_s17 <- analysis_data %>%
  dplyr::group_by(
    Gene,
    TreatmentCode,
    Treatment,
    Time
  ) %>%
  dplyr::summarise(
    n =
      dplyr::n(),
    Mean_DeltaCt =
      mean(DeltaCt),
    SD =
      stats::sd(DeltaCt),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Mean_DeltaCt_rounded =
      round(
        Mean_DeltaCt,
        3
      )
  )

nc_means_s17 <- delta_summary_s17 %>%
  dplyr::filter(
    TreatmentCode == "NC"
  ) %>%
  dplyr::select(
    Gene,
    Time,
    NC_mean_DeltaCt_rounded =
      Mean_DeltaCt_rounded
  )

table_s17 <- delta_summary_s17 %>%
  dplyr::filter(
    TreatmentCode != "NC"
  ) %>%
  dplyr::left_join(
    nc_means_s17,
    by = c(
      "Gene",
      "Time"
    )
  ) %>%
  dplyr::mutate(
    `log2FC (-DeltaDeltaCt)` =
      NC_mean_DeltaCt_rounded -
      Mean_DeltaCt_rounded,
    Time =
      as.integer(
        as.character(Time)
      )
  ) %>%
  dplyr::select(
    Gene,
    Treatment,
    Time,
    n,
    `log2FC (-DeltaDeltaCt)`,
    SD
  )

# -------------------------
# Table S18
# Mean DeltaCt ± SD
# -------------------------

table_s18 <- analysis_data %>%
  dplyr::group_by(
    Gene,
    Treatment,
    Time
  ) %>%
  dplyr::summarise(
    n =
      dplyr::n(),
    Mean_DeltaCt =
      mean(DeltaCt),
    SD_DeltaCt =
      stats::sd(DeltaCt),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Time =
      paste0(
        as.character(Time),
        " dpi"
      ),
    `Mean DeltaCt ± SD` =
      sprintf(
        "%.3f ± %.3f",
        Mean_DeltaCt,
        SD_DeltaCt
      )
  ) %>%
  dplyr::select(
    Gene,
    Treatment,
    Time,
    n,
    `Mean DeltaCt ± SD`
  )

# -------------------------
# Model contrasts
# -------------------------

directed_methods <- list(
  "NC+BC - NC" =
    c(1, -1, 0, 0, 0, 0),
  "1013 - NC" =
    c(1, 0, -1, 0, 0, 0),
  "1013+BC - NC" =
    c(1, 0, 0, -1, 0, 0),
  "1017 - NC" =
    c(1, 0, 0, 0, -1, 0),
  "1017+BC - NC" =
    c(1, 0, 0, 0, 0, -1),
  "1013+BC - NC+BC" =
    c(0, 1, 0, -1, 0, 0),
  "1017+BC - NC+BC" =
    c(0, 1, 0, 0, 0, -1),
  "1013+BC - 1013" =
    c(0, 0, 1, -1, 0, 0),
  "1017+BC - 1017" =
    c(0, 0, 0, 0, 1, -1)
)

interaction_methods <- list(
  "Ps. putida 1013 × B. cinerea" =
    c(-1, 1, 1, -1, 0, 0),
  "Pa. agglomerans 1017 × B. cinerea" =
    c(-1, 1, 0, 0, 1, -1)
)

genes <- levels(
  droplevels(
    analysis_data$Gene
  )
)

anova_results <- list()
directed_results <- list()
interaction_results <- list()
diagnostics <- list()

for (gene_i in genes) {

  dat_g <- analysis_data %>%
    dplyr::filter(
      Gene == gene_i
    ) %>%
    droplevels()

  model_g <- lmerTest::lmer(
    DeltaCt ~
      TreatmentCode *
      Time +
      (1 | Plant_ID),
    data = dat_g,
    REML = TRUE,
    control =
      lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl =
          list(
            maxfun = 200000
          )
      )
  )

  a <- stats::anova(
    model_g,
    type = 3,
    ddf = "Satterthwaite"
  )

  anova_results[[gene_i]] <- tibble::tibble(
    Gene = gene_i,
    Effect = rownames(a),
    `Numerator df` =
      a$NumDF,
    `Denominator df` =
      a$DenDF,
    `F value` =
      a$`F value`,
    p_numeric =
      a$`Pr(>F)`
  ) %>%
    dplyr::filter(
      Effect != "(Intercept)"
    ) %>%
    dplyr::mutate(
      Effect = dplyr::recode(
        Effect,
        "TreatmentCode" =
          "Treatment",
        "Time" =
          "Time",
        "TreatmentCode:Time" =
          "Treatment × Time"
      ),
      `p-value` =
        dplyr::if_else(
          p_numeric < 0.001,
          "<0.001",
          sprintf(
            "%.4f",
            p_numeric
          )
        ),
      Significance =
        sig_code(
          p_numeric
        )
    ) %>%
    dplyr::select(
      Gene,
      Effect,
      `Numerator df`,
      `Denominator df`,
      `F value`,
      `p-value`,
      Significance
    )

  emm_g <- emmeans::emmeans(
    model_g,
    ~ TreatmentCode | Time,
    weights = "equal",
    lmer.df = "satterthwaite"
  )

  dc <- emmeans::contrast(
    emm_g,
    method = directed_methods,
    adjust = "none"
  ) %>%
    summary(
      infer = c(
        TRUE,
        TRUE
      )
    ) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      Gene = gene_i,
      Contrast = contrast,
      Time =
        as.character(Time),
      log2FC = estimate,
      SE = SE,
      df = df,
      t_value = t.ratio,
      Raw_p = p.value
    ) %>%
    dplyr::group_by(
      Gene,
      Time
    ) %>%
    dplyr::mutate(
      BH_adjusted_p =
        stats::p.adjust(
          Raw_p,
          method = "BH"
        ),
      Significance =
        sig_code(
          BH_adjusted_p
        )
    ) %>%
    dplyr::ungroup()

  directed_results[[gene_i]] <- dc

  ic <- emmeans::contrast(
    emm_g,
    method = interaction_methods,
    adjust = "none"
  ) %>%
    summary(
      infer = c(
        TRUE,
        TRUE
      )
    ) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      Gene = gene_i,
      Interaction_tested = contrast,
      Time =
        as.character(Time),
      Interaction_log2FC =
        estimate,
      SE = SE,
      df = df,
      t_value =
        t.ratio,
      Raw_p =
        p.value
    ) %>%
    dplyr::mutate(
      BH_adjusted_p =
        stats::p.adjust(
          Raw_p,
          method = "BH"
        ),
      Significance =
        sig_code(
          BH_adjusted_p
        )
    )

  interaction_results[[gene_i]] <- ic

  vc <- as.data.frame(
    lme4::VarCorr(
      model_g
    )
  )

  plant_var <- vc %>%
    dplyr::filter(
      grp == "Plant_ID"
    ) %>%
    dplyr::pull(vcov)

  residual_var <-
    stats::sigma(
      model_g
    )^2

  diagnostics[[gene_i]] <- tibble::tibble(
    Gene = gene_i,
    Observed_DeltaCt =
      nrow(dat_g),
    Plants =
      dplyr::n_distinct(
        dat_g$Plant_ID
      ),
    Residual_variance =
      residual_var,
    Plant_random_intercept_variance =
      plant_var,
    ICC =
      plant_var /
      (
        plant_var +
        residual_var
      ),
    Boundary_fit =
      dplyr::if_else(
        lme4::isSingular(
          model_g,
          tol = 1e-5
        ),
        "YES",
        "NO"
      ),
    Converged =
      dplyr::if_else(
        is.null(
          model_g@optinfo$conv$lme4$messages
        ),
        "YES",
        "NO"
      )
  )
}

typeIII <- dplyr::bind_rows(
  anova_results
)

directed_detailed <- dplyr::bind_rows(
  directed_results
)

interaction_detailed <- dplyr::bind_rows(
  interaction_results
)

model_diagnostics <- dplyr::bind_rows(
  diagnostics
)

# -------------------------
# Table S20
# -------------------------

table_s20 <- directed_detailed %>%
  dplyr::mutate(
    Display =
      sprintf(
        "%.3f%s",
        log2FC,
        Significance
      ),
    Time =
      paste0(
        Time,
        " dpi"
      )
  ) %>%
  dplyr::select(
    Gene,
    Contrast,
    Time,
    Display
  ) %>%
  tidyr::pivot_wider(
    names_from = Time,
    values_from = Display
  )

# -------------------------
# Tables S21 and S22
# -------------------------

table_s21 <- interaction_detailed %>%
  dplyr::transmute(
    Gene,
    `Interaction tested` =
      Interaction_tested,
    `Time (dpi)` =
      as.integer(Time),
    `Interaction log2FC` =
      Interaction_log2FC,
    SE,
    df,
    `t value` =
      t_value,
    `Raw p-value` =
      Raw_p,
    `Adjusted p-value` =
      BH_adjusted_p,
    Significance
  )

table_s22 <- interaction_detailed %>%
  dplyr::mutate(
    Display =
      sprintf(
        "%.3f%s",
        Interaction_log2FC,
        Significance
      ),
    Time =
      paste0(
        Time,
        " dpi"
      )
  ) %>%
  dplyr::select(
    Gene,
    Interaction_tested,
    Time,
    Display
  ) %>%
  tidyr::pivot_wider(
    names_from = Time,
    values_from = Display
  ) %>%
  dplyr::rename(
    `Interaction tested` =
      Interaction_tested
  )

# -------------------------
# Figure 5 data
# -------------------------

figure5_data <- interaction_detailed %>%
  dplyr::mutate(
    Strain =
      dplyr::if_else(
        grepl(
          "1013",
          Interaction_tested
        ),
        "1013",
        "1017"
      ),
    Time_dpi =
      as.integer(Time)
  ) %>%
  dplyr::select(
    Strain,
    Gene,
    Time_dpi,
    Interaction_log2FC,
    Significance,
    Adjusted_p =
      BH_adjusted_p
  )

# -------------------------
# Export
# -------------------------

readr::write_csv(
  table_s17,
  file.path(
    output_dir,
    "Table_S17_Gene_Expression.csv"
  )
)

readr::write_csv(
  table_s18,
  file.path(
    output_dir,
    "Table_S18_qPCR_Summary_Statistics.csv"
  )
)

readr::write_csv(
  typeIII,
  file.path(
    output_dir,
    "Table_S19_qPCR_TypeIII_Tests.csv"
  )
)

readr::write_csv(
  table_s20,
  file.path(
    output_dir,
    "Table_S20_qPCR_Directed_Contrasts.csv"
  )
)

readr::write_csv(
  table_s21,
  file.path(
    output_dir,
    "Table_S21_Interaction_Values.csv"
  )
)

readr::write_csv(
  table_s22,
  file.path(
    output_dir,
    "Table_S22_qPCR_Interaction_Effect.csv"
  )
)

readr::write_csv(
  figure5_data,
  file.path(
    output_dir,
    "Figure_5_data_from_qPCR_analysis.csv"
  )
)

readr::write_csv(
  model_diagnostics,
  file.path(
    output_dir,
    "qPCR_model_diagnostics.csv"
  )
)


cat("\nqPCR analysis completed successfully.\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Generated Tables S17-S22 and Figure 5 source data.\n")



# -------------------------
# Figure S7 source matrix
# -------------------------
figure_s7_source <- directed_detailed %>%
  dplyr::filter(
    Contrast %in% c(
      "NC+BC - NC",
      "1013 - NC",
      "1013+BC - NC",
      "1017 - NC",
      "1017+BC - NC"
    )
  ) %>%
  dplyr::mutate(
    Treatment_label = dplyr::recode(
      Contrast,
      "NC+BC - NC" = "Control+B.cinerea",
      "1013 - NC" = "1013",
      "1013+BC - NC" = "1013+B.cinerea",
      "1017 - NC" = "1017",
      "1017+BC - NC" = "1017+B.cinerea"
    ),
    Column = paste0(Time, "_", Treatment_label)
  ) %>%
  dplyr::select(Gene, Column, log2FC) %>%
  tidyr::pivot_wider(names_from = Column, values_from = log2FC)

column_order <- as.vector(
  unlist(
    lapply(
      c(3, 5, 7, 10, 14),
      function(tt) paste0(
        tt, "_",
        c(
          "Control+B.cinerea",
          "1013",
          "1013+B.cinerea",
          "1017",
          "1017+B.cinerea"
        )
      )
    )
  )
)

figure_s7_source <- figure_s7_source %>%
  dplyr::select(Gene, dplyr::all_of(column_order))

readr::write_delim(
  figure_s7_source,
  file.path(output_dir, "Figure_S7_expression_data.csv"),
  delim = ";",
  na = "",
  quote = "needed"
)

})

# ---- Figure 5 ----
local({
# ==========================================================
# Figure 5 - Endophyte-dependent interaction heatmap
# Reads the FINAL Satterthwaite reanalysis directly from Excel
# ==========================================================

required_packages <- c(
  "readr", "dplyr", "ggplot2", "showtext",
  "sysfonts", "scales", "ragg"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

library(readr)
library(dplyr)
library(ggplot2)
library(showtext)
library(sysfonts)
library(scales)
library(ragg)

# ----------------------------------------------------------
# 1. Input and output paths
# ----------------------------------------------------------

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

input_file <- file.path(project_dir, "output", "09_qpcr_analysis", "Figure_5_data_from_qPCR_analysis.csv")

output_dir <- file.path(project_dir, "output", "09_qpcr_analysis")

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


# ----------------------------------------------------------
# 2. Read and prepare Figure 5 data
# ----------------------------------------------------------

df <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(decimal_mark = ","),
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
) %>%
  transmute(
    Strain = as.character(Strain),
    Gene = as.character(Gene),
    Time_dpi = as.integer(Time_dpi),
    log2FC_interaction = as.numeric(Interaction_log2FC),
    Significance = as.character(Significance),
    Adjusted_p = as.numeric(Adjusted_p)
  ) %>%
  mutate(
    # Do not print "ns" inside non-significant cells.
    Significance = if_else(
      is.na(Significance) | Significance == "ns",
      "",
      Significance
    ),

    # White text is used only on strongly coloured cells.
    TextColour = if_else(
      abs(log2FC_interaction) >= 4.5,
      "white",
      "black"
    ),

    Strain = factor(
      Strain,
      levels = c("1013", "1017")
    ),

    Time = factor(
      paste0(Time_dpi, " dpi"),
      levels = c(
        "3 dpi", "5 dpi", "7 dpi",
        "10 dpi", "14 dpi"
      )
    ),

    Gene = factor(
      Gene,
      levels = c(
        "VvJAZ1", "VvPAL1",
        "VvRBOHD", "VvSTS1.2"
      )
    )
  )

expected_rows <- 2 * 4 * 5

if (nrow(df) != expected_rows) {
  warning(
    "Expected ", expected_rows,
    " rows (2 strains × 4 genes × 5 times), but found ",
    nrow(df), "."
  )
}

if (
  any(is.na(df$Strain)) ||
  any(is.na(df$Gene)) ||
  any(is.na(df$Time)) ||
  any(is.na(df$log2FC_interaction))
) {
  stop(
    "Figure_5.csv contains missing or unrecognized values ",
    "in Strain, Gene, Time_dpi, or Interaction_log2FC."
  )
}


# ----------------------------------------------------------
# 3. Font
# ----------------------------------------------------------

if ("Times New Roman" %in% sysfonts::font_families()) {
  font_add(
    family = "TNR",
    regular = "Times New Roman"
  )
} else {
  # Same fallback used in the original script.
  font_add_google(
    "Tinos",
    "TNR"
  )
}

showtext_auto()

# ----------------------------------------------------------
# 4. Heatmap
# ----------------------------------------------------------

p_heat <- ggplot(
  df,
  aes(
    x = Time,
    y = Gene,
    fill = log2FC_interaction
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = Significance,
      color = TextColour
    ),
    fontface = "bold",
    size = 5.2,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~Strain,
    nrow = 1
  ) +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-7, 7),
    breaks = c(-7, -3.5, 0, 3.5, 7),
    oob = scales::squish,
    na.value = "white",
    name = expression(
      "Interaction " * log[2] * "FC"
    )
  ) +
  scale_color_identity() +
  scale_y_discrete(
    labels = function(x) {
      parse(
        text = paste0(
          "italic('", x, "')"
        )
      )
    }
  ) +
  labs(
    title = expression(
      "Endophyte-dependent modulation of the grapevine response to " *
        italic("Botrytis cinerea")
    ),
    x = "Time (days post-inoculation, dpi)",
    y = "Gene"
  ) +
  theme_minimal(
    base_family = "TNR",
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

# ----------------------------------------------------------
# 5. Save outputs in the same folder as the Excel workbook
# ----------------------------------------------------------

output_pdf <- file.path(
  output_dir,
  "Figure_5.pdf"
)

output_png <- file.path(
  output_dir,
  "Figure_5.png"
)

ggsave(
  filename = output_pdf,
  plot = p_heat,
  width = 8.5,
  height = 3.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ragg::agg_png(
  filename = output_png,
  width = 180,
  height = 75,
  units = "mm",
  res = 600,
  background = "white"
)

print(p_heat)
dev.off()

cat("\nFigure 5 generated successfully.\n")
cat("Files saved in:\n", output_dir, "\n\n")
cat("- ", basename(output_pdf), "\n", sep = "")
cat("- ", basename(output_png), "\n", sep = "")

})

# ---- Supplementary Figure S7 ----
local({
# ============================================================
# SUPPLEMENTARY FIGURE S7 – FINAL PUBLICATION VERSION
#
# - Datos: FigureS7_expression_data_STYLE_MATCH.csv
# - Sin título dentro de la figura
# - Fuente Times New Roman
# - Nombres de los genes en cursiva
# - Se conserva el diseño y la escala de colores original
# - PNG y TIFF a 600 dpi
# - PDF de alta resolución
# ============================================================


# ------------------------------------------------------------
# 1. INSTALAR Y CARGAR PAQUETES
# ------------------------------------------------------------

paquetes <- c(
  "pheatmap",
  "RColorBrewer",
  "ragg",
  "png",
  "systemfonts"
)

missing_packages <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

library(pheatmap)
library(RColorBrewer)
library(grid)
library(ragg)


# ------------------------------------------------------------
# 2. CERRAR DISPOSITIVOS GRÁFICOS ANTERIORES
# ------------------------------------------------------------

graphics.off()

if (requireNamespace("showtext", quietly = TRUE)) {
  showtext::showtext_auto(enable = FALSE)
}


# ------------------------------------------------------------
# 3. INPUT AND OUTPUT PATHS
# ------------------------------------------------------------

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

input_file <- file.path(project_dir, "output", "09_qpcr_analysis", "Figure_S7_expression_data.csv")

output_dir <- file.path(project_dir, "output", "09_qpcr_analysis")

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


# ------------------------------------------------------------
# 4. COMPROBAR TIMES NEW ROMAN
# ------------------------------------------------------------

fuentes_disponibles <- systemfonts::system_fonts()

if (!any(
  tolower(fuentes_disponibles$family) ==
  "times new roman"
)) {
  stop(
    "No se ha encontrado Times New Roman entre las fuentes instaladas."
  )
}

font_family <- "Times New Roman"


# ------------------------------------------------------------
# 5. LEER LOS DATOS
# ------------------------------------------------------------

# El archivo utiliza punto y coma como separador
# y coma como separador decimal.

df <- read.csv2(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

if (!"Gene" %in% names(df)) {
  stop("El archivo no contiene la columna Gene.")
}

rownames(df) <- df$Gene

df$Gene <- NULL


# ------------------------------------------------------------
# 6. CONSTRUIR LA MATRIZ
# ------------------------------------------------------------

# Convertir los valores a formato numérico.
# No se recalculan los valores de expresión.

to_num <- function(x) {
  
  as.numeric(
    gsub(",", ".", as.character(x))
  )
  
}

mat <- as.matrix(
  data.frame(
    lapply(df, to_num),
    check.names = FALSE
  )
)

rownames(mat) <- rownames(df)


if (any(!is.finite(mat))) {
  
  stop(
    "El CSV contiene valores ausentes o no numéricos. Revisa el archivo seleccionado."
  )
  
}


# ------------------------------------------------------------
# 7. IDENTIFICAR TIEMPOS Y TRATAMIENTOS
# ------------------------------------------------------------

cn <- colnames(mat)

Time_num <- as.integer(
  sub("^([0-9]+)_.*", "\\1", cn)
)

Time_dpi <- paste0(
  Time_num,
  " dpi"
)


Treatment <- sub(
  "^[0-9]+_",
  "",
  cn
)

Treatment <- gsub(
  "Control\\+B\\.cinerea",
  "NC+Bc",
  Treatment
)

Treatment <- gsub(
  "B\\.cinerea",
  "Bc",
  Treatment
)


# Mantener el orden original de los tratamientos.

treat_levels <- c(
  "NC+Bc",
  "1013",
  "1013+Bc",
  "1017",
  "1017+Bc"
)

Treatment <- factor(
  Treatment,
  levels = treat_levels
)


# Orden cronológico.

time_levels <- c(
  "3 dpi",
  "5 dpi",
  "7 dpi",
  "10 dpi",
  "14 dpi"
)

Time_dpi <- factor(
  Time_dpi,
  levels = time_levels
)


if (
  any(is.na(Treatment)) ||
  any(is.na(Time_dpi))
) {
  
  stop(
    "Hay tiempos o tratamientos que no coinciden con el diseño esperado."
  )
  
}


# ------------------------------------------------------------
# 8. ANOTACIONES SUPERIORES
# ------------------------------------------------------------

ann_col <- data.frame(
  
  Time = Time_dpi,
  
  Treatment = Treatment,
  
  row.names = cn,
  
  check.names = FALSE
  
)


# Ordenar las columnas por tiempo y tratamiento.

ord <- order(
  Time_num,
  ann_col$Treatment
)

mat <- mat[
  ,
  ord,
  drop = FALSE
]

ann_col <- ann_col[
  ord,
  ,
  drop = FALSE
]


# Comprobar que hay 25 columnas:
# 5 tiempos x 5 tratamientos.

if (ncol(mat) != 25L) {
  
  stop(
    "Se esperaban 25 combinaciones de tiempo y tratamiento."
  )
  
}


# ------------------------------------------------------------
# 9. ORDEN DE LOS GENES
# ------------------------------------------------------------

gene_order <- c(
  "VvJAZ1",
  "VvRBOHD",
  "VvPAL1",
  "VvSTS1.2"
)


if (!setequal(rownames(mat), gene_order)) {
  
  stop(
    "Los genes del CSV no coinciden con los cuatro genes esperados."
  )
  
}


mat <- mat[
  gene_order,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 10. ESCALA DE EXPRESIÓN Y COLORES
# ------------------------------------------------------------

# Se mantiene la escala visual original de -2 a 2.
#
# Los valores superiores a 2 o inferiores a -2
# se limitan únicamente para asignarles un color.
#
# No se modifican los datos originales del CSV.

cap <- 2

mat_cap <- mat

mat_cap[mat_cap > cap] <- cap

mat_cap[mat_cap < -cap] <- -cap


breaks <- seq(
  -cap,
  cap,
  length.out = 101
)


pal <- colorRampPalette(
  rev(
    RColorBrewer::brewer.pal(
      11,
      "RdBu"
    )
  )
)(100)


# ------------------------------------------------------------
# 11. COLORES DE LAS ANOTACIONES
# ------------------------------------------------------------

time_palette <- c(
  
  "3 dpi" = "#DEEBF7",
  
  "5 dpi" = "#C6DBEF",
  
  "7 dpi" = "#9ECAE1",
  
  "10 dpi" = "#6BAED6",
  
  "14 dpi" = "#3182BD"
  
)


ann_colors <- list(
  
  Time = time_palette,
  
  Treatment = c(
    
    "NC+Bc" = "#BDBDBD",
    
    "1013" = "#1B9E77",
    
    "1013+Bc" = "#66C2A5",
    
    "1017" = "#7570B3",
    
    "1017+Bc" = "#BC80BD"
    
  )
  
)


# ------------------------------------------------------------
# 12. PREPARAR LOS NOMBRES DE LOS ARCHIVOS
# ------------------------------------------------------------

archivo_png <- file.path(
  output_dir,
  "Figure_S7.png"
)

archivo_tiff <- file.path(
  output_dir,
  "Figure_S7.tiff"
)

archivo_pdf <- file.path(
  output_dir,
  "Figure_S7.pdf"
)


# ------------------------------------------------------------
# 13. ABRIR EL PNG ANTES DE CREAR EL GRÁFICO
# ------------------------------------------------------------

# De esta forma, R calcula el tamaño de las etiquetas
# utilizando un dispositivo compatible con Times New Roman.

ragg::agg_png(
  
  filename = archivo_png,
  
  width = 12,
  
  height = 4.8,
  
  units = "in",
  
  res = 600
  
)


# ------------------------------------------------------------
# 14. CONSTRUIR EL MAPA DE CALOR
# ------------------------------------------------------------

# silent = TRUE permite recuperar el gráfico y modificar
# las etiquetas de los genes antes de dibujarlo.

p <- pheatmap::pheatmap(
  
  mat_cap,
  
  color = pal,
  
  breaks = breaks,
  
  border_color = NA,
  
  na_col = "grey90",
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  gaps_col = c(5, 10, 15, 20),
  
  annotation_col = ann_col,
  
  annotation_colors = ann_colors,
  
  show_colnames = FALSE,
  
  labels_row = gene_order,
  
  fontsize = 10,
  
  fontsize_row = 11,
  
  fontsize_col = 8,
  
  fontfamily = font_family,
  
  # NO se incluye main: así desaparece el título.
  
  silent = TRUE
  
)


# ------------------------------------------------------------
# 15. PONER LOS GENES EN CURSIVA
# ------------------------------------------------------------

# Buscamos el elemento gráfico que contiene los nombres
# de los genes y aplicamos cursiva únicamente a esos nombres.

indice_genes <- which(
  p$gtable$layout$name == "row_names"
)


if (length(indice_genes) != 1L) {
  
  stop(
    "No se ha encontrado el elemento gráfico de los nombres de los genes."
  )
  
}


p$gtable$grobs[[indice_genes]]$gp <-
  grid::gpar(
    
    fontsize = 11,
    
    fontfamily = font_family,
    
    fontface = "italic",
    
    col = "black"
    
  )


# ------------------------------------------------------------
# 16. DIBUJAR Y GUARDAR EL PNG
# ------------------------------------------------------------

grid::grid.newpage()

grid::grid.draw(
  p$gtable
)

grDevices::dev.off()


# ------------------------------------------------------------
# 17. GUARDAR TAMBIÉN EN TIFF A 600 DPI
# ------------------------------------------------------------

ragg::agg_tiff(
  
  filename = archivo_tiff,
  
  width = 12,
  
  height = 4.8,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw"
  
)


grid::grid.newpage()

grid::grid.draw(
  p$gtable
)

grDevices::dev.off()


# ------------------------------------------------------------
# 18. GUARDAR TAMBIÉN EN PDF
# ------------------------------------------------------------

# Para evitar el error PDF CID font, incorporamos al PDF
# el PNG que acabamos de generar a 600 dpi.
#
# El PDF mantiene exactamente el mismo aspecto que el PNG
# y el TIFF, aunque su contenido es una imagen rasterizada.

imagen_png <- png::readPNG(
  archivo_png
)


grDevices::cairo_pdf(
  
  filename = archivo_pdf,
  
  width = 12,
  
  height = 4.8,
  
  family = "sans"
  
)


grid::grid.newpage()

grid::grid.raster(
  
  imagen_png,
  
  x = 0.5,
  
  y = 0.5,
  
  width = grid::unit(
    1,
    "npc"
  ),
  
  height = grid::unit(
    1,
    "npc"
  ),
  
  interpolate = FALSE
  
)

grDevices::dev.off()


# ------------------------------------------------------------
# 19. CONFIRMAR LOS ARCHIVOS GENERADOS
# ------------------------------------------------------------

cat(
  
  "\nFIGURA S7 GENERADA CORRECTAMENTE\n\n",
  
  "PNG:\n",
  archivo_png,
  "\n\n",
  
  "TIFF:\n",
  archivo_tiff,
  "\n\n",
  
  "PDF:\n",
  archivo_pdf,
  "\n\n"
  
)
})
