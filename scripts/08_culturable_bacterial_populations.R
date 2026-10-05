# ============================================================
# 08 – Culturable bacterial populations
# Clean Zenodo reproducibility script
# ============================================================

# ============================================================
# CFU analysis according to the manuscript
# Supplementary Tables S11–S13
# ============================================================

# Input:
#   data/08_culturable_bacterial_populations.csv
#
# Descriptive reporting:
#   log10 CFU g^-1 fresh weight, mean ± SD at T0–T4.
#
# Inferential analysis (as described in the manuscript):
#   1) For each plant and compartment, calculate within-plant
#      baseline-normalized change:
#
#         Delta_log10_CFU = log10_CFU_t - log10_CFU_T0
#
#   2) Exclude T0 from model fitting.
#   3) Fit separate linear mixed-effects models for Leaf and Soil:
#
#         Delta_log10_CFU ~ Treatment * Time + (1 | PlantID)
#
#      where PlantID = Treatment:Plant.
#
#   4) Evaluate Treatment, Time and Treatment × Time using
#      Type III Wald chi-square tests with sum-to-zero contrasts.
#
#   5) Compare treatments within each post-baseline time using
#      estimated marginal means with Tukey adjustment.
#
# This script intentionally follows the statistical method stated
# in the manuscript, even if older archived outputs differ.

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "lme4",
  "lmerTest",
  "emmeans",
  "car",
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
  "08_culturable_bacterial_populations.csv"
)

output_dir <- file.path(project_dir, "output", "08_culturable_bacterial_populations")

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

df <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(
    decimal_mark = ","
  ),
  show_col_types = FALSE,
  trim_ws = TRUE,
  na = c("", "NA", "NaN")
)

required_columns <- c(
  "Treatment",
  "Plant",
  "Time",
  "log10_CFU",
  "Compartment"
)

missing_columns <- setdiff(
  required_columns,
  names(df)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

treatment_levels <- c(
  "NC",
  "B.C",
  "1013",
  "1017",
  "1013_plus_B.C",
  "1017_plus_B.C"
)

time_levels_all <- c(
  "T0",
  "T1",
  "T2",
  "T3",
  "T4"
)

time_levels_model <- c(
  "T1",
  "T2",
  "T3",
  "T4"
)

df <- df %>%
  dplyr::mutate(
    Treatment = factor(
      as.character(Treatment),
      levels = treatment_levels
    ),
    Plant = factor(
      as.character(Plant)
    ),
    Time = factor(
      as.character(Time),
      levels = time_levels_all,
      ordered = TRUE
    ),
    Compartment = factor(
      as.character(Compartment),
      levels = c("Leaf", "Soil")
    ),
    log10_CFU = as.numeric(log10_CFU),
    PlantID = interaction(
      Treatment,
      Plant,
      drop = TRUE
    )
  )

if (anyNA(df[, required_columns])) {
  stop(
    "Missing or unrecognized values were detected in 08_culturable_bacterial_populations.csv."
  )
}

# Check expected balanced design:
# 2 compartments × 6 treatments × 5 times × 3 plants = 180
expected_n <- 2 * 6 * 5 * 3

if (nrow(df) != expected_n) {
  warning(
    "Expected ",
    expected_n,
    " observations, but found ",
    nrow(df),
    "."
  )
}

balance_check <- df %>%
  dplyr::count(
    Compartment,
    Treatment,
    Time,
    Plant,
    name = "n"
  )

if (any(balance_check$n != 1L)) {
  stop(
    "The experimental design is not balanced as expected."
  )
}

# -------------------------
# Table S11
# Descriptive log10 CFU mean ± SD
# -------------------------

summary_tbl <- df %>%
  dplyr::group_by(
    Compartment,
    Treatment,
    Time
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_log10 = mean(log10_CFU),
    sd_log10 = stats::sd(log10_CFU),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_sd = sprintf(
      "%.2f ± %.2f",
      mean_log10,
      sd_log10
    )
  )

treatment_labels <- c(
  "NC" = "Negative control",
  "B.C" = "B. cinerea",
  "1013" = "Ps. putida 1013",
  "1017" = "Pa. agglomerans 1017",
  "1013_plus_B.C" =
    "Ps. putida 1013 + B. cinerea",
  "1017_plus_B.C" =
    "Pa. agglomerans 1017 + B. cinerea"
)

table_s11 <- summary_tbl %>%
  dplyr::mutate(
    Treatment_label =
      unname(
        treatment_labels[
          as.character(Treatment)
        ]
      ),
    Time = as.character(Time)
  ) %>%
  dplyr::select(
    Compartment,
    Treatment = Treatment_label,
    Time,
    mean_sd
  ) %>%
  tidyr::pivot_wider(
    names_from = Time,
    values_from = mean_sd
  )

# -------------------------
# Baseline normalization
# -------------------------

baseline <- df %>%
  dplyr::filter(
    Time == "T0"
  ) %>%
  dplyr::select(
    Compartment,
    Treatment,
    Plant,
    PlantID,
    baseline_log10_CFU = log10_CFU
  )

analysis_df <- df %>%
  dplyr::filter(
    Time != "T0"
  ) %>%
  dplyr::left_join(
    baseline,
    by = c(
      "Compartment",
      "Treatment",
      "Plant",
      "PlantID"
    )
  ) %>%
  dplyr::mutate(
    Delta_log10_CFU =
      log10_CFU -
      baseline_log10_CFU,
    Time = factor(
      as.character(Time),
      levels = time_levels_model,
      ordered = TRUE
    )
  )

if (anyNA(analysis_df$baseline_log10_CFU)) {
  stop(
    "At least one post-baseline observation has no matching T0 baseline."
  )
}

# -------------------------
# Mixed model by compartment
# -------------------------

fit_compartment <- function(
  data,
  compartment_name
) {

  d <- data %>%
    dplyr::filter(
      Compartment == compartment_name
    ) %>%
    droplevels()

  model <- lme4::lmer(
    Delta_log10_CFU ~
      Treatment *
      Time +
      (1 | PlantID),
    data = d,
    REML = FALSE
  )

  type3 <- as.data.frame(
    car::Anova(
      model,
      type = 3
    )
  ) %>%
    tibble::rownames_to_column(
      "Effect"
    ) %>%
    dplyr::filter(
      Effect != "(Intercept)"
    ) %>%
    dplyr::transmute(
      Compartment = compartment_name,
      Effect = dplyr::recode(
        Effect,
        "Treatment" = "Treatment",
        "Time" = "Time",
        "Treatment:Time" =
          "Treatment × Time"
      ),
      Chi_square = Chisq,
      df = Df,
      p_value = `Pr(>Chisq)`
    )

  emm <- emmeans::emmeans(
    model,
    ~ Treatment | Time
  )

  # Full Tukey-adjusted pairwise comparisons
  pairwise_tukey <- emmeans::contrast(
    emm,
    method = "pairwise",
    adjust = "tukey"
  ) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      Compartment =
        compartment_name
    )

  list(
    model = model,
    type3 = type3,
    emmeans = emm,
    pairwise_tukey = pairwise_tukey
  )
}

leaf_res <- fit_compartment(
  analysis_df,
  "Leaf"
)

soil_res <- fit_compartment(
  analysis_df,
  "Soil"
)

type3_all <- dplyr::bind_rows(
  leaf_res$type3,
  soil_res$type3
)

pairwise_all <- dplyr::bind_rows(
  leaf_res$pairwise_tukey,
  soil_res$pairwise_tukey
)

# -------------------------
# Table S12
# -------------------------

table_s12 <- type3_all %>%
  dplyr::mutate(
    `χ²` =
      round(
        Chi_square,
        2
      ),
    `p-value` =
      dplyr::if_else(
        p_value < 0.001,
        "< 0.001",
        sprintf(
          "%.4f",
          p_value
        )
      )
  ) %>%
  dplyr::select(
    Compartment,
    Effect,
    `χ²`,
    df,
    `p-value`
  )

# -------------------------
# Table S13
# Comparisons versus NC extracted from the
# full Tukey-adjusted pairwise family
# -------------------------

sig_code <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# emmeans may output either "A - B" or "B - A".
# Re-orient comparisons so that all reported effects are
# treatment minus NC, matching the manuscript table.
vs_nc <- pairwise_all %>%
  dplyr::filter(
    grepl("NC", contrast, fixed = TRUE)
  ) %>%
  dplyr::mutate(
    contrast = as.character(contrast),
    treatment_name = dplyr::case_when(
      grepl(" - NC$", contrast) ~
        sub(" - NC$", "", contrast),
      grepl("^NC - ", contrast) ~
        sub("^NC - ", "", contrast),
      TRUE ~ NA_character_
    ),
    estimate_oriented = dplyr::case_when(
      grepl(" - NC$", contrast) ~
        estimate,
      grepl("^NC - ", contrast) ~
        -estimate,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(
    !is.na(treatment_name),
    treatment_name != "NC"
  )

contrast_labels <- c(
  "1013" =
    "Pseudomonas putida 1013 - NC",
  "1013_plus_B.C" =
    "Pseudomonas putida 1013 + B. cinerea - NC",
  "1017" =
    "Pantoea agglomerans 1017 - NC",
  "1017_plus_B.C" =
    "Pantoea agglomerans 1017 + B. cinerea - NC",
  "B.C" =
    "B. Cinerea - NC"
)

table_s13_long <- vs_nc %>%
  dplyr::mutate(
    Time =
      as.character(Time),
    Comparison =
      unname(
        contrast_labels[
          treatment_name
        ]
      ),
    Significance =
      sig_code(
        p.value
      ),
    Display =
      paste0(
        sprintf(
          "%.3f",
          estimate_oriented
        ),
        Significance
      )
  ) %>%
  dplyr::select(
    Compartment,
    Comparison,
    Time,
    Display,
    estimate_oriented,
    SE,
    df,
    t.ratio,
    p.value
  )

table_s13 <- table_s13_long %>%
  dplyr::select(
    Compartment,
    Comparison,
    Time,
    Display
  ) %>%
  tidyr::pivot_wider(
    names_from = Time,
    values_from = Display
  )

# -------------------------
# Model diagnostics
# -------------------------

extract_diagnostics <- function(
  model,
  compartment_name
) {
  vc <- as.data.frame(
    lme4::VarCorr(model)
  )

  plant_var <- vc %>%
    dplyr::filter(
      grp == "PlantID"
    ) %>%
    dplyr::pull(vcov)

  residual_var <-
    stats::sigma(model)^2

  tibble::tibble(
    Compartment =
      compartment_name,
    Observations =
      stats::nobs(model),
    Plant_random_intercept_variance =
      plant_var,
    Residual_variance =
      residual_var,
    ICC =
      plant_var /
      (
        plant_var +
        residual_var
      ),
    Singular_fit =
      lme4::isSingular(
        model,
        tol = 1e-5
      )
  )
}

diagnostics <- dplyr::bind_rows(
  extract_diagnostics(
    leaf_res$model,
    "Leaf"
  ),
  extract_diagnostics(
    soil_res$model,
    "Soil"
  )
)

# -------------------------
# Export
# -------------------------

readr::write_csv(
  table_s11,
  file.path(
    output_dir,
    "Table_S11_Bacterial_Populations.csv"
  )
)

readr::write_csv(
  table_s12,
  file.path(
    output_dir,
    "Table_S12_Statistical_analysis.csv"
  )
)

readr::write_csv(
  table_s13,
  file.path(
    output_dir,
    "Table_S13_Pairwise_comparisons.csv"
  )
)

readr::write_csv(
  analysis_df,
  file.path(
    output_dir,
    "CFU_baseline_normalized_data.csv"
  )
)

readr::write_csv(
  type3_all,
  file.path(
    output_dir,
    "CFU_TypeIII_detailed.csv"
  )
)

readr::write_csv(
  table_s13_long,
  file.path(
    output_dir,
    "CFU_Tukey_vs_NC_detailed.csv"
  )
)

readr::write_csv(
  pairwise_all,
  file.path(
    output_dir,
    "CFU_Tukey_all_pairwise.csv"
  )
)

readr::write_csv(
  diagnostics,
  file.path(
    output_dir,
    "CFU_model_diagnostics.csv"
  )
)


cat(
  "\nCFU manuscript-method analysis completed successfully.\n"
)
cat(
  "Input file:",
  input_file,
  "\n"
)
cat(
  "Model: Delta_log10_CFU ~ Treatment * Time + (1 | PlantID)\n"
)
cat(
  "T0 used only for baseline normalization; models fitted to T1-T4.\n"
)
cat(
  "Post hoc comparisons: EMMs with Tukey adjustment.\n"
)
