# ============================================================
# 02 – Mineral solubilization
# Clean Zenodo reproducibility script
# ============================================================

# ============================================================
# Table S3 – Solubilization indices
# Reproducible analysis for phosphate, potassium and zinc
# solubilization
# ============================================================

# Input:
#   data/02_mineral_solubilization.csv
#
# SI = (colony diameter + halo diameter) / colony diameter
# SI = 1  -> no detectable solubilization
# SI > 1  -> positive solubilization
#
# Statistics:
#   Kruskal-Wallis followed by Dunn tests with BH correction.
#   Compact letters are calculated only among strains with
#   detectable activity (SI > 1) for each mineral.

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "rstatix",
  "multcompView",
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

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

input_file <- file.path(project_dir, "data", "02_mineral_solubilization.csv")
output_dir <- file.path(project_dir, "output", "02_mineral_solubilization")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# -------------------------
# Read data
# -------------------------

raw <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(decimal_mark = ","),
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
)

required_columns <- c("Strain", "P", "K", "Zn")
missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

strain_key <- tibble::tribble(
  ~Strain_ID, ~Species,
  "1013", "Pseudomonas putida",
  "1014", "Stenotrophomonas maltophilia",
  "1015", "Erwinia billingiae",
  "1016", "Erwinia billingiae",
  "1017", "Pantoea agglomerans",
  "1018", "Erwinia billingiae",
  "1090", "Stenotrophomonas maltophilia",
  "890",  "Penicillium chrysogenum",
  "932",  "Nothophoma quercina",
  "950",  "Fusarium avenaceum",
  "1107", "Aureobasidium pullulans"
)

si_long <- raw %>%
  dplyr::mutate(
    Strain = as.character(Strain),
    Strain_ID = stringr::str_extract(Strain, "^[0-9]+"),
    Replicate = stringr::str_extract(Strain, "(?<=\\.)[0-9]+")
  ) %>%
  dplyr::left_join(strain_key, by = "Strain_ID") %>%
  dplyr::mutate(
    Species = dplyr::if_else(is.na(Species), Strain_ID, Species),
    Strain_label = paste(Species, Strain_ID)
  ) %>%
  tidyr::pivot_longer(
    cols = c(P, K, Zn),
    names_to = "Nutrient",
    values_to = "SI"
  ) %>%
  dplyr::mutate(
    SI = as.numeric(SI),
    Positive = SI > 1
  )

# -------------------------
# Summary statistics
# -------------------------

si_summary <- si_long %>%
  dplyr::group_by(
    Species,
    Strain_ID,
    Strain_label,
    Nutrient
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    Mean_SI = mean(SI, na.rm = TRUE),
    SD_SI = stats::sd(SI, na.rm = TRUE),
    SEM_SI = SD_SI / sqrt(n),
    Positive_replicates = sum(Positive, na.rm = TRUE),
    Result = dplyr::case_when(
      Positive_replicates == 0 ~ "Negative",
      Positive_replicates == n ~ "Positive",
      TRUE ~ "Variable"
    ),
    Mean_SD = sprintf("%.2f ± %.2f", Mean_SI, SD_SI),
    .groups = "drop"
  )

# -------------------------
# Statistics
# -------------------------

variable_nutrients <- si_long %>%
  dplyr::group_by(Nutrient) %>%
  dplyr::summarise(
    Variation = dplyr::n_distinct(SI) > 1,
    .groups = "drop"
  ) %>%
  dplyr::filter(Variation) %>%
  dplyr::pull(Nutrient)

kruskal_results <- si_long %>%
  dplyr::filter(Nutrient %in% variable_nutrients) %>%
  dplyr::group_by(Nutrient) %>%
  rstatix::kruskal_test(SI ~ Strain_label) %>%
  dplyr::ungroup()

dunn_results <- si_long %>%
  dplyr::filter(Nutrient %in% variable_nutrients) %>%
  dplyr::group_by(Nutrient) %>%
  rstatix::dunn_test(
    SI ~ Strain_label,
    p.adjust.method = "BH"
  ) %>%
  dplyr::ungroup()

get_dunn_letters <- function(
  data,
  nutrient_name,
  positive_only = TRUE
) {

  df <- data %>%
    dplyr::filter(Nutrient == nutrient_name)

  if (positive_only) {
    df <- df %>%
      dplyr::filter(SI > 1)
  }

  if (
    dplyr::n_distinct(df$Strain_label) < 2 ||
    dplyr::n_distinct(df$SI) < 2
  ) {
    return(
      df %>%
        dplyr::distinct(Nutrient, Strain_label) %>%
        dplyr::mutate(Letters = "a")
    )
  }

  kw <- stats::kruskal.test(
    SI ~ Strain_label,
    data = df
  )

  if (kw$p.value >= 0.05) {
    return(
      df %>%
        dplyr::distinct(Nutrient, Strain_label) %>%
        dplyr::mutate(Letters = "a")
    )
  }

  dunn <- df %>%
    rstatix::dunn_test(
      SI ~ Strain_label,
      p.adjust.method = "BH"
    ) %>%
    dplyr::filter(!is.na(p.adj))

  pvals <- dunn$p.adj
  names(pvals) <- paste(
    dunn$group1,
    dunn$group2,
    sep = "-"
  )

  cld <- multcompView::multcompLetters(
    pvals,
    compare = "<",
    threshold = 0.05
  )$Letters

  tibble::tibble(
    Nutrient = nutrient_name,
    Strain_label = names(cld),
    Letters = unname(cld)
  )
}

solub_letters <- purrr::map_dfr(
  c("P", "K"),
  ~ get_dunn_letters(
    si_long,
    .x,
    positive_only = TRUE
  )
)

si_summary_letters <- si_summary %>%
  dplyr::left_join(
    solub_letters,
    by = c(
      "Nutrient",
      "Strain_label"
    )
  ) %>%
  dplyr::mutate(
    Mean_SD_letters = dplyr::case_when(
      Result == "Negative" ~ "1.00 ± 0.00",
      is.na(Letters) ~ Mean_SD,
      TRUE ~ paste0(Mean_SD, Letters)
    )
  )

# -------------------------
# Table S3
# -------------------------

table_s3 <- si_summary_letters %>%
  dplyr::select(
    Species,
    Strain_ID,
    Strain_label,
    Nutrient,
    Mean_SD_letters,
    Mean_SI
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(
      Species,
      Strain_ID,
      Strain_label
    ),
    names_from = Nutrient,
    values_from = c(
      Mean_SD_letters,
      Mean_SI
    ),
    names_glue = "{Nutrient}_{.value}"
  ) %>%
  dplyr::mutate(
    Solubilization_profile = dplyr::case_when(
      P_Mean_SI > 1 & K_Mean_SI > 1 & Zn_Mean_SI > 1 ~
        "P/K/Zn solubilizer",
      P_Mean_SI > 1 & K_Mean_SI > 1 & Zn_Mean_SI <= 1 ~
        "P/K solubilizer",
      P_Mean_SI > 1 & K_Mean_SI <= 1 & Zn_Mean_SI <= 1 ~
        "P solubilizer",
      P_Mean_SI <= 1 & K_Mean_SI > 1 & Zn_Mean_SI <= 1 ~
        "K solubilizer",
      P_Mean_SI <= 1 & K_Mean_SI <= 1 & Zn_Mean_SI > 1 ~
        "Zn solubilizer",
      TRUE ~
        "No detectable activity"
    )
  ) %>%
  dplyr::select(
    Strain = Strain_label,
    `P SI` = P_Mean_SD_letters,
    `K SI` = K_Mean_SD_letters,
    `Zn SI` = Zn_Mean_SD_letters,
    Solubilization_profile
  )

# -------------------------
# Export
# -------------------------

readr::write_csv(
  table_s3,
  file.path(
    output_dir,
    "Table_S3_Summary_SI.csv"
  )
)

readr::write_csv(
  kruskal_results,
  file.path(
    output_dir,
    "Table_S3_Kruskal_Wallis.csv"
  )
)

readr::write_csv(
  dunn_results,
  file.path(
    output_dir,
    "Table_S3_Dunn_BH.csv"
  )
)


cat("\nTable S3 analysis completed successfully.\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
