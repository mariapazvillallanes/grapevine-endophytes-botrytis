# ============================================================
# 01 – Indole-related compound production
# Clean Zenodo reproducibility script
# ============================================================

local({
# ============================================================
# Table S2 – Time-course production of indole-related compounds
# Reproducible analysis from the original biological replicates
# ============================================================

# Input:
#   data/01_indole_production.csv
#
# Experimental structure:
#   11 strains
#   3 incubation times: 72, 120 and 168 h
#   2 tryptophan conditions: -Trp and +Trp
#   3 independent replicates per cell
#
# Statistics:
#   Kruskal-Wallis followed by Dunn tests with BH correction
#   within each time × tryptophan condition.
#
# Compact letters are calculated from the adjusted pairwise
# p-values for each condition independently.

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "rstatix",
  "multcompView",
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
  "01_indole_production.csv"
)

output_dir <- file.path(project_dir, "output", "01_indole_production")

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

iaa <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(
    decimal_mark = ","
  ),
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
)

required_columns <- c(
  "Species_strain",
  "Species",
  "Strain",
  "Time_h",
  "Trp",
  "Replicate",
  "IAA"
)

missing_columns <- setdiff(
  required_columns,
  names(iaa)
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

iaa <- iaa %>%
  dplyr::mutate(
    Species_strain =
      as.character(Species_strain),
    Species =
      as.character(Species),
    Strain =
      as.character(Strain),
    Time_h =
      as.numeric(Time_h),
    Trp =
      as.character(Trp),
    Replicate =
      as.character(Replicate),
    IAA =
      as.numeric(IAA)
  )

expected_n <- 11 * 3 * 2 * 3

if (nrow(iaa) != expected_n) {
  warning(
    "Expected ",
    expected_n,
    " observations, but found ",
    nrow(iaa),
    "."
  )
}

rep_check <- iaa %>%
  dplyr::count(
    Species_strain,
    Time_h,
    Trp,
    name = "n"
  )

if (any(rep_check$n != 3L)) {
  stop(
    "Each strain × time × Trp condition must contain exactly three replicates."
  )
}

# -------------------------
# Descriptive statistics
# -------------------------

summary_tbl <- iaa %>%
  dplyr::group_by(
    Species_strain,
    Species,
    Strain,
    Time_h,
    Trp
  ) %>%
  dplyr::summarise(
    n =
      dplyr::n(),
    Mean =
      mean(IAA),
    SD =
      stats::sd(IAA),
    .groups = "drop"
  )

# -------------------------
# Kruskal-Wallis and Dunn-BH
# -------------------------

kruskal_results <- iaa %>%
  dplyr::group_by(
    Time_h,
    Trp
  ) %>%
  rstatix::kruskal_test(
    IAA ~ Species_strain
  ) %>%
  dplyr::ungroup()

dunn_results <- iaa %>%
  dplyr::group_by(
    Time_h,
    Trp
  ) %>%
  rstatix::dunn_test(
    IAA ~ Species_strain,
    p.adjust.method = "BH"
  ) %>%
  dplyr::ungroup()

# -------------------------
# Compact letter display
# -------------------------

get_letters <- function(
  data,
  time_value,
  trp_value
) {

  d <- data %>%
    dplyr::filter(
      Time_h == time_value,
      Trp == trp_value
    )

  kw <- stats::kruskal.test(
    IAA ~ Species_strain,
    data = d
  )

  if (kw$p.value >= 0.05) {
    return(
      d %>%
        dplyr::distinct(
          Species_strain
        ) %>%
        dplyr::mutate(
          Letters = "a",
          Time_h = time_value,
          Trp = trp_value
        )
    )
  }

  dunn <- d %>%
    rstatix::dunn_test(
      IAA ~ Species_strain,
      p.adjust.method = "BH"
    ) %>%
    dplyr::filter(
      !is.na(p.adj)
    )

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
    Species_strain = names(cld),
    Letters = unname(cld),
    Time_h = time_value,
    Trp = trp_value
  )
}

letter_grid <- tidyr::expand_grid(
  Time_h = c(
    72,
    120,
    168
  ),
  Trp = c(
    "-Trp",
    "+Trp"
  )
)

letters_tbl <- purrr::map2_dfr(
  letter_grid$Time_h,
  letter_grid$Trp,
  ~ get_letters(
    iaa,
    .x,
    .y
  )
)

summary_letters <- summary_tbl %>%
  dplyr::left_join(
    letters_tbl,
    by = c(
      "Species_strain",
      "Time_h",
      "Trp"
    )
  ) %>%
  dplyr::mutate(
    Mean_SD =
      sprintf(
        "%.2f ± %.2f",
        Mean,
        SD
      ),
    Mean_SD_letters =
      paste0(
        Mean_SD,
        Letters
      ),
    Condition =
      paste0(
        Time_h,
        " h ",
        Trp
      )
  )

# -------------------------
# Derived 168 h parameters
# -------------------------

derived_168 <- summary_tbl %>%
  dplyr::filter(
    Time_h == 168
  ) %>%
  dplyr::select(
    Species_strain,
    Trp,
    Mean
  ) %>%
  tidyr::pivot_wider(
    names_from = Trp,
    values_from = Mean
  ) %>%
  dplyr::mutate(
    `Delta indole-related compounds at 168 h` =
      round(`+Trp` - `-Trp`, 2),
    `Induction at 168 h (%)` =
      round(
        100 * (`+Trp` - `-Trp`) / `-Trp`,
        1
      ),
    `Fold change at 168 h (+Trp/-Trp)` =
      round(`+Trp` / `-Trp`, 2)
  )

# -------------------------
# Final Table S2
# -------------------------

condition_order <- c(
  "72 h -Trp",
  "120 h -Trp",
  "168 h -Trp",
  "72 h +Trp",
  "120 h +Trp",
  "168 h +Trp"
)

table_s2 <- summary_letters %>%
  dplyr::select(
    Species_strain,
    Condition,
    Mean_SD_letters
  ) %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = Mean_SD_letters
  ) %>%
  dplyr::left_join(
    derived_168 %>%
      dplyr::select(
        Species_strain,
        `Delta indole-related compounds at 168 h`,
        `Induction at 168 h (%)`,
        `Fold change at 168 h (+Trp/-Trp)`
      ),
    by = "Species_strain"
  ) %>%
  dplyr::select(
    Strain = Species_strain,
    dplyr::all_of(
      condition_order
    ),
    `Delta indole-related compounds at 168 h`,
    `Induction at 168 h (%)`,
    `Fold change at 168 h (+Trp/-Trp)`
  )

# -------------------------
# Export
# -------------------------

readr::write_csv(
  table_s2,
  file.path(
    output_dir,
    "Table_S2_Indole_Compounds.csv"
  )
)

readr::write_csv(
  kruskal_results,
  file.path(
    output_dir,
    "Table_S2_Kruskal_Wallis.csv"
  )
)

readr::write_csv(
  dunn_results,
  file.path(
    output_dir,
    "Table_S2_Dunn_BH.csv"
  )
)

readr::write_csv(
  letters_tbl,
  file.path(
    output_dir,
    "Table_S2_letters.csv"
  )
)


cat(
  "\nTable S2 analysis completed successfully.\n"
)

cat(
  "Input file:",
  input_file,
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)

})

# ---- Supplementary Figure S1 ----
local({
# ============================================================
# Figure S1 – IAA production heatmap
# Reproducible analysis for Supplementary Figure S1
# ============================================================

# Input:
#   data/01_indole_production.csv
#
# The input contains the three original biological replicates
# for each strain, incubation time and tryptophan condition.
#
# Workflow:
#   1. Calculate the mean IAA value from the three replicates.
#   2. Build the strain × condition matrix.
#   3. Calculate a Z-score independently for each strain.
#   4. Generate Figure S1.


# -------------------------
# 1) Required packages
# -------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "ComplexHeatmap",
  "circlize",
  "ragg",
  "systemfonts",
  "gridtext"
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

graphics.off()

if (requireNamespace("showtext", quietly = TRUE)) {
  showtext::showtext_auto(enable = FALSE)
}


# -------------------------
# 2) Input and output paths
# -------------------------

# The script can be run either from the project root
# or from the /scripts directory.

if (dir.exists("data")) {
  project_dir <- "."
} else if (dir.exists("../data")) {
  project_dir <- ".."
} else {
  stop(
    "Project directory not found. ",
    "Run the script from the project root or from the /scripts directory."
  )
}

input_file <- file.path(
  project_dir,
  "data",
  "01_indole_production.csv"
)

output_dir <- file.path(project_dir, "output", "01_indole_production")

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
# 3) Read and validate data
# -------------------------

iaa_long <- readr::read_delim(
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
)

required_columns <- c(
  "Species_strain",
  "Species",
  "Strain",
  "Time_h",
  "Trp",
  "Replicate",
  "IAA"
)

missing_columns <- setdiff(
  required_columns,
  names(iaa_long)
)

if (length(missing_columns) > 0) {
  stop(
    "The input file is missing the following required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

iaa_long <- iaa_long %>%
  dplyr::mutate(
    Species_strain =
      trimws(
        as.character(Species_strain)
      ),
    Species =
      trimws(
        as.character(Species)
      ),
    Strain =
      as.character(Strain),
    Time_h =
      as.numeric(Time_h),
    Trp =
      as.character(Trp),
    Replicate =
      as.character(Replicate),
    IAA =
      as.numeric(IAA)
  )

expected_times <- c(
  72,
  120,
  168
)

expected_trp <- c(
  "-Trp",
  "+Trp"
)

expected_replicates <- c(
  "R1",
  "R2",
  "R3"
)

if (
  anyNA(
    iaa_long[
      ,
      required_columns
    ]
  )
) {
  stop(
    "Missing values were detected in 01_indole_production.csv."
  )
}

replicate_check <- iaa_long %>%
  dplyr::count(
    Species_strain,
    Time_h,
    Trp,
    name = "n_replicates"
  )

if (
  any(
    replicate_check$n_replicates != 3L
  )
) {
  print(
    replicate_check %>%
      dplyr::filter(
        n_replicates != 3L
      )
  )

  stop(
    "Every strain × time × tryptophan condition must contain ",
    "exactly three biological replicates."
  )
}

if (
  !setequal(
    unique(iaa_long$Time_h),
    expected_times
  )
) {
  stop(
    "Unexpected incubation times were detected."
  )
}

if (
  !setequal(
    unique(iaa_long$Trp),
    expected_trp
  )
) {
  stop(
    "Unexpected tryptophan conditions were detected."
  )
}

if (
  !setequal(
    unique(iaa_long$Replicate),
    expected_replicates
  )
) {
  stop(
    "Unexpected replicate identifiers were detected."
  )
}


# -------------------------
# 4) Mean IAA values
# -------------------------

iaa_summary <- iaa_long %>%
  dplyr::group_by(
    Species,
    Strain,
    Time_h,
    Trp
  ) %>%
  dplyr::summarise(
    n =
      dplyr::n(),
    IAA_mean =
      mean(IAA),
    IAA_SD =
      stats::sd(IAA),
    .groups = "drop"
  )


# -------------------------
# 5) Heatmap matrix
# -------------------------

heatmap_matrix <- iaa_summary %>%
  dplyr::mutate(
    Condition = paste0(
      Time_h,
      " h ",
      Trp
    ),
    Strain_label = paste0(
      Species,
      " ",
      Strain
    )
  ) %>%
  dplyr::select(
    Strain_label,
    Condition,
    IAA_mean
  ) %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = IAA_mean
  ) %>%
  tibble::column_to_rownames(
    "Strain_label"
  ) %>%
  as.matrix()

condition_order <- c(
  "72 h -Trp",
  "120 h -Trp",
  "168 h -Trp",
  "72 h +Trp",
  "120 h +Trp",
  "168 h +Trp"
)

if (
  !all(
    condition_order %in%
      colnames(heatmap_matrix)
  )
) {
  stop(
    "Not all expected experimental conditions were found."
  )
}

heatmap_matrix <- heatmap_matrix[
  ,
  condition_order,
  drop = FALSE
]


# -------------------------
# 6) Z-score by strain
# -------------------------

# Each row is standardized independently, matching
# the original analysis.

heatmap_scaled <- t(
  scale(
    t(heatmap_matrix)
  )
)

if (
  any(
    !is.finite(
      heatmap_scaled
    )
  )
) {
  stop(
    "A Z-score could not be calculated for at least one strain. ",
    "Check for missing data or a row without variation."
  )
}


# -------------------------
# 7) Scientific-name labels
# -------------------------

row_label_html <- vapply(
  rownames(heatmap_scaled),
  function(x) {

    species <- sub(
      "\\s+[0-9]+$",
      "",
      x
    )

    strain <- sub(
      "^.*\\s+",
      "",
      x
    )

    genus <- strsplit(
      species,
      " ",
      fixed = TRUE
    )[[1]][1]

    epithet <- sub(
      "^[^ ]+\\s+",
      "",
      species
    )

    genus_abbrev <- if (
      genus == "Pseudomonas"
    ) {

      "Ps."

    } else if (
      genus == "Pantoea"
    ) {

      "Pa."

    } else {

      paste0(
        substr(
          genus,
          1,
          1
        ),
        "."
      )
    }

    paste0(
      "<i>",
      genus_abbrev,
      " ",
      epithet,
      "</i> ",
      strain
    )
  },
  character(1)
)


# -------------------------
# 8) Font configuration
# -------------------------

font_matches <- systemfonts::system_fonts()

if (
  !any(
    tolower(
      font_matches$family
    ) == "times new roman"
  )
) {
  stop(
    "Times New Roman was not found on this computer."
  )
}

font_family <- "Times New Roman"


# -------------------------
# 9) Figure S1
# -------------------------

col_fun <- circlize::colorRamp2(
  c(
    -1.5,
    0,
    1.5
  ),
  c(
    "#9ecae1",
    "#ffffcc",
    "#d7301f"
  )
)

draw_figure_s1 <- function() {

  figure_s1 <- ComplexHeatmap::Heatmap(
    heatmap_scaled,

    name = "Z-score",
    col = col_fun,

    cluster_rows = TRUE,
    cluster_columns = FALSE,

    row_dend_side = "left",
    row_dend_width =
      grid::unit(
        10,
        "mm"
      ),

    row_names_side = "right",

    row_labels =
      ComplexHeatmap::gt_render(
        row_label_html,
        align_widths = TRUE,
        padding =
          grid::unit(
            c(
              0,
              0,
              0,
              0
            ),
            "pt"
          )
      ),

    row_names_gp =
      grid::gpar(
        fontsize = 9.5,
        fontfamily =
          font_family
      ),

    row_names_max_width =
      grid::unit(
        48,
        "mm"
      ),

    column_names_gp =
      grid::gpar(
        fontsize = 9.5,
        fontfamily =
          font_family
      ),

    column_names_rot = 45,
    column_names_centered = TRUE,

    column_title = NULL,
    row_title = NULL,

    width =
      grid::unit(
        100,
        "mm"
      ),

    height =
      grid::unit(
        78,
        "mm"
      ),

    border = FALSE,

    rect_gp =
      grid::gpar(
        col = "white",
        lwd = 0.35
      ),

    heatmap_legend_param =
      list(
        title = "Z-score",
        direction = "horizontal",
        title_position = "leftcenter",

        legend_width =
          grid::unit(
            65,
            "mm"
          ),

        at = c(
          -1.5,
          -1,
          -0.5,
          0,
          0.5,
          1,
          1.5
        ),

        title_gp =
          grid::gpar(
            fontsize = 9,
            fontfamily =
              font_family
          ),

        labels_gp =
          grid::gpar(
            fontsize = 8,
            fontfamily =
              font_family
          )
      )
  )

  ComplexHeatmap::draw(
    figure_s1,
    heatmap_legend_side = "bottom",
    padding =
      grid::unit(
        c(
          4,
          4,
          4,
          4
        ),
        "mm"
      )
  )
}


# TIFF
ragg::agg_tiff(
  filename = file.path(
    output_dir,
    "Figure_S1.tiff"
  ),
  width = 8.2,
  height = 5.2,
  units = "in",
  res = 600,
  compression = "lzw"
)

draw_figure_s1()

grDevices::dev.off()


# PDF
grDevices::cairo_pdf(
  filename = file.path(
    output_dir,
    "Figure_S1.pdf"
  ),
  width = 8.2,
  height = 5.2,
  family = "sans"
)

draw_figure_s1()

grDevices::dev.off()


# -------------------------
# 10) Export associated data
# -------------------------

readr::write_csv(
  iaa_summary,
  file.path(
    output_dir,
    "Figure_S1_mean_IAA.csv"
  )
)

z_scores_output <- as.data.frame(
  heatmap_scaled,
  check.names = FALSE
) %>%
  tibble::rownames_to_column(
    "Strain"
  )

readr::write_csv(
  z_scores_output,
  file.path(
    output_dir,
    "Figure_S1_Z_scores.csv"
  )
)


# -------------------------
# 11) Completion message
# -------------------------

cat(
  "\nFigure S1 analysis completed successfully.\n"
)

cat(
  "Input file:",
  input_file,
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)

})
