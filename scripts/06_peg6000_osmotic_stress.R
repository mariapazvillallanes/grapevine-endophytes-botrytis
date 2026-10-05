# ============================================================
# 06 – PEG 6000 osmotic stress
# Clean Zenodo reproducibility script
# ============================================================

# ============================================================
# Figure 3 – Growth response to PEG 6000
# Reproducible analysis for Figure 3, Figure S4 and associated
# supplementary statistical outputs
# ============================================================

# Experimental structure:
#   8 non-filamentous isolates
#   5 PEG 6000 concentrations: 0, 5, 10, 20 and 40% (w/v)
#   3 independent biological cultures per isolate
#
# Important:
#   OD600 values in 06_peg6000_osmotic_stress.csv are already background-corrected
#   using the corresponding uninoculated medium + PEG blank.
#   No additional blank subtraction is performed in this script.
#
# Main inferential analyses:
#   - Table S9: bacteria-only split-plot repeated-measures ANOVA
#   - Table S10 / Figure 3: bacteria-only relative AUC at 40% PEG
#     followed by Tukey-adjusted comparisons
#
# Aureobasidium pullulans 1107 was cultured in YPD and is retained
# descriptively but excluded from bacteria-only inferential tests.


# -------------------------
# 1) Required packages
# -------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "multcompView",
  "systemfonts",
  "sysfonts",
  "showtext"
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
  "06_peg6000_osmotic_stress.csv"
)

output_dir <- file.path(project_dir, "output", "06_peg6000_osmotic_stress")

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
# 3) Helper functions
# -------------------------

to_numeric_safely <- function(x) {
  suppressWarnings(
    as.numeric(
      stringr::str_replace_all(
        as.character(x),
        ",",
        "."
      )
    )
  )
}

trapz_auc <- function(
  time,
  response
) {

  ord <- order(time)
  time <- time[ord]
  response <- response[ord]

  sum(
    diff(time) *
      (
        head(response, -1) +
          tail(response, -1)
      ) / 2
  )
}

p_display <- function(p) {
  ifelse(
    p < 0.001,
    "<0.001",
    sprintf("%.4f", p)
  )
}

sig_label <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

normalise_letters <- function(
  raw_letters,
  ordered_groups
) {

  mapping <- character(0)
  next_letter_index <- 1L

  for (group_name in ordered_groups) {

    symbols <- strsplit(
      raw_letters[[group_name]],
      split = ""
    )[[1]]

    for (symbol in symbols) {

      if (!(symbol %in% names(mapping))) {
        mapping[symbol] <- letters[next_letter_index]
        next_letter_index <- next_letter_index + 1L
      }
    }
  }

  vapply(
    raw_letters,
    function(letter_string) {

      symbols <- strsplit(
        letter_string,
        split = ""
      )[[1]]

      mapped <- unname(
        mapping[symbols]
      )

      mapped <- mapped[
        order(
          match(
            mapped,
            letters
          )
        )
      ]

      paste0(
        mapped,
        collapse = ""
      )
    },
    character(1)
  )
}


# -------------------------
# 4) Read and prepare data
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
  dplyr::rename_with(
    ~ tolower(
      trimws(.x)
    )
  )

required_columns <- c(
  "strain",
  "peg",
  "time",
  "n1",
  "n2",
  "n3"
)

missing_columns <- setdiff(
  required_columns,
  names(raw)
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

dat <- raw %>%
  dplyr::transmute(
    strain = stringr::str_trim(
      as.character(strain)
    ),
    peg = to_numeric_safely(peg),
    time = to_numeric_safely(time),
    n1 = to_numeric_safely(n1),
    n2 = to_numeric_safely(n2),
    n3 = to_numeric_safely(n3)
  ) %>%
  dplyr::mutate(
    strain = stringr::str_replace(
      strain,
      "\\.0$",
      ""
    )
  )

if (
  any(is.na(dat$peg)) ||
  any(is.na(dat$time))
) {
  stop(
    "Non-numeric PEG concentration or time values were detected."
  )
}

if (
  any(is.na(dat$n1)) ||
  any(is.na(dat$n2)) ||
  any(is.na(dat$n3))
) {
  stop(
    "Missing or non-numeric OD600 measurements were detected."
  )
}


# -------------------------
# 5) Experimental design
# -------------------------

bacterial_strains <- c(
  "1013",
  "1014",
  "1015",
  "1016",
  "1017",
  "1018",
  "1090"
)

yeast_strain <- "1107"

microbial_strains <- c(
  bacterial_strains,
  yeast_strain
)

expected_pegs <- c(
  0,
  5,
  10,
  20,
  40
)

microbial_wide <- dat %>%
  dplyr::filter(
    strain %in% microbial_strains
  )

control_wide <- dat %>%
  dplyr::filter(
    strain %in% c(
      "TSB",
      "YPD"
    )
  )

missing_strains <- setdiff(
  microbial_strains,
  unique(
    microbial_wide$strain
  )
)

if (length(missing_strains) > 0) {
  stop(
    "Missing microbial strains: ",
    paste(
      missing_strains,
      collapse = ", "
    )
  )
}

observed_pegs <- sort(
  unique(
    microbial_wide$peg
  )
)

if (
  !identical(
    observed_pegs,
    expected_pegs
  )
) {
  stop(
    "Microbial PEG concentrations must be exactly: ",
    paste(
      expected_pegs,
      collapse = ", "
    ),
    "%."
  )
}

duplicate_check <- microbial_wide %>%
  dplyr::count(
    strain,
    peg,
    time,
    name = "n_rows"
  )

if (
  any(
    duplicate_check$n_rows != 1L
  )
) {

  print(
    duplicate_check %>%
      dplyr::filter(
        n_rows != 1L
      )
  )

  stop(
    "Duplicated or missing strain × PEG × time rows were detected."
  )
}

microbial_long <- microbial_wide %>%
  tidyr::pivot_longer(
    cols = c(
      n1,
      n2,
      n3
    ),
    names_to = "replicate_source",
    values_to = "OD600_corrected"
  ) %>%
  dplyr::mutate(
    biological_replicate =
      dplyr::case_when(
        replicate_source == "n1" ~ 1L,
        replicate_source == "n2" ~ 2L,
        replicate_source == "n3" ~ 3L,
        TRUE ~ NA_integer_
      ),
    medium = dplyr::if_else(
      strain == yeast_strain,
      "YPD",
      "TSB"
    ),
    culture_id = interaction(
      strain,
      biological_replicate,
      drop = TRUE
    )
  )

expected_times <- sort(
  unique(
    microbial_long$time
  )
)

trajectory_check <- microbial_long %>%
  dplyr::group_by(
    strain,
    peg,
    biological_replicate
  ) %>%
  dplyr::summarise(
    n_times = dplyr::n_distinct(
      time
    ),
    .groups = "drop"
  )

if (
  any(
    trajectory_check$n_times !=
      length(expected_times)
  )
) {

  print(
    trajectory_check %>%
      dplyr::filter(
        n_times !=
          length(expected_times)
      )
  )

  stop(
    "At least one microbial trajectory is incomplete."
  )
}

expected_microbe_observations <-
  length(microbial_strains) *
  length(expected_pegs) *
  length(expected_times) *
  3L

if (
  nrow(microbial_long) !=
    expected_microbe_observations
) {
  stop(
    "Unexpected number of microbial observations. Expected ",
    expected_microbe_observations,
    " but found ",
    nrow(microbial_long),
    "."
  )
}

# OD600 values are already background-corrected.
# TSB and YPD control rows are retained only for audit
# and are not subtracted again.

initial_check <- microbial_long %>%
  dplyr::filter(
    time == min(expected_times)
  ) %>%
  dplyr::group_by(
    strain,
    peg
  ) %>%
  dplyr::summarise(
    Mean_initial_corrected_OD =
      mean(OD600_corrected),
    SD_initial_corrected_OD =
      stats::sd(OD600_corrected),
    Minimum =
      min(OD600_corrected),
    Maximum =
      max(OD600_corrected),
    .groups = "drop"
  )

control_summary <- if (
  nrow(control_wide) > 0
) {

  control_wide %>%
    tidyr::pivot_longer(
      cols = c(
        n1,
        n2,
        n3
      ),
      names_to = "blank_well",
      values_to = "Control_OD600"
    ) %>%
    dplyr::group_by(
      strain,
      peg,
      time
    ) %>%
    dplyr::summarise(
      Mean_control_OD600 =
        mean(Control_OD600),
      SD_control_OD600 =
        stats::sd(Control_OD600),
      n_control_wells =
        dplyr::n(),
      .groups = "drop"
    )

} else {

  tibble::tibble(
    strain = character(0),
    peg = numeric(0),
    time = numeric(0),
    Mean_control_OD600 = numeric(0),
    SD_control_OD600 = numeric(0),
    n_control_wells = integer(0)
  )
}


# -------------------------
# 6) AUC calculations
# -------------------------

auc_data <- microbial_long %>%
  dplyr::arrange(
    strain,
    peg,
    biological_replicate,
    time
  ) %>%
  dplyr::group_by(
    strain,
    medium,
    peg,
    biological_replicate,
    culture_id
  ) %>%
  dplyr::summarise(
    AUC_corrected = trapz_auc(
      time,
      OD600_corrected
    ),
    .groups = "drop"
  )

auc0 <- auc_data %>%
  dplyr::filter(
    peg == 0
  ) %>%
  dplyr::select(
    strain,
    biological_replicate,
    AUC_PEG0 = AUC_corrected
  )

if (
  any(
    auc0$AUC_PEG0 <= 0
  )
) {

  print(
    auc0 %>%
      dplyr::filter(
        AUC_PEG0 <= 0
      )
  )

  stop(
    "At least one PEG-free corrected AUC is zero or negative; ",
    "relative AUC cannot be calculated reliably."
  )
}

auc_data <- auc_data %>%
  dplyr::left_join(
    auc0,
    by = c(
      "strain",
      "biological_replicate"
    )
  ) %>%
  dplyr::mutate(
    Relative_AUC =
      100 *
      AUC_corrected /
      AUC_PEG0
  )


# -------------------------
# 7) Table S9
# Split-plot repeated-measures ANOVA
# -------------------------

auc_bacteria <- auc_data %>%
  dplyr::filter(
    strain %in%
      bacterial_strains
  )

a <- dplyr::n_distinct(
  auc_bacteria$strain
)

n <- dplyr::n_distinct(
  auc_bacteria$biological_replicate
)

b <- dplyr::n_distinct(
  auc_bacteria$peg
)

if (
  a != 7L ||
  n != 3L ||
  b != 5L
) {
  stop(
    "The bacteria-only split-plot analysis requires ",
    "7 strains, 3 cultures and 5 PEG conditions."
  )
}

grand_mean <- mean(
  auc_bacteria$AUC_corrected
)

strain_means <- auc_bacteria %>%
  dplyr::group_by(
    strain
  ) %>%
  dplyr::summarise(
    strain_mean =
      mean(AUC_corrected),
    .groups = "drop"
  )

culture_means <- auc_bacteria %>%
  dplyr::group_by(
    strain,
    biological_replicate
  ) %>%
  dplyr::summarise(
    culture_mean =
      mean(AUC_corrected),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    strain_means,
    by = "strain"
  )

peg_means <- auc_bacteria %>%
  dplyr::group_by(
    peg
  ) %>%
  dplyr::summarise(
    peg_mean =
      mean(AUC_corrected),
    .groups = "drop"
  )

strain_peg_means <- auc_bacteria %>%
  dplyr::group_by(
    strain,
    peg
  ) %>%
  dplyr::summarise(
    strain_peg_mean =
      mean(AUC_corrected),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    strain_means,
    by = "strain"
  ) %>%
  dplyr::left_join(
    peg_means,
    by = "peg"
  )

residual_data <- auc_bacteria %>%
  dplyr::left_join(
    culture_means %>%
      dplyr::select(
        strain,
        biological_replicate,
        culture_mean
      ),
    by = c(
      "strain",
      "biological_replicate"
    )
  ) %>%
  dplyr::left_join(
    strain_peg_means %>%
      dplyr::select(
        strain,
        peg,
        strain_peg_mean
      ),
    by = c(
      "strain",
      "peg"
    )
  ) %>%
  dplyr::left_join(
    strain_means,
    by = "strain"
  ) %>%
  dplyr::mutate(
    within_residual =
      AUC_corrected -
      culture_mean -
      strain_peg_mean +
      strain_mean
  )

SS_strain <- b * n * sum(
  (
    strain_means$strain_mean -
      grand_mean
  )^2
)

SS_culture_strain <- b * sum(
  (
    culture_means$culture_mean -
      culture_means$strain_mean
  )^2
)

SS_peg <- a * n * sum(
  (
    peg_means$peg_mean -
      grand_mean
  )^2
)

SS_interaction <- n * sum(
  (
    strain_peg_means$strain_peg_mean -
      strain_peg_means$strain_mean -
      strain_peg_means$peg_mean +
      grand_mean
  )^2
)

SS_peg_culture_strain <- sum(
  residual_data$within_residual^2
)

df_strain <- a - 1
df_culture_strain <- a * (n - 1)
df_peg <- b - 1
df_interaction <- (a - 1) * (b - 1)
df_peg_culture_strain <-
  a * (n - 1) * (b - 1)

MS_strain <-
  SS_strain /
  df_strain

MS_culture_strain <-
  SS_culture_strain /
  df_culture_strain

MS_peg <-
  SS_peg /
  df_peg

MS_interaction <-
  SS_interaction /
  df_interaction

MS_peg_culture_strain <-
  SS_peg_culture_strain /
  df_peg_culture_strain

F_strain <-
  MS_strain /
  MS_culture_strain

F_peg <-
  MS_peg /
  MS_peg_culture_strain

F_interaction <-
  MS_interaction /
  MS_peg_culture_strain

p_strain <- stats::pf(
  F_strain,
  df_strain,
  df_culture_strain,
  lower.tail = FALSE
)

p_peg <- stats::pf(
  F_peg,
  df_peg,
  df_peg_culture_strain,
  lower.tail = FALSE
)

p_interaction <- stats::pf(
  F_interaction,
  df_interaction,
  df_peg_culture_strain,
  lower.tail = FALSE
)

table_s9 <- tibble::tibble(
  Effect = c(
    "Bacterial strain",
    "PEG 6000 concentration",
    "Bacterial strain × PEG 6000 concentration"
  ),
  `Numerator df` = c(
    df_strain,
    df_peg,
    df_interaction
  ),
  `Denominator df` = c(
    df_culture_strain,
    df_peg_culture_strain,
    df_peg_culture_strain
  ),
  `F value` = c(
    F_strain,
    F_peg,
    F_interaction
  ),
  `p value` = c(
    p_display(p_strain),
    p_display(p_peg),
    p_display(p_interaction)
  ),
  Significance = c(
    sig_label(p_strain),
    sig_label(p_peg),
    sig_label(p_interaction)
  ),
  `Error term used for the F test` = c(
    "Culture nested within strain",
    "PEG × culture nested within strain",
    "PEG × culture nested within strain"
  )
)


# -------------------------
# 8) Table S10 and Figure 3
# Relative AUC at 40% PEG
# -------------------------

auc_40_bacteria <- auc_bacteria %>%
  dplyr::filter(
    peg == 40
  ) %>%
  dplyr::mutate(
    strain = factor(
      strain,
      levels = bacterial_strains
    )
  )

model_40 <- stats::aov(
  Relative_AUC ~ strain,
  data = auc_40_bacteria
)

anova_40 <- summary(
  model_40
)[[1]]

F_40 <- anova_40[
  "strain",
  "F value"
]

p_40 <- anova_40[
  "strain",
  "Pr(>F)"
]

tukey_40 <- stats::TukeyHSD(
  model_40,
  "strain"
)

tukey_matrix <- tukey_40[
  ["strain"]
]

p_values <- tukey_matrix[
  ,
  "p adj"
]

names(p_values) <- rownames(
  tukey_matrix
)

raw_letters <- multcompView::multcompLetters(
  p_values,
  threshold = 0.05
)$Letters

summary_40_bacteria <- auc_40_bacteria %>%
  dplyr::mutate(
    strain = as.character(
      strain
    )
  ) %>%
  dplyr::group_by(
    strain
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    Mean =
      mean(Relative_AUC),
    SD =
      stats::sd(Relative_AUC),
    Minimum =
      min(Relative_AUC),
    Maximum =
      max(Relative_AUC),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean)
  )

display_letters <- normalise_letters(
  raw_letters,
  summary_40_bacteria$strain
)

letters_table <- tibble::tibble(
  strain =
    names(display_letters),
  Tukey_group =
    unname(display_letters)
)

summary_40_bacteria <- summary_40_bacteria %>%
  dplyr::left_join(
    letters_table,
    by = "strain"
  )

summary_40_yeast <- auc_data %>%
  dplyr::filter(
    strain == yeast_strain,
    peg == 40
  ) %>%
  dplyr::summarise(
    strain = yeast_strain,
    n = dplyr::n(),
    Mean =
      mean(Relative_AUC),
    SD =
      stats::sd(Relative_AUC),
    Minimum =
      min(Relative_AUC),
    Maximum =
      max(Relative_AUC),
    Tukey_group = "—"
  )

species_names <- c(
  "1013" = "Pseudomonas putida 1013",
  "1014" = "Stenotrophomonas maltophilia 1014",
  "1015" = "Erwinia billingiae 1015",
  "1016" = "Erwinia billingiae 1016",
  "1017" = "Pantoea agglomerans 1017",
  "1018" = "Erwinia billingiae 1018",
  "1090" = "Stenotrophomonas maltophilia 1090",
  "1107" = "Aureobasidium pullulans 1107"
)

table_s10_numeric <- dplyr::bind_rows(
  summary_40_yeast,
  summary_40_bacteria
) %>%
  dplyr::mutate(
    Species =
      unname(
        species_names[strain]
      )
  ) %>%
  dplyr::select(
    strain,
    Species,
    n,
    Mean,
    SD,
    Minimum,
    Maximum,
    Tukey_group
  )

table_s10_publication <- table_s10_numeric %>%
  dplyr::mutate(
    `Relative AUC (% of PEG 0 control)` =
      sprintf(
        "%.2f ± %.2f",
        Mean,
        SD
      )
  ) %>%
  dplyr::select(
    Strain = Species,
    n,
    `Relative AUC (% of PEG 0 control)`,
    `Tukey group` =
      Tukey_group
  )

tukey_pairs <- as.data.frame(
  tukey_matrix
) %>%
  tibble::rownames_to_column(
    "Comparison"
  )


# -------------------------
# 9) Model diagnostics
# -------------------------

shapiro_split <- stats::shapiro.test(
  residual_data$within_residual
)

shapiro_40 <- stats::shapiro.test(
  stats::residuals(model_40)
)

bartlett_40 <- stats::bartlett.test(
  Relative_AUC ~ strain,
  data = auc_40_bacteria
)

diagnostics <- tibble::tibble(
  Check = c(
    "Split-plot within-stratum residual normality",
    "PEG 40% ANOVA residual normality",
    "PEG 40% homogeneity of variances"
  ),
  Test = c(
    "Shapiro-Wilk",
    "Shapiro-Wilk",
    "Bartlett"
  ),
  Statistic = c(
    unname(
      shapiro_split$statistic
    ),
    unname(
      shapiro_40$statistic
    ),
    unname(
      bartlett_40$statistic
    )
  ),
  `p value` = c(
    shapiro_split$p.value,
    shapiro_40$p.value,
    bartlett_40$p.value
  )
)


# -------------------------
# 10) Font configuration
# -------------------------

font_family <- "serif"

times_regular <- systemfonts::match_font(
  "Times New Roman",
  bold = FALSE,
  italic = FALSE
)$path

times_bold <- systemfonts::match_font(
  "Times New Roman",
  bold = TRUE,
  italic = FALSE
)$path

times_italic <- systemfonts::match_font(
  "Times New Roman",
  bold = FALSE,
  italic = TRUE
)$path

times_bolditalic <- systemfonts::match_font(
  "Times New Roman",
  bold = TRUE,
  italic = TRUE
)$path

font_paths <- c(
  times_regular,
  times_bold,
  times_italic,
  times_bolditalic
)

if (
  all(!is.na(font_paths)) &&
  all(file.exists(font_paths))
) {

  sysfonts::font_add(
    family = "TimesNRPEG",
    regular = times_regular,
    bold = times_bold,
    italic = times_italic,
    bolditalic = times_bolditalic
  )

  showtext::showtext_auto(
    enable = TRUE
  )

  font_family <- "TimesNRPEG"
}

publication_theme <- function(
  base_size = 12
) {

  ggplot2::theme_bw(
    base_family = font_family,
    base_size = base_size
  ) +
    ggplot2::theme(
      panel.grid.minor =
        ggplot2::element_blank(),
      panel.grid.major =
        ggplot2::element_line(
          linewidth = 0.30,
          colour = "grey88"
        ),
      axis.title =
        ggplot2::element_text(
          size = base_size + 1,
          colour = "black"
        ),
      axis.text =
        ggplot2::element_text(
          size = base_size - 1,
          colour = "black"
        ),
      strip.text =
        ggplot2::element_text(
          size = base_size - 1,
          colour = "black"
        ),
      plot.title =
        ggplot2::element_text(
          hjust = 0.5,
          size = base_size + 1,
          face = "plain",
          colour = "black"
        ),
      plot.margin =
        ggplot2::margin(
          8,
          8,
          8,
          8
        )
    )
}

save_publication_plot <- function(
  plot_object,
  filename,
  width,
  height
) {

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        filename,
        ".pdf"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        filename,
        ".tiff"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw"
  )
}


# -------------------------
# 11) Figure 3
# -------------------------

strain_labels_plotmath <- c(
  "1013" =
    "italic(Ps.~putida)~'1013'",
  "1014" =
    "italic(S.~maltophilia)~'1014'",
  "1015" =
    "italic(E.~billingiae)~'1015'",
  "1016" =
    "italic(E.~billingiae)~'1016'",
  "1017" =
    "italic(Pa.~agglomerans)~'1017'",
  "1018" =
    "italic(E.~billingiae)~'1018'",
  "1090" =
    "italic(S.~maltophilia)~'1090'"
)

strain_order <-
  summary_40_bacteria$strain

figure3_points <- auc_40_bacteria %>%
  dplyr::mutate(
    strain = factor(
      as.character(strain),
      levels = strain_order
    )
  )

plot_range <- diff(
  range(
    figure3_points$Relative_AUC
  )
)

figure3_summary <- summary_40_bacteria %>%
  dplyr::mutate(
    strain = factor(
      strain,
      levels = strain_order
    ),
    label_y =
      Mean +
      SD +
      0.07 *
      plot_range
  )

figure3 <- ggplot2::ggplot(
  figure3_points,
  ggplot2::aes(
    x = strain,
    y = Relative_AUC
  )
) +
  ggplot2::geom_point(
    position =
      ggplot2::position_jitter(
        width = 0.07,
        height = 0,
        seed = 123
      ),
    size = 2.4
  ) +
  ggplot2::geom_errorbar(
    data = figure3_summary,
    ggplot2::aes(
      x = strain,
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    data = figure3_summary,
    ggplot2::aes(
      x = strain,
      y = Mean
    ),
    inherit.aes = FALSE,
    size = 3.3
  ) +
  ggplot2::geom_text(
    data = figure3_summary,
    ggplot2::aes(
      x = strain,
      y = label_y,
      label = Tukey_group
    ),
    inherit.aes = FALSE,
    family = font_family,
    size = 4.2
  ) +
  ggplot2::scale_x_discrete(
    labels = function(x) {
      parse(
        text = unname(
          strain_labels_plotmath[x]
        )
      )
    }
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.05,
        0.14
      )
    )
  ) +
  ggplot2::labs(
    title =
      "Strain-specific growth retention at 40% PEG 6000",
    x =
      "Strain",
    y =
      "Relative AUC (% of PEG 0 control)"
  ) +
  publication_theme(
    13
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 55,
        hjust = 1
      ),
    panel.grid.major.x =
      ggplot2::element_blank()
  )

save_publication_plot(
  figure3,
  "Figure_3",
  width = 7.4,
  height = 6.1
)


# -------------------------
# 12) Figure S4
# -------------------------

facet_labels <- c(
  "1013" =
    "italic(Ps.~putida)~'1013'",
  "1014" =
    "italic(S.~maltophilia)~'1014'",
  "1015" =
    "italic(E.~billingiae)~'1015'",
  "1016" =
    "italic(E.~billingiae)~'1016'",
  "1017" =
    "italic(Pa.~agglomerans)~'1017'",
  "1018" =
    "italic(E.~billingiae)~'1018'",
  "1090" =
    "italic(S.~maltophilia)~'1090'",
  "1107" =
    "italic(A.~pullulans)~'1107'"
)

figure_s4_data <- microbial_long %>%
  dplyr::group_by(
    strain,
    peg,
    time
  ) %>%
  dplyr::summarise(
    Mean_corrected_OD =
      mean(OD600_corrected),
    SD =
      stats::sd(OD600_corrected),
    SE =
      SD /
      sqrt(dplyr::n()),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    strain = factor(
      strain,
      levels = microbial_strains
    ),
    PEG = factor(
      peg,
      levels = expected_pegs
    )
  )

figure_s4 <- ggplot2::ggplot(
  figure_s4_data,
  ggplot2::aes(
    x = time,
    y = Mean_corrected_OD,
    group = PEG,
    linetype = PEG
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin =
        Mean_corrected_OD -
        SE,
      ymax =
        Mean_corrected_OD +
        SE
    ),
    width = 0.28,
    linewidth = 0.35
  ) +
  ggplot2::geom_line(
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    size = 1.6
  ) +
  ggplot2::facet_wrap(
    ~ strain,
    ncol = 3,
    labeller = ggplot2::as_labeller(
      facet_labels,
      ggplot2::label_parsed
    ),
    scales = "free_y"
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "0" = "solid",
      "5" = "longdash",
      "10" = "dashed",
      "20" = "dotdash",
      "40" = "dotted"
    )
  ) +
  ggplot2::labs(
    title =
      "Growth kinetics under PEG 6000-amended conditions",
    x =
      "Time (h)",
    y =
      expression(
        "Background-corrected " *
          OD[600]
      ),
    linetype =
      "PEG 6000 (%)"
  ) +
  publication_theme(
    11
  ) +
  ggplot2::theme(
    legend.position =
      "right"
  )

save_publication_plot(
  figure_s4,
  "Figure_S4",
  width = 8.2,
  height = 6.8
)


# -------------------------
# 13) Export statistical outputs
# -------------------------

readr::write_csv(
  microbial_long %>%
    dplyr::select(
      strain,
      medium,
      peg,
      time,
      biological_replicate,
      OD600_corrected
    ),
  file.path(
    output_dir,
    "Figure_3_background_corrected_OD.csv"
  )
)

readr::write_csv(
  auc_data,
  file.path(
    output_dir,
    "Figure_3_AUC_by_biological_culture.csv"
  )
)

readr::write_csv(
  table_s9,
  file.path(
    output_dir,
    "Figure_3_Table_S9_ANOVA.csv"
  )
)

readr::write_csv(
  table_s10_publication,
  file.path(
    output_dir,
    "Figure_3_Table_S10_PEG40_Tukey.csv"
  )
)

readr::write_csv(
  table_s10_numeric,
  file.path(
    output_dir,
    "Figure_3_Table_S10_numeric.csv"
  )
)

readr::write_csv(
  tukey_pairs,
  file.path(
    output_dir,
    "Figure_3_Tukey_pairwise.csv"
  )
)

readr::write_csv(
  diagnostics,
  file.path(
    output_dir,
    "Figure_3_model_diagnostics.csv"
  )
)

readr::write_csv(
  figure_s4_data %>%
    dplyr::mutate(
      strain =
        as.character(strain),
      PEG =
        as.character(PEG)
    ),
  file.path(
    output_dir,
    "Figure_S4_data.csv"
  )
)



# -------------------------
# 14) Completion message
# -------------------------

cat(
  "\nFigure 3 analysis completed successfully.\n"
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

cat(
  "OD600 input values were already background-corrected; ",
  "no additional blank subtraction was performed.\n"
)

cat(
  "Figure 3, Figure S4, Table S9 and Table S10 were generated.\n"
)
