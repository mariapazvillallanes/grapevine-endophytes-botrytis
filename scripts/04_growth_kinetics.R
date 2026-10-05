# ============================================================
# 04 – Growth kinetics
# Clean Zenodo reproducibility script
# ============================================================

# ============================================================
# Figure 2 – Growth kinetics and derived kinetic parameters
# Reproducible analysis for Figure 2, Figure S2 and associated
# statistical outputs
# ============================================================

# Experimental structure:
#   8 strains
#   3 independent biological cultures per strain
#   3 technical wells per biological culture
#
# Statistical unit:
#   biological culture (n = 3), not the individual technical well.
#
# Inferential comparisons:
#   The seven bacterial strains cultured in TSB are compared using
#   Kruskal-Wallis tests followed by Dunn tests with BH correction.
#   Aureobasidium pullulans 1107, cultured in YPDB, is retained
#   descriptively but excluded from inferential comparisons.


# -------------------------
# 1) Required packages
# -------------------------

required_packages <- c(
  "tidyverse",
  "zoo",
  "pracma",
  "rstatix",
  "multcompView",
  "ragg",
  "sysfonts",
  "showtext",
  "systemfonts",
  "ggtext"
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
  "04_growth_kinetics.csv"
)

output_dir <- file.path(project_dir, "output", "04_growth_kinetics")

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
# 3) Strain definitions
# -------------------------

yeast_strains <- c("1107")

bacterial_strains <- c(
  "1013", "1014", "1015", "1016",
  "1017", "1018", "1090"
)

strain_levels <- c(
  "1013", "1014", "1015", "1016",
  "1017", "1018", "1090", "1107"
)

strain_species <- c(
  "1013" = "Ps. putida (1013)",
  "1014" = "S. maltophilia (1014)",
  "1015" = "E. billingiae (1015)",
  "1016" = "E. billingiae (1016)",
  "1017" = "Pa. agglomerans (1017)",
  "1018" = "E. billingiae (1018)",
  "1090" = "S. maltophilia (1090)",
  "1107" = "A. pullulans (1107)"
)

strain_colors <- c(
  "1013" = "#E69F00",
  "1014" = "#56B4E9",
  "1015" = "#009E73",
  "1016" = "#999933",
  "1017" = "#0072B2",
  "1018" = "#D55E00",
  "1090" = "#CC79A7",
  "1107" = "#000000"
)

strain_labels_plotmath <- c(
  "1013" = "italic(P.~putida)~'(1013)'",
  "1014" = "italic(S.~maltophilia)~'(1014)'",
  "1015" = "italic(E.~billingiae)~'(1015)'",
  "1016" = "italic(E.~billingiae)~'(1016)'",
  "1017" = "italic(P.~agglomerans)~'(1017)'",
  "1018" = "italic(E.~billingiae)~'(1018)'",
  "1090" = "italic(S.~maltophilia)~'(1090)'",
  "1107" = "italic(A.~pullulans)~'(1107)'"
)


# -------------------------
# 4) Plot helpers
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
    family = "TimesNRFullGrowth",
    regular = times_regular,
    bold = times_bold,
    italic = times_italic,
    bolditalic = times_bolditalic
  )

  showtext::showtext_auto(enable = TRUE)
  font_family <- "TimesNRFullGrowth"
}

publication_theme <- function(base_size = 13) {

  ggplot2::theme_classic(
    base_family = font_family,
    base_size = base_size
  ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = base_size + 2,
        colour = "black"
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 1,
        colour = "black"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.7,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.7,
        colour = "black"
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = base_size + 1,
        face = "bold",
        colour = "black"
      ),
      legend.title = ggplot2::element_text(
        size = base_size,
        colour = "black"
      ),
      legend.text = ggplot2::element_text(
        size = base_size - 2,
        colour = "black"
      ),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

save_pub_plot <- function(
  plot,
  filename,
  width,
  height
) {

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(filename, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(filename, ".tiff")
    ),
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw"
  )
}


# -------------------------
# 5) Read input data
# -------------------------

raw <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(decimal_mark = ","),
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
)

required_columns <- c(
  "Type",
  "Strain_or_control",
  "Biological_replicate",
  "Technical_replicate",
  "Time_h",
  "Medium",
  "OD600_raw",
  "OD600_corrected"
)

missing_columns <- setdiff(
  required_columns,
  names(raw)
)

if (length(missing_columns) > 0) {
  stop(
    "The input file is missing the following required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# -------------------------
# 6) Blank correction
# -------------------------

blank_data <- raw %>%
  dplyr::filter(Type == "Control") %>%
  dplyr::transmute(
    Medium = as.character(Medium),
    Time_h = as.numeric(Time_h),
    Blank_OD = as.numeric(OD600_raw)
  ) %>%
  dplyr::group_by(
    Medium,
    Time_h
  ) %>%
  dplyr::summarise(
    Blank_OD = mean(
      Blank_OD,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

sample_data <- raw %>%
  dplyr::filter(Type == "Sample") %>%
  dplyr::transmute(
    Strain = as.character(Strain_or_control),
    Bio_rep = as.integer(Biological_replicate),
    Tech_rep = as.integer(Technical_replicate),
    Time_h = as.numeric(Time_h),
    Medium = as.character(Medium),
    OD_raw = as.numeric(OD600_raw),
    OD_corrected_input = as.numeric(OD600_corrected)
  )

technical_corrected <- sample_data %>%
  dplyr::left_join(
    blank_data,
    by = c(
      "Medium",
      "Time_h"
    )
  ) %>%
  dplyr::mutate(
    Strain = factor(
      Strain,
      levels = strain_levels
    ),
    OD_corrected_recalculated = OD_raw - Blank_OD,
    correction_difference =
      OD_corrected_recalculated -
      OD_corrected_input,
    OD_corrected = OD_corrected_recalculated
  )

if (any(is.na(technical_corrected$Strain))) {

  unknown_strains <- sample_data %>%
    dplyr::filter(
      !Strain %in% strain_levels
    ) %>%
    dplyr::distinct(Strain) %>%
    dplyr::pull(Strain)

  stop(
    "Unexpected strain identifiers found in the input file: ",
    paste(unknown_strains, collapse = ", ")
  )
}

if (any(is.na(technical_corrected$Blank_OD))) {
  stop(
    "At least one sample row could not be matched to a ",
    "blank control with the same medium and time point."
  )
}

max_correction_difference <- max(
  abs(
    technical_corrected$correction_difference
  ),
  na.rm = TRUE
)

if (
  is.finite(max_correction_difference) &&
  max_correction_difference > 1e-4
) {

  warning(
    "The recalculated blank-corrected OD differs from the ",
    "OD600_corrected input column by as much as ",
    signif(max_correction_difference, 4),
    ". The analysis uses the recalculated values."
  )

} else {

  message(
    "Blank-correction check passed: recalculated values ",
    "agree with the input values within rounding tolerance."
  )
}


# -------------------------
# 7) Experimental design checks
# -------------------------

bio_check <- technical_corrected %>%
  dplyr::distinct(
    Strain,
    Bio_rep
  ) %>%
  dplyr::count(
    Strain,
    name = "n_biological_replicates"
  )

tech_check <- technical_corrected %>%
  dplyr::count(
    Strain,
    Bio_rep,
    Time_h,
    name = "n_technical_wells"
  )

if (
  any(
    bio_check$n_biological_replicates != 3L
  )
) {

  print(bio_check)

  stop(
    "Every strain must contain exactly three independent ",
    "biological cultures."
  )
}

if (
  any(
    tech_check$n_technical_wells != 3L
  )
) {

  print(
    tech_check %>%
      dplyr::filter(
        n_technical_wells != 3L
      )
  )

  stop(
    "Every biological culture and time point must contain ",
    "exactly three technical wells."
  )
}

observed_times <- sort(
  unique(
    technical_corrected$Time_h
  )
)

if (
  !identical(
    observed_times,
    as.numeric(0:36)
  )
) {

  warning(
    "The input file does not contain exactly the expected ",
    "hourly time points from 0 to 36 h. The script will ",
    "continue using the observed Time_h values."
  )
}

message(
  "Design check passed: every strain contains 3 biological ",
  "cultures × 3 technical wells."
)


# -------------------------
# 8) Biological growth curves
# -------------------------

# Technical wells are averaged before growth descriptors
# are calculated.

biological_curves <- technical_corrected %>%
  dplyr::group_by(
    Strain,
    Bio_rep,
    Medium,
    Time_h
  ) %>%
  dplyr::summarise(
    OD_corrected = mean(
      OD_corrected,
      na.rm = TRUE
    ),
    Technical_SD = sd(
      OD_corrected,
      na.rm = TRUE
    ),
    n_technical = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Biological_curve_ID = paste(
      Strain,
      Bio_rep,
      sep = "_bio"
    )
  )

if (
  dplyr::n_distinct(
    biological_curves$Biological_curve_ID
  ) != 24L
) {

  stop(
    "Expected 24 biological curves ",
    "(8 strains × 3 biological replicates)."
  )
}


# -------------------------
# 9) Growth descriptors
# -------------------------

calc_growth_params <- function(
  df,
  window_size = 4,
  baseline_points = 3
) {

  df <- df %>%
    dplyr::arrange(Time_h) %>%
    dplyr::mutate(
      OD_positive = dplyr::if_else(
        OD_corrected > 0,
        OD_corrected,
        NA_real_
      ),
      lnOD = log(OD_positive)
    )

  Ymax <- max(
    df$OD_corrected,
    na.rm = TRUE
  )

  AUC <- pracma::trapz(
    df$Time_h,
    df$OD_corrected
  )

  possible_windows <- nrow(df) - window_size + 1

  if (possible_windows < 1) {
    return(
      tibble::tibble(
        Ymax = Ymax,
        AUC = AUC,
        mu_max = NA_real_,
        lag_h = NA_real_,
        exp_phase_start_h = NA_real_,
        exp_phase_end_h = NA_real_,
        exp_phase_r2 = NA_real_
      )
    )
  }

  window_results <- purrr::map_dfr(
    seq_len(possible_windows),
    function(i) {

      sub <- df[
        i:(i + window_size - 1),
      ]

      valid <- sub %>%
        dplyr::filter(
          !is.na(lnOD)
        )

      if (nrow(valid) < 3) {
        return(
          tibble::tibble(
            i = i,
            slope = NA_real_,
            r2 = NA_real_
          )
        )
      }

      fit <- stats::lm(
        lnOD ~ Time_h,
        data = valid
      )

      tibble::tibble(
        i = i,
        slope = unname(
          stats::coef(fit)[2]
        ),
        r2 = summary(fit)$r.squared
      )
    }
  ) %>%
    dplyr::filter(
      !is.na(slope),
      slope > 0
    )

  if (nrow(window_results) == 0) {
    return(
      tibble::tibble(
        Ymax = Ymax,
        AUC = AUC,
        mu_max = NA_real_,
        lag_h = NA_real_,
        exp_phase_start_h = NA_real_,
        exp_phase_end_h = NA_real_,
        exp_phase_r2 = NA_real_
      )
    )
  }

  best <- window_results %>%
    dplyr::arrange(
      dplyr::desc(slope),
      dplyr::desc(r2)
    ) %>%
    dplyr::slice(1)

  best_i <- best$i

  exp_phase <- df[
    best_i:(best_i + window_size - 1),
  ] %>%
    dplyr::filter(
      !is.na(lnOD)
    )

  fit <- stats::lm(
    lnOD ~ Time_h,
    data = exp_phase
  )

  intercept <- unname(
    stats::coef(fit)[1]
  )

  slope <- unname(
    stats::coef(fit)[2]
  )

  baseline_values <- df$OD_positive[
    seq_len(
      min(
        baseline_points,
        nrow(df)
      )
    )
  ]

  OD0_used <- mean(
    baseline_values,
    na.rm = TRUE
  )

  lag_h <- if (
    is.na(OD0_used) ||
    is.na(slope) ||
    slope <= 0
  ) {

    NA_real_

  } else {

    (
      log(OD0_used) -
      intercept
    ) / slope
  }

  if (
    !is.na(lag_h) &&
    lag_h < 0
  ) {
    lag_h <- 0
  }

  tibble::tibble(
    Ymax = Ymax,
    AUC = AUC,
    mu_max = slope,
    lag_h = lag_h,
    exp_phase_start_h = min(
      exp_phase$Time_h
    ),
    exp_phase_end_h = max(
      exp_phase$Time_h
    ),
    exp_phase_r2 = summary(fit)$r.squared
  )
}

growth_results <- biological_curves %>%
  dplyr::group_by(
    Strain,
    Bio_rep,
    Medium,
    Biological_curve_ID
  ) %>%
  dplyr::group_modify(
    ~ calc_growth_params(
      .x,
      window_size = 4,
      baseline_points = 3
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Species = unname(
      strain_species[
        as.character(Strain)
      ]
    )
  )

n_check <- growth_results %>%
  dplyr::count(
    Strain,
    name = "n"
  )

if (any(n_check$n != 3L)) {

  print(n_check)

  stop(
    "Kinetic descriptors must contain n = 3 biological ",
    "values per strain."
  )
}

growth_results_inference <- growth_results %>%
  dplyr::filter(
    as.character(Strain) %in%
      bacterial_strains
  ) %>%
  dplyr::mutate(
    Strain = droplevels(Strain)
  )


# -------------------------
# 10) Descriptive summaries
# -------------------------

summary_results <- growth_results %>%
  dplyr::group_by(
    Strain,
    Species,
    Medium
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    Ymax_mean = mean(
      Ymax,
      na.rm = TRUE
    ),
    Ymax_sd = stats::sd(
      Ymax,
      na.rm = TRUE
    ),
    AUC_mean = mean(
      AUC,
      na.rm = TRUE
    ),
    AUC_sd = stats::sd(
      AUC,
      na.rm = TRUE
    ),
    mu_max_mean = mean(
      mu_max,
      na.rm = TRUE
    ),
    mu_max_sd = stats::sd(
      mu_max,
      na.rm = TRUE
    ),
    lag_mean = mean(
      lag_h,
      na.rm = TRUE
    ),
    lag_sd = stats::sd(
      lag_h,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# -------------------------
# 11) Kruskal-Wallis and Dunn-BH tests
# -------------------------

run_kw_dunn <- function(
  data,
  parameter_name
) {

  df <- data %>%
    dplyr::select(
      Strain,
      dplyr::all_of(parameter_name)
    ) %>%
    dplyr::rename(
      Value =
        dplyr::all_of(parameter_name)
    ) %>%
    dplyr::filter(
      !is.na(Value)
    )

  kw <- rstatix::kruskal_test(
    df,
    Value ~ Strain
  ) %>%
    dplyr::mutate(
      Parameter = parameter_name
    ) %>%
    dplyr::select(
      Parameter,
      dplyr::everything()
    )

  dunn <- rstatix::dunn_test(
    df,
    Value ~ Strain,
    p.adjust.method = "BH"
  ) %>%
    dplyr::mutate(
      Parameter = parameter_name
    ) %>%
    dplyr::select(
      Parameter,
      dplyr::everything()
    )

  list(
    kw = kw,
    dunn = dunn
  )
}

parameters <- c(
  "mu_max",
  "Ymax",
  "lag_h",
  "AUC"
)

stat_results <- lapply(
  parameters,
  function(parameter_name) {
    run_kw_dunn(
      growth_results_inference,
      parameter_name
    )
  }
)

kw_results <- dplyr::bind_rows(
  lapply(
    stat_results,
    `[[`,
    "kw"
  )
)

dunn_results <- dplyr::bind_rows(
  lapply(
    stat_results,
    `[[`,
    "dunn"
  )
)

make_letters_from_dunn <- function(
  dunn_df,
  parameter_name,
  alpha = 0.05
) {

  d <- dunn_df %>%
    dplyr::filter(
      Parameter == parameter_name
    ) %>%
    dplyr::mutate(
      comparison = paste(
        group1,
        group2,
        sep = "-"
      )
    )

  p_values <- d$p.adj
  names(p_values) <- d$comparison

  letters <- multcompView::multcompLetters(
    p_values,
    threshold = alpha
  )$Letters

  tibble::tibble(
    Strain = names(letters),
    Letters = unname(letters),
    Parameter_raw = parameter_name
  )
}

letters_table <- dplyr::bind_rows(
  lapply(
    parameters,
    function(parameter_name) {
      make_letters_from_dunn(
        dunn_results,
        parameter_name
      )
    }
  )
)


# -------------------------
# 12) Figure S2 – Growth curves
# -------------------------

growth_plot_data <- biological_curves %>%
  dplyr::group_by(
    Strain,
    Time_h
  ) %>%
  dplyr::summarise(
    OD_mean = mean(
      OD_corrected,
      na.rm = TRUE
    ),
    OD_sd = stats::sd(
      OD_corrected,
      na.rm = TRUE
    ),
    n_biological = dplyr::n(),
    OD_sem =
      OD_sd /
      sqrt(n_biological),
    .groups = "drop"
  )

p_growth <- ggplot2::ggplot(
  growth_plot_data,
  ggplot2::aes(
    x = Time_h,
    y = OD_mean,
    color = Strain,
    fill = Strain,
    group = Strain
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = OD_mean - OD_sem,
      ymax = OD_mean + OD_sem
    ),
    alpha = 0.12,
    color = NA,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(
    linewidth = 1.15
  ) +
  ggplot2::scale_color_manual(
    values = strain_colors,
    labels = function(x) {
      parse(
        text = unname(
          strain_labels_plotmath[x]
        )
      )
    }
  ) +
  ggplot2::scale_fill_manual(
    values = strain_colors,
    labels = function(x) {
      parse(
        text = unname(
          strain_labels_plotmath[x]
        )
      )
    }
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(
      0,
      36,
      by = 6
    ),
    limits = c(
      0,
      36
    ),
    expand = c(
      0.01,
      0
    )
  ) +
  ggplot2::labs(
    x = "Time (h)",
    y = expression(
      "Corrected " * OD[600]
    ),
    color = "Microorganism"
  ) +
  publication_theme(13) +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(
      family = font_family,
      size = 10
    )
  )

save_pub_plot(
  p_growth,
  "Figure_S2",
  width = 8.5,
  height = 5.5
)


# -------------------------
# 13) Diagnostic growth curves
# -------------------------

# One curve per independent biological culture after
# averaging its three technical wells.

p_growth_biological <- ggplot2::ggplot(
  biological_curves,
  ggplot2::aes(
    x = Time_h,
    y = OD_corrected,
    group = Bio_rep,
    linetype = factor(Bio_rep)
  )
) +
  ggplot2::geom_line(
    linewidth = 0.85,
    colour = "black",
    alpha = 0.85
  ) +
  ggplot2::facet_wrap(
    ~ Strain,
    ncol = 2,
    labeller = ggplot2::as_labeller(
      strain_labels_plotmath,
      ggplot2::label_parsed
    )
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "1" = "solid",
      "2" = "dashed",
      "3" = "dotted"
    ),
    name = "Biological culture"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(
      0,
      36,
      by = 6
    ),
    limits = range(
      biological_curves$Time_h,
      na.rm = TRUE
    ),
    expand = c(
      0.01,
      0
    )
  ) +
  ggplot2::labs(
    x = "Time (h)",
    y = expression(
      "Corrected " * OD[600]
    )
  ) +
  publication_theme(12) +
  ggplot2::theme(
    legend.position = "bottom"
  )

save_pub_plot(
  p_growth_biological,
  "Figure_2_diagnostic_growth_curves",
  width = 8.5,
  height = 9.5
)


# -------------------------
# 14) Figure 2 – Kinetic parameters
# -------------------------

parameter_labels <- c(
  "mu_max" = "mu[max]~'(h'^-1*')'",
  "Ymax" = "Y[max]",
  "lag_h" = "lambda~'(h)'",
  "AUC" = "AUC"
)

parameter_order <- c(
  "mu_max",
  "Ymax",
  "lag_h",
  "AUC"
)

panel_tags <- tibble::tibble(
  Parameter_raw = factor(
    parameter_order,
    levels = parameter_order
  ),
  panel_tag = c(
    "A",
    "B",
    "C",
    "D"
  )
)

growth_results_bacteria <- growth_results %>%
  dplyr::filter(
    as.character(Strain) %in%
      bacterial_strains
  )

growth_params_long <- growth_results_bacteria %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      parameter_order
    ),
    names_to = "Parameter_raw",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Parameter_raw = factor(
      Parameter_raw,
      levels = parameter_order
    ),
    Parameter = factor(
      unname(
        parameter_labels[
          as.character(Parameter_raw)
        ]
      ),
      levels = unname(
        parameter_labels[
          parameter_order
        ]
      )
    ),
    Strain = factor(
      Strain,
      levels = bacterial_strains
    )
  )

panel_tags_for_plot <- panel_tags %>%
  dplyr::mutate(
    Parameter = factor(
      unname(
        parameter_labels[
          as.character(Parameter_raw)
        ]
      ),
      levels = unname(
        parameter_labels[
          parameter_order
        ]
      )
    )
  )

letters_for_plot <- letters_table %>%
  dplyr::filter(
    Strain %in% bacterial_strains
  ) %>%
  dplyr::mutate(
    Parameter = factor(
      unname(
        parameter_labels[
          Parameter_raw
        ]
      ),
      levels = unname(
        parameter_labels[
          parameter_order
        ]
      )
    ),
    Strain = factor(
      Strain,
      levels = bacterial_strains
    )
  )

letter_positions <- growth_params_long %>%
  dplyr::group_by(
    Parameter,
    Strain
  ) %>%
  dplyr::summarise(
    y_max = max(
      Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    letters_for_plot,
    by = c(
      "Parameter",
      "Strain"
    )
  ) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::mutate(
    panel_range = diff(
      range(
        growth_params_long$Value[
          growth_params_long$Parameter ==
            dplyr::first(Parameter)
        ],
        na.rm = TRUE
      )
    ),
    panel_range = dplyr::if_else(
      panel_range > 0,
      panel_range,
      1
    ),
    label_y =
      y_max +
      0.08 * panel_range
  ) %>%
  dplyr::ungroup()

p_parameters <- ggplot2::ggplot(
  growth_params_long,
  ggplot2::aes(
    x = Strain,
    y = Value,
    fill = Strain
  )
) +
  ggplot2::geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.5,
    colour = "black",
    alpha = 0.78
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(
      color = Strain
    ),
    width = 0.10,
    height = 0,
    size = 2.4,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = letter_positions,
    ggplot2::aes(
      x = Strain,
      y = label_y,
      label = Letters
    ),
    inherit.aes = FALSE,
    family = font_family,
    size = 4.2
  ) +
  ggplot2::geom_text(
    data = panel_tags_for_plot,
    ggplot2::aes(
      x = -Inf,
      y = Inf,
      label = panel_tag
    ),
    inherit.aes = FALSE,
    hjust = -0.35,
    vjust = 1.15,
    family = font_family,
    fontface = "bold",
    size = 5.0
  ) +
  ggplot2::facet_wrap(
    ~ Parameter,
    scales = "free_y",
    ncol = 2,
    labeller = ggplot2::label_parsed
  ) +
  ggplot2::scale_fill_manual(
    values = strain_colors[
      bacterial_strains
    ]
  ) +
  ggplot2::scale_color_manual(
    values = strain_colors[
      bacterial_strains
    ]
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.05,
        0.18
      )
    )
  ) +
  ggplot2::scale_x_discrete(
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "Strain",
    y = NULL
  ) +
  publication_theme(13) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

save_pub_plot(
  p_parameters,
  "Figure_2",
  width = 8.5,
  height = 6.8
)


# -------------------------
# 15) Export statistical outputs
# -------------------------

readr::write_csv(
  technical_corrected,
  file.path(
    output_dir,
    "Figure_2_technical_wells_corrected.csv"
  )
)

readr::write_csv(
  biological_curves,
  file.path(
    output_dir,
    "Figure_2_biological_growth_curves.csv"
  )
)

readr::write_csv(
  growth_results,
  file.path(
    output_dir,
    "Figure_2_growth_parameters.csv"
  )
)

readr::write_csv(
  summary_results,
  file.path(
    output_dir,
    "Figure_2_Table_S7_summary.csv"
  )
)

readr::write_csv(
  kw_results,
  file.path(
    output_dir,
    "Figure_2_Kruskal_Wallis.csv"
  )
)

readr::write_csv(
  dunn_results,
  file.path(
    output_dir,
    "Figure_2_Table_S8_Dunn_pairwise_BH.csv"
  )
)

readr::write_csv(
  letters_table,
  file.path(
    output_dir,
    "Figure_2_letters.csv"
  )
)



# -------------------------
# 16) Completion message
# -------------------------

cat(
  "\nFigure 2 analysis completed successfully.\n"
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
  "Statistical unit: 3 independent biological cultures per strain.\n"
)

cat(
  "Inferential comparisons: seven TSB-grown bacterial strains only.\n"
)

cat(
  "Aureobasidium pullulans 1107 is retained descriptively.\n"
)
