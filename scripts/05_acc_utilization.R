# ============================================================
# 05 – ACC utilization
# Clean Zenodo reproducibility script
# ============================================================

# ============================================================
# Figure S3 – ACC utilisation screening
# Reproducible analysis for Supplementary Figure S3
# ============================================================

# Input:
#   data/05_acc_utilization.csv
#
# Experimental structure:
#   8 strains
#   3 independent biological replicates
#   3 sampling days: 3, 5 and 7
#   3 conditions: ACC-, ACC+ and (NH4)2SO4
#
# The figure is descriptive. No statistical significance
# is inferred or displayed.


# -------------------------
# 1) Required packages
# -------------------------

required_packages <- c(
  "readr",
  "ggplot2",
  "dplyr",
  "ragg",
  "systemfonts",
  "png"
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
  "05_acc_utilization.csv"
)

output_dir <- file.path(project_dir, "output", "05_acc_utilization")

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

acc <- readr::read_delim(
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
  "Strain",
  "Replicate",
  "Time",
  "Condition",
  "Value"
)

missing_columns <- setdiff(
  required_columns,
  names(acc)
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

acc <- acc %>%
  dplyr::mutate(
    Strain =
      as.character(Strain),
    Replicate =
      as.character(Replicate),
    Time =
      as.numeric(Time),
    Condition =
      as.character(Condition),
    Value =
      as.numeric(Value)
  )

strains <- c(
  "1013",
  "1014",
  "1015",
  "1016",
  "1017",
  "1018",
  "1090",
  "1107"
)

conditions <- c(
  "ACC-",
  "ACC+",
  "(NH4)2SO4"
)

times <- c(
  3,
  5,
  7
)

replicates <- c(
  "R1",
  "R2",
  "R3"
)

expected <- expand.grid(
  Strain = strains,
  Replicate = replicates,
  Time = times,
  Condition = conditions,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

row_id <- function(d) {
  paste(
    d$Strain,
    d$Replicate,
    d$Time,
    d$Condition,
    sep = "|"
  )
}

if (
  nrow(acc) !=
    nrow(expected) ||
  anyNA(
    acc[
      ,
      required_columns
    ]
  ) ||
  any(
    !is.finite(
      acc$Value
    )
  ) ||
  anyDuplicated(
    row_id(acc)
  ) ||
  !setequal(
    row_id(acc),
    row_id(expected)
  )
) {
  stop(
    "Incomplete or unexpected dataset. ",
    "05_acc_utilization.csv should contain exactly 216 biological values."
  )
}


# -------------------------
# 4) Summary statistics
# -------------------------

acc_summary <- acc %>%
  dplyr::group_by(
    Strain,
    Time,
    Condition
  ) %>%
  dplyr::summarise(
    n =
      dplyr::n(),
    Mean_OD600 =
      mean(Value),
    SD_OD600 =
      stats::sd(Value),
    SE_OD600 =
      stats::sd(Value) /
      sqrt(
        dplyr::n()
      ),
    .groups = "drop"
  )


# -------------------------
# 5) Font configuration
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
# 6) Figure S3
# -------------------------

strain_labels <- c(
  "1013" =
    "italic('Ps. putida')~'1013'",
  "1014" =
    "italic('S. maltophilia')~'1014'",
  "1015" =
    "italic('E. billingiae')~'1015'",
  "1016" =
    "italic('E. billingiae')~'1016'",
  "1017" =
    "italic('Pa. agglomerans')~'1017'",
  "1018" =
    "italic('E. billingiae')~'1018'",
  "1090" =
    "italic('S. maltophilia')~'1090'",
  "1107" =
    "italic('A. pullulans')~'1107'"
)

acc_summary <- acc_summary %>%
  dplyr::mutate(
    Strain = factor(
      Strain,
      levels = strains
    ),
    Condition = factor(
      Condition,
      levels = conditions
    )
  )

p_s3 <- ggplot2::ggplot(
  acc_summary,
  ggplot2::aes(
    x = Time,
    y = Mean_OD600,
    group = Condition,
    linetype = Condition,
    shape = Condition
  )
) +
  ggplot2::geom_line(
    linewidth = 0.78,
    colour = "black"
  ) +
  ggplot2::geom_point(
    size = 1.8,
    stroke = 0.4,
    colour = "black"
  ) +
  ggplot2::facet_wrap(
    ~ Strain,
    ncol = 4,
    labeller =
      ggplot2::labeller(
        Strain =
          ggplot2::as_labeller(
            strain_labels,
            default =
              ggplot2::label_parsed
          )
      )
  ) +
  ggplot2::scale_x_continuous(
    breaks = times,
    limits = range(times)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(
      0.1,
      0.4,
      by = 0.1
    ),
    limits = c(
      0.08,
      0.46
    ),
    labels = function(x) {
      sprintf(
        "%.1f",
        x
      )
    }
  ) +
  ggplot2::scale_linetype_manual(
    name = "Condition",
    breaks = conditions,
    values = c(
      "ACC-" = "solid",
      "ACC+" = "dotted",
      "(NH4)2SO4" = "longdash"
    ),
    labels = c(
      "ACC−",
      "ACC+",
      "(NH₄)₂SO₄"
    )
  ) +
  ggplot2::scale_shape_manual(
    name = "Condition",
    breaks = conditions,
    values = c(
      "ACC-" = 16,
      "ACC+" = 17,
      "(NH4)2SO4" = 15
    ),
    labels = c(
      "ACC−",
      "ACC+",
      "(NH₄)₂SO₄"
    )
  ) +
  ggplot2::labs(
    x = "Incubation time (days)",
    y = expression(
      OD[600]
    )
  ) +
  ggplot2::theme_bw(
    base_family = font_family,
    base_size = 11
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_blank(),
    panel.grid =
      ggplot2::element_blank(),
    panel.border =
      ggplot2::element_rect(
        fill = NA,
        colour = "grey45",
        linewidth = 0.35
      ),
    strip.background =
      ggplot2::element_rect(
        fill = "grey95",
        colour = "grey45",
        linewidth = 0.35
      ),
    strip.text =
      ggplot2::element_text(
        size = 10.5,
        colour = "black"
      ),
    axis.text =
      ggplot2::element_text(
        size = 9.5,
        colour = "black"
      ),
    axis.title =
      ggplot2::element_text(
        size = 11.5,
        colour = "black"
      ),
    legend.position = "bottom",
    legend.title =
      ggplot2::element_text(
        size = 10.5
      ),
    legend.text =
      ggplot2::element_text(
        size = 10
      ),
    legend.key.width =
      grid::unit(
        12,
        "mm"
      ),
    panel.spacing =
      grid::unit(
        4,
        "mm"
      ),
    plot.margin =
      ggplot2::margin(
        8,
        10,
        7,
        8
      )
  )


# -------------------------
# 7) Export Figure S3
# -------------------------

tiff_file <- file.path(
  output_dir,
  "Figure_S3.tiff"
)

png_file <- file.path(
  output_dir,
  "Figure_S3.png"
)

pdf_file <- file.path(
  output_dir,
  "Figure_S3.pdf"
)

ggplot2::ggsave(
  tiff_file,
  plot = p_s3,
  device = ragg::agg_tiff,
  width = 11,
  height = 6.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggplot2::ggsave(
  png_file,
  plot = p_s3,
  device = ragg::agg_png,
  width = 11,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

# The original workflow creates the PDF from the high-resolution
# PNG to avoid PDF/PostScript font issues.

image <- png::readPNG(
  png_file,
  native = TRUE
)

grDevices::pdf(
  pdf_file,
  width = 11,
  height = 6.2,
  family = "sans",
  useDingbats = FALSE,
  paper = "special"
)

grid::grid.newpage()

grid::grid.raster(
  image,
  x = 0.5,
  y = 0.5,
  width =
    grid::unit(
      1,
      "npc"
    ),
  height =
    grid::unit(
      1,
      "npc"
    ),
  interpolate = FALSE
)

grDevices::dev.off()


# -------------------------
# 8) Export associated data
# -------------------------

readr::write_csv(
  acc_summary %>%
    dplyr::mutate(
      Strain =
        as.character(Strain),
      Condition =
        as.character(Condition)
    ),
  file.path(
    output_dir,
    "Figure_S3_summary.csv"
  )
)


# -------------------------
# 9) Completion message
# -------------------------

cat(
  "\nFigure S3 analysis completed successfully.\n"
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
