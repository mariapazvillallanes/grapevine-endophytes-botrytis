# ============================================================
# 03 – In vitro antagonism against Botrytis cinerea
# Clean Zenodo reproducibility script
# ============================================================

# ============================================================
# Figure 1 – In vitro antagonism against Botrytis cinerea
# VERSION 6: full font family + plotmath true italics
# Correct two-factor analysis, simple effects and publication figure
# ============================================================

# install.packages(c(
#   "tidyverse", "car", "emmeans", "multcomp", "multcompView",
#   "writexl", "ragg", "showtext", "sysfonts"
# ))

library(tidyverse)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
if (!requireNamespace("ggtext", quietly = TRUE)) stop("Missing required R package: ggtext")
library(ggtext)

# IMPORTANT: all dplyr::recode() calls are explicitly written as dplyr::recode()
# to prevent masking by functions from other loaded packages.
library(ragg)
library(showtext)
library(sysfonts)

# -------------------------
# 1) Input and output paths
# -------------------------
if (dir.exists("data")) {
  project_dir <- "."
} else if (dir.exists("../data")) {
  project_dir <- ".."
} else {
  stop("Project directory not found. Run from the project root or /scripts.")
}

input_file <- file.path(project_dir, "data", "03_in_vitro_antagonism.csv")
output_dir <- file.path(project_dir, "output", "03_in_vitro_antagonism")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

# Times New Roman with regular, bold, italic, and bold-italic faces.
# Registering only the regular face prevents scientific names from
# appearing in true italics.
font_family <- "serif"

times_regular <- systemfonts::match_font(
  "Times New Roman", bold = FALSE, italic = FALSE
)$path
times_bold <- systemfonts::match_font(
  "Times New Roman", bold = TRUE, italic = FALSE
)$path
times_italic <- systemfonts::match_font(
  "Times New Roman", bold = FALSE, italic = TRUE
)$path
times_bolditalic <- systemfonts::match_font(
  "Times New Roman", bold = TRUE, italic = TRUE
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
    family = "TimesNRFull",
    regular = times_regular,
    bold = times_bold,
    italic = times_italic,
    bolditalic = times_bolditalic
  )
  showtext::showtext_auto(enable = TRUE)
  font_family <- "TimesNRFull"
}

# -------------------------
# 2) Read and clean data
# -------------------------
raw <- readr::read_delim(
  input_file,
  delim = ";",
  locale = readr::locale(decimal_mark = ","),
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
)

strain_levels <- c(
  "1013", "1014", "1015", "1016", "1017", "1018",
  "1090", "890", "932", "950", "1107"
)

antag <- raw %>%
  dplyr::transmute(
    Strain = as.character(Strain),
    Trait = as.character(Trait),
    Replicate = as.integer(Rep),
    Inhibition = as.numeric(Value)
  ) %>%
  dplyr::filter(
    Trait %in% c(
      "Inhib_B0510", "Inhib_B05.10", "Inhib_UCA992"
    )
  ) %>%
  dplyr::mutate(
    Isolate = dplyr::case_when(
      Trait %in% c("Inhib_B0510", "Inhib_B05.10") ~ "B05.10",
      Trait == "Inhib_UCA992" ~ "UCA992",
      TRUE ~ NA_character_
    ),
    Strain = factor(Strain, levels = strain_levels),
    Isolate = factor(Isolate, levels = c("B05.10", "UCA992"))
  ) %>%
  drop_na(Inhibition)

# Three independent observations are expected for every cell.
design_check <- antag %>%
  dplyr::count(Strain, Isolate, name = "n")

if (any(design_check$n != 3L)) {
  print(design_check)
  stop(
    "Every Strain × Isolate combination must contain exactly ",
    "three independent biological replicates."
  )
}

# -------------------------
# 3) Two-way ANOVA
# -------------------------
old_contrasts <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
on.exit(options(old_contrasts), add = TRUE)

model <- lm(
  Inhibition ~ Strain * Isolate,
  data = antag
)

anova_table <- car::Anova(
  model,
  type = 3
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  dplyr::filter(Effect != "(Intercept)") %>%
  dplyr::transmute(
    Effect = dplyr::recode(
      Effect,
      "Strain" = "Antagonist strain",
      "Isolate" = "B. cinerea isolate",
      "Strain:Isolate" = "Antagonist strain × B. cinerea isolate",
      "Residuals" = "Residual"
    ),
    df = Df,
    `F value` = `F value`,
    `p value` = `Pr(>F)`
  )

# Basic diagnostic checks
diagnostics <- tibble(
  Test = c(
    "Shapiro-Wilk test of model residuals",
    "Levene test across Strain × Isolate cells"
  ),
  Statistic = c(
    unname(shapiro.test(residuals(model))$statistic),
    unname(car::leveneTest(
      Inhibition ~ interaction(Strain, Isolate),
      data = antag
    )[1, "F value"])
  ),
  `p value` = c(
    shapiro.test(residuals(model))$p.value,
    car::leveneTest(
      Inhibition ~ interaction(Strain, Isolate),
      data = antag
    )[1, "Pr(>F)"]
  )
)

# -------------------------
# 4) Post hoc comparisons
# -------------------------
# A) Antagonists compared within each B. cinerea isolate
emm_strain <- emmeans(
  model,
  ~ Strain | Isolate
)

tukey_pairs <- pairs(
  emm_strain,
  adjust = "tukey"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  as_tibble()

letters_table <- multcomp::cld(
  emm_strain,
  adjust = "tukey",
  Letters = letters,
  sort = FALSE
) %>%
  as.data.frame() %>%
  as_tibble() %>%
  dplyr::transmute(
    Strain = as.character(Strain),
    Isolate = as.character(Isolate),
    Tukey_group = stringr::str_remove_all(.group, "\\s+")
  )

# B) B05.10 versus UCA992 within each antagonist
# The raw p-values are adjusted jointly across the 11 strain-specific tests.
emm_isolate <- emmeans(
  model,
  ~ Isolate | Strain
)

isolate_contrasts <- contrast(
  emm_isolate,
  method = "revpairwise",
  adjust = "none"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  as_tibble() %>%
  dplyr::mutate(
    `Holm-adjusted p value` = p.adjust(p.value, method = "holm"),
    Significance = dplyr::case_when(
      `Holm-adjusted p value` < 0.001 ~ "***",
      `Holm-adjusted p value` < 0.01 ~ "**",
      `Holm-adjusted p value` < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# -------------------------
# 5) Summary tables
# -------------------------
species_names <- c(
  "1013" = "Pseudomonas putida 1013",
  "1014" = "Stenotrophomonas maltophilia 1014",
  "1015" = "Erwinia billingiae 1015",
  "1016" = "Erwinia billingiae 1016",
  "1017" = "Pantoea agglomerans 1017",
  "1018" = "Erwinia billingiae 1018",
  "1090" = "Stenotrophomonas maltophilia 1090",
  "890"  = "Penicillium chrysogenum 890",
  "932"  = "Nothophoma quercina 932",
  "950"  = "Fusarium avenaceum 950",
  "1107" = "Aureobasidium pullulans 1107"
)

summary_table <- antag %>%
  group_by(Strain, Isolate) %>%
  dplyr::summarise(
    Mean = mean(Inhibition),
    SD = sd(Inhibition),
    n = n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Strain_code = as.character(Strain)) %>%
  dplyr::left_join(
    letters_table,
    by = c("Strain_code" = "Strain", "Isolate")
  ) %>%
  dplyr::mutate(
    Microorganism = unname(species_names[Strain_code])
  )

table_s4 <- summary_table %>%
  dplyr::select(Microorganism, Isolate, n, Mean, SD) %>%
  pivot_wider(
    names_from = Isolate,
    values_from = c(Mean, SD)
  )

table_s6 <- summary_table %>%
  dplyr::select(Strain_code, Microorganism, Isolate, Tukey_group) %>%
  pivot_wider(
    names_from = Isolate,
    values_from = Tukey_group,
    names_prefix = "Tukey group – "
  ) %>%
  dplyr::left_join(
    isolate_contrasts %>%
      dplyr::transmute(
        Strain_code = as.character(Strain),
        `UCA992 − B05.10 estimate` = estimate,
        `Holm-adjusted p value` = `Holm-adjusted p value`,
        Significance
      ),
    by = "Strain_code"
  ) %>%
  dplyr::arrange(factor(Strain_code, levels = strain_levels)) %>%
  dplyr::select(-Strain_code)

# -------------------------
# 6) Publication figure
# -------------------------

# Order microorganisms according to mean inhibition against B05.10.
plot_order <- summary_table %>%
  dplyr::filter(Isolate == "B05.10") %>%
  dplyr::arrange(Mean) %>%
  dplyr::pull(Strain_code)

# Plotmath labels: scientific names italic and strain numbers upright.
species_labels_plotmath <- c(
  "1013" = "italic(P.~putida)~'1013'",
  "1014" = "italic(S.~maltophilia)~'1014'",
  "1015" = "italic(E.~billingiae)~'1015'",
  "1016" = "italic(E.~billingiae)~'1016'",
  "1017" = "italic(P.~agglomerans)~'1017'",
  "1018" = "italic(E.~billingiae)~'1018'",
  "1090" = "italic(S.~maltophilia)~'1090'",
  "890"  = "italic(P.~chrysogenum)~'890'",
  "932"  = "italic(N.~quercina)~'932'",
  "950"  = "italic(F.~avenaceum)~'950'",
  "1107" = "italic(A.~pullulans)~'1107'"
)

plot_df <- summary_table %>%
  dplyr::mutate(
    Strain_code = factor(Strain_code, levels = plot_order),
    Isolate = factor(Isolate, levels = c("B05.10", "UCA992")),
    label_x = Mean + SD + 2.2
  )

p <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = Mean,
    y = Strain_code,
    fill = Isolate
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.25
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = Mean - SD,
      xmax = Mean + SD
    ),
    width = 0.18,
    linewidth = 0.55
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = label_x,
      label = Tukey_group
    ),
    hjust = 0,
    size = 4.0,
    family = font_family
  ) +
  ggplot2::facet_wrap(
    ~ Isolate,
    nrow = 1
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "B05.10" = "grey42",
      "UCA992" = "grey72"
    ),
    guide = "none"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      parse(
        text = unname(
          species_labels_plotmath[x]
        )
      )
    }
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 85),
    breaks = seq(0, 80, by = 20),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    x = "Inhibition (%)",
    y = NULL
  ) +
  ggplot2::theme_bw(
    base_family = font_family,
    base_size = 11
  ) +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(
      size = 9
    ),
    panel.spacing = grid::unit(1.0, "cm")
  )
    
    
    
ggplot2::ggsave(
  file.path(output_dir, "Figure_1_antagonism_FINAL_publication_v6.pdf"),
  p,
  device = cairo_pdf,
  width = 190,
  height = 120,
  units = "mm"
)

ggplot2::ggsave(
  file.path(output_dir, "Figure_1_antagonism_FINAL_publication_v6.tiff"),
  p,
  device = ragg::agg_tiff,
  width = 190,
  height = 120,
  units = "mm",
  res = 600,
  compression = "lzw"
)

# -------------------------
# 7) Export results
# -------------------------
readr::write_csv(
  anova_table,
  file.path(output_dir, "Table_S5_ANOVA_antagonism.csv")
)
readr::write_csv(
  table_s4,
  file.path(output_dir, "Table_S4_antagonism_summary.csv")
)
readr::write_csv(
  table_s6,
  file.path(output_dir, "Table_S6_posthoc_antagonism.csv")
)
readr::write_csv(
  tukey_pairs,
  file.path(output_dir, "Antagonists_within_isolate_Tukey.csv")
)
readr::write_csv(
  isolate_contrasts,
  file.path(output_dir, "Isolates_within_antagonist_Holm.csv")
)
readr::write_csv(
  diagnostics,
  file.path(output_dir, "Model_diagnostics.csv")
)


cat("\nCompleted successfully.\n")
cat("Output folder:", output_dir, "\n")
cat(
  "Important interpretation: because the interaction is significant, ",
  "strain-specific and isolate-specific contrasts should be used.\n"
)
