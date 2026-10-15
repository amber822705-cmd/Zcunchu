# ============================================================
# EMF_region_difference_raincloud_FINAL_v5.R
# Publication-oriented final script
# Comparison of FT vs GQ for:
# EMF-carbon, EMF-nutrient, EMF-enzyme, EMF-vegetation, EMF-avg
# Outputs:
# 1) Detailed statistics workbook
# 2) Publication-style three-line table workbook
# 3) Five TIFF raincloud plots
# 4) Editable PPTX: 4-panel functional dimensions figure
# 5) Editable PPTX: EMF-avg figure
# 6) Text summary of descriptive statistics and group comparisons
# 7) Editable Word three-line table and panel-specific figure legends
# 8) Combined four-panel TIFF and PNG previews
# ============================================================

# ---------------------------
# 0. Package management
# ---------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "ggdist", "car",
  "openxlsx", "officer", "rvg", "patchwork", "ggsignif", "stringr", "flextable"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(car)
library(openxlsx)
library(officer)
library(rvg)
library(patchwork)
library(ggsignif)
library(stringr)

# ---------------------------
# 1. Paths
# ---------------------------
input_file <- "D:/Users/Lenovo/Desktop/原始/EMFclean.xlsx"
input_sheet <- "EMF"
# Optional environment overrides allow a separate output folder for each run.
input_file <- Sys.getenv("EMF_INPUT_FILE", unset = input_file)
output_dir <- Sys.getenv(
  "EMF_OUTPUT_DIR", unset = "D:/Users/Lenovo/Desktop/原始/3.2/3.2.2/revised"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---------------------------
# 2. Helpers
# ---------------------------
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) {
    formatC(p, format = "e", digits = 2)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}

sig_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

safe_shapiro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || length(unique(x)) < 3) return(NA_real_)
  shapiro.test(x)$p.value
}

cohens_d <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  pooled_sd <- sqrt((((length(x) - 1) * sd(x)^2) + ((length(y) - 1) * sd(y)^2)) /
                      (length(x) + length(y) - 2))
  if (is.na(pooled_sd) || pooled_sd == 0) return(NA_real_)
  (mean(x) - mean(y)) / pooled_sd
}

rank_biserial <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  wt <- wilcox.test(x, y, exact = FALSE)
  U <- as.numeric(wt$statistic)
  n1 <- length(x)
  n2 <- length(y)
  # Positive ES consistently means FT > GQ, as for Cohen's d.
  (2 * U) / (n1 * n2) - 1
}

mean_ci <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  se <- s / sqrt(n)
  tval <- qt(0.975, df = n - 1)
  c(mean = m, sd = s, se = se, ci_low = m - tval * se, ci_high = m + tval * se)
}

# ---------------------------
# 3. Read and validate data
# ---------------------------
raw_dat <- read_excel(input_file, sheet = input_sheet)

required_cols <- c(
  "Plot ID", "Region",
  "EMF-carbon", "EMF-nutrient", "EMF-enzyme", "EMF-vegetation", "EMF-avg"
)

missing_cols <- setdiff(required_cols, colnames(raw_dat))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

dat <- raw_dat %>%
  select(all_of(required_cols)) %>%
  filter(!is.na(Region))

group_levels <- c("FT", "GQ")
if (!all(group_levels %in% unique(dat$Region))) {
  stop("Region column must contain both FT and GQ.")
}
dat$Region <- factor(dat$Region, levels = group_levels)

indices <- c(
  "EMF-carbon",
  "EMF-nutrient",
  "EMF-enzyme",
  "EMF-vegetation",
  "EMF-avg"
)

if (!all(vapply(dat[indices], is.numeric, logical(1)))) {
  stop("All EMF columns must be numeric.")
}
dat <- dat %>% filter(!is.na(Region))
if (any(vapply(dat[indices], function(x) any(!is.finite(x) & !is.na(x)), logical(1)))) {
  stop("EMF columns contain infinite values. Please check the input data.")
}

# ---------------------------
# 4. Statistical analysis
# ---------------------------
desc_list <- list()
assump_list <- list()
comp_list <- list()

for (ind in indices) {
  tmp <- dat %>%
    select(`Plot ID`, Region, all_of(ind)) %>%
    rename(Value = all_of(ind)) %>%
    filter(!is.na(Value))

  # descriptive stats
  desc_ind <- tmp %>%
    group_by(Region) %>%
    summarise(
      n = n(),
      Mean = mean(Value),
      SD = sd(Value),
      SE = SD / sqrt(n),
      CI_low = Mean - qt(0.975, df = n - 1) * SE,
      CI_high = Mean + qt(0.975, df = n - 1) * SE,
      Min = min(Value),
      Max = max(Value),
      .groups = "drop"
    ) %>%
    mutate(Index = ind) %>%
    relocate(Index)

  desc_list[[ind]] <- desc_ind

  ft <- tmp %>% filter(Region == "FT") %>% pull(Value)
  gq <- tmp %>% filter(Region == "GQ") %>% pull(Value)

  if (length(ft) < 2 || length(gq) < 2) {
    stop(paste(ind, "requires at least two observations per region."))
  }

  p_shapiro_ft <- safe_shapiro(ft)
  p_shapiro_gq <- safe_shapiro(gq)
  p_levene <- car::leveneTest(Value ~ Region, data = tmp)$`Pr(>F)`[1]

  normal_condition <- !is.na(p_shapiro_ft) && !is.na(p_shapiro_gq) &&
    p_shapiro_ft > 0.05 && p_shapiro_gq > 0.05
  variance_condition <- !is.na(p_levene) && p_levene > 0.05

  if (normal_condition && variance_condition) {
    test_obj <- t.test(Value ~ Region, data = tmp, var.equal = TRUE)
    method <- "Student's independent-samples t-test"
    stat_name <- "t"
    effect <- cohens_d(ft, gq)
    effect_type <- "Cohen's d"
  } else {
    test_obj <- wilcox.test(Value ~ Region, data = tmp, exact = FALSE)
    method <- "Mann-Whitney U test"
    # R's two-sample Wilcoxon W is the FT Mann-Whitney U statistic.
    stat_name <- "U"
    effect <- rank_biserial(ft, gq)
    effect_type <- "Rank-biserial correlation"
  }

  assump_list[[ind]] <- data.frame(
    Index = ind,
    Shapiro_FT_P = p_shapiro_ft,
    Shapiro_GQ_P = p_shapiro_gq,
    Levene_P = p_levene,
    Normality_met = normal_condition,
    Homogeneity_met = variance_condition
  )

  comp_list[[ind]] <- data.frame(
    Index = ind,
    Method = method,
    Statistic_name = stat_name,
    Statistic = unname(test_obj$statistic),
    P_value = test_obj$p.value,
    Significance = sig_stars(test_obj$p.value),
    Effect_size = effect,
    Effect_type = effect_type
  )
}

desc_table <- bind_rows(desc_list)
assump_table <- bind_rows(assump_list)
comp_table <- bind_rows(comp_list)

# Results-specific captions: no test labels, P values or ES inside plot panels.
caption_common <- paste(
  "Futian (FT, blue) and Gaoqiao (GQ, vermilion) mangrove forests.",
  "Points show individual observations; boxes show the median and interquartile",
  "range, with whiskers extending to observations within 1.5 times the IQR.",
  "Right half-violins show untrimmed Gaussian kernel density estimates.",
  "All comparisons are two-sided. Student's independent-samples t-test with pooled",
  "variance was used when both Shapiro-Wilk normality tests and the median-centered",
  "Levene test gave P > 0.05; otherwise, the Mann-Whitney U test was used",
  "(normal approximation with continuity correction; exact = FALSE).",
  "ES denotes Cohen's d for the t-test or rank-biserial correlation (r_rb) for the",
  "Mann-Whitney test; positive ES indicates higher values in FT.",
  "*P < 0.05, **P < 0.01, ***P < 0.001. No asterisk indicates P >= 0.05.",
  "P values are unadjusted for multiple comparisons."
)

panel_result <- function(ind, panel = NULL) {
  st <- comp_table %>% filter(Index == ind)
  ns <- desc_table %>% filter(Index == ind)
  paste0(
    if (!is.null(panel)) paste0("(", panel, ") ") else "", ind, ": ",
    st$Method, "; ES (", if (st$Effect_type == "Cohen's d") "Cohen's d" else "r_rb",
    ") = ", fmt_num(st$Effect_size), "; n(FT/GQ) = ",
    ns$n[ns$Region == "FT"], "/", ns$n[ns$Region == "GQ"], "."
  )
}

legend_functions <- paste(
  "Figure X. Regional differences in four ecosystem multifunctionality dimensions.",
  caption_common,
  paste(vapply(seq_len(4), function(i) panel_result(indices[i], LETTERS[i]), character(1)),
        collapse = " ")
)
legend_avg <- paste(
  "Figure Y. Regional difference in overall ecosystem multifunctionality.",
  caption_common, panel_result("EMF-avg")
)
writeLines(c(legend_functions, "", legend_avg),
           file.path(output_dir, "EMF_figure_legends.txt"), useBytes = TRUE)

# Publication table
pub_desc <- desc_table %>%
  mutate(`Mean ± SD` = paste0(fmt_num(Mean), " ± ", fmt_num(SD))) %>%
  select(Index, Region, `Mean ± SD`) %>%
  pivot_wider(names_from = Region, values_from = `Mean ± SD`)

pub_table <- pub_desc %>%
  left_join(comp_table, by = "Index") %>%
  transmute(
    Function = Index,
    `FT (mean ± SD)` = FT,
    `GQ (mean ± SD)` = GQ,
    `Statistical test` = Method,
    Statistic = paste0(Statistic_name, " = ", fmt_num(Statistic)),
    `P value` = vapply(P_value, fmt_p, character(1)),
    `Effect size` = paste0(fmt_num(Effect_size), " (", Effect_type, ")"),
    Significance = Significance
  )

# Narrative summary
summary_lines <- c(
  "Regional difference summary for ecosystem multifunctionality (FT vs GQ)",
  "======================================================================",
  ""
)

for (ind in indices) {
  ds <- desc_table %>% filter(Index == ind)
  cp <- comp_table %>% filter(Index == ind)
  ft_row <- ds %>% filter(Region == "FT")
  gq_row <- ds %>% filter(Region == "GQ")

  one_line <- paste0(
    ind, ": FT = ", fmt_num(ft_row$Mean), " ± ", fmt_num(ft_row$SD),
    " (n = ", ft_row$n, "); GQ = ", fmt_num(gq_row$Mean), " ± ", fmt_num(gq_row$SD),
    " (n = ", gq_row$n, "). ",
    cp$Method, ": ", cp$Statistic_name, " = ", fmt_num(cp$Statistic),
    ", P = ", fmt_p(cp$P_value), ", effect size = ", fmt_num(cp$Effect_size),
    " [", cp$Effect_type, "], significance = ", cp$Significance, "."
  )
  summary_lines <- c(summary_lines, one_line)
}

writeLines(summary_lines, con = file.path(output_dir, "EMF_region_difference_summary.txt"))

# ---------------------------
# 5. Export Excel workbooks
# ---------------------------
# 5.1 Detailed statistics workbook
wb_detail <- createWorkbook()

addWorksheet(wb_detail, "Descriptive_statistics")
writeData(wb_detail, "Descriptive_statistics", desc_table)

addWorksheet(wb_detail, "Assumption_tests")
writeData(wb_detail, "Assumption_tests", assump_table)

addWorksheet(wb_detail, "Statistical_comparison")
writeData(wb_detail, "Statistical_comparison", comp_table)

addWorksheet(wb_detail, "Publication_table")
writeData(wb_detail, "Publication_table", pub_table)
note_row <- nrow(pub_table) + 4
writeData(
  wb_detail, "Publication_table",
  x = "Table note: Values are presented as mean ± SD. Differences between Futian (FT) and Gaoqiao (GQ) were tested using independent-sample t-tests when assumptions of normality and homogeneity of variance were satisfied; otherwise, Mann–Whitney U tests were applied. Effect sizes were calculated using Cohen's d for parametric tests and rank-biserial correlation for non-parametric tests.",
  startRow = note_row, startCol = 1
)

for (sh in names(wb_detail)) {
  setColWidths(wb_detail, sh, cols = 1:20, widths = "auto")
  freezePane(wb_detail, sh, firstRow = TRUE)
  showGridLines(wb_detail, sh, showGridLines = FALSE)
}

saveWorkbook(
  wb_detail,
  file = file.path(output_dir, "EMF_region_difference_statistics.xlsx"),
  overwrite = TRUE
)

# 5.2 Editable publication-style three-line table workbook
wb_pub <- createWorkbook()
addWorksheet(wb_pub, "Three_line_table")
showGridLines(wb_pub, "Three_line_table", showGridLines = FALSE)

# title
writeData(wb_pub, "Three_line_table",
          x = "Table X. Differences in ecosystem multifunctionality indices between Futian (FT) and Gaoqiao (GQ) mangrove forests.",
          startRow = 1, startCol = 1)
mergeCells(wb_pub, "Three_line_table", cols = 1:ncol(pub_table), rows = 1)

# table body
start_row <- 3
writeData(wb_pub, "Three_line_table", pub_table, startRow = start_row, startCol = 1)

# styles
title_style <- createStyle(fontName = "Times New Roman", fontSize = 11, textDecoration = "bold")
header_style <- createStyle(
  fontName = "Times New Roman", fontSize = 10, textDecoration = "bold",
  halign = "center", valign = "center",
  border = c("top", "bottom"), borderStyle = "thin"
)
body_style <- createStyle(
  fontName = "Times New Roman", fontSize = 10,
  halign = "center", valign = "center"
)
last_row_style <- createStyle(
  fontName = "Times New Roman", fontSize = 10,
  halign = "center", valign = "center",
  border = "bottom", borderStyle = "thin"
)
note_style <- createStyle(fontName = "Times New Roman", fontSize = 9, wrapText = TRUE)

addStyle(wb_pub, "Three_line_table", title_style,
         rows = 1, cols = 1, gridExpand = TRUE, stack = TRUE)

addStyle(wb_pub, "Three_line_table", header_style,
         rows = start_row, cols = 1:ncol(pub_table), gridExpand = TRUE, stack = TRUE)

if (nrow(pub_table) > 1) {
  addStyle(wb_pub, "Three_line_table", body_style,
           rows = (start_row + 1):(start_row + nrow(pub_table) - 1),
           cols = 1:ncol(pub_table), gridExpand = TRUE, stack = TRUE)
}

addStyle(wb_pub, "Three_line_table", last_row_style,
         rows = start_row + nrow(pub_table),
         cols = 1:ncol(pub_table), gridExpand = TRUE, stack = TRUE)

setColWidths(wb_pub, "Three_line_table", cols = 1:ncol(pub_table), widths = "auto")

note_text <- paste(
  "Table note: Values are presented as mean ± SD. Differences between Futian (FT) and Gaoqiao (GQ)",
  "were tested using independent-sample t-tests when assumptions of normality and homogeneity of variance",
  "were satisfied; otherwise, Mann–Whitney U tests were applied. Effect sizes were calculated using",
  "Cohen's d for parametric tests and rank-biserial correlation for non-parametric tests.",
  "All tests are two-sided; Mann-Whitney P values use a normal approximation with",
  "continuity correction. Positive effect sizes indicate FT > GQ.",
  "P values are unadjusted for multiple comparisons.",
  "*P < 0.05, **P < 0.01, ***P < 0.001; ns, P >= 0.05."
)

note_start <- start_row + nrow(pub_table) + 3
writeData(wb_pub, "Three_line_table", note_text, startRow = note_start, startCol = 1)
mergeCells(wb_pub, "Three_line_table", cols = 1:ncol(pub_table), rows = note_start)
addStyle(wb_pub, "Three_line_table", note_style,
         rows = note_start, cols = 1, gridExpand = TRUE, stack = TRUE)
setRowHeights(wb_pub, "Three_line_table", rows = note_start, heights = 45)

saveWorkbook(
  wb_pub,
  file = file.path(output_dir, "EMF_difference_three_line_table.xlsx"),
  overwrite = TRUE
)

# 5.3 Native editable Word table: exactly three horizontal rules, no vertical rules.
# Abbreviations keep the same statistical content readable on a journal table page.
word_table <- pub_table
word_table$`Statistical test` <- ifelse(
  comp_table$Method[match(word_table$Function, comp_table$Index)] ==
    "Student's independent-samples t-test", "Student's t-test", "Mann-Whitney U"
)
word_table$`Effect size` <- paste0(
  ifelse(comp_table$Effect_type == "Cohen's d", "d = ", "r_rb = "),
  fmt_num(comp_table$Effect_size)
)
word_table <- as.data.frame(word_table, check.names = FALSE)
ft_word <- flextable::flextable(word_table)
ft_word <- flextable::set_header_labels(
  ft_word, Function = "Function", `FT (mean ± SD)` = "FT\nMean ± SD",
  `GQ (mean ± SD)` = "GQ\nMean ± SD", `Statistical test` = "Test",
  Statistic = "Statistic", `P value` = "P value", `Effect size` = "ES",
  Significance = "Sig."
)
ft_word <- flextable::font(ft_word, fontname = "Times New Roman", part = "all")
ft_word <- flextable::fontsize(ft_word, size = 10, part = "all")
ft_word <- flextable::color(ft_word, color = "black", part = "all")
ft_word <- flextable::bg(ft_word, bg = "white", part = "all")
ft_word <- flextable::bold(ft_word, part = "header")
ft_word <- flextable::align(ft_word, align = "center", part = "all")
ft_word <- flextable::align(ft_word, j = 1, align = "left", part = "all")
ft_word <- flextable::valign(ft_word, valign = "center", part = "all")
ft_word <- flextable::padding(ft_word, padding = 5, part = "all")
ft_word <- flextable::width(
  ft_word, j = seq_len(ncol(word_table)),
  width = c(1.35, 1.2, 1.2, 1.4, 1.0, 0.7, 1.0, 0.55)
)
ft_word <- flextable::border_remove(ft_word)
ft_word <- flextable::hline_top(
  ft_word, border = officer::fp_border(color = "black", width = 1.2), part = "header"
)
ft_word <- flextable::hline_bottom(
  ft_word, border = officer::fp_border(color = "black", width = 0.6), part = "header"
)
ft_word <- flextable::hline_bottom(
  ft_word, border = officer::fp_border(color = "black", width = 1.2), part = "body"
)
ft_word <- flextable::set_table_properties(
  ft_word, layout = "fixed", align = "center",
  opts_word = list(split = FALSE, repeat_headers = TRUE)
)
doc <- officer::read_docx()
doc <- officer::body_set_default_section(
  doc, officer::prop_section(
    page_size = officer::page_size(width = 8.5, height = 11, orient = "landscape"),
    page_margins = officer::page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75)
  )
)
doc <- officer::body_add_fpar(doc, officer::fpar(
  officer::ftext(
    "Table X. Differences in ecosystem multifunctionality between Futian and Gaoqiao mangrove forests",
    officer::fp_text(font.family = "Times New Roman", font.size = 11, bold = TRUE)
  ), fp_p = officer::fp_par(keep_with_next = TRUE, padding.bottom = 8)
))
doc <- flextable::body_add_flextable(doc, value = ft_word)
doc <- officer::body_add_fpar(doc, officer::fpar(
  officer::ftext(paste(note_text, "ES: effect size; d: Cohen's d; r_rb: rank-biserial correlation.",
                     "U is the Mann-Whitney statistic for FT; FT and GQ denote Futian and Gaoqiao."),
                officer::fp_text(font.family = "Times New Roman", font.size = 10)),
  fp_p = officer::fp_par(padding.top = 8)
))
print(doc, target = file.path(output_dir, "EMF_difference_three_line_table.docx"))

# ---------------------------
# 6. Plot functions
# ---------------------------
plot_labels <- c(
  "EMF-carbon" = "EMF-carbon",
  "EMF-nutrient" = "EMF-nutrient",
  "EMF-enzyme" = "EMF-enzyme",
  "EMF-vegetation" = "EMF-vegetation",
  "EMF-avg" = "EMF-avg"
)

# Okabe-Ito colors, kept identical in every panel and export.
region_palette <- c(FT = "#0072B2", GQ = "#D55E00")

make_raincloud <- function(ind, show_caption = FALSE) {
  df <- dat %>%
    select(Region, all_of(ind)) %>%
    rename(Value = all_of(ind)) %>%
    filter(!is.na(Value)) %>%
    mutate(Group_x = as.numeric(Region))
  st <- comp_table %>% filter(Index == ind)
  star_lab <- if (is.finite(st$P_value[1]) && st$P_value[1] < 0.05) {
    st$Significance[1]
  } else ""

  # Separate lanes: jittered points on the left, box at center, half-violin on right.
  # Explicit unbounded KDE and trim = FALSE preserve both density tails.
  p <- ggplot(df, aes(x = Group_x, y = Value, fill = Region, group = Region)) +
    ggdist::stat_halfeye(
      density = "unbounded", trim = FALSE, adjust = 0.9,
      width = 0.34, justification = 0, side = "right",
      position = position_nudge(x = 0.12), orientation = "vertical",
      .width = 0, show_point = FALSE, show_interval = FALSE,
      slab_colour = NA, alpha = 0.65
    ) +
    geom_boxplot(width = 0.14, alpha = 0.8, outlier.shape = NA,
                 colour = "black", linewidth = 0.45) +
    geom_point(
      aes(x = Group_x - 0.20),
      position = position_jitter(width = 0.065, height = 0, seed = 20261014),
      size = 1.9, shape = 21, stroke = 0.3, colour = "black", alpha = 0.85
    ) +
    scale_fill_manual(values = region_palette, drop = FALSE) +
    scale_x_continuous(breaks = seq_along(group_levels), labels = group_levels,
                       expand = expansion(mult = 0)) +
    labs(x = NULL, y = plot_labels[ind],
         caption = if (show_caption) paste(caption_common, panel_result(ind)) else NULL) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      plot.caption = element_text(size = 9, hjust = 0),
      plot.margin = margin(10, 12, 10, 10)
    )

  # Get actual untrimmed KDE endpoints before setting the view limits.
  # Using the raw data range here would cut off the violin tails.
  density_y <- ggplot_build(p)$data[[1]]$y
  y_bounds <- range(c(df$Value, density_y), finite = TRUE)
  y_range <- diff(y_bounds)
  if (!is.finite(y_range) || y_range <= 0) y_range <- 1
  bracket_y <- y_bounds[2] + y_range * 0.055
  upper_pad <- if (nzchar(star_lab)) 0.15 else 0.045
  p <- p + coord_cartesian(
    xlim = c(0.62, 2.55),
    ylim = c(y_bounds[1] - y_range * 0.045, y_bounds[2] + y_range * upper_pad),
    expand = FALSE, clip = "on"
  )
  if (nzchar(star_lab)) {
    # Display only the significance asterisks; hide non-significant labels.
    p <- p + annotate("text", x = 1.5, y = bracket_y,
                      label = star_lab, size = 5.2, colour = "black")
  }
  p
}

plots <- lapply(indices, make_raincloud, show_caption = FALSE)
names(plots) <- indices

# ---------------------------
# 7. Export TIFF figures
# ---------------------------
for (ind in indices) {
  ggsave(
    filename = file.path(output_dir, paste0(ind, ".tiff")),
    plot = plots[[ind]],
    width = 6,
    height = 5,
    dpi = 600,
    compression = "lzw",
    units = "in"
  )
}

# ---------------------------
# 8. Combined 4-panel figure and editable PPTX
# ---------------------------
function_fig <- (plots[["EMF-carbon"]] + plots[["EMF-nutrient"]]) /
  (plots[["EMF-enzyme"]] + plots[["EMF-vegetation"]]) +
  plot_annotation(
    tag_levels = "A",
    title = "Regional differences in ecosystem multifunctionality functional dimensions"
  )

ggsave(file.path(output_dir, "Figure_EMF_function_raincloud.tiff"), function_fig,
       width = 10, height = 8, dpi = 600, compression = "lzw", units = "in")
ggsave(file.path(output_dir, "Figure_EMF_function_raincloud_preview.png"), function_fig,
       width = 10, height = 8, dpi = 160, units = "in")

# 9. EMF-avg and editable PPTX figures with editable caption slides.
avg_plot <- plots[["EMF-avg"]] +
  labs(title = "Overall ecosystem multifunctionality")
ggsave(file.path(output_dir, "Figure_EMF_avg_raincloud_preview.png"), avg_plot,
       width = 6, height = 5, dpi = 160, units = "in")

export_editable_pptx <- function(fig, legend_text, filename) {
  ppt <- read_pptx()
  slide_dims <- slide_size(ppt)
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  ppt <- ph_with(ppt, value = dml(ggobj = fig),
                 location = ph_location(left = 0.3, top = 0.3,
                                        width = slide_dims$width - 0.6,
                                        height = slide_dims$height - 0.6))
  # A dedicated legend slide keeps the complete panel-specific caption readable.
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  ppt <- ph_with(ppt, value = fpar(ftext("Figure legend", fp_text(font.size = 18, bold = TRUE))),
                 location = ph_location(left = 0.45, top = 0.35,
                                        width = slide_dims$width - 0.9, height = 0.5))
  ppt <- ph_with(ppt,
                 value = fpar(ftext(legend_text, fp_text(font.size = 14, font.family = "Arial"))),
                 location = ph_location(left = 0.45, top = 1.05,
                                        width = slide_dims$width - 0.9,
                                        height = slide_dims$height - 1.4))
  print(ppt, target = file.path(output_dir, filename))
}
export_editable_pptx(function_fig, legend_functions, "Figure_EMF_function_raincloud_editable.pptx")
export_editable_pptx(avg_plot, legend_avg, "Figure_EMF_avg_raincloud_editable.pptx")

cat("All analyses and outputs completed successfully.\n")
