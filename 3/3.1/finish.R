# =============================================================================
# Mangrove carbon indicators: descriptive statistics, assumption checks,
# between-site comparisons, and publication-style boxplots
# Input : D:/Users/Lenovo/Desktop/原始/EMFclean.xlsx, worksheet "计算"
# Output: D:/Users/Lenovo/Desktop/原始/二
# Sites : Futian (福田 / SZ) vs Gaoqiao (高桥 / GQ)
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ---- 1. User settings --------------------------------------------------------
input_file <- "D:/Users/Lenovo/Desktop/原始/EMFclean.xlsx"
sheet_name <- "计算"
output_dir <- "D:/Users/Lenovo/Desktop/原始/二"

# Leave these as NULL to identify columns from their headers automatically.
# If automatic detection fails, set an override to the exact Excel header text,
# e.g. id_col_override <- "Plot ID".
id_col_override     <- NULL
region_col_override <- NULL
ocd_col_override    <- NULL
mbc_col_override    <- NULL

site_levels <- c("Futian", "Gaoqiao")
site_colors <- c("Futian" = "#287C8E", "Gaoqiao" = "#D49A38")

# ---- 2. Install/load packages -----------------------------------------------
required_packages <- c(
  "readxl", "dplyr", "ggplot2", "car", "openxlsx",
  "flextable", "officer", "readr", "ggbeeswarm",
  "patchwork", "rvg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

# ---- 3. Check paths and read Excel ------------------------------------------
if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPlease check the path in Section 1. Use forward slashes (/) in Windows paths."
  )
}
if (!sheet_name %in% readxl::excel_sheets(input_file)) {
  stop(
    "Worksheet '", sheet_name, "' was not found. Available worksheets: ",
    paste(readxl::excel_sheets(input_file), collapse = ", ")
  )
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- readxl::read_excel(input_file, sheet = sheet_name, .name_repair = "minimal")
header_names <- names(raw)

# Find the first column matching one of the header patterns.
find_column <- function(headers, patterns, label, override = NULL) {
  if (!is.null(override)) {
    if (!override %in% headers) {
      stop("The override for ", label, " was not found: '", override,
           "'. Available headers are: ", paste(headers, collapse = " | "))
    }
    return(override)
  }
  for (pat in patterns) {
    idx <- which(grepl(pat, headers, ignore.case = TRUE, perl = TRUE))
    if (length(idx) > 0) return(headers[idx[1]])
  }
  stop(
    "Cannot identify the column for ", label, ".\nAvailable headers are:\n- ",
    paste(headers, collapse = "\n- "),
    "\nPlease set the corresponding *_col_override in Section 1 to the exact header."
  )
}

id_col <- find_column(
  header_names,
  c("^\\s*plot\\s*id\\s*$", "plot.?id", "样点编号|样方编号|样点ID|样方ID|样点名称|样点"),
  "Plot ID", id_col_override
)
region_col <- find_column(
  header_names,
  c("^\\s*region\\s*$", "样地|地区|区域|研究区|site.?name|site.?group"),
  "Region", region_col_override
)
ocd_col <- find_column(
  header_names,
  c("有机碳密度", "organic\\s*carbon\\s*density", "organiccarbon.*density", "soc.?density"),
  "organic carbon density", ocd_col_override
)
mbc_col <- find_column(
  header_names,
  c("微生物.*生物量.*碳", "\\bMBC\\b", "microbial.*biomass.*carbon"),
  "microbial biomass carbon (MBC)", mbc_col_override
)

# Convert Excel numeric or numeric-looking text cells to numeric safely.
to_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "—", "-", "未测", "缺失")] <- NA_character_
  x <- gsub(",", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

classify_site <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("福田|futian|^SZ", x, ignore.case = TRUE, perl = TRUE) ~ "Futian",
    grepl("高桥|高橋|gao\\s*qiao|gao.?qiao|^GQ", x, ignore.case = TRUE, perl = TRUE) ~ "Gaoqiao",
    TRUE ~ NA_character_
  )
}

selected <- data.frame(
  Plot_ID = trimws(as.character(raw[[id_col]])),
  Region_raw = trimws(as.character(raw[[region_col]])),
  Organic_carbon_density_kg_m2 = to_numeric(raw[[ocd_col]]),
  MBC_mg_kg = to_numeric(raw[[mbc_col]]),
  stringsAsFactors = FALSE
)
selected$Site <- classify_site(selected$Region_raw)
selected$Plot_ID[is.na(selected$Plot_ID) | selected$Plot_ID == ""] <- NA_character_
selected$Region_raw[is.na(selected$Region_raw) | selected$Region_raw == ""] <- "[missing Region]"

# Record the Region-to-site mapping explicitly so unmapped labels cannot be hidden.
site_mapping <- selected |>
  dplyr::count(Region_raw, Site, name = "Rows") |>
  dplyr::arrange(is.na(Site), Region_raw)

analysis_wide <- selected |>
  dplyr::filter(!is.na(Site)) |>
  dplyr::mutate(Site = factor(Site, levels = site_levels))

if (nrow(analysis_wide) == 0) {
  stop(
    "No site labels were recognized. The Region values need to include 福田/Futian/SZ or 高桥/Gaoqiao/GQ.\n",
    "Check the '分组映射' table logic in the script and adapt classify_site() to your actual labels."
  )
}
if (length(unique(stats::na.omit(as.character(analysis_wide$Site)))) < 2) {
  stop("Only one study site was recognized. Check the Region column and classify_site() mapping.")
}

# Long data: one row per site x indicator x plot.
dat_long <- dplyr::bind_rows(
  data.frame(
    Plot_ID = analysis_wide$Plot_ID,
    Site = analysis_wide$Site,
    Indicator = "Organic carbon density",
    Value = analysis_wide$Organic_carbon_density_kg_m2
  ),
  data.frame(
    Plot_ID = analysis_wide$Plot_ID,
    Site = analysis_wide$Site,
    Indicator = "Microbial biomass carbon",
    Value = analysis_wide$MBC_mg_kg
  )
) |>
  dplyr::filter(!is.na(Value), !is.na(Site)) |>
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_levels),
    Indicator = factor(Indicator, levels = c("Organic carbon density", "Microbial biomass carbon"))
  )

if (nrow(dat_long) == 0) stop("The selected indicator columns contain no usable numeric values.")

# ---- 4. Descriptive statistics ----------------------------------------------
summarise_values <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  avg <- if (n > 0) mean(x) else NA_real_
  sdev <- if (n > 1) stats::sd(x) else NA_real_
  minv <- if (n > 0) min(x) else NA_real_
  maxv <- if (n > 0) max(x) else NA_real_
  data.frame(
    N = n,
    Mean = avg,
    SD = sdev,
    CV_percent = if (is.finite(avg) && avg != 0 && is.finite(sdev)) 100 * sdev / abs(avg) else NA_real_,
    Min = minv,
    Max = maxv,
    Range = if (n > 0) paste0(format(minv, digits = 5, trim = TRUE), "–", format(maxv, digits = 5, trim = TRUE)) else NA_character_,
    Median = if (n > 0) stats::median(x) else NA_real_,
    Q1 = if (n > 0) as.numeric(stats::quantile(x, 0.25, names = FALSE)) else NA_real_,
    Q3 = if (n > 0) as.numeric(stats::quantile(x, 0.75, names = FALSE)) else NA_real_,
    IQR = if (n > 0) stats::IQR(x) else NA_real_
  )
}

indicator_names <- levels(dat_long$Indicator)
site_names <- site_levels

desc_by_site <- dat_long |>
  dplyr::group_by(Indicator, Site) |>
  dplyr::group_modify(~ summarise_values(.x$Value)) |>
  dplyr::ungroup() |>
  dplyr::mutate(Site = as.character(Site))

desc_all <- dat_long |>
  dplyr::group_by(Indicator) |>
  dplyr::group_modify(~ summarise_values(.x$Value)) |>
  dplyr::ungroup() |>
  dplyr::mutate(Site = "All sites") |>
  dplyr::select(Indicator, Site, dplyr::everything())

descriptive_stats <- dplyr::bind_rows(desc_by_site, desc_all) |>
  dplyr::select(Indicator, Site, N, Mean, SD, CV_percent, Min, Max, Range, Median, Q1, Q3, IQR) |>
  dplyr::arrange(match(Indicator, indicator_names), match(Site, c(site_levels, "All sites")))

# ---- 5. Assumption tests and group comparison -------------------------------
shapiro_rows <- list()
levene_rows <- list()
test_rows <- list()

for (ind in indicator_names) {
  d <- dat_long |>
    dplyr::filter(as.character(Indicator) == ind) |>
    dplyr::mutate(Site = factor(as.character(Site), levels = site_levels))

  # Shapiro-Wilk normality test separately in each site (valid for 3-5000 values).
  for (site in site_levels) {
    x <- d$Value[as.character(d$Site) == site]
    x <- x[is.finite(x)]
    if (length(x) >= 3 && length(x) <= 5000 && length(unique(x)) >= 2) {
      st <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
      w_val <- if (!is.null(st)) unname(st$statistic) else NA_real_
      p_val <- if (!is.null(st)) st$p.value else NA_real_
      status <- if (is.na(p_val)) "Test unavailable" else if (p_val >= 0.05) "No significant departure detected (p >= 0.05)" else "Departure from normality detected (p < 0.05)"
    } else {
      w_val <- NA_real_
      p_val <- NA_real_
      status <- "Not tested: requires 3-5000 values and at least 2 unique values"
    }
    shapiro_rows[[length(shapiro_rows) + 1]] <- data.frame(
      Indicator = ind, Site = site, N = length(x), W = w_val, P_value = p_val,
      Interpretation = status, stringsAsFactors = FALSE
    )
  }

  # Brown-Forsythe / median-centred Levene test for homogeneity of variance.
  lev <- tryCatch(car::leveneTest(Value ~ Site, data = d, center = median), error = function(e) NULL)
  if (!is.null(lev)) {
    lev_p <- as.numeric(lev[1, "Pr(>F)"])
    lev_f <- as.numeric(lev[1, "F value"])
    lev_df1 <- as.numeric(lev[1, "Df"])
    lev_df2 <- if (nrow(lev) >= 2) as.numeric(lev[2, "Df"]) else NA_real_
    lev_interpretation <- if (is.na(lev_p)) "Test unavailable" else if (lev_p >= 0.05) "No significant variance heterogeneity detected (p >= 0.05)" else "Variance heterogeneity detected (p < 0.05)"
  } else {
    lev_p <- lev_f <- lev_df1 <- lev_df2 <- NA_real_
    lev_interpretation <- "Test unavailable"
  }
  levene_rows[[length(levene_rows) + 1]] <- data.frame(
    Indicator = ind, F = lev_f, df1 = lev_df1, df2 = lev_df2,
    P_value = lev_p, Interpretation = lev_interpretation, stringsAsFactors = FALSE
  )

  this_shapiro <- dplyr::bind_rows(shapiro_rows) |>
    dplyr::filter(Indicator == ind, Site %in% site_levels)
  normality_ok <- nrow(this_shapiro) == 2 && all(!is.na(this_shapiro$P_value)) && all(this_shapiro$P_value >= 0.05)
  variance_ok <- !is.na(lev_p) && lev_p >= 0.05

  # Pre-specified decision rule for this analysis:
  # normality in both groups AND equal variances -> pooled independent-samples t-test;
  # otherwise -> Kruskal-Wallis rank-sum test, as requested.
  if (normality_ok && variance_ok) {
    comparison <- tryCatch(stats::t.test(Value ~ Site, data = d, var.equal = TRUE), error = function(e) NULL)
    if (!is.null(comparison)) {
      method <- "Student's independent-samples t-test"
      statistic <- unname(comparison$statistic)
      df_value <- unname(comparison$parameter)
      p_raw <- comparison$p.value
      futian_mean <- mean(d$Value[as.character(d$Site) == "Futian"], na.rm = TRUE)
      gaoqiao_mean <- mean(d$Value[as.character(d$Site) == "Gaoqiao"], na.rm = TRUE)
      estimate <- futian_mean - gaoqiao_mean
      ci_low <- unname(comparison$conf.int[1])
      ci_high <- unname(comparison$conf.int[2])
      estimate_type <- "Mean difference (Futian - Gaoqiao)"
    } else {
      comparison <- stats::kruskal.test(Value ~ Site, data = d)
      method <- "Kruskal-Wallis test (fallback after t-test error)"
      statistic <- unname(comparison$statistic)
      df_value <- unname(comparison$parameter)
      p_raw <- comparison$p.value
      estimate <- stats::median(d$Value[as.character(d$Site) == "Futian"], na.rm = TRUE) -
        stats::median(d$Value[as.character(d$Site) == "Gaoqiao"], na.rm = TRUE)
      ci_low <- ci_high <- NA_real_
      estimate_type <- "Median difference (Futian - Gaoqiao; descriptive only)"
    }
  } else {
    comparison <- stats::kruskal.test(Value ~ Site, data = d)
    method <- "Kruskal-Wallis rank-sum test"
    statistic <- unname(comparison$statistic)
    df_value <- unname(comparison$parameter)
    p_raw <- comparison$p.value
    estimate <- stats::median(d$Value[as.character(d$Site) == "Futian"], na.rm = TRUE) -
      stats::median(d$Value[as.character(d$Site) == "Gaoqiao"], na.rm = TRUE)
    ci_low <- ci_high <- NA_real_
    estimate_type <- "Median difference (Futian - Gaoqiao; descriptive only)"
  }

  test_rows[[length(test_rows) + 1]] <- data.frame(
    Indicator = ind,
    N_Futian = sum(as.character(d$Site) == "Futian"),
    N_Gaoqiao = sum(as.character(d$Site) == "Gaoqiao"),
    Shapiro_normality_passed_both_sites = normality_ok,
    Levene_p = lev_p,
    Variance_homogeneity_passed = variance_ok,
    Test_method = method,
    Statistic = statistic,
    df = df_value,
    Estimate = estimate,
    Estimate_definition = estimate_type,
    CI95_lower = ci_low,
    CI95_upper = ci_high,
    P_raw = p_raw,
    stringsAsFactors = FALSE
  )
}

shapiro_results <- dplyr::bind_rows(shapiro_rows)
levene_results <- dplyr::bind_rows(levene_rows)
p_to_stars <- function(p) {
  if (length(p) == 0 || is.na(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

test_results <- dplyr::bind_rows(test_rows) |>
  dplyr::mutate(
    P_BH_adjusted = stats::p.adjust(P_raw, method = "BH")
  )

# Effect size for manuscript reporting.
# Parametric comparison: Hedges' g; non-parametric comparison: Cliff's delta.
calc_effect_size <- function(indicator_name) {
  d <- dat_long |> dplyr::filter(Indicator == indicator_name)
  x <- d$Value[as.character(d$Site) == "Futian"]
  y <- d$Value[as.character(d$Site) == "Gaoqiao"]
  method <- test_results$Test_method[test_results$Indicator == indicator_name][1]

  if (grepl("Student", method, fixed = TRUE)) {
    nx <- length(x); ny <- length(y)
    sp <- sqrt(((nx - 1) * stats::var(x) + (ny - 1) * stats::var(y)) / (nx + ny - 2))
    g <- ((mean(x) - mean(y)) / sp) * (1 - 3 / (4 * (nx + ny) - 9))
    data.frame(Effect_size = g, Effect_size_method = "Hedges' g")
  } else {
    delta <- (sum(outer(x, y, "-") > 0) - sum(outer(x, y, "-") < 0)) / (length(x) * length(y))
    data.frame(Effect_size = delta, Effect_size_method = "Cliff's delta")
  }
}

effect_sizes <- dplyr::bind_rows(lapply(levels(dat_long$Indicator), calc_effect_size)) |>
  dplyr::mutate(Indicator = levels(dat_long$Indicator))

test_results <- dplyr::left_join(test_results, effect_sizes, by = "Indicator") |>
  dplyr::mutate(Significance = vapply(P_BH_adjusted, p_to_stars, character(1)))

# ---- 6. Final submission-style figures ---------------------------------------
# This section fixes the figure size, font size, line width, colour palette,
# point style, and caption layout for a consistent submission-ready style.
format_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

format_p_sentence <- function(p, prefix = "P") {
  if (length(p) == 0 || is.na(p)) return(paste0(prefix, " = NA"))
  if (p < 0.001) return(paste0(prefix, " = ", formatC(p, format = "e", digits = 2)))
  paste0(prefix, " = ", formatC(p, format = "f", digits = 3))
}

get_y_label <- function(ind) {
  if (ind == "Organic carbon density") return("Organic carbon density (kg m\u207B\u00B2)")
  "Microbial biomass carbon (mg kg\u207B\u00B9)"
}

get_short_test_label <- function(test_method) {
  if (grepl("Student", test_method, fixed = TRUE)) {
    "Student's t-test"
  } else {
    "Kruskal\u2013Wallis test"
  }
}

get_full_test_label <- function(test_method) {
  if (grepl("Student", test_method, fixed = TRUE)) {
    "pooled independent-samples Student's t-test"
  } else {
    "Kruskal\u2013Wallis rank-sum test"
  }
}

# Final fixed graphical parameters. Keep all figure geometry here for easy tuning.
fig_theme <- list(
  font_family = "Arial",
  axis_text_pt = 8.5,
  axis_title_pt = 9.5,
  caption_pt = 8.0,
  panel_tag_pt = 10.0,
  box_width = 0.52,
  box_alpha = 0.60,       # 40% transparent fill
  box_linewidth = 0.70,
  axis_linewidth = 0.42,
  tick_linewidth = 0.36,
  bracket_linewidth = 0.46,
  point_size = 1.75,
  point_stroke = 0.30,
  beeswarm_width = 0.18,
  single_width_mm = 88,
  single_height_mm = 96,
  combined_width_mm = 178,
  combined_height_mm = 138,
  figure_dpi = 600
)

# Muted Okabe-Ito-inspired, colour-blind-friendly blue and orange.
site_colors <- c("Futian" = "#79AFC7", "Gaoqiao" = "#D6A24E")

make_plot <- function(ind, panel_tag = NULL) {
  d <- dat_long |>
    dplyr::filter(as.character(Indicator) == ind) |>
    dplyr::mutate(Site = factor(as.character(Site), levels = site_levels))
  tr <- test_results |>
    dplyr::filter(Indicator == ind) |>
    dplyr::slice(1)

  vals <- d$Value[is.finite(d$Value)]
  vmin <- min(vals, na.rm = TRUE)
  vmax <- max(vals, na.rm = TRUE)
  spread <- vmax - vmin
  if (!is.finite(spread) || spread <= 0) spread <- max(abs(vmax), 1) * 0.20

  # A fixed rule for balanced y limits and a higher significance bracket.
  y_lower <- vmin - 0.08 * spread
  bracket_y <- vmax + 0.17 * spread
  star_y <- bracket_y + 0.030 * spread
  y_upper <- star_y + 0.085 * spread

  if (!is.finite(y_lower)) y_lower <- vmin
  if (!is.finite(y_upper)) y_upper <- vmax + 0.30 * spread

  star_label <- tr$Significance[1]

  box_args <- list(
    width = fig_theme$box_width,
    alpha = fig_theme$box_alpha,
    outlier.shape = NA,
    colour = "#1F2933",
    show.legend = FALSE
  )
  if (utils::packageVersion("ggplot2") >= "3.4.0") {
    box_args$linewidth <- fig_theme$box_linewidth
  } else {
    box_args$size <- fig_theme$box_linewidth
  }
  if ("staplewidth" %in% names(formals(ggplot2::geom_boxplot))) {
    box_args$staplewidth <- 0.48
  }
  box_layer <- do.call(ggplot2::geom_boxplot, box_args)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = Site, y = Value, fill = Site)) +
    box_layer +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(fill = Site),
      shape = 21,
      size = fig_theme$point_size,
      stroke = fig_theme$point_stroke,
      colour = "#111111",
      alpha = 0.96,
      width = fig_theme$beeswarm_width,
      varwidth = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = site_colors, drop = FALSE) +
    ggplot2::annotate(
      "segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
      linewidth = fig_theme$bracket_linewidth, colour = "#1F2933"
    ) +
    ggplot2::annotate(
      "segment", x = 1, xend = 1, y = bracket_y,
      yend = bracket_y - 0.027 * spread,
      linewidth = fig_theme$bracket_linewidth, colour = "#1F2933"
    ) +
    ggplot2::annotate(
      "segment", x = 2, xend = 2, y = bracket_y,
      yend = bracket_y - 0.027 * spread,
      linewidth = fig_theme$bracket_linewidth, colour = "#1F2933"
    ) +
    ggplot2::annotate(
      "text", x = 1.5, y = star_y, label = star_label,
      size = 3.65, fontface = "bold", colour = "#111111", vjust = 0
    ) +
    ggplot2::scale_y_continuous(
      limits = c(y_lower, y_upper),
      breaks = pretty(c(y_lower, y_upper), n = 5),
      labels = function(x) format(x, trim = TRUE, scientific = FALSE)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = get_y_label(ind)) +
    ggplot2::theme_classic(base_size = fig_theme$axis_title_pt, base_family = fig_theme$font_family) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        family = fig_theme$font_family, colour = "#263238",
        size = fig_theme$axis_text_pt, margin = ggplot2::margin(t = 4)
      ),
      axis.text.y = ggplot2::element_text(
        family = fig_theme$font_family, colour = "#263238",
        size = fig_theme$axis_text_pt
      ),
      axis.title.x = ggplot2::element_text(
        family = fig_theme$font_family, size = fig_theme$axis_title_pt,
        margin = ggplot2::margin(t = 7)
      ),
      axis.title.y = ggplot2::element_text(
        family = fig_theme$font_family, face = "plain",
        size = fig_theme$axis_title_pt, margin = ggplot2::margin(r = 7)
      ),
      axis.line = ggplot2::element_line(colour = "#30363B", linewidth = fig_theme$axis_linewidth),
      axis.ticks = ggplot2::element_line(colour = "#30363B", linewidth = fig_theme$tick_linewidth),
      axis.ticks.length = grid::unit(2.0, "mm"),
      plot.margin = ggplot2::margin(t = 10, r = 8, b = 6, l = 7),
      legend.position = "none"
    )

  if (!is.null(panel_tag)) {
    p <- p +
      ggplot2::annotate(
        "text", x = -Inf, y = Inf, label = panel_tag,
        hjust = -0.25, vjust = 1.25,
        family = fig_theme$font_family,
        size = 3.7, fontface = "bold", colour = "#111111"
      )
  }
  p
}

p_ocd <- make_plot("Organic carbon density", panel_tag = NULL)
p_mbc <- make_plot("Microbial biomass carbon", panel_tag = NULL)
p_ocd_tagged <- make_plot("Organic carbon density", panel_tag = "A")
p_mbc_tagged <- make_plot("Microbial biomass carbon", panel_tag = "B")

# Sample sizes for the caption are derived from the non-missing data used in each panel.
count_caption_n <- function(ind, site) {
  dat_long |>
    dplyr::filter(as.character(Indicator) == ind, as.character(Site) == site) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::pull(n)
}

n_ocd_f <- count_caption_n("Organic carbon density", "Futian")
n_ocd_g <- count_caption_n("Organic carbon density", "Gaoqiao")
n_mbc_f <- count_caption_n("Microbial biomass carbon", "Futian")
n_mbc_g <- count_caption_n("Microbial biomass carbon", "Gaoqiao")

tr_ocd <- test_results |> dplyr::filter(Indicator == "Organic carbon density") |> dplyr::slice(1)
tr_mbc <- test_results |> dplyr::filter(Indicator == "Microbial biomass carbon") |> dplyr::slice(1)

# Universal number formatter used for manuscript tables and figure captions.
format_number <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

caption_text <- paste0(
  "Figure 1. Organic carbon density (A) and microbial biomass carbon (B) across mangrove sites (Futian and Gaoqiao). ",
  "Boxes indicate median and interquartile range; whiskers represent 1.5 × IQR; points show individual plots. ",
  "Differences were tested using the Shapiro–Wilk and median-centred Levene tests followed by Student's t-test or Kruskal–Wallis test as appropriate. ",
  "P values shown are Benjamini–Hochberg adjusted. ",
  "Organic carbon density: P = ", format_p(tr_ocd$P_BH_adjusted[1]),
  ", effect size = ", format_number(tr_ocd$Effect_size[1], 3), " (", tr_ocd$Effect_size_method[1], "); ",
  "microbial biomass carbon: P = ", format_p(tr_mbc$P_BH_adjusted[1]),
  ", effect size = ", format_number(tr_mbc$Effect_size[1], 3), " (", tr_mbc$Effect_size_method[1], ")."
)

wrapped_caption <- paste(strwrap(caption_text, width = 150), collapse = "\n")

combined_panels <- make_plot("Organic carbon density", panel_tag = "A") +
  make_plot("Microbial biomass carbon", panel_tag = "B") +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1))

caption_plot <- ggplot2::ggplot() +
  ggplot2::theme_void() +
  ggplot2::annotate(
    "text", x = 0, y = 1, label = wrapped_caption,
    hjust = 0, vjust = 1,
    family = fig_theme$font_family,
    size = 2.75, lineheight = 1.08, colour = "#20252A"
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 8, b = 0, l = 8))

final_figure <- patchwork::plot_spacer() / combined_panels / caption_plot +
  patchwork::plot_layout(heights = c(1, 0.28))

# Output files: no combined PNG and no SVG are generated in this final version.
fig_ocd_tiff <- file.path(output_dir, "Organic_Carbon_Density_Boxplot.tiff")
fig_mbc_tiff <- file.path(output_dir, "MBC_Boxplot.tiff")
fig_combined_pptx <- file.path(output_dir, "Figure1_Carbon_Indicators_with_Caption_editable.pptx")
fig_caption_txt <- file.path(output_dir, "Figure1_Caption.txt")

ggplot2::ggsave(
  filename = fig_ocd_tiff,
  plot = p_ocd,
  width = fig_theme$single_width_mm,
  height = fig_theme$single_height_mm,
  units = "mm",
  dpi = fig_theme$figure_dpi,
  compression = "lzw",
  limitsize = FALSE,
  bg = "white"
)

ggplot2::ggsave(
  filename = fig_mbc_tiff,
  plot = p_mbc,
  width = fig_theme$single_width_mm,
  height = fig_theme$single_height_mm,
  units = "mm",
  dpi = fig_theme$figure_dpi,
  compression = "lzw",
  limitsize = FALSE,
  bg = "white"
)

readr::write_lines(caption_text, fig_caption_txt)

# Editable composite figure: vector panels plus editable caption text in PowerPoint.
tryCatch({
  ppt <- officer::read_pptx()
  ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
  ppt <- officer::ph_with(
    ppt,
    value = rvg::dml(ggobj = combined_panels),
    location = officer::ph_location(left = 0.48, top = 0.15, width = 12.38, height = 4.80)
  )
  ppt <- officer::ph_with(
    ppt,
    value = wrapped_caption,
    location = officer::ph_location(left = 0.60, top = 5.10, width = 12.10, height = 1.55)
  )
  print(ppt, target = fig_combined_pptx)
}, error = function(e) {
  warning("Editable PPTX figure was not written: ", conditionMessage(e))
})

# ---- 7. Data quality checks --------------------------------------------------
site_indicator_n <- dat_long |>
  dplyr::group_by(Indicator, Site) |>
  dplyr::summarise(N_non_missing = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(Site = as.character(Site))

qc_info <- data.frame(
  Item = c(
    "Expected plot count (as specified by user)",
    "Rows in source worksheet",
    "Rows with non-missing Plot ID",
    "Unique non-missing Plot IDs",
    "Duplicated Plot ID rows (duplicate instances after first occurrence)",
    "Rows with recognized site label",
    "Recognized Futian rows",
    "Recognized Gaoqiao rows",
    "Rows excluded from site comparison because Region was unmapped",
    "Organic carbon density non-missing values (recognized sites)",
    "MBC non-missing values (recognized sites)"
  ),
  Value = c(
    143,
    nrow(raw),
    sum(!is.na(selected$Plot_ID)),
    dplyr::n_distinct(selected$Plot_ID[!is.na(selected$Plot_ID)]),
    sum(duplicated(selected$Plot_ID[!is.na(selected$Plot_ID)])),
    sum(!is.na(selected$Site)),
    sum(selected$Site == "Futian", na.rm = TRUE),
    sum(selected$Site == "Gaoqiao", na.rm = TRUE),
    sum(is.na(selected$Site)),
    sum(!is.na(analysis_wide$Organic_carbon_density_kg_m2), na.rm = TRUE),
    sum(!is.na(analysis_wide$MBC_mg_kg), na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

method_notes <- data.frame(
  Topic = c(
    "Input and grouping", "Mean", "Standard deviation (SD)", "Coefficient of variation (CV)",
    "Range", "Shapiro-Wilk", "Levene test", "Comparison decision rule", "Multiple testing",
    "Interpretation note", "GitHub / reproducibility", "Figures and tables"
  ),
  Details = c(
    paste0("Worksheet: ", sheet_name, "; sites mapped from Region as Futian/福田/SZ and Gaoqiao/高桥/GQ."),
    "Arithmetic mean calculated from non-missing plot-level values.",
    "Sample standard deviation (denominator n - 1).",
    "CV (%) = SD / abs(mean) × 100; reported as NA if the mean is zero or cannot be calculated.",
    "Minimum and maximum values are reported separately; Range is a display field combining them.",
    "Run separately within each site. P >= 0.05 means no significant departure from normality was detected; it does not prove normality.",
    "Median-centred Levene (Brown-Forsythe) test; P >= 0.05 means no significant variance heterogeneity was detected.",
    "Use pooled independent-samples Student's t-test only if both site-specific Shapiro tests have P >= 0.05 and Levene P >= 0.05; otherwise use Kruskal-Wallis, following the requested rule.",
    "Benjamini-Hochberg adjusted P values are provided across the two indicator-level comparisons; the figure displays significance codes based on adjusted P.",
    "Kruskal-Wallis tests rank distributions. Interpret as a difference in location/median only when group distribution shapes are reasonably similar. Report effect sizes and distributions, not P values alone.",
    "The analysis uses standard CRAN packages and does not require GitHub. For reproducibility, commit this script and package/session information to a private or approved GitHub repository; do not upload sensitive/raw data publicly.",
    "Final submission-style figures are exported as two separate 600-dpi TIFF panels, one captioned combined 600-dpi TIFF figure, and one editable PowerPoint composite figure; combined PNG and SVG files are intentionally not generated. Editable results are also saved to XLSX, UTF-8-BOM CSV, and a Word document with three-line tables."
  ),
  stringsAsFactors = FALSE
)

# ---- 8. Export results -------------------------------------------------------
output_xlsx <- file.path(output_dir, "Carbon_Indicators_Statistical_Results.xlsx")
result_sheets <- list(
  "Descriptive statistics" = descriptive_stats,
  "Shapiro-Wilk" = shapiro_results,
  "Levene" = levene_results,
  "Site comparisons" = test_results,
  "Region mapping" = site_mapping,
  "Data QC" = qc_info,
  "N by site and indicator" = site_indicator_n,
  "Analysis data" = analysis_wide |>
    dplyr::select(Plot_ID, Region_raw, Site, Organic_carbon_density_kg_m2, MBC_mg_kg),
  "Method notes" = method_notes
)

wb <- openxlsx::createWorkbook()
header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF", fgFill = "#173F4F", textDecoration = "bold",
  halign = "left", valign = "center", border = "Bottom", borderColour = "#173F4F"
)
body_style <- openxlsx::createStyle(valign = "top", wrapText = TRUE)
for (sheet_nm in names(result_sheets)) {
  tab <- result_sheets[[sheet_nm]]
  openxlsx::addWorksheet(wb, sheet_nm, gridLines = FALSE)
  openxlsx::writeData(wb, sheet_nm, tab, withFilter = nrow(tab) > 0)
  if (ncol(tab) > 0) {
    openxlsx::addStyle(wb, sheet_nm, header_style, rows = 1, cols = seq_len(ncol(tab)),
                       gridExpand = TRUE, stack = TRUE)
    if (nrow(tab) > 0) {
      openxlsx::addStyle(wb, sheet_nm, body_style, rows = 2:(nrow(tab) + 1),
                         cols = seq_len(ncol(tab)), gridExpand = TRUE, stack = TRUE)
    }
    widths <- rep(18, ncol(tab))
    names_tab <- names(tab)
    if ("Indicator" %in% names_tab) widths[match("Indicator", names_tab)] <- 30
    if ("Region_raw" %in% names_tab) widths[match("Region_raw", names_tab)] <- 24
    if ("Item" %in% names_tab) widths[match("Item", names_tab)] <- 48
    if ("Details" %in% names_tab) widths[match("Details", names_tab)] <- 88
    if ("Interpretation" %in% names_tab) widths[match("Interpretation", names_tab)] <- 58
    if ("Estimate_definition" %in% names_tab) widths[match("Estimate_definition", names_tab)] <- 48
    openxlsx::setColWidths(wb, sheet_nm, cols = seq_len(ncol(tab)), widths = widths)
    openxlsx::setRowHeights(wb, sheet_nm, rows = 1, heights = 24)
    openxlsx::freezePane(wb, sheet_nm, firstRow = TRUE)
  }
}
openxlsx::saveWorkbook(wb, output_xlsx, overwrite = TRUE)


# ---- 8A. Editable manuscript tables (Word, three-line/booktabs style) --------
format_p_table <- function(x) {
  ifelse(is.na(x), "—", ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3)))
}

indicator_label <- function(x) {
  dplyr::recode(
    x,
    "Organic carbon density" = "Organic carbon density (kg m\u207B\u00B2)",
    "Microbial biomass carbon" = "Microbial biomass carbon (mg kg\u207B\u00B9)"
  )
}

# Avoid Unicode escape sequences in backtick-quoted R names: those can trigger
# "\uxxxx sequences not supported inside backticks" on some R parsers.
# Build the columns with ASCII-safe names first, then assign publication labels.
paper_desc <- descriptive_stats |>
  dplyr::filter(Site %in% site_levels) |>
  dplyr::mutate(
    Indicator = indicator_label(Indicator),
    Mean_SD = paste0(format_number(Mean), " \u00B1 ", format_number(SD)),
    CV_percent_display = format_number(CV_percent, 1),
    Range_min_max = paste0(format_number(Min), "\u2013", format_number(Max))
  ) |>
  dplyr::transmute(
    Indicator, Site, n = N, Mean_SD, CV_percent_display, Range_min_max
  )
names(paper_desc)[names(paper_desc) == "Mean_SD"] <- "Mean ± SD"
names(paper_desc)[names(paper_desc) == "CV_percent_display"] <- "CV (%)"
names(paper_desc)[names(paper_desc) == "Range_min_max"] <- "Range (min–max)"

paper_tests <- test_results |>
  dplyr::mutate(
    Indicator = indicator_label(Indicator),
    Test = ifelse(
      Test_method == "Student's independent-samples t-test",
      "Student's t-test", "Kruskal\u2013Wallis"
    ),
    Statistic = ifelse(
      Test_method == "Student's independent-samples t-test",
      paste0("t = ", format_number(Statistic), "; df = ", format_number(df, 0)),
      paste0("H = ", format_number(Statistic), "; df = ", format_number(df, 0))
    ),
    P = format_p_table(P_raw),
    `BH-adjusted P` = format_p_table(P_BH_adjusted),
    Effect_size = format_number(Effect_size, 3),
    Effect_method = Effect_size_method,
    Significance = Significance
  ) |>
  dplyr::select(Indicator, Test, Statistic, P, `BH-adjusted P`,
                Effect_size, Effect_method, Significance)

make_booktabs_table <- function(x, font_size = 9) {
  ft <- flextable::flextable(x)
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = font_size, part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::align(ft, align = "center", part = "all")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  ft <- flextable::padding(ft, padding = 4, part = "all")
  ft <- flextable::autofit(ft)
  ft <- flextable::set_table_properties(ft, layout = "autofit", width = 1)
  ft
}

ft_desc <- make_booktabs_table(paper_desc, font_size = 8.5)
ft_tests <- make_booktabs_table(paper_tests, font_size = 8.5)

add_table_caption <- function(doc, text) {
  officer::body_add_fpar(
    doc,
    value = officer::fpar(officer::ftext(
      text,
      prop = officer::fp_text(font.family = "Arial", font.size = 10, bold = TRUE)
    ))
  )
}
add_table_note <- function(doc, text) {
  officer::body_add_fpar(
    doc,
    value = officer::fpar(officer::ftext(
      text,
      prop = officer::fp_text(font.family = "Arial", font.size = 8.5)
    ))
  )
}

table_docx <- officer::read_docx()
table_docx <- add_table_caption(
  table_docx,
  "Table 1. Descriptive statistics of carbon indicators by study site"
)
table_docx <- flextable::body_add_flextable(table_docx, value = ft_desc)
table_docx <- add_table_note(
  table_docx,
  "Note. Values are mean \u00B1 SD. CV = coefficient of variation; range is minimum\u2013maximum. Units are shown in the indicator names."
)
table_docx <- officer::body_add_par(table_docx, "", style = "Normal")
table_docx <- add_table_caption(
  table_docx,
  "Table 2. Statistical comparisons between Futian and Gaoqiao"
)
table_docx <- flextable::body_add_flextable(table_docx, value = ft_tests)
table_docx <- add_table_note(
  table_docx,
  "Note. The comparison method was selected using the prespecified Shapiro\u2013Wilk and median-centred Levene tests. BH-adjusted P values control the false-discovery rate across the two indicator comparisons. ns, not significant; * P < 0.05; ** P < 0.01; *** P < 0.001 (based on BH-adjusted P)."
)

table_docx_path <- file.path(output_dir, "Carbon_Indicators_Three_Line_Tables.docx")
print(table_docx, target = table_docx_path)

# UTF-8-with-BOM CSV copies: Excel on Windows can open these without mojibake.
readr::write_excel_csv(descriptive_stats, file.path(output_dir, "Descriptive_Statistics.csv"), na = "")
readr::write_excel_csv(shapiro_results, file.path(output_dir, "Shapiro_Wilk_Results.csv"), na = "")
readr::write_excel_csv(levene_results, file.path(output_dir, "Levene_Results.csv"), na = "")
readr::write_excel_csv(test_results, file.path(output_dir, "Site_Comparison_Results.csv"), na = "")

# Save the execution environment for reproducibility.
sink(file.path(output_dir, "R_sessionInfo.txt"))
cat("Analysis completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", sep = "")
cat("Detected columns:\n")
cat("Plot ID: ", id_col, "\nRegion: ", region_col,
    "\nOrganic carbon density: ", ocd_col,
    "\nMBC: ", mbc_col, "\n\n", sep = "")
cat("R session information:\n")
print(utils::sessionInfo())
sink()

message("Analysis complete. Results saved to: ", output_dir)
message("Individual TIFF: ", fig_ocd_tiff)
message("Individual TIFF: ", fig_mbc_tiff)
message("Editable combined PPTX: ", fig_combined_pptx)
message("Figure caption text file: ", fig_caption_txt)
message("Editable statistics workbook: ", output_xlsx)
message("Editable three-line Word tables: ", table_docx_path)

# =============================================================================
# Interpretation checklist before writing the manuscript:
# 1) Check 'Region mapping' and 'Data QC' sheets to confirm all 143 plots are
#    assigned to the intended sites and no unexpected rows were excluded.
# 2) Confirm indicator units and outliers against the field/lab records.
# 3) Report mean +/- SD, CV, min-max, the chosen test, exact P and BH-adjusted P.
# 4) Shapiro P >= 0.05 should be phrased as 'no significant departure detected',
#    not as proof that the data are normally distributed.
# =============================================================================
