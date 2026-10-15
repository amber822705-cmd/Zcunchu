# EMF_nutrient 环境因子随机森林分析与期刊图件
# 输入：模块.xlsx / 自变量 / A:X
# 响应：X 列 EMF-vegetation；预测：C:S 的 17 个环境变量。
# 排除 Plot ID、Region 和 T:W 的四个 EMF 功能分数。
# 输出 TIFF、可编辑 PPTX 和绘图数据汇总 XLSX。模型参数与软件版本也输出至 R 控制台。
# 运行：在 RStudio 打开本文件，点击 Source；或 Rscript EMF_RF_analysis.R。
# 可选命令行：Rscript EMF_RF_analysis.R 输入文件 输出目录
# 首次运行若缺少包：install.packages(c("readxl", "randomForest", "ragg",
#   "officer", "rvg", "xml2", "zip", "openxlsx"), repos = "https://cloud.r-project.org")
# 方法：区域分层 5 折外层、4 折内层嵌套交叉验证；100 次验证集置换/变量/折。
# 图中 SD 为折均值之间的样本标准差，不是置信区间，也不是显著性检验。
# 重要性反映预测关联。负值表示打乱变量后预测略有改善，不代表负向生态效应。
# R randomForest 和其他随机森林实现的数值与排序可能不同。
# 文档：https://cran.r-project.org/web/packages/randomForest/randomForest.pdf
# 编辑图形：https://davidgohel.github.io/rvg/reference/ph_with.dml.html

if (.Platform$OS.type == "windows") {
  user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library",
                       paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = "."))
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
}
required <- c("readxl", "randomForest", "ragg", "officer", "rvg", "xml2", "zip", "openxlsx")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "),
       "。请先运行 install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), ")).")
}
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[1L] else
  "D:/Users/Lenovo/Desktop/原始/3.3/3.3.1/模块.xlsx"
output_dir <- if (length(args) >= 2L) args[2L] else
  "D:/Users/Lenovo/Desktop/原始/3.3/3.3.4"
SEED <- 20261014L
OUTER_K <- 5L
INNER_K <- 4L
TUNE_TREES <- 300L
FINAL_TREES <- 1000L
N_PERM <- 100L
FIG_WIDTH <- 7.3
FIG_HEIGHT <- 7.6
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(SEED)
if (!file.exists(input_file)) stop("找不到输入文件：", input_file)
dat <- as.data.frame(readxl::read_excel(input_file, sheet = "自变量", range = readxl::cell_cols("A:X"),
                                      .name_repair = "minimal"), check.names = FALSE)
names(dat) <- trimws(names(dat))
dat <- dat[!is.na(dat[[1L]]) & nzchar(trimws(as.character(dat[[1L]]))), , drop = FALSE]
if (ncol(dat) != 24L || !identical(names(dat)[c(1L, 2L, 24L)],
                                c("Plot ID", "Region", "EMF-avg"))) {
  stop("工作表应为 A:X 共24列，A=Plot ID、B=Region、X=EMF-avg。")
}
if (anyDuplicated(dat[[1L]])) stop("Plot ID 存在重复，请先确认样点对应关系。")
if (anyNA(dat[[2L]])) stop("Region 存在缺失，无法进行区域分层验证。")
# 固定使用12个环境预测因子（根据A:X表结构）
# E 土壤含水率
# H Salinity
# I pH
# J C:N比
# L N:P比
# M SOC
# N TN
# O TP
# P MAT_C
# Q MAP_mm_yr
# R Distance_coast_m
# S Distance_estuary_model_m

selected_cols <- c(5, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19)

feature_names <- names(dat)[selected_cols]

X <- as.data.frame(dat[, selected_cols, drop = FALSE],
                   check.names = FALSE)
for (j in seq_along(X)) {
  if (!is.numeric(X[[j]])) stop("自变量不是数值型：", feature_names[j])
}
y <- dat[["EMF-vegetation"]]
region <- as.character(dat[["Region"]])
if (!is.numeric(y) || any(!is.finite(y)) || any(!is.finite(as.matrix(X)))) {
  stop("响应或自变量存在缺失、非数值或非有限值；脚本不会自动填补或删除。")
}
if (any(table(region) < OUTER_K)) stop("每个 Region 至少需要5个样点。")
if (var(y) <= 0) stop("EMF-vegetation 无变异，无法计算 R²。")
n <- nrow(X)
p <- ncol(X)
# 使用合法内部变量名，避免中文、冒号与单位符号被公式解析。
# 图件和结果表仍保留原始变量名，映射顺序完全一致。
names(X) <- sprintf("X%02d", seq_len(p))
labels_en <- c("Soil water content", "Salinity", "pH",
  "C:N ratio", "N:P ratio", "SOC", "TN", "TP",
  "Mean annual temperature", "Mean annual precipitation",
  "Distance to coast", "Distance to estuary")

short_labels <- c("Water content", "Salinity", "pH",
  "C:N", "N:P", "SOC", "TN", "TP",
  "MAT", "MAP", "Coast distance", "Estuary distance")

score <- function(observed, predicted) {
  c(R2 = 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2),
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted)))
}
stratified_folds <- function(groups, k, seed) {
  set.seed(seed)
  folds <- integer(length(groups))
  for (r in sort(unique(groups))) {
    ids <- which(groups == r)
    if (length(ids) < k) stop("分层组样本量少于折数：", r)
    folds[ids] <- sample(rep(seq_len(k), length.out = length(ids)))
  }
  folds
}
fit_rf <- function(x, response, parameter, trees, seed) {
  set.seed(seed)
  randomForest::randomForest(x = x, y = response, ntree = trees,
    mtry = as.integer(parameter$mtry), nodesize = as.integer(parameter$nodesize),
    replace = TRUE, sampsize = nrow(x), importance = FALSE, keep.forest = TRUE)
}
grid_parameters <- expand.grid(mtry = unique(c(floor(p / 3), floor(p / 2), p)),
                               nodesize = c(1L, 3L, 5L))
tune_rf <- function(x, response, groups, k, seed) {
  folds <- stratified_folds(groups, k, seed)
  trial_mse <- matrix(NA_real_, nrow(grid_parameters), k)
  for (g in seq_len(nrow(grid_parameters))) {
    for (f in seq_len(k)) {
      tr <- folds != f
      te <- folds == f
      model <- fit_rf(x[tr, , drop = FALSE], response[tr], grid_parameters[g, ],
                      TUNE_TREES, seed + g * 100L + f)
      pr <- predict(model, x[te, , drop = FALSE])
      trial_mse[g, f] <- mean((response[te] - pr)^2)
    }
  }
  result <- cbind(grid_parameters, mean_CV_MSE = rowMeans(trial_mse),
                  fold_MSE_SD = apply(trial_mse, 1L, sd))
  result$rank <- rank(result$mean_CV_MSE, ties.method = "min")
  list(best = grid_parameters[which.min(result$mean_CV_MSE), , drop = FALSE],
       results = result)
}
permutation_scores <- function(model, x, response, columns, seed, within = NULL) {
  # 每次置换仅打乱指定列；联合组中的列同步打乱以保留组内关系。
  # 批量 predict 与逐次预测完全等价，用于加快计算。
  set.seed(seed)
  nv <- nrow(x)
  permutations <- replicate(N_PERM, {
    ids <- seq_len(nv)
    if (is.null(within)) ids <- sample(ids) else {
      for (r in sort(unique(within))) {
        ix <- which(within == r)
        ids[ix] <- sample(ix)
      }
    }
    ids
  })
  baseline_mse <- mean((response - predict(model, x))^2)
  response_variance <- mean((response - mean(response))^2)
  values <- matrix(NA_real_, length(columns), N_PERM)
  for (j in seq_along(columns)) {
    shuffled <- x[rep(seq_len(nv), times = N_PERM), , drop = FALSE]
    for (col in columns[[j]]) shuffled[[col]] <- x[[col]][as.vector(permutations)]
    predictions <- matrix(predict(model, shuffled), nrow = nv, ncol = N_PERM)
    values[j, ] <- (colMeans((predictions - response)^2) - baseline_mse) / response_variance
  }
  single <- x
  for (col in columns[[1L]]) single[[col]] <- x[[col]][permutations[, 1L]]
  check <- (mean((predict(model, single) - response)^2) - baseline_mse) / response_variance
  stopifnot(abs(check - values[1L, 1L]) < 1e-10)
  values
}

outer_fold <- stratified_folds(region, OUTER_K, SEED)
oof_prediction <- rep(NA_real_, n)
baseline_oof <- rep(NA_real_, n)
fold_importance <- matrix(NA_real_, p, OUTER_K)
within_importance <- matrix(NA_real_, p, OUTER_K)
joint_importance <- matrix(NA_real_, nrow = length(joint_columns), ncol = OUTER_K)
fold_performance <- vector("list", OUTER_K)
tuning_records <- vector("list", OUTER_K + 1L)
permutation_repeats <- vector("list", OUTER_K)
joint_columns <- list(
  c(1L, 2L, 3L),
  c(4L, 5L, 6L, 7L, 8L),
  c(9L, 10L),
  c(11L, 12L)
)
cat("有效样点：", n, "；预测变量：", p, "\n")
print(table(region))
for (f in seq_len(OUTER_K)) {
  train <- outer_fold != f
  test <- outer_fold == f
  tuning <- tune_rf(X[train, , drop = FALSE], y[train], region[train], INNER_K, SEED + f)
  model <- fit_rf(X[train, , drop = FALSE], y[train], tuning$best, FINAL_TREES, SEED + f)
  prediction <- predict(model, X[test, , drop = FALSE])
  oof_prediction[test] <- prediction
  baseline_oof[test] <- mean(y[train])
  perf <- score(y[test], prediction)
  fold_performance[[f]] <- data.frame(fold = f, n_train = sum(train), n_test = sum(test),
    R2 = perf["R2"], RMSE = perf["RMSE"], MAE = perf["MAE"],
    mtry = tuning$best$mtry, nodesize = tuning$best$nodesize, row.names = NULL)
  tuning_records[[f]] <- cbind(stage = paste0("outer_", f), tuning$results)
  values <- permutation_scores(model, X[test, , drop = FALSE], y[test],
    as.list(seq_len(p)), SEED + 1000L + f)
  fold_importance[, f] <- rowMeans(values)
  permutation_repeats[[f]] <- values
  within <- permutation_scores(model, X[test, , drop = FALSE], y[test],
    as.list(seq_len(p)), SEED + 2000L + f, region[test])
  within_importance[, f] <- rowMeans(within)
  joint_importance[, f] <- rowMeans(permutation_scores(model, X[test, , drop = FALSE],
    y[test], joint_columns, SEED + 3000L + f))
  cat(sprintf("外层第 %d/%d 折完成，验证 R² = %.4f\n", f, OUTER_K, perf["R2"]))
  flush.console()
}
stopifnot(all(is.finite(oof_prediction)))
final_tuning <- tune_rf(X, y, region, OUTER_K, SEED + 100L)
tuning_records[[OUTER_K + 1L]] <- cbind(stage = "full_data_final", final_tuning$results)
final_model <- fit_rf(X, y, final_tuning$best, FINAL_TREES, SEED)
cv_metrics <- score(y, oof_prediction)
oob_metrics <- score(y, final_model$predicted)
train_metrics <- score(y, predict(final_model, X))
fold_ranks <- apply(fold_importance, 2L, function(v) rank(-v, ties.method = "min"))
mdi <- as.numeric(randomForest::importance(final_model, type = 2L, scale = FALSE))
importance_table <- data.frame(variable = feature_names, label = labels_en, short_label = short_labels,
  mean_delta_R2 = rowMeans(fold_importance), fold_SD = apply(fold_importance, 1L, sd),
  top3_folds = rowSums(fold_ranks <= 3L),
  within_region_delta_R2 = rowMeans(within_importance),
  within_region_SD = apply(within_importance, 1L, sd), MDI_fraction = mdi / sum(mdi),
  stringsAsFactors = FALSE)
importance_table <- importance_table[order(-importance_table$mean_delta_R2), ]
rownames(importance_table) <- NULL
importance_table$rank <- seq_len(p)
stopifnot(abs(sum(importance_table$MDI_fraction) - 1) < 1e-10)
joint_table <- data.frame(group = c("土壤环境（水分/盐度/pH）",
  "资源化学计量（C:N/N:P/SOC/TN/TP）",
  "气候（MAT/MAP）",
  "空间距离（距海/距河口）"),
  mean_delta_R2 = rowMeans(joint_importance), fold_SD = apply(joint_importance, 1L, sd))

# 结果留在 R 会话供进一步分析；绘图数据另存为 XLSX，不另存 CSV、模型或说明文件。
performance_table <- rbind(Nested_OOF = cv_metrics, OOB = oob_metrics,
                          Training = train_metrics, Mean_baseline_OOF = score(y, baseline_oof))
description_table <- do.call(rbind, lapply(c("All", sort(unique(region))), function(r) {
  ids <- if (r == "All") seq_len(n) else which(region == r)
  z <- cbind(EMF_nutrient = y[ids], dat[ids, selected_cols, drop = FALSE])
  do.call(rbind, lapply(seq_along(z), function(j) data.frame(region = r,
    variable = names(z)[j], N = length(z[[j]]), mean = mean(z[[j]]), SD = sd(z[[j]]),
    min = min(z[[j]]), median = median(z[[j]]), max = max(z[[j]]))))
}))
transfer_table <- do.call(rbind, lapply(sort(unique(region)), function(held) {
  tr <- region != held
  te <- region == held
  transfer <- fit_rf(X[tr, , drop = FALSE], y[tr], data.frame(mtry = floor(p/2), nodesize = 3L),
                     FINAL_TREES, SEED + 4000L)
  s <- score(y[te], predict(transfer, X[te, , drop = FALSE]))
  data.frame(train_region = paste(unique(region[tr]), collapse = "+"), test_region = held,
             N = sum(te), R2 = s["R2"], RMSE = s["RMSE"], MAE = s["MAE"], row.names = NULL)
}))
cat("\n模型性能\n"); print(round(performance_table, 5))
cat("\n最终参数\n"); print(final_tuning$best)
cat("\n变量重要性（降序，前三名为候选核心因子）\n"); print(importance_table, digits = 5, row.names = FALSE)
cat("\n逐折模型性能\n"); print(do.call(rbind, fold_performance), row.names = FALSE)
cat("\n联合组敏感性（不具有可加性）\n"); print(joint_table, row.names = FALSE)
cat("\n跨区域预测（未在验证区域调参）\n"); print(transfer_table, row.names = FALSE)
cat("\n描述统计\n"); print(description_table, row.names = FALSE)

draw_importance <- function() {
  # 全部使用 grid 文字和几何对象；rvg 将其转为 PowerPoint 可编辑对象。
  grid::grid.newpage()
  ink <- "#172932"; muted <- "#536978"; dark <- "#1F4E67"
  text <- function(label, x, y, size = 9, col = ink, face = "plain", just = "left") {
    grid::grid.text(label, x, y, just = just,
      gp = grid::gpar(fontfamily = "Arial", fontsize = size, col = col, fontface = face))
  }
  segment <- function(x0, y0, x1, y1, col = "#40535F", width = 0.7) {
    # 用极细的原生矩形表示水平/垂直线，便于 PowerPoint 和预览器一致渲染。
    if (abs(y1 - y0) < 1e-12) {
      grid::grid.rect(x = (x0 + x1)/2, y = y0, width = abs(x1 - x0),
        height = grid::unit(width * 0.75, "points"), gp = grid::gpar(fill = col, col = NA))
    } else if (abs(x1 - x0) < 1e-12) {
      grid::grid.rect(x = x0, y = (y0 + y1)/2, width = grid::unit(width * 0.75, "points"),
        height = abs(y1 - y0), gp = grid::gpar(fill = col, col = NA))
    } else {
      grid::grid.lines(x = c(x0, x1), y = c(y0, y1), gp = grid::gpar(col = col, lwd = width))
    }
  }
  means <- importance_table$mean_delta_R2
  sd <- importance_table$fold_SD
  low <- min(0, means - sd)
  high <- max(0, means + sd)
  span <- high - low
  limits <- c(low - 0.02*span, high + 0.21*span)
  plot_left <- 0.35; plot_right <- 0.94
  map_x <- function(x) plot_left + (x - limits[1]) / diff(limits) * (plot_right - plot_left)
  ys <- seq(0.803, 0.190, length.out = p)
  row_height <- abs(diff(ys)[1])
  bar_height <- row_height * 0.67
  zero <- map_x(0)
  palette <- grDevices::colorRampPalette(c(dark, "#C1D1DA"))(p)
  # 前三行浅色背景；同一蓝灰色系的降序渐变。
  grid::grid.rect(x = (plot_left + plot_right)/2, y = mean(ys[1:3]),
    width = plot_right - plot_left, height = row_height * 3,
    gp = grid::gpar(fill = "#F1F5F7", col = NA))
  ticks <- pretty(limits, n = 5)
  ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
  for (t in ticks) {
    segment(map_x(t), 0.164, map_x(t), ys[1] + row_height/2, "#E3EAEE", 0.6)
    segment(map_x(t), 0.164, map_x(t), 0.159, "#AAB5BD", 0.6)
    text(sprintf("%.2f", t), map_x(t), 0.144, size = 8, col = muted, just = "center")
  }
  segment(zero, 0.164, zero, ys[1] + row_height/2, "#8B9BA5", 0.7)
  for (i in seq_len(p)) {
    endpoint <- map_x(means[i])
    grid::grid.rect(x = (zero + endpoint)/2, y = ys[i],
      width = max(abs(endpoint - zero), 0.00001), height = bar_height,
      gp = grid::gpar(fill = palette[i], col = NA))
    a <- map_x(means[i] - sd[i]); b <- map_x(means[i] + sd[i])
    segment(a, ys[i], b, ys[i])
    segment(a, ys[i] - 0.003, a, ys[i] + 0.003)
    segment(b, ys[i] - 0.003, b, ys[i] + 0.003)
    font <- if (i <= 3L) "bold" else "plain"
    label_col <- if (i <= 3L) dark else ink
    text(sprintf("%02d", i), 0.067, ys[i], 8.8, if (i <= 3L) dark else "#8A9CA7", font, "center")
    text(importance_table$label[i], 0.334, ys[i], 9, label_col, font, "right")
    text(sprintf("%.4f", means[i]), b + 0.010, ys[i], 8.2,
         if (i <= 3L) dark else muted, font)
  }
  segment(plot_left, 0.164, plot_right, 0.164, "#AAB5BD", 0.6)
  text("Environmental predictors of EMF_nutrient", 0.055, 0.957, 14.2, face = "bold")
  text("Random forest regression", 0.055, 0.921, 10.3, col = muted)
  text(sprintf("n = %d     %d predictors     Nested CV R² = %.3f     RMSE = %.3f",
               n, p, cv_metrics["R2"], cv_metrics["RMSE"]), 0.055, 0.882, 8.8, muted)
  text("RANK", 0.055, 0.840, 7.6, "#768B97")
  text("Validation permutation importance (mean ΔR²)",
       (plot_left + plot_right)/2, 0.118, 9, ink, just = "center")
  leading <- paste(sprintf("%d %s", 1:3, importance_table$short_label[1:3]), collapse = ";  ")
  text(paste0("Leading predictors: ", leading), 0.055, 0.082, 8.7, dark, "bold")
  text("Bars: mean of five validation folds. Whiskers: ±1 fold SD. 100 permutations per feature per fold.",
       0.055, 0.048, 7.1, "#6C7E88")
  text("Importance describes predictive association. Correlated variables can share information.",
       0.055, 0.028, 7.1, "#6C7E88")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
tiff_path <- file.path(output_dir, "EMF-vegetation_RF_importance.tiff")
pptx_path <- file.path(output_dir, "EMF-vegetation_RF_importance_editable.pptx")
ragg::agg_tiff(tiff_path, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in",
               res = 600, background = "white", compression = "lzw")
tryCatch(draw_importance(), finally = grDevices::dev.off())

# rvg 的 DrawingML 为原生文字、条形、线段；编辑时在 PowerPoint 中取消组合。
deck <- officer::read_pptx()
layouts <- officer::layout_summary(deck)
blank_id <- which(layouts$layout == "Blank")[1L]
if (is.na(blank_id)) stop("PowerPoint 模板中未找到 Blank 版式。")
deck <- officer::add_slide(deck, layout = layouts$layout[blank_id], master = layouts$master[blank_id])
deck <- officer::ph_with(deck, value = rvg::dml(code = draw_importance(),
  bg = "white", fonts = list(sans = "Arial", serif = "Arial", mono = "Arial"), editable = TRUE),
  location = officer::ph_location(left = 0, top = 0, width = FIG_WIDTH, height = FIG_HEIGHT))
notes <- paste("Input:", input_file, "/ 自变量 / A:X.",
  "Response = EMF-vegetation. Predictors = C:S only. Region stratifies validation.",
  "Outer 5 folds, inner 4 folds, 300 tuning trees, 1000 validation/final trees.",
  "100 validation-set permutations per feature per fold. Error bars are fold SD.",
  "https://cran.r-project.org/web/packages/randomForest/randomForest.pdf",
  "https://davidgohel.github.io/rvg/reference/ph_with.dml.html")
deck <- officer::set_notes(deck, value = notes, location = officer::notes_location_type("body"))
print(deck, target = pptx_path)

# 将幻灯片画布设为与 TIFF 相同的 7.3 x 7.6 英寸，无缩放、无裁切。
# officer::slide_size() 仅用于读取尺寸，因此在输出包中修改标准 OOXML 尺寸。
unpack_dir <- tempfile("emf_pptx_")
dir.create(unpack_dir)
utils::unzip(pptx_path, exdir = unpack_dir)
presentation_xml <- file.path(unpack_dir, "ppt", "presentation.xml")
doc <- xml2::read_xml(presentation_xml)
ns <- xml2::xml_ns(doc)
size_node <- xml2::xml_find_first(doc, ".//p:sldSz", ns)
xml2::xml_set_attr(size_node, "cx", as.character(round(FIG_WIDTH * 914400)))
xml2::xml_set_attr(size_node, "cy", as.character(round(FIG_HEIGHT * 914400)))
xml2::xml_set_attr(size_node, "type", "custom")
xml2::write_xml(doc, presentation_xml)
zip_files <- list.files(unpack_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
zip_files <- zip_files[!file.info(file.path(unpack_dir, zip_files))$isdir]
zip::zipr(zipfile = pptx_path, files = zip_files, root = unpack_dir,
          mode = "mirror", include_directories = FALSE)
unlink(unpack_dir, recursive = TRUE)

# 导出实际绘图使用的均值、SD、误差线端点、标签，以及其逐折/逐次来源。
export_plot_tables <- function(path) {
  ids <- match(importance_table$variable, feature_names)
  plot_data <- data.frame(
    排名 = importance_table$rank, 原始变量名 = importance_table$variable,
    图中英文标签 = importance_table$label, 重要性均值 = importance_table$mean_delta_R2,
    折间SD = importance_table$fold_SD,
    误差线下限 = importance_table$mean_delta_R2 - importance_table$fold_SD,
    误差线上限 = importance_table$mean_delta_R2 + importance_table$fold_SD,
    图中显示数值 = sprintf("%.4f", importance_table$mean_delta_R2),
    前三名折数 = importance_table$top3_folds,
    图中突出显示 = ifelse(importance_table$rank <= 3L, "前三候选因子", ""),
    stringsAsFactors = FALSE, check.names = FALSE)
  fold_data <- data.frame(排名 = importance_table$rank, 原始变量名 = importance_table$variable,
                         fold_importance[ids, , drop = FALSE], check.names = FALSE)
  names(fold_data)[3:7] <- paste0("外层折", 1:5, "均值")
  fold_data$五折平均 <- importance_table$mean_delta_R2
  fold_data$五折SD <- importance_table$fold_SD
  repeat_data <- do.call(rbind, lapply(seq_along(permutation_repeats), function(f) {
    do.call(rbind, lapply(seq_len(nrow(importance_table)), function(i) {
      values <- permutation_repeats[[f]][ids[i], ]
      data.frame(外层折 = f, 绘图排名 = i, 原始变量名 = importance_table$variable[i],
                 置换序号 = seq_along(values), delta_R2 = as.numeric(values),
                 stringsAsFactors = FALSE, check.names = FALSE)
    }))
  }))
  openxlsx::write.xlsx(list(绘图汇总 = plot_data, 逐折重要性 = fold_data, 置换重复 = repeat_data),
    file = path, overwrite = TRUE, colWidths = "auto", rowNames = FALSE,
    headerStyle = openxlsx::createStyle(fontName = "Microsoft YaHei", fontSize = 10,
      fontColour = "#FFFFFF", fgFill = "#214E66", textDecoration = "bold"))
  # 清除部分openxlsx版本遗留的空绘图/打印关系，避免Excel修复提示。
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  table_dir <- tempfile("emf_xlsx_")
  dir.create(table_dir)
  utils::unzip(path, exdir = table_dir)
  sheets <- list.files(file.path(table_dir, "xl", "worksheets"),
                       pattern = "^sheet[0-9]+\\.xml$", full.names = TRUE)
  for (sheet in sheets) {
    xml <- xml2::read_xml(sheet)
    xml2::xml_remove(xml2::xml_find_all(xml,
      "//*[local-name()='drawing' or local-name()='legacyDrawing' or local-name()='legacyDrawingHF' or local-name()='pageSetup']"))
    xml2::write_xml(xml, sheet)
  }
  relations <- list.files(table_dir, recursive = TRUE, all.files = TRUE,
                          pattern = "\\.rels$", full.names = TRUE)
  for (relation in relations) {
    xml <- xml2::read_xml(relation)
    nodes <- xml2::xml_find_all(xml, "//*[local-name()='Relationship']")
    for (node in nodes) {
      mode <- xml2::xml_attr(node, "TargetMode")
      if (!is.na(mode) && mode == "External") next
      target <- xml2::xml_attr(node, "Target")
      physical <- if (startsWith(target, "/")) file.path(table_dir, sub("^/", "", target)) else
        file.path(dirname(dirname(relation)), target)
      if (!file.exists(physical)) xml2::xml_remove(node)
    }
    xml2::write_xml(xml, relation)
  }
  types <- file.path(table_dir, "[Content_Types].xml")
  xml <- xml2::read_xml(types)
  for (node in xml2::xml_find_all(xml, "//*[local-name()='Override']")) {
    part <- sub("^/", "", xml2::xml_attr(node, "PartName"))
    if (!file.exists(file.path(table_dir, part))) xml2::xml_remove(node)
  }
  xml2::write_xml(xml, types)
  files <- list.files(table_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(file.path(table_dir, files))$isdir]
  zip::zipr(path, files = files, root = table_dir, mode = "mirror", include_directories = FALSE)
  unlink(table_dir, recursive = TRUE)
}
xlsx_path <- file.path(output_dir, "EMF-vegetation_RF_绘图数据汇总.xlsx")
export_plot_tables(xlsx_path)
cat("\n已输出：\n", tiff_path, "\n", pptx_path, "\n", xlsx_path, "\n")
cat("PPTX 编辑：选择图形，右键 → 组合 → 取消组合；必要时再次取消组合。\n")
cat("\n软件版本\n")
for (pkg in required) cat(pkg, as.character(getNamespaceVersion(pkg)), "\n")
print(sessionInfo())
