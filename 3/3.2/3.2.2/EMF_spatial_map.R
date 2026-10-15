# Spatial distribution of ecosystem multifunctionality (EMF_avg).
# Source: the supplied nine-column CSV. No values or coordinates are altered.
# GitHub references reviewed on 2026-10-14:
# https://github.com/tidyverse/ggplot2
# https://github.com/davidgohel/rvg
# https://github.com/davidgohel/officer
#
# Run: Rscript EMF_spatial_map.R [input_csv] [output_directory]
# RStudio: source("EMF_spatial_map.R") after installing the packages below.
# install.packages(c("ggplot2", "rvg", "officer", "patchwork", "ragg"))
# Optional EMF_R_LIB environment variable selects an existing package library.

if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)))
}
args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1L) args[1] else
  "D:/Users/Lenovo/Desktop/\u539f\u59cb/3.2/3.2.2/\u6837\u70b9+EMF.csv"
out_dir <- if (length(args) >= 2L) args[2] else
  "D:/Users/Lenovo/Desktop/\u539f\u59cb/3.2/3.2.2/EMF_spatial_20261014"
local_lib <- Sys.getenv("EMF_R_LIB", "")
if (nzchar(local_lib)) .libPaths(c(local_lib, .libPaths()))
required <- c("ggplot2", "rvg", "officer", "patchwork", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these R packages before running: ", paste(missing, collapse = ", "))
library(ggplot2)
library(patchwork)
library(grid)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input_csv, check.names = FALSE, fileEncoding = "UTF-8-BOM")
expected <- c("Plot ID", "Region", "Longitude", "Latitude", "EMF-carbon",
              "EMF-nutrient", "EMF-enzyme", "EMF-vegetation", "EMF-avg")
if (!identical(names(d), expected)) stop("The input must have the original nine columns in the documented order.")
names(d) <- c("Plot_ID", "Region", "Longitude", "Latitude", "EMF_carbon",
              "EMF_nutrient", "EMF_enzyme", "EMF_vegetation", "EMF_avg")
stopifnot(nrow(d) == 143L, all(complete.cases(d)), !anyDuplicated(d$Plot_ID),
          all(d$Region %in% c("FT", "GQ")), all(d$Longitude >= -180 & d$Longitude <= 180),
          all(d$Latitude >= -90 & d$Latitude <= 90))
stopifnot(all(vapply(d[3:9], is.numeric, logical(1))))
if (anyDuplicated(d[c("Longitude", "Latitude")])) warning("Some sample plots share coordinates.")
average_error <- max(abs(d$EMF_avg - rowMeans(d[5:8])))
if (average_error > 1e-6) warning("The supplied EMF_avg differs from the mean of its four components.")

# Local equirectangular coordinates are used only for distances and diagnostics.
# Maps retain the original longitude and latitude. The coordinate datum was
# not supplied. No datum transformation, basemap or boundary is inferred.
earth_radius_m <- 6371008.8
local_xy <- function(z) {
  cbind((z$Longitude - mean(z$Longitude)) * pi / 180 * earth_radius_m *
          cos(mean(z$Latitude) * pi / 180),
        (z$Latitude - mean(z$Latitude)) * pi / 180 * earth_radius_m)
}
haversine <- function(lon1, lat1, lon2, lat2) {
  p <- pi / 180
  a <- sin((lat2-lat1)*p/2)^2 + cos(lat1*p)*cos(lat2*p)*sin((lon2-lon1)*p/2)^2
  2 * earth_radius_m * atan2(sqrt(a), sqrt(pmax(0, 1-a)))
}

# Exploratory IDW diagnostic only, not a fitted/validated interpolation model.
# Leave each plot out, use up to eight neighbours and inverse-square weights.
# The regional-mean baseline also excludes the held-out plot. Random/LOO
# validation can be optimistic for clustered data, so no inferential claim
# about all possible interpolation methods follows from this diagnostic.
summary_rows <- lapply(c("FT", "GQ"), function(region) {
  z <- d[d$Region == region, ]
  xy <- local_xy(z)
  distance <- as.matrix(dist(xy)); diag(distance) <- Inf
  eigenvalues <- prcomp(xy)$sdev^2
  idw <- vapply(seq_len(nrow(z)), function(i) {
    j <- order(distance[i, ])[seq_len(min(8L, nrow(z)-1L))]
    weighted.mean(z$EMF_avg[j], 1 / distance[i, j]^2)
  }, numeric(1))
  baseline <- vapply(seq_len(nrow(z)), function(i) mean(z$EMF_avg[-i]), numeric(1))
  data.frame(Region = region, n = nrow(z), mean = mean(z$EMF_avg),
             sd = sd(z$EMF_avg), median = median(z$EMF_avg),
             min = min(z$EMF_avg), max = max(z$EMF_avg),
             lon_min = min(z$Longitude), lon_max = max(z$Longitude),
             lat_min = min(z$Latitude), lat_max = max(z$Latitude),
             east_west_extent_m = diff(range(xy[, 1])),
             north_south_extent_m = diff(range(xy[, 2])),
             median_nearest_neighbour_m = median(apply(distance, 1, min)),
             principal_axis_variance_fraction = eigenvalues[1] / sum(eigenvalues),
             IDW_LOOCV_RMSE = sqrt(mean((idw-z$EMF_avg)^2)),
             regional_mean_LOOCV_RMSE = sqrt(mean((baseline-z$EMF_avg)^2)))
})
regional_summary <- do.call(rbind, summary_rows)
write.csv(regional_summary, file.path(out_dir, "EMF_regional_summary.csv"), row.names = FALSE)

# Identify spatial clusters by connected components: join Gaoqiao plots
# within 1.5 km. These labels are display groups, not ecological boundaries.
gq <- d[d$Region == "GQ", ]
distance <- as.matrix(dist(local_xy(gq)))
membership <- integer(nrow(gq)); k <- 0L
for (i in seq_len(nrow(gq))) {
  if (membership[i] != 0L) next
  k <- k + 1L; membership[i] <- k; queue <- i
  while (length(queue)) {
    j <- queue[1]; queue <- queue[-1]
    neighbours <- which(distance[j, ] <= 1500 & membership == 0L)
    membership[neighbours] <- k; queue <- c(queue, neighbours)
  }
}
cluster_centres <- aggregate(gq[c("Longitude", "Latitude")], list(old = membership), mean)
order_west_to_east <- cluster_centres$old[order(cluster_centres$Longitude)]
gq$Cluster <- match(membership, order_west_to_east)
cluster_centres$Cluster <- match(cluster_centres$old, order_west_to_east)
cluster_centres <- cluster_centres[order(cluster_centres$Cluster), ]
write.csv(gq[c("Plot_ID", "Region", "Longitude", "Latitude", "EMF_avg", "Cluster")],
          file.path(out_dir, "EMF_Gaoqiao_cluster_membership.csv"), row.names = FALSE)

# One signed, symmetric colour scale in every panel. Equal-size circles avoid
# mapping negative EMF values to bubble area or concealing sign by shifting.
colour_limit <- ceiling(max(abs(d$EMF_avg)) * 2) / 2
emf_scale <- function() scale_fill_gradient2(
  name = "EMF_avg", low = "#2166AC", mid = "#F7F7F7", high = "#D66000",
  midpoint = 0, limits = c(-colour_limit, colour_limit),
  breaks = seq(-colour_limit, colour_limit, by = 0.5),
  labels = function(x) sprintf("%.1f", x),
  # Rectangles remain individual editable objects in PowerPoint.
  guide = guide_colourbar(direction = "horizontal", title.position = "top",
                         title.hjust = 0.5, label.position = "bottom",
                         barwidth = unit(72, "mm"), barheight = unit(3.6, "mm"),
                         display = "rectangles", nbin = 80,
                         frame.colour = "#555555", ticks.colour = "#555555"))

map_bounds <- function(z, aspect = 1) {
  lat0 <- mean(z$Latitude)
  mx <- earth_radius_m*pi/180*cos(lat0*pi/180)
  my <- earth_radius_m*pi/180
  side <- max(diff(range(z$Longitude))*mx, diff(range(z$Latitude))*my*aspect, 180) * 1.22
  lon0 <- mean(range(z$Longitude)); latcentre <- mean(range(z$Latitude))
  list(x = lon0+c(-0.5,0.5)*side/mx, y = latcentre+c(-0.5,0.5)*side/(my*aspect),
       mx = mx, my = my, side = side, lat0 = lat0)
}
nice_scale <- function(max_length) {
  steps <- sort(as.vector(outer(c(1, 2, 5), 10^(-1:6))))
  max(steps[steps <= max_length])
}
map_panel <- function(z, title, subtitle, point_size = 3.2, detail = FALSE,
                      label_clusters = FALSE) {
  b <- map_bounds(z, aspect = if (detail) 1.7 else 1); dx <- diff(b$x); dy <- diff(b$y)
  scale_m <- nice_scale(b$side * 0.24)
  xs <- b$x[1]+dx*0.075; ys <- b$y[1]+dy*0.07
  # Put map furniture in the emptier corner so it does not cover observations.
  scaled_points <- cbind((z$Longitude-b$x[1])/dx, (z$Latitude-b$y[1])/dy)
  nearest_point_distance <- function(x, y) min((scaled_points[,1]-x)^2+(scaled_points[,2]-y)^2)
  north_fraction <- 0.93
  if (detail) {
    north_options <- c(0.07, 0.93)
    north_fraction <- north_options[which.max(vapply(north_options, function(x) nearest_point_distance(x,0.12+0.75), numeric(1)))]
    scale_fraction <- scale_m/(b$mx*dx)
    scale_options <- c(0.075, 0.925-scale_fraction)
    scale_start <- scale_options[which.max(vapply(scale_options, function(x) nearest_point_distance(x+scale_fraction/2,0.10), numeric(1)))]
    xs <- b$x[1]+dx*scale_start
  }
  xe <- xs+scale_m/b$mx
  label <- if (scale_m >= 1000) paste0(scale_m/1000, " km") else paste0(scale_m, " m")
  precision <- if (dx < 0.006) 4L else if (dx < 0.1) 3L else 2L
  axis_font <- if (detail) 7.2 else 9
  p <- ggplot(z, aes(Longitude, Latitude)) +
    geom_point(aes(fill = EMF_avg), shape = 21, size = point_size,
               colour = "#333333", stroke = if (detail) 0.20 else 0.28) +
    emf_scale() +
    annotate("segment", x = xs, xend = xe, y = ys, yend = ys, linewidth = 0.5) +
    annotate("segment", x = c(xs,xe), xend = c(xs,xe), y = ys-dy*0.006,
             yend = ys+dy*0.006, linewidth = 0.4) +
    annotate("text", x = (xs+xe)/2, y = ys+dy*(if (detail) 0.060 else 0.026), label = label,
             size = if (detail) 2.2 else 2.7, family = "Arial", colour = "#333333") +
    annotate("segment", x = b$x[1]+dx*north_fraction, xend = b$x[1]+dx*north_fraction,
             y = b$y[2]-dy*0.18, yend = b$y[2]-dy*0.075,
             arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = 0.35) +
    annotate("text", x = b$x[1]+dx*north_fraction, y = b$y[2]-dy*0.045, label = "N",
             size = if (detail) 2.3 else 2.8, family = "Arial") +
    coord_fixed(ratio = 1/cos(b$lat0*pi/180), xlim = b$x, ylim = b$y, expand = FALSE) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = if (detail) 2 else 3),
      labels = function(x) paste0(formatC(x, format = "f", digits = precision), if (detail) "" else "\u00b0E")) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = if (detail) 2 else 3),
      labels = function(x) paste0(formatC(x, format = "f", digits = precision), if (detail) "" else "\u00b0N")) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_family = "Arial", base_size = 10) +
    theme(panel.border = element_rect(fill = NA, colour = "#7B848A", linewidth = 0.45),
          panel.grid.major = element_line(colour = "#E3E7EA", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = if (detail) 10 else 14, face = "bold", colour = "#182F3A"),
          plot.subtitle = element_text(size = if (detail) 8 else 9.5, colour = "#49545B", margin = margin(b = 7)),
          axis.text = element_text(size = axis_font, colour = "#3F4850"),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          plot.margin = margin(5, 8, 5, 5),
          legend.position = "bottom",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          legend.margin = margin(1, 4, 1, 4))
  if (label_clusters) {
    lab <- cluster_centres
    lab$label_x <- lab$Longitude+dx*0.040
    lab$label_y <- lab$Latitude+dy*ifelse(lab$Cluster == 4, -0.035, 0.035)
    p <- p + geom_segment(data = lab,
        aes(x = Longitude, y = Latitude, xend = label_x-dx*0.006, yend = label_y),
        inherit.aes = FALSE, linewidth = 0.2, colour = "#87939A") +
      geom_text(data = lab,
        aes(x = label_x, y = label_y, label = paste0("C", Cluster)),
        inherit.aes = FALSE, size = 2.7, family = "Arial", fontface = "bold", colour = "#152E3D")
  }
  p
}

ft <- d[d$Region == "FT", ]
panel_subtitle <- function(z) sprintf("n = %d     Mean = %.3f     Range = %.3f to %.3f",
                                      nrow(z), mean(z$EMF_avg), min(z$EMF_avg), max(z$EMF_avg))
main_figure <- (map_panel(ft, "a  Futian", panel_subtitle(ft), 3.2) +
               map_panel(gq, "b  Gaoqiao", panel_subtitle(gq), 3.2, label_clusters = TRUE) +
               plot_layout(ncol = 2, guides = "collect")) +
  plot_annotation(
    title = "Spatial distribution of ecosystem multifunctionality",
    subtitle = "143 sample plots in Futian and Gaoqiao",
    caption = paste0("Circles show observed EMF_avg. Blue indicates lower values and orange indicates higher values.\n",
                     "Panels use different map scales. C1\u2013C", k, " identify the Gaoqiao detail panels.\n",
                     "Coordinates as supplied, datum unspecified. Scale bars are approximate. Overlapping plots retain their locations."),
    theme = theme(plot.title = element_text(family = "Arial", size = 18, face = "bold", colour = "#182F3A"),
                  plot.subtitle = element_text(family = "Arial", size = 11, colour = "#53616B"),
                  plot.caption = element_text(family = "Arial", size = 8, colour = "#53616B", hjust = 0),
                  plot.margin = margin(13, 15, 11, 15)))
main_figure <- main_figure & theme(legend.position = "bottom")

detail_panels <- lapply(seq_len(k), function(i) {
  z <- gq[gq$Cluster == i, ]
  map_panel(z, paste0("C", i), sprintf("n = %d   Mean = %.3f", nrow(z), mean(z$EMF_avg)),
            point_size = 2.7, detail = TRUE)
})
detail_figure <- wrap_plots(detail_panels, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Gaoqiao: spatial variation within sampling clusters",
    subtitle = "120 plots. Cluster labels match the regional overview",
    caption = "The same EMF_avg colour scale applies throughout. Each panel has its own map scale. Axes show longitude (\u00b0E) and latitude (\u00b0N).\nGroups link plots within 1.5 km for display only. Coordinates as supplied, datum unspecified. Scale bars are approximate.",
    theme = theme(plot.title = element_text(family = "Arial", size = 17, face = "bold", colour = "#182F3A"),
                  plot.subtitle = element_text(family = "Arial", size = 10, colour = "#53616B"),
                  plot.caption = element_text(family = "Arial", size = 7.8, colour = "#53616B", hjust = 0),
                  plot.margin = margin(11, 13, 10, 13)))
detail_figure <- detail_figure & theme(legend.position = "bottom")

# 600-dpi, losslessly compressed TIFF, 6000 x 4500 pixels.
# The corresponding PNG is an optional convenient viewing copy.
export_figure <- function(plot, stem) {
  ragg::agg_tiff(file.path(out_dir, paste0(stem, ".tiff")), width = 10, height = 7.5,
                 units = "in", res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
  ragg::agg_png(file.path(out_dir, paste0(stem, ".png")), width = 10, height = 7.5,
                units = "in", res = 180, background = "white")
  print(plot); dev.off()
}
export_figure(main_figure, "EMF_spatial_distribution")
export_figure(detail_figure, "EMF_Gaoqiao_details")

# DrawingML vectors: points, text, gridlines and colour-bar rectangles are
# editable PowerPoint objects. Ungroup the figure to edit individual elements.
# These are editable vector graphics, not linked PowerPoint chart workbooks.
doc <- officer::read_pptx()
doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
doc <- officer::ph_with(doc, rvg::dml(ggobj = main_figure, editable = TRUE),
                       location = officer::ph_location_fullsize())
doc <- officer::set_notes(doc, value = paste0(
  "Source: supplied sample-plot CSV, 143 records, 23 Futian and 120 Gaoqiao. ",
  "Equal-size circles preserve the signed EMF_avg values in a shared colour scale. ",
  "No interpolation, jitter, coordinate transformation or inferred boundary. ",
  "Scale bars use a spherical local distance approximation. ",
  "References: https://github.com/tidyverse/ggplot2 ; https://github.com/davidgohel/rvg ; https://github.com/davidgohel/officer ."),
  location = officer::notes_location_type("body"))
doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
doc <- officer::ph_with(doc, rvg::dml(ggobj = detail_figure, editable = TRUE),
                       location = officer::ph_location_fullsize())
doc <- officer::set_notes(doc, value = paste0(
  "Source: the 120 Gaoqiao records in the supplied CSV. ",
  "Display clusters are connected components with a 1500-m neighbour threshold, ",
  "ordered by centroid longitude from west to east. They are not ecological or administrative boundaries. ",
  "All 120 records appear once in the detail atlas. Every panel has its own scale bar. ",
  "References: https://github.com/tidyverse/ggplot2 ; https://github.com/davidgohel/rvg ; https://github.com/davidgohel/officer ."),
  location = officer::notes_location_type("body"))
print(doc, target = file.path(out_dir, "EMF_spatial_distribution_editable.pptx"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "EMF_R_sessionInfo.txt"))
message("Finished. TIFF and editable PPTX files are in: ", normalizePath(out_dir, winslash = "/"))
