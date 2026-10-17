#!/usr/bin/env Rscript

# Plot the five lowest and the highest snowfall winters as 3D perspective maps:
# the unmasked snowfall field is draped over the topography.
# Usage: Rscript 04_plot_field_shaded_multimodel_3D.R [snowfall_sorted.txt] [out.png]

suppressPackageStartupMessages(library(ncdf4))

BASEDIR <- "/highres_nobackup/svenk/SEKEstorylines_snow"
args <- commandArgs(trailingOnly = TRUE)
ranking_file <- if (length(args) >= 1L) args[1] else file.path(
  BASEDIR, "RESULTS", "snowfall_sorted.txt"
)
out_file <- if (length(args) >= 2L) args[2] else file.path(
  BASEDIR, "PLOTS", "snowfall_lowest5_unmasked_3D.png"
)
topo_file <- file.path(BASEDIR, "META", "topo.swiss1_ch01r.swiss.lv95.nc")

SHADE_STRENGTH <- 0.7
SHADE_GREY <- 0.35
GRID_STEP <- 2L        # use every n-th grid cell in the 3D surface
VERTICAL_EXAG <- 4     # vertical exaggeration of the terrain
VIEW_THETA <- 0        # azimuth of the viewer (0 = looking from the south, north at the back)
VIEW_PHI <- 45         # elevation angle of the viewer

results <- read.delim(ranking_file, stringsAsFactors = FALSE)
required_columns <- c("snowfall", "year", "model")
if (!all(required_columns %in% names(results))) {
  stop("Ranking file must contain snowfall, year, and model columns.", call. = FALSE)
}
if (nrow(results) < 6L) {
  stop("Ranking file contains fewer than six winters.", call. = FALSE)
}
results$snowfall <- as.numeric(results$snowfall)
results$year <- as.integer(as.character(results$year))
if (any(!is.finite(results$snowfall)) || anyNA(results$year) ||
    anyNA(results$model) || any(!nzchar(results$model))) {
  stop("Ranking file contains invalid snowfall, year, or model values.", call. = FALSE)
}
lowest5 <- results[order(results$snowfall), ][seq_len(5L), ]
maximum <- results[which.max(results$snowfall), ]
selected <- rbind(lowest5, maximum)

topo_nc <- nc_open(topo_file)
height <- ncvar_get(topo_nc, "height", start = c(1, 1, 1), count = c(-1, -1, 1))
easting <- ncvar_get(topo_nc, "E")
northing <- ncvar_get(topo_nc, "N")
nc_close(topo_nc)
height[!is.finite(height)] <- 0

read_ranked_field <- function(model, year) {
  field_file <- file.path(
    BASEDIR, "RESULTS",
    sprintf("result_field_%s_gwl3.0_nonmasked.nc", model)
  )
  if (!file.exists(field_file)) {
    stop(sprintf("Missing unmasked field file: %s", field_file), call. = FALSE)
  }

  nc <- nc_open(field_file)
  on.exit(nc_close(nc), add = TRUE)
  time <- ncvar_get(nc, "time")
  units <- ncatt_get(nc, "time", "units")
  calendar <- ncatt_get(nc, "time", "calendar")
  if (!units$hasatt || !grepl("^days since 1-1-1", units$value) ||
      !calendar$hasatt || !calendar$value %in% c("365_day", "noleap")) {
    stop(
      sprintf("Unsupported time coordinate in %s.", field_file),
      call. = FALSE
    )
  }
  time_years <- floor(time / 365) + 1L
  time_index <- which(time_years == year)
  if (length(time_index) != 1L) {
    stop(
      sprintf("Expected one timestep for year %s in %s.", year, field_file),
      call. = FALSE
    )
  }

  field_easting <- ncvar_get(nc, "E")
  field_northing <- ncvar_get(nc, "N")
  if (!identical(dim(height), c(length(field_easting), length(field_northing))) ||
      !isTRUE(all.equal(easting, field_easting)) ||
      !isTRUE(all.equal(northing, field_northing))) {
    stop(sprintf("Grid does not match topography in %s.", field_file), call. = FALSE)
  }
  field <- ncvar_get(
    nc, "pr",
    start = c(1L, 1L, time_index),
    count = c(-1L, -1L, 1L)
  )
  list(field = field, model = model, year = year, file = field_file)
}

fields <- lapply(seq_len(nrow(selected)), function(i) {
  read_ranked_field(selected$model[i], selected$year[i])
})

res <- mean(diff(easting))
nx <- nrow(height)
ny <- ncol(height)
x_prev <- c(1, seq_len(nx - 1L))
x_next <- c(seq_len(nx)[-1L], nx)
y_prev <- c(1, seq_len(ny - 1L))
y_next <- c(seq_len(ny)[-1L], ny)
dzdx <- (height[x_next, ] - height[x_prev, ]) / (2 * res)
dzdy <- (height[, y_next] - height[, y_prev]) / (2 * res)
slope <- atan(sqrt(dzdx^2 + dzdy^2))
aspect <- atan2(dzdy, -dzdx)

# Sun from the north-west, 40 degrees above the horizon.
azimuth <- (360 - 315 + 90) * pi / 180
altitude <- 40 * pi / 180
shade <- cos(pi / 2 - altitude) * cos(slope) +
  sin(pi / 2 - altitude) * sin(slope) * cos(azimuth - aspect)
shade <- pmin(pmax(shade, 0), 1)
flat_shade <- cos(pi / 2 - altitude)
grey_weight <- SHADE_STRENGTH *
  pmin(pmax((flat_shade - shade) / flat_shade, 0), 1)
if (max(grey_weight) > 0) {
  grey_weight <- grey_weight / max(grey_weight) * SHADE_STRENGTH
}

field_values <- unlist(lapply(fields, `[[`, "field"), use.names = FALSE)
field_range <- range(field_values, finite = TRUE)
if (diff(field_range) == 0) field_range[2] <- field_range[1] + 1
breaks <- seq(field_range[1], field_range[2], length.out = 65L)
palette <- grDevices::colorRampPalette(
  c("#EBEBEB", "#40E0D0", "#00008B")
)(64)
half_cell <- res / 2000

# Colours per grid node (snowfall palette blended with the hillshade)
make_colors <- function(field) {
  outside <- !is.finite(field)
  class_index <- findInterval(field, breaks, all.inside = TRUE)
  class_index[outside] <- 1L
  base_rgb <- col2rgb(palette[class_index]) / 255
  blend <- function(channel) {
    matrix(channel * (1 - grey_weight) + SHADE_GREY * grey_weight, nx, ny)
  }
  list(
    colors = matrix(
      rgb(blend(base_rgb[1, ]), blend(base_rgb[2, ]), blend(base_rgb[3, ])),
      nx, ny
    ),
    outside = outside
  )
}

xi <- seq(1L, nx, by = GRID_STEP)
yi <- seq(1L, ny, by = GRID_STEP)
x_km <- easting[xi] / 1000
y_km <- northing[yi] / 1000
height_km <- height / 1000
z_max <- max(height_km) * VERTICAL_EXAG

plot_fields <- function() {
  grDevices::png(out_file, width = 1800, height = 1300, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(
    matrix(c(1, 2, 0, 3, 4, 7, 5, 6, 0), nrow = 3L, byrow = TRUE),
    widths = c(1, 1, 0.45)
  )
  par(oma = c(0, 0, 3, 0))

  for (i in seq_along(fields)) {
    item <- fields[[i]]
    rendered <- make_colors(item$field)
    # Cells outside Switzerland are dropped from the surface
    z <- height_km[xi, yi] * VERTICAL_EXAG
    z[rendered$outside[xi, yi]] <- NA
    facet_colors <- rendered$colors[xi, yi][-length(xi), -length(yi)]
    par(mar = c(0, 0, 3, if (i == 6L) 0 else 1))
    panel_title <- if (i == 6L) {
      sprintf("Highest: %s | winter %d | %.2f mm/day",
              item$model, item$year, selected$snowfall[i])
    } else {
      sprintf(
        "%d. %s | winter %d | %.2f mm/day",
        i, item$model, item$year, selected$snowfall[i]
      )
    }
    persp(
      x_km, y_km, z, zlim = c(0, z_max), scale = FALSE, expand = 1,
      theta = VIEW_THETA, phi = VIEW_PHI, r = 3, col = facet_colors,
      border = NA, box = FALSE, axes = FALSE, shade = NA, main = panel_title
    )
    if (i == 6L) box(which = "plot", col = "grey75", lwd = 2)
  }

  par(mar = c(1, 0, 3, 4))
  plot(
    NA, xlim = c(0, 6), ylim = range(breaks), axes = FALSE, xaxs = "i",
    yaxs = "i", xlab = "", ylab = ""
  )
  rect(0, breaks[-65], 1, breaks[-1], col = palette, border = NA)
  axis(4, pos = 1)
  rect(0, min(breaks), 1, max(breaks), border = "black")
  text(
    x = 3.4, y = mean(range(breaks)), labels = "Mean winter snowfall [mm/day]",
    srt = 90, xpd = NA
  )
  mtext(
    "Five lowest and single highest snowfall winters (unmasked distribution), GWL3.0",
    3, 1, outer = TRUE, cex = 1.4
  )
}

plot_fields()
message("Wrote ", out_file)
