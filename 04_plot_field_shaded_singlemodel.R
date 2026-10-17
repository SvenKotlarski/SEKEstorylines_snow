#!/usr/bin/env Rscript

# Plot the first time step of an annual-mean snowfall field with a subtle grey
# topographic shading. Usage: Rscript 03_plot_field_shaded.R [field.nc] [out.png]

suppressPackageStartupMessages(library(ncdf4))

BASEDIR <- "/highres/svenk/SEKEstorylines_snow"
args <- commandArgs(trailingOnly = TRUE)
field_file <- if (length(args) >= 1L) args[1] else file.path(
  BASEDIR, "RESULTS", "result_field_clmcom-cclm4-miroc_gwl3.0_nonmasked.nc"
)
out_file <- if (length(args) >= 2L) args[2] else file.path(
  BASEDIR, "WORK", "plot_nonmasked_miroc_field1_shaded.png"
)
topo_file <- file.path(BASEDIR, "META", "topo.swiss1_ch01r.swiss.lv95.nc")

SHADE_STRENGTH <- 0.3 # max share of grey blended into the snowfall colours
SHADE_GREY <- 0.35

topo_nc <- nc_open(topo_file)
height <- ncvar_get(topo_nc, "height", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(topo_nc)

nc <- nc_open(field_file)
easting <- ncvar_get(nc, "E")
northing <- ncvar_get(nc, "N")
field <- ncvar_get(nc, "pr", start = c(1, 1, 1), count = c(-1, -1, 1))
nc_close(nc)
stopifnot(identical(dim(height), dim(field)))

res <- mean(diff(easting))
nx <- nrow(height)
ny <- ncol(height)
height[!is.finite(height)] <- 0
x_prev <- c(1, seq_len(nx - 1L))
x_next <- c(seq_len(nx)[-1L], nx)
y_prev <- c(1, seq_len(ny - 1L))
y_next <- c(seq_len(ny)[-1L], ny)
dzdx <- (height[x_next, ] - height[x_prev, ]) / (2 * res)
dzdy <- (height[, y_next] - height[, y_prev]) / (2 * res)
slope <- atan(sqrt(dzdx^2 + dzdy^2))
aspect <- atan2(dzdy, -dzdx)

# Sun from the north-west, 40 degrees above the horizon
azimuth <- (360 - 315 + 90) * pi / 180
altitude <- 40 * pi / 180
shade <- cos(pi / 2 - altitude) * cos(slope) +
  sin(pi / 2 - altitude) * sin(slope) * cos(azimuth - aspect)
shade <- pmin(pmax(shade, 0), 1)
# Flat terrain gives no shading; only slopes darker than flat ground are greyed
flat_shade <- cos(pi / 2 - altitude)
grey_weight <- SHADE_STRENGTH * pmin(pmax((flat_shade - shade) / flat_shade, 0), 1)
grey_weight <- grey_weight / max(grey_weight) * SHADE_STRENGTH

breaks <- seq(0, max(field, na.rm = TRUE), length.out = 65)
palette <- hcl.colors(64, "YlGnBu", rev = TRUE)
class_index <- pmin(pmax(findInterval(field, breaks, all.inside = TRUE), 1), 64)
base_rgb <- col2rgb(palette[class_index]) / 255
blend <- function(channel) {
  matrix(channel * (1 - grey_weight) + SHADE_GREY * grey_weight, nx, ny)
}
red <- blend(base_rgb[1, ])
green <- blend(base_rgb[2, ])
blue <- blend(base_rgb[3, ])
outside <- is.na(field)
red[outside] <- green[outside] <- blue[outside] <- 1
image_matrix <- t(matrix(rgb(red, green, blue), nx, ny))[ny:1, ]

png(out_file, 1100, 760, res = 120)
layout(matrix(1:2, 1), widths = c(5, 1))
par(mar = c(1, 1, 3, 1))
plot(
  NA,
  xlim = range(easting) / 1000, ylim = range(northing) / 1000, asp = 1,
  axes = FALSE, xlab = "", ylab = "",
  main = paste0(basename(field_file), "\nfirst time step")
)
half_cell <- res / 2000
rasterImage(
  as.raster(image_matrix),
  min(easting) / 1000 - half_cell, min(northing) / 1000 - half_cell,
  max(easting) / 1000 + half_cell, max(northing) / 1000 + half_cell,
  interpolate = FALSE
)
# Border of Switzerland: edge of the valid-data area of the grid
contour(
  easting / 1000, northing / 1000, matrix(as.numeric(!outside), nx, ny),
  levels = 0.5, drawlabels = FALSE, add = TRUE, lwd = 1.2, col = "grey20"
)
par(mar = c(1, 2.6, 3, 4))
plot(NA, xlim = 0:1, ylim = range(breaks), axes = FALSE, xaxs = "i",
     yaxs = "i", xlab = "", ylab = "")
rect(0, breaks[-65], 1, breaks[-1], col = palette, border = NA)
axis(4)
box()
mtext("Mean winter snowfall [mm/day]", 4, 2.5)
invisible(dev.off())
message("Wrote ", out_file)
