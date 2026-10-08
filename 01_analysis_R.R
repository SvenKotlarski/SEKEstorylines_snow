#!/usr/bin/env Rscript

# Compute the Nov-Apr mean snowfall proxy above 1000 m, then average over
# the Swiss grid. Requires the ncdf4 package.

suppressPackageStartupMessages(library(ncdf4))

DATADIR <- "/net/stratus/c2sm-data/CH2025/ogd-climate-scenarios-ch2025-grid"
METADIR <- "/highres/svenk/SEKEstorylines_snow/META"
RESDIR <- "/highres/svenk/SEKEstorylines_snow/RESULTS"
FILE_BASE <- "ogd-climate-scenarios-ch2025-grid_ch"
PERIOD <- "gwl3.0"
TEMPERATURE_THRESHOLD <- 2

get_data_variable <- function(nc, preferred_name) {
  if (preferred_name %in% names(nc$var)) {
    return(preferred_name)
  }
  stop(
    sprintf("Variable '%s' is missing from %s", preferred_name, nc$filename),
    call. = FALSE
  )
}

get_time_axis <- function(variable) {
  matches <- which(vapply(
    variable$dim,
    function(dimension) tolower(dimension$name) == "time",
    logical(1)
  ))
  if (length(matches) != 1L) {
    stop("Expected exactly one time dimension.", call. = FALSE)
  }
  matches
}

decode_time <- function(time_values, units, calendar) {
  if (is.null(calendar) || !calendar %in%
      c("standard", "gregorian", "proleptic_gregorian")) {
    stop(
      sprintf("Unsupported CF calendar: %s", calendar %||% "<missing>"),
      call. = FALSE
    )
  }

  unit_match <- regexec(
    "^(seconds|minutes|hours|days) since (.+)$",
    units,
    ignore.case = TRUE
  )
  unit_parts <- regmatches(units, unit_match)[[1]]
  if (length(unit_parts) != 3L) {
    stop(sprintf("Unsupported CF time units: %s", units), call. = FALSE)
  }

  origin_text <- sub(" (UTC|Z)$", "", unit_parts[3], ignore.case = TRUE)
  origin <- as.POSIXct(
    origin_text,
    tz = "UTC",
    tryFormats = c(
      "%Y-%m-%d %H:%M:%OS",
      "%Y-%m-%d %H:%M",
      "%Y-%m-%d"
    )
  )
  if (is.na(origin)) {
    stop(sprintf("Could not parse time origin: %s", origin_text), call. = FALSE)
  }

  seconds_per_unit <- switch(
    tolower(unit_parts[2]),
    seconds = 1,
    minutes = 60,
    hours = 3600,
    days = 86400
  )
  as.Date(origin + time_values * seconds_per_unit, tz = "UTC")
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

read_time_block <- function(nc, variable_name, variable, time_index,
                            first_step, number_of_steps) {
  start <- rep(1L, length(variable$dim))
  count <- vapply(variable$dim, function(dimension) dimension$len, integer(1))
  start[time_index] <- first_step
  count[time_index] <- number_of_steps

  values <- ncvar_get(
    nc,
    variable_name,
    start = start,
    count = count,
    collapse_degen = FALSE
  )
  permutation <- c(time_index, setdiff(seq_along(start), time_index))
  values <- aperm(values, permutation)
  matrix(values, nrow = number_of_steps)
}

read_topography <- function(path, spatial_names, spatial_lengths) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))

  variable_name <- if ("height" %in% names(nc$var)) {
    "height"
  } else {
    stop(sprintf("Topography variable 'height' is missing from %s", path),
         call. = FALSE)
  }
  variable <- nc$var[[variable_name]]
  time_index <- get_time_axis(variable)
  topo_spatial_names <- vapply(
    variable$dim[-time_index],
    function(dimension) dimension$name,
    character(1)
  )
  topo_spatial_lengths <- vapply(
    variable$dim[-time_index],
    function(dimension) dimension$len,
    integer(1)
  )

  if (!setequal(topo_spatial_names, spatial_names) ||
      !identical(
        unname(topo_spatial_lengths[match(spatial_names, topo_spatial_names)]),
        unname(spatial_lengths)
      )) {
    stop("Topography and climate data grids do not match.", call. = FALSE)
  }

  values <- ncvar_get(
    nc,
    variable_name,
    start = rep(1L, length(variable$dim)),
    count = vapply(variable$dim, function(dimension) dimension$len, integer(1)),
    collapse_degen = FALSE
  )
  spatial_indices <- setdiff(seq_along(variable$dim), time_index)
  permutation <- c(time_index, spatial_indices[
    match(spatial_names, topo_spatial_names)
  ])
  values <- aperm(values, permutation)
  as.vector(matrix(values, nrow = 1L))
}

process_model <- function(model, topo_path) {
  tas_path <- file.path(
    DATADIR, "tas", sprintf("%s_tas_%s_%s.nc", FILE_BASE, model, PERIOD)
  )
  pr_path <- file.path(
    DATADIR, "pr", sprintf("%s_pr_%s_%s.nc", FILE_BASE, model, PERIOD)
  )
  if (!file.exists(tas_path) || !file.exists(pr_path)) {
    stop(
      sprintf("Missing input for %s; expected:\n%s\n%s", model, tas_path, pr_path),
      call. = FALSE
    )
  }

  tas_nc <- nc_open(tas_path)
  on.exit(nc_close(tas_nc), add = TRUE)
  pr_nc <- nc_open(pr_path)
  on.exit(nc_close(pr_nc), add = TRUE)

  tas_name <- get_data_variable(tas_nc, "tas")
  pr_name <- get_data_variable(pr_nc, "pr")
  tas_var <- tas_nc$var[[tas_name]]
  pr_var <- pr_nc$var[[pr_name]]
  tas_time_index <- get_time_axis(tas_var)
  pr_time_index <- get_time_axis(pr_var)

  tas_dim_names <- vapply(tas_var$dim, function(dimension) dimension$name, "")
  pr_dim_names <- vapply(pr_var$dim, function(dimension) dimension$name, "")
  if (!identical(tas_dim_names, pr_dim_names) ||
      !identical(
        vapply(tas_var$dim, function(dimension) dimension$len, integer(1)),
        vapply(pr_var$dim, function(dimension) dimension$len, integer(1))
      ) ||
      tas_time_index != pr_time_index) {
    stop(sprintf("Temperature and precipitation grids differ for %s.", model),
         call. = FALSE)
  }
  if (length(tas_var$dim) != 3L) {
    stop("Expected time plus two spatial dimensions.", call. = FALSE)
  }

  time_values <- tas_var$dim[[tas_time_index]]$vals
  pr_time_values <- pr_var$dim[[pr_time_index]]$vals
  if (length(time_values) != length(pr_time_values) ||
      !isTRUE(all.equal(time_values, pr_time_values))) {
    stop(sprintf("Temperature and precipitation times differ for %s.", model),
         call. = FALSE)
  }
  time_name <- tas_var$dim[[tas_time_index]]$name
  calendar_attribute <- ncatt_get(tas_nc, time_name, "calendar")
  time_calendar <- if (calendar_attribute$hasatt) {
    calendar_attribute$value
  } else {
    "standard"
  }
  dates <- decode_time(
    time_values,
    tas_var$dim[[tas_time_index]]$units,
    time_calendar
  )
  if (anyNA(dates)) {
    stop(sprintf("Could not decode all dates for %s.", model), call. = FALSE)
  }

  spatial_dims <- tas_var$dim[-tas_time_index]
  spatial_names <- vapply(spatial_dims, function(dimension) dimension$name, "")
  spatial_lengths <- vapply(spatial_dims, function(dimension) dimension$len, integer(1))
  number_of_cells <- prod(spatial_lengths)
  month_key <- format(dates, "%Y-%m")
  month_groups <- split(seq_along(dates), factor(
    month_key,
    levels = unique(month_key)
  ))
  month_dates <- as.Date(paste0(names(month_groups), "-15"))

  monthly_values <- matrix(
    NA_real_,
    nrow = length(month_groups),
    ncol = number_of_cells
  )
  for (month_index in seq_along(month_groups)) {
    steps <- month_groups[[month_index]]
    tas <- read_time_block(
      tas_nc, tas_name, tas_var, tas_time_index, min(steps), length(steps)
    )
    pr <- read_time_block(
      pr_nc, pr_name, pr_var, pr_time_index, min(steps), length(steps)
    )

    snow <- pr
    snow[is.na(tas) | tas > TEMPERATURE_THRESHOLD] <- NA_real_
    snow[is.na(pr)] <- NA_real_
    valid_counts <- colSums(!is.na(snow))
    monthly_sum <- colSums(snow, na.rm = TRUE)
    valid <- valid_counts > 0L
    monthly_values[month_index, valid] <-
      monthly_sum[valid] / valid_counts[valid]
  }

  winter_months <- as.integer(format(month_dates, "%m")) %in%
    c(11L, 12L, 1L, 2L, 3L, 4L)
  selected <- which(winter_months)
  if (length(selected) <= 6L) {
    stop(sprintf("Not enough Nov-Apr data to trim for %s.", model),
         call. = FALSE)
  }
  selected <- selected[seq.int(5L, length(selected) - 2L)]
  selected_dates <- month_dates[selected]
  selected_values <- monthly_values[selected, , drop = FALSE]

  month_numbers <- as.integer(format(selected_dates, "%m"))
  winter_years <- as.integer(format(selected_dates, "%Y")) +
    as.integer(month_numbers >= 11L)
  output_years <- sort(unique(winter_years))
  annual_values <- t(vapply(output_years, function(year) {
    values <- colMeans(
      selected_values[winter_years == year, , drop = FALSE],
      na.rm = TRUE
    )
    values[!is.finite(values)] <- NA_real_
    values
  }, numeric(number_of_cells)))

  topo <- read_topography(topo_path, spatial_names, spatial_lengths)
  # The supplied *_gec1000.nc file contains 1 for included cells and missing
  # values elsewhere; it is not the original elevation field in metres.
  annual_values[, !(is.finite(topo) & topo > 0)] <- NA_real_
  area_mean <- rowMeans(annual_values, na.rm = TRUE)
  area_mean[!is.finite(area_mean)] <- NA_real_

  result_txt <- file.path(RESDIR, sprintf("result_%s_%s.txt", model, PERIOD))
  result_field <- file.path(
    RESDIR, sprintf("result_field_%s_%s.rds", model, PERIOD)
  )
  write.table(
    data.frame(year = output_years, value = area_mean),
    file = result_txt,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  dimension_values <- lapply(spatial_dims, function(dimension) dimension$vals)
  names(dimension_values) <- spatial_names
  dimension_labels <- lapply(dimension_values, as.character)
  field <- lapply(seq_along(output_years), function(index) {
    array(annual_values[index, ], dim = spatial_lengths,
          dimnames = setNames(dimension_labels, spatial_names))
  })
  names(field) <- as.character(output_years)
  saveRDS(
    list(
      years = output_years,
      field = field,
      spatial_dimensions = dimension_values,
      spatial_dimension_names = spatial_names,
      units = ncatt_get(pr_nc, pr_name, "units")$value
    ),
    result_field
  )

  message(sprintf("Finished %s: wrote %s and %s", model, result_txt, result_field))
}

model_file <- file.path(METADIR, "model_list.txt")
if (!file.exists(model_file)) {
  stop(sprintf("Model list not found: %s", model_file), call. = FALSE)
}
if (!dir.exists(RESDIR)) {
  stop(sprintf("Results directory not found: %s", RESDIR), call. = FALSE)
}

topo_path <- file.path(METADIR, "topo.swiss1_ch01r.swiss.lv95_gec1000.nc")
if (!file.exists(topo_path)) {
  stop(sprintf("Elevation mask not found: %s", topo_path), call. = FALSE)
}
models <- trimws(readLines(model_file, warn = FALSE))
models <- models[nzchar(models)]
if (length(models) == 0L) {
  stop(sprintf("Model list is empty: %s", model_file), call. = FALSE)
}
for (model in models) {
  process_model(model, topo_path)
}
