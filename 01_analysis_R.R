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
      c("standard", "gregorian", "proleptic_gregorian",
        "365_day", "noleap", "360_day")) {
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
  if (calendar %in% c("standard", "gregorian", "proleptic_gregorian")) {
    return(as.Date(origin + time_values * seconds_per_unit, tz = "UTC"))
  }

  # Fixed-length-year calendars: decode arithmetically. Only year and month
  # are used downstream; 360_day days are capped at 28 to stay valid dates.
  origin_parts <- as.integer(strsplit(format(origin, "%Y-%m-%d"), "-")[[1]])
  origin_fraction <- (as.numeric(origin) %% 86400) / 86400
  elapsed_days <- time_values * seconds_per_unit / 86400
  if (calendar == "360_day") {
    year_length <- 360
    origin_day <- origin_parts[1] * year_length +
      (origin_parts[2] - 1L) * 30 + origin_parts[3] - 1
  } else {
    year_length <- 365
    cumulative <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
    origin_day <- origin_parts[1] * year_length +
      cumulative[origin_parts[2]] + origin_parts[3] - 1
  }
  total_days <- floor(origin_day + origin_fraction + elapsed_days + 1e-9)
  year <- total_days %/% year_length
  day_of_year <- total_days %% year_length
  if (calendar == "360_day") {
    month <- day_of_year %/% 30 + 1
    day <- pmin(day_of_year %% 30 + 1, 28)
  } else {
    month <- findInterval(day_of_year, cumulative)
    day <- day_of_year - cumulative[month] + 1
  }
  as.Date(sprintf("%04d-%02d-%02d", as.integer(year), as.integer(month),
                  as.integer(day)))
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

get_result_paths <- function(model) {
  result_base <- file.path(RESDIR, sprintf("result_%s_%s", model, PERIOD))
  suffix_number <- 0L
  repeat {
    suffix <- if (suffix_number == 0L) {
      ""
    } else {
      sprintf("_%d", suffix_number)
    }
    paths <- list(
      nc = sprintf("%s%s_R.nc", result_base, suffix),
      txt = sprintf("%s%s_R.txt", result_base, suffix),
      field = file.path(RESDIR, sprintf(
        "result_field_%s_%s%s_R.nc", model, PERIOD, suffix
      )),
      field_nonmasked = file.path(RESDIR, sprintf(
        "result_field_%s_%s%s_nonmasked_R.nc", model, PERIOD, suffix
      ))
    )
    if (!any(file.exists(unlist(paths)))) {
      return(paths)
    }
    suffix_number <- suffix_number + 1L
  }
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

  # month_groups/month_dates are only used to classify each day's calendar
  # month (for the Nov-Apr selection, winter-year assignment, and trimming
  # of the incomplete leading/trailing winter below). No monthly mean or
  # sum is ever computed; daily values are accumulated directly into
  # winter-year bins, consistent with cdo's yearmean applied directly to
  # daily values (01_analysis_cdo_v02.sh).
  winter_months <- as.integer(format(month_dates, "%m")) %in%
    c(11L, 12L, 1L, 2L, 3L, 4L)
  selected <- which(winter_months)
  if (length(selected) <= 6L) {
    stop(sprintf("Not enough Nov-Apr data to trim for %s.", model),
         call. = FALSE)
  }
  # Drop the first four and last two Nov-Apr months: they are always an
  # incomplete leading/trailing winter season.
  selected <- selected[seq.int(5L, length(selected) - 2L)]
  selected_dates <- month_dates[selected]
  month_numbers <- as.integer(format(selected_dates, "%m"))
  winter_years <- as.integer(format(selected_dates, "%Y")) +
    as.integer(month_numbers >= 11L)
  output_years <- sort(unique(winter_years))
  year_row <- match(winter_years, output_years)

  # Day-weighted annual mean: accumulate daily sums and valid-day counts
  # directly into the winter-year they belong to, then divide once per
  # year (equivalent to cdo's yearmean on daily data).
  annual_sums <- matrix(0, nrow = length(output_years), ncol = number_of_cells)
  annual_counts <- matrix(0L, nrow = length(output_years), ncol = number_of_cells)
  for (i in seq_along(selected)) {
    steps <- month_groups[[selected[i]]]
    tas <- read_time_block(
      tas_nc, tas_name, tas_var, tas_time_index, min(steps), length(steps)
    )
    pr <- read_time_block(
      pr_nc, pr_name, pr_var, pr_time_index, min(steps), length(steps)
    )

    snow <- pr
    snow[!is.na(tas) & tas > TEMPERATURE_THRESHOLD] <- 0
    snow[is.na(tas) | is.na(pr)] <- NA_real_

    row <- year_row[i]
    annual_counts[row, ] <- annual_counts[row, ] + colSums(!is.na(snow))
    annual_sums[row, ] <- annual_sums[row, ] + colSums(snow, na.rm = TRUE)
  }

  annual_values <- matrix(
    NA_real_, nrow = length(output_years), ncol = number_of_cells
  )
  valid <- annual_counts > 0L
  annual_values[valid] <- annual_sums[valid] / annual_counts[valid]

  topo <- read_topography(topo_path, spatial_names, spatial_lengths)
  # The supplied *_gec1000.nc file contains 1 for included cells and missing
  # values elsewhere; it is not the original elevation field in metres.
  annual_values_nonmasked <- annual_values
  annual_values[, !(is.finite(topo) & topo > 0)] <- NA_real_
  area_mean <- rowMeans(annual_values, na.rm = TRUE)
  area_mean[!is.finite(area_mean)] <- NA_real_

  result_paths <- get_result_paths(model)
  write.table(
    data.frame(year = output_years, value = area_mean),
    file = result_paths$txt,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  dimension_values <- lapply(spatial_dims, function(dimension) dimension$vals)
  names(dimension_values) <- spatial_names
  output_dimensions <- c(
    list(ncdim_def("year", "year", as.double(output_years))),
    lapply(seq_along(spatial_dims), function(index) {
      ncdim_def(
        spatial_names[[index]],
        spatial_dims[[index]]$units %||% "1",
        dimension_values[[index]]
      )
    })
  )

  missing_value <- -9999
  units <- ncatt_get(pr_nc, pr_name, "units")$value
  write_nc <- function(path, values, dimensions, dim_lengths, longname) {
    variable <- ncvar_def(
      "snowfall_proxy", units, dimensions,
      missval = missing_value, longname = longname
    )
    nc <- nc_create(path, variable)
    on.exit(nc_close(nc), add = TRUE)
    output_values <- array(
      as.vector(t(values)),
      dim = c(length(output_years), dim_lengths)
    )
    output_values[!is.finite(output_values)] <- missing_value
    ncvar_put(nc, variable, output_values)
  }

  # Area mean as year x 1 x 1, like the CDO fldmean result
  mean_dimensions <- c(
    list(ncdim_def("year", "year", as.double(output_years))),
    lapply(seq_along(spatial_dims), function(index) {
      ncdim_def(
        spatial_names[[index]],
        spatial_dims[[index]]$units %||% "1",
        mean(dimension_values[[index]])
      )
    })
  )
  write_nc(
    result_paths$nc, matrix(area_mean, ncol = 1L), mean_dimensions,
    rep(1L, length(spatial_dims)),
    "Annual mean snowfall proxy above 1000 m, spatial mean"
  )
  write_nc(
    result_paths$field, annual_values, output_dimensions, spatial_lengths,
    "Annual mean snowfall proxy above 1000 m"
  )
  write_nc(
    result_paths$field_nonmasked, annual_values_nonmasked,
    output_dimensions, spatial_lengths,
    "Annual mean snowfall proxy, all elevations"
  )

  message(sprintf(
    "Finished %s: wrote %s",
    model, paste(unlist(result_paths), collapse = ", ")
  ))
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
