# Functions

library(sf)
library(dplyr)

calc_angle <- function(x, y) {
  dst_diff <- as.numeric(x - y)
  return(atan2(dst_diff[1], dst_diff[2]) + pi)
}

parallel_lines <- function(line, width, field, offset_path = 0, passes = 50, clip = FALSE) {
  
  cl <- sf::st_coordinates(line)[, 1:2]
  angle <- calc_angle(cl[1, ], cl[nrow(cl), ]) + pi / 2
  #n <- ceiling(max_dist / width)
  nl <- map(seq(0, n), simplify = FALSE, ~sf::st_geometry(line) + (offset_path + .) * width * c(sin(angle), cos(angle)))
  
  path_lines <- sf::st_sfc(do.call(rbind, nl), crs = sf::st_crs(line))
  
  path_lines <- st_transform(path_lines, 4326)
  
  # Remove the lines outside the field:
  fieldb <- sf::st_union(sf::st_buffer(field, width / 2))
  crit <- !is.na(as.numeric(sf::st_intersects(path_lines, fieldb)))
  path_lines <- path_lines[crit]
  df_id <- data.frame(id = 1:length(path_lines))
  path_lines <- sf::st_as_sf(df_id, geometry = path_lines)
  if (clip) {
    fieldb <- sf::st_buffer(field, width)
    sf::st_agr(path_lines) <- "constant"
    path_lines <- sf::st_intersection(path_lines, sf::st_geometry(fieldb))
    sf::st_agr(path_lines) <- "constant"
    path_lines <- sf::st_cast(sf::st_cast(path_lines, "MULTILINESTRING"), "LINESTRING")
  }
  return(path_lines)
}