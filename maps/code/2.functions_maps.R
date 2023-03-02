# Functions

library(sf)
library(dplyr)


#' Function to get the heading angle between two points:
#'
#' This take in any two point features of class 'sf' and returns
#'  heading angle between the two points:
#'
#' @param x point feature of class 'sf'
#' @param y point feature of class 'sf'
#' @keywords Angle, heading, Simple Features, sf
#' @export
#' @examples
#' \dontrun{
#'
#' points <- sf::read_sf("data/Points.gpkg")
#' x <- sf::st_coordinates(points[1, 4])
#' y <- sf::st_coordinates(points[2, 4])
#' }
#' calc_angle(x, y)
calc_angle <- function(x, y) {
  dst_diff <- as.numeric(x - y)
  return(atan2(dst_diff[1], dst_diff[2]) + pi)
}


#' Function to generate parallel straight lines:
#'
#' This takes in a line and draws parallel lines to it
#' across the field, separated by the distance given in width
#' @param line of class 'sfc' to draw lines parallel to
#' @param width of class 'numeric' to separate parallel lines by
#' @param field of class 'sf'
#' @param max_dist of class 'numeric'
#' @keywords Over, Intersect, Simple Features, sf
#' @export
#' @examples
#' \dontrun{
#'
#' }
parallel_lines <- function(line, width, field, offset_path = 0, p, clip = FALSE) {
  
  
  cl <- sf::st_coordinates(line)[, 1:2]
  angle <- calc_angle(cl[1, ], cl[nrow(cl), ]) + pi / 2

  nl <- map(seq(0, # might be 1 or 0 need to check with field data
                p-1), simplify = FALSE, ~sf::st_geometry(line) + (offset_path + .) * width * c(sin(angle), cos(angle)))
  
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

#' Function to get distances between consecutive points:
#'
#' This goes through every point in sf_obj and sums the distances between each consecutive point. If parameter max is FALSE sf_obj should
#' be a point feature. If max = TRUE, sf_obj should be a polygon and it will return the diagonal distance of the field.
#' @param sf_obj object of class 'sf'
#' @param max conditional variable - default = FALSE
#' @keywords Distance, Simple Features, sf
#' @export
#' @examples
#' #Default
#' points =  sf::read_sf('./data/Points.gpkg')
#' DIFMR::calc_distance(points)
#' #max = TRUE
#' field = DIFMR::fields[4,]
#' DIFMR::calc_distance(field, max = TRUE)
#'

calc_distance <- function(sf_obj, max = FALSE) {
  if (max == TRUE) {
    pts <- sf::st_cast(sf::st_as_sfc(sf::st_bbox(sf_obj)), "POINT")
    return(as.numeric(sf::st_distance(pts[1], pts[3])))
  } else {
    ab_dist <- function(a, b) {
      return(sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2))
    }
    coords <- sf::st_coordinates(sf_obj)
    result <- as.numeric(sapply(c(2:nrow(coords)), function(x) {
      ab_dist(coords[x - 1, ], coords[x, ])
    }))
    return(c(result[1], result))
  }
}
