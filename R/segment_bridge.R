#' Find the bridges in a point cloud
#'
#' Find bridges in a point cloud and return their bounding box. The function does not scan the entire
#' point cloud but rather relies on an existing map that is supposed to be approximately correct. Then,
#' the function searches at this location within a buffer.
#'
#' About planarity: we use a different definition than the regular definition. \eqn{a_1, a_2, a_3}{a1, a2, a3}
#' denote the eigenvalues of the covariance matrix of the neighboring points in ascending order. The planarity
#' is \eqn{(a_2 - a_1) / a_3}{(a2 - a1) / a3} and ranges in  \eqn{[0, 1]}. This definition is the correct definition
#' of the planarity but, in practice, may lead to weird behaviors in some point cloud sampling patterns. Instead,
#' we use \eqn{a_3 / a_1}{a3/a1}. This value ranges between \eqn{[0, +infinity]}. This is not the correct definition
#' of the planarity but appears to work effectively better without being affected by some sampling patterns.
#' `th_planarity` = 15 is very tolerant, 40 is tolerant; 75 will detect only very coplanar points. A
#' planar object may easily have a planarity threshold above 100, but with noisy data, 40 or even 15 is more likely.
#' If `th_planarity` < 1 then the algorithm uses the true definition of the planarity.
#'
#' @param las A LAS object from lidR
#' @param location numeric. The approximate coordinates of the bridges to look for (in the CRS of the point cloud).
#' Default is NULL: the method searches from the centroid of the point cloud.
#' @param buffer numeric. Bigger than the bridge. This is used to clip the point cloud and reduce the search area,
#' but it is not used to find the bridge itself. If the buffer is too small and the bridge is partially
#' missing, the bridge could be missed. If the buffer is too large, the computation time is longer.
#' @param th_planarity numeric. A threshold to estimate the planarity of the local geometry. More than
#' this value is considered planar (see details).
#' @param th_angle numeric. A threshold to estimate the horizontality of the local geometry. A bridge is
#' planar AND horizontal. An angle of 5 degrees is enough in most cases, but with noisy data being more tolerant
#' may be suitable.
#' @param th_hag numeric. A threshold on height above ground. A bridge is a horizontal plane above the ground.
#' Above the ground is defined by this parameter. Default 90 cm.
#' @param query_stream a function capable of querying from a database of streams. It is used for robust
#' filtering of false positive bridges. It can be omitted. To know how to create your own function, look
#' at the code of `ALSbridges:::jr_internal_stream_query`. Copy-paste the code and adjust to your needs.
#' @param force_stream boolean. By default, the streams are queried only if there are more than one bridge
#' found. If force_stream = TRUE, a stream is always queried to check the bridges.
#' @param ... Internals for debugging, e.g., display = TRUE.
#' @export
#' @returns NULL if no bridge found. An `sf` object with the bounding boxes of the bridges found.
#' @md
#' @seealso [classify_bridges()]
#' @examples
#' library(lidR)
#' library(ALSbridges)
#'
#' LASfile <- system.file("extdata", "16_2175571f05_dc_000045.laz", package="ALSbridges")
#' las = readLAS(LASfile)
#'
#' bridge = find_bridges(las, th_planarity = 50, th_angle = 5, buffer = 130, th_hag = 0.9)
#' bridge
#'
#' plot(header(las))
#' plot(bridge, add = TRUE)
#'
#' plot(las) |> add_hulls3d(bridge, col = "red")
#'
#' las = classify_bridges(las, bridge, upper_buffer = 25)
#' plot(las, color = "Classification")
find_bridges = function(las, location = NULL,  buffer = 75, th_planarity = 40, th_angle = 8, th_hag = 0.9, query_stream = NULL, force_stream = FALSE,...)
{
  Classification <- Z <- Shape <- cluster <- NULL

  debug_options = list(...)
  display = isTRUE(debug_options$display)
  display_dtm = isTRUE(debug_options$display_dtm)
  interactive = isTRUE(debug_options$interactive)

  # Pulse density
  area <- as.numeric(sf::st_area(las))
  npulses <- las@header[["Number of points by return"]][1]
  d <- npulses/area

  # Extract a smaller region of interest
  if (is.null(location))
  {
    centerx = mean(las$X)
    centery = mean(las$Y)
  }
  else
  {
    centerx = location[1]
    centery = location[2]
  }
  las = lidR::clip_circle(las, centerx, centery, buffer+5)

  # Compute the original DTM
  #cat("Computing DTM from original ground classification...\n")
  dtm1 = lidR::rasterize_terrain(las, 0.5)

  # Recompute the ground points from original ground points + extra water
  #cat("Reclassifying ground points...\n")
  ground = lidR::filter_poi(las, Classification %in% c(lidR::LASGROUND, lidR::LASWATER))
  #if (!no_water)
  #{
  #  ground@data = ground@data[, .(X,Y,Z, Classification)]
  #  ground@data = rbind(ground@data, extra_water@data)
  #}
  reground = lidR::classify_ground(ground, lidR::mcc(2, 0.1), last_returns = F) |> suppressMessages()

  #plot(reground, color = "Classification")

  # Compute the reground DTM
  #cat("Computing DTM from new ground classification...\n")
  dtm2 = lidR::rasterize_terrain(reground, 0.5)

  # Combine the two dtm such as the miss-classified bridge are no longer ground
  dtm = combine_dtm(dtm1, dtm2)

  if (display_dtm)
  {
    lidR::plot_dtm3d(dtm1)
    lidR::plot_dtm3d(dtm2)
    lidR::plot_dtm3d(dtm)
  }

  # Compute HAG
  las = lidR::merge_spatial(las, dtm, "hag")
  las$hag = las$Z - las$hag
  las$hag[las$hag < 0] = 0

  # Segment horizontal planes
  #cat("Segmenting horizontal planes...\n")
  #las = segment_shapes(las, shp_hplane(th1 = 10, th2 = 6, th3 = 0.99, k = 40))
  u = lidR::point_eigenvalues(las, k = 30, r = 2, coeffs = T)
  u$eigen_smallest[u$eigen_smallest == 0] = 1e-8
  angle = (acos(abs(u$coeff22))%%pi)*180/pi

  if (th_planarity < 1) # true planarity
  {
    planarity = (u$eigen_medium-u$eigen_smallest)/u$eigen_largest
  }
  else # works better
  {
    planarity = (u$eigen_largest/u$eigen_smallest)
  }

  planarity[is.infinite(planarity)] = 0
  planarity[is.nan(planarity)] = 0

  rm(u)

  if (display)
  {
    las@data$planarity = planarity
    las@data$angle = angle
    if (th_planarity > 1)
      lidR::plot(las, color = "planarity", breaks = "quantile", legend = T)
    else
      lidR::plot(las, color = "planarity", legend = T)
    lidR::plot(las, color = "angle", legend = T)
    lidR::plot(las, color = "hag", legend = T)
  }

  # A bridge is a plane above the ground
  las@data$Shape = planarity > th_planarity & angle < th_angle
  las@data$bridge = las$Shape == TRUE & las$hag > th_hag

  if (display)
  {
    x = lidR::plot(lidR::filter_poi(las, Shape == FALSE), color = "hag", size = 2)
    lidR::plot(lidR::filter_poi(las, Shape == TRUE), add = x, pal = "pink", size = 2)
    if (sum(las$bridge) > 0 ) lidR::plot(lidR::filter_poi(las, bridge == TRUE), add = x, pal = "red", size = 4)
    #if (!is.null(extra_water) && !is.empty(extra_water)) plot(extra_water, add = x, pal = "cyan")
    #lidR::add_dtm3d(x, dtm)
  }

  # Keep only the bridge points
  bridge = lidR::filter_poi(las, bridge == TRUE)

  if (lidR::is.empty(bridge))
  {
    warning("No bridge found...\n")
    return(NULL)
  }

  # Cluster the bridge points
  #cat("Clustering potential bridges...\n")
  coord = sf::st_coordinates(bridge)
  ps = 4*1/sqrt(d)
  Dbscan_cl <- dbscan::dbscan(coord, eps = ps, minPts = 5)

  bridge@data$cluster = Dbscan_cl$cluster
  sizes = table(Dbscan_cl$cluster)
  min_cluster_size = 10*d
  sizes = sizes[sizes > min_cluster_size]
  valid_cluster = as.numeric(names(sizes))
  valid_cluster = valid_cluster[valid_cluster > 0]

  if (length(valid_cluster) == 0)
  {
    warning("No bridge found...\n")
    return(NULL)
  }

  # Remove false positive clusters
  bridge = lidR::filter_poi(bridge, cluster %in% valid_cluster)

  # Compute the hull of each bridge
  #cat("Computing the hull of the bridge...\n")
  polys = compute_cluster_hulls(bridge)
  bbox  = lapply(polys, MBR)
  bbox = do.call(c, bbox)
  sf::st_crs(bbox) = sf::st_crs(polys)
  zpont = sapply(polys, function(x) { xyz = sf::st_coordinates(x) ; xyz[1,3]})

  # Buffer and merge
  buffered_bbox = sf::st_buffer(bbox, 1.75)
  intersections = sf::st_intersects(buffered_bbox)
  intersections = unique(intersections)

  z = lapply(intersections, function(x)
  {
    dz = diff(range(zpont[x]))
    z = if (dz < 2) mean(zpont[x]) else zpont[x]
    return(z)
  })

  buffered_bbox = lapply(intersections, function(x)
  {
    dz = diff(range(zpont[x]))
    z = if (dz < 2) sf::st_union(buffered_bbox[x]) else buffered_bbox[x]
    return(z)
  })

  buffered_bbox = do.call(c, buffered_bbox)
  buffered_bbox = lapply(buffered_bbox, MBR)
  buffered_bbox = do.call(c, buffered_bbox)
  buffered_bbox = sf::st_zm(buffered_bbox)
  sf::st_crs(buffered_bbox) = sf::st_crs(las)
  zpont = do.call(c, z)

  # Classify using the polygon
  keep = lapply(seq_along(buffered_bbox), function(i)
  {
    hull = buffered_bbox[i]
    z = zpont[i]

    # Check the connection to ground points to filter roofs that are not connected to ground
    valid = las[hull]
    pont = lidR::filter_poi(valid, Z > z-0.74, Z < z+4)
    ngnd = sum(pont$Classification == lidR::LASGROUND)

    # Check if the bridge is in the forest. It might be a FP
    above_pont = lidR::filter_poi(valid, Z > z-0.75)
    ratio = lidR::npoints(pont)/lidR::npoints(above_pont)

    if(ngnd > 0 & ratio > 0.7)
    {
      hull = sf::st_as_sf(hull)
      hull$Z = z
      return(data.frame(keep = TRUE, Z = z, ngnd = ngnd, ratio = ratio))
    }
    else
    {
      return(data.frame(keep = FALSE, Z = z, ngnd = ngnd, ratio = ratio))
    }
  })

  keep = do.call(rbind, keep)

  bpoly = sf::st_as_sf(buffered_bbox)
  bpoly$Z = keep$Z

  if (display)
  {
    add_hulls3d(x, bpoly, col = ifelse(keep$keep, "green", "red"))
  }

  if (!any(keep$keep))
  {
    warning("No bridge found...\n")
    return(NULL)
  }

  hull = bpoly[keep$keep,]

  # More than 1 bridge? Use the stream to validate the false positive
  if (nrow(hull) > 0 & (nrow(hull) > 1 | force_stream))
  {
    if (!is.null(query_stream))
    {
      if (sf::st_crs(las) != sf::NA_crs_)
      {
        cat("Using stream database to remove false positive\n")
        stream = query_stream(sf::st_bbox(las))
        if (nrow(stream) > 0)
        {
          u = sf::st_intersects(sf::st_buffer(hull, -1.5), stream)
          u = sapply(u, length) > 0
          hull = hull[u,]
        }
        else
        {
          warning("No water in this area. Returning all the objects found for manual debugging")
        }
      }
      else
      {
        warning("Point cloud without CRS: cannot query the streams")
      }
    }
  }

  return(hull)
}

#' Classify bridges using the output of \link{find_bridges}
#'
#' For an example see \link{find_bridges}.?
#'
#' @param las A LAS object.
#' @param bbox The bounding box found by \link{find_bridges}.
#' @param bottom_buffer Points below the bounding box are classified as bridge.
#' @param upper_buffer Points above the bounding box are classified as bridge. 2 m
#' in forest, 20 m in urban looks good.
#' @export
classify_bridges = function(las, bbox, bottom_buffer = 0.75, upper_buffer = 2)
{
  if (!is.null(bbox))
  {
    z = bbox$Z
    geom = sf::st_geometry(bbox)

    for (i in seq_along(geom))
    {
      las = lidR::classify_poi(las, lidR::LASBRIGDE, ~Z > z[i]-bottom_buffer & Z < z[i]+upper_buffer, geom[i])
    }
  }

  return(las)
}
