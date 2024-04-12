
# Internal
combine_dtm = function(dtm1, dtm2)
{
  terra::mosaic(dtm1, dtm2, fun = "min")
}

# Internal
compute_cluster_hulls = function(x)
{
  X <- Y <- Z <- .N <- cluster <- NULL
  hull = function(x,y,z)
  {
    z = stats::quantile(z, probs = 0.05)
    xy = cbind(x,y)
    h = grDevices::chull(xy)
    h = xy[h,]
    h = rbind(h, h[1,])
    h = cbind(h,z)
    list(h)
  }

  y = x@data[, list(hull = hull(X,Y,Z), n = .N), by = cluster]

  polys = lapply(y$hull, function(x) sf::st_sfc(sf::st_polygon(list(x))))
  polys = do.call(c, polys)
  sf::st_crs(polys) = sf::st_crs(x)
  polys
}

# Internal: minimum bounding rectangle
MBR <- function(poly)
{
  z = sf::st_coordinates(poly)[,3][1]
  p = sf::st_coordinates(poly)[,1:2]
  p = p[-1,]

  # Analyze the convex hull edges
  a <- grDevices::chull(p)                         # Indexes of extremal points
  a <- c(a, a[1])                                 # Close the loop
  e <- p[a[-1],] - p[a[-length(a)], ]             # Edge directions
  norms <- sqrt(rowSums(e^2))                     # Edge lengths
  v <- e / norms                                  # Unit edge directions
  w <- cbind(-v[,2], v[,1])                       # Normal directions to the edges

  # Find the MBR
  vertices <- p[a, ]                              # Convex hull vertices
  x <- apply(vertices %*% t(v), 2, range)         # Extremes along edges
  y <- apply(vertices %*% t(w), 2, range)         # Extremes normal to edges
  areas <- (y[1,]-y[2,])*(x[1,]-x[2,])            # Areas
  k <- which.min(areas)                           # Index of the best edge (smallest area)

  # Form a rectangle from the extremes of the best edge
  bb = cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,])

  bb = cbind(bb, z)
  bb = sf::st_polygon(list(bb))
  bb = sf::st_sfc(bb)
  sf::st_crs(bb) = sf::st_crs(poly)
  return(bb)
}
