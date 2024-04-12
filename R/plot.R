#' Plot hull 3d
#'
#' Function to display the bridges on the point cloud
#'
#' @param x The output of the function plot used with a LAS object.
#' @param polys sf of sfc polygons returned by \link{find_bridges}
#' @param col color
#' @export
add_hulls3d = function(x, polys, col = "red")
{
  if (length(col) == 1) col = rep(col, length(polys))

  if (methods::is(polys, "sf"))
  {
    z = polys$Z
    polys = sf::st_geometry(polys)

    for (i in seq_along(polys))
    {
      xyz = sf::st_coordinates(polys[i])
      rgl::lines3d(xyz[,1]-x[1], xyz[,2]-x[2], z = z[i], col = col[i], lwd = 4)
    }
  }
  else
  {
    for (i in seq_along(polys))
    {
      xyz = sf::st_coordinates(polys[i])
      rgl::lines3d(xyz[,1]-x[1], xyz[,2]-x[2], z = xyz[,3], col = col[i], lwd = 4)
    }
  }

  return(x)
}

add_stream3d = function(x, lines, z)
{
  lines = sf::st_geometry(lines)
  for (i in seq_along(lines))
  {
    xyz = sf::st_coordinates(lines[i])
    rgl::lines3d(xyz[,1]-x[1], xyz[,2]-x[2], z = z, col = "red", lwd = 4)
  }

  return(x)
}
