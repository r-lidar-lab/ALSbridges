jr_internal_stream_query = function(bbox)
{
  db = "/home/jr/Documents/Ulaval/2024 PDR/Projet ponts/Water stream/streams.gpkg"
  sql = "SELECT * FROM stream WHERE ST_Intersects(geom, ST_GeomFromText('{wkt}'));"

  from = sf::st_crs(bbox)
  if (is.na(from)) stop("The bounding box does not have any CRS")

  bbox = sf::st_as_sfc(bbox)
  target = sf::st_crs(32198)

  bb = bbox
  if (from != target) bb = sf::st_transform(bb, target)

  wkt = sf::st_as_text(bb)

  query = glue::glue(sql)

  u = sf::st_read(db, query = query, quiet = TRUE)
  if (from != target) u = sf::st_transform(u, from)

  u = sf::st_crop(u, bbox)
  return(u)
}
