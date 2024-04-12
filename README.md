# ALSbridges

Segment the bridges in an ALS dataset. Project created for Québec Ministry and provided under GPL-3 license

## Installation

``` r
remotes::install_github("metafor-ulaval/ALSbridges")
```

## Example

``` r
library(lidR)
library(ALSbridges)

LASfile <- system.file("extdata", "16_2175571f05_dc_000045.laz", package="ALSbridges")
las = readLAS(LASfile)

bridge = find_bridges(las, th_planarity = 50, th_angle = 5, buffer = 130, th_hag = 0.9)
bridge

plot(header(las))
plot(bridge, add = TRUE)

plot(las) |> add_hulls3d(bridge, col = "red")

las = classify_bridges(las, bridge, upper_buffer = 25)
plot(las, color = "Classification")
```
![](https://github.com/metafor-ulaval/ALSbridges/blob/main/inst/extdata/screenshot.png?raw=true)
