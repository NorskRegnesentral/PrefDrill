# Script with miscellaneous useful functions

logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  res = 1/(1+exp(-x))
  res[x >= 37] = 1
  res[x <= -101] = 0
  res
}

# Use bilinear interpolation on a grid, where gridDat has long format
# 
# inputs:
# pts: easting/northing coordinates of pts to extract bilinearly interpolated 
#      value for based on gridDat
# gridDat: 3 columns, with first two being names "east" and "north", third being 
#          what to interpolate
# 
# outputs: vector of interpolated values at input pts
bilinearInterp = function(pts, gridDat) {
  require(akima)
  
  # use bilinear interpolation to get seismic estimate at that point
  # First must convert to akima format into matrix grid form for z
  eastGrid = sort(unique(gridDat[,2]))
  northGrid = eastGrid
  sortI1 = order(gridDat[,1]) # east
  gridDat = gridDat[sortI1,]
  sortI2 = order(gridDat[,2]) # north
  gridDat = gridDat[sortI2,]
  
  z = matrix(gridDat[,3], nrow=length(eastGrid), ncol=length(northGrid))
  
  seismicEst = akima::bilinear(x=eastGrid, y=northGrid, z=z, x0 = pts[,1], y0=pts[,2])$z
  seismicEst
}




