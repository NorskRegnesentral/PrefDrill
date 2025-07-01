# script with functions for reading RMS surfaces


readSurfaceRMS = function(filename) {
  # << "  -996" << " "
  # << std::setw(5)  << grid.GetNY()     << "  "
  # << std::setw(12) << grid.GetXInc()   << "  "
  # << std::setw(12) << grid.GetYInc()   << "\n"
  # << std::setprecision(4)
  # << std::setw(12) << grid.GetXStart() << "  "
  # << std::setw(12) << xend             << "  "
  # << std::setw(12) << grid.GetYStart() << "  "
  # << std::setw(12) << yend             << "\n "
  # << std::setw(5)  << grid.GetNX()     << "  "
  # << std::setw(12) << Conversion::RadiansToDegrees(grid.GetRotation()) << "  " // rotation origin
  # << std::setw(12) << grid.GetXStart() << "  "
  # << std::setw(12) << grid.GetYStart() << "\n"
  # << "   0   0   0   0   0   0   0\n";
  
  # read file, separate header and values, remove NAs from values
  tab = read.table(filename, fill=TRUE)
  tab = as.matrix(tab)
  header = tab[1:3,]
  vals = tab[-(1:4),]
  vals = c(t(vals))
  vals = vals[!is.na(vals)]
  
  # extract information from header
  ny = header[1,2]
  xinc = header[1,3]
  yinc = header[1,4]
  xstart = header[2,1]
  xend = header[2,2]
  ystart = header[2,3]
  yend = header[2,4]
  nx = header[3,1]
  rotation = header[3,2]
  xstart = header[3,2] # not sure why there are multiple xstarts
  ystart = header[3,3] # not sure why there are multiple ystarts
  
  if(rotation != 0) {
    stop("rotation is nonzero, which isn't yet supported")
  }
  
  xgrid = seq(from=xstart, to=xend, l=nx)
  ygrid = seq(from=ystart, to=yend, l=ny)
  surfMat = matrix(vals, nrow=nx, ncol=ny)
  
  # convert to long format
  surfFrame = data.frame(east = rep(xgrid, ny), north=rep(ygrid, each=nx), 
                         seismicEst=c(surfMat))
  
  list(xgrid=xgrid, ygrid=ygrid, surfMat=surfMat, surfFrame=surfFrame, 
       nx=nx, ny=ny, xstart=xstart, ystart=ystart, xend=xend, yend=yend)
  
  
}





