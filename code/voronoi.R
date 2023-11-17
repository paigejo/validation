# Using Voronoi cells to estimate distributional weights

getVCellAreas = function(pts, domainPoly=NULL) {
  require(deldir)
  
  if(is.null(domainPoly)) {
    domainPolyInds = chull(pts)
    domainPoly = pts[c(domainPolyInds, domainPolyInds[1]),]
  }
  
  out = deldir(pts)
  tiles = tile.list(out, clipp=list(x=domainPoly[,1], y=domainPoly[,2]))
  area = sapply(tiles, function(x) {x$area})
  vcellPolys = sapply(tiles, function(x) {cbind(c(x$x, x$x[1]), c(x$y, x$y[1]))})
  
  list(vcellPolys=vcellPolys, area=area)
}

getVCellRates = function(pts, popDistn, domainPoly=NULL) {
  # get Voronoi cells around each point
  out = getVCellAreas(pts, domainPoly)
  vcells = out$vcellPolys
  areas = out$area
  
  # get info from population distribution
  pixelArea = popDistn$pixelArea
  locs = popDistn$locs
  pixelRates = popDistn$truth
  
  # for each cell:
  popWeights = numeric(length(vcells))
  for(i in 1:length(vcells)) {
    # integrate population distribution in the cell
    thisCellPoly = vcells[[i]]
    inCell = in.poly(locs, thisCellPoly)
    thisPts = pts[inCell,]
    thisRates = pixelRates[inCell]
    
    popWeights[i] = sum(thisRates)
  }
  popProps = popWeights/sum(popWeights)
  
  # return results
  list(vcellPolys=vcellPolys, area=area, cellRate=popProps)
}


