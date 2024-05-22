# Using Voronoi cells to estimate distributional weights

getVCells = function(pts) {
  require(deldir)
  deldir(pts)
}

getCVCells = function(pts, domainPoly=NULL, maxit=1, tol=1e-4, stopcrit = c("maxit", "change"), ...) {
  require(deldir)
  stopcrit = match.arg(stopcrit)
  VCells = deldir(pts)
  
  if(is.null(domainPoly)) {
    domainPolyInds = chull(pts)
    domainPoly = pts[c(domainPolyInds, domainPolyInds[1]),]
  }
  
  # clip Voronoi cells to be within the domain polygon
  tiles = tile.list(VCells, clipp=list(x=domainPoly[,1], y=domainPoly[,2]))
  
  # calculate centroidal Voronoi cells and put them in the typical format
  CVCells = cvt(tiles, tol=tol, maxit=maxit, stopcrit=stopcrit, ...)
  deldir(CVCells$centroids)
}

getVCellAreas = function(pts=NULL, domainPoly=NULL, VCells=NULL, useCVC=FALSE, ...) {
  require(deldir)
  
  if(is.null(domainPoly) && !is.null(pts)) {
    domainPolyInds = chull(pts)
    domainPoly = pts[c(domainPolyInds, domainPolyInds[1]),]
  }
  
  if(is.null(VCells)) {
    if(useCVC) {
      VCells = getCVCells(pts, domainPoly=domainPoly, ...)
    } else {
      VCells = deldir(pts)
    }
  }
  
  if(is.null(domainPoly)) {
    tiles = tile.list(VCells)
  } else {
    tiles = tile.list(VCells, clipp=list(x=domainPoly[,1], y=domainPoly[,2]))
  }
  
  area = sapply(tiles, function(x) {x$area})
  vcellPolys = sapply(tiles, function(x) {cbind(c(x$x, x$x[1]), c(x$y, x$y[1]))})
  
  if(useCVC && !is.null(pts)) {
    # match each point with respective tile/centroid
    centroids = do.call("rbind", lapply(tiles, function(x) {x$pt}))
    dists = rdist(pts, centroids)
    area2ptI = apply(dists, 1, which.min)
    ptArea = area[area2ptI]
    nPerArea = sapply(1:length(area), function(i) {sum(area2ptI == i)})
    ptNPerArea = nPerArea[area2ptI]
  } else {
    centroids = pts
    ptArea = area
    area2ptI = NULL
    nPerArea = NULL
    ptNPerArea = NULL
  }
  
  list(vcellPolys=vcellPolys, area=area, ptArea=ptArea, centroids=centroids, 
       area2ptI=area2ptI, nPerArea=nPerArea, ptNPerArea=ptNPerArea)
}

plotVCells = function(VCells, add=TRUE, xlim=NULL, ylim=NULL, main="", 
                      xlab="", ylab="", ...) {
  polys = getVCellAreas(VCells=VCells)$vcellPolys
  
  if(!add) {
    allXs = sapply(VCells, rbind)
    xlim = range(allXs[,1])
    ylim = range(allXs[,2])
  }
  
  for(i in 1:length(VCells)) {
    thisPoly = polys[[i]]
    
    if(i == 1 && !add) {
      plot(allXs, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
    }
    
    polygon(thisPoly, ...)
  }
}

# avg: whether to average or just integrate within each cell
integrateVCells = function(pts=NULL, gridLocs, valsAtGrid=1, domainPoly=NULL, 
                           VCells=NULL, adjustIntegralsFromArea=(length(valsAtGrid) != 1)) {
  
  if(length(valsAtGrid) == 1) {
    valsAtGrid = rep(valsAtGrid, nrow(gridLocs))
  }
  
  # get Voronoi cells around each point
  out = getVCellAreas(pts, domainPoly, VCells=VCells)
  vcells = out$vcellPolys
  trueAreas = out$area
  
  # get info from population distribution
  locs = gridLocs
  deltas = apply(locs, 2, function(x) {diff(sort(unique(x))[1:2])})
  pixelArea = prod(deltas)
  
  # for each cell:
  ints = numeric(length(vcells))
  avgs = numeric(length(vcells))
  allInCell = matrix(FALSE, nrow=nrow(locs), ncol=length(vcells))
  for(i in 1:length(vcells)) {
    # integrate values in the cell
    thisCellPoly = vcells[[i]]
    inCell = in.poly(locs, thisCellPoly)
    allInCell[,i] = inCell
    theseVals = valsAtGrid[inCell]
    
    # integrate values in the Voronoi cell (later adjusted)
    ints[i] = sum(theseVals)
    
    # get average value in the Voronoi cell
    approxArea = sum(inCell)
    avgs[i] = ints[i] / approxArea
    
    # adjust integral based on the true versus approximate area of the cell. This 
    # may work poorly if extreme values are at the edges of the Voronoi cell
    if(adjustIntegralsFromArea) {
      ints[i] = ints[i] * trueAreas[i]/approxArea
    }
  }
  
  # return results
  list(ints=ints, avgs=avgs, vcellPolys=vcells, area=area, allInCell=allInCell)
}

# integrates one VCell density estimate over another set of VCells
# avg: whether to average or just integrate within each cell
integrateVCells2 = function(pts1=NULL, VCells1=NULL,
                            pts2=NULL, VCells2=NULL, domainPoly=NULL) {
  
  if(is.null(domainPoly) && (!is.null(pts1) && !is.null(pts2))) {
    domainPoly = chull(rbind(pts1, pts2))
  }
  
  # get Voronoi cell polygons and areas
  out1 = getVCellAreas(pts1, domainPoly, VCells=VCells1)
  vcells1 = out1$vcellPolys
  trueAreas1 = out1$area
  
  out2 = getVCellAreas(pts2, domainPoly, VCells=VCells2)
  vcells2 = out2$vcellPolys
  trueAreas2 = out2$area
  rates2 = 1/trueAreas2
  
  # convert polygons to sf objects
  vcellsSF1 = lapply(vcells1, function(x) {st_polygon(list(x))})
  vcellsSF2 = lapply(vcells2, function(x) {st_polygon(list(x))})
  vcellsSF1List = st_sfc(vcellsSF1)
  vcellsSF2List = st_sfc(vcellsSF2)
  # each element of intersects is a list of VCells2 intersecting with the given VCells1
  intersects = st_intersects(vcellsSF1List, vcellsSF2List)
  
  # for each VCell1:
  ints = numeric(length(vcells1)) # initialize at 0
  avgs = numeric(length(vcells1))
  for(i in 1:length(vcells1)) {
    # get which VCells2 intersect with this VCell1
    thisAllIntersectsI = intersects[[i]]
    
    # for each VCell2 intersection, get proportion of area in this VCell1
    for(j in 1:length(thisIntersects)) {
      # get intersection
      vcell2I = thisAllIntersectsI[j]
      thisIntersection = st_intersection(vcellsSF1[[i]], vcellsSF2[[vcell2I]])
      
      # get area of the intersection
      intersectArea = st_area(thisIntersection)
      
      # add contribution to the rate integral in this VCell1
      ints[i] = ints[i] + rates2[vcell2I] * intersectArea
    }
    
    # divide integral by area of VCell1 to get average in VCell1
    avgs[i] = ints[i] / trueAreas[i]
  }
  
  # return results
  list(ints=ints, avgs=avgs)
}

getVCellRatesOnGrid = function(gridLocs, VCells, domainPoly) {
  # get Voronoi cell polygons
  out = getVCellAreas(pts=NULL, domainPoly, VCells=VCells)
  vcells = out$vcellPolys
  
  # get which grid points (and how many) are associated with each cell
  out = integrateVCells(pts=NULL, gridLocs=gridLocs, domainPoly=domainPoly, 
                        VCells=VCells)$ints
  nPerCell = out$ints
  inCell = out$allInCell
  
  # calculate rates
  rates = numeric(nrow(gridLocs))
  for(i in 1:length(nPerCell)) {
    thisI = inCell[,i]
    rates[thisI] = 1/nPerCell[i]
  }
  
  rates
}

