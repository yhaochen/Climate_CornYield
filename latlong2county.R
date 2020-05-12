# Found from https://stackoverflow.com/questions/13316185/r-convert-zipcode-or-lat-long-to-county

library(sp)
library(maps)
library(maptools)
# The single argument to this function, pointsDF, is a data.frame in which:
#   - column 1 contains the longitude in degrees (negative in the US)
#   - column 2 contains the latitude in degrees

latlong2county <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per county (plus DC, minus HI & AK)
  countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
  countys_sp <- map2SpatialPolygons(countys, IDs=IDs,
                                    proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, countys_sp)
  
  # Return the county names of the Polygons object containing each point
  countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
  countyNames[indices]
}

# Test the function using points in Wisconsin and Oregon.
testPoints <- data.frame(x = c(-77.4805), y = c(40.85))

latlong2county(testPoints)

