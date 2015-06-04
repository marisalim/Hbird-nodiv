#Code for making distribution inputs for nodiv analysis with hummingbird range maps

MLwd = "D:/MarisaLimfiles/Make_sitexsp_matrix/for sitexsp"
setwd(MLwd)

library(raster)
library(maptools)

myvars <- "bio1.bil"
myvars <- stack(myvars)
plot(myvars)
myvars

myvars2 <- aggregate(myvars, fact=6) #this makes the res=c(1,1)
myvars2
plot(myvars2)

# crop extent to just N. America and S. America
myex<-extent(c(-165, -30, -60, 75))
myvarscrop <- crop(myvars2, myex)
class(myvarscrop) # this is a rasterbrick
res(myvarscrop)
myvarscrop
plot(myvarscrop)

# make presence/absence 1/0 matrix
siteXsp <- function(){
  #lists all NatureServe files with .shp extension name (n=338 files)
  NSfolder = "Trochilidae"
  shpfilenames <- list.files(paste(NSfolder,sep="/"),pattern=".shp",full.name=T,recursive=T)
  
  # use this statement when testing subset of the data, otherwise comment out
#   shpfilenames <- shpfilenames[1:5]
  
  # this gets the sitexsp matrix of 0s and 1s
  rasters <- list()
  d <- data.frame("Species"="species")
  for (i in 1:length(shpfilenames)){
    speciesrange <- readShapePoly(shpfilenames[i]) 
    rasters[[i]] <- rasterize(speciesrange, myvarscrop, field=1, background=0) 
    keepname <- data.frame(Species=speciesrange$SCINAME[1])
    d <- rbind(d, keepname)
    print(i)
  }
  s <- stack(rasters)
  x <- as.data.frame(s)
  d <- d[-1,]
  drename <- strsplit(x=as.character(d), split = " ")
  drename2 <- sapply(drename, function(x){
    paste(x[[1]], x[[2]], sep="_")
  })
  colnames(x) <- drename2
  
  # #   use these to check that the ranges are properly mapping onto the extent
  #   plot(s) # note: this plots the map of the species on the overall extent
  #   plot(x) # note: this is really slow, it just shows that there are 0s and 1s
  write.csv(x, "hbird_sitexsp_matrix.csv")
}
siteXsp()

mysitexspmat <- read.csv("hbird_sitexsp_matrix.csv") # really big file!! will take some time to load
head(mysitexspmat)
dim(mysitexspmat)

#check that there are presences!
length(mysitexspmat[mysitexspmat$Urosticte_ruficrissa == "1", 339])
length(mysitexspmat[mysitexspmat$Metallura_tyrianthina == "1", 251])
length(mysitexspmat[mysitexspmat$Selasphorus_rufus == "1", 312])


# Get cell coords
cell_coords <- function(){
  # get cell IDs
  cellcoords <- xyFromCell(myvarscrop, cell=1:ncell(myvarscrop)) 
  # Then, use the xy coords where species is present to get the cell position IDs in the var raster
  cellIDs <- cellFromXY(myvarscrop, cellcoords)
  
  # this gets the cellIDs and XY coordinates
  mycoords <- data.frame("cellIDs"=cellIDs, "X"=cellcoords[,1], "Y"=cellcoords[,2])
  
  write.csv(mycoords, "cell_coords.csv")
}
cell_coords()

mycellcoords <- read.csv("cell_coords.csv")
head(mycellcoords)

# Make shape file of extent
writeRaster(myvarscrop, "myvarscrop.grd")



