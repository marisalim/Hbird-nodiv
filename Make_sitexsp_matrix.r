#Code for making distribution inputs for nodiv analysis with hummingbird range maps

MLwd = "D:/MarisaLimfiles/Make_sitexsp_matrix/for sitexsp"
setwd(MLwd)

library(raster)
library(maptools)
library(reshape2)

myvars <- "bio1.bil"
myvars <- stack(myvars)
plot(myvars)
myvars

myvars2 <- aggregate(myvars, fact=6) #fact=6 makes the res=c(1,1); fact=3 makes the res=c(0.5,0.5)
myvars2
plot(myvars2)

# crop extent to just N. America and S. America
myex<-extent(c(-165, -30, -60, 75))
myvarscrop <- crop(myvars, myex)
class(myvarscrop) # this is a rasterlayer
res(myvarscrop)
myvarscrop
plot(myvarscrop)

# make presence/absence 1/0 matrix
# this includes all the space where hummingbirds don't occur - e.g., ocean, so a lot of absences that I don't need
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

head(mysitexspmat[mysitexspmat$Loddigesia_mirabilis == "1",230]) #range too small? has 0 presences...the mean of the space must be 0

#try Ben's code:
siteXsp_Ben <- function(){
  #create blank raster
  r <- raster(myvarscrop)
  
  #Turn each shapefile into a raster map
  species_rasterize <- list()
  
  #lists all NatureServe files with .shp extension name (n=338 files)
  NSfolder = "Trochilidae"
  shpfilenames <- list.files(paste(NSfolder,sep="/"),pattern=".shp",full.name=T,recursive=T)
  d <- data.frame("Species"="species")
  #loop through all files
  for(x in 1:length(shpfilenames)){
    #load file
    try(sp_shp <- readShapePoly(shpfilenames[[x]], delete_null_obj=TRUE))
    
    #turn to raster
    sp_ras <- rasterize(sp_shp,r)
    
    #which cells are 1's?
    presence_cells <- Which(sp_ras, cell=TRUE)
    
    #check extent, return NA if no presences
    if(length(presence_cells)==0){
      species_rasterize[[x]] <- (NA); next}
    
    #get xy of presence cells
    species_rasterize[[x]] <- presence_cells
    
    #keep species names for naming columns
    keepname <- data.frame(Species=sp_shp$SCINAME[1])
    d <- rbind(d, keepname)
    
    print(x)
  }
  
  d <- d[-1,]
  drename <- strsplit(x=as.character(d), split = " ")
  drename2 <- sapply(drename, function(x){
    paste(x[[1]], x[[2]], sep="_")
  })
  names(species_rasterize) <- drename2
  
  #remove species with just NA (length==1)
  species_full <- species_rasterize[!sapply(species_rasterize,length)==1]
  
  #yields how many species?
  length(species_full)
  
  #melt list into dataframe
  m.species <- melt(species_full)
  
  sitexsppmat <- t(as.data.frame.array(table(m.species$L1, m.species$value)))
  write.csv(sitexsppmat, "sitexspmat.csv")
}
siteXsp_Ben()

mysitexspmat <- read.csv("sitexspmat.csv") 
head(mysitexspmat)
dim(mysitexspmat)
rownames(mysitexspmat) <- mysitexspmat$X
colnames(mysitexspmat)[1] <- "cellIDs"
length(mysitexspmat[mysitexspmat$Calypte_anna == "1", 65])
head(mysitexspmat[mysitexspmat$Urosticte_ruficrissa == "1", 324])
head(mysitexspmat[mysitexspmat$Loddigesia_mirabilis == "1",217])

# Get cell coords
cell_coords <- function(){
  # get cell IDs
  cellcoords <- xyFromCell(myvarscrop, cell=1:ncell(myvarscrop)) 
  # Then, use the xy coords where species is present to get the cell position IDs in the var raster
  cellIDs <- cellFromXY(myvarscrop, cellcoords)
  
  # this gets the cellIDs and XY coordinates
  mycoords <- data.frame("cellIDs"=cellIDs, "X"=cellcoords[,1], "Y"=cellcoords[,2])
  
  # make it match the cellIDs in the sitexsp matrix
  mycells <- data.frame(cellIDs=mysitexspmat$cellIDs)
  dim(mycells)
  dim(mycoords)
  mycoords2 <- merge(mycoords, mycells, by.y="cellIDs")
  dim(mycoords2)
  # head(mycoords2)
  
  write.csv(mycoords, "cell_coords.csv")
  write.csv(mycoords2, "cell_coords2.csv")
}
cell_coords()

mycellcoords <- read.csv("cell_coords.csv")
head(mycellcoords)

# Make shape file of extent
writeRaster(myvarscrop, "myvarscrop.grd", overwrite=TRUE)



