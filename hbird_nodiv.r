MLwd = "D:/MarisaLimfiles/Make_sitexsp_matrix/"
setwd(MLwd)

# # Installing the package
# install.packages("nodiv") 
#OR
library(devtools)
install_github("mkborregaard/nodiv") #Gives the newest version - is probably LESS buggy than CRAN version

# # load the package
library(nodiv)
# help(package = nodiv)

library(ape)
library(raster)

# Building distrib_data and nodiv_data object
hb_commat <- read.csv("for sitexsp/hbird_sitexsp_matrix.csv") #a 0/1 site-by-species matrix
dim(hb_commat)
hb_commat2 <- hb_commat[,-1]
hb_coords <- read.csv("for sitexsp/cell_coords.csv")[,-1] # a data.frame with one row for each site, and lat/long coordinates in 2 of the columns
head(hb_coords)
colnames(hb_coords) <- c("cellIDs", "X", "Y")

# hb_commat2test <- hb_commat2[250613:258933,]
# hb_coordstest <- hb_coords[250613:258933,]
# hb_test <- distrib_data(hb_commat2test,hb_coordstest)

# The coords do not need to have the same length or be ordered the same way
hummers <- distrib_data(hb_commat2, hb_coords) #FIXME: CRASHES with 0.1degree raster! I think the data is too big; it's fine when I change the res to 1degree

# Describing the data set
hummers
summary(hummers)
names(hummers)

# create a species richness plot as diagnostic
plot(hummers)  

myshape <- raster("for sitexsp/myvarscrop.grd")
plot(myshape)
hummers <- add_shape(hummers, myshape)
plot(hummers)

# We can also add the phylogenetic information
tree <- read.tree("humtree280.tre") 
hummers <- nodiv_data(tree, hummers) 
summary(hummers)

# Now the default plot also shows the phylogeny
plot(hummers)
jpeg("Hbird_richness.jpg", height=6, width=6, units="in", res=600)
plot_richness(hummers, main="Species richness")
dev.off()

# The most basic summary metrics are then available:
head(species(hummers))
Nspecies(hummers)
head(sites(hummers))
Nsites(hummers)

hist(occupancy(hummers))

plot_species(hummers, "Chalcostigma_herrani")
plot_species(hummers, 16) # species, sites and nodes can all be identified as numbers or names
i <- 16
plot_species(hummers, i <- i + 1)

occurrences(hummers, "Chalcostigma_herrani")
assemblage(hummers, site = 6)
assemblage(hummers, site = 6, value = "names")

# Tools for working with nodes
plot_nodes_phylo(Node_occupancy(hummers), hummers$phylo, cex = 0.8)

Node_species(hummers, node = 150)
Sister(150, hummers)
Parent(150, hummers)
Descendants(150, hummers)

# An example of the nodiv procedure

#Optionally, you can restrict the analysis to only sites with at least 4 species
# hummers <- subsample(hummers, sites = richness(hummers) > 3)

humm_nodiv <- Node_analysis(hummers, repeats = 200, method = "rdtable") #TODO: how many repeats?
summary(humm_nodiv)
plot(humm_nodiv, label=nodenumbers(humm_nodiv), col="white")

jpeg("Hbird_GND_geography.jpg", height=8, width=10, units="in", res=600)
plot(humm_nodiv, label=humm_nodiv$GND, col=heat.colors(10), direction="upwards")
dev.off()
#plot SOS values for a given node
jpeg("Hbird_SOS_Andeancladevsrest.jpg", height=6, width=6, units="in", res=600)
plotSOS(humm_nodiv, 304, main="node 304: Andean clade vs. rest") 
dev.off()
jpeg("Hbird_SOS_Basevsrest.jpg", height=6, width=6, units="in", res=600)
plotSOS(humm_nodiv, 251, main="node 251: HermitsTopazes vs. rest")
dev.off()

hist(SOS(humm_nodiv, 304))
GND(humm_nodiv)


