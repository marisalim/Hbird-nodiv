# Code to run nodiv on species range map data (site x sp matrix)
# Marisa Lim (c)2015

# load libraries
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

# set working directory
#MLwd = "D:/MarisaLimfiles/Make_sitexsp_matrix/"
MLwd = "C:/Users/mcwlim/Dropbox/resultsfromSarahcomputer/nodesigfiles"
setwd(MLwd)

# ------------ Building distrib_data and nodiv_data object ------------
# Read in site by species matrix
#hb_commat <- read.csv("for sitexsp/hbird_sitexspmat_sub.csv") #a 0/1 site-by-species matrix
hb_commat <- read.csv("hbird_sitexspmat_sub.csv")
dim(hb_commat)
hb_commat[1:5,1:4]
rownames(hb_commat) <- hb_commat$X
hb_commat2 <- hb_commat[,-1]

# Read in cell coordinate information
#hb_coords <- read.csv("for sitexsp/cell_coords_rounded.csv")[,-1] # a data.frame with one row for each site, and lat/long coordinates in 2 of the columns
hb_coords <- read.csv("cell_coords_rounded.csv")[,-1]
head(hb_coords)
#colnames(hb_coords) <- c("cellIDs", "X", "Y")

# The coords do not need to have the same length or be ordered the same way
hummers <- distrib_data(hb_commat2, hb_coords) 

# Describing the data set 
hummers
summary(hummers)
names(hummers)

# create a species richness plot as diagnostic
plot(hummers)  

# Read in geographic extent shape file
#myshape <- raster("for sitexsp/myvarscrop.grd")
myshape<- raster("myvarscrop.grd")
plot(myshape,col="grey")
hummers <- add_shape(hummers, myshape)
plot(hummers)

# ------------ Add the phylogenetic information ------------
#tree <- read.tree("humtree280.tre") 
tree <- read.tree("C:/Users/mcwlim/Desktop/StonyBrook/GrahamLab/Dissertation idea materials/THESIS PROJECTS/Hummingbird_Genomics/humtree280.tre")
hummers <- nodiv_data(tree, hummers) 
summary(hummers)

#look at node numbers
plot(hummers$phylo)
nodelabels(bg="white")

# Now the default plot also shows the phylogeny
jpeg("Hbird_richness.jpg", height=6, width=6, units="in", res=600)
plot_richness(hummers, main="Species richness")
dev.off()
jpeg("Hbird_noderichness.jpg", height=6, width=10, units="in", res=600)
par(mfrow=c(1,2))
plot_node(hummers, node=329, main="Species richness: Andean clade vs. rest")
plot_node(hummers, node=272, main="Species richness at node: Topazes,Hermits vs. rest")
dev.off()

jpeg("Hbird_noderichness_outofSA.jpg", height=6, width=6, units="in", res=600)
plot_node(hummers,node=416, main="Species richness; out of SA")
dev.off()
jpeg("Hbird_noderichness_Coqvs.Brill.jpg", height=6, width=6, units="in", res=600)
plot_node(hummers,node=330, main="Species richness; Coquettes vs. Brilliants")
dev.off()

# ------------ The most basic summary metrics are then available: ------------
head(species(hummers))
Nspecies(hummers)
head(sites(hummers))
Nsites(hummers)

hist(occupancy(hummers))

plot_species(hummers, "Archilochus_colubris")
plot_species(hummers, "Calypte_anna")
plot_species(hummers, "Metallura_phoebe")

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

#Optionally, you can restrict the analysis to only sites with at least 4 species
# hummers <- subsample(hummers, sites = richness(hummers) > 3)

# ------------ Run nodiv! ------------
humm_nodiv <- Node_analysis(hummers, repeats = 200, method = "rdtable") #TODO: how many repeats?
summary(humm_nodiv)
par(mfrow=c(1,1))
plot(humm_nodiv, label=nodenumbers(humm_nodiv), col="HMrainbow")

jpeg("Hbird_GND_geography.jpg", height=8, width=15, units="in", res=600)
plot(humm_nodiv, label=humm_nodiv$GND, col=heat.colors(10), direction="upwards")
dev.off()
#plot SOS values for a given node
# check the node numbers if you change resolution of data
jpeg("Hbird_SOS_Andeancladevsrest.jpg", height=6, width=6, units="in", res=600)
plotSOS(humm_nodiv, 329, main="Andean clade vs. rest", col="HMrainbow") 
dev.off()
jpeg("Hbird_SOS_Basevsrest.jpg", height=6, width=6, units="in", res=600)
plotSOS(humm_nodiv, 272, main="HermitsTopazes vs. rest", col="HMrainbow")
dev.off()
jpeg("Hbird_SOS_BeesGemsvsEms.jpg", height=6, width=6, units="in", res=600)
plotSOS(humm_nodiv, 416, main="BeesGems vs. Emeralds", col="HMrainbow")
dev.off()

jpeg("Hbird_SOS_hists.jpg", height=6, width=8, units="in", res=600)
par(mfrow=c(1,2))
hist(SOS(humm_nodiv, 329))
hist(SOS(humm_nodiv, 272))
dev.off()

GND(humm_nodiv)

# ------------ code for running without Archilochus colubris -----------
#hb_commat2 without A. colubris
# hb_commat_noAc <- hb_commat2[,!colnames(hb_commat2)=="Archilochus_colubris"]
# hummers_noAc <- distrib_data(hb_commat_noAc, hb_coords)
# hummers_noAc <- add_shape(hummers_noAc, myshape)

# hummers_noAc <- nodiv_data(tree, hummers_noAc)
# summary(hummers_noAc)
# plot_richness(hummers_noAc, main="Species richness")
# plot_node(hummers_noAc, node=303)
# plot_node(hummers_noAc, node=250)

# humm_nodiv_noAc <- Node_analysis(hummers_noAc, repeats = 200, method = "rdtable") #TODO: how many repeats?
# summary(humm_nodiv_noAc)
# plot(humm_nodiv_noAc, label=nodenumbers(humm_nodiv_noAc), col="white")
# 
# jpeg("Hbird_GND_geographynoAc.jpg", height=8, width=10, units="in", res=600)
# plot(humm_nodiv_noAc, label=humm_nodiv_noAc$GND, col=heat.colors(10), direction="upwards")
# dev.off()
# #plot SOS values for a given node
# jpeg("Hbird_SOS_AndeancladevsrestnoAc.jpg", height=6, width=6, units="in", res=600)
# plotSOS(humm_nodiv_noAc, 303, main="node 303: Andean clade vs. rest") 
# dev.off()
# jpeg("Hbird_SOS_BasevsrestnoAc.jpg", height=6, width=6, units="in", res=600)
# plotSOS(humm_nodiv_noAc, 250, main="node 250: HermitsTopazes vs. rest")
# dev.off()
# 
# jpeg("Hbird_SOS_histsnoAc.jpg", height=6, width=8, units="in", res=600)
# par(mfrow=c(1,2))
# hist(SOS(humm_nodiv_noAc, 303))
# hist(SOS(humm_nodiv_noAc, 250))
# dev.off()



