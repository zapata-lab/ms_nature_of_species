##############################
# State of taxonomic species #
##############################

# This R script is written by members of the Zapata Lab

---------------------------------------------------------------------------------
# This script runs the analysis on the state of taxonomic species using phenotypic data. For details on the analysis and interpretation of results, please read the Methods and Results in our manuscript

# We use n-cubes (i.e., hypervolumes) to describe geometrically each taxonomic species in phenospace. To build each n-cube, we use the as vertices the minimum and maxium values of all traits used to delimit each taxonomic species based on the taxonomic monograph.

# Once we generate the n-cubes for each species, we ask two main questions:
# Do n-cubes (i.e., taxonomic species) occupy distinct areas of phenospace?
# Do specimens hypothesized to belong to a given taxonomic species fall inside or outside of its corresponding n-cube (i.e., taxonomic species)?

# There are two examples provided below.
# The first runs an example using mock-data, simple enough to calculate and draw-out by hand to gain an intuition as to what the code does and this analysis works.
# The second uses CladeI as an example of how this analysis was performed using our empirical data. For other clades, please modify the necessary lines acordingly

--------------------------------------------------------------------------------



# Preliminaries
--------------------------------------------------------------------
library(tidyverse)
library(geometry)
library(reshape2)
library(RColorBrewer)
library(plot3D)

#################################
#
# MOCK DATA
#
#################################
# This example is intended to be drawn by hand, as a complementary visualization of what's happening at each step

# Given one species with the following min/max measurements:
# spp1: trait1 = 2-6, trait2 = 3-13, trait3 = 1-7
# we can build a volume (3 dimensions, or 3-cube) that uses the ranges of values to create the sides of the shape
# the number of vertices for a 3-dimensional n-cube would be equal to 2^3
# So, 8 vertices in total to describe this species; each of these sets of three numbers are the x,y,z coordinates of one corner of the cube
testVertices <- rbind(c(2,3,1),
                      c(6,3,1),
                      c(2,13,1),
                      c(6,13,1),
                      c(2,3,7),
                      c(6,3,7),
                      c(2,13,7),
                      c(6,13,7))

# In the same graphing area, draw out the cube corresponding to these three trait values
# Similarly, spp2 has min/max measurements for the same traits:
# spp2: trait1 = 5-14, trait2 = 5-8, trait3 = 1-3
# The 8 vertices to describe this species would be
testVertices2 <- rbind(c(5,5,1),
                       c(14,5,1),
                       c(5,8,1),
                       c(14,8,1),
                       c(5,5,3),
                       c(14,5,3),
                       c(5,8,3),
                       c(14,8,3))

# to compute the n-cubes, we use the function convhulln from the geometry package
testHull <- convhulln(testVertices)
testHull2 <- convhulln(testVertices2)

# Next we calculate the intersection of these n-cubes when they are plotted in the same morphospace
# we use the function intersectn from the geometry package
# this is where we will extract the volumes of the n-cubes, and the volume of the intersection of the n-cubes
# the volume is calculated as V = l*w*h
testHullIntersects <- intersectn(as.matrix(testVertices), as.matrix(testVertices2))

testHullIntersects$ch1$vol # the volume of the n-cube for spp1 = 240
testHullIntersects$ch2$vol # the volume of the n-cube for spp2 = 54
testHullIntersects$ch$vol # the volume of the interesection of spp1 and spp2 = 6


# For 10 specimens, we take point measurements for the three traits
# here we list the point measurements in order as (trait1, trait2, trait3) and
# given these x,y,z coordinates, we can see whether they fall 'in' or 'out' of the n-cubes for spp1 and spp2
# for each specimen, we provide its position, relative to both spp's n-cubes
testPoints <- rbind(
  c(3,4,2), #specimen 1: in spp1, out spp2
  c(8,8,3), #specimen 2: out spp1, in spp2
  c(13.5,13.5,3), #specimen 3: out spp1, out spp2
  c(5,8,3), #specimen 4: in spp1, in spp2
  c(8,4,5), #specimen 5: out spp1, out spp2
  c(5,14,4), #specimen 6: out spp1, out spp2
  c(14,11,6), #specimen 7: out spp1, out spp2
  c(5,8,1), #specimen 8: in spp1, in spp2
  c(5,8,7),  #specimen 9: in spp1, out spp2
  c(3,7,8)) #specimen 10: out spp1, out spp2

# We use the function inhulln to test if each specimen lies within (T or F) of n-cube for spp1, followed by spp2
# these results match those listed next to each specimen above
testInHulln <- inhulln(testHull, testPoints)
testInHulln2 <- inhulln(testHull2, testPoints)

# to test any set of coordinates, plug in the x,y,z values here
testPoints <- testPoints[,c(3,2,1)]
testInHulln <- inhulln(testHull, testPoints)



#################################
#
# EMPIRICAL DATA
#
#################################

# Preliminaries
--------------------------------------------------------------------

# define the function to calculate the intersection of multiple species' hull and record the information we want - only the volume of sp1, sp2, and the volume of the intersection of the two species n-cubes
intersection_by_pair = function(i,j){
  inthull <- data.frame(NA)
  chullsp1 <- data.frame(NA)
  chullsp2 <- data.frame(NA)
  result = intersectn(i,j)
  chullsp1[1,1] = result$ch1$vol #store the volume of the first species n-cube
  chullsp2[1,1] = result$ch2$vol #store the volume of the second species n-cube
  inthull[1,1] = result$ch$vol #store the volume of the intersection hull
}


# Read in data
## MONOGRAPH translations: these are pre-formatted long, tidy data, that record the min and max measurements for all quantitative traits recorded in species descriptions. Change the input file for other clades, accordingly
indataTemp <- read_csv("../../data/CladeI/CladeI_sppHull_rangeAdj.csv")

## SPECIMEN data: Read in specimen sample data (quantitative measures taken from specimens)
# filter samples by clade membership - keeping only samples that belong to species of the same clade; filter out sterile samples (those with no floral measurements)
samples <- read_csv("../../data/morphometrics_mspub.csv")
samples_CladeI_temp = samples %>% 
  dplyr::filter(CladeID == "I", 
                PETLEN != 0) %>% #ensures the retention of only those specimens with both leaf AND floral measurements
  dplyr::arrange(Taxon)


## Filter datasets
# For the monograph translations and the specimen data, filter to retain traits to be used for the analysis from the larger set 
# NOTE: The trait order of the monograph traits MUST MATCH the order of specimen traits (test for the match is included below)

# filter the monograph translations
indataAll <- indataTemp %>% 
  dplyr::filter(morpho_trait %in% c("CALLOBLEN",
                                    "FILLEN",
                                    "LAMLEN",
                                    "LAMWID",
                                    "OVALEN",
                                    "PEDLEN",
                                    "PETLEN",
                                    "PETWID",
                                    "STYLEN"))
# filter the specimen data
samples_CladeI_All = samples_CladeI_temp %>% 
  dplyr::select(CALLOBLEN,
                FILLEN,
                LAMLEN,
                LAMWID,
                OVALEN,
                PEDLEN,
                PETLEN,
                PETWID,
                STYLEN) 


# Order the trait columns of the specimens to match that of the monograph translations
samples_CladeI_All <- samples_CladeI_All[,order(match(colnames(samples_CladeI_All), unique(indataAll$morpho_trait)))]

# Confirm the order - this should result in TRUE 
all(colnames(samples_CladeI_All) == unique(indataAll$morpho_trait))


# CALCULATE N-CUBES
--------------------------------------------------------------------

# parse trait measurements, find VERTICES of species n-cubes, and create a list of vertices for each species based on those traits
# format the data from the monograph translation - parse traits by species and trait
# find the VERTICES of the n-cube for each species (these increase in size with the number of dimensions used - 2^number of dimensions)
all_sp_list_All = lapply(split(indataAll, indataAll$species), function(x) lapply(split(x[3:4], x$morpho_trait), unlist, use.names = FALSE))
list_of_all_vertices_All = lapply(all_sp_list_All, function(x) as.matrix(expand.grid(x)))

# calculate the N-CUBE for each species using the convhulln function in package geometry - depending on the number of vertices, this can take a while to calculate (on the order of tens of minutes)
chull_by_group_All <- lapply(list_of_all_vertices_All, function(x) convhulln(x))

# find and calculate the volume of the INTERSECTION of all pairwise combinations of species n-cubes - similarly, can take a while to calculate, depending on the number of vertices and the number of n-cubes considered
# uses the function defined at the beginning 
intersectionHulls_All <- lapply(list_of_all_vertices_All, function(x) lapply(list_of_all_vertices_All, function(y) intersection_by_pair(x, y)))

# the above results in a list of volumes of intersections
# $micrantha
# $micrantha$micrantha --- the volume of the intersection of the "micrantha" n-cube within the "micrantha" n-cube. Because this is a self-comparison, the volume of the intersection is equal to the volume of the micrantha n-cube
# [1] 45.5625
# 
# $micrantha$millegrana --- the volume of the intersection of the "millegrana" n-cube within the "micrantha" n-cube. These n-cubes do not overlap (i.e., there is no intersection), so this volume is zero
# [1] 0
# 
# 
# $millegrana
# $millegrana$micrantha --- the volume of the intersection of the "micrantha" n-cube within the "millegrana" n-cube. These n-cubes do not overlap (i.e., there is no intersection), so this volume is zero
# [1] 0
# 
# $millegrana$millegrana --- the volume of the intersection of the "millegrana" n-cube within the "millegrana" n-cube. Because this is a self-comparison, the volume of the intersection is equal to the volume of the millegrana n-cube
# [1] 2250


# ANALYSIS 
--------------------------------------------------------------------
#

#############################################################
# PROPORTION OF OVERLAP: Do n-cubes of different species overlap? 
# By how much?
############################################################# 
# convert the information about the intersection of all combinations of species n-cubes (a list of lists) into a dataframe and store the diagonal element (the species specific n-cube volumes)
intersectnDF_All <- data.frame(t(sapply(intersectionHulls_All,c)))
DF_diag_All <- unlist(diag((as.matrix(intersectnDF_All)), names=FALSE))


# calculate the proportion of overlap, per column, based on the diagonal element residing in that column 
# perform the column-wise calculation on the first column, using the first element in the diagonal list
overlapProp_All <- data.frame()

for(i in 1:length(DF_diag_All)){
  overlapProp_All <- rbind(overlapProp_All, lapply(intersectnDF_All[,i], function(x) x/DF_diag_All[[i]]))
}
rownames(overlapProp_All) <- rownames(intersectnDF_All)


# make the results tidy and long; as.matrix part is essential to get var1 and var2, as well as the value
melted_overlapProp_All <- melt(as.matrix(overlapProp_All))


#############################################################
#
# Matching-prediction tests: Are specimes inside their respective
# n-cube?  True or False ###
#
#############################################################


# for each species' n-cube, ask if each specimen occurs within the n-cube (True) or outside of the n-cube (False)
# will need one 'inHull' assignment for every entity being examined
inHull_All_1 <- inhulln(chull_by_group_All[[1]], as.matrix(samples_CladeI_All))
inHull_All_2 <- inhulln(chull_by_group_All[[2]], as.matrix(samples_CladeI_All))

# bind the True/False results together, and add two new columns that give a unique identifier and indicate assignment expectations (i.e., taxonomic identification of each specimen)
CladeI_hullResults_All <- as.data.frame(cbind(samples_CladeI_temp$SpecimenID,inHull_All_1,inHull_All_2,samples_CladeI_temp$Taxon))
CladeI_hullResults_All <- CladeI_hullResults_All %>% arrange(V4, V1) # order by taxonomic identification, then by specimenID


#############################################################
#
# Visualizations
#
############################################################# 

# melt the data for visualization 
# the second id column will need to be adjusted, depending on how many species are being examined
melted_CladeI_hullResults_All <- melt(CladeI_hullResults_All, id=c("V1","V4"))

# create a grouping variable for plotting on the x-axis this assigns a number to each sample, that is repeated for each hull assignment test creates a specimen-associated variable that can be used to group assignment tests of the same specimen for all species n-cubes
melted_CladeI_hullResults_All$uniqueID <- rep(rep(1:nrow(CladeI_hullResults_All)),length(unique(melted_CladeI_hullResults_All$variable)))

# create a grouping variable for rapidly visualizing (by color) the results of assignment test if the sample IS not assigned to a hull, the value of the "value" column will be FALSE, and we will assign it the white color if a sample IS assigned to a hull, the value will be true, and we will assign a color value here we use black and white, but see below to plot in color
melted_CladeI_hullResults_All$col <- ifelse( melted_CladeI_hullResults_All$value == 'FALSE', "no assignment", "assignment") 


# sometimes, it is the case that no specimens are assigned to any n-cube, in which case no color will be present on the plot however, we want to be able to illustrate (in the legend) the possible results of the assignment test to do so, we need to inform ggplot of the possible levels of the 'col' column, so that it can plot both assigned and unassigned possibilities in the legend keys, regardless of outcome of assignment
# explicitly provide the assignment possibilities for the n-cube overlap proportions and the matching-prediction tests
melted_CladeI_hullResults_All$col <- factor(melted_CladeI_hullResults_All$col, levels = c("assignment", "no assignment"))
colnames(melted_overlapProp_All) <- c("X1","X2","value")

# set the color for the assignment
pal <- c('#636363')
names(pal) <- c('assignment')

# plot the overlap proportions of the n-cubes
CladeI_propOverlap_All <- ggplot(melted_overlapProp_All, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) +
  geom_text(aes(label=round(value,digits=1)),color='darkgrey') +
  scale_fill_gradient(low='white', high="black") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=6),
        axis.text.x = element_text(angle=-45, hjust=0.1)) +
  theme(legend.text = element_text(size=4),
        legend.margin = margin(2,2,6,2),
        legend.title = element_text(size=6),
        legend.background = element_rect(color="#cccccc")) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 2))


# plot the results of the matching-prediction tests
# for aesthetics:
# the following is useful for drawing line segments along the bottom number ranges 
# indicate the uniqueID's associated with the taxonomic identification of each sample these numbers correspond to row numbers in each test dataframe
# use: CladeI_hullResults_All$V4 to retrieve row numbers 
# 1:11 corresponds to specimens identified as E. micrantha
# 12:33 corresponds to specimens identified as E. millegrana
# you will need to fiddle with margins, coord_cartesian (the y-limits), annotations (text and segments) and positioning of text 

CladeI_inHullAssignments_All <- ggplot(melted_CladeI_hullResults_All, 
                                       aes(x=uniqueID, y=variable, fill=col)) +
  geom_tile(aes(width=0.99, height=0.99), 
            color='black') +
  scale_fill_manual(values = c("assignment" = '#636363', "no assignment" = "white"),
                    name=NULL, 
                    labels=c('yes','no'),
                    drop = FALSE) +
  scale_y_discrete(labels=c("In E. micrantha",
                            "In E. millegrana")) +
  scale_x_discrete(limits=factor(1:33)) + # limits equal to number of specimens
  theme_bw() +
  theme(panel.ontop = TRUE,
        panel.background = element_rect(fill='transparent'),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'top',
        text = element_text(size=6),
        plot.margin = margin(1,20,20,10),
        legend.key.height = unit(0.25,'cm'),
        legend.key.width = unit(0.25,'cm')) +
  coord_cartesian(ylim=c(0.5,2),clip = 'off') +
  annotate(geom='text',x=(11/2), y=0.2,label='As E. micrantha',size=2, angle=320, hjust=0) +
  annotate(geom='text',x=(12 + (33-11)/2), y=0.2,label='As E. millegrana',size=2, angle=320, hjust=0) +
  annotate(geom='segment',x=1, xend=11, y=0.3, yend=0.3, color='black', size=0.75) +
  annotate(geom='segment',x=12, xend=33, y=0.3, yend=0.3, color='black', size=0.75) 