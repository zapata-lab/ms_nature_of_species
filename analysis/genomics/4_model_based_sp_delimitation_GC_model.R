########################################################################
# GC Model-based species discovery and delimitation using genomic data #
########################################################################

# This is a script developed and written by members of the Zapata Lab

---------------------------------------------------------------------------------

# This script models the number of species and the assignment of specimens to species using 
# genomic data without prior information about taxonomy by applying Gaussian finite mixture 
# modeling (GFMM). For details on the analysis and interpretation of results, please read 
# the Methods and Results in our manuscript

# There are three steps to the analysis: 
# First, using a data matrix with a single SNP per locus, we estimate shared 
# allele distances using alleleinit from prabclus.

# Next, we apply non-metric multidimensional scaling (NMDS) to reduce dimensionality.

# Finally, Gaussian finite mixture modeling is employed to search for genomic clusters. 
# This is performed by mclust, through the prabclus wrapper function prabclust


---------------------------------------------------------------------------------

# Preliminaries
--------------------------------------------------------------------

library(tidyverse)
library(prabclus)
library(RColorBrewer)
library(patchwork)


# This approach takes as input one data file:
# A STRUCTURE-formatted file separated by commas. This comma-delimited file is ONLY provided
# for CladeI. You will need to generate it for the other clades.


######################################################################
######################################################################

# import SNP data (STRUCTURE file format, ".str") and calculate distance matrix.
# Change input file for other clades accordingly. 

# NOTE: THIS COMMA-DELIMITED FILE IS NOT PROVIDED FOR ALL CLADES, ONLY CLADE_I
# CREATE A FILE FOR EACH CLADE
infile.tmp <- read.delim("../..data/CladeI/CladeI_filtered_str_prabclusFormat.csv", sep = ',', header=FALSE)
ncol(infile.tmp) # records the number of SNPs in your file - put this number in the in the next line

# convert the STRUCTURE-formatted file into format needed for prabclus --- a matrix in prabclus format
infile <- alleleconvert(strmatrix=infile.tmp[,c(1,8:5624)], # total number of SNPs goes here
                        format.in="structure",
                        firstcolname = TRUE, 
                        orig.nachar = "-9")

# calculate allele distances
infile_init <- alleleinit(allelematrix=infile, distance="alleledist")


# read in collection information: list of each individual that we have locality for (likely more than the individuals in the SNP dataset)
# Change input file for other clades accordingly
cladeI.collinfo <- read.csv(file="../../data/CladeI/CladeI_geneticSamples.csv", header=TRUE)

# preps a dataframe to bind clustering results
cladeI.master <- cladeI.collinfo

# get an idea about number of dimensions needed to explain variation
# per Hausdorf and Hennig, results are reliable when we can keep stress as low as possible (at least below 20%, ideally lower, if possible)
# at the same time we want to keep the number of dimensions as low as possible
# our strategy was to keep as few dimensions as possible and stay below 15-20%
cladeI.stressvalues <- stressvals(infile_init, mdsdim=1:5)
cladeI.stressvalues$MDSstress
plot(cladeI.stressvalues$MDSstress)

# prepare mclust for clustering
# confirm that mclust (accessed by prabclust) is using the correct variable ('VARS') for hcUse
mclust.options()
mclust.options(hcUse="VARS")



#########################################################
###################  PRABCLUST  #########################
#########################################################
#### RUN the series of prabclust analyses

# 2 dimensions
# 100 permutations
infile.pclust.dim2.n0 <- prabclust(infile_init, mdsmethod = "kruskal", nnk=0, nclus=1:3, modelid="all", permutations=100, mdsdim=2)

# 3 dimensions
# 100 permutations
infile.pclust.dim3.n0 <- prabclust(infile_init, mdsmethod = "kruskal", nnk=0, nclus=1:3, modelid="all", permutations=100, mdsdim=3)


#########################################################
############  SAVE THE RESULTS  #########################
#########################################################


# we want the classification and associated probability of assignment of each individual to a genetic group
# we also want the NMDS points from dimensions 1 and 2, for each analysis performed above

# the classifications are stored under $clustering
# the probability of assignment is stored under $clustsummary$z --- one column for each of the groups found
# the NMDS scores are stored under $points




##
# dimensions = 2, noise = 0
##
infile.pclust.dim2.n0$clustsummary 
# record classification of each individual - what group do they belong to
cladeI.master <- cladeI.master %>% mutate(class.dim2.n0=infile.pclust.dim2.n0$clustering[match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim2.n0$clustering)))])
# record probability of assignment for each individual for each group
cladeI.master <- cladeI.master %>% mutate(probass.dim2.n0.g1=infile.pclust.dim2.n0$clustsummary$z[,1][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim2.n0$clustering)))])
cladeI.master <- cladeI.master %>% mutate(probass.dim2.n0.g2=infile.pclust.dim2.n0$clustsummary$z[,2][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim2.n0$clustering)))])
cladeI.master <- cladeI.master %>% mutate(probass.dim2.n0.g3=infile.pclust.dim2.n0$clustsummary$z[,3][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim2.n0$clustering)))])
# record position in modeled dimensions for each individual
cladeI.master <- cladeI.master %>% mutate(NMDS1.dim2.n0=infile.pclust.dim2.n0$points[,1][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim2.n0$points)))])
cladeI.master <- cladeI.master %>% mutate(NMDS2.dim2.n0=infile.pclust.dim2.n0$points[,2][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim2.n0$points)))])
infile.pclust.dim2.n0$bicsummary 

##
# dimensions = 3, noise = 0
##
infile.pclust.dim3.n0$clustsummary 
# record classification of each individual - what group do they belong to
cladeI.master <- cladeI.master %>% mutate(class.dim3.n0=infile.pclust.dim3.n0$clustering[match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim3.n0$clustering)))])
# record probability of assignment for each individual for each group
cladeI.master <- cladeI.master %>% mutate(probass.dim3.n0.g1=infile.pclust.dim3.n0$clustsummary$z[,1][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim3.n0$clustering)))])
cladeI.master <- cladeI.master %>% mutate(probass.dim3.n0.g2=infile.pclust.dim3.n0$clustsummary$z[,2][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim3.n0$clustering)))])
cladeI.master <- cladeI.master %>% mutate(probass.dim3.n0.g3=infile.pclust.dim3.n0$clustsummary$z[,3][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim3.n0$clustering)))])
# record position in modeled dimensions for each individual
cladeI.master <- cladeI.master %>% mutate(NMDS1.dim3.n0=infile.pclust.dim3.n0$points[,1][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim3.n0$points)))])
cladeI.master <- cladeI.master %>% mutate(NMDS2.dim3.n0=infile.pclust.dim3.n0$points[,2][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim3.n0$points)))])
cladeI.master <- cladeI.master %>% mutate(NMDS3.dim3.n0=infile.pclust.dim3.n0$points[,2][match(rownames(cladeI.master), rownames(as.data.frame(infile.pclust.dim3.n0$points)))])
infile.pclust.dim3.n0$bicsummary 

cladeI.master
dim(cladeI.master)
colnames(cladeI.master)

# cladeI.master now holds all results of the analysis


#########################################################
############    VISUALIZATION   #########################
#########################################################


mycols=brewer.pal(3, 'Greys')
p1 <- ggplot(cladeI.master, aes(x=NMDS1.dim2.n0, y=NMDS2.dim2.n0, fill=as.factor(class.dim2.n0))) +
  geom_jitter(pch=21, cex=5, width=0.008, height=0.008) +
  scale_fill_manual(values=mycols, labels=c('group1','group2','group3'), name='Classification') +
  theme_bw() +
  theme(text = element_text(size=15)) +
  xlab('Dimension 1') + ylab('Dimension 2')

p2 <- ggplot(cladeI.master, aes(x=NMDS1.dim2.n0, y=Latitude, fill=as.factor(class.dim2.n0))) +
  geom_jitter(pch=21, cex=5, width=0.008, height=0.008) +
  scale_fill_manual(values=mycols, labels=c('group1','group2','group3'), name='Classification') +
  theme_bw() +
  theme(text = element_text(size=15)) +
  xlab('Dimension 1') + ylab('Latitude')

p3 <- ggplot(cladeI.master, aes(x=NMDS1.dim2.n0, y=Elevation, fill=as.factor(class.dim2.n0))) +
  geom_jitter(pch=21, cex=5, width=0.008, height=0.008) +
  scale_fill_manual(values=mycols, labels=c('group1','group2','group3'), name='Classification') +
  theme_bw() +
  theme(text = element_text(size=15)) +
  xlab('Dimension 1') + ylab('Elevation')

p1+p2+p3+plot_layout(guides='collect')



