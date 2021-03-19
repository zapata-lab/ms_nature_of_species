########################################################################
# Model-based species discovery and delimitation using phenotypic data #
########################################################################

# This is a script developed and written by members of the Zapata Lab

---------------------------------------------------------------------------------

# This script models the number of species and the assignment of specimens to species using phenotypic data without prior information about taxonomy by applying Gaussian finite mixture modeling (GFMM). For details on the analysis and interpretation of results, please read the Methods and Results in our manuscript

# There are four steps to the analysis: 
# First, we rotate the original, log-transformed data matrix into orthogonal axes using robust PCA using package pcaPP
# Next, we reduce the dimensionality of the orthogonal axes to only those that optimize shape, orientation, and the number of clusters (i.e., phenotypic-based species) using packages clustvarsel and mclust
# Then, we identify the best Gaussian Mixture Model (the naive model) using a Bayesian information criterion (BIC) and an integrated complete-data likelihood (ICL) approach using mclust
# Finally, we assess support for alternative models where specimens are assigned to groups a priori that correspond to taxonomy (taxonomy model) and phenogroups identified independent of taxonomy (taxonomy unaware model) using mclust

---------------------------------------------------------------------------------

# Preliminaries
--------------------------------------------------------------------
library(tidyverse)
library(mclust)
library(clustvarsel)
library(parallel)
library(pcaPP)
library(patchwork)

# read data for morphologicalDelimitation
indata = read_csv("../../data/morphometrics_mspub.csv")
CladeI = indata %>%
  dplyr::filter(CladeID == "I", PETLEN != 0) # filter for CladeI and keep only specimens with measurements for both leaves and flowers. Change for other clades, accordingly
  
# filter traits as needed
CladeI_log = CladeI %>%
  dplyr::select( LAMLEN,
          LAMWID,
          PEDLEN,
          OVALEN,
          # CALYTUBLEN, # Taxonomic species in Clade I do not have this trait (calyx tubes)
          CALLOBLEN,
          PETLEN,
          PETWID,
          FILLEN,
          STYLEN ) %>% 
  mutate_all( list( log = log ) ) # log-transform traits
   
CladeI_all = inner_join(CladeI, CladeI_log)

data_analysis = CladeI_all %>%
  dplyr::select(LAMLEN_log,
         LAMWID_log,
         PEDLEN_log,
         OVALEN_log,
         #CALYTUBLEN_log,
         CALLOBLEN_log,
         PETLEN_log,
         PETWID_log,
         FILLEN_log,
         STYLEN_log) 



#
# 1. Robust PCA to rotate the original data into orthogonal axes
--------------------------------------------------------------------
#

# perform robust PCA on the nine, log-transformed variables using function PCAgrid
PCAgrid_qn <- PCAgrid(data_analysis, k=9, maxiter = 100, method='qn', scale='mad', center='median')


#
# 2. reduce dimensionality of the orthogonal axes
--------------------------------------------------------------------
#

# use the scores of the data points on the principal components calculated above to perform clustering and dimensionality reduction in clustvarsel. 
# Make sure to set hcUse=VARS (no variable transformation to initiate clustering.
# Use forward variable selection procedure 
# Model all possible cluster shapes

mclust.options(hcUse = "VARS")
PCAgrid_qn_clustvarsel <- clustvarsel(PCAgrid_qn$scores,
                                      set.seed(4567),
                                      G = 1:3, # set the limit for the number of clusters modeled
                                      search = "greedy",
                                      direction = "forward",
                                      emModels2 = mclust.options("emModelNames"),
                                      samp = FALSE,
                                      hcModel = "VVV",
                                      fit = TRUE)


PCAgrid_qn_clustvarsel
# selected variables include components 7, 5, 1 and 3

#
# 3. Identify the best Gaussian mixture model using BIC and ICL
--------------------------------------------------------------------
#

# as part of the variable selection process, model-based clustering is performed and the BIC score of the models is used to help decide whether to retain or omit a variable. Therefore, the model associated with the best set of variables (and its associated BIC and ICL scores) can be found in the clustvarsel object

PCAgrid_qn_clustvarsel$model$bic
PCAgrid_qn_clustvarsel$model$icl


# but, if you want to run mclust independently (which would allow you to explore the results of fitting other types of models), followed by an independent ICL analysis, do the following:

# model choice using BIC criterion in mclust
cladeI_mclustBIC <- mclustBIC(PCAgrid_qn$scores[,PCAgrid_qn_clustvarsel$subset],
                              set.seed(4567),
                              G = 1:3, # set the limit for the number of clusters modeled
                              search = "greedy",
                              direction = "forward",
                              modelNames = mclust.options("emModelNames"),
                              hcModel = "VVV",
                              fit = TRUE)

cladeI_mclustBIC


# model choice using ICL criterion (as opposed to BIC) in mclust
cladeI_mclustICL <- mclustICL(PCAgrid_qn$scores[,PCAgrid_qn_clustvarsel$subset],
                              set.seed(4567),
                              G = 1:3,
                              search = "greedy",
                              direction = "forward",
                              modelNames = mclust.options("emModelNames"),
                              hcModel = "VVV",
                              fit = TRUE)

cladeI_mclustICL


#
# 4a. Assess support for alternative models 
--------------------------------------------------------------------
#
# Perform a likelihood ratio test to assess whether models with a nested number of phenogroups of equal shape and orientation as the best fit model were appropriate, and approximate the significance with bootstrap for the LRT statistic
# Here you'll need to update the modelName (currently set at "EEI") -- use PCAgrid_qn_clustvarsel$model$modelName to retrieve best model name
# and set the value for G one unit above the number of groups modeled in mclust
cladeI_bootstrapLRT <- mclustBootstrapLRT(PCAgrid_qn$scores[,PCAgrid_qn_clustvarsel$subset], 
                                          modelName = "EEI", 
                                          maxG = 4, 
                                          nboot = 2000)
cladeI_bootstrapLRT

# to plot the test statistic against the bootstraps use the following, changing the value for G to explore distributions that modeled different components
plot(cladeI_bootstrapLRT, G=1)
plot(cladeI_bootstrapLRT, G=2)
plot(cladeI_bootstrapLRT, G=3)
plot(cladeI_bootstrapLRT, G=4)


# save the classification and bind it to the dataset for plotting

CladeI <- CladeI %>% add_column(PCAgrid_qn_clustvarsel$model$classification)

# the position of individuals in morphospace are sometimes best represented using another function that maximizes the differences between points (this is ONLY for visualization); this is a dimension reduction for the features that capture the groups modeled in our winning model

winner_model_2dr <- MclustDR(PCAgrid_qn_clustvarsel$model)
# the 'dir' value stores the estimated directions of each point, from the dimension reduction
# we'll store this for use in plotting later
morphDR <- as.data.frame(winner_model_2dr$dir)


#
# 4b. Examine alternative models that correspond to taxonomy and taxonomy unaware groups 
--------------------------------------------------------------------
#

# in our morphometrics_mspub.csv file, the variable coding current taxonmy is the 'Taxon' column and the variable coding our taxonomy unaware groups is the 'MorphogroupID' column

# use a discriminant analysis (DA) based on Gaussian finite mixture modeling where specimens are assigned to groups a priori
# DA of taxonomic groups
taxonomy =
  as_tibble(PCAgrid_qn$scores) %>%
  dplyr::select( names(PCAgrid_qn_clustvarsel$subset) ) %>% 
  MclustDA( class=as.vector(CladeI_all$Taxon), G = 1 )

# DA of taxonomy unaware groups
morphogroups =
  as_tibble(PCAgrid_qn$scores) %>%
  dplyr::select( names(PCAgrid_qn_clustvarsel$subset) ) %>%
  MclustDA( class=as.vector(CladeI_all$MorphogroupID), G=1)


#
# Visualizations
--------------------------------------------------------------------
#

# bind the position of each specimen in the winning naive model (note: this is using the positions from the MclustDR analysis performed on the winning mclustBIC model) to its classification, and also bind the probability of its classification to each of the identified clusters, as well as any specimen-specific information that you might want to plot --- elevation, latitude, longitude, etc.
morphDR <- cbind(morphDR, 
                 PCAgrid_qn_clustvarsel$model$classification, 
                 PCAgrid_qn_clustvarsel$model$z,
                 CladeI_all[,c(4,6,7,8)]) # identified which columns from cladeI_all corresponded with SpecimenID(4), elevation(6), latitude(7) and longitude(8)

# set some helpful column names for the tibble
colnames(morphDR) <- c("MclustDR_Dir1", 
                       "BestModel_Classification", 
                       "BestModel_ProbAssign1","BestModel_ProbAssign2",
                       "SpecimenID","Elevation","Latitude","Longitude")



# now, we'll build a dataframe to help us plot the BIC rankings of all the models we've considered: the naive, taxonomy, and taxonomy-unaware models
# retrieve the BIC score for each number of naive groups modeled; identify the 'best' BIC score and its corresponding number of groups
BIC.Best.Model.Per.G <- apply(PCAgrid_qn_clustvarsel$model$BIC, 1, max, na.rm=T)
max.BIC <- max(BIC.Best.Model.Per.G)
# store the ranking of the naive, taxonomy, and taxonomy unaware models in vectors of the same length as the number of groups modeled
numMorphGroups <- c(1,2,3) #total number of morphogroups modeled
taxonBIC <- c(NA,taxonomy$bic,NA) # taxonomy model BIC, stored in the vector position that corresponds to the number of currently recognized taxa
morphoGroupsBIC <- c(morphogroups$bic,NA,NA) # taxonomy unaware model BIC, stored in the vector position that corresponds to the number of morphogroups

# convert vectors into dataframe (for plotting) and calculate delta-BIC
CladeI_modelsRanked <- as.data.frame(cbind(numMorphGroups, 
                                           BIC.Best.Model.Per.G, 
                                           taxonBIC, 
                                           morphoGroupsBIC))
CladeI_modelsRanked$deltaBIC_naive <- unlist(lapply(CladeI_modelsRanked$BIC.Best.Model.Per.G, function(x) max.BIC-(x)))
CladeI_modelsRanked$deltaBIC_taxon <- unlist(lapply(CladeI_modelsRanked$taxonBIC, function(x) max.BIC-(x)))
CladeI_modelsRanked$deltaBIC_morphoGroups <- unlist(lapply(CladeI_modelsRanked$morphoGroupsBIC, function(x) max.BIC-(x)))

# Plot
modelsRanked = ggplot(CladeI_modelsRanked, 
                      aes(x=numMorphGroups, 
                          y=deltaBIC_naive)) +  
  geom_point(aes(x=numMorphGroups, 
                 y=deltaBIC_naive, 
                 shape="A_naive_models"), 
             cex=3, 
             col="black", 
             bg='white') +
  geom_point(aes(x=numMorphGroups[order(deltaBIC_naive)[1]], 
                 y=deltaBIC_naive[order(deltaBIC_naive)[1]], 
                 shape="B_best_naive_model"), 
             cex=3, 
             col='black', 
             bg = 'black') +
  geom_point(aes(x=numMorphGroups[taxonBIC != 'NA'], 
                 y=max.BIC-taxonBIC[taxonBIC != 'NA'], 
                 shape="D_currentTax"), 
             cex=3, 
             col='black', 
             bg='red') +
  geom_point(aes(x=numMorphGroups[morphoGroupsBIC != 'NA'], 
                 y=max.BIC-morphoGroupsBIC[morphoGroupsBIC != 'NA'], 
                 shape="E_morphoGroups"), 
             cex=3, 
             col='black', 
             bg='red') +
  scale_shape_manual(values=c(21,21,25,24),
                     labels=c("Best model overall",
                              "Second best model overall",
                              "Naive models", 
                              "Current taxonomy",
                              "Taxonomy-unaware morphogroups"), 
                     name=NULL) +
  guides(shape = guide_legend(override.aes = list(shape = c(21,21,25,24),
                                                  bg=c('black','white','red','red'),
                                                  color=c('black','black','black','black'),
                                                  size=3))) +
  scale_x_discrete(limits=factor(1:3)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab('Number of pheno-groups') +
  ylab(expression(paste("Empirical support (", Delta, "BIC)", sep="")))



# plot ordination of best naive model
myShapes = c(21,22)
Q2 = ggplot(morphDR) +
  geom_boxplot(aes(x = MclustDR_Dir1,
                   y = factor(1)),
               width = 0.001,
               lwd = 0) +
  geom_jitter(aes(x = MclustDR_Dir1,
                  y = factor(1),
                  shape = as.factor(BestModel_Classification)),
              color = 'black',
              bg = '#cccccc',
              size = 2,
              stroke = 0.25,
              width = 0,
              height = 0.03) +
  scale_shape_manual(values = myShapes,
                     name = 'Classification',
                     labels = c('group1', 'group2')) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  theme(legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,
                                    color = 'black'),
        text = element_text(size=10),
        axis.ticks.y = element_blank()) +
  guides(shape = guide_legend(override.aes = list(size=2)))


Q3 = ggplot(morphDR) +
  geom_point(aes(x = MclustDR_Dir1,
                 y = Latitude,
                 shape = as.factor(BestModel_Classification)),
             color = 'black',
             bg = '#cccccc',
             size = 2,
             stroke = 0.25) +
  scale_shape_manual(values = myShapes,
                     name = 'Classification',
                     labels = c('group1', 'group2')) +
  xlab("Dimension 1") +
  ylab("Latitude") +
  theme(legend.key = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,
                                    color='black'),
        text = element_text(size=10)) +
  guides(shape = guide_legend(override.aes = list(size=2)))

Q4 = ggplot(morphDR) +
  geom_point(aes(x = MclustDR_Dir1,
                 y = Elevation,
                 shape = as.factor(BestModel_Classification)),
             color = 'black',
             bg = '#cccccc',
             size = 2,
             stroke = 0.25) +
  scale_shape_manual(values = myShapes,
                     name = 'Classification',
                     labels = c('group1', 'group2')) +
  xlab("Dimension 1") +
  ylab("Elevation") +
  theme(legend.key = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,
                                    color = 'black'),
        text = element_text(size=10)) +
  guides(shape = guide_legend(override.aes = list(size=2)))


# plot everything together - requires patchwork library
(modelsRanked + plot_spacer() + plot_layout(widths=c(3,1))) /
  (Q2 + Q3 + Q4 + plot_layout(guides='collect'))
