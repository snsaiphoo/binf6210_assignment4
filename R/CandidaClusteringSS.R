#### Install Required Packages ####
# These lines of code are meant to be un-commented, in case the user needs to 
# install the packages 
#install.packages(c("rentrez", "cluster", "fpc", "vegan"))
#BiocManager::install(c("Biostrings", "DECIPHER"))

#### Load Packages ####
# Initializing the libraries necessary to run the clustering script
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(cluster)
library(fpc)
library(vegan)

#### Data Collection ####

# This is calling the file "functions.R" that is located in the same folder as 
# this R file. This file houses all the functions seen in this script

source("functions.R")

# These are a list of the 6 Candida species we will be investigating,
# 3 are of high resistance concern -> Candida auris, Candida glabrata, Candida krusei
# 3 are of low resistance concern -> Candida albicans, Candida parapsilosis, Candida tropicalis
# Here we are generating the starting part of the search term, this will be used 
# for both ITS and ERG11 Genes

candida6_species_term <- paste0(
  '("Candida auris"[ORGN] OR ',
  '"Candida glabrata"[ORGN] OR ',
  '"Candida krusei"[ORGN] OR ',
  '"Candida albicans"[ORGN] OR ',
  '"Candida tropicalis"[ORGN] OR ',
  '"Candida parapsilosis"[ORGN])'
)

##### ITS Gene #####
# Searching for the ITS region, this is commonly referred to as the fungal barcode gene.
# The ITS sequence is typically 400–800 bp, so we filter by length.
# This is building the NCBI Search Term

candida6_ITS_term <- paste0(
  candida6_species_term,
  ' AND "internal transcribed spacer"[All Fields] AND 400:800[SLEN]'
)

# This is a verification step, that the search term looks as we want it.

candida6_ITS_term

# This is a variable we will use as an argument for the functions in order to 
# label graphs and files accordingly.

geneITS = "ITS"

# Calling the rentrez search and fetch function to get the top 500 results.
# Commenting about the functions can be seen in the functions.R file

search_fetch(candida6_ITS_term, geneITS)

# The previous step, creates the FASTA file of the 500 ITS matches and stores them
# in a file in our raw folder. We then call on the readDNAStringSet function to 
# change the FASTA file into a format we can work with and analyze.

candida6_ITS_seqs <- readDNAStringSet("../data/raw/candida6_ITS_500.fasta")

# This function takes the output from the previous step and turns it into a dataframe
# More information on this function can be seen in the functions.R file 
# This returns a dataframe of 500 ITS samples from the 6 Candida species

dfITS <- gene_dataframe(candida6_ITS_seqs)

# This is used to ensure there are only 6 Candida species names 

table(dfITS$Species_Name)

# This step is creating exploratory histograms to identify how our data is distributed between
# the species of high and low concern. These plots are saved in the figs folder, but won't be 
# used for analysis.

createhistograms(dfITS, geneITS)

# This next step creates a subset of the ITS dataframe so that it can be used for clustering.
# This function determines the smallest number of samples from one species, and re-samples
# the others to match the minimum number. This is to ensure there is minimal bias when clustering

dfITS_subset <- subset_species(dfITS)

# The subsetted dataframe is then written to a csv file for later use. 

write.csv(dfITS_subset,
          file = "../data/clean/candida6_ITS_minperSpecies.csv",
          row.names = FALSE)

##### ERG111 Gene #####
# The steps are reproduced from the ITS gene, 
# ERG11 is a drug-target gene involved in ergosterol synthesis. Mutations in ERG11
# contribute to azole resistance, making it important for comparing species-level
# clustering patterns.

erg11_gene_term <- paste0(
  '(',
  'ERG11[Gene] OR ',
  '"lanosterol 14-alpha demethylase"[Title] OR ',
  '"lanosterol 14 alpha demethylase"[Title] OR ',
  'CYP51[Gene]',
  ')'
)

length_term <- '1000:2500[SLEN]'

candida6_ERG11_term <- paste(
  candida6_species_term,
  "AND",
  erg11_gene_term,
  "AND",
  length_term
)

candida6_ERG11_term
geneERG11 = "ERG11"

search_fetch(candida6_ERG11_term, geneERG11)

candida6_ERG11_seqs <- readDNAStringSet("../data/raw/candida6_ERG11_500.fasta")

dfERG11 <- gene_dataframe(candida6_ERG11_seqs)
table(dfERG11$Species_Name)
createhistograms(dfERG11, geneERG11)
dfERG11_subset <- subset_species(dfERG11)

write.csv(dfERG11_subset,
          file = "../data/clean/candida6_ERG11_minperSpecies.csv",
          row.names = FALSE)

#### Pairwise Distances ####
# Computed alignment-based pairwise distance matrices for ITS and ERG11.
# These distances quantify how different each sequence is from every other sequence
# and are used as input for hierarchical clustering, NMDS, and silhouette analysis.

ITS_dist   <- pairwise_distance(dfITS_subset,  geneITS)
ERG11_dist <- pairwise_distance(dfERG11_subset, geneERG11)

#### Hierarchical Clustering ####
# Perform hierarchical clustering on the ITS and ERG11 distance matrices.
# These dendrograms visualize how sequences group based on alignment-based
# genetic similarity across the six high & low resistance Candida species.

h_clustering(ITS_dist, geneITS)
h_clustering(ERG11_dist, geneERG11)


#### NMDS - Alignment-Based Clustering ####
# Non-metric multidimensional scaling 
# alignment-based distance matrices. NMDS projects the high-dimensional
# distances into 2D space, allowing visualization of species-level clustering
# patterns across high and low resistance Candida.

aligncluster(dfITS_subset, ITS_dist, geneITS, "topright")
aligncluster(dfERG11_subset, ERG11_dist, geneERG11, "right")

#### Feature Based Clustering  + PCA #####
# basic characteristics like length, GC content, nucleotide frequencies are calculated
# for each sequence, and these features are used to generate both a 
# hierarchical clustering dendrogram and a PCA ordination plot. The ordination plot
# will be used in the analysis. The dendrogram is for exploratory analysis.
# This provides an alternative view of species separation independent of sequence alignment.

gene_features(dfITS_subset, geneITS, "topright")
gene_features(dfERG11_subset, geneERG11, "topright")


#### Silhouette Index for Alignment-Based Clustering ####
# Silhouette values quantify how well each sequence fits within its
# species cluster relative to other species, providing a numerical measure of
# clustering quality for the alignment-based distance matrices.

dindex_ITS <- silhouetteplot(dfITS_subset, ITS_dist, geneITS)
dindex_ERG11 <- silhouetteplot(dfERG11_subset, ERG11_dist, geneERG11)

#### Dunn Index ####
# The Dunn index evaluates the ratio
# between inter-cluster separation and intra-cluster compactness; higher values
# indicate stronger, more well-defined species clusters.

dindex_ITS$dunn
dindex_ERG11$dunn








