### Install Required Packages ###
#install.packages(c("rentrez", "cluster", "fpc", "vegan"))
#BiocManager::install(c("Biostrings", "DECIPHER"))

### Load Packages ###
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(cluster)
library(fpc)
library(vegan)

# Data Collection
# These are a list of the species we will be using

source("functions.R")

candida6_species_term <- paste0(
  '("Candida auris"[ORGN] OR ',
  '"Candida glabrata"[ORGN] OR ',
  '"Candida krusei"[ORGN] OR ',
  '"Candida albicans"[ORGN] OR ',
  '"Candida tropicalis"[ORGN] OR ',
  '"Candida parapsilosis"[ORGN])'
)

##### ITS Gene #####

candida6_ITS_term <- paste0(
  candida6_species_term,
  ' AND "internal transcribed spacer"[All Fields] AND 400:800[SLEN]'
)

candida6_ITS_term
geneITS = "ITS"

search_fetch(candida6_ITS_term, geneITS)

library(Biostrings)

candida6_ITS_seqs <- readDNAStringSet("../data/raw/candida6_ITS_500.fasta")

dfITS <- gene_dataframe(candida6_ITS_seqs)

table(dfITS$Species_Name)

createhistograms(dfITS, geneITS)

dfITS_subset <- subset_species(dfITS)

write.csv(dfITS_subset,
          file = "../data/clean/candida6_ITS_minperSpecies.csv",
          row.names = FALSE)

#### ERG111 Gene #####

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

## Part 2.1: Pairwise Distance Calculated for Both Genes

ITS_dist   <- pairwise_distance(dfITS_subset,  geneITS)
ERG11_dist <- pairwise_distance(dfERG11_subset, geneERG11)

## Part 3.1: Hierarchical Clustering for ITS and ERG11

h_clustering(ITS_dist, geneITS)
h_clustering(ERG11_dist, geneERG11)


## Part 7.1 NMDS for alignment-based clustering 


aligncluster(dfITS_subset, ITS_dist, geneITS, "topright")
aligncluster(dfERG11_subset, ERG11_dist, geneERG11, "right")

## Part 4.1: Feature Based Clustering ITS and ERG11

# Feature Extraction + PCA

gene_features(dfITS_subset, geneITS, "topright")
gene_features(dfERG11_subset, geneERG11, "topright")


## Part 5.1 Silhouette Index for Alignment-Based Clustering

dindex_ITS <- silhouetteplot(dfITS_subset, ITS_dist, geneITS)
dindex_ERG11 <- silhouetteplot(dfERG11_subset, ERG11_dist, geneERG11)

## Part 6.1 Dunn Index

dindex_ITS$dunn
dindex_ERG11$dunn








