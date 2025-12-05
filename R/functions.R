# Write the searched sequences to a FASTA file using base R. 
# If seq_vector has names, they are used as FASTA headers.

write_fasta <- function(seq_vector, filepath) {
  write(seq_vector, file = filepath, sep = "\n")
}

# The search_fetch function takes the term and the gene name,
# uses the rentrez library to search for 500 matches of the species,
# and saves them into a FASTA file.

search_fetch <- function(term, gene_name) {
  
  library(rentrez)
  
  search_res <- entrez_search(
    db = "nuccore",
    term = term,
    retmax = 500,
    use_history = TRUE
  )
  
  fasta_res <- entrez_fetch(
    db = "nuccore",
    web_history = search_res$web_history,
    rettype = "fasta",
    retmax = 500
  )
  
  outfile <- paste0("../data/raw/candida6_", gene_name, "_500.fasta")
  write_fasta(fasta_res, outfile)
  
}

# The gene_dataframe function takes a named sequence vector and creates a
# clean data frame by extracting species names from the sequence titles,
# standardizing species labels, and returning a formatted table for analysis.

gene_dataframe <- function(geneframe) {
  
  # Make initial data frame with consistent column names
  dfgene <- data.frame(
    Title    = names(geneframe),
    Sequence = as.character(geneframe),
    stringsAsFactors = FALSE
  )
  
  # Split titles by spaces
  title_split <- strsplit(dfgene$Title, " ")
  
  # Extract species name (2nd and 3rd word)
  species_name <- sapply(
    title_split,
    function(x) {
      if (length(x) >= 3) {
        paste(x[2], x[3])
      } else {
        NA_character_
      }
    }
  )
  
  dfgene$Species_Name <- species_name
  
  # Reorder columns
  dfgene <- dfgene[, c("Title", "Species_Name", "Sequence")]
  
  print(table(dfgene$Species_Name))
  
  # Fix species naming
  dfgene$Species_Name[dfgene$Species_Name == "Candidozyma auris"]    <- "Candida auris"
  dfgene$Species_Name[dfgene$Species_Name == "Pichia kudriavzevii"]  <- "Candida krusei"
  dfgene$Species_Name[dfgene$Species_Name == "Nakaseomyces glabratus"] <- "Candida glabrata"
  
  return(dfgene)
}

# The createhistograms function groups species by resistance concern,
# calculates sequence lengths, and saves two histogram plots (high vs. low
# resistance concern) for the given gene.

createhistograms <- function(geneframe, gene_name) {
  
  geneframe$Resistance_group <- NA
  
  high_res <- c("Candida auris", "Candida glabrata", "Candida krusei")
  low_res <- c("Candida albicans", "Candida parapsilosis", "Candida tropicalis")
  
  geneframe$Resistance_group[geneframe$Species_Name %in% high_res] <- "High resistance concern"
  geneframe$Resistance_group[geneframe$Species_Name %in% low_res]  <- "Lower resistance concern"
  
  geneframe$Seq_Length <- nchar(geneframe$Sequence)
  
  outfile1 <- paste0("../figs/", gene_name, "_high_res_hist.png")
  
  png(outfile1, width = 800, height = 600)
  hist(
    geneframe$Seq_Length[geneframe$Resistance_group == "High resistance concern"],
    main = paste0(gene_name, ": Sequence Lengths (High Resistance Concern)"),
    xlab = "Sequence Length (bp)",
    col = "tomato",
    breaks = 20
  )
  dev.off()
  
  outfile2 <- paste0("../figs/", gene_name, "_low_res_hist.png")
  
  png(outfile2, width = 800, height = 600)
  hist(
    geneframe$Seq_Length[geneframe$Resistance_group == "Lower resistance concern"],
    main = paste0(gene_name,": ITS Sequence Lengths (Lower Resistance Concern Species)"),
    xlab = "Sequence Length (bp)",
    col = "skyblue",
    breaks = 20
  )
  dev.off()
  
}

# The subset_species function takes an unequal dataset and randomly selects
# an equal number of sequences per species, returning a balanced subset.

subset_species <- function(geneframe) {
  
  table(geneframe$Species_Name)
  min_n <- min(table(geneframe$Species_Name))  
  idx_by_species <- split(seq_len(nrow(geneframe)), geneframe$Species_Name)
  
  set.seed(123)  
  idx_sampled <- unlist(
    lapply(idx_by_species, function(idx) sample(idx, min_n))
  )
  
  subset_frame <- geneframe[idx_sampled, ]
  table(subset_frame$Species_Name)
  sum(is.na(subset_frame$Sequence))
  
  return(subset_frame)
}

# The pairwise_distance function aligns the sequences, computes a pairwise
# distance matrix, saves it as a CSV file, and returns the distance object.

pairwise_distance <- function(genesubset, gene_name) {
  library(DECIPHER)
  library(Biostrings)
  
  gene <- DNAStringSet(genesubset$Sequence)
  names(gene) <- paste(genesubset$Species_Name,
                       seq_len(nrow(genesubset)),
                       sep = "_")
  
  gene_aln <- AlignSeqs(gene, processors = NULL)
  
  gene_dist <- DistanceMatrix(gene_aln, type = "dist")
  
  gene_dist_mat <- as.matrix(gene_dist)
  write.csv(
    gene_dist_mat,
    file = paste0("../data/clean/", gene_name, "_pairwise_distance_matrix.csv"),
    quote = FALSE
  )
  
  return(gene_dist)
}

# The h_clustering function performs hierarchical clustering on the distance
# matrix and saves the resulting dendrogram as a PNG file.

h_clustering <- function(dist_mat, gene_name) {
  
  hc <- hclust(dist_mat, method = "average")
  
  outfile <- paste0("../figs/", gene_name, "_hierarchal_clustering.png")
  
  png(outfile, width = 800, height = 600)
  
  plot(
    hc,
    main = paste0("Hierarchical Clustering of ", gene_name, " Sequences"),
    xlab = "",
    sub = "",
    cex = 0.6
  )
  dev.off()
}

# The compute_features_biostring function calculates basic sequence features
# (length, GC content, and nucleotide frequencies) using Biostrings.

compute_features_biostring <- function(seq_vec) {
  
  library(Biostrings)
  
  dna <- DNAStringSet(seq_vec)
  
  # Length
  length_vec <- width(dna)
  
  # Nucleotide counts (A, C, G, T)
  nuc_counts <- letterFrequency(dna, letters = c("A","C","G","T"))
  
  # GC content
  GC <- (nuc_counts[, "G"] + nuc_counts[, "C"]) / length_vec
  
  # Convert counts to frequencies
  nuc_freq <- nuc_counts / length_vec
  
  # Combine into a dataframe
  data.frame(
    Seq_Length = length_vec,
    GC_content = GC,
    A_freq = nuc_freq[, "A"],
    C_freq = nuc_freq[, "C"],
    G_freq = nuc_freq[, "G"],
    T_freq = nuc_freq[, "T"]
  )
}

# The gene_features function computes sequence features, scales them,
# performs feature-based clustering and PCA, and saves both

gene_features <- function(genesubset, gene_name, leg) {
  
  features <- compute_features_biostring(genesubset$Sequence)
  
  features$Species <- genesubset$Species_Name
  features$Resistance_group <- genesubset$Resistance_group
  
  head(features)
  
  f_numeric <- features[, c("Seq_Length", "GC_content",
                            "A_freq", "C_freq", "G_freq", "T_freq")]
  
  scaled <- scale(f_numeric)
  
  dist_feat <- dist(scaled, method = "euclidean")
  hc_feat <- hclust(dist_feat, method = "average")
  
  outfile <- paste0("../figs/", gene_name, "_featurebased_clustering.png")
  
  png(outfile, width = 800, height = 600)
  
  plot(
    hc_feat,
    labels = features$Species,
    main = paste0("Feature-based Hierarchical Clustering of ", gene_name, " Sequences"),
    xlab = "",
    sub = "",
    cex = 0.6
  )
  dev.off()
  
  pca <- prcomp(scaled)
  
  outfile <- paste0("../figs/", gene_name, "_PCA_featurebasedclustering.png")
  
  png(outfile, width = 800, height = 600)
  
  plot(
    pca$x[,1], pca$x[,2],
    col = as.factor(features$Species),
    pch = 19,
    main = paste0("PCA of ", gene_name, " Features"),
    xlab = "PC1", ylab = "PC2"
  )
  legend(leg, legend = levels(as.factor(features$Species)),
         col = 1:6, pch = 19, cex = 1.2)
  dev.off()

}

# The silhouetteplot function computes silhouette scores for the clustering
# results, saves the silhouette plot, and returns clustering statistics.

silhouetteplot <- function(genesubset, dist, gene_name) {
  
  library(cluster)
  library(fpc)
  
  # Convert species labels to numeric clusters
  species <- as.factor(genesubset$Species_Name)
  cluster_labels <- as.numeric(species)
  
  # Compute silhouette scores
  sil <- silhouette(cluster_labels, dist)
  
  # Average silhouette width
  mean_sil <- mean(sil[, 3])

  stats <- cluster.stats(
    d = as.matrix(dist),
    clustering = cluster_labels
  )
  
  outfile <- paste0("../figs/", gene_name, "_silhouette.png")
  
  png(outfile, width = 800, height = 600)
  
  # Silhouette plot
  plot(
    sil,
    border = NA,
    main = paste0("Silhouette Plot – ", gene_name, " Alignment-Based Clustering")
  )
  dev.off()
  
  return(stats)
}

# The aligncluster function performs NMDS on the distance matrix and saves
# a plot showing species clustering.

aligncluster <- function(genesubset, dist, gene_name, leg) {
  
  library(vegan)
  nmds <- metaMDS(as.matrix(dist))
  
  outfile <- paste0("../figs/", gene_name, "_NMDS.png")
  png(outfile, width = 800, height = 600)
  
  plot(
    nmds$points,
    col = as.factor(genesubset$Species_Name),
    pch = 19,
    main = paste0("NMDS – ", gene_name, " Distances")
  )
  legend(leg, legend = levels(as.factor(genesubset$Species_Name)),
         col = 1:6, pch = 19, cex = 1.2)
  dev.off()
  
}
