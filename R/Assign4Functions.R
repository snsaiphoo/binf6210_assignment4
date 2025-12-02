# Data Collection
# These are a list of the species we will be using
candida6_species_term <- paste0(
  '("Candida auris"[ORGN] OR ',
  '"Candida glabrata"[ORGN] OR ',
  '"Candida krusei"[ORGN] OR ',
  '"Candida albicans"[ORGN] OR ',
  '"Candida tropicalis"[ORGN] OR ',
  '"Candida parapsilosis"[ORGN])'
)

candida6_ITS_term <- paste0(
  candida6_species_term,
  ' AND "internal transcribed spacer"[All Fields] AND 400:800[SLEN]'
)

candida6_ITS_term
geneITS = "ITS"

write_fasta <- function(seq_vector, filepath) {
  write(seq_vector, file = filepath, sep = "\n")
}

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

search_fetch(candida6_ITS_term, geneITS)

library(Biostrings)

candida6_ITS_seqs <- readDNAStringSet("../data/raw/candida6_ITS_500.fasta")

# function to create the dataframe
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
  
  # Optional: show species counts
  print(table(dfgene$Species_Name))
  
  # Fix species naming
  dfgene$Species_Name[dfgene$Species_Name == "Candidozyma auris"]    <- "Candida auris"
  dfgene$Species_Name[dfgene$Species_Name == "Pichia kudriavzevii"]  <- "Candida krusei"
  dfgene$Species_Name[dfgene$Species_Name == "Nakaseomyces glabratus"] <- "Candida glabrata"
  
  return(dfgene)
}

dfITS <- gene_dataframe(candida6_ITS_seqs)

table(dfITS$Species_Name)

#function for histograms

dfITS$Resistance_group <- NA

high_res <- c("Candida auris", "Candida glabrata", "Candida krusei")
low_res  <- c("Candida albicans", "Candida parapsilosis", "Candida tropicalis")

dfITS$Resistance_group[dfITS$Species_Name %in% high_res] <- "High resistance concern"
dfITS$Resistance_group[dfITS$Species_Name %in% low_res]  <- "Lower resistance concern"

dfITS$Seq_Length <- nchar(dfITS$ITS_Sequence)

hist(
  dfITS$Seq_Length[dfITS$Resistance_group == "High resistance concern"],
  main = "ITS Sequence Lengths: High Resistance Concern Species",
  xlab = "Sequence Length (bp)",
  col = "tomato",
  breaks = 20
)

hist(
  dfITS$Seq_Length[dfITS$Resistance_group == "Lower resistance concern"],
  main = "ITS Sequence Lengths: Lower Resistance Concern Species",
  xlab = "Sequence Length (bp)",
  col = "skyblue",
  breaks = 20
)


#minimum species

## Sample min ITS sequences per species ----

table(dfITS$Species_Name)

min_n <- min(table(dfITS$Species_Name))  

idx_by_species <- split(seq_len(nrow(dfITS)), dfITS$Species_Name)

set.seed(123)  

# Sample min_n indices per species
idx_sampled <- unlist(
  lapply(idx_by_species, function(idx) sample(idx, min_n))
)

# Subset dataframe
dfITS_subset <- dfITS[idx_sampled, ]

# Check that we have 19 per species
table(dfITS_subset$Species_Name)

# Making sure no missing sequences
sum(is.na(dfITS$ITS_Sequence))


write.csv(dfITS_subset,
          file = "../data/clean/candida6_ITS_minperSpecies.csv",
          row.names = FALSE)

#### ERG111

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


candida6_ERG11_search <- entrez_search(
  db          = "nuccore",
  term        = candida6_ERG11_term,
  retmax      = 500,      
  use_history = TRUE
)

candida6_ERG11_search$count
length(candida6_ERG11_search$ids)

erg11_fasta <- entrez_fetch(
  db          = "nuccore",
  web_history = candida6_ERG11_search$web_history,
  rettype     = "fasta",
  retmax = 500,
  retmode     = "text"
)

write_fasta(erg11_fasta, "../data/raw/candida6_ERG11_raw.fasta")

candida6_ERG11_seqs <- readDNAStringSet("../data/raw/candida6_ERG11_raw.fasta")
length(candida6_ERG11_seqs)

seq_lengths <- width(candida6_ERG11_seqs)
summary(seq_lengths)


dfERG11 <- data.frame(
  ERG11_Title = names(candida6_ERG11_seqs),
  ERG11_Sequence = paste(candida6_ERG11_seqs)
)

# Base R species name extraction 
title_split <- strsplit(dfERG11$ERG11_Title, " ")

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

dfERG11$Species_Name <- species_name
dfERG11 <- dfERG11[, c("ERG11_Title", "Species_Name", "ERG11_Sequence")]
table(dfERG11$Species_Name)

dfERG11$Species_Name[dfERG11$Species_Name == "Candidozyma auris"] <- "Candida auris"

dfERG11$Species_Name[dfERG11$Species_Name == "Pichia kudriavzevii"] <- "Candida krusei"

dfERG11$Species_Name[dfERG11$Species_Name == "Nakaseomyces glabratus"] <- "Candida glabrata"

table(dfERG11$Species_Name)


dfERG11$Resistance_group <- NA

high_res <- c("Candida auris", "Candida glabrata", "Candida krusei")
low_res  <- c("Candida albicans", "Candida parapsilosis", "Candida tropicalis")

dfERG11$Resistance_group[dfERG11$Species_Name %in% high_res] <- "High resistance concern"
dfERG11$Resistance_group[dfERG11$Species_Name %in% low_res]  <- "Lower resistance concern"

dfERG11$Seq_Length <- nchar(dfERG11$ERG11_Sequence)

hist(
  dfERG11$Seq_Length[dfERG11$Resistance_group == "High resistance concern"],
  main = "ERG11 Sequence Lengths: High Resistance Concern Species",
  xlab = "Sequence Length (bp)",
  col = "tomato",
  breaks = 20
)

hist(
  dfERG11$Seq_Length[dfERG11$Resistance_group == "Lower resistance concern"],
  main = "ERG11 Sequence Lengths: Lower Resistance Concern Species",
  xlab = "Sequence Length (bp)",
  col = "skyblue",
  breaks = 20
)


## Sample min ERG11 sequences per species ----

table(dfERG11$Species_Name)

min_n <- min(table(dfERG11$Species_Name))  

idx_by_species <- split(seq_len(nrow(dfERG11)), dfERG11$Species_Name)

set.seed(123)  

# Sample min_n indices per species
idx_sampled <- unlist(
  lapply(idx_by_species, function(idx) sample(idx, min_n))
)

# Subset dataframe
dfERG11_subset <- dfERG11[idx_sampled, ]

# Check that we have 19 per species
table(dfERG11_subset$Species_Name)

# Making sure no missing sequences
sum(is.na(dfERG11$ERG11_Sequence))


write.csv(dfERG11_subset,
          file = "../data/clean/candida6_ERG11_minperSpecies.csv",
          row.names = FALSE)

