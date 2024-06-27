# Load necessary libraries
library(readxl)
library(WGCNA)

# expression data of excel file
gene_data <- read_excel("C:/Users/13046/Desktop/data_set.xlsx")

# path to cluster info
file_path <- "C:/Users/13046/Desktop/Work Paper/trials/david_file.txt"

lines <- readLines(file_path)

# init cluster list
clusters <- list()

# extract cluster info
i <- 1
while (i <= length(lines)) {
  if (grepl("^Block", lines[i])) {
    cluster_id <- sub(".*?(\\d+).*", "\\1", lines[i])  # Extract the cluster number
    i <- i + 2  # Skip the block info line
    
    # Extract genes from cluster info line
    cluster_genes <- unlist(strsplit(lines[i], ",\\s*"))
    
    # Store the genes in the clusters list if the cluster has more than 100 genes
    if (length(cluster_genes) > 100) {
      clusters[[cluster_id]] <- cluster_genes
    }
  }
  i <- i + 1
}

# verify no empty
if (any(sapply(clusters, length) == 0)) {
  cat("Warning: There are empty clusters.\n")
}

# remove if necessary
clusters <- clusters[sapply(clusters, length) > 0]

# check for no clusters
if (length(clusters) == 0) {
  stop("Error: No non-empty clusters found.")
}

# main func, calcs eigengene for each defined cluster (according to text file), finds hub-genes
calculateModuleEigengenes <- function(gene_data, clusters) {
  hub_genes <- list()  # Initialize list to store hub genes
  
  for (cluster_id in names(clusters)) {
    cluster_genes <- clusters[[cluster_id]]
    
    cluster_rows <- gene_data[gene_data[[1]] %in% cluster_genes, ]
    
    expression_matrix <- as.matrix(cluster_rows[, -1])  # exclude gene_id
    
    expression_matrix <- t(expression_matrix)
    
    colors <- rep(cluster_id, ncol(expression_matrix))
    
    MEList <- moduleEigengenes(expression_matrix, colors = colors)
    
    cat("Eigengene values for Cluster", cluster_id, ":\n")
    print(MEList$eigengenes)
    cat("\n")
    
    # Calc pearson
    cor_coef <- cor(expression_matrix, MEList$eigengenes[, 1])
    
    # Sort indicies
    sorted_indices <- order(cor_coef, decreasing = TRUE)
    
    # find top 15 hub_genes
    top_hub_genes <- cluster_genes[sorted_indices[1:15]]
    
    # store with ranks
    hub_genes[[cluster_id]] <- data.frame(
      Gene_ID = top_hub_genes,
      Correlation_Coefficient = cor_coef[sorted_indices[1:15]],
      Rank_within_Tricluster = 1:15
    )
  }
  
  return(hub_genes)
}

# grab
hub_genes <- calculateModuleEigengenes(gene_data, clusters)

# output of hub-genes
for (cluster_id in names(hub_genes)) {
  cat("Top 15 Hub Genes for Cluster", cluster_id, ":\n")
  print(hub_genes[[cluster_id]])
  cat("\n")
}