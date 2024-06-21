# Load necessary libraries
library(readxl)
library(WGCNA)

# Read gene expression data from Excel file
gene_data <- read_excel("C:/Users/13046/Desktop/data_set.xlsx")

# Define the path to the text file containing cluster information
file_path <- "C:/Users/13046/Desktop/Work Paper/trials/david_file.txt"

# Read the file into a character vector, each element is a line from the file
lines <- readLines(file_path)

# Initialize a list to store clusters
clusters <- list()

# Iterate through the lines to extract cluster info
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

# Check if there are any empty clusters
if (any(sapply(clusters, length) == 0)) {
  cat("Warning: There are empty clusters.\n")
}

# Remove empty clusters
clusters <- clusters[sapply(clusters, length) > 0]

# Check if there are any clusters left
if (length(clusters) == 0) {
  stop("Error: No non-empty clusters found.")
}

# For each cluster, extract the relevant rows from gene_data, transpose, and calculate eigengene
for (cluster_id in names(clusters)) {
  cluster_genes <- clusters[[cluster_id]]
  
  # Find rows in gene_data that match the cluster genes
  cluster_rows <- gene_data[gene_data[[1]] %in% cluster_genes, ]
  
  # Convert cluster data to matrix format (rows as genes, columns as samples)
  expression_matrix <- as.matrix(cluster_rows[, -1])  # Exclude the first column (gene ID)
  
  # Transpose expression matrix
  expression_matrix <- t(expression_matrix)
  
  # Assign a unique color for each gene (each column in the transposed matrix)
  complex_colors_genes <- paste0("color", 1:ncol(expression_matrix))
  
  # Calculate eigengene for the current cluster
  MEList_genes <- moduleEigengenes(expression_matrix, colors = complex_colors_genes)
  
  # Print eigengene values for the current cluster
  cat("Eigengene values for Cluster", cluster_id, ":\n")
  print(MEList_genes$eigengenes)
  cat("\n")
}