# Load necessary libraries
library(readxl)
library(WGCNA)
library(limma)

# Expression data from an Excel file
gene_data <- read_excel("C:/Users/13046/Desktop/data_set.xlsx")

# Path to cluster info
file_path <- "C:/Users/13046/Desktop/Work Paper/trials/david_file.txt"
lines <- readLines(file_path)

# Initialize cluster list
clusters <- list()

# Extract cluster info
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

# Verify non-empty clusters
if (any(sapply(clusters, length) == 0)) {
  cat("Warning: There are empty clusters.\n")
}

# Remove empty clusters
clusters <- clusters[sapply(clusters, length) > 0]

# Check if clusters exist
if (length(clusters) == 0) {
  stop("Error: No non-empty clusters found.")
}

# Initialize lists to store row names and original matrices
row_names_list <- list()
original_matrices <- list()

# Function to calculate module eigengenes and hub genes
calculateModuleEigengenes <- function(gene_data, clusters, row_names_list, original_matrices) {
  hub_genes <- list()  # Initialize list to store hub genes
  eigengenes <- list()  # Initialize list to store eigengenes
  
  for (cluster_id in names(clusters)) {
    cluster_genes <- clusters[[cluster_id]]
    
    cluster_rows <- gene_data[gene_data[[1]] %in% cluster_genes, ]
    
    # Store row names
    row_names_list[[cluster_id]] <- cluster_rows[[1]]  # Assuming the first column has gene IDs
    
    # Store original untransposed matrix
    original_matrices[[cluster_id]] <- as.matrix(cluster_rows[, -1])  # Exclude gene_id
    
    expression_matrix <- t(original_matrices[[cluster_id]])
    
    colors <- rep(cluster_id, ncol(expression_matrix))
    
    MEList <- moduleEigengenes(expression_matrix, colors = colors)
    
    # Store eigengene values
    eigengenes[[cluster_id]] <- MEList$eigengenes
    
    cat("Cluster", cluster_id, ":\n")
    cat("Eigengene values:\n")
    print(MEList$eigengenes)
    cat("\n")
    
    # Calculate Pearson correlation coefficients
    cor_coef <- cor(expression_matrix, MEList$eigengenes[, 1])
    
    # Sort indices based on correlation coefficients
    sorted_indices <- order(cor_coef, decreasing = TRUE)
    
    # Find top 15 hub genes
    top_hub_genes <- row_names_list[[cluster_id]][sorted_indices[1:15]]
    
    # Store hub genes with ranks
    hub_genes[[cluster_id]] <- data.frame(
      Gene_ID = top_hub_genes,
      Correlation_Coefficient = cor_coef[sorted_indices[1:15]],
      Rank_within_Tricluster = 1:15
    )
  }
  
  return(list(hub_genes = hub_genes, eigengenes = eigengenes, original_matrices = original_matrices))
}

# Function to replace row numbers with gene IDs
replaceRowNumbersWithGeneIDs <- function(row_numbers, row_names) {
  return(row_names[row_numbers])
}

# Calculate module eigengenes and hub genes for clusters
result <- calculateModuleEigengenes(gene_data, clusters, row_names_list, original_matrices)

# Perform differential expression analysis for each cluster and redirect output to file
sink("C:/Users/13046/Desktop/master_file.txt")

for (cluster_id in names(result$hub_genes)) {
  cat("Cluster", cluster_id, ":\n")
  cat("Eigengene values:\n")
  print(result$eigengenes[[cluster_id]])
  cat("\n")
  
  cat("Hub genes:\n")
  print(result$hub_genes[[cluster_id]])
  cat("\n")
  
  # Perform differential expression analysis
  original_matrix <- result$original_matrices[[cluster_id]]
  row_names <- row_names_list[[cluster_id]]
  
  n_rows <- nrow(original_matrix)
  expression_matrix <- matrix(as.vector(t(original_matrix)), nrow = n_rows, byrow = TRUE)
  
  group <- factor(rep(c("Group1", "Group2", "Group3", "Group4"), each = 3))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  contrast_matrix <- makeContrasts(
    Group2 - Group1,
    Group3 - Group1,
    Group4 - Group1,
    Group3 - Group2,
    Group4 - Group2,
    Group4 - Group3,
    levels = design
  )
  
  fit <- lmFit(expression_matrix, design)
  fit_contrast <- contrasts.fit(fit, contrast_matrix)
  fit_contrast <- eBayes(fit_contrast)
  
  for (i in 1:ncol(contrast_matrix)) {
    contrast_name <- colnames(contrast_matrix)[i]
    top_genes <- topTable(fit_contrast, coef = i, adjust = "BH")
    filtered_genes <- top_genes[abs(top_genes$logFC) > 2 & top_genes$P.Value < 0.001, ]
    
    # Replace row numbers with gene IDs
    filtered_genes$Gene_ID <- replaceRowNumbersWithGeneIDs(filtered_genes$Gene_ID, row_names)
    
    cat("Top Differentially Expressed Genes for", contrast_name, ":\n")
    print(filtered_genes)
    cat("\n")
  }
}

sink()  # Close the sink