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

# Function to calculate module eigengenes and hub genes
calculateModuleEigengenes <- function(gene_data, clusters) {
  hub_genes <- list()  # Initialize list to store hub genes
  results <- list()    # Initialize list to store all results
  
  for (cluster_id in names(clusters)) {
    cluster_genes <- clusters[[cluster_id]]
    
    cluster_rows <- gene_data[gene_data[[1]] %in% cluster_genes, ]
    
    # Create expression matrix and transpose
    expression_matrix <- as.matrix(cluster_rows[, -1])  # Exclude gene_id
    expression_matrix <- t(expression_matrix)
    
    # Assign row names from the original data
    row_names <- cluster_rows[[1]]  # Assuming the first column is gene IDs
    
    colors <- rep(cluster_id, ncol(expression_matrix))
    
    # Calculate eigengenes
    MEList <- moduleEigengenes(expression_matrix, colors = colors)
    
    cat("Eigengene values for Cluster", cluster_id, ":\n")
    print(MEList$eigengenes)
    cat("\n")
    
    # Calculate Pearson correlation coefficients
    cor_coef <- cor(expression_matrix, MEList$eigengenes[, 1])
    
    # Sort indices based on correlation coefficients
    sorted_indices <- order(cor_coef, decreasing = TRUE)
    
    # Find top 15 hub genes
    top_hub_genes <- row_names[sorted_indices[1:15]]
    
    # Store hub genes with ranks
    hub_genes[[cluster_id]] <- data.frame(
      Gene_ID = top_hub_genes,
      Correlation_Coefficient = cor_coef[sorted_indices[1:15]],
      Rank_within_Tricluster = 1:15
    )
    
    # Perform differential expression analysis
    # Transpose back the expression matrix
    expression_matrix <- t(expression_matrix)
    
    # Define group factor and design matrix
    group <- factor(rep(c("Group1", "Group2", "Group3", "Group4"), each = 3))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    # Check dimensions
    if (ncol(design) != ncol(expression_matrix)) {
      stop("Error: Dimensions of design matrix do not match expression matrix.")
    }
    
    # Define contrasts
    contrast_matrix <- makeContrasts(
      Group2 - Group1,
      Group3 - Group1,
      Group4 - Group1,
      Group3 - Group2,
      Group4 - Group2,
      Group4 - Group3,
      levels = design
    )
    
    # Fit linear model
    fit <- lmFit(expression_matrix, design)
    
    # Fit contrasts
    fit_contrast <- contrasts.fit(fit, contrast_matrix)
    
    # Perform empirical Bayes moderation
    fit_contrast <- eBayes(fit_contrast)
    
    # Get top differentially expressed genes for each contrast
    top_genes <- list()
    for (i in 1:ncol(contrast_matrix)) {
      contrast_name <- colnames(contrast_matrix)[i]
      top_genes[[contrast_name]] <- topTable(fit_contrast, coef = i, adjust = "BH")
      
      # Filter genes based on logFC and p-value criteria
      filtered_genes <- top_genes[[contrast_name]][abs(top_genes[[contrast_name]]$logFC) > 2 & top_genes[[contrast_name]]$P.Value < 0.001, ]
      
      # Store results in list
      results[[paste("Cluster", cluster_id, "Contrast", contrast_name)]] <- filtered_genes
    }
  }
  
  # Write results to file
  write_results_to_file(results)
  
  return(hub_genes)
}

# Function to write results to a file
write_results_to_file <- function(results) {
  file_path <- file.path("~/Desktop", "master_file.txt")
  cat("Gene Expression Analysis Results\n", file = file_path, append = TRUE)
  
  for (key in names(results)) {
    cat("\n", key, ":\n", file = file_path, append = TRUE)
    cat("Top Differentially Expressed Genes:\n", file = file_path, append = TRUE)
    print(results[[key]], file = file_path, append = TRUE)
    cat("\n", file = file_path, append = TRUE)
  }
}

# Calculate module eigengenes and hub genes for clusters
hub_genes <- calculateModuleEigengenes(gene_data, clusters)

# Output hub genes
for (cluster_id in names(hub_genes)) {
  cat("Top 15 Hub Genes for Cluster", cluster_id, ":\n")
  print(hub_genes[[cluster_id]])
  cat("\n")
}

cat("Results written to 'master_file.txt' on your desktop.\n")