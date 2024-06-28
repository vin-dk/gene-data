library(readxl)
library(WGCNA)
library(limma)

# Raw data
gene_data <- read_excel("C:/Users/13046/Desktop/data_set.xlsx")

options(max.print = 10000, width = 10000)

# Path to cluster info
file_path <- "C:/Users/13046/Desktop/Work Paper/trials/steps/Trial 1_1/david_file_1.txt"
lines <- readLines(file_path)

# Initialize cluster list
clusters <- list()

# Extract cluster info, only if > 100 genes
i <- 1
while (i <= length(lines)) {
  if (grepl("^Block", lines[i])) {
    cluster_id <- sub(".*?(\\d+).*", "\\1", lines[i])  # Extract the cluster number
    i <- i + 2  # Skip the block info line
    
    
    cluster_genes <- unlist(strsplit(lines[i], ",\\s*"))
    
    
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

# Initialize lists to store original matrices
original_matrices <- list()

# Function to calculate module eigengenes and hub genes
calculateModuleEigengenes <- function(gene_data, clusters, original_matrices) {
  hub_genes <- list()  
  eigengenes <- list()  
  
  for (cluster_id in names(clusters)) {
    cluster_genes <- clusters[[cluster_id]]
    
    cluster_rows <- gene_data[gene_data[[1]] %in% cluster_genes, ]
    
    # Store original untransposed matrix
    original_matrices[[cluster_id]] <- as.matrix(cluster_rows[, -1])  
    
    expression_matrix <- t(original_matrices[[cluster_id]])
    
    colors <- rep(cluster_id, ncol(expression_matrix))
    
    MEList <- moduleEigengenes(expression_matrix, colors = colors)
    
    # Store eigengene values
    eigengenes[[cluster_id]] <- MEList$eigengenes
    
    cat("Cluster", cluster_id, ":\n")
    cat("Eigengene values:\n")
    print(MEList$eigengenes)
    cat("\n")
    
    # Write eigengene values to master_file_2
    sink("C:/Users/13046/Desktop/master_file_2.txt", append = TRUE)
    cat(sprintf("Cluster %s:\n", cluster_id))
    cat("Eigengene values:\n")
    print(MEList$eigengenes)
    cat("\n")
    sink()  
    
    # Calculate Pearson correlation coefficients
    cor_coef <- cor(expression_matrix, MEList$eigengenes[, 1])
    
    # Sort indices based on correlation coefficients
    sorted_indices <- order(cor_coef, decreasing = TRUE)
    
    # Find top 15 hub genes
    top_hub_genes <- cluster_genes[sorted_indices[1:15]]
    
    # Store with ranks
    hub_genes[[cluster_id]] <- data.frame(
      Gene_ID = top_hub_genes,
      Correlation_Coefficient = cor_coef[sorted_indices[1:15]],
      Rank_within_Tricluster = 1:15
    )
    
    # Write hub genes to master_file_2
    sink("C:/Users/13046/Desktop/master_file_2.txt", append = TRUE)
    cat(sprintf("Cluster %s:\n", cluster_id))
    cat("Hub genes:\n")
    print(hub_genes[[cluster_id]])
    cat("\n")
    sink()  
  }
  
  return(list(hub_genes = hub_genes, eigengenes = eigengenes, original_matrices = original_matrices))
}

# Calculate module eigengenes and hub genes for clusters
result <- calculateModuleEigengenes(gene_data, clusters, original_matrices)

# Differential analysis
sink("C:/Users/13046/Desktop/master_file.txt")


diff_exp_genes <- list()

for (cluster_id in names(result$hub_genes)) {
  cat("Cluster", cluster_id, ":\n")
  cat("Eigengene values:\n")
  print(result$eigengenes[[cluster_id]])
  cat("\n")
  
  cat("Hub genes:\n")
  print(result$hub_genes[[cluster_id]])
  cat("\n")
  
  # Analysis portion
  original_matrix <- result$original_matrices[[cluster_id]]
  row_names <- clusters[[cluster_id]]  
  
  n_rows <- nrow(original_matrix)
  expression_matrix <- matrix(as.vector(t(original_matrix)), nrow = n_rows, byrow = TRUE)
  
  while (length(row_names) != nrow(expression_matrix)) {
    row_names <- head(row_names, -1)  
    rownames(expression_matrix) <- row_names
    n_rows <- nrow(original_matrix)
    expression_matrix <- matrix(as.vector(t(original_matrix)), nrow = n_rows, byrow = TRUE)
  }
  
  rownames(expression_matrix) <- row_names
  
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
  
  diff_exp_genes[[cluster_id]] <- c()  
  
  for (i in 1:ncol(contrast_matrix)) {
    contrast_name <- colnames(contrast_matrix)[i]
    top_genes <- topTable(fit_contrast, coef = i, adjust = "BH", number = Inf)
    
    cat("Top Differentially Expressed Genes for", contrast_name, ":\n")
    print(top_genes)
    cat("\n")
    
    sig_genes <- rownames(top_genes)[top_genes$adj.P.Val < 0.05]
    diff_exp_genes[[cluster_id]] <- unique(c(diff_exp_genes[[cluster_id]], sig_genes))
  }
  
  cat("\n")
}

sink()


all_diff_exp_genes <- unique(unlist(diff_exp_genes))

# make the differential expression summary
sink("C:/Users/13046/Desktop/diff_exp_summary.txt")
for (cluster_id in names(diff_exp_genes)) {
  total_genes <- length(clusters[[cluster_id]])
  num_unique_genes <- length(diff_exp_genes[[cluster_id]])
  
  percentage <- (num_unique_genes / total_genes) * 100
  
  cat("Cluster", cluster_id, ":\n")
  cat("Total number of genes in the cluster:", total_genes, "\n")
  cat("Number of unique differentially expressed genes:", num_unique_genes, "\n")
  cat("Percentage of differentially expressed genes:", percentage, "%\n\n")
}
sink()

# EOF