library(readxl)
library(WGCNA)
library(limma)

#This is a modified version of the main file to do eigengene analysis

# Raw data
gene_data <- read_excel("C:/Users/13046/Desktop/data_set.xlsx")

options(max.print = 10000, width = 10000)

# Path to cluster info
file_path <- "C:/Users/13046/Desktop/Work Paper/trials/steps/Trial 3_10/david_file_30.txt"
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

original_matrices <- list()

# Function to calculate module eigengenes and hub genes
calculateModuleEigengenes <- function(gene_data, clusters, original_matrices) {
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
  }
  
  return(list(eigengenes = eigengenes, original_matrices = original_matrices))
}

# Calculate module eigengenes and hub genes for clusters
result <- calculateModuleEigengenes(gene_data, clusters, original_matrices)

# Initialize counters for differential expression
de_counts <- list(
  "T3h-T0h" = 0,
  "T6h-T0h" = 0,
  "T12h-T0h" = 0,
  "T6h-T3h" = 0,
  "T12h-T3h" = 0,
  "T12h-T6h" = 0
)
total_eigengenes <- 0

# Differential analysis
sink("C:/Users/13046/Desktop/eigengene_summary.txt")

for (i in seq_along(result$eigengenes)) {
  eigengene_id <- names(result$eigengenes)[i]
  eigengene <- result$eigengenes[[eigengene_id]]
  total_eigengenes <- total_eigengenes + 1
  
  # transpose
  eigengene_matrix <- t(eigengene)
  
  # design matrix
  time_points <- rep(c("T0h", "T3h", "T6h", "T12h"), each = 3)
  design <- model.matrix(~ 0 + factor(time_points))
  colnames(design) <- c("T0h", "T3h", "T6h", "T12h")
  
  
  fit <- lmFit(eigengene_matrix, design)
  
  contrast.matrix <- makeContrasts(
    "T3h-T0h", "T6h-T0h", "T12h-T0h",
    "T6h-T3h", "T12h-T3h",
    "T12h-T6h",
    levels = design
  )
  
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  
  cat(sprintf("Eigengene %s:\n", eigengene_id))
  print(eigengene)
  
  comparisons <- c("T3h-T0h", "T6h-T0h", "T12h-T0h", "T6h-T3h", "T12h-T3h", "T12h-T6h")
  
  for (comp in comparisons) {
    res <- topTable(fit2, coef = comp, adjust.method = "BH", number = Inf)
    
    de <- ifelse(any(res$adj.P.Val < 0.05), "Yes", "No")
    
    cat(sprintf("%s: %s\n", comp, de))
    
    if (de == "Yes") {
      de_counts[[comp]] <- de_counts[[comp]] + 1
    }
    
    
    cat("Numeric Information:\n")
    print(res[, c("logFC", "adj.P.Val")])
    cat("\n")
  }
  cat("\n")
}


cat("Summary of differential expression:\n")
for (comp in names(de_counts)) {
  cat(sprintf("%s: %d/%d (%.2f%%)\n", 
              comp, de_counts[[comp]], total_eigengenes, 
              100 * de_counts[[comp]] / total_eigengenes))
}

sink()  # Close the sink