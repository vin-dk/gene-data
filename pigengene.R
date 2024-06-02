# Load necessary library
library(readxl)
library(Pigengene)
library(HDO.db)

# Read gene expression data from Excel file
gene_data <- read_excel("C:/Users/13046/Desktop/data_set.xlsx")

# Define the path to the text file
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
    
    # Store the genes in the clusters list
    clusters[[cluster_id]] <- cluster_genes
  }
  i <- i + 1
}

# Print the clusters to verify
print(clusters)

# Initialize a list to store the clustered data
clustered_data <- list()

# For each cluster, extract the relevant rows from gene_data
for (cluster_id in names(clusters)) {
  cluster_genes <- clusters[[cluster_id]]
  
  # Find rows in gene_data that match the cluster genes
  cluster_rows <- gene_data[gene_data[[1]] %in% cluster_genes, ]
  
  # Store the cluster data in the list
  clustered_data[[cluster_id]] <- cluster_rows
}

# Print the clustered data to verify
print(clustered_data)