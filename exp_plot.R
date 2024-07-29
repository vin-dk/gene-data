library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

file_path <- "C:/Users/13046/Desktop/data_set.xlsx"

# Read in the gene expression data from the Excel file
gene_data <- read_excel(file_path)

new_gene_list <- c("236306_at", "238999_at", "243255_at", "211835_at", 
                   "1563117_at", "237128_at", "217974_at", "223706_at", 
                   "241270_at", "230536_at", "228981_at", "233138_at", 
                   "244176_at", "1564050_at", "216928_at")

# Filter the gene_data to include only the genes in the new_gene_list
filtered_gene_data <- gene_data %>% filter(gene_ref %in% new_gene_list)

# Compute average expression values for each time point
averaged_data <- filtered_gene_data %>%
  mutate(
    zero = rowMeans(select(., starts_with("zero"))),
    three = rowMeans(select(., starts_with("three"))),
    six = rowMeans(select(., starts_with("six"))),
    twelve = rowMeans(select(., starts_with("twelve")))
  ) %>%
  select(gene_ref, zero, three, six, twelve) %>%
  pivot_longer(cols = zero:twelve, names_to = "time_point", values_to = "expression")

averaged_data$time_point <- factor(averaged_data$time_point, levels = c("zero", "three", "six", "twelve"), labels = c("0", "3", "6", "12"))

# Plot the expression profiles
plot <- ggplot(averaged_data, aes(x = time_point, y = expression, group = gene_ref, color = gene_ref)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene Expression Profiles", x = "Time Points (hours)", y = "Expression Level") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid.major = element_line(color = "gray"),
    panel.grid.minor = element_line(color = "gray"),
    text = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    legend.background = element_rect(fill = "black"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    plot.title = element_text(color = "white")
  )

# Save
ggsave("C:/Users/13046/Desktop/cluster7_exp_graph.png", plot)
