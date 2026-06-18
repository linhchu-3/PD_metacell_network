#iterative pruning procedure to remove genes and nuclei in a data matrix

library(Seurat)
library(SeuratObject)
library(MEGENA)
library(MAST)

obj <-readRDS("/MC_all_clusters_gamma15_avg_PC10_knn5_n2000_unweighted.seurat.RDS")

logfc.threshold <- log2(1)

#Get the cluster index from the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
clusterid <- as.numeric(args[1])

# Get the expression data for the specified cluster
expr_dat <- GetAssayData(obj[, Idents(obj) == clusterid], slot = 'data')
expr_dat = expr_dat[rowSums(expr_dat > 0) >= 0.2 * ncol(expr_dat), ]

# Define the range of T_g and T_m values
T_g_values <- seq(0.1, 0.5, by = 0.05)  # Example: T_g from 0.1 to 0.5
T_m_values <- seq(0.1, 0.5, by = 0.05)  # Example: T_m from 0.1 to 0.5

# Initialize a list to store DEG results for each combination of T_g and T_m
deg_results_list <- list()
filtered_expr_data_list <- list()

# Loop over cutoff values from 0.1 to 0.5 with an interval of 0.05
for (T_g in T_g_values) {
  # Loop over T_m values
  for (T_m in T_m_values) {
  # Filter genes: Remove genes with exceed % missing rate (zero counts)
  filtered_expr_dat <- expr_dat[rowMeans(expr_dat == 0) <= T_g, ]
  
  # Filter metacells: remove metacells with zero counts for all genes
  filtered_expr_dat <- filtered_expr_dat[, colMeans(filtered_expr_dat == 0) <= T_m]
  
  # Subset the original Seurat object to retain only the filtered cells and genes
  subset_obj <- subset(
    obj,  # Original Seurat object
    cells = colnames(filtered_expr_dat),  # Retain only the filtered cells
    features = rownames(filtered_expr_dat)  # Retain only the filtered genes
  )
  
  #Store the filtered expression data in the list
  filtered_expr_data_list[[paste0("T_g_", T_g, "_T_m_", T_m)]] <- filtered_expr_dat
  
  # Perform DEG analysis using MAST for the current cutoff
  r1 <- try(
    FindMarkers(
      subset_obj,
      slot = 'data',
      assay = "RNA",
      ident.1 = "PD",
      ident.2 = 'Control',
      group.by = 'Dx',
      test.use = "MAST",
      subset.ident = clusterid,
      logfc.threshold = logfc.threshold
    ),
    silent = TRUE
  )
  
  # If DEG analysis succeeds, save the results
  if (!inherits(r1, 'try-error')) {
    # Add metadata columns to the DEG results
    r1 <- data.frame(
      Geneid = rownames(r1),
      Cluster = clusterid,
      Contrast = 'PD-vs-Control',
      T_g = T_g,
      T_m = T_m,
      r1,
      stringsAsFactors = FALSE
    )
  }
  
  # Store the DEG results in the list
  deg_results_list[[paste0("T_g_", T_g, "_T_m_", T_m)]] <- r1
  
  # Save the DEG results to a CSV file
  filename <- paste0("DEG_Results_T_g_", T_g, "_T_m_", T_m, "_Cluster_", clusterid, ".csv")
  write.csv(r1, filename, row.names = FALSE)
  
  #Save the matrix and deg list
  saveRDS(filtered_expr_data_list, "filtered_expr_data_list.RDS")
  saveRDS(deg_results_list, "deg_results_list.RDS")
  }
}


