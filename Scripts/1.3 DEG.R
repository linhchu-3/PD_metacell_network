library(Seurat)
library(SeuratObject)
library(MEGENA)
library(MAST)

MC_all_clusters_unweighted.seurat <-readRDS("/MC_all_clusters_gamma15_avg_PC10_knn5_n2000_unweighted.seurat.RDS")

logfc.threshold <- log2(1)

cluster_ids <- (0:11)

#unweighted

#loop through different clusters

for (clusterid in cluster_ids){
  cat('Do Cluster', clusterid, '...\n')
  result_all_clusters <- NULL
  #Find DEGs
  r1 <- try(FindMarkers(MC_all_clusters_unweighted.seurat, slot='data', assay = "RNA", ident.1 = "PD", ident.2 = 'Control', group.by = 'Dx', test.use = "MAST", subset.ident = clusterid, logfc.threshold = logfc.threshold), silent = TRUE)
  if(!inherits(r1, 'try-error'))  
  result_all_clusters  = rbind(result_all_clusters, data.frame(Geneid= rownames(r1), Cluster = clusterid, Contrast = 'PD-vs-Control', r1, stringsAsFactors =  FALSE))
  
  #write results
  filename <- paste0("Cluster", clusterid, "_unweighted_DEGs_MAST.csv",sep ="")
  write.csv(result_all_clusters, filename)
}

##filter based on the threshold
#unweighted
results_mc_avg_unw_df <- data.frame()

clusterid_values <- c(0:11)
for (clusterid in clusterid_values) {
      # Construct the file name
      file_name <- paste0("/Cluster", clusterid, "_unweighted_DEGs_MAST.csv")
      # Read the data from the file
      data <- read.csv(file_name)  # Modify according to your file structure
      # Create a new column in the results dataframe
      data$module  <- paste("Cluster", clusterid)
      data$Condition <- "Avg_unw"
      # Append the new row to the results dataframe
      results_mc_avg_unw_df <- rbind(results_mc_avg_unw_df, data)
}


#choose the threshold, here, we use results with FC = 1.1, padj < 0.05
results_mc_avg_unw_df_padj <- filter(results_mc_avg_unw_df, p_val_adj < 0.05)
results_mc_avg_unw_df_up <- filter(results_mc_avg_unw_df_padj, avg_log2FC > log2(1.1))
results_mc_avg_unw_df_dn <- filter(results_mc_avg_unw_df_padj, avg_log2FC < log2(1/1.1))
results_mc_avg_unw_df_sig <- rbind(results_mc_avg_unw_df_up, results_mc_avg_unw_df_dn)

for (i in 0:11) {
  # Filter the data for the current module
  result_filtered <- dplyr::filter(results_mc_avg_unw_df_sig, module == paste("Cluster", i))
  result_filtered$FC <- 2^(result_filtered$avg_log2FC)
  result_filtered <- result_filtered %>%
    select(-X)
  # Define the filename
  filename <- paste0("result_cluster_", i, ".csv")
  assign(paste0("result_mc_avg_unw_C", i, "_sig"), result_filtered)
}
