#load library
library(dplyr)

#choose cluster, use cluster 1 as example

cluster_id <- 1
cluster_dir <- paste0("Cluster", cluster_id)

#Extract the signed cor with FDR < 0.05 
#example, mincell cutoff = 40 for the cluster

for (cutoff in 40) {
  
  # Construct file name
  rds_file <- file.path(cluster_dir, paste0("ijw_FDR_0.05_pos_neg_mincell", cutoff, ".rds")
  
  # Load the RDS file
  ijw_FDR <- readRDS(rds_file)
  
  # Create a subset of 'signif.ijw'
  ijw_pos_subset <- ijw_FDR[[1]] 
  ijw_neg_subset <- ijw_FDR[[2]] 
  ijw_both_subset <- rbind(ijw_pos_subset, ijw_neg_subset)

  # Create a dynamic variable name
  var_name <- paste0("ijw_FDR_0.05_subset_mincell", cutoff)
  
  # Save the subset to the R environment with the dynamic name
  assign(var_name, ijw_both_subset)
  }

#Extract mean value
Mean_list <- readRDS(paste0(cluster_dir,"PDControl_avg_signed_Mean_list.rds"))
MSE_list <- readRDS(paste0(cluster_dir, "PDControl_avg_signed_MSE_list.rds"))

Mean <- Mean_list[["40"]]

#Exclude the (0, mean)
ijw_FDR_0.05_substract <- ijw_FDR_0.05_subset_mincell40 %>%
  filter(rho <= 0 | rho >= Mean)

#Calculate new rho - BMC
BMC <- ijw_FDR_0.05_substract

saveRDS(BMC, paste0(cluster_dir, "BMC.RDS")
