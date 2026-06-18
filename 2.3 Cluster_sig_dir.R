#rm(list = ls()) #delete environment
#Using MEGENA v.1.6 - need to uninstall prev MEGENA and install new one

library(Seurat)
library(SeuratObject)
library(MEGENA)

#Get the cluster index from the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
clusterid <- as.numeric(args[1])

# Construct file paths based on clusterid
cluster_dir <- paste0("/PD/SnRNAseq/MEGENA_no0_Metacells_corrected/Gamma15_avg_PC10_knn5_n2000_mincell_opt/unweighted/PDControl/Cluster",clusterid)
markers_file <- paste0("/PD/SnRNAseq/PD_Metacells_corrected/Gamma15_avg_PC10_knn5_n2000/Conserved markers/Unweighted/Cluster", clusterid, "_unweighted_conserved_markers.csv")

#load file
obj <- readRDS("/PD/SnRNAseq/PD_Metacells_corrected/Gamma15_avg_PC10_knn5_n2000/MC_all_clusters_gamma15_avg_PC10_knn5_n2000_unweighted.seurat.RDS")

#set some parameters

n.cores <- 14; #number of cores/threads to call for PCP
doPar <- TRUE; #do we want to parallelize?
method = "spearman" #method for correlation. either pearson or spearman.
FDR.cutoff = 0.05 #FDR threshold to define significant correlation upon shuffling samples
module.pval = 0.05 #module significance p-value. Recommend is 0.05
hub.pval = 0.05 #module significance p-value. Recommend is 0.05
#cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
#hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 

# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2

#Find conserved markers
Cluster_conserved_markers <- read.csv(markers_file)
colnames(Cluster_conserved_markers)[1] <- "Geneid"

obj.cls = obj[,(Idents(obj) == clusterid)]

#loop through differen min.cell cutoffs (from 0.1 to 0.5)

for (cutoff in seq(0.1, 0.4, by = 0.1)) {
  min.cells <- ncol(obj.cls)* cutoff #Calculate the minimum number of cells expressed
  cnt <- GetAssayData(obj.cls, slot = "counts") 
  
  #filter expressed genes
  nc <- rowSums(cnt >= 1E-320,na.rm = TRUE)
  ii = which(nc >= min.cells) # index of expressed genes
  ii2 = intersect(ii, which(is.element(rownames(obj@assays$RNA@data), Cluster_conserved_markers$Geneid)))
  print(length(ii2))
  
  
  #remove noisy cells/nuclei
  min.genes = 500;
  ng = colSums(cnt[ii2,] > 1E-320,na.rm = TRUE)
  jj = which(ng >= min.genes)
  print(length(jj))

  #Extract the expression data matrix
   
  datExpr = GetAssayData(obj.cls,slot = "data")[ii2,jj]

  print(class(datExpr))
  print(dim(datExpr))
  print(datExpr[1:5, 1:5])

  #replace all 0 to NA 
  datExpr[datExpr == 0] <- NA
  print(datExpr[1:5, 1:5])

  #change it to matrix

  datExpr <- as.matrix(datExpr)
  print(class(datExpr))
  print(dim(datExpr))

 # Perform pairwise correlation without applying the FDR cutoff
 
  ijw_pos <- calculate.rho.signed(datExpr, 
                              n.perm = 10, 
                              FDR.cutoff = FDR.cutoff, 
                              estimator = method,
                              use.obs = "pairwise.complete.obs",
                              direction = "positive",
                              rho.thresh = NULL,
                              sort.el = TRUE)

  ijw_neg <- calculate.rho.signed(datExpr, 
                              n.perm = 10, 
                              FDR.cutoff = FDR.cutoff, 
                              estimator = method,
                              use.obs = "pairwise.complete.obs",
                              direction = "negative",
                              rho.thresh = NULL,
                              sort.el = TRUE)

  ijw <- rbind(ijw_pos, ijw_neg)

  # Save the correlation results (ijw) as an RDS file with a unique name based on the cutoff
  ijw_filename <- paste0(cluster_dir, "/ijw_FDR_0.05_pos_neg_mincell", as.integer(cutoff * 100), ".rds")
  saveRDS(ijw, file = ijw_filename)
}


