#rm(list = ls()) #delete environment

#load library
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(scran)
library(cowplot)
library(Matrix)
library(SuperCell)

options(future.globals.maxSize = 4000 * 1024^2)

#load file
obj <- readRDS("dbl.classification.SeuratObject.clean.Cluster_final.RDS")
obj.cor <- readRDS("obj.cor.RDS")

## Simplify single-cell data at

gamma <- 15 # graining level
k.knn <- 5
hvg <- VariableFeatures(obj) #number of features = 2000
pca <- 10 #PC to generate metacells
clusterid_values  <- 0:11 

#avg method

#Loop through different clustersid and Dx 
for (clusterid in clusterid_values){
  for (Dx in c("Control", "PD")) {
    #subset data based on cluster ID and Dx
    obj_subset <- obj.cor[,Idents(obj.cor) == clusterid & obj.cor@meta.data$Dx == Dx]
    
    #compute metacells using SuperCell Package
    MC <- SCimplify(
      X = GetAssayData(obj_subset),
      gene.use = hvg,
      gamma = gamma,
      n.pc = pca
    )
    
    #Annotate metacells to cluster and Dx
    MC$Cluster <- clusterid
    MC$Dx <- Dx
   
    # Compute gene expression of metacells by simply averaging gene expression within each metacell
    MC.ge <- supercell_GE(
      ge = GetAssayData(obj_subset),
      groups = MC$membership
    )
    
    #convert to Seurat object
    MC.seurat <- supercell_2_Seurat(
      SC.GE = MC.ge,
      SC = MC, 
      fields = c("Cluster", "Dx"),
      var.genes = MC$genes.use, 
      N.comp = pca
    )
    
    #Set Cluster identities
    Idents(MC.seurat) <- "Cluster"
    
    #Update Cluster and Dx in Seurat object
    MC.seurat@meta.data$Cluster <- as.factor(clusterid)
    MC.seurat@meta.data$Dx <- Dx
    
    #Save Seurat object
    filename <- paste("MC_Cluster", clusterid, "_", Dx, "_avg.seurat.RDS", sep="")
    saveRDS(MC.seurat, filename)
  }
}


