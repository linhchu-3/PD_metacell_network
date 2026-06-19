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
vars.to.regress = c('Sex', 'Age', 'PMI')

all.genes <- rownames(obj)

obj.cor <- ScaleData(obj, features = all.genes, vars.to.regress = vars.to.regress)

saveRDS(obj.cor, "obj.cor.RDS")

