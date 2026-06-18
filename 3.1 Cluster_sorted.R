#rm(list = ls()) #delete environment
#Using MEGENA v.1.6 - installed from https://github.com/songw01/MEGENA/blob/master/MEGENA_1.6.tar.gz

library(Seurat)
library(SeuratObject)
library(MEGENA)

#Get the cluster index from the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
clusterid <- as.numeric(args[1])

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

# Construct file paths based on clusterid
cluster_dir <- paste0("/Cluster",clusterid)
output_dir <- paste0("/Cluster",clusterid)

#load pairwise biologically meaningful correlation (BMC) matrix
ijw_filename <- paste0(cluster_dir, "/BMC.RDS")
ijw <- readRDS(file = ijw_filename)

#make the weight an absolute value and sort them in descending order

ijw <- ijw %>%
  mutate(across(3, abs)) %>%
  arrange(desc(.data[[3]]))

#calculate PFN
if (doPar & getDoParWorkers() == 1)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}

el <- calculate.PFN(ijw, doPar = doPar, num.cores = n.cores, keep.track = FALSE)

output_file <- paste0(output_dir, "/MEGENA_network.txt")
write.table(el, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#Do multiscale clustering
g <- graph.data.frame(el,directed = F)
if (doPar & getDoParWorkers() == 1)
{
  set.parallel.backend(num.cores = n.cores)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}

MEGENA.output <- do.MEGENA(g,
                           mod.pval = module.pval, hub.pval = hub.pval, remove.unsig = TRUE,
                           #min.size = 20,#max.size = vcount(g)/2,
                           doPar = TRUE, num.cores = n.cores, n.perm = 100,
                           save.output = TRUE)

save_file <- paste0(output_dir, "/MEGENA_ouput.RData")
save(MEGENA.output, file = save_file)

quit(save = "no",status = 0)

