#update prior path
prior_path_new <-  "/TF_target_BN_priors.tsv"

#use example cluster 1
clusterid = 1

filtered_expr_data_list <- readRDS(paste0("Cluster", clusterid, "/filtered_expr_data_list.RDS"))

expr_data_filt <- filtered_expr_data_list[["T_g_0.4_T_m_0.4"]]

expr_dat_new = expr_data_filt[apply(expr_data_filt, 1, function(x) length(unique(x))) > 3, ]

#discretization fnx
discretizeK = function(x, k = 3){
    y = kmeans(x, k)
    order(y$centers)[y$cluster] - 1
}

expr_dat_3 <- t(apply(expr_dat_new, 1, discretizeK))

#prepare output directory #use Tg = 0.4 and Tm = 0.4
                                    
outdir = paste0('Cluster', clusterid, '/BN_Tg_0.4_Tm_0.4/')
if(!file.exists(outdir)) dir.create(outdir, recursive = TRUE)

write.table(expr_dat_3, file = paste0(outdir, 'data.discretized.txt'), col.names = FALSE, row.names = TRUE, quote = FALSE, sep = "\t")

##Step 2 - update ban files with addition priors
#here I use Transcription factor -> target information downloaded from TFLink database
                                    
# Add row and colnames to make life easier
rownames(identity) <- rownames(expr_dat_3)
colnames(identity) <- rownames(expr_dat_3)
identity <- as.data.frame(identity)

for(file in prior_path_new){
  if(file.exists(file)){
    print(paste0("Adding priors from ", basename(file), " to bn.banned.txt"))
    
    # Load up priors file
    prior <- read.table(file, sep = "\t")
    
    # Rename to have consistent column headers
    colnames(prior) <- c("Gene1", "Gene2", "Prior")
    
    # Filter just in case
    prior <- prior[prior$Prior > 0.9,]
    
    # Filter to only gene pairs that both genes are actually in our simplified expression dataset
    prior <- prior[(prior$Gene1 %in% colnames(identity)) & (prior$Gene2 %in% colnames(identity)),]
    
    # Get indices of both lists of genes
    gene1 <- match(prior$Gene1, colnames(identity))
    gene2 <- match(prior$Gene2, colnames(identity))
    
    index_matrix <- unique(cbind(gene2, gene1))
    identity[index_matrix] <- 1
    
    print(paste0(sum(identity > 0) - nrow(identity), " banned gene pairs (", (sum(identity > 0) - nrow(identity))*100/(nrow(identity)*nrow(identity)), "% of pairs)"))
  }
}


# Convert back to a simple matrix
identity <- as.matrix(identity)
rownames(identity) <- NULL
colnames(identity) <- NULL

# And write to file
write.table(identity, paste0(outdir, "bn.banned.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
