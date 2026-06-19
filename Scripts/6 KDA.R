##download KDA package from here: https://github.com/mw201608/mnml-public/blob/master/pkgs/KDA_0.2.2.zip 
###install.packages("KDA_0.2.2.tar.gz", repos=NULL, type = "source")

##note run the KDA for L1 (undirected), L2 (undirected), and L3 (directed) then combine them to get the final list

directed <- FALSE ##choose FALSE for L1, L2, and TRUE for L3
layer <- 2 # choose 1,2, or 3
minDsCut <- -1
fgeneinfo <- NULL

library( class )
library( cluster )
library( rpart )
#library( sma ) # this is needed for plot.mat below
library( lattice ) # require is design for use inside functions 

# read in network

cnet <- as.matrix(cnet)
dim(cnet)

totalnodes <- union(cnet[,1] , cnet[,2] )
totalnodes <- as.data.frame(totalnodes)
colnames(totalnodes) <- "Geneid"

#read in gene lists
listMatrix <- totalnodes

listMatrix <- merge(listMatrix, result_padj_0.05)
listMatrix <- listMatrix %>% select(1,2)
dim( listMatrix )
listMatrix <- as.matrix( listMatrix )
listMatrix5[1:2,]
ncols <- dim( listMatrix)[2]

modules <- names(table(listMatrix[,ncols]))

xkdrMatrix <- NULL
paraMatrix <- NULL

# process each gene list

for ( em in modules )
{

	print( paste( "*****************" , em , "********************" ) )
	
	esel <- listMatrix[,ncols] == em
	
	genes <- union( listMatrix[esel,1] , NULL )
	
	if(layer >=1 )
	{
		# expand network by K-hop nearest neighbors layers
		expandNet <- findNLayerNeighborsLinkPairs( linkpairs = cnet , subnetNodes = genes ,
				nlayers = layer , directed = FALSE )
	} 
	else
	{
		# no expansion
		expandNet <- getSubnetworkLinkPairs( linkpairs = cnet , subnetNodes = genes )
	}
	dim( expandNet )
	
	allnodes <- union( expandNet[,1] , expandNet[,2] )

	################################################################################################
# keydriver for a given network
#
	
  if (directed)
  {
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes, background=NULL, directed = directed ,
			         nlayers = 6 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
			         boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }
  else
  {
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes , directed = directed ,
			         nlayers = 2 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
				boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }
	
	if ( is.null( ret ) )
	{
		next
	}
	
	fkd <- ret[[1]]
	parameters <- ret[[2]]
	
	fkd2<- cbind( rep( em , dim( fkd )[1] ) , fkd )
	xkdrMatrix <- rbind( xkdrMatrix , fkd2 )
	paraMatrix<- rbind( paraMatrix , parameters)
	
}


colnames(xkdrMatrix)[2] <- "Geneid"
xkdrMatrix <- merge(xkdrMatrix, annotation, by = "Geneid")


cluster_KDA_L1 <- xkdrMatrix_cluster
cluster_KDA_L1 <- merge(cluster_KDA_L1, result_padj_0.05, by = "Geneid", all = TRUE)
cluster_KDA_L1$Cell_type <- "OPC"
cluster_KDA_L1 <- mutate(cluster_KDA_L1, "log10Padj_network" = (-log10(as.numeric(cluster_KDA_L1$pvalue_corrected_subnet))))
cluster_KDA_L1 <- cluster_KDA_L1[complete.cases(cluster_KDA_L1$log10Padj_network),]
cluster_KDA_L1 <- merge(cluster_KDA_L1, annotation, by = "Geneid")
cluster_KDA_L1$KDA <- "KDA"


#write.csv(cluster_KDA_L3, "cluster_KDA_L3.csv")
#write.csv(cluster_KDA_L1, "cluster_KDA_L1.csv")
#write.csv(cluster_KDA_L2, "cluster_KDA_L2.csv")

cluster_KDA_L1$KDA <- "KDA_L1"
cluster_KDA_L2$KDA <- "KDA_L2"
cluster_KDA_L1_L2_int <- merge(cluster_KDA_L1, cluster_KDA_L2, by = "Symbol")
cluster_KDA_L1_L2_L3 <- merge(cluster_KDA_L1_L2_int, cluster_KDA_L3, by = "Symbol", all = TRUE)
cluster_KDA_L1_L2_L3_clean <- dplyr::select(cluster_KDA_L1_L2_L3, "Symbol", "L1_padj_network" = "pvalue_corrected_subnet.x", "L1_FE" = "fold_change_subnet.x", "L1_log10_padj_network" = "log10Padj_network.x", "KDA_L1" = "KDA.x", "L2_padj_network" = "pvalue_corrected_subnet.y", "L2_FE" = "fold_change_subnet.y", "L2_log10_padj_network" = "log10Padj_network.y", "KDA_L2" = "KDA.y", "L3_FE" = "fold_change_subnet", "L3_padj_network" = "pvalue_corrected_subnet", "L3_log10_padj_network" = "log10Padj_network", "KDA_L3" = "KDA")
cluster_KDA_L1_L2_L3_sig <- merge(cluster_KDA_L1_L2_L3_clean, result5_padj_0.05, by = "Symbol")
cluster_KDA_L1_L2_L3_final <- cluster_KDA_L1_L2_L3_clean 

cluster_KDA_L1_L2_L3_final$maxlog10_network <- apply(cluster_KDA_L1_L2_L3_final[, c("L1_log10_padj_network", "L2_log10_padj_network", "L3_log10_padj_network")], 1, function(row) {
  max(ifelse(is.na(row), 0, row))
})

# Replace NA with 0 and find the column (A, B, or C) for max
cluster_KDA_L1_L2_L3_final$KDAmax <- apply(cluster_KDA_L1_L2_L3_final[, c("L1_log10_padj_network", "L2_log10_padj_network", "L3_log10_padj_network")], 1, function(row) {
  names(which.max(ifelse(is.na(row), 0, row)))
})


cluster_KDA_L1_L2_L3_final$FE <- apply(cluster_KDA_L1_L2_L3_final[, c("L1_FE", "L2_FE", "L3_FE")], 1, function(row) {
  max(ifelse(is.na(row), 0, row))
})

cluster_KDA_L1_L2_L3_final$Cell_type <- "Oli"
cluster_KDA_L1_L2_L3_final <- merge(cluster_KDA_L1_L2_L3_final, result_padj_0.05, by = "Symbol", all = TRUE)
cluster_KDA_L1_L2_L3_final <- cluster_KDA_L1_L2_L3_final[complete.cases(cluster_KDA_L1_L2_L3_final$maxlog10_network),]


write.csv(cluster_KDA_L1_L2_L3_final, "cluster_0.4_KDA_L1_L2_L3_final_02.csv")
