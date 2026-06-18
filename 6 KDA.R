install.packages("KDA_0.2.2.tar.gz", repos=NULL, type = "source")

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
	
	esel_c5 <- listMatrix[,ncols] == em
	
	genes_c5 <- union( listMatrix[esel,1] , NULL )
	
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
	dim( expandNet_c5 )
	
	allnodes_c5 <- union( expandNet_c5[,1] , expandNet_c5[,2] )

	################################################################################################
# keydriver for a given network
#
	
  if (directed)
  {
    ret_c5 <- keydriverInSubnetwork( linkpairs = expandNet_c5 , signature = genes_c5, background=NULL, directed = directed ,
			         nlayers = 6 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
			         boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }
  else
  {
    ret_c5 <- keydriverInSubnetwork( linkpairs = expandNet_c5 , signature = genes_c5 , directed = directed ,
			         nlayers = 2 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
				boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }
	
	if ( is.null( ret_c5 ) )
	{
		next
	}
	
	fkd_c5 <- ret_c5[[1]]
	parameters_c5 <- ret_c5[[2]]
	
	fkd2_c5 <- cbind( rep( em , dim( fkd_c5 )[1] ) , fkd_c5 )
	xkdrMatrix_c5 <- rbind( xkdrMatrix_c5 , fkd2_c5 )
	paraMatrix_c5 <- rbind( paraMatrix_c5 , parameters_c5)
	
}


colnames(xkdrMatrix_c5)[2] <- "Geneid"
xkdrMatrix_c5 <- merge(xkdrMatrix_c5, annotation, by = "Geneid")

c5_KDA_L1 <- xkdrMatrix_c5
c5_KDA_L1 <- merge(c5_KDA_L1, result5_padj_0.05, by = "Geneid", all = TRUE)
c5_KDA_L1$Cell_type <- "OPC"
c5_KDA_L1 <- mutate(c5_KDA_L1, "log10Padj_network" = (-log10(as.numeric(c5_KDA_L1$pvalue_corrected_subnet))))
c5_KDA_L1 <- c5_KDA_L1[complete.cases(c5_KDA_L1$log10Padj_network),]
c5_KDA_L1 <- merge(c5_KDA_L1, annotation, by = "Geneid")
c5_KDA_L1$KDA <- "KDA"


#write.csv(c6_KDA_L3, "c5_KDA_L3.csv")
#write.csv(c3_KDA_L1, "c3_KDA_L1.csv")
#write.csv(c3_KDA_L2, "c3_KDA_L2.csv")

c5_KDA_L1$KDA <- "KDA_L1"
c5_KDA_L2$KDA <- "KDA_L2"
c5_KDA_L1_L2_int <- merge(c5_KDA_L1, c5_KDA_L2, by = "Symbol")
c5_KDA_L1_L2_L3 <- merge(c5_KDA_L1_L2_int, c5_KDA_L3, by = "Symbol", all = TRUE)
c5_KDA_L1_L2_L3_clean <- dplyr::select(c5_KDA_L1_L2_L3, "Symbol", "L1_padj_network" = "pvalue_corrected_subnet.x", "L1_FE" = "fold_change_subnet.x", "L1_log10_padj_network" = "log10Padj_network.x", "KDA_L1" = "KDA.x", "L2_padj_network" = "pvalue_corrected_subnet.y", "L2_FE" = "fold_change_subnet.y", "L2_log10_padj_network" = "log10Padj_network.y", "KDA_L2" = "KDA.y", "L3_FE" = "fold_change_subnet", "L3_padj_network" = "pvalue_corrected_subnet", "L3_log10_padj_network" = "log10Padj_network", "KDA_L3" = "KDA")
c5_KDA_L1_L2_L3_sig <- merge(c5_KDA_L1_L2_L3_clean, result5_padj_0.05, by = "Symbol")
c5_KDA_L1_L2_L3_final <- c5_KDA_L1_L2_L3_clean 

c5_KDA_L1_L2_L3_final$maxlog10_network <- apply(c5_KDA_L1_L2_L3_final[, c("L1_log10_padj_network", "L2_log10_padj_network", "L3_log10_padj_network")], 1, function(row) {
  max(ifelse(is.na(row), 0, row))
})

# Replace NA with 0 and find the column (A, B, or C) for max
c5_KDA_L1_L2_L3_final$KDAmax <- apply(c5_KDA_L1_L2_L3_final[, c("L1_log10_padj_network", "L2_log10_padj_network", "L3_log10_padj_network")], 1, function(row) {
  names(which.max(ifelse(is.na(row), 0, row)))
})


c5_KDA_L1_L2_L3_final$FE <- apply(c5_KDA_L1_L2_L3_final[, c("L1_FE", "L2_FE", "L3_FE")], 1, function(row) {
  max(ifelse(is.na(row), 0, row))
})

c5_KDA_L1_L2_L3_final$Cell_type <- "OPC"
c5_KDA_L1_L2_L3_final <- merge(c5_KDA_L1_L2_L3_final, result5_padj_0.05, by = "Symbol", all = TRUE)
c5_KDA_L1_L2_L3_final <- c5_KDA_L1_L2_L3_final[complete.cases(c5_KDA_L1_L2_L3_final$maxlog10_network),]


write.csv(c5_KDA_L1_L2_L3_final, "c5_0.4_KDA_L1_L2_L3_final_02.csv")
