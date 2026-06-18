#library

#Signed
#Plot correlation distribution with different cutoff

# Loop through different RDS files
for (cutoff in seq(10, 80, by = 10)) {
  
  # Construct file name
  rds_file <- paste0("Cluster1/ijw_FDR_1_no_perm_pos_neg_mincell", cutoff, ".rds")
  
  # Load the RDS file
  ijw_FDR <- readRDS(rds_file)
  
  # Create a subset of 'signif.ijw'
  ijw_pos_subset <- ijw_FDR[[1]] 
  ijw_neg_subset <- ijw_FDR[[2]] 
  ijw_both_subset <- rbind(ijw_pos_subset, ijw_neg_subset)

  # Create a dynamic variable name
  var_name <- paste0("ijw_subset_mincell", cutoff)
  
  # Save the subset to the R environment with the dynamic name
  assign(var_name, ijw_both_subset)
  }


# Define a color palette for the histograms
colors <- brewer.pal(n = 8, name = "Set1")  # Adjust the number if needed

# Initialize a list to store individual plots
plot_list <- list()
MSE_list <-  list()
Mean_list <- list()

for (cutoff in seq(10, 80, by = 10)) {
  
  # Create the variable name for the subset
  var_name <- paste0("ijw_subset_mincell", cutoff)
  
  # Access the ijw_subset variable from the environment
  ijw_subset <- get(var_name)
  
  # Calculate the mean and variance of 'rho' in the current ijw_subset
  mean_val <- mean(ijw_subset$rho)
  var_val <- var(ijw_subset$rho)
  sd_val <- sqrt(var_val)  # Standard deviation
  
  # Generate truncated normal distribution using the calculated mean and sd
  x_vals <- seq(min(ijw_subset$rho), max(ijw_subset$rho), length.out = 100)
  truncated_norm <- dtruncnorm(x_vals, a = min(ijw_subset$rho), b = max(ijw_subset$rho),  mean = mean_val, sd = sd_val)
  truncated_norm_df <- data.frame(x_vals = x_vals, expected_val = truncated_norm)
  
  # Calculate the histogram for the current ijw_subset
  hist_data <- hist(ijw_subset$rho, plot = FALSE)

  
  # Interpolate the histogram data to match the x-values of the truncated normal distribution
  observed_density_interpolated <- approx(hist_data$mids, hist_data$density, xout = x_vals)$y
  observed_density_interpolated[is.na(observed_density_interpolated)] <- 0
  

  # Create a data frame for plotting the expected truncated normal distribution
  observed_data <- data.frame(x_vals = x_vals, observed_val = observed_density_interpolated)
  
  # Calculate the mean squared error (MSE) for every value of x
  MSE_df <- merge(observed_data, truncated_norm_df, by = "x_vals")
  MSE_df$MSE_val <- (MSE_df$observed_val - MSE_df$expected_val)^2
  MSE_value <- sum(MSE_df$MSE_val)/nrow(MSE_df)
  
  #store MSE values by appending 
  MSE_list[[as.character(cutoff)]] <- MSE_value
  Mean_list[[as.character(cutoff)]] <- mean_val
  
  # Create a histogram for the current ijw_subset with a line of the truncated normal distribution
  hist_df <- data.frame(mids = hist_data$mids, density = hist_data$density)
  truncated_norm <- data.frame(x_vals = x_vals, density = truncated_norm)
  
  p <- ggplot()+
    geom_histogram(data = hist_df, aes(x = mids, y = density), stat = "identity", binwidth = 0.05, fill = colors[cutoff / 10], color = "black", alpha = 0.7) +
    geom_line(data = truncated_norm, aes(x = x_vals, y = density), 
            color = "red", size = 1) +  # Add the expected line
    labs(title = paste("mincell", cutoff),
         x = "Signed Rho",
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"))+
    annotate("text", x = Inf, y = Inf, label = paste("MSE:", round(MSE_value, 4), "\nMean:", round(mean_val, 4),
      "\nVar:", round(var_val, 4), "\nStd:", round(sd_val, 4)),hjust = 1.1, vjust = 1.1, size = 3, color = "red", fontface = "italic") 
  
  # Add the plot to the list
  plot_list[[as.character(cutoff)]] <- p
}

main_title <- textGrob("C1 PD+Control All gene pairs Avg no 0 Signed", gp = gpar(fontsize = 10, fontface = "bold"))


# Arrange all plots in a 3x4 grid
png("Cluster1/PDControl_C1_FDR_1_signed_correlation_density_MSE.png", width=3000,height=2000,res=300)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 2, top = main_title)
dev.off()

saveRDS(MSE_list, file = "Cluster1/C1_PDControl_avg_signed_MSE_list.rds")
saveRDS(Mean_list, file = "Cluster1/C1_PDControl_avg_signed_Mean_list.rds")


