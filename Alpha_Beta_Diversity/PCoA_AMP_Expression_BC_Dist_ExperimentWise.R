df <- read.table("AMP_exp_data.txt", sep = '\t', header = TRUE, row.names = 1)

# Select the unique values in the column 'experiment'
experiment_names <- unique(df$experiment)
adonis_results <- data.frame()
# Loop through each element of 'experiment_names'
for (experiment_name in experiment_names) {
  # Subset the dataframe if the element is present in the 'experiment' column
  subset_df <- df[df$experiment == experiment_name, ]
  # Transpose the dataframe
  KO_test <- t(subset_df)
  
  # Subset rows 7:12 of the transposed dataframe
  KO_test_subset <- as.data.frame(KO_test[7:12, ])
  
  # Convert the dataframe into numeric
  df_numeric <- KO_test_subset
  for (i in 1:ncol(KO_test_subset)) {
     df_numeric[, i] <- as.numeric(KO_test_subset[, i])
  }
  KO <- df_numeric
  
  KO_test_T <- as.data.frame(t(KO_test))
  KO_test_1 <- as.data.frame(KO_test_T[7:12])
  rownames(KO) <- colnames(KO_test_1)
  colnames(KO) <- row.names(KO_test_1)
  
  #KO_proportions <- KO/colSums(KO)[col(KO)]
  #KO_proportions1 <- as.data.frame(t(KO_proportions))
  
  # use these lines only when you need to plot the PCoA from raw data without calculating relative abundance
  #-------------------------------------------------
  df_numeric <- KO_test_1
  for (i in 1:ncol(KO_test_1)) {
    df_numeric[, i] <- as.numeric(KO_test_1[, i])
  }
  KO_proportions1 <- df_numeric
  #-------------------------------------------------
  
  covar <- as.vector(colnames(KO_test_T[1:6]))
  # Elements to remove
  elements_to_remove <- c("mouse.number", "sample")
  
  # Remove specific elements
  covar_upd <- covar[!(covar %in% elements_to_remove)]
  for (i in 1:length(covar_upd)) {
    element <- covar_upd[i]
    if (element == "diet") {
      Colors <- c("darkolivegreen4", "red", "orange")
      Colors1 <- c("darkolivegreen","red4", "orangered")
    } else if (element == "experiment"){
      Colors <- c("mediumpurple3", "snow4", "aquamarine")
      Colors1 <- c("violet", "gray30", "deepskyblue")
    }else if (element == "Mouse_provider"){
      Colors <- c("goldenrod1", "salmon4")
      Colors1 <- c("goldenrod4", "salmon")
    }else if (element == "sex"){
      Colors <- c("lightblue", "gray")
      Colors1 <- c("midnightblue", "darkgray")
    }else if (element == "diet_sex"){
      Colors <- c("darkolivegreen", "darkolivegreen4", "red", "red4")
      Colors1 <- c("darkolivegreen4", "darkolivegreen", "red4", "red")
    }
    
      # Check if unique elements in the vector are greater than or equal to 2
    if (length(unique(KO_test_T[[element]])) >= 2) {
        # Print the vector
      print(KO_test_T[[element]])
      class <- KO_test_T[[element]]
      Bray_pcoa <-pcoa(vegdist(KO_proportions1, "bray"))
      Bray_pcoa$values[1:2,]
      mds.var.per = round(Bray_pcoa$values$Eigenvalues/sum(Bray_pcoa$values$Eigenvalues)*100, 1)
      Bray_PCoA_MATRIX <- Bray_pcoa$vectors[,1:2]
      Bray_PCoA_MATRIX <- as.data.frame(Bray_PCoA_MATRIX)
        
      Bray_distances <-vegdist(KO_proportions1, "bray")
      adonis2(Bray_distances ~ as.factor(KO_test_T[[element]]))
      file_name1 <- paste(experiment_name, element, "Bray_Curtis_Distance.txt", sep = "_")
        
      write.table(as.matrix(Bray_distances), file = file_name1, quote = FALSE, sep = '\t')
        
      Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, class)
      file_name2 <- paste(experiment_name, element, "PCA_data_proportions.txt", sep = "_")
      write.table(Bray_PCoA_MATRIX_New, file = file_name2, quote = FALSE, sep = '\t')
      BC_adonis_results <- data.frame()
      out <- as.data.frame(adonis2(Bray_distances ~ as.factor(KO_test_T[[element]])))
      out1 <- bind_cols(out, status = element)
      out1 <- bind_cols(out1, Dist = "BC")
      out2 <- bind_cols(out1, Experiment = experiment_name)
      BC_adonis_results  <- bind_rows(out2, BC_adonis_results)
      adonis_results <- bind_rows(adonis_results, BC_adonis_results)
      # Print the combined data frame
      print(BC_adonis_results )
      pc <- c(1,2)
      file_name3 <- paste(experiment_name, element, "Bray-Curtis.jpg", sep = "_")
      jpeg(file_name3, height = 10, width = 10, units = 'in', res = 600)
      
      plot(Bray_pcoa$vectors[,1:2], bg=Colors[as.factor(Bray_PCoA_MATRIX_New$class)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
      ordiellipse(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 70, col = Colors)
      #ordispider(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, lty=3, spider ="centroid", lwd=1, col="black")
      legend("topright", legend = sort(unique(Bray_PCoA_MATRIX_New$class)), col = Colors,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
      text(Bray_pcoa$vectors[,1:2], labels=as.factor(rownames(KO_proportions1)), cex=0.6, font=1, pos=1)
      #abline(h=0, v=0, col = "gray60")
      dev.off ()
    }
  }
} 
write.table (adonis_results, file = "Adonis_Test_Results.txt", sep = "\t", col.names = NA, row.names = TRUE)
