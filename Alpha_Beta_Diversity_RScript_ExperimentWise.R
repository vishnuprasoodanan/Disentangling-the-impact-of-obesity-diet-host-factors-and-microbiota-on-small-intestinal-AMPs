#This code in written by Vishnu Prasoodanan P K
#The input data for this 
library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(phyloseq)
library (ape)
library (ggplot2)
library(vegan)

#--------------------------------------------------------
Colors <- c("darkolivegreen4", "red") #colors used for different diets
#Colors <- c("lightblue", "gray") #colors used for gender

Colors1 <- c("darkolivegreen","red4") #colors used for different diets
#Colors1 <- c("midnightblue", "darkgray") #colors used for gender
#-------------------------------------------------------

otu_table_in <- read.csv("Selected_FeatureTable_rarefied_edited.txt", sep = "\t")
otu_table_t <- setNames(data.frame(t(otu_table_in[,-1])), otu_table_in[,1])
otu_table_t <- tibble::rownames_to_column(otu_table_t)

names(otu_table_t)[names(otu_table_t) == "rowname"] <- "SampleID"

otu_table_t$SampleID <- gsub('^[X]', '', otu_table_t$SampleID) 
otu_table_t$SampleID <- gsub('[.]', '-', otu_table_t$SampleID)

metadata_all <- read.table("selected_Metadata.tsv", sep="\t", row.names = 1, header=T)
metadata_all <- tibble::rownames_to_column(metadata_all)
names(metadata_all)[names(metadata_all) == "rowname"] <- "SampleID"

merge_otu_table_t <- merge(otu_table_t, metadata_all, by.x = "SampleID")
merge_otu_table_t$experiment <- as.factor(merge_otu_table_t$experiment)
merge_otu_table_t$source <- as.factor(merge_otu_table_t$source)
# Get the factor column name
factor_column <- "experiment"

# Iterate through each factor variable
factor_levels <- levels(merge_otu_table_t[[factor_column]])
for (level in factor_levels) {
  # Filter the data frame for the current factor level
  Si8_FP001 <- merge_otu_table_t[merge_otu_table_t[[factor_column]] == level, ]
  numeric_columns <- Si8_FP001[, sapply(Si8_FP001, is.numeric)]
  numeric_columns <- numeric_columns[, colSums(numeric_columns) > 0]
  
  non_numeric_columns <- Si8_FP001[, !sapply(Si8_FP001, is.numeric)]
  
  # Create a new data frame
  filtered_df <- data.frame(numeric_columns, non_numeric_columns)
  row.names(filtered_df) <- Si8_FP001$SampleID
  
  Si8_FP001_Final1 <- as.data.frame(t(filtered_df[,1:(ncol(filtered_df)-8)]))
  Si8_FP001_Metadata1 <- filtered_df[,(ncol(filtered_df)-7):ncol(filtered_df)]
  Si8_FP001_Metadata <- rownames_to_column(Si8_FP001_Metadata1, var = "RowID")
  colnames(Si8_FP001_Metadata)[1] <- "SampleID"
  Si8_FP001_Final <- rownames_to_column(Si8_FP001_Final1, var = "RowID")
  colnames(Si8_FP001_Final)[1] <- "SampleID"
  file_name1 <- paste(level, "_feature-table.txt", sep = "_")
  file_name2 <- paste(level, "_Metadata.txt", sep = "_")
  write.table(Si8_FP001_Final, file = file_name1, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(Si8_FP001_Metadata, file = file_name2, sep = "\t", quote = FALSE, row.names = FALSE)
  filtered_df$source <- as.factor(filtered_df$source)
  factor_column2 <- "source"
  factor_levels2 <- levels(filtered_df[[factor_column2]])
  for (level2 in factor_levels2) {
    gut <- filtered_df[filtered_df[[factor_column2]] == level2, ]
    numeric_columns2 <- gut[, sapply(gut, is.numeric)]
    numeric_columns2 <- numeric_columns2[, colSums(numeric_columns2) > 0]
    
    non_numeric_columns2 <- gut[, !sapply(gut, is.numeric)]
    # Create a new data frame
    filtered_df2 <- data.frame(numeric_columns2, non_numeric_columns2)
    row.names(filtered_df2) <- gut$SampleID
    gut_Final1 <- as.data.frame(t(filtered_df2[,1:(ncol(filtered_df2)-8)]))
    gut_Metadata1 <- filtered_df2[,(ncol(filtered_df2)-7):ncol(filtered_df2)]
    gut_Metadata <- rownames_to_column(gut_Metadata1, var = "RowID")
    colnames(gut_Metadata)[1] <- "SampleID"
    gut_Final <- rownames_to_column(gut_Final1, var = "RowID")
    colnames(gut_Final)[1] <- "SampleID"
    file_name3 <- paste(level, level2, "_feature-table.txt", sep = "_")
    file_name4 <- paste(level, level2, "_Metadata.txt", sep = "_")
    write.table(gut_Final, file = file_name3, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(gut_Metadata, file = file_name4, sep = "\t", quote = FALSE, row.names = FALSE)
    
    otu_table_final <- gut_Final
    row.names(otu_table_final) <- otu_table_final$SampleID
    otu_table_final <- otu_table_final[, -1]
    metadata <- gut_Metadata
    row.names(metadata) <- metadata$SampleID
    metadata <- metadata[, -1]
    taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
    taxonomy <- as.matrix(taxonomy)
    # Read in metadata
    # Read in tree
    phy_tree <- read_tree("tree.nwk")
    # Import all as phyloseq objects
    OTU <- otu_table(otu_table_final, taxa_are_rows = TRUE)
    TAX <- tax_table(taxonomy)
    META <- sample_data(metadata)
    # Sanity checks for consistent OTU names
    head(taxa_names(TAX))
    head(taxa_names(OTU))
    head(taxa_names(phy_tree))
    # Same sample names
    sample_names(OTU) <- gsub('[.]', '-', sample_names(OTU))
    sample_names(OTU) <- gsub('^[X]', '', sample_names(OTU)) 
    head(sample_names(OTU))
    head(sample_names(META))
    file_name5 <- paste(level, level2, "_alpha_diversity_estimates_diet.pdf", sep = "_")
    file_name6 <- paste(level, level2, "_alpha_diversity_estimates_diet.jpg", sep = "_")
    # Finally merge to create Phyloseq object!
    ps <- phyloseq(OTU, TAX, META, phy_tree)
    #mat <- t(otu_table(otu_table_final, taxa_are_rows = TRUE))
    #class(mat) <- "matrix"
    #raremax <- min(rowSums(mat))
    #system.time(rarecurve(mat, step = 100, sample = raremax, col = "blue", label = TRUE ))
    #ps <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 711)

    p <- plot_richness(ps, "diet", measures = c("Observed", "Chao1", "Shannon"), color = NULL, shape = NULL)
    p <- p + geom_boxplot(data = p$data, aes(x= diet, y = value, fill = diet)) + 
      labs(x="",y="Alpha Diversity Measure") + 
      theme_bw() +
      scale_color_manual(values = Colors1)+
      scale_fill_manual(values = Colors)+
      geom_point(aes(colour = factor(diet)), position=position_jitterdodge(jitter.width = 1), size = 5) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24)) 
    ggsave(file_name5, p, width = 10, height = 10, limitsize = FALSE)
    jpeg(file_name6, height = 10, width = 15, units = 'in', res = 600)
    p
    dev.off ()
    #---------------------------PCoA Analysis
    
    Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )
    
    #---------------------- weighted UniFrac PCoA
    file_name7 <- paste(level, level2, "W-UniFrac_distances.txt", sep = "_")
    UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
    UniFrac_distances <- as.matrix(UniFrac_distances)
    UniFrac_dist_column <- melt(UniFrac_distances)
    write.table (UniFrac_dist_column, file = file_name7, sep = "\t")
    
    Uni_pcoa <- pcoa(UniFrac_distances)
    Uni_pcoa$values[1:2,]
    mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
    pc <- c(1,2)
    file_name8 <- paste(level, level2, "W-Unifrac_PCOA.jpg", sep = "_")
    jpeg(file_name8, height = 10, width = 10, units = 'in', res = 600)
    plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen","red4")[as.factor(Bushman2@sam_data$diet)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
    text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
    ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("darkolivegreen","red4"))
    ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, lty=3, spider ="centroid", lwd=1, col="black")
    legend("bottomleft", legend = c("Chow", "WSD"), col = c("darkolivegreen","red4"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
    abline(h=0, v=0, col = "gray60")
    dev.off ()
    adonis2(UniFrac_distances ~ Bushman2@sam_data$diet)
    
    #--------------------- unweighted UniFrac PCoA
    
    UWUniFrac_distances <- UniFrac(Bushman2, weighted=FALSE)
    UWUniFrac_distances <- as.matrix(UWUniFrac_distances)
    UniFrac_dist_column <- melt(UWUniFrac_distances)
    file_name9 <- paste(level, level2, "unweighted_UWUniFrac_distances.txt", sep = "_")
    write.table (UniFrac_dist_column, file = file_name9, sep = "\t")
    Uni_pcoa <- pcoa(UWUniFrac_distances)
    Uni_pcoa$values[1:2,]
    mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
    pc <- c(1,2)
    
    file_name10 <- paste(level, level2, "unweighted_UWUniFrac_distances.jpg", sep = "_")
    jpeg(file_name10, height = 10, width = 10, units = 'in', res = 600)
    plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen","red4")[as.factor(Bushman2@sam_data$diet)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
    text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
    ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("darkolivegreen","red4"))
    ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, lty=3, spider ="centroid", lwd=1, col="black")
    legend("bottomleft", legend = c("Chow", "WSD"), col = c("darkolivegreen","red4"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
    abline(h=0, v=0, col = "gray60")
    dev.off ()
    adonis2(UWUniFrac_distances ~ Bushman2@sam_data$diet)
    
    #-------------------------- Bray-Curtis PCoA
    BC_distances <- distance(Bushman2, method="bray")
    BC_distances <- as.matrix(BC_distances)
    UniFrac_dist_column <- melt(BC_distances)
    
    file_name11 <- paste(level, level2, "BC_distances.txt", sep = "_")
    write.table (UniFrac_dist_column, file = file_name11, sep = "\t")
    Uni_pcoa <- pcoa(BC_distances)
    Uni_pcoa$values[1:2,]
    mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
    pc <- c(1,2)
    
    file_name11 <- paste(level, level2, "BC_distances.jpg", sep = "_")
    jpeg(file_name11, height = 10, width = 10, units = 'in', res = 600)
    plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen","red4")[as.factor(Bushman2@sam_data$diet)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
    text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
    ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("darkolivegreen","red4"))
    ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, lty=3, spider ="centroid", lwd=1, col="black")
    legend("topleft", legend = c("Chow", "WSD"), col = c("darkolivegreen","red4"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
    abline(h=0, v=0, col = "gray60")
    dev.off ()
    adonis2(BC_distances ~ Bushman2@sam_data$diet)
    
  }
}
