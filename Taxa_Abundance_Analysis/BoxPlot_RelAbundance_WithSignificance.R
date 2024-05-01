#This Script is created by Vishnu Prasoodanan P K
#It will take Genus-level feature table (i.e, ASV count at phylum level) and metadata file. Then it will convert the count values to relative abundance. 
#Then separate the relative abundance data of each experiment, then separate the relative abundance data of each source (MAB and content).
#Then it will plot stacked barplot for each relative abundance data for experiment-source combinations

library(tidyr)
library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(phyloseq)
library (ape)
library (ggplot2)
library(vegan)
library(RColorBrewer)
library(viridis)
library(colorspace)

p_data <- read.table(file = "feature-table.txt", sep = "\t", header = T, row.names = 1)
p_data_filtered <- p_data[rowSums(p_data) >= 1, ]
#p_data_filtered <- p_data_filtered[colSums(p_data_filtered) != 0, ]
p_data_filtered  <- p_data_filtered [, colSums(p_data_filtered ) != 0]

df <- p_data_filtered
names(df) <- gsub("^X", "", names(df))
names(df) <- gsub("\\.", "-", names(df))
# calculate the relative abundance
df <- as.data.frame(t(df))
df_relative <- t(apply(df, 1, function(x) x/sum(x)))

# transpose the result back to its original orientation
df_relative <- as.data.frame(t(df_relative))
write.table(df_relative, file = "sel_feature_table_abundance.txt", sep = '\t', quote = FALSE, row.names = TRUE)

#----------------------------------------------------------------------------------------
#experimwnt-wise separation of data

#otu_table_in <- read.csv("sel_feature_table_abundance.txt", sep = "\t")
otu_table_in <- df_relative
rownames(otu_table_in) <- gsub("[^[:alnum:]_]", "_", rownames(otu_table_in))
otu_table_t <- as.data.frame(t(otu_table_in ))

#otu_table_t <- setNames(data.frame(t(otu_table_in[,-1])), otu_table_in[,1])
otu_table_t <- tibble::rownames_to_column(otu_table_t)

names(otu_table_t)[names(otu_table_t) == "rowname"] <- "SampleID"

otu_table_t$SampleID <- gsub('^[X]', '', otu_table_t$SampleID) 
otu_table_t$SampleID <- gsub('[.]', '-', otu_table_t$SampleID)

metadata_all <- read.table("Metadata_corrected_analysis.txt", sep="\t", row.names = 1, header=T)
metadata_all <- tibble::rownames_to_column(metadata_all)
names(metadata_all)[names(metadata_all) == "rowname"] <- "SampleID"

merge_otu_table_t <- merge(otu_table_t, metadata_all, by.x = "SampleID")
merge_otu_table_t$experiment <- as.factor(merge_otu_table_t$experiment)
merge_otu_table_t$source <- as.factor(merge_otu_table_t$source)
# Get the factor column name
factor_column <- "experiment"

# Iterate through each factor variable in experiment
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
  # Iterate through each factor variable in source (content and MAB)
  for (level2 in factor_levels2) {
    gut <- filtered_df[filtered_df[[factor_column2]] == level2, ]
    numeric_columns2 <- gut[, sapply(gut, is.numeric)]
    numeric_columns2 <- numeric_columns2[, colSums(numeric_columns2) > 0]
    
    non_numeric_columns2 <- gut[, !sapply(gut, is.numeric)]
    # Create a new data frame
    filtered_df2 <- data.frame(numeric_columns2, non_numeric_columns2)
    row.names(filtered_df2) <- gut$SampleID
    gut_Final1 <- as.data.frame(t(filtered_df2[,1:(ncol(filtered_df2)-9)]))
    gut_Metadata1 <- filtered_df2[,(ncol(filtered_df2)-7):ncol(filtered_df2)]
    gut_Metadata <- rownames_to_column(gut_Metadata1, var = "RowID")
    colnames(gut_Metadata)[1] <- "SampleID"
    gut_Final <- rownames_to_column(gut_Final1, var = "RowID")
    colnames(gut_Final)[1] <- "SampleID"
    file_name3 <- paste(level, level2, "_feature-table.txt", sep = "_")
    file_name4 <- paste(level, level2, "_Metadata.txt", sep = "_")
    write.table(gut_Final, file = file_name3, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(gut_Metadata, file = file_name4, sep = "\t", quote = FALSE, row.names = FALSE)
    genus_rel <-as.data.frame(t(gut_Final))
    # Store the first row as the header
    header <- unname(as.character(genus_rel[1, ]))
    # Remove the first row from the dataframe
    genus_rel <- genus_rel[-1, ]
    # Set the column names of the dataframe using the stored header
    colnames(genus_rel) <- header
    genus_rel$SampleID <- rownames(genus_rel)
    rownames(gut_Metadata) <- gut_Metadata$SampleID
    # Set the "SampleID" column as the first column
    genus_rel <- genus_rel[, c(ncol(genus_rel), 1:(ncol(genus_rel)-1))]
    # Compare row names
    if(identical(rownames(genus_rel), rownames(gut_Metadata))) {
      print("Matching")
    } else {
      print("Not matching")
    }
    boxplot_input <- merge(genus_rel, gut_Metadata[, "diet_sex", drop = FALSE], by = "row.names")
    boxplot_input$Row.names <- NULL
    boxplot_input$diet_sex <- factor(boxplot_input$diet_sex)
    rownames(boxplot_input) <- boxplot_input$SampleID
    
    # Move the diet_sex column to the second position
    boxplot_input <- boxplot_input %>%
      select(SampleID, diet_sex, everything())
    boxplot_input[, -c(1, 2)] <- apply(boxplot_input[, -c(1, 2)], 2, as.numeric)
    
    # Print the modified dataframe
    #print(boxplot_input)
    # Extract levels of the factor column diet_sex
    lev_group <- levels(boxplot_input$diet_sex)
    print(kruskalmc(boxplot_input$d__Bacteria_p__Firmicutes_c__Clostridia_o__Lachnospirales_f__Lachnospiraceae_g__Roseburia ~ boxplot_input$diet_sex, data = boxplot_input, probs = 0.05, alpha = 0.05))
    ####boxplot script
    #print(lev_group)
    file_name2 <- paste(level, level2, "_boxplot.pdf", sep = "_")
    pdf(file = file_name2)
    colnames(boxplot_input[,1:ncol(boxplot_input)]) -> Species_name
    
    Colors <- c("darkgreen", "darkolivegreen4", "red", "red4","#1F77B4", "#FF7F0E","#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4")
    Colors1 <- c("darkolivegreen4", "darkgreen", "red4","red", "midnightblue", "orange", "green", "salmon4", "purple", "blue")
    
    plot_lst <- vector("list", length = ncol(boxplot_input))
    for (i in 3:ncol(boxplot_input)) {    
      
      species = Species_name[i]
      data1 = boxplot_input[,c(1:2,i)]
      colnames(data1) <- c("ID", "Status", "name")
      
      print(species)
      species <- gsub("_", " ", species)
      
      #CHECK SPECIFIC TAXA
      if (species == "d Bacteria p Firmicutes c Clostridia o Peptostreptococcales Tissierellales f Peptostreptococcaceae g Romboutsia" && level == "Si8_FP001" && level2 == "Content")
      {
        print(kruskalmc(name ~ Status, data = data1, probs = 0.05, alpha = 0.05))
      }
      # P<- ggplot(data1, aes(Status, name, fill=Status))+
      #   ggtitle(species)+
      #   labs(y = "Relative-Abundance")+
      #   geom_boxplot(outlier.shape=NA)+ 
      #   geom_point(aes(colour = factor(Status)), position=position_jitterdodge(jitter.width = 0.5))+
      #   scale_color_manual(values = Colors1)+
      #   theme(plot.title = element_text(size = 8, face = "bold")) +
      #   scale_fill_manual(values = Colors)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
      # plot_lst[[i]] <- P
      
      # Check the number of unique values in the 'Status' column
      unique_values <- unique(data1$Status)
      num_unique <- length(unique_values)
      
      # Initialize a vector to store significance stars
      significance <- rep("", length(data1$Status))
      
      # Perform the appropriate test based on the number of unique values
      if (num_unique == 2) {
        # Wilcoxon test
        p_val <- wilcox.test(name ~ Status, data = data1, paired = FALSE, p.adjust.method = "fdr", alternative = "two.sided", conf.level = 0.95)$p.value
        # Check significance
        if (p_val < 0.05) {
          significant_levels <- levels(data1$Status)
          # Combine levels into pairs
          significant_pairs <- combn(significant_levels, 2, simplify = FALSE)
          }
        } else if (num_unique > 2) {
        # Kruskal-Wallis test
        kruskal_test_df <- as.data.frame(kruskalmc(name ~ Status, data = data1, probs = 0.05, alpha = 0.05))
        kruskal_test_df$Comparisons <- rownames(kruskal_test_df)

        
        # Move the new column to the first position
        kruskal_test_df <- kruskal_test_df[, c(ncol(kruskal_test_df), 1:(ncol(kruskal_test_df)-1))]
        significant_pairs <- kruskal_test_df %>%
          filter(dif.com.stat.signif) %>%
          select(Comparisons)
        
        # Extract significant pairs
        significant_pairs <- strsplit(significant_pairs$Comparisons, "-")
        significant_pairs <- lapply(significant_pairs, trimws)
        #print(significant_pairs)
      }
      p <- ggplot(data1, aes(Status, name, fill = Status)) +
        ggtitle(species) +
        labs(y = "Relative-Abundance") +
        geom_boxplot(outlier.shape = NA) + 
        #geom_jitter(width = 0.2, size = 4, alpha = 0.6) +
        #geom_point(aes(colour = factor(Status)), size = 4, alpha = 0.6) +
        geom_point(aes(colour = factor(Status)), position = position_jitterdodge(jitter.width = 0.5),  size = 4, alpha = 0.6) +
        scale_color_manual(values = Colors1) +
        scale_fill_manual(values = Colors) +
        theme(plot.title = element_text(size = 8, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      # Add stars between significant pairs
      for (pair in significant_pairs) {
        #print(pair)
        p <- p + geom_signif(comparisons = list(pair), map_signif_level = TRUE, textsize = 6)
        }
      plot_lst[[i]] <- p
    }
    ml <- marrangeGrob(plot_lst, nrow = 1, ncol = 1)
    print(ml)
    dev.off()
  }
}
