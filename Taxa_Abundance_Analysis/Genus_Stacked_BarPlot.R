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

# Calculate the row-wise mean abundance
row_means <- rowMeans(otu_table_in)

# Create a dataframe with rownames and their corresponding row means
mean_df <- data.frame(Genera = rownames(otu_table_in), Mean_Abundance = row_means)

# Sort the dataframe by Mean_Abundance in descending order
mean_df <- mean_df[order(-mean_df$Mean_Abundance), ]

# Print the top 50 rownames based on the mean abundance
selected_51 <- c(mean_df$Genera[1:42], "Others")
write.table(selected_51, file = "top_genera.txt", sep = '/t', quote = FALSE)
# Generate a palette of 215 colors from the viridis color map
colors_g <- c("rosybrown3", "gold4", "blue1", "blueviolet", "brown1", "burlywood1", "yellowgreen", 
                     "chocolate1", "red", "cornflowerblue", "orange", "cyan", 
                     "darkblue", "darkcyan","indianred4", "sienna1", "darkgreen", "darkkhaki", 
                     "darkmagenta", "tomato2", "plum1","firebrick3", "midnightblue",
                     "darkred", "darkseagreen1", "darkslateblue", "lightcoral",
                     "darkslategray1", "darkturquoise","mediumorchid1", "deeppink1", 
                     "deepskyblue1", "palevioletred1", "dodgerblue1", "peru", 
                     "forestgreen", "deeppink4", "lemonchiffon4", "lightgreen", "lightcyan4",
                     "lightgoldenrod4","lightpink4","gray")

# Map the genera to colors
mapped_colors <- setNames(colors_g, selected_51)
#names(colors_g) = levels(selected_51)


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
    
    gut_Final_v1 <- as.data.frame(t(gut_Final1))
    gut_Final_v1 <- tibble::rownames_to_column( gut_Final_v1)
    names(gut_Final_v1)[names(gut_Final_v1) == "rowname"] <- "SampleID"
    gut_Metadata_v1 <- gut_Metadata
    gut_Metadata_v2 <- gut_Metadata_v1[, -c( 2:6,9)]
    gut_Final_v2 <- merge(gut_Final_v1, gut_Metadata_v2, by.x = "SampleID")
    
    gut_Final_diet <- gut_Final_v2 %>% select(SampleID, diet_sex, everything())%>% select(-sex)
    row.names(gut_Final_diet) <- gut_Final_diet$SampleID
    gut_Final_diet1 <- gut_Final_diet[, -1]
    combined <- paste(gut_Final_diet1$diet_sex, rownames(gut_Final_diet1), sep = "_")
    row.names(gut_Final_diet1) <- combined
    gut_Final_diet1$diet_sex <- NULL
    gut_Final_diet_t <- as.data.frame(t(gut_Final_diet1))
    gut_Final_diet_t <- tibble::rownames_to_column(gut_Final_diet_t)
    names(gut_Final_diet_t)[names(gut_Final_diet_t) == "rowname"] <- "Name"
    # Calculate row sums excluding first column
    row_sums <- rowSums(gut_Final_diet_t[, -1])
    
    # Find the indices of the rows with the 5 highest row sums
    top_rows <- order(row_sums, decreasing = TRUE)[1:20]
    
    # Select the rows with the 5 highest row sums
    df_top <- gut_Final_diet_t[top_rows, ]
    
    # Calculate column sums
    column_sums <- 1-(colSums(df_top[, -1]))
    
    # Create a new row with column sums and label "Others"
    sum_row <- data.frame(Name = "Others", t(column_sums))
    names(sum_row) <- gsub("\\.", "-", names(sum_row))
    # Add the new row as the last row of the dataframe
    df_with_sum_row <- rbind(df_top, sum_row)
    
    df <- df_with_sum_row 
    data_order <- df[order(apply(df, 1, min)),]
    # convert data from wide to long format using tidyr
    df_long <- pivot_longer(df, cols = -Name, names_to = "Sample", values_to = "Abundance")
    # create stacked barplot using ggplot2
    
    header1 <- paste(level, level2, "_Genus_relative_abundance", sep = "_")
    file_name5 <- paste(level, level2, "_Genus_abundance_diet.pdf", sep = "_")
    p <- ggplot(df_long, aes(x = Sample, y = Abundance, fill = Name)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = mapped_colors) +  # Use the color mapping
      labs(title = header1, x = "Sample", y = "Relative Abundance") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggsave(file_name5, p, width = 25, height = 10, limitsize = FALSE)
    
    data <- gut_Final_diet
    rownames(data) <- NULL
    # Compute average abundance by group
    group_averages<- data %>% group_by(diet_sex) %>% dplyr::summarize(across(starts_with("d__"), mean))
    
    file_name6 <- paste(level, level2, "_average_group_abundance.txt", sep = "_")
    write.table(group_averages, file = file_name6, sep = '\t', quote = FALSE, row.names = FALSE)
    
    group_averages_cs <- colSums(group_averages[, -1])
    
    # Find the names of the columns with the 20 highest column sums
    group_averages_top_columns <- names(group_averages_cs)[order(group_averages_cs, decreasing = TRUE)[1:20]]
    
    # Select the columns with the 20 highest column sums
    group_averages_top_20 <- group_averages[, c("diet_sex", group_averages_top_columns)]
    group_averages_top_21 <- group_averages_top_20 %>% mutate(Others = 1-(rowSums(select(., -1))))
    
    #group_averages_top_21 <- group_averages_top_20 %>% mutate(Others = 1-(rowSums(.)))
    group_averages_long <- pivot_longer(group_averages_top_21, cols = -diet_sex, names_to = "SampleID", values_to = "Abundance")
    # create stacked barplot using ggplot2
    file_name7 <- paste(level, level2, "_average_group_abundance.pdf", sep = "_")
    header2 <- paste(level, level2, "_Genus_group_average_relative_abundance", sep = "_")
    q <- ggplot(group_averages_long, aes(x = diet_sex, y = Abundance, fill = SampleID)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = mapped_colors) +
      labs(title = header2, x = "diet", y = "Average Relative Abundance") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(file_name7, q, width = 20, height = 10, limitsize = FALSE)
  }
