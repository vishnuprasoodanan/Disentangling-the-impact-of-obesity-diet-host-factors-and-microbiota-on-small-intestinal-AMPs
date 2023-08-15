folder_path <- "/Users/vishnu/WORK/FABIOLA_ANALYSIS/FP_STUDY_12062023/QIIME_ANALYSIS/STATS_ANALYSIS/ALPHA_BETA_DIVERSITY/Intersample_distance/INTER_SAMPLE_DISTANCE_FILES"
library(dplyr)
# Get the list of file names ending with "xyz.txt"
file_names <- list.files(folder_path, pattern = "*.txt_gender", full.names = TRUE)
Colors <- c("lightblue", "gray")
Colors1 <- c("midnightblue", "darkgray")
plot_lst <- vector("list", 0)
pdf(file = "Abundance_boxplots.pdf")
# Read each file iteratively
for (file in file_names) {
  data <- read.table(file = file, sep = "\t", header = T)
  selected_rows <- data %>% filter((Gender1 == "Male" & Gender2 == "Male") | (Gender1 == "Female" & Gender2 == "Female"))
  subset_data <- selected_rows[, c(4, 3)]
  updated_variable <- sub(folder_path, "", file)
  updated_variable <- sub("/", "", updated_variable)
  updated_variable <- sub("_distances.txt_gender", "", updated_variable)
  P<- ggplot(subset_data, aes(Gender1, Value, fill=Gender1))+
    ggtitle(updated_variable)+
    labs(y = "inter-sample distance")+
    geom_boxplot(outlier.shape=NA)+ 
    geom_point(aes_string(colour = factor(subset_data$Gender1)), position=position_jitterdodge(jitter.width = 0.5))+
    scale_color_manual(values = Colors1)+
    scale_fill_manual(values = Colors)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_lst[[file]] <- P
}
ml <- marrangeGrob(plot_lst, nrow = 2, ncol = 2)
print(ml)
dev.off()

#----------------------------------------------------------------------------

folder_path <- "/Users/vishnu/WORK/FABIOLA_ANALYSIS/FP_STUDY_12062023/QIIME_ANALYSIS/STATS_ANALYSIS/ALPHA_BETA_DIVERSITY/Intersample_distance/INTER_SAMPLE_DISTANCE_FILES"
library(dplyr)
# Get the list of file names ending with "xyz.txt"
file_names <- list.files(folder_path, pattern = "*.txt_diet", full.names = TRUE)

Colors <- c("darkolivegreen4", "red") #colors used for different diets
Colors1 <- c("darkolivegreen","red4") #colors used for different diets

plot_lst <- vector("list", 0)
pdf(file = "Abundance_boxplots.pdf")
# Read each file iteratively
for (file in file_names) {
  data <- read.table(file = file, sep = "\t", header = T)
  selected_rows <- data %>% filter((Diet1 == "Chow" & Diet2 == "Chow") | (Diet1 == "WSD" & Diet2 == "WSD")  | (Diet1 == "HFD" & Diet2 == "HFD"))
  subset_data <- selected_rows[, c(4, 3)]
  updated_variable <- sub(folder_path, "", file)
  updated_variable <- sub("/", "", updated_variable)
  updated_variable <- sub("_distances.txt_diet", "", updated_variable)
  P<- ggplot(subset_data, aes(Diet1, Value, fill=Diet1))+
    ggtitle(updated_variable)+
    labs(y = "inter-sample distance")+
    geom_boxplot(outlier.shape=NA)+ 
    geom_point(aes_string(colour = factor(subset_data$Diet1)), position=position_jitterdodge(jitter.width = 0.5))+
    scale_color_manual(values = Colors1)+
    scale_fill_manual(values = Colors)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_lst[[file]] <- P
}
ml <- marrangeGrob(plot_lst, nrow = 2, ncol = 2)
print(ml)
dev.off()
