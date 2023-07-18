#------------------------------------------------------------------
#Boruta on different diets
folder_path <- "/Users/vishnu/WORK/FABIOLA_ANALYSIS/FP_STUDY_12062023/QIIME_ANALYSIS/STATS_ANALYSIS/TAXONOMIC_ANALYSIS/Files_for_LEFSe_and_Boruta/Boruta_Analysis"
file_names <- list.files(folder_path, pattern = "*_diet.txt", full.names = TRUE)
for (file in file_names) {
  data <- read.table(file = file, header = TRUE, row.names = 1, sep = '\t')

  boruta.train <- Boruta(data[,1:(ncol(data)-1)], as.factor(data$diet), pValue = 0.01, mcAdj = TRUE, maxRuns = 100,doTrace = 0, holdHistory = TRUE, getImp = getImpRfZ)

  file_name1 <- paste(file, "Boruta_finaldecision.txt", sep = "_")
  file_name2 <- paste(file, "Boruta_Impotance.txt", sep = "_")
  write.table(boruta.train$finalDecision, file = file_name1, sep = '\t', quote = FALSE)
  write.table(boruta.train$ImpHistory, file = file_name2, sep = '\t', quote = FALSE)
}
#-------------------------------------------------------------------
# Boruta on different genders

folder_path <- "/Users/vishnu/WORK/FABIOLA_ANALYSIS/FP_STUDY_12062023/QIIME_ANALYSIS/STATS_ANALYSIS/TAXONOMIC_ANALYSIS/Files_for_LEFSe_and_Boruta/Boruta_Analysis"
file_names <- list.files(folder_path, pattern = "*_sex.txt", full.names = TRUE)
for (file in file_names) {
  data <- read.table(file = file, header = TRUE, row.names = 1, sep = '\t')
  
  boruta.train <- Boruta(data[,1:(ncol(data)-1)], as.factor(data$sex), pValue = 0.01, mcAdj = TRUE, maxRuns = 100,doTrace = 0, holdHistory = TRUE, getImp = getImpRfZ)
  
  file_name1 <- paste(file, "Boruta_finaldecision.txt", sep = "_")
  file_name2 <- paste(file, "Boruta_Impotance.txt", sep = "_")
  write.table(boruta.train$finalDecision, file = file_name1, sep = '\t', quote = FALSE)
  write.table(boruta.train$ImpHistory, file = file_name2, sep = '\t', quote = FALSE)
}
#----------------------------------------------------------------------
#labdsv on different diets

folder_path <- "/Users/vishnu/WORK/FABIOLA_ANALYSIS/FP_STUDY_12062023/QIIME_ANALYSIS/STATS_ANALYSIS/TAXONOMIC_ANALYSIS/Files_for_LEFSe_and_Boruta/Boruta_Analysis"
file_names <- list.files(folder_path, pattern = "*_diet.txt", full.names = TRUE)
for (file in file_names) {
  famdata <- read.table(file = file, sep = '\t', header = TRUE, row.names = 1)
  iva <- indval(famdata[,1:(ncol(famdata)-1)], famdata$diet)
  gr <- iva$maxcls[iva$pval<=0.01]
  iv <- iva$indcls[iva$pval<=0.01]
  pv <- iva$pval[iva$pval<=0.01]
  fr <- apply(famdata[,1:(ncol(famdata)-1)] >0, 2, sum)[iva$pval<=0.01]
  indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
  indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
  file_name3 <- paste(file, "indvalsummary_groups.txt", sep = "_")
  write.table (indvalsummary, file=file_name3 , sep = "\t")
}

#-----------------------------------------------------------------------
#labdsv on different gender
folder_path <- "/Users/vishnu/WORK/FABIOLA_ANALYSIS/FP_STUDY_12062023/QIIME_ANALYSIS/STATS_ANALYSIS/TAXONOMIC_ANALYSIS/Files_for_LEFSe_and_Boruta/Boruta_Analysis"
file_names <- list.files(folder_path, pattern = "*_sex.txt", full.names = TRUE)
for (file in file_names) {
  famdata <- read.table(file = file, sep = '\t', header = TRUE, row.names = 1)
  iva <- indval(famdata[,1:(ncol(famdata)-1)], famdata$sex)
  gr <- iva$maxcls[iva$pval<=0.01]
  iv <- iva$indcls[iva$pval<=0.01]
  pv <- iva$pval[iva$pval<=0.01]
  fr <- apply(famdata[,1:(ncol(famdata)-1)] >0, 2, sum)[iva$pval<=0.01]
  indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
  indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
  file_name3 <- paste(file, "indvalsummary_groups.txt", sep = "_")
  write.table (indvalsummary, file=file_name3 , sep = "\t")
}
