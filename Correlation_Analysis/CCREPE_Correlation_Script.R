library("ccrepe")
library("reshape")
library("reshape2")
library("dplyr")
library("knitr")
Rel_abun_asv <- read.table(file="FP004_AMP_exp.txt", sep = "\t", header = T, row.names = 1)

df <- as.data.frame(t(Rel_abun_asv)) 
names(df) <- gsub("^X", "", names(df))
names(df) <- gsub("\\.", "-", names(df))
# calculate the relative abundance
df <- as.data.frame(t(df))
df_relative <- apply(df, 1, function(x) x/sum(x))

Rel_abun_asv <- as.matrix(t(df_relative))

Rel_abun_asv2 <- read.table(file="Si8_FP004_MAB__Genus-table.txt", sep = "\t", header = T, row.names = 1)
#Rel_abun_asv2 <- as.matrix(Rel_abun_asv2)
df2 <- as.data.frame(Rel_abun_asv2) 
names(df2) <- gsub("^X", "", names(df2))
names(df2) <- gsub("\\.", "-", names(df2))
names(df2) <- gsub("MAB", "C", names(df2))
names(df2) <- gsub("_.*$", "", colnames(df2))
df2 <- as.data.frame(t(df2))
df2 <- subset(df2, rownames(df2) != "53-Si8-MAB")
Rel_abun_asv2 <- df2 %>% select(where(~ sum(. > 0) >= 3))
Rel_abun_asv2 <- as.matrix(Rel_abun_asv2)

correlation.table<-ccrepe(x = Rel_abun_asv, y = Rel_abun_asv2, sim.score = cor, sim.score.args = list(method="spearman", use="complete.obs"), min.subj = 6)

r_ccrepe <- correlation.table$sim.score
p_ccrepe <- correlation.table$p.values
q_ccrepe <- correlation.table$q.values
r_ccrepe <- melt(r_ccrepe)
p_ccrepe <- melt(p_ccrepe)
q_ccrepe <- melt (q_ccrepe)

r_p_q_ccrepe_results <- cbind(r_ccrepe,p_ccrepe,q_ccrepe)
r_p_q_ccrepe_results_Final <- r_p_q_ccrepe_results[,c(1:3,6,9)]

colnames(r_p_q_ccrepe_results_Final) <- c("genus1", "genus2", "r_value", "p_value", "q_value")
write.table (r_p_q_ccrepe_results_Final, file = "Combined_corr.txt", sep = "\t", row.names = FALSE, quote = FALSE)
