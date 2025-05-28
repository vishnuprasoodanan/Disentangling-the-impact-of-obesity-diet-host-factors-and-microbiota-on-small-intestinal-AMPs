# Load required libraries
library(tidyverse)
library(ggpubr)
library(robustbase)
library(compositions)
library(stringr)

# Function to clean taxonomic labels
clean_taxa_names <- function(name) {
  if (grepl(";g__", name)) {
    return(str_split(name, ";g__")[[1]][2])
  } else if (grepl(";f__", name)) {
    return(str_split(name, ";f__")[[1]][2])
  } else if (grepl(";o__", name)) {
    return(str_split(name, ";o__")[[1]][2])
  } else if (grepl(";c__", name)) {
    return(str_split(name, ";c__")[[1]][2])
  } else if (grepl(";p__", name)) {
    return(str_split(name, ";p__")[[1]][2])
  } else if (grepl(";d__", name)) {
    return(str_split(name, ";d__")[[1]][2])
  } else {
    return(name)
  }
}

# Read and clean microbiome data
mb_data <- read.delim("mb_data_content.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
new_names <- sapply(rownames(mb_data), clean_taxa_names)
rownames(mb_data) <- make.unique(new_names)

# Read AMP expression and metadata
amp_data <- read.delim("amp_data_content.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata <- read.delim("metadata_all.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure sample consistency
all_samples <- intersect(colnames(mb_data), colnames(amp_data))
mb_data <- mb_data[, all_samples]
amp_data <- amp_data[, all_samples]
metadata <- metadata[metadata$sample_id %in% all_samples, ]
rownames(metadata) <- metadata$sample_id

# CLR transformation of microbiome data (add pseudocount 1)
clr_mb_data <- clr(mb_data + 1)

# Helper to wrap axis labels to two lines if > 30 characters
split_label <- function(label) {
  if (nchar(label) > 30) return(str_wrap(label, width = 30)) else return(label)
}

# Function to remove up to 2 Mahalanobis outliers
remove_outliers <- function(df, sample_names) {
  if (nrow(df) <= 3) return(list(data = df, outliers = character(0)))
  tryCatch({
    cov_robust <- covMcd(df, alpha = 0.75)
    d2 <- mahalanobis(df, cov_robust$center, cov_robust$cov)
    cutoff <- quantile(d2, probs = 0.975)
    outlier_idx <- which(d2 > cutoff)
    if (length(outlier_idx) > 2) outlier_idx <- head(order(d2, decreasing = TRUE), 2)
    list(data = df[-outlier_idx, ], outliers = sample_names[outlier_idx])
  }, error = function(e) list(data = df, outliers = character(0)))
}

# Loop over metadata grouping columns
grouping_cols <- setdiff(names(metadata), "sample_id")

for (group_col in grouping_cols) {
  unique_values <- unique(metadata[[group_col]])
  
  for (val in unique_values) {
    selected_samples <- rownames(metadata[metadata[[group_col]] == val, ])
    mb_sub <- mb_data[, selected_samples, drop = FALSE]
    amp_sub <- amp_data[, selected_samples, drop = FALSE]
    clr_mb_sub <- clr_mb_data[, selected_samples, drop = FALSE]
    
    mb_t <- as.data.frame(t(mb_sub))
    clr_mb_t <- as.data.frame(t(clr_mb_sub))
    amp_t <- as.data.frame(t(amp_sub))
    
    # Loop for each transformation type
    for (type in c("RAW", "LOG", "CLR")) {
      plots <- list()
      outlier_log <- c()
      plot_count <- 0
      
      pdf_file <- paste0("ScatterPlots_", type, "_", group_col, "_", val, ".pdf")
      pdf(pdf_file, width = 10, height = 8)
      
      for (genus in rownames(mb_data)) {
        for (amp in rownames(amp_data)) {
          if ((type == "CLR" && genus %in% colnames(clr_mb_t) || genus %in% colnames(mb_t)) &&
              amp %in% colnames(amp_t)) {
            
            if (type == "RAW") {
              df <- data.frame(Genus = mb_t[[genus]], AMP = amp_t[[amp]])
            } else if (type == "LOG") {
              df <- data.frame(Genus = mb_t[[genus]], AMP = log1p(amp_t[[amp]]))
            } else if (type == "CLR") {
              df <- data.frame(Genus = clr_mb_t[[genus]], AMP = log1p(amp_t[[amp]]))
            }
            
            df <- df[complete.cases(df), ]
            sample_names <- rownames(df)
            
            if (nrow(df) >= 4) {
              res <- remove_outliers(df, sample_names)
              df_clean <- res$data
              outliers <- res$outliers
              
              if (length(outliers)) {
                outlier_log <- c(outlier_log,
                                 paste(genus, amp, ":", paste(outliers, collapse = ", "))
                )
              }
              
              xlab_txt <- split_label(ifelse(type == "CLR", paste("CLR(Abundance):", genus), paste("Abundance:", genus)))
              ylab_txt <- split_label(ifelse(type == "RAW",
                                             paste("Expression:", amp),
                                             paste("log(AMP + 1):", amp)))
              
              # Publication-quality plot
              p <- ggplot(df_clean, aes(x = Genus, y = AMP)) +
                geom_point(size = 3, shape = 21, fill = "steelblue", color = "black", alpha = 0.5) +
                geom_smooth(method = "lm", color = "gray40", size = 0.8, se = FALSE) +
                ggtitle(paste(type, "- Genus:", genus, "vs AMP:", amp)) +
                xlab(xlab_txt) + ylab(ylab_txt) +
                theme_minimal(base_size = 12) +
                theme(
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
                  axis.title = element_text(face = "bold", size = 10),
                  axis.text = element_text(size = 9)
                )
              
              plots[[length(plots) + 1]] <- p
              plot_count <- plot_count + 1
              
              # Print 2x2 plots per page
              if (plot_count %% 4 == 0 || (genus == tail(rownames(mb_data), 1) && amp == tail(rownames(amp_data), 1))) {
                print(ggarrange(plotlist = plots, ncol = 2, nrow = 2))
                plots <- list()
              }
            }
          }
        }
      }
      
      dev.off()
      writeLines(outlier_log, paste0("Outliers_", type, "_", group_col, "_", val, ".txt"))
    }
  }
}
