# Load required packages
library(tidyverse)
library(lme4)
library(broom)
library(patchwork)

# Read the dataset
df <- read.csv("lean_Ob_table.csv")

# Set 'diet_duration' as a factor with '6w' as the reference level
df$diet_duration <- factor(df$diet_duration, levels = c("Lean", "Ob"))

# Identify gene expression columns
gene_cols <- grep("^(Si1|Si5|Si8)_", names(df), value = TRUE)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = all_of(gene_cols),
    names_to = c("cell_type", "gene_name"),
    names_sep = "_",
    values_to = "expression_value"
  )

# Re-affirm factor levels
df_long$diet_duration <- factor(df_long$diet_duration, levels = c("Lean", "Ob"))

# Define metabolic traits
metabolic_vars <- c("Bodyweight", "Insulin", "Glucose", "HOMA_IR")

# Loop through each cell type
for (ct in unique(df_long$cell_type)) {
  df_cell <- df_long %>% filter(cell_type == ct)
  
  models <- list()
  stats_combined <- list()
  
  for (met in metabolic_vars) {
    all_fixed_plots <- list()
    all_slope_plots <- list()
    
    for (gene in unique(df_cell$gene_name)) {
      df_gene <- df_cell %>%
        filter(gene_name == gene) %>%
        mutate(expression_log = log10(expression_value + 1))
      
      # Fit the model
      model_formula <- as.formula(paste(met, "~ expression_log + diet_duration"))
      mod <- lm(model_formula, data = df_gene)
      models[[paste(gene, met, sep = "_")]] <- mod
      
      # Extract fixed effects using broom
      fixed_effects <- tidy(mod, conf.int = TRUE) %>%
        mutate(
          gene = gene,
          metabolic_trait = met,
          cell_type = ct
        )
      stats_combined[[paste(gene, met, sep = "_")]] <- fixed_effects
      
      ## --- FIXED EFFECTS PLOT ---
      p_fixed <- ggplot(fixed_effects, aes(x = estimate, y = term)) +
        geom_point(color = "blue") +
        geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        labs(
          title = paste("Fixed Effects:", gene, "-", ct),
          x = "Estimate ± 95% CI", y = NULL
        ) +
        theme_minimal()
      all_fixed_plots[[gene]] <- p_fixed
      
      ## --- SLOPE PLOT FOR EACH DIET DURATION ---
      new_expr_log <- seq(min(df_gene$expression_log), max(df_gene$expression_log), length.out = 100)
      pred_data <- expand.grid(
        expression_log = new_expr_log,
        diet_duration = levels(df_gene$diet_duration)
      )
      pred_data[[met]] <- predict(mod, newdata = pred_data)
      
      df_points <- df_gene %>%
        select(expression_log, diet_duration, !!sym(met))
      
      gene_stats <- fixed_effects %>%
        filter(term == "expression_log") %>%
        mutate(
          label = paste0(
            "Gene effect: estimate = ", round(estimate, 2),
            #", SE = ", round(std.error, 2),
            #", t = ", round(statistic, 2),
            #", 95% CI = [", round(conf.low, 2), ", ", round(conf.high, 2), "]",
            ", p = ", format.pval(p.value, digits = 2)
          )
        ) %>%
        pull(label)
      
      wrapped_stats <- strwrap(gene_stats, width = 50)
      gene_stats <- paste(wrapped_stats, collapse = "\n")
      
      p_slope <- ggplot(pred_data, aes(x = expression_log, y = .data[[met]], color = diet_duration)) +
        geom_line(size = 1) +
        geom_point(
          data = df_points,
          aes(x = expression_log, y = .data[[met]], fill = diet_duration),
          shape = 21, color = "black", alpha = 0.6, size = 2, stroke = 0.4
        ) +
        labs(
          title = paste("Slope Plot:", gene, "-", ct),
          subtitle = gene_stats,
          x = paste("log10(", gene, " Expression + 1)", sep = ""),
          y = met
        ) +
        theme_minimal() +
        theme(plot.subtitle = element_text(size = 9, color = "gray30"))
      
      all_slope_plots[[gene]] <- p_slope
    }
    
    ## --- SAVE PLOTS FOR THIS CELL TYPE & TRAIT ---
    pdf(paste0("Fixed_Effects_", ct, "_", met, ".pdf"), width = 10, height = 8)
    print(wrap_plots(all_fixed_plots) + plot_annotation(title = paste("Fixed Effects for", met, "in", ct)))
    dev.off()
    
    pdf(paste0("Slope_Plots_", ct, "_", met, ".pdf"), width = 10, height = 8)
    print(wrap_plots(all_slope_plots) + plot_annotation(title = paste("Slope Plots for", met, "in", ct)))
    dev.off()
  }
  
  ## --- SAVE STATS FOR THIS CELL TYPE ---
  all_stats_df <- bind_rows(stats_combined)
  write_tsv(all_stats_df, paste0("Model_Stats_", ct, ".tsv"))
}
