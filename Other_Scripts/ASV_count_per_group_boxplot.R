# Load the necessary library
library(ggplot2)
df <- read.table("boxplot_input2.txt", sep="\t", row.names = 1, header=T) 

#OTU_ID	SUM	Count	Mouse_provider	experiment	source	sample	sex	diet_sex	diet
#121-Si8-C_S25	2334	7	Taconic	Si8_FP002	Content	121_Si8_C	Male	WSD_Male	WSD
#121-Si8-MAB_S59	17520	29	Taconic	Si8_FP002	MAB	121_Si8_MAB	Male	WSD_Male	WSD
#122-Si8-MAB_S60	27480	140	Taconic	Si8_FP002	MAB	122_Si8_MAB	Male	WSD_Male	WSD
#123-Si8-C_S27	876	8	Taconic	Si8_FP002	Content	123_Si8_C	Male	WSD_Male	WSD
#123-Si8-MAB_S61	14863	91	Taconic	Si8_FP002	MAB	123_Si8_MAB	Male	WSD_Male	WSD


# Define the multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Open a PDF device
pdf("boxplot_output2.pdf")

# Convert 'SUM' column to numeric
df$SUM <- as.numeric(df$SUM)

# Convert categorical columns to factors
df$Mouse_provider <- as.factor(df$Mouse_provider)
df$experiment <- as.factor(df$experiment)
df$source <- as.factor(df$source)
df$sex <- as.factor(df$sex)
df$diet_sex <- as.factor(df$diet_sex)
df$diet <- as.factor(df$diet)

# Define the factor columns
factor_columns <- c("Mouse_provider", "experiment", "source", "sex", "diet_sex", "diet")

# Define colors for the boxes
box_colors <- c("red", "blue", "green", "orange", "purple", "yellow")

# Create a list to store the plots
plot_list <- list()

# Loop through factor columns and create boxplots
for (factor_col in factor_columns) {
  p <- ggplot(df, aes(x = !!as.symbol(factor_col), y = SUM, fill = !!as.symbol(factor_col))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(title = paste("Boxplot of SUM by", factor_col)) +
    scale_fill_manual(values = box_colors) +
    theme(legend.position = "top")
  
  # Add the plot to the list
  plot_list[[factor_col]] <- p
}

# Arrange all the plots in a single page
multiplot(plotlist = plot_list, cols = 3)

# Close the PDF device
dev.off()
