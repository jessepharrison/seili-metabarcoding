### PCA for environmental parameters (April 2018) - Jesse Harrison ###
# based on: https://tgmstat.wordpress.com/2013/11/28/computing-and-visualizing-pca-in-r/

# Seili data set

#  --- Clean up R ---

rm(list=ls())

### --- Load data ---

require(vegan)
require(data.table)
require(ggplot2)
require(gridExtra)

envdata <- read.table('/home/jesse/jpharris/micca/Seili/SeiliEnv/SeiliEnv_General_PorCorr.csv', sep = ',', header = TRUE)

# old file without porosity-corrected data:
# envdata <- read.table('/home/jesse/jpharris/micca/Seili/SeiliEnv/SeiliEnv_General.csv', sep = ',', header = TRUE)

theme_set(theme_light())

### --- Split data from site assignment ---

envdata.values <- envdata[, 1:18]
envdata.site <- envdata[, 19] # 19 is the column that gives different sites (we won't use it further here)

### --- PCA ---

envdata.pca <- prcomp(envdata.values,
                       center = TRUE,
                       scale. = TRUE) 

print(envdata.pca)
plot(envdata.pca, type = "l")
summary(envdata.pca) # here it's a good idea to check how much variation (%) is explained by the first two PCs

theta <- seq(0,2*pi,length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle, aes(x,y)) + geom_path()

loadings <- data.frame(envdata.pca$rotation, 
                        .names = row.names(envdata.pca$rotation))
plot <- p + geom_text(data=loadings, 
                        mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
  coord_fixed(ratio=1) +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position="none")
plot

tiff(filename = 'PCA_PorCorr.tiff', 
    width = 1500, height = 1500, res = 300)
plot
dev.off()
