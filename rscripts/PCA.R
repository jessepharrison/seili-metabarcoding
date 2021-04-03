# Seili metabarcoding study - PCA for environmental variables
# jesse harrison 2020-2021
# using seili-r Singularity container (based on seili-r.def)

# additional libpath ####
# (see extra_RPackages.R for extra package installs)

.libPaths(c("/home/jharriso/seili-singularity/rpackages", .libPaths()))

# packages ####

packages <- c("vegan", "data.table", "ggplot2", 
              "Cairo", "ggrepel")

lapply(packages, require, character.only = TRUE)

# sessioninfo ####

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS

# Matrix products: default
# BLAS/LAPACK: /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so

# locale:
#  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8   
#  [6] LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C        
# [11] LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

# attached base packages:
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] multcompView_0.1-8   Cairo_1.5-12.2       RColorBrewer_1.1-2   knitr_1.30           microbiome_1.12.0    RVAideMemoire_0.9-78
#  [7] cowplot_1.1.0        plyr_1.8.6           data.table_1.13.2    gridExtra_2.3        vegan_2.5-7          lattice_0.20-41     
# [13] permute_0.9-5        ggplot2_3.3.2        phyloseq_1.34.0     

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.5          ape_5.4-1           tidyr_1.1.2         prettyunits_1.1.1   Biostrings_2.58.0   digest_0.6.27      
#  [7] foreach_1.5.1       R6_2.5.0            stats4_4.0.2        pillar_1.4.7        zlibbioc_1.36.0     rlang_0.4.9        
# [13] progress_1.2.2      rstudioapi_0.13     S4Vectors_0.28.0    Matrix_1.2-18       labeling_0.4.2      splines_4.0.2      
# [19] Rtsne_0.15          stringr_1.4.0       igraph_1.2.6        munsell_0.5.0       tinytex_0.27        compiler_4.0.2     
# [25] xfun_0.19           pkgconfig_2.0.3     BiocGenerics_0.36.0 multtest_2.46.0     mgcv_1.8-33         biomformat_1.18.0  
# [31] tidyselect_1.1.0    tibble_3.0.4        IRanges_2.24.0      codetools_0.2-18    crayon_1.3.4        dplyr_1.0.2        
# [37] withr_2.3.0         MASS_7.3-53         rhdf5filters_1.2.0  nlme_3.1-150        jsonlite_1.7.1      gtable_0.3.0       
# [43] lifecycle_0.2.0     magrittr_2.0.1      scales_1.1.1        stringi_1.5.3       farver_2.0.3        XVector_0.30.0     
# [49] reshape2_1.4.4      ellipsis_0.3.1      generics_0.1.0      vctrs_0.3.5         Rhdf5lib_1.12.0     iterators_1.0.13   
# [55] tools_4.0.2         ade4_1.7-16         Biobase_2.50.0      glue_1.4.2          purrr_0.3.4         hms_0.5.3          
# [61] parallel_4.0.2      survival_3.2-7      colorspace_2.0-0    rhdf5_2.34.0        cluster_2.1.0

# ggplot2 theme ####

theme_set(theme_classic())

# working directory ####

setwd("/home/jharriso/git/seili-metabarcoding/")

# load porosity-corrected environmental data ####

envdata <- read.table('envdata/SeiliEnv_General_PorCorr.csv', sep = ',', header = TRUE)

# split data from site assignments ####

# sites are listed column 19 (we won't use it further here)
envdata.values <- envdata[, 1:18]
envdata.site <- envdata[, 19]

# PCA (centered and scaled) ####

envdata.pca <- prcomp(envdata.values,
                       center = TRUE,
                       scale. = TRUE) 

# PCA loadings and line plot (commented out, kept for completeness)
# print(envdata.pca)
# plot(envdata.pca, type = "l")

# PCA summary (for checking % variation explained by first two PCs) ####

pca.summary <- summary(envdata.pca)

cat(" -----------------------------------",
    "\n", "POROSITY-CORRECTED ENV DATA - PCA SUMMARY", "\n",
    "------------------------------------",
    "\n",
    file = "stats/pca_envdata.txt")

capture.output(pca.summary, 
               file = "stats/pca_envdata.txt",
               append = TRUE)

# convert PCA loadings with variable names into a data frame ####

pca.loadings <- data.frame(envdata.pca$rotation,
                           names = row.names(envdata.pca$rotation))

# plot the PCA loadings ####

Cairo(file = "figures/r_output/Fig2.png", 
      type = "png", 
      units = "cm", 
      width = 20, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg = "white")

pca.plot <- ggplot(pca.loadings, 
                   aes(x = PC1,
                       y = PC2)) +
  geom_point(color = "red") +
  geom_text_repel(aes(label = names),
                  color = "gray45",
                  min.segment.length = 0, 
                  seed = 42, 
                  box.padding = 0.6,
                  size = 5) + 
  scale_colour_manual(values = "black.palette") +
  coord_fixed(ratio = 1) +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5) +
  theme(legend.title = element_blank()) +
  theme(text = element_text(size = 20)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, 
                                                    r = 5, 
                                                    b = 0, 
                                                    l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, 
                                                    r = 0, 
                                                    b = 0, 
                                                    l = 0)))

pca.plot
dev.off()