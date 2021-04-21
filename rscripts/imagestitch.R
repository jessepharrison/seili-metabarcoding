# Seili 16S analysis - jesse harrison 2020-2021
# using seili-r Singularity container (based on seili-r.def)
# additional libpath ####
# (see extra_RPackages.R for extra package installs)

.libPaths(c("~/seili-singularity/rpackages", .libPaths()))

# packages ####

packages <- c("magick")
lapply(packages, require, character.only = TRUE)

# sessioninfo ####

# > sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS

# Matrix products: default
# BLAS/LAPACK: /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so

# locale:
# [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
# [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] magick_2.5.2

# loaded via a namespace (and not attached):
#  [1] Biobase_2.50.0      pkgload_1.1.0       tidyr_1.1.2         jsonlite_1.7.1      splines_4.0.2       foreach_1.5.1       assertthat_0.2.1   
#  [8] BiocManager_1.30.10 stats4_4.0.2        phyloseq_1.34.0     remotes_2.2.0       progress_1.2.2      sessioninfo_1.1.1   pillar_1.4.7       
# [15] lattice_0.20-41     glue_1.4.2          digest_0.6.27       XVector_0.30.0      colorspace_2.0-0    Matrix_1.2-18       plyr_1.8.6         
# [22] pkgconfig_2.0.3     devtools_2.3.2      microbiome_1.12.0   zlibbioc_1.36.0     purrr_0.3.4         scales_1.1.1        processx_3.4.4     
# [29] Rtsne_0.15          tibble_3.0.4        mgcv_1.8-33         generics_0.1.0      IRanges_2.24.0      ggplot2_3.3.2       usethis_1.6.3      
# [36] ellipsis_0.3.1      withr_2.3.0         BiocGenerics_0.36.0 cli_2.2.0           survival_3.2-7      magrittr_2.0.1      crayon_1.3.4       
# [43] memoise_1.1.0       ps_1.4.0            fs_1.5.0            fansi_0.4.1         nlme_3.1-150        MASS_7.3-53         pkgbuild_1.1.0     
# [50] vegan_2.5-7         tools_4.0.2         data.table_1.13.2   prettyunits_1.1.1   hms_0.5.3           lifecycle_0.2.0     stringr_1.4.0      
# [57] Rhdf5lib_1.12.0     S4Vectors_0.28.0    munsell_0.5.0       cluster_2.1.0       callr_3.5.1         Biostrings_2.58.0   ade4_1.7-16        
# [64] compiler_4.0.2      tinytex_0.27        rlang_0.4.9         rhdf5_2.34.0        grid_4.0.2          iterators_1.0.13    rhdf5filters_1.2.0 
# [71] biomformat_1.18.0   rstudioapi_0.13     igraph_1.2.6        testthat_3.0.0      gtable_0.3.0        codetools_0.2-18    multtest_2.46.0    
# [78] reshape2_1.4.4      R6_2.5.0            dplyr_1.0.2         rprojroot_2.0.2     permute_0.9-5       desc_1.2.0          ape_5.4-1          
# [85] stringi_1.5.3       parallel_4.0.2      Rcpp_1.0.5          vctrs_0.3.5         tidyselect_1.1.0    xfun_0.19 

# ggplot2 theme ####

theme_set(theme_classic())

# working directory ####

setwd("~/git/seili-metabarcoding/")

# note: running gc() several times to clear cache ####
# fig 3 ####

# read in the files
fig3a <- image_read('figures/r_output/Fig3a_16S.tiff')
fig3b <- image_read('figures/r_output/Fig3b_18S.tiff')

# combine the elements
fig3 <- c(fig3a, fig3b)
fig3 <- image_append(fig3)

# add some white space
fig3 <- image_border(fig3, color = "white", geometry = "00x120")

# add letters
fig3 <- image_annotate(fig3, text = "a)", 
                       size = 140, color = "black",
                       location = "0+0")

fig3 <- image_annotate(fig3, text = "b)", 
                       size = 140, color = "black",
                       location = "+1772+0")

# crop the bottom white space
fig3 <- image_crop(fig3, "3544x1900")

# export the figure
image_write(fig3, path = "figures/r_output/combined_figs/fig3.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(fig3, fig3a, fig3b)

# clear up the env
gc()

# fig 4 ####

# read in the files
fig4a <- image_read('figures/r_output/Fig4a_16S.tiff')
fig4b <- image_read('figures/r_output/Fig4b_18S.tiff')

# combine the elements
fig4 <- c(fig4a, fig4b)
fig4 <- image_append(fig4, stack = TRUE)

# add some white space
fig4 <- image_border(fig4, color = "white", geometry = "00x120")

# add letters
fig4 <- image_annotate(fig4, text = "a)", 
                       size = 140, color = "black",
                       location = "0+0")

fig4 <- image_annotate(fig4, text = "b)", 
                       size = 140, color = "black",
                       location = "+0+2362")

# clear up the env
gc()

# crop the bottom white space
fig4 <- image_crop(fig4, "3189x4840")

# export the figure
image_write(fig4, path = "figures/r_output/combined_figs/fig4.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(fig4, fig4a, fig4b)

# clear up the env
gc()

# fig 5 ####

# read in the files
fig5a <- image_read('figures/r_output/Fig5a_16S.tiff')
fig5b <- image_read('figures/r_output/Fig5b_18S.tiff')

# add some white space to fig5a for composite centering
fig5a <- image_border(fig5a, color = "white", geometry = "472x00")

# combine the elements
fig5 <- c(fig5a, fig5b)
fig5 <- image_append(fig5, stack = TRUE)

# add some white space to the composite
fig5 <- image_border(fig5, color = "white", geometry = "00x120")

# add letters
fig5 <- image_annotate(fig5, text = "a)", 
                       size = 140, color = "black",
                       location = "0+0")

fig5 <- image_annotate(fig5, text = "b)", 
                       size = 140, color = "black",
                       location = "+0+2362")

# clear up the env
gc()

# crop the bottom white space
fig5 <- image_crop(fig5, "3189x4750")

# export the figure
image_write(fig5, path = "figures/r_output/combined_figs/fig5.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(fig5, fig5a, fig5b)

# clear up the env
gc()

# fig 6 ####

#read in the file
fig6 <- image_read('figures/r_output/Fig6_18S.tiff')

# trim the image
fig6 <- image_trim(fig6)

# export the figure
image_write(fig6, path = "figures/r_output/Fig6_18s.tiff", 
            format = "tiff", quality = 100)

# clear up the env
gc()

# fig 7 ####

# read in the files
fig7a <- image_read('figures/r_output/Fig7a_16S.tiff')
fig7b <- image_read('figures/r_output/Fig7b_18S.tiff')
fig7c <- image_read('figures/r_output/Fig7c_16S.tiff')
fig7d <- image_read('figures/r_output/Fig7d_18S.tiff')

# clear up the env
gc()

# combine top elements
fig7.up <- c(fig7a, fig7b)
fig7.up <- image_append(fig7.up)

# add some white space to top element
fig7.up <- image_border(fig7.up, color = "white", geometry = "00x120")

# add letters to top element
fig7.up <- image_annotate(fig7.up, text = "a)", 
                          size = 140, color = "black",
                          location = "0+0")

fig7.up <- image_annotate(fig7.up, text = "b)", 
                          size = 140, color = "black",
                          location = "+2362+0")

# crop the bottom white space of the top element
fig7.up <- image_crop(fig7.up, "4724x2480")

# combine the bottom elements
fig7.down <- c(fig7c, fig7d)
fig7.down <- image_append(fig7.down)

# add some white space to bottom element
fig7.down <- image_border(fig7.down, color = "white", geometry = "00x120")

# clear up the env
gc()

# add letters to top element
fig7.down <- image_annotate(fig7.down, text = "c)", 
                            size = 140, color = "black",
                            location = "0+0")

fig7.down <- image_annotate(fig7.down, text = "d)", 
                            size = 140, color = "black",
                            location = "+2362+0")

# crop the bottom white space of the bottom element
fig7.down <- image_crop(fig7.down, "4724x2480")

# stack the top + bottom elements

fig7 <- c(fig7.up, fig7.down)
fig7 <- image_append(fig7, stack = TRUE)

# export the figure
image_write(fig7, path = "figures/r_output/combined_figs/fig7.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(fig7, fig7.down, fig7.up, fig7a, fig7b, fig7c, fig7d)

# clear up the env
gc()

# fig S1 ####

# read in the files
figS1a <- image_read('figures/r_output/FigS1a_16S.tiff')
figS1b <- image_read('figures/r_output/FigS1b_18S.tiff')

# combine the elements
figS1 <- c(figS1a, figS1b)
figS1 <- image_append(figS1)

# adjust the white space
figS1 <- image_crop(figS1, "5906x2900")
figS1 <- image_trim(figS1)
figS1 <- image_border(figS1, color = "white", geometry = "00x120")

# add letters
figS1 <- image_annotate(figS1, text = "a)", 
                        size = 140, color = "black",
                        location = "0+0")

figS1 <- image_annotate(figS1, text = "b)", 
                        size = 140, color = "black",
                        location = "+2953+0")

# crop the bottom white space
figS1 <- image_crop(figS1, "5757x2750")

# export the figure
image_write(figS1, path = "figures/r_output/combined_figs/figS1.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(figS1, figS1a, figS1b)

# clear up the env
gc()

# fig S2 ####

# read in the files
figS2a <- image_read('figures/r_output/FigS2a_16S.tiff')
figS2b <- image_read('figures/r_output/FigS2b_18S.tiff')

# combine the elements
figS2 <- c(figS2a, figS2b)
figS2 <- image_append(figS2, stack = TRUE)

# add some white space
figS2 <- image_border(figS2, color = "white", geometry = "00x120")

# add letters
figS2 <- image_annotate(figS2, text = "a)", 
                       size = 140, color = "black",
                       location = "0+0")

figS2 <- image_annotate(figS2, text = "b)", 
                       size = 140, color = "black",
                       location = "+0+2362")

# crop the bottom white space
figS2 <- image_crop(figS2, "3189x4850")

# export the figure
image_write(figS2, path = "figures/r_output/combined_figs/figS2.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(figS2, figS2a, figS2b)

# clear up the env
gc()

# fig S3 ####

# read in the files
figS3a <- image_read('figures/r_output/FigS3a_16S.tiff')
figS3b <- image_read('figures/r_output/FigS3b_16S.tiff')
figS3c <- image_read('figures/r_output/FigS3c_16S.tiff')
figS3d <- image_read('figures/r_output/FigS3d_18S.tiff')
figS3e <- image_read('figures/r_output/FigS3e_18S.tiff')

# combine the elements (upper section)
figS3.up <- c(figS3a, figS3b, figS3c)
figS3.up <- image_append(figS3.up)

# add some white space
figS3.up <- image_border(figS3.up, color = "white", geometry = "00x120")

# add letters
figS3.up <- image_annotate(figS3.up, text = "a)", 
                           size = 140, color = "black",
                           location = "0+0")

figS3.up <- image_annotate(figS3.up, text = "b)", 
                           size = 140, color = "black",
                           location = "+1772+0")

figS3.up <- image_annotate(figS3.up, text = "c)", 
                           size = 140, color = "black",
                           location = "+3544+0")

# crop the bottom white space
figS3.up <- image_crop(figS3.up, "5316x1900")

# combine the elements (lower section)
figS3.down <- c(figS3d, figS3e)
figS3.down <- image_append(figS3.down)

# add some white space
figS3.down <- image_border(figS3.down, color = "white", geometry = "00x120")

# add letters
figS3.down <- image_annotate(figS3.down, text = "d)", 
                           size = 140, color = "black",
                           location = "0+0")

figS3.down <- image_annotate(figS3.down, text = "e)", 
                           size = 140, color = "black",
                           location = "+1772+0")

# crop the bottom white space
figS3.down <- image_crop(figS3.down, "3544x1900")

# clear up the env
gc()

# combine the upper and lower sections
figS3 <- c(figS3.up, figS3.down)
figS3 <- image_append(figS3, stack = TRUE)

# export the figure
image_write(figS3, path = "figures/r_output/combined_figs/figS3.tiff", 
            format = "tiff", quality = 100)

# remove objects
rm(figS3, figS3.down, figS3.up, figS3a, figS3b, figS3c, figS3d, figS3e)

# fig S4 ####

# read in the files
figS4a <- image_read('figures/r_output/FigS4a_16S.tiff')
figS4b <- image_read('figures/r_output/FigS4b_18S.tiff')

# combine the elements
figS4 <- c(figS4a, figS4b)
figS4 <- image_append(figS4, stack = TRUE)

# clear up the env
gc()

# add some white space
figS4 <- image_border(figS4, color = "white", geometry = "00x120")

# clear up the env
gc()

# add letters
figS4 <- image_annotate(figS4, text = "a)", 
                       size = 140, color = "black",
                       location = "0+0")

# clear up the env
gc()

figS4 <- image_annotate(figS4, text = "b)", 
                       size = 140, color = "black",
                       location = "+0+4138")

# clear up the env
gc()

# crop the bottom white space
figS4 <- image_crop(figS4, "5315x8000")

# export the figure
image_write(figS4, path = "figures/r_output/combined_figs/figS4.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(figS4, figS4a, figS4b)

# clear up the env
gc()