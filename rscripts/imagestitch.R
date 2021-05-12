# haverö metabarcoding study - image processing
# jesse harrison 2020-2021

# additional libpath ####
# (see extra_RPackages.R for extra package installs)

.libPaths(c("~/seili-singularity/rpackages", .libPaths()))

# packages ####

packages <- c("magick")
lapply(packages, require, character.only = TRUE)

# ggplot2 theme ####

theme_set(theme_classic())

# working directory ####

setwd("~/git/seili-metabarcoding/")

# note: running gc() several times to clear cache ####
# fig 2 ####

# read in the file and trim it
fig2 <- image_read('figures/r_output/Fig2.tiff')
fig2 <- image_trim(fig2)

# export the figure
image_write(fig2, path = "figures/r_output/Fig2.tiff", 
            format = "tiff", quality = 100)

# clear up the env
gc()

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
image_write(fig6, path = "figures/r_output/Fig6_18S.tiff", 
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
figS2b <- image_read('figures/r_output/FigS2b_16S.tiff')
figS2c <- image_read('figures/r_output/FigS2c_16S.tiff')
figS2d <- image_read('figures/r_output/FigS2d_18S.tiff')
figS2e <- image_read('figures/r_output/FigS2e_18S.tiff')

# combine the elements (upper section)
figS2.up <- c(figS2a, figS2b, figS2c)
figS2.up <- image_append(figS2.up)

# add some white space
figS2.up <- image_border(figS2.up, color = "white", geometry = "00x120")

# add letters
figS2.up <- image_annotate(figS2.up, text = "a)", 
                           size = 140, color = "black",
                           location = "0+0")

figS2.up <- image_annotate(figS2.up, text = "b)", 
                           size = 140, color = "black",
                           location = "+1772+0")

figS2.up <- image_annotate(figS2.up, text = "c)", 
                           size = 140, color = "black",
                           location = "+3544+0")

# crop the bottom white space
figS2.up <- image_crop(figS2.up, "5316x1900")

# combine the elements (lower section)
figS2.down <- c(figS2d, figS2e)
figS2.down <- image_append(figS2.down)

# add some white space
figS2.down <- image_border(figS2.down, color = "white", geometry = "00x120")

# add letters
figS2.down <- image_annotate(figS2.down, text = "d)", 
                           size = 140, color = "black",
                           location = "0+0")

figS2.down <- image_annotate(figS2.down, text = "e)", 
                           size = 140, color = "black",
                           location = "+1772+0")

# crop the bottom white space
figS2.down <- image_crop(figS2.down, "3544x1900")

# clear up the env
gc()

# combine the upper and lower sections
figS2 <- c(figS2.up, figS2.down)
figS2 <- image_append(figS2, stack = TRUE)

# export the figure
image_write(figS2, path = "figures/r_output/combined_figs/figS2.tiff", 
            format = "tiff", quality = 100)

# remove objects
rm(figS2, figS2.down, figS2.up, figS2a, figS2b, figS2c, figS2d, figS2e)

# fig S3 ####

# read in the files
figS3a <- image_read('figures/r_output/FigS3a_16S.tiff')
figS3b <- image_read('figures/r_output/FigS3b_18S.tiff')

# combine the elements
figS3 <- c(figS3a, figS3b)
figS3 <- image_append(figS3, stack = TRUE)

# clear up the env
gc()

# add some white space
figS3 <- image_border(figS3, color = "white", geometry = "00x120")

# clear up the env
gc()

# add letters
figS3 <- image_annotate(figS3, text = "a)", 
                       size = 140, color = "black",
                       location = "0+0")

# clear up the env
gc()

figS3 <- image_annotate(figS3, text = "b)", 
                       size = 140, color = "black",
                       location = "+0+4138")

# clear up the env
gc()

# crop the bottom white space
figS3 <- image_crop(figS3, "5315x8000")

# export the figure
image_write(figS3, path = "figures/r_output/combined_figs/figS3.tiff", 
            format = "tiff", quality = 100)

# remove files
rm(figS3, figS3a, figS3b)

# clear up the env
gc()