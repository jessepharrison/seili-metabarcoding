# haverö metabarcoding study - PCA for environmental variables
# jesse harrison 2020-2022

# make sure your working directory is the github repository root!
# (i.e. path/to/seili-metabarcoding)

# add singularity folder to libpaths

.libPaths(c("singularity", .libPaths()))

# packages ####

packages <- c("vegan", "data.table", "ggplot2", 
              "Cairo", "ggrepel", "colorspace")

lapply(packages, require, character.only = TRUE)

# ggplot2 theme ####

theme_set(theme_classic())

# disable scientific notation
options(scipen=10000)

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


# add a column with with model groupings ####

groups <- c("Included in Model 1", 
            "No model", 
            "No model", 
            "No model", 
            "Included in Model 1", 
            "No model", 
            "Included in Model 2", 
            "Included in Model 2", 
            "No model", 
            "Included in Model 1", 
            "No model", 
            "Included in Model 1", 
            "No model", 
            "No model", 
            "No model", 
            "Included in Model 1", 
            "Included in Model 2", 
            "Included in Model 2")

pca.loadings$group <- as.factor(groups)

# plot the PCA loadings ####

# Create the plot

Cairo(file = "figures/r_output/Fig2.tiff", 
      type = "tiff", 
      units = "cm", 
      width = 35, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg = "white")

scale_colour_discrete <- scale_colour_viridis_d

pca.plot <- ggplot(pca.loadings, 
                   aes(x = PC1,
                       y = PC2,
                       fill = names,
                       colour = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = names),
                  color = "gray45",
                  min.segment.length = 0, 
                  seed = 42, 
                  box.padding = 0.6,
                  size = 5) + 
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
                                                    l = 0))) +
  scale_fill_discrete(name = "names",
                      labels = c(expression(bold("C_1cm:") ~ "sediment C"[{"org"}]* " content, top 1 cm (% dry mass)"),
                                 expression(bold("CN_1cm:") ~ "sediment C"[{"org"}]* ":N"[{"tot"}]* " ratio, top 1 cm (mol:mol)"),
                                 expression(bold("Depth:") ~ "water column depth (m)"),
                                 expression(bold("Farm:") ~ "distance to the farm (km)"),
                                 expression(bold("Grain:") ~ "sediment median grain size (μm)"),
                                 expression(bold("HS_Inv:") ~ "porewater H"[2]* "S inventory (μmol cm"^{"-2"}*")"),
                                 expression(bold("NH4_Inv:") ~ "porewater NH"[4]^{"+"}* " inventory (μmol cm"^{"-2"}*")"),
                                 expression(bold("NO2_Inv:") ~ "porewater NO"[2]^{"-"}* " inventory (μmol cm"^{"-2"}*")"),
                                 expression(bold("NOX_BW:") ~ "bottom water NO"[{"x"}]* " concentration (μmol L"^{"-1"}*")"),
                                 expression(bold("NOX_Inv:") ~ "porewater NO"[{"x"}]* " inventory (μmol cm"^{"-2"}*")"),
                                 expression(bold("O2_BW:") ~ "bottom water O"[2]* " concentration (µmol L"^{"-1"}*")"),
                                 expression(bold("O2_Int:") ~ "sediment depth-integrated O"[2]* " consumption rate (μmol"^{"-1"}* " cm"^{"-2"}* " s"^{"-1"}*")"),
                                 expression(bold("O2_Pen:") ~ "sediment O"[2]* " penetration depth (mm)"),
                                 expression(bold("P_BW:") ~ "bottom water P concentration (μmol L"^{"-1"}*")"),
                                 expression(bold("P_Inv:") ~ "porewater P inventory (μmol cm"^{"-2"}*")"),
                                 expression(bold("Sal_BW:") ~ "water column salinity (no unit)"),
                                 expression(bold("Temp_BW:") ~ "water column temperature (°C)"),
                                 expression(bold("X13C_1cm:") ~ "sediment δ"^13* "C"[{"org"}]* ", top 1 cm (Δ PDB)")
                      )) + 
  scale_colour_discrete_qualitative(breaks = c("Included in Model 1", 
                                               "Included in Model 2"),
                                    palette = "Harmonic") +
  guides(fill = guide_legend(override.aes = list(fill = "black"))) +
  theme(legend.text=element_text(size = 12)) +
  theme(legend.text.align = 0)

pca.plot

dev.off()
