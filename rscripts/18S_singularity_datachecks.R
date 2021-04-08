# Seili 18S analysis - jesse harrison 2020-2021
# using seili-r Singularity container (based on seili-r.def)
# additional libpath ####
# (see extra_RPackages.R for extra package installs)

.libPaths(c("/home/jharriso/seili-singularity/rpackages", .libPaths()))

# packages ####

packages <- c("phyloseq", "ggplot2", "vegan", "grid", "gridExtra", "data.table", "plyr", "cowplot", 
              "RVAideMemoire", "microbiome", "knitr", "RColorBrewer", "Cairo", "multcompView", 
              "QsRutils", "dplyr", "ggvegan")

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

# load 18S RData ####

load("rdata/Seili18s.RData")
# based on: import_biom("tables_sampledata_export.biom", treefilename="tree_rooted.tree", refseqfilename="otus.fasta")
# 'rawdata' is a version of the dataset where no step has been taken to remove NAs


# dataset with NAs removed at phylum level ####

rawdata.noNA <- subset_taxa(rawdata, 
                                !is.na(Rank2) & 
                                  !Rank2 %in% c("", "uncharacterized"))

table(tax_table(rawdata.noNA)[, "Rank2"], exclude = NULL)
# the above shows we also have one phylum called "Ambiguous_taxa", treat it as NA

filterPhyla.raw <- c("Ambiguous_taxa")
rawdata.noNA <- subset_taxa(rawdata.noNA, !Rank2 %in% filterPhyla.raw)

# further pre-processing (prevalences) ####

# source script with fast_melt() function
# ("melts" OTU table into format with three main columns: SampleID, TaxaID and count)
# taxa_summary.R can be downloaded from: http://evomics.org/phyloseq/taxa_summary-r/

source("rscripts/taxa_summary.R", local = TRUE)

# pre-processing following steps in analysis script

# fast_melt
mdt.raw <- fast_melt(rawdata)
mdt.raw2 <- fast_melt(rawdata.noNA)

# calculating relative abundances
mdt.raw[, RelativeAbundance := count / sum(count), by = SampleID]
prevdt.raw = mdt.raw[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

mdt.raw2[, RelativeAbundance := count / sum(count), by = SampleID]
prevdt.raw2 = mdt.raw2[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

# prevalence calculations
prevdf.raw = apply(X = otu_table(rawdata),
                   MARGIN = ifelse(taxa_are_rows(rawdata), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

prevdf.raw = data.frame(Prevalence = prevdf.raw,
                        TotalAbundance = taxa_sums(rawdata),
                        tax_table(rawdata))

prevdf.raw2 = apply(X = otu_table(rawdata.noNA),
                    MARGIN = ifelse(taxa_are_rows(rawdata.noNA), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

prevdf.raw2 = data.frame(Prevalence = prevdf.raw2,
                         TotalAbundance = taxa_sums(rawdata.noNA),
                         tax_table(rawdata.noNA))

# phylum-level prevalences
prevdf.raw.phyla = subset(prevdf.raw, Rank2 %in% get_taxa_unique(rawdata, "Rank2"))
prevdf.raw.phyla2 = subset(prevdf.raw2, Rank2 %in% get_taxa_unique(rawdata.noNA, "Rank2"))

# 5% prevalence filtering ####

# defining 5% prevalence thresholds for filtering
prevalenceThreshold.raw <- 0.05 * nsamples(rawdata)
prevalenceThreshold.raw2 <- 0.05 * nsamples(rawdata.noNA)

# execute prevalence filter, using `prune_taxa()` function
keepTaxa.raw <- rownames(prevdf.raw.phyla)[(prevdf.raw.phyla$Prevalence >= prevalenceThreshold.raw)]
rawdata.2 <- prune_taxa(keepTaxa.raw, rawdata)

keepTaxa.raw2 = rownames(prevdf.raw.phyla2)[(prevdf.raw.phyla2$Prevalence >= prevalenceThreshold.raw2)]
rawdata.noNA.2 = prune_taxa(keepTaxa.raw2, rawdata.noNA)

# Hellinger transformation ####

# taxon composition plots showed as proportions (%) for readability
# note: here not using CLR-transformed data but Hellinger (for visualisation)

raw2.ra <- microbiome::transform(rawdata.2, 'hellinger')
raw2.ra <- transform_sample_counts(raw2.ra, 
                                   function(x) x/sum(x))

raw2.noNA.ra <- microbiome::transform(rawdata.noNA.2, 'hellinger')
raw2.noNA.ra <- transform_sample_counts(raw2.noNA.ra, 
                                        function(x) x/sum(x))
# stacked bar plot using Hellinger-transformed data ####

colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
             "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
             "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
             "#771122", "#AA4455", "#DD7788", "#DD776655")

# subset the data to domain level (= rank1)
raw2.ra.r1 <- tax_glom(raw2.ra, taxrank = rank_names(raw2.ra)[1])
raw2.ra.r1 <- psmelt(raw2.ra.r1)

Cairo(file = "figures/r_output/FigS2b_18S.png", 
      type = "png", 
      units = "cm", 
      width = 27, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

r1.plot <- ggplot(raw2.ra.r1 , 
                  aes(x = Replicate, y = Abundance, fill = Rank1)) + 
  geom_bar(stat = "identity", 
           position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))

r1.plot
dev.off()
# CLR transformation ####

# CLR transformation (microbiome package) for data using 5% prevalence filter

rawdata.2.clr <- transform(rawdata.2, 'clr')
rawdata.noNA.2.clr <- transform(rawdata.noNA.2, 'clr')


# nMDS using CLR-transformed data ####

set.seed(1)
# using Aitchison distance (CLR + Euclidean distances)
rawdata.2.nmds <- ordinate(physeq = rawdata.2.clr, 
                           method = "NMDS", 
                           distance = "euclidean") # stress = 0.16

set.seed(1)
rawdata.noNA.2.nmds <- ordinate(physeq = rawdata.noNA.2.clr, 
                                method = "NMDS", 
                                distance = "euclidean") # stress = 0.16

# nMDS plots ####

Cairo(file = "figures/r_output/FigS3d_18S.png", 
      type = "png", 
      units = "cm", 
      width = 15, 
      height = 15, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

nMDS.plot.rawdata.2.clr <- plot_ordination(rawdata.2.clr, 
                                           rawdata.2.nmds, 
                                           color = "Site") + 
  theme(aspect.ratio = 1) + 
  geom_point(size = 6.5, alpha = 0.75) + 
  scale_colour_brewer(type = "qual", 
                      palette = "Dark2") + 
  theme(axis.title = element_blank()) + 
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

nMDS.plot.rawdata.2.clr
dev.off()

Cairo(file = "figures/r_output/FigS3e_18S.png", 
      type = "png", 
      units = "cm", 
      width = 15, 
      height = 15, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

nMDS.plot.rawdata.noNA.2.clr <- plot_ordination(rawdata.noNA.2.clr, 
                                                rawdata.noNA.2.nmds, 
                                                color = "Site") + 
  theme(aspect.ratio = 1) + 
  geom_point(size = 6.5, alpha = 0.75) + 
  scale_colour_brewer(type = "qual", 
                      palette = "Dark2") + 
  theme(axis.title = element_blank()) + 
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

nMDS.plot.rawdata.noNA.2.clr
dev.off()
