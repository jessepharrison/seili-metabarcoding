# Seili 16S analysis - jesse harrison 2020-2021
# using seili-r Singularity container (based on seili-r.def)

# additional libpath ####
# (see extra_RPackages.R for extra package installs)
.libPaths(c("/home/jharriso/seili-singularity/rpackages", .libPaths()))

# packages ####

packages <- c("phyloseq", "ggplot2", "vegan", "grid", "gridExtra", "data.table", "plyr", "cowplot", 
              "RVAideMemoire", "microbiome", "knitr", "RColorBrewer", "Cairo", "multcompView")

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

setwd("/home/jharriso/GitHub/seili-metabarcoding/")

# load 16S RData ####

load("rdata/Seili16s.RData")
# based on: import_biom("tables_sampledata_export.biom", treefilename="tree_rooted.tree", refseqfilename="otus.fasta")

# rarefaction curve (without rarefying to even depth) ####

Cairo(file = "figures/r_output/FigS1a_16S.png", 
      type = "png", 
      units = "cm", 
      width = 25, 
      height = 25, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

rarecurve(t(otu_table(rawdata)), 
          step = 100, 
          cex.lab = 1.5, 
          cex.axis = 1.5, 
          label = FALSE,
          xlab = "No. of sequences",
          ylab = "OTUs")

dev.off()

# alpha diversity and evenness (without rarefying to even depth) ####

richness <- estimate_richness(rawdata, measures = c("Observed", "Chao1", "Shannon"))
richness <- cbind(richness, rawdata@sam_data[,1:2]) # add replicate and site
write.csv(richness, "tables/richness_rawdata_16S.csv")

# alpha diversity and evenness (rarefied to even depth - for completeness, but not used) ####

# set.seed(1)
# alpharare <- rarefy_even_depth(rawdata)

# estimate_richness(alpharare, measures=c("Observed", "Chao1", "Shannon"))
# evenness(alpharare, 'pielou')

# data inspection + filtering (non-rarefied) ####

# no. of samples and taxa
nsamples(rawdata)
ntaxa(rawdata)

# per-sample no.s of sequences
sample_sums(rawdata)

# sample variables
sample_variables(rawdata)

# taxonomy table (phyla)
table(tax_table(rawdata)[, "Rank2"], exclude = NULL)

# a few phyla that occur only once: candidate division WPS-2, Chlamydiae, Thermotogae
# there are also 6035 NAs that are likely to be artifacts

# first, take out NAs, mitochondrial + chloroplast sequences + make new tax table
rawdata.bac <- subset_taxa(rawdata, 
                           Rank1 == "Bacteria" & 
                             Rank3 != "Chloroplast")

table(tax_table(rawdata.bac)[, "Rank2"], exclude = NULL)

# no. of taxa after removing mitochondrial + chloroplast sequences
ntaxa(rawdata.bac)

# workflow based on: https://f1000research.com/articles/5-1492/v2
# this step does nothing as the previous filter step already removed NAs, but keeping it for sake of
# consistency between scripts.

rawdata.bac.0 <- subset_taxa(rawdata.bac, 
                             !is.na(Rank2) & 
                               !Rank2 %in% c("", "uncharacterized"))

table(tax_table(rawdata.bac.0)[, "Rank2"], exclude = NULL)

# next, inspect the data in more detail to see if certain phyla are comprised mostly of low-prevalence
# features (prevalence = no. of times an OTU is observed at least once)

# compute prevalence of each feature, store as data.frame
prevdf.raw <- apply(X = otu_table(rawdata.bac.0),
                   MARGIN = ifelse(taxa_are_rows(rawdata.bac.0), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

# add taxonomy and total read counts to data.frame
prevdf.raw <- data.frame(Prevalence = prevdf.raw,
                        TotalAbundance = taxa_sums(rawdata.bac.0),
                        tax_table(rawdata.bac.0))

plyr::ddply(prevdf.raw, "Rank2", 
            function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# the phyla that occurred only once are also of very low prevalence (Chlamydiae: 1, Thermotogae: 4, WPS-2: 23).
# out of these, would filter Chlamydiae and Thermotogae.

# define phyla to filter and remove them from the data set
filterPhyla.raw <- c("Chlamydiae", "Thermotogae")
rawdata.1 <- subset_taxa(rawdata.bac.0, !Rank2 %in% filterPhyla.raw)

table(tax_table(rawdata.1)[, "Rank2"], exclude = NULL)

# additional ways to inspect OTU counts ####

# source script with fast_melt() function
# ("melts" OTU table into format with three main columns: SampleID, TaxaID and count)
# taxa_summary.R can be downloaded from: http://evomics.org/phyloseq/taxa_summary-r/

source("rscripts/taxa_summary.R", local = TRUE)

# fast_melt rawdata1
mdt.raw <- fast_melt(rawdata.1)

# omit NAs
mdt.raw <- mdt.raw[!is.na(count)]

# calculate relative abundance
mdt.raw[, RelativeAbundance := count / sum(count), by = SampleID]

prevdt.raw <- mdt.raw[, list(Prevalence = sum(count > 0), 
                             TotalCounts = sum(count)), by = TaxaID]

# draw prevalence plot
ggplot(prevdt.raw, aes(Prevalence)) + 
  geom_histogram()

# how many with 0 seqs?
prevdt.raw[(Prevalence <= 0), .N] ## [1] 0

# how many singletons?
prevdt.raw[(Prevalence <= 1), .N] ## [1] 298

# how many doubletons?
prevdt.raw[(Prevalence <= 2), .N] ## [1] 725

# taxa cumulative sum with prevalence on x axis
prevcumsum.raw <- prevdt.raw[, .N, by = Prevalence]
setkey(prevcumsum.raw, Prevalence)
prevcumsum.raw[, CumSum := cumsum(N)]
pPrevCumSum.raw = ggplot(prevcumsum.raw, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering threshold (prevalence)") +
  ylab("OTUs filtered")
pPrevCumSum.raw

# prevalence vs. total count scatter plot
ggplot(prevdt.raw, 
       aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

# prevalence plot for phyla

Cairo(file = "figures/r_output/FigS4a_16S.png", 
      type = "png", 
      units = "cm", 
      width = 45, 
      height = 33, 
      pointsize = 14, 
      dpi = 300, 
      bg = "white")

prevdf.raw.phyla <- subset(prevdf.raw, 
                           Rank2 %in% get_taxa_unique(rawdata.1, "Rank2"))

ggplot(prevdf.raw.phyla, 
       aes(TotalAbundance, 
           Prevalence / nsamples(rawdata.bac.0),
           color = Rank2)) +
  # Include a guess for parameter (5 % of total samples)
  geom_hline(yintercept = 0.05, 
             alpha = 0.5, 
             linetype = 2) + 
  geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +
  facet_wrap(~Rank2) + 
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20)) +
  theme(strip.text.x = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  xlab("\nTotal abundance") +
  ylab("Prevalence [frac. samples]\n")

dev.off()

# 5% prevalence filter ####

# define 5% prevalence threshold for filtering
prevalenceThreshold.raw <- 0.05 * nsamples(rawdata.bac.0)

# execute prevalence filter, using `prune_taxa()` function
keepTaxa.raw <- rownames(prevdf.raw.phyla)[(prevdf.raw.phyla$Prevalence >= prevalenceThreshold.raw)]
rawdata.2 <- prune_taxa(keepTaxa.raw, rawdata.bac.0)

# data transformation ####

# CLR transformation (microbiome package) for data using 5% prevalence filter (rawdata.2)

rawdata.2.clr <- transform(rawdata.2, 'clr')

# Hellinger transformation
# (Included for completeness in case required, commented out)

# rawdata.2.hellinger <- rawdata.2
# otu_table(rawdata.2.hellinger) <- otu_table(decostand(otu_table(rawdata.2.hellinger), 
#                                                      method = "hellinger"), taxa_are_rows = TRUE)

# % transformation (commented out)
# rawdata.2.hellinger <- transform_sample_counts(rawdata.2, function(x) x/sum(x))

# stacked taxon composition plots for %RA data (no CLR) ####

# taxon composition plots showed as proportions (%) for readability
# note: here not using CLR-transformed data

# relative abundance conversion
raw2.ra <- transform_sample_counts(rawdata.2, 
                                   function(x) x/sum(x))

# subset the data to phylum level (= rank2), cut out low-abundance taxa
raw2.ra.phylum <- tax_glom(raw2.ra, 
                           taxrank = rank_names(raw2.ra)[2], 
                           NArm = T,
                           bad_empty = c(NA,""," ","\t"))

# melt to long format (also converts to data frame)
raw2.ra.phylum <- psmelt(raw2.ra.phylum)

# cut out OTUs less than 2%
raw2.ra.phylum.02 <- subset(raw2.ra.phylum, 
                            Abundance > 0.02)

# % retained after removing OTUs < 2%
sum(raw2.ra.phylum.02$Abundance)/sum(raw2.ra.phylum$Abundance) # 91.4%

# define colour palette
colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
             "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
             "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", 
             "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", 
             "#DD7788", "#CBD588")

# phylum-level bar plot
Cairo(file = "figures/r_output/Fig4a_16S.png", 
      type = "png", 
      units = "cm", 
      width = 27, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

phylum.plot <- ggplot(raw2.ra.phylum.02 , 
                      aes(x = Replicate, 
                          y = Abundance, fill = Rank2)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >2%)\n") +
  theme(axis.text.y = element_text(size = 12, 
                                   hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, 
                                   angle = 30, 
                                   vjust = 0.5, 
                                   hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))

phylum.plot
dev.off()

# class level plot (= rank3)
raw2.ra.class <- tax_glom(raw2.ra, 
                          taxrank = rank_names(rawdata.2)[3], 
                          NArm = T,
                          bad_empty = c(NA,""," ","\t"))

# psmelt raw2.ra.class
raw2.ra.class <- psmelt(raw2.ra.class)

# cut out OTUs <3%
raw2.ra.class.03 <- subset(raw2.ra.class, 
                           Abundance > 0.03)

# % remaining after removing OTUs <10%
sum(raw2.ra.class.03$Abundance)/sum(raw2.ra.class$Abundance) # retains 75% of original data

Cairo(file = "figures/r_output/FigS6_16S.png", 
      type = "png", units = "cm", 
      width = 27, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

class.plot <- ggplot(raw2.ra.class.03, 
                        aes(x = Replicate, 
                            y = Abundance, 
                            fill = Rank3)) + 
  geom_bar(stat = "identity", 
           position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >3%)\n") +
  theme(axis.text.y = element_text(size = 12, 
                                   hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, 
                                   angle = 30, 
                                   vjust = 0.5, 
                                   hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))

class.plot
dev.off()

# nMDS for CLR-transformed data ####

set.seed(1)
# using Aitchison distance (CLR + Euclidean distances)
rawdata.2.nmds <- ordinate(physeq = rawdata.2.clr, 
                           method = "NMDS", 
                           distance = "euclidean") # stress = 0.04

Cairo(file = "figures/r_output/Fig3a_16S.png", 
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
  theme(legend.position = "bottom")

nMDS.plot.rawdata.2.clr
dev.off()

# global + pairwise PERMANOVA (~ site) and PERMDISP ####

# make a data frame from the sample_data
rawdata.2.clr.df <- data.frame(sample_data(rawdata.2.clr))

# extract Aitchison distances
set.seed(1)
rawdata.2.clr.euc <- phyloseq::distance(rawdata.2.clr, 
                                        method = 'euclidean')

# global PERMANOVA (after setting seed = 1)
set.seed(1)
global.permanova <- adonis(rawdata.2.clr.euc ~ Site,
                           data = rawdata.2.clr.df)

# save global PERMANOVA output

cat(" ------------------------------------",
    "\n", "16S - GLOBAL PERMANOVA", "\n",
    "------------------------------------",
    "\n",
    file = "stats/permanova_permdisp_16S.txt")

capture.output(global.permanova, 
               file = "stats/permanova_permdisp_16S.txt",
               append = TRUE)

# post-hoc pairwise PERMANOVA (after setting seed = 1)
set.seed(1)
pairwise.permanova <- pairwise.perm.manova(rawdata.2.clr.euc,
                                           rawdata.2.clr.df$Site,
                                           nperm = 999,
                                           p.method = "BH")

# save pairwise PERMANOVA output

cat("\n",
    "------------------------------------",
    "\n", "16S - PAIRWISE PERMANOVA", "\n",
    "------------------------------------",
    "\n",
    file = "stats/permanova_permdisp_16S.txt",
    append = TRUE)

capture.output(pairwise.permanova, 
               file = "stats/permanova_permdisp_16S.txt", 
               append = TRUE)

# betadisper / PERMDISP (after setting seed = 1)
set.seed(1)
beta <- betadisper(rawdata.2.clr.euc, 
                   rawdata.2.clr.df$Site)
set.seed(1)
permdisp <- permutest(beta)

# save PERMDISP output

cat("\n",
    "------------------------------------",
    "\n", "16S - PERMDISP", "\n",
    "------------------------------------",
    "\n",
    file = "stats/permanova_permdisp_16S.txt",
    append = TRUE)

capture.output(permdisp, 
               file = "stats/permanova_permdisp_16S.txt", 
               append = TRUE)
