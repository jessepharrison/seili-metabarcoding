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

# disable scientific notation
options(scipen=10000)

# working directory ####

setwd("/home/jharriso/git/seili-metabarcoding/")

# load 18S RData ####

load("rdata/Seili18s.RData")
# based on: import_biom("tables_sampledata_export.biom", treefilename="tree_rooted.tree", refseqfilename="otus.fasta")

# rarefaction curve (without rarefying to even depth) ####

Cairo(file = "figures/r_output/FigS1b_18S.tiff", 
      type = "tiff", 
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
          ylab = "")

dev.off()

# alpha diversity and Good's coverage (without rarefying to even depth) ####

richness <- estimate_richness(rawdata, measures = c("Observed", "Chao1", "Shannon"))
richness <- cbind(richness, rawdata@sam_data[,1:2]) # add replicate and site
write.csv(richness, "tables/richness_rawdata_18S.csv")

# Good's coverage

psotu2veg <- function(physeq) {
   OTU <- otu_table(physeq)
   if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
   }
   return(as(OTU, "matrix"))
}

rawdata.df <- psotu2veg(rawdata)
goods(rawdata.df)

# alpha diversity (rarefied to even depth - for completeness, but not used) ####

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

# one phylum that occurs only once: "Ambiguous_taxa"
# there are also 225 NAs that are likely to be artifacts

# first, take out ambiguous sequences and make new tax table
rawdata.0 <- subset_taxa(rawdata, !is.na(Rank2) & !Rank2 %in% c("", "uncharacterized"))
table(tax_table(rawdata.0)[, "Rank2"], exclude = NULL) # exclude NAs

# no. of taxa after removing ambiguous sequences
ntaxa(rawdata.0)

# next, inspect the data in more detail to see if certain phyla are comprised mostly of low-prevalence
# features (prevalence = no. of times an OTU is observed at least once)

# compute prevalence of each feature, store as data.frame
prevdf.raw = apply(X = otu_table(rawdata.0),
                   MARGIN = ifelse(taxa_are_rows(rawdata.0), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

# add taxonomy and total read counts to data.frame
prevdf.raw = data.frame(Prevalence = prevdf.raw,
                        TotalAbundance = taxa_sums(rawdata.0),
                        tax_table(rawdata.0))

plyr::ddply(prevdf.raw, "Rank2", 
            function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# "Ambiguous_taxa" has low prevalence (2), so it should be safe to filter.

# define phyla to filter and remove them from the data set
filterPhyla.raw <- c("Ambiguous_taxa")
rawdata.1 <- subset_taxa(rawdata.0, !Rank2 %in% filterPhyla.raw)

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
prevdt.raw[(Prevalence <= 0), .N]

# how many singletons?
prevdt.raw[(Prevalence <= 1), .N]

# how many doubletons?
prevdt.raw[(Prevalence <= 2), .N]

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

Cairo(file = "figures/r_output/FigS4b_18S.tiff", 
      type = "tiff", 
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
           Prevalence / nsamples(rawdata.0),
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
prevalenceThreshold.raw <- 0.05 * nsamples(rawdata.1)

# execute prevalence filter, using `prune_taxa()` function
keepTaxa.raw <- rownames(prevdf.raw.phyla)[(prevdf.raw.phyla$Prevalence >= prevalenceThreshold.raw)]
rawdata.2 <- prune_taxa(keepTaxa.raw, rawdata.1)

# summarise sequence and OTU no.s ####

# no. of sequences
sample_sums(rawdata.2)

# overall sum of sequences
sum(sample_sums(rawdata.2))

# no. of taxa
ntaxa(rawdata.2)

# CLR transformation ####

# CLR transformation (microbiome package) for data using 5% prevalence filter (rawdata.2)

rawdata.2.clr <- transform(rawdata.2, 'clr')

# stacked taxon composition plot for %RA data (Hellinger transformation) ####

# taxon composition plots showed as proportions (%) for readability
# note: here not using CLR-transformed data but Hellinger (for visualisation)

# Hellinger transformation + relative abundance conversion
raw2.ra <- microbiome::transform(rawdata.2, 'hellinger')
raw2.ra <- transform_sample_counts(raw2.ra, 
                                   function(x) x/sum(x))

# subset the data to Rank 4 of the SILVA majority taxonomy mapping file, cut out low-abundance taxa
# (note that ranks comprise multiple levels of taxonomic classification)

raw2.ra.r4 <- tax_glom(raw2.ra, 
                       taxrank = rank_names(raw2.ra)[4], 
                       NArm = T,
                       bad_empty = c(NA,""," ","\t"))

# melt to long format (also converts to data frame)
raw2.ra.r4 <- psmelt(raw2.ra.r4)

# cut out OTUs less than 2%
raw2.ra.r4.02 <- subset(raw2.ra.r4, 
                        Abundance > 0.02)

# % retained after removing OTUs < 2%
sum(raw2.ra.r4.02$Abundance)/sum(raw2.ra.r4$Abundance) # 93.9%

# define colour palette
colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
             "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
             "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", 
             "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", 
             "#DD7788", "#CBD588")

# rank 4 bar plot
Cairo(file = "figures/r_output/Fig4b_18S.tiff", 
      type = "tiff", 
      units = "cm", 
      width = 27, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

r4.plot <- ggplot(raw2.ra.r4.02 , 
                      aes(x = Replicate, 
                          y = Abundance, fill = Rank4)) + 
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

r4.plot
dev.off()

# Metazoa could be of further interest.

# create Metazoa subset ####

# for visualising subset of the data including Metazoa only,
# use earlier Hellinger-transformed %RA data set

metazoa.ra <- subset_taxa(raw2.ra, Rank4 == "Metazoa (Animalia)")

# stacked taxon composition plot for %RA data, focusing on Metazoa (Hellinger transformation) ####

# subset the data to Rank 6 of the SILVA majority taxonomy mapping file, cut out low-abundance taxa
# (note that ranks comprise multiple levels of taxonomic classification)

metazoa.ra.r6 <- tax_glom(metazoa.ra, 
                       taxrank = rank_names(metazoa.ra)[6], 
                       NArm = T,
                       bad_empty = c(NA,""," ","\t"))

# melt to long format (also converts to data frame)
metazoa.ra.r6 <- psmelt(metazoa.ra.r6)

# cut out OTUs less than 0.5%
metazoa.ra.r6.005 <- subset(metazoa.ra.r6, 
                            Abundance > 0.005)

# % retained after removing OTUs <2%
sum(metazoa.ra.r6.005$Abundance)/sum(metazoa.ra.r6$Abundance) # retains 92.3% of original data

# rank 6 bar plot
Cairo(file = "figures/r_output/Fig5b_18S.tiff", 
      type = "tiff", 
      units = "cm", 
      width = 27, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

meta.plot <- ggplot(metazoa.ra.r6.005, 
                  aes(x = Replicate, 
                      y = Abundance, fill = Rank6)) + 
   geom_bar(stat = "identity", position = "stack") +
   facet_wrap(~ Site) +
   scale_fill_manual(values = colours) +
   scale_y_continuous(name = "Relative abundance (taxa >0.5%)\n") +
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

meta.plot
dev.off()

# Mann-Whitney U test for comparing Shannon diversity (S1-3 vs S4-6) of Metazoa ####

# the taxon composition bar plot for Metazoa suggests that sites S1-S3 may be 
# less diverse than S4-S6.

# first let's calculate diversity indices for Metazoa.
# it's important to note that, for this, we should use the original untransformed data set.

# subset the original raw data
metazoa.raw <- subset_taxa(rawdata, Rank4 == "Metazoa (Animalia)")

# calculate diversity indices
richness.meta <- estimate_richness(metazoa.raw, measures = c("Observed", "Chao1", "Shannon"))
richness.meta <- cbind(richness.meta, metazoa.raw@sam_data[,1:2]) # add replicate and site
write.csv(richness.meta, "tables/richness_rawdata_metazoaonly_18S.csv")

# add a grouping column to the data frame
groupvec <- rep(c("S1-S3", "S4-S6"), each = nrow(richness.meta) / 2)
richness.meta <- cbind(richness.meta, Grouping = groupvec)

# quick distribution check
ggplot(richness.meta, aes(x = Shannon)) +
   geom_histogram(color = "black", fill = "white") +
   facet_wrap(facets = vars(Grouping)) +
   theme_bw() +
   theme(text = element_text(size = 16))

# likely best to use a non-parametric test to be on the safe side.

# perform Mann-Whitney U test
richness.meta.mann <- wilcox.test(Shannon ~ Grouping, richness.meta)

# save Mann-Whitney U test results

cat(" -------------------------------------------------------------------",
    "\n", "18S - Shannon index (S1-S3 vs S4-S6) - Mann-Whitney U", "\n",
    "-------------------------------------------------------------------",
    "\n",
    file = "stats/mannwhitney_shannon_metazoa_18S.txt")

capture.output(richness.meta.mann, 
               file = "stats/mannwhitney_shannon_metazoa_18S.txt",
               append = TRUE)

# nMDS for CLR-transformed data ####

set.seed(1)
# using Aitchison distance (CLR + Euclidean distances)
rawdata.2.nmds <- ordinate(physeq = rawdata.2.clr, 
                           method = "NMDS", 
                           distance = "euclidean") # stress = 0.16

Cairo(file = "figures/r_output/Fig3b_18S.tiff", 
      type = "tiff", 
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
    "\n", "18S - GLOBAL PERMANOVA", "\n",
    "------------------------------------",
    "\n",
    file = "stats/permanova_permdisp_18S.txt")

capture.output(global.permanova, 
               file = "stats/permanova_permdisp_18S.txt",
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
    "\n", "18S - PAIRWISE PERMANOVA", "\n",
    "------------------------------------",
    "\n",
    file = "stats/permanova_permdisp_18S.txt",
    append = TRUE)

capture.output(pairwise.permanova, 
               file = "stats/permanova_permdisp_18S.txt", 
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
    "\n", "18S - PERMDISP", "\n",
    "------------------------------------",
    "\n",
    file = "stats/permanova_permdisp_18S.txt",
    append = TRUE)

capture.output(permdisp, 
               file = "stats/permanova_permdisp_18S.txt", 
               append = TRUE)

# contrary to 16S data set, for 18S we see a significant PERMDISP result.

# pairwise Tukey HSD comparisons for sites (following significant PERMDISP) ####

betaHSD <- TukeyHSD(beta)

# save Tukey HSD output

cat(" ------------------------------------",
    "\n", "18S - TUKEY HSD", "\n",
    "------------------------------------",
    "\n",
    file = "stats/tukey_18S.txt")

capture.output(betaHSD, 
               file = "stats/tukey_18S.txt",
               append = TRUE)

# pairwise Tukey HSD plot ####

# compact letter display
# (not saved but used to inform plot construction)
multcompLetters(extract_p(betaHSD$group))

# manual box plot

Cairo(file = "figures/r_output/Fig6_18S.tiff", 
      type = "tiff", 
      units = "cm", 
      width = 15, 
      height = 15, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

boxplot(beta, xlab = "Site", ylim = c(18, 37))

text(x = 1, y = 27.9, labels= "c")
text(x = 2, y = 33.0, labels= "ab")
text(x = 3, y = 29.5, labels= "ac")
text(x = 4, y = 36.5, labels= "b")
text(x = 5, y = 25.3, labels= "c")
text(x = 6, y = 30.9, labels= "abc")

dev.off()


# db-RDA pt 1: load porosity-corrected environmental data and scale it (mean = 0, sd = 1) ####

# load the environmental data
# note: same values as in SeiliEnv_General_PorCorr.csv 
# (but with phyloseq-friendly formatting)

envdata <- read.csv("envdata/envfit/16S_envdata_envfit.csv")

# scale to mean zero and unit variance

envdata <- as.data.frame(scale(envdata))

# trim envdata from five replicates per site to four (to match 18S dataset)

nth.delete <- function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
envdata <- nth.delete(envdata, 5) # deleting every 5th row

# db-RDA pt 2: merge scaled env data with phyloseq object ####

# add sequencing sample IDs to env data row names
# (required for correct merging with phyloseq object)
rownames(envdata) <- rownames(sample_data(rawdata.2.clr))

# merge env data with phyloseq object
envdata <- sample_data(envdata)
rawdata.2.clr <- merge_phyloseq(rawdata.2.clr, envdata)









# db-RDA pt 3: capscale for 18S models ####

# db-RDA test 3
# (keeping annotation the same here as what is used in the 16S script)
# Farm + CN_1cm + IntO2 + NH4_Inv + grain

# convert phyloseq OTU table into vegan format
rawdata.2.clr.otus <- veganotu(rawdata.2.clr)

# extract sample data from phyloseq object
rawdata.2.clr.samdata <- data.frame(sample_data(rawdata.2.clr))

# select continuous variables (i.e. drop Replicate and Site)
rawdata.2.clr.samdata.numeric <- select_if(rawdata.2.clr.samdata, is.numeric)

set.seed(1)
viftest.res3 <- capscale(rawdata.2.clr.otus ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain, 
                         data = rawdata.2.clr.samdata.numeric, 
                         distance = 'euclidean')

# calculate VIFs
viftest.res3.vifs <- sort(vif.cca(viftest.res3))

# drop1() for db-RDA test 3
set.seed(1)
viftest.res3.drop1 <- drop1(viftest.res3, test = "perm")

# save db-RDA test 3 output

cat(" ------------------------------------",
    "\n", "18S - test model 3", "\n",
    "------------------------------------",
    "\n", "\n",
    file = "stats/vifs_testmodels_18S.txt")

capture.output(viftest.res3, 
               file = "stats/vifs_testmodels_18S.txt",
               append = TRUE)

cat(" Model 3 VIFs (lowest to highest)", "\n", "\n",
    file = "stats/vifs_testmodels_18S.txt",
    append = TRUE)

capture.output(viftest.res3.vifs, 
               file = "stats/vifs_testmodels_18S.txt",
               append = TRUE)

cat("\n", "Model 3 drop1() analysis", "\n", "\n",
    file = "stats/vifs_testmodels_18S.txt",
    append = TRUE)

capture.output(viftest.res3.drop1, 
               file = "stats/vifs_testmodels_18S.txt",
               append = TRUE)

# db-RDA test 5
# C_1cm + O2_BW + NOX_BW + P_BW

set.seed(1)
viftest.res5 <- capscale(rawdata.2.clr.otus ~ C_1cm + O2_BW + NOX_BW + P_BW, 
                         data = rawdata.2.clr.samdata.numeric, 
                         distance = 'euclidean')

# calculate VIFs
viftest.res5.vifs <- sort(vif.cca(viftest.res5))

# drop1() for db-RDA test 5
set.seed(1)
viftest.res5.drop1 <- drop1(viftest.res5, test = "perm")

# save db-RDA test 5 output

cat("\n",
    "------------------------------------",
    "\n", "18S - test model 5", "\n",
    "------------------------------------",
    "\n", "\n",
    file = "stats/vifs_testmodels_18S.txt",
    append = TRUE)

capture.output(viftest.res5, 
               file = "stats/vifs_testmodels_18S.txt",
               append = TRUE)

cat(" Model 5 VIFs (lowest to highest)", "\n", "\n",
    file = "stats/vifs_testmodels_18S.txt",
    append = TRUE)

capture.output(viftest.res5.vifs, 
               file = "stats/vifs_testmodels_18S.txt",
               append = TRUE)

cat("\n", "Model 5 drop1() analysis", "\n", "\n",
    file = "stats/vifs_testmodels_18S.txt",
    append = TRUE)

capture.output(viftest.res5.drop1, 
               file = "stats/vifs_testmodels_18S.txt",
               append = TRUE)

# db-RDA pt 4: perform db-RDA + envfit for VIF-based model ####

# See 16S script for VIF calculations

# extract sample data from phyloseq object
rawdata.2.clr.samdata <- data.frame(sample_data(rawdata.2.clr))

# db-RDA
set.seed(1)
vifmod.dbrda <- ordinate(rawdata.2.clr, 
                         method = "CAP", 
                         distance = "euclidean",
                         formula = ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain)

# envfit (as a data frame using fortify)
# calculated on linear combinations of scores
# see: https://www.davidzeleny.net/anadat-r/doku.php/en:confusions

set.seed(1)
envfit.lc <- fortify(
   envfit(vifmod.dbrda ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain,
          data = rawdata.2.clr.samdata,
          permutations = 999, display = 'lc')
)





# db-RDA pt 5: plot db-RDA + envfit for VIF-based model ####

# Note: arbitrary multiplier addded to xend + yend
# for plotting purposes (could also be handled through scaling)

# Layout for Fig. 6:
# a) 16S, VIF-based model (16S script)
# b) 18S, VIF-based model (this script)
# c) 16S, Monitoring-based model (16S script)
# d) 18S, Monitoring-based model (this script)

Cairo(file = "figures/r_output/Fig7b_18S.tiff", 
      type = "tiff", 
      units = "cm", 
      width = 20, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

vifmod.dbrda.plot <- plot_ordination(rawdata.2.clr, 
                                     vifmod.dbrda, 
                                     color = "Site") +
   theme(aspect.ratio = 1) +
   scale_colour_brewer(type = "qual", 
                       palette = "Dark2") + 
   geom_point(size = 6.5, alpha = 0.75) +
   geom_segment(data = envfit.lc,
                aes(x = 0, 
                    y = 0, 
                    xend = CAP1*6, 
                    yend = CAP2*6), 
                size = 1, 
                alpha = 0.5, 
                colour = "grey30") +
   geom_label(data = envfit.lc, 
              aes(x = CAP1*6, 
                  y = CAP2*6),
              label = envfit.lc$Label, 
              colour = "black", 
              fontface = "bold",
              size = 3.5) + 
   theme(axis.title = element_blank()) +
   theme(axis.ticks = element_blank()) +
   theme(axis.text = element_blank()) +
   theme(legend.position = "none")

vifmod.dbrda.plot
dev.off()





# db-RDA pt 6: permutation test for db-RDA explanatory variables (VIF-based model) ####

# global test (999 permutations)
dbrda.glob <- anova(vifmod.dbrda, permutations = 999)

# save test results

cat(" -------------------------------------------------------------------",
    "\n", "18S - model 1 - global test", "\n",
    "-------------------------------------------------------------------",
    "\n", "\n",
    file = "tables/permtest_global_mod1_18S.txt")

capture.output(dbrda.glob, 
               file = "stats/permtest_global_mod1_18S.txt",
               append = TRUE)

# test done using 999 permutations with remaining variables as covariates
# + Benjamini-Hochberg correction

set.seed(1)
margin.vifmod.dbrda <- anova(vifmod.dbrda, by = 'margin', parallel = 4)
margin.vifmod.dbrda$`Pr(>F)` <- p.adjust(margin.vifmod.dbrda$`Pr(>F)`, method = 'BH')

# save permutation test results

cat(" -------------------------------------------------------------------",
    "\n", "18S - model 1 - permutation test with remaining vars as covariates", "\n",
    "-------------------------------------------------------------------",
    "\n", "\n",
    file = "tables/permtest_covariates_mod1_18S.txt")

capture.output(margin.vifmod.dbrda, 
               file = "tables/permtest_covariates_mod1_18S.txt",
               append = TRUE)

# db-RDA pt 7: perform db-RDA + envit for monitoring-based model ####

# db-RDA
set.seed(1)
vifmod.dbrda.mon <- ordinate(rawdata.2.clr, 
                             method = "CAP", 
                             distance = "euclidean",
                             formula = ~ C_1cm + O2_BW + NOX_BW + P_BW)

# envfit (as a data frame using fortify)
# calculated on linear combinations of scores
# see: https://www.davidzeleny.net/anadat-r/doku.php/en:confusions

set.seed(1)
envfit.lc.mon <- fortify(
   envfit(vifmod.dbrda.mon ~ C_1cm + O2_BW + NOX_BW + P_BW,
          data = rawdata.2.clr.samdata,
          permutations = 999, display = 'lc')
)

# db-RDA pt 8: plot db-RDA + envit for monitoring-based model ####

# Note: arbitrary multiplier addded to xend + yend
# for plotting purposes (could also be handled through scaling)

Cairo(file = "figures/r_output/Fig7d_18S.tiff", 
      type = "tiff", 
      units = "cm", 
      width = 20, 
      height = 20, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

vifmod.dbrda.mon.plot <- plot_ordination(rawdata.2.clr, 
                                         vifmod.dbrda.mon, 
                                         color = "Site") +
   theme(aspect.ratio = 1) +
   scale_colour_brewer(type = "qual", 
                       palette = "Dark2") + 
   geom_point(size = 6.5, alpha = 0.75) +
   geom_segment(data = envfit.lc.mon,
                aes(x = 0, 
                    y = 0, 
                    xend = CAP1*6, 
                    yend = CAP2*6), 
                size = 1, 
                alpha = 0.5, 
                colour = "grey30") +
   geom_label(data = envfit.lc.mon, 
              aes(x = CAP1*6, 
                  y = CAP2*6),
              label = envfit.lc.mon$Label, 
              colour = "black", 
              fontface = "bold",
              size = 3.5) + 
   theme(axis.title = element_blank()) +
   theme(axis.ticks = element_blank()) +
   theme(axis.text = element_blank()) +
   theme(legend.position = "none") +
   expand_limits(x = 5.7) # extra line to make sure 02_BW label included in full

vifmod.dbrda.mon.plot
dev.off()

# db-RDA pt 9: permutation test for db-RDA explanatory variables (monitoring-based model) ####

# global test (999 permutations)
dbrda.glob.mon <- anova(vifmod.dbrda.mon, permutations = 999)

# save test results

cat(" -------------------------------------------------------------------",
    "\n", "18S - model 2 - global test", "\n",
    "-------------------------------------------------------------------",
    "\n", "\n",
    file = "tables/permtest_global_mod2_18S.txt")

capture.output(dbrda.glob.mon, 
               file = "stats/permtest_global_mod2_18S.txt",
               append = TRUE)

set.seed(1)
margin.vifmod.dbrda.mon <- anova(vifmod.dbrda.mon, by = 'margin', parallel = 4)
margin.vifmod.dbrda.mon$`Pr(>F)` <- p.adjust(margin.vifmod.dbrda.mon$`Pr(>F)`, method = 'BH')

# save permutation test results

cat(" -------------------------------------------------------------------",
    "\n", "18S - model 2 - permutation test with remaining vars as covariates", "\n",
    "-------------------------------------------------------------------",
    "\n", "\n",
    file = "tables/permtest_covariates_mod2_18S.txt")

capture.output(margin.vifmod.dbrda.mon, 
               file = "tables/permtest_covariates_mod2_18S.txt",
               append = TRUE)
