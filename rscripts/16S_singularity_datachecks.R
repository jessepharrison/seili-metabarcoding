# haverö metabarcoding study - 16S analysis (data checks)
# jesse harrison 2020-2021

# additional libpath ####
# (see extra_RPackages.R for extra package installs)

.libPaths(c("/home/jharriso/seili-singularity/rpackages", .libPaths()))

# packages ####

packages <- c("phyloseq", "ggplot2", "vegan", "grid", "gridExtra", "data.table", "plyr", "cowplot", 
              "RVAideMemoire", "microbiome", "knitr", "RColorBrewer", "Cairo", "multcompView", 
              "QsRutils", "dplyr", "ggvegan")

lapply(packages, require, character.only = TRUE)

# ggplot2 theme ####

theme_set(theme_classic())

# disable scientific notation
options(scipen=10000)

# working directory ####

setwd("/home/jharriso/git/seili-metabarcoding/")

# load 16S RData ####

load("rdata/Seili16s.RData")
# based on: import_biom("tables_sampledata_export.biom", treefilename="tree_rooted.tree", refseqfilename="otus.fasta")
# 'rawdata' is a version of the dataset where no step has been taken to remove NAs + chloroplast / mitochondrial sequences

# dataset with NAs removed at phylum level, retaining potential chloroplast / mitochondrial sequencies ####

rawdata.bac.noNA <- subset_taxa(rawdata, 
                            !is.na(Rank2) & 
                              !Rank2 %in% c("", "uncharacterized"))

# bacterial dataset with NAs removed at phylum level and with class Chloroplast / family Mitochondria removed ####

rawdata.bac.noNA.noCyano <- subset_taxa(rawdata.bac.noNA, 
                                    Rank3 != "Chloroplast" & # Class Chloroplast removed
                                      Rank5 != "Mitochondria") # Family Mitochondria removed

# removal of Thermotogae and Chlamydiae from all datasets ####

# phyla Thermotogae + Chlamydiae were removed from analyses due to low prevalences (see 16S script).

filterPhyla.raw = c("Chlamydiae", "Thermotogae")

rawdata.bac.1 <- subset_taxa(rawdata, !Rank2 %in% filterPhyla.raw)
rawdata.bac.noNA.1 <- subset_taxa(rawdata.bac.noNA, !Rank2 %in% filterPhyla.raw)
rawdata.bac.noNA.noCyano.1 <- subset_taxa(rawdata.bac.noNA.noCyano, !Rank2 %in% filterPhyla.raw)

# further pre-processing (prevalences) ####

# source script with fast_melt() function
# ("melts" OTU table into format with three main columns: SampleID, TaxaID and count)
# taxa_summary.R can be downloaded from: http://evomics.org/phyloseq/taxa_summary-r/

source("rscripts/taxa_summary.R", local = TRUE)

# pre-processing following steps in analysis script

# fast_melt + simplifying names
mdt.raw <- fast_melt(rawdata.bac.1)
mdt.raw2 <- fast_melt(rawdata.bac.noNA.1)
mdt.raw3 <- fast_melt(rawdata.bac.noNA.noCyano.1)

# calculating relative abundances
mdt.raw[, RelativeAbundance := count / sum(count), by = SampleID]
prevdt.raw = mdt.raw[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

mdt.raw2[, RelativeAbundance := count / sum(count), by = SampleID]
prevdt.raw2 = mdt.raw2[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

mdt.raw3[, RelativeAbundance := count / sum(count), by = SampleID]
prevdt.raw3 = mdt.raw3[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

# prevalence calculations
prevdf.raw = apply(X = otu_table(rawdata.bac.1),
                   MARGIN = ifelse(taxa_are_rows(rawdata.bac.1), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

prevdf.raw = data.frame(Prevalence = prevdf.raw,
                        TotalAbundance = taxa_sums(rawdata.bac.1),
                        tax_table(rawdata.bac.1))

prevdf.raw2 = apply(X = otu_table(rawdata.bac.noNA.1),
                    MARGIN = ifelse(taxa_are_rows(rawdata.bac.noNA.1), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

prevdf.raw2 = data.frame(Prevalence = prevdf.raw2,
                         TotalAbundance = taxa_sums(rawdata.bac.noNA.1),
                         tax_table(rawdata.bac.noNA.1))

prevdf.raw3 = apply(X = otu_table(rawdata.bac.noNA.noCyano.1),
                    MARGIN = ifelse(taxa_are_rows(rawdata.bac.noNA.noCyano.1), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

prevdf.raw3 = data.frame(Prevalence = prevdf.raw3,
                         TotalAbundance = taxa_sums(rawdata.bac.noNA.noCyano.1),
                         tax_table(rawdata.bac.noNA.noCyano.1))

# phylum-level prevalences
prevdf.raw.phyla = subset(prevdf.raw, Rank2 %in% get_taxa_unique(rawdata.bac.1, "Rank2"))
prevdf.raw.phyla2 = subset(prevdf.raw2, Rank2 %in% get_taxa_unique(rawdata.bac.noNA.1, "Rank2"))
prevdf.raw.phyla3 = subset(prevdf.raw3, Rank2 %in% get_taxa_unique(rawdata.bac.noNA.noCyano.1, "Rank2"))

# 5% prevalence filtering ####

# defining 5% prevalence thresholds for filtering
prevalenceThreshold.raw <- 0.05 * nsamples(rawdata.bac.1)
prevalenceThreshold.raw2 <- 0.05 * nsamples(rawdata.bac.noNA.1)
prevalenceThreshold.raw3 <- 0.05 * nsamples(rawdata.bac.noNA.noCyano.1)

# execute prevalence filter using `prune_taxa()` function
keepTaxa.raw <- rownames(prevdf.raw.phyla)[(prevdf.raw.phyla$Prevalence >= prevalenceThreshold.raw)]
rawdata.2 <- prune_taxa(keepTaxa.raw, rawdata.bac.1)

keepTaxa.raw2 <- rownames(prevdf.raw.phyla2)[(prevdf.raw.phyla2$Prevalence >= prevalenceThreshold.raw2)]
rawdata.noNA.2 <- prune_taxa(keepTaxa.raw2, rawdata.bac.noNA.1)

keepTaxa.raw3 <- rownames(prevdf.raw.phyla3)[(prevdf.raw.phyla3$Prevalence >= prevalenceThreshold.raw3)]
rawdata.noNA.noCyano.2 <- prune_taxa(keepTaxa.raw3, rawdata.bac.noNA.noCyano.1)

# %RA transformation ####

# taxon composition plots showed as proportions (%) for readability
# note: here not using CLR-transformed data

raw2.ra <- transform_sample_counts(rawdata.2, 
                                   function(x) x/sum(x))

raw2.noNA.ra <- transform_sample_counts(rawdata.noNA.2, 
                                        function(x) x/sum(x))

raw2.noNA.noCyano.ra <- transform_sample_counts(rawdata.noNA.noCyano.2, 
                                        function(x) x/sum(x))

# stacked bar plot using %RA data ####

colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
             "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
             "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
             "#771122", "#AA4455", "#DD7788", "#DD776655")

# subset the data to domain level (= rank1)
raw2.ra.r1 <- tax_glom(raw2.ra, taxrank = rank_names(raw2.ra)[1])
raw2.ra.r1 <- psmelt(raw2.ra.r1)

Cairo(file = "figures/r_output/FigS2a_16S.tiff", 
      type = "tiff", 
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
  theme(axis.title.y = element_text(size = 14)) +
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
rawdata.noNA.noCyano.2.clr <- transform(rawdata.noNA.noCyano.2, 'clr')
# nMDS using CLR-transformed data ####

set.seed(1)
# using Aitchison distance (CLR + Euclidean distances)
rawdata.2.nmds <- ordinate(physeq = rawdata.2.clr, 
                           method = "NMDS", 
                           distance = "euclidean") # stress = 0.04

set.seed(1)
rawdata.noNA.2.nmds <- ordinate(physeq = rawdata.noNA.2.clr, 
                               method = "NMDS", 
                               distance = "euclidean") # stress = 0.04

set.seed(1)
rawdata.noNA.noCyano.2.nmds <- ordinate(physeq = rawdata.noNA.noCyano.2.clr, 
                                        method = "NMDS", 
                                        distance = "euclidean") # stress = 0.04


# nMDS plots ####

Cairo(file = "figures/r_output/FigS3a_16S.tiff", 
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
  theme(legend.position = c(0.1, 0.8),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 13))  

nMDS.plot.rawdata.2.clr
dev.off()

Cairo(file = "figures/r_output/FigS3b_16S.tiff", 
      type = "tiff", 
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

Cairo(file = "figures/r_output/FigS3c_16S.tiff", 
      type = "tiff", 
      units = "cm", 
      width = 15, 
      height = 15, 
      pointsize = 14, 
      dpi = 300, 
      bg= "white")

nMDS.plot.rawdata.noNA.noCyano.2.clr <- plot_ordination(rawdata.noNA.noCyano.2.clr, 
                                                        rawdata.noNA.noCyano.2.nmds, 
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

nMDS.plot.rawdata.noNA.noCyano.2.clr
dev.off()
