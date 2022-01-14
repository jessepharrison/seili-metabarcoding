# haverö metabarcoding study - 18S analysis (data checks)
# jesse harrison 2020-2022

# make sure your working directory is the github repository root!
# (i.e. path/to/seili-metabarcoding)

# add singularity folder to libpaths

.libPaths(c("singularity", .libPaths()))

# packages ####

packages <- c("phyloseq", "ggplot2", "vegan", "grid", "gridExtra", "data.table", "plyr", "cowplot", 
              "RVAideMemoire", "microbiome", "knitr", "RColorBrewer", "Cairo", "multcompView", 
              "QsRutils", "dplyr", "ggvegan")

lapply(packages, require, character.only = TRUE)

# ggplot2 theme ####

theme_set(theme_classic())

# disable scientific notation
options(scipen=10000)

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

# %RA transformation ####

# taxon composition plots showed as proportions (%) for readability
# note: here not using CLR-transformed data

raw2.ra <- transform_sample_counts(rawdata.2, 
                                   function(x) x/sum(x))

raw2.noNA.ra <- transform_sample_counts(rawdata.noNA.2, 
                                        function(x) x/sum(x))

# fractions of classified vs unclassified ####

# subset the data to domain level (= rank1)
raw2.ra.r1 <- tax_glom(raw2.ra, taxrank = rank_names(raw2.ra)[1])
raw2.ra.r1 <- psmelt(raw2.ra.r1)
raw2.ra.r1

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

Cairo(file = "figures/r_output/FigS2d_18S.tiff", 
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

Cairo(file = "figures/r_output/FigS2e_18S.tiff", 
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
