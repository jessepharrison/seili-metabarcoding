## Seili Data - Phyloseq Analysis (Bacteria) ######################
## Jesse Harrison (June 2018) ####################################

## CLEAN-UP, LIBRARIES, PLOT THEME ################################

# phyloseq installation:

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("phyloseq")

# microbiome installation:
# BiocManager::install("microbiome")

# clean up R
rm(list=ls())

# libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(grid)
library(gridExtra)
library(data.table)
library(plyr)
library(cowplot)
library(RVAideMemoire)
library(microbiome)
library(knitr)
library(RColorBrewer)
library(Cairo)
library(multcompView)

# other attached packages:
# [1] multcompView_0.1-7   Cairo_1.5-10         RColorBrewer_1.1-2   knitr_1.25           microbiome_1.6.0     RVAideMemoire_0.9-73 cowplot_1.0.0       
# [8] plyr_1.8.4           data.table_1.12.6    gridExtra_2.3        vegan_2.5-6          lattice_0.20-38      permute_0.9-5        ggplot2_3.2.1       
# [15] phyloseq_1.28.0 

# theme for plots
theme_set(theme_classic())

## WORKING DIR + DATA IMPORT ######################################

# import from rdata files

setwd("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/RData/")
load("Seili16s.RData")

# data import from source files

# set working directory (based on location of micca output)
# setwd("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharr") # primer sequences removed in micca

# use import_biom() to import BIOM file, phylogenetic tree file and reference sequence file
# rawdata <- import_biom("tables_sampledata_export.biom", treefilename="tree_rooted.tree", refseqfilename="otus.fasta")
# the 'export' version of the BIOM contains minimal metadata - if want more, remove '_export' from file name

## RAREFACTION WITHOUT EQUAL SAMPLING DEPTH #######################

# Cairo(file = "FigS1_Bac_Winter20.png", type = "png", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")
rarecurve(t(otu_table(rawdata)), step = 100, 
          cex.lab = 1.5, cex.axis = 1.5, label = FALSE, ylab = "Phylotype")
# dev.off()

## ALPHA DIVERSITY ###

# Cairo(file = "Seili_16S_alpha_rawdata_Winter20.png", type = "png", units = "cm", width = 22, height = 12, pointsize = 14, dpi = 300, bg= "white")

alphaplot <- plot_richness(rawdata, x = "Site", measures = c("Observed", "Chao1", "Shannon"))
alphaplot <- alphaplot +  
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust = 0.5)) +
  xlab("\nSite") +
  ylab("Alpha Diversity Measure\n")

alphaplot + geom_boxplot(data = alphaplot$data, 
                         aes(x = Site, color = NULL), alpha = 0.1)
# dev.off()

estimate_richness(rawdata, measures=c("Observed", "Chao1", "Shannon"))
evenness(rawdata, 'pielou')

# alpha div for rarefied raw data

set.seed(1)
alpharare <- rarefy_even_depth(rawdata)

estimate_richness(alpharare, measures=c("Observed", "Chao1", "Shannon"))
evenrare <- evenness(alpharare, 'pielou')

## INSPECTING THE DATA AND FILTERING OUT ARTIFACTS / LOW-PREVALENCE PHYLA ######

# basic commands
nsamples(rawdata) ## [1] 30
ntaxa(rawdata) ## [1] 15394

sample_sums(rawdata)

sample_variables(rawdata)

table(tax_table(rawdata)[, "Rank2"], exclude = NULL)

#          Acidobacteria              Actinobacteria               Aminicenantes             Armatimonadetes 
#                    415                         304                           3                          54 
#          Bacteroidetes                        BRC1    candidate division WPS-1    candidate division WPS-2 
#                   1540                          47                          11                           1 
# candidate division ZB3 Candidatus Saccharibacteria                  Chlamydiae                 Chloroflexi 
#                     16                          19                           1                         399 
#          Cloacimonetes   Cyanobacteria/Chloroplast             Deferribacteres         Deinococcus-Thermus 
#                     14                          72                          23                           4 
#          Elusimicrobia               Fibrobacteres                  Firmicutes                Fusobacteria 
#                     27                          48                         224                          15 
#       Gemmatimonadetes             Hydrogenedentes             Ignavibacteriae             Latescibacteria 
#                     20                         126                          39                          55 
#          Lentisphaerae              Marinimicrobia              Microgenomates                 Nitrospinae 
#                    234                           2                           2                          10 
#            Nitrospirae               Parcubacteria              Planctomycetes              Proteobacteria 
#                     15                         114                         386                        4654 
#           Spirochaetes                         SR1               Synergistetes                 Tenericutes 
#                    173                          47                           2                           5 
#            Thermotogae             Verrucomicrobia                        <NA> 
#                      1                         237                        6035 

# a few phyla that occur only once: candidate division WPS-2, Chlamydiae, Thermotogae
# there are also 6035 NAs that are likely to be artifacts and should be removed.

# rawdata.psmelt <- psmelt(rawdata)
# View(rawdata.psmelt)

# first let's filter out any mitochondrial + chloroplast sequences (focusing on bacteria only)
rawdata.bac <- subset_taxa(rawdata, Rank1 == "Bacteria" & Rank3 != "Chloroplast")
table(tax_table(rawdata.bac)[, "Rank2"], exclude = NULL)

ntaxa(rawdata.bac) # [1] 7645

# rawdata.bac.psmelt <- psmelt(rawdata.bac)
# View(rawdata.bac.psmelt)

# following workflow based on: https://f1000research.com/articles/5-1492/v2

rawdata.bac.0 <- subset_taxa(rawdata.bac, !is.na(Rank2) & !Rank2 %in% c("", "uncharacterized"))
table(tax_table(rawdata.bac.0)[, "Rank2"], exclude = NULL)

# this doesn't do anything as the previous filter step already removed NAs, but keeping it for sake of
# consistency between scripts.

# next, we can inspect the data in more detail to see if certain phyla are comprised mostly of low-prevalence
# features (prevalence = no. of times an OTU is observed at least once).

# compute prevalence of each feature, store as data.frame
prevdf.raw = apply(X = otu_table(rawdata.bac.0),
                   MARGIN = ifelse(taxa_are_rows(rawdata.bac.0), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

# add taxonomy and total read counts to this data.frame
prevdf.raw = data.frame(Prevalence = prevdf.raw,
                        TotalAbundance = taxa_sums(rawdata.bac.0),
                        tax_table(rawdata.bac.0))

plyr::ddply(prevdf.raw, "Rank2", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#                          Rank2         1     2
# 1                Acidobacteria 15.965517  6482
# 2               Actinobacteria 16.347518  4610
# 3                Aminicenantes 11.666667    35
# 4              Armatimonadetes 11.640000   582
# 5                Bacteroidetes 14.757534 10773
# 6                         BRC1  9.085106   427
# 7     candidate division WPS-1  6.000000    66
# 8     candidate division WPS-2 23.000000    23
# 9       candidate division ZB3  5.000000    80
# 10 Candidatus Saccharibacteria  7.631579   145
# 11                  Chlamydiae  1.000000     1
# 12                 Chloroflexi 16.614706  5649
# 13               Cloacimonetes 13.571429   190
# 14   Cyanobacteria/Chloroplast 10.250000    41
# 15             Deferribacteres 21.304348   490
# 16         Deinococcus-Thermus 18.750000    75
# 17               Elusimicrobia  8.280000   207
# 18               Fibrobacteres 10.583333   508
# 19                  Firmicutes 10.475962  2179
# 20                Fusobacteria 11.733333   176
# 21            Gemmatimonadetes 22.950000   459
# 22             Hydrogenedentes 10.666667  1344
# 23             Ignavibacteriae 19.666667   767
# 24             Latescibacteria 16.036364   882
# 25               Lentisphaerae 12.027027  1335
# 26              Marinimicrobia 25.000000    50
# 27              Microgenomates  9.000000    18
# 28                 Nitrospinae 14.700000   147
# 29                 Nitrospirae 23.933333   359
# 30               Parcubacteria  4.657895   531
# 31              Planctomycetes  9.500000  3591
# 32              Proteobacteria 11.974371 49059
# 33                Spirochaetes 11.658960  2017
# 34                         SR1  6.659574   313
# 35               Synergistetes  5.500000    11
# 36                 Tenericutes 15.000000    75
# 37                 Thermotogae  4.000000     4
# 38             Verrucomicrobia 11.432836  2298

# the phyla that occurred only once are also of very low prevalence (Chlamydiae: 1, Thermotogae: 4, WPS-2: 23).
# out of these, would filter Chlamydiae and Thermotogae.

# define phyla to filter and remove them from the data set
filterPhyla.raw = c("Chlamydiae", "Thermotogae")
rawdata.1 = subset_taxa(rawdata.bac.0, !Rank2 %in% filterPhyla.raw)
rawdata.1

table(tax_table(rawdata.1)[, "Rank2"], exclude = NULL)

## ADDITIONAL WAYS TO INSPECT PREVALENCE / OTU COUNTS #############################################################

# source script with fast_melt() function
# ("melts" OTU table into format with three main columns: SampleID, TaxaID and count)
# taxa_summary.R can be downloaded from: http://evomics.org/phyloseq/taxa_summary-r/

source("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/Scripts/taxa_summary.R", local = TRUE)

mdt.raw = fast_melt(rawdata.1)

# omit NAs
mdt.raw <- mdt.raw[!is.na(count)]

# calculate relative abundance
mdt.raw[, RelativeAbundance := count / sum(count), by = SampleID]
prevdt.raw = mdt.raw[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

# draw prevalence plot
ggplot(prevdt.raw, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Taxon prevalence (uneven)")

# how many with 0 seqs?
prevdt.raw[(Prevalence <= 0), .N] ## [1] 0

# how many singletons?
prevdt.raw[(Prevalence <= 1), .N] ## [1] 298

# how many doubletons?
prevdt.raw[(Prevalence <= 2), .N] ## [1] 725

# taxa cumulative sum with prevalence on x axis
prevcumsum.raw = prevdt.raw[, .N, by = Prevalence]
setkey(prevcumsum.raw, Prevalence)
prevcumsum.raw[, CumSum := cumsum(N)]
pPrevCumSum.raw = ggplot(prevcumsum.raw, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold (uneven)")
pPrevCumSum.raw

# prevalence vs. total count scatter plot
ggplot(prevdt.raw, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

# prevalence plot for phyla

# Cairo(file = "Fig_S4_Bac_Winter20.png", type = "png", units = "cm", width = 42, height = 30, pointsize = 14, dpi = 300, bg= "white")
prevdf.raw.phyla = subset(prevdf.raw, Rank2 %in% get_taxa_unique(rawdata.1, "Rank2"))

ggplot(prevdf.raw.phyla, aes(TotalAbundance, Prevalence / nsamples(rawdata.bac.0),color = Rank2)) +
  # Include a guess for parameter (5 % of total samples)
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Rank2) + theme(legend.position="none") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1)) +
  xlab("\nTotal Abundance") +
  ylab("Prevalence [Frac. Samples]\n")

# dev.off()

# BEFORE RUNNING PREVALENCE FILTER: ALSO CHECK DATA WITHOUT SINGLETONS / DOUBLETONS

# defining 5% prevalence threshold for filtering
prevalenceThreshold.raw = 0.05 * nsamples(rawdata.bac.0)

# execute prevalence filter, using `prune_taxa()` function
keepTaxa.raw = rownames(prevdf.raw.phyla)[(prevdf.raw.phyla$Prevalence >= prevalenceThreshold.raw)]
rawdata.2 = prune_taxa(keepTaxa.raw, rawdata.bac.0)

# Cairo(file = "Seili_16S_rare_afterprevfiltNEW_Winter20.png", type = "png", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg = "white")
rarecurve(t(otu_table(rawdata.2)), step = 100,
          cex.lab = 1.5, cex.axis = 1.5, label = FALSE, ylab = "Phylotype")
# dev.off()

## EXAMINING DATA WITHOUT SINGLETONS / DOUBLETONS

rawdata.f <- filter_taxa(rawdata.bac, function (x) {sum(x > 0) > 2}, prune = TRUE) # removes singletons and doubletons

rarecurve(t(otu_table(rawdata.f)), step = 100)

table(tax_table(rawdata.f)[, "Rank2"], exclude = NULL)
# single-occurrence phyla: candidate division WPS-2, Synergistetes, Thermotogae
# NAs: none!
ntaxa(rawdata.f) # 6919; no need to create rawdata.f.0 because no NAs

prevdf.raw.f = apply(X = otu_table(rawdata.f),
                     MARGIN = ifelse(taxa_are_rows(rawdata.f), yes = 1, no = 2),
                     FUN = function(x){sum(x > 0)})

prevdf.raw.f = data.frame(Prevalence = prevdf.raw.f,
                          TotalAbundance = taxa_sums(rawdata.f),
                          tax_table(rawdata.f))

plyr::ddply(prevdf.raw.f, "Rank2", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# WPS-2: 23, Synergistetes: 10, Thermotogae: 4

# define phyla to filter and remove them from the data set. based on above, would remove Thermotogae
filterPhyla.raw.f = c("Thermotogae")
rawdata.1.f = subset_taxa(rawdata.f, !Rank2 %in% filterPhyla.raw.f)
rawdata.1.f

table(tax_table(rawdata.1.f)[, "Rank2"], exclude = NULL)

mdt.raw.f = fast_melt(rawdata.1.f)
mdt.raw.f <- mdt.raw.f[!is.na(count)]
mdt.raw.f[, RelativeAbundance := count / sum(count), by = SampleID]
prevdt.raw.f = mdt.raw.f[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

ggplot(prevdt.raw.f, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Taxon prevalence (uneven)")

prevcumsum.raw.f = prevdt.raw.f[, .N, by = Prevalence]
setkey(prevcumsum.raw.f, Prevalence)
prevcumsum.raw.f[, CumSum := cumsum(N)]
pPrevCumSum.raw.f = ggplot(prevcumsum.raw.f, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold (uneven)")
pPrevCumSum.raw.f

ggplot(prevdt.raw.f, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

prevdf.raw.phyla.f = subset(prevdf.raw.f, Rank2 %in% get_taxa_unique(rawdata.1.f, "Rank2"))
ggplot(prevdf.raw.phyla.f, aes(TotalAbundance, Prevalence / nsamples(rawdata.1.f),color = Rank2)) +
  # Include a guess for parameter (5 % of total samples)
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Rank2) + theme(legend.position="none")

# Cairo(file = "Seili_16S_rare_rawdata_nosingleNEW_winter20.png", type = "png", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")
rarecurve(t(otu_table(rawdata.1.f)), step = 100,
          cex.lab = 1.5, cex.axis = 1.5, label = FALSE, ylab = "Phylotype")
# dev.off()

# results are pretty much identical to what get with a 5% prevalence filter.
# however, prevalence filter is slightly less stringent with the 18S data.
# using prevalence filter for both data sets.

## HELLINGER TRANSFORMATION + RA CONVERSION ############################################

# Hellinger transformation for data set using 5% prevalence filter (rawdata.2)
# (= square root of relative abundance)

# sample_sums(rawdata.2)

rawdata.2.hellinger <- rawdata.2
otu_table(rawdata.2.hellinger) <-otu_table(decostand(otu_table(rawdata.2.hellinger), 
                                                     method = "hellinger"), taxa_are_rows = TRUE)

rawdata.2.hellinger <- transform_sample_counts(rawdata.2.hellinger, function(x) x/sum(x)) # convert to proportions

# ntaxa(rawdata.2.hellinger)

# export the data as a CSV file
# rawdata.2.hellinger.psmelt <- psmelt(rawdata.2.hellinger) # melt to long format (also converts to data frame)
# fwrite(rawdata.2.hellinger.psmelt, "Seili_16S_Hellinger_RA_noprim_Aug18.csv")

# if want to view in table format inside R: 
# View(rawdata.2.hellinger.psmelt)

# RA without Hellinger transformation? 
# rawdata.2.nohellinger <- transform_sample_counts(rawdata.2, function(x) x/sum(x)) # convert to proportions

## TAXON BAR PLOTS #####################################################################

# subset the data to phylum level (= rank2), cut out low-abundance taxa

raw2.ra.phylum <- tax_glom(rawdata.2.hellinger, taxrank = rank_names(rawdata.2.hellinger)[2], NArm = T,
                           bad_empty = c(NA,""," ","\t")) # filters out any junk (not necessarily required but keeping it there anyway)

raw2.ra.phylum <- psmelt(raw2.ra.phylum) # melt to long format (also converts to data frame)
raw2.ra.phylum.02 <- subset(raw2.ra.phylum, Abundance > 0.02) # cut out OTUs less than 2%
sum(raw2.ra.phylum.02$Abundance)/sum(raw2.ra.phylum$Abundance) # retains 88.7% of original data

# View(raw2.ra.phylum)

# choose a colour palette
colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
             "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",
             "#CBD588")

# Cairo(file = "Fig5_Bac_Phyla_Winter20.png", type = "png", units = "cm", width = 23, height = 20, pointsize = 14, dpi = 300, bg= "white")
phylum.plot <- ggplot(raw2.ra.phylum.02 , aes(x = Replicate, y = Abundance, fill = Rank2)) + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >2%)\n") +
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))
phylum.plot
# dev.off()

# proteobacteria are the majority everywhere
# ... but it also looks like there's an increase in bacteroidetes in S1 (and maybe S2+3)

# class level plot (= rank3) with 3% cut-off

raw2.ra.class <- tax_glom(rawdata.2.hellinger, taxrank = rank_names(rawdata.2.hellinger)[3], NArm = T,
                          bad_empty = c(NA,""," ","\t"))

raw2.ra.class <- psmelt(raw2.ra.class)
raw2.ra.class.03 <- subset(raw2.ra.class, Abundance > 0.03) # cut out OTUs less than 3%
sum(raw2.ra.class.03$Abundance)/sum(raw2.ra.class$Abundance) # retains 65% of original data

# Cairo(file = "FigS5_Bac_Class_Winter20.png", type = "png", units = "cm", width = 23, height = 20, pointsize = 14, dpi = 300, bg= "white")
class.plot.03 <- ggplot(raw2.ra.class.03, aes(x = Replicate, y = Abundance, fill = Rank3)) + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >3%)\n") +
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))
class.plot.03
# dev.off()

## nMDS ORDINATION ######################################################################

set.seed(1)

rawdata.2.nmds <- ordinate(physeq = rawdata.2.hellinger, method = "NMDS", distance = "bray") 
# stress = 0.05

# Cairo(file = "Fig3_Bac_Winter20.png", type = "png", units = "cm", width = 15, height = 15, pointsize = 14, dpi = 300, bg= "white")
nMDS.plot.rawdata.2.hellinger <- plot_ordination(rawdata.2.hellinger, rawdata.2.nmds, color = "Site") + theme(aspect.ratio = 1)
nMDS.plot.rawdata.2.hellinger <- nMDS.plot.rawdata.2.hellinger + geom_point(size = 6.5, alpha = 0.75)
nMDS.plot.rawdata.2.hellinger <- nMDS.plot.rawdata.2.hellinger + scale_colour_brewer(type = "qual", palette = "Dark2")
nMDS.plot.rawdata.2.hellinger <- nMDS.plot.rawdata.2.hellinger + theme(axis.title = element_blank())
nMDS.plot.rawdata.2.hellinger <- nMDS.plot.rawdata.2.hellinger + theme(axis.ticks = element_blank())
nMDS.plot.rawdata.2.hellinger <- nMDS.plot.rawdata.2.hellinger + theme(axis.text = element_blank())
nMDS.plot.rawdata.2.hellinger <- nMDS.plot.rawdata.2.hellinger + theme(legend.title = element_blank())
nMDS.plot.rawdata.2.hellinger <- nMDS.plot.rawdata.2.hellinger + theme(legend.position = "bottom")
nMDS.plot.rawdata.2.hellinger
# dev.off()

# set.seed(1)

# rawdata.2.nmds.nohellinger <- ordinate(physeq = rawdata.2.nohellinger, method = "NMDS", distance = "bray") 

# nMDS.plot.rawdata.2.nohellinger <- plot_ordination(rawdata.2.nohellinger, rawdata.2.nmds.nohellinger, color = "Site") + theme(aspect.ratio = 1)
# nMDS.plot.rawdata.2.nohellinger <- nMDS.plot.rawdata.2.nohellinger + geom_point(size = 6.5, alpha = 0.75)
# nMDS.plot.rawdata.2.nohellinger <- nMDS.plot.rawdata.2.nohellinger + scale_colour_brewer(type = "qual", palette = "Dark2")
# nMDS.plot.rawdata.2.nohellinger <- nMDS.plot.rawdata.2.nohellinger + theme(axis.title = element_blank())
# nMDS.plot.rawdata.2.nohellinger <- nMDS.plot.rawdata.2.nohellinger + theme(axis.ticks = element_blank())
# nMDS.plot.rawdata.2.nohellinger <- nMDS.plot.rawdata.2.nohellinger + theme(axis.text = element_blank())
# nMDS.plot.rawdata.2.nohellinger <- nMDS.plot.rawdata.2.nohellinger + theme(legend.title = element_blank())
# nMDS.plot.rawdata.2.nohellinger <- nMDS.plot.rawdata.2.nohellinger + theme(legend.position = "bottom")
# nMDS.plot.rawdata.2.nohellinger

# using Hellinger transformation gives a slight improvement in terms of separating S2 and S3

## PERMANOVA FOR FACTOR SITE ############################################################

## ADONIS() (PERMANOVA global test, Bray-Curtis distance matrix)

# make a data frame from the sample_data
rawdata.2.hellinger.df <- data.frame(sample_data(rawdata.2.hellinger))

# extract bray-curtis distances

### NOTE MAY 2018: occasional compatibility issues between vegan and phyloseq.
### if need an older version of vegan, then: 
### install_version("vegan", version = "2.4-6", repos = "http://cran.us.r-project.org")

set.seed(1)
rawdata.2.hellinger.bray <- phyloseq::distance(rawdata.2.hellinger, method = 'bray')

# do the test
set.seed(1)
adonis(rawdata.2.hellinger.bray ~ Site, data = rawdata.2.hellinger.df)

# Call:
# adonis(formula = rawdata.2.hellinger.bray ~ Site, data = rawdata.2.hellinger.df) 

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs  F.Model     R2  Pr(>F)    
# Site       5    2.6750 0.53500  5.5209 0.53492  0.001 ***
# Residuals 24    2.3257 0.09691         0.46508           
# Total     29    5.0008                 1.00000          
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# post-hoc pairwise comparison

set.seed(1)
pairwise.perm.manova(rawdata.2.hellinger.bray, rawdata.2.hellinger.df$Site, nperm = 999, p.method = "BH")

# Pairwise comparisons using permutation MANOVAs on a distance matrix 

# data:  rawdata.2.hellinger.bray by rawdata.2.hellinger.df$Site
# 999 permutations 

#    S1    S2    S3    S4    S5   
# S2 0.014 -     -     -     -    
# S3 0.014 0.014 -     -     -    
# S4 0.014 0.014 0.014 -     -    
# S5 0.014 0.014 0.014 0.014 -    
# S6 0.014 0.014 0.014 0.014 0.014

# P value adjustment method: BH 

## Betadisper test

set.seed(1)
beta <- betadisper(rawdata.2.hellinger.bray, rawdata.2.hellinger.df$Site)
set.seed(1)
permutest(beta)

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999

# Response: Distances
# Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)   
# Groups     5 0.0083936 0.00167871 3.9908    999  0.008 **
# Residuals 24 0.0100955 0.00042065                        
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# as with 18S, there is a significant multivariate dispersion effect that needs to be taken into
# account when interpreting the PERMANOVA output.

# pairwise Tukey HSD

betaHSD <- TukeyHSD(beta)

# Cairo(file = "Seili_16S_Tukey1NEW_Winter20.png", type = "png", units = "cm", width = 15, height = 15, pointsize = 12, dpi = 300, bg= "white")
psig = as.numeric(apply(betaHSD$`group`[,2:3],1,prod)>=0)+1
op = par(mar=c(4.2,5,3.8,2))
plot(betaHSD, col = psig, yaxt = "n")
for (j in 1:length(psig)){
  axis(2,at=j,labels=rownames(betaHSD$`group`)[length(psig)-j+1],
       las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
}
par(op)
# dev.off()

# we can see from here that the differences in dispersion occur when comparing 
# S6 with S1, S3 or S4.

# or can extract letters like this:
multcompLetters(extract_p(betaHSD$group))
#  S2   S3   S4   S5   S6   S1 
# "ab"  "a"  "a" "ab"  "b"  "a"

# box plot with letters

# Cairo(file = "Fig4_Bac_Winter20.png", type = "png", units = "cm", width = 15, height = 15, pointsize = 12, dpi = 300, bg= "white")
boxplot(beta, xlab = "Site", ylim=c(0.23, 0.39))
text(x= 1, y= 0.284, labels= "a")
text(x= 2, y= 0.326, labels= "ab")
text(x= 3, y= 0.284, labels= "a")
text(x= 4, y= 0.285, labels= "a")
text(x= 5, y= 0.305, labels= "ab")
text(x= 6, y= 0.380, labels= "b")
# dev.off()

## SAVE THE DATA ########################################################################

savedatafile = "Seili_16S_Winter20"
datestamp = gsub(":", "_", gsub("[[:space:]]+", "-", date()), fixed = TRUE)
save.image(paste0(savedatafile, "-", datestamp, ".RData"))

## CCA ##################################################################################

# Load saved data?
# load("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/RData/16S_Winter20_UntilCCA/Seili_16S_Winter20-Thu-Jan-16-13_55_08-2020.RData")

## First CCA ("sediment" parameters)

# version of "sediment" CCA with HS data included and grain size excluded

# all CCA metadata files: updated for porosity-corrected data (added "_porrcorr" in file name)

# import metadata - values are the same as for BIOM construction

Env.CCA <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_ccatrim_porcorr.txt", header = TRUE)
Env.vectors <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_vectorstrim_porcorr.txt", header = TRUE)

# old commands:
# Env.CCA <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris/micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_ccatrim_porcorr.txt", header = TRUE)
# Env.vectors <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris//micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_vectorstrim_porcorr.txt", header = TRUE)
# separate "vectors" file without site + replicate data - want to plot only environmental variables as vectors!

# extract abundance matrix from the phyloseq object
Abundance = as(otu_table(rawdata.2.hellinger), "matrix")

# transpose the data and coerce to data frame
Abundance_tr <- t(Abundance)
Abundance_df <- as.data.frame(Abundance_tr)

# create ordination
CCA.vegan <- cca(Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv + HS_Inv, data = Env.CCA)

# check variance inflation factors

vif.cca(CCA.vegan)

# OLD RESULTS
##       Farm     CN_1cm      IntO2    NH4_Inv     HS_Inv 
##   5.538346   3.537974   3.101337 107.278162  88.622290

# NEW RESULTS
##       Farm     CN_1cm      IntO2    NH4_Inv     HS_Inv 
##   5.579770   4.952682   3.068927 177.646171 149.028974 

set.seed(1)
drop1(CCA.vegan, test="perm")

#         Df    AIC      F Pr(>F)   
# <none>     25.881                 
# Farm     1 26.164 1.8975  0.005 **
# CN_1cm   1 25.694 1.4953  0.010 **
# IntO2    1 26.183 1.9136  0.005 **
# NH4_Inv  1 25.552 1.3744  0.075 . 
# HS_Inv   1 25.547 1.3709  0.075 . 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# very high VIFs for NH4_Inv and HS_Inv, these are also p > 0.05 in drop1.
# based on the output, remove HS_Inv and re-run the CCA.

# version of "sediment" CCA where HS excluded

# import metadata v2 - values are the same as for BIOM construction

Env.CCA2 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_ccatrim2_porcorr.txt", header = TRUE)
Env.vectors2 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_vectorstrim2_porcorr.txt", header = TRUE)

# old commands
# Env.CCA2 <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris/micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_ccatrim2_porcorr.txt", header = TRUE)
# Env.vectors2 <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris/micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_vectorstrim2_porcorr.txt", header = TRUE)

# create ordination v2
CCA.vegan2 <- cca(Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv, data = Env.CCA2)

# check variance inflation factors for v2
vif.cca(CCA.vegan2) # now everything looks OK!

# OLD RESULTS:
##     Farm   CN_1cm    IntO2  NH4_Inv 
## 2.168320 1.726793 3.063382 1.481779

# NEW RESULTS:
##     Farm   CN_1cm    IntO2  NH4_Inv 
## 2.170685 1.726942 3.067182 1.491056 

set.seed(1)
drop1(CCA.vegan2, test="perm")

#         Df    AIC      F Pr(>F)   
# <none>     25.547                 
# Farm     1 26.803 2.8659  0.005 **
# CN_1cm   1 25.160 1.3805  0.030 * 
# IntO2    1 25.729 1.8853  0.005 **
# NH4_Inv  1 27.481 3.5026  0.005 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(summary(CCA.vegan2)) # 31.3% constrained (in older version was 38%)

with(Env.CCA2, levels(Site))
## [1] "S1" "S2" "S3" "S4" "S5" "S6"

colvec <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD")

# Cairo(file = "Seili_CCA_Trim2.png", type = "png", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")
# set.seed(1337)

par(mar=c(5,5,2,2))
plot(CCA.vegan2, display = "sites", type= "n", xlim = c(-1.25, 0.75), scaling = "sites", cex.lab = 1.5, cex.axis = 1.5)
with(Env.CCA2, points(CCA.vegan2, display = "sites", col = colvec[Site],
                      scaling = "sites", pch = 21, bg = colvec[Site], cex = 1.5))
with(Env.CCA2, legend("topleft", legend = levels(Site), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

set.seed(1)
fit2 <- envfit(CCA.vegan2, Env.vectors2, permutations = 999) # fit vectors
fit2 # print

# ***VECTORS

#             CCA1     CCA2     r2 Pr(>r)    
# Farm     0.68588 -0.72772 0.9182  0.001 ***
# CN_1cm  -0.73139  0.68196 0.1534  0.120    
# NH4_Inv -0.98013 -0.19835 0.9251  0.001 ***
# IntO2   -0.87417  0.48562 0.7273  0.001 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  Permutation: free
#  Number of permutations: 999

plot(fit2, p.max = 0.05) # limit plotting to most significant variables with argument p.max

# dev.off()

# version of "sediment" CCA with HS excluded and grain size added (cca3)

# import metadata v3 - values are the same as for BIOM construction

Env.CCA3 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_ccatrim3_porcorr.txt", header = TRUE)
Env.vectors3 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_vectorstrim3_porcorr.txt", header = TRUE)

# old commands
# Env.CCA3 <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris/micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_ccatrim3_porcorr.txt", header = TRUE)
# Env.vectors3 <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris/micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_vectorstrim3_porcorr.txt", header = TRUE)

# Env.CCA3.scaled <- as.data.frame(scale(Env.CCA3[2:6])) # scale to mean zero and unit variance
# Env.vectors3.scaled <- as.data.frame(scale(Env.vectors3))
  
# create ordination v3
CCA.vegan3 <- cca(Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain, data = Env.CCA3)

anova(CCA.vegan3, permutations = 999)

# Model: cca(formula = Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain, data = Env.CCA3)
#           Df ChiSquare      F Pr(>F)    
# Model      5   0.85662 2.5887  0.001 ***
#  Residual 24   1.58836                  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# CCA.vegan3.scaled <- cca(Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain, data = Env.CCA3.scaled)

# check variance inflation factors for v3
vif.cca(CCA.vegan3) # looks reasonable - max is 6.04

# OLD RESULTS:
##     Farm   CN_1cm    IntO2  NH4_Inv    Grain 
## 5.356383 2.305965 6.050600 2.077094 5.494467

# NEW RESULTS:
##    Farm   CN_1cm    IntO2  NH4_Inv    Grain 
## 5.350545 2.305570 6.036453 2.083147 5.476191 

set.seed(1)
drop1(CCA.vegan3, test="perm") # although this suggests that grain size doesn't "add much" ...

# NEW RESULTS

#         Df    AIC      F Pr(>F)   
# <none>     25.881                 
# Farm     1 26.193 1.9227  0.010 **
# CN_1cm   1 25.451 1.2897  0.050 * 
# IntO2    1 25.855 1.6327  0.020 * 
# NH4_Inv  1 27.259 2.8607  0.005 **
# Grain    1 25.547 1.3709  0.070 . 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# OLD RESULTS

#         Df    AIC      F Pr(>F)   
# <none>     25.881                 
# Farm     1 26.193 1.9227  0.005 **
# CN_1cm   1 25.451 1.2897  0.045 * 
# IntO2    1 25.855 1.6327  0.010 **
# NH4_Inv  1 27.259 2.8607  0.005 **
# Grain    1 25.547 1.3709  0.075 . 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(summary(CCA.vegan3)) # 35% constrained (in older data set, 41.5%)

# with(Env.CCA3, levels(Site))
## [1] "S1" "S2" "S3" "S4" "S5" "S6"

# Cairo(file = "Fig6_Bac_Winter20.png", type = "png", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")

set.seed(1337)

par(mar=c(5,5,2,2))
plot(CCA.vegan3, display = "sites", type= "n", xlim = c(-1.25, 0.75), scaling = "sites", cex.lab = 1.5, cex.axis = 1.5)
with(Env.CCA3, points(CCA.vegan3, display = "sites", col = colvec[Site],
                      scaling = "sites", pch = 21, bg = colvec[Site], cex = 2))
with(Env.CCA3, legend("topleft", legend = levels(Site), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

fit3 <- envfit(CCA.vegan3, Env.vectors3, permutations = 999, display = 'lc') # fit vectors
fit3 # print

# ***VECTORS

#  CCA1     CCA2     r2 Pr(>r)    
#  Farm     0.68417 -0.72932 0.9245  0.001 ***
#  CN_1cm  -0.73544  0.67759 0.1539  0.096 .  
#  NH4_Inv -0.98055 -0.19629 0.9289  0.001 ***
#  IntO2   -0.87405  0.48583 0.7316  0.001 ***
#  Grain    0.90902  0.41674 0.2619  0.021 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999

plot(fit3, cex = 1.3)
dev.off()
dev.off()

# fit3.scaled <- envfit(CCA.vegan3.scaled, Env.vectors3.scaled, permutations = 999) # fit vectors
# plot(fit3.scaled, p.max = 0.05)
# scaling the data produces identical results, so not using it

# dev.off()

test_margin_CCAVegan3 <- anova(CCA.vegan3, by = 'margin', parallel = 4)
test_margin_CCAVegan3$`Pr(>F)` <- p.adjust (test_margin_CCAVegan3$`Pr(>F)`, method = 'BH')
test_margin_CCAVegan3

# Permutation test for cca under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999

# NEW RESULTS (Jan 20)

# Model: cca(formula = Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain, data = Env.CCA3)
#          Df ChiSquare      F  Pr(>F)   
# Farm      1   0.12724 1.9227 0.00750 **
# CN_1cm    1   0.08536 1.2897 0.08000 . 
# IntO2     1   0.10806 1.6327 0.01833 * 
# NH4_Inv   1   0.18932 2.8607 0.00500 **
# Grain     1   0.09073 1.3709 0.04250 * 
# Residual 24   1.58836                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# OLD RESULTS

# Model: cca(formula = Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain, data = Env.CCA3)
#           Df ChiSquare      F  Pr(>F)   
# Farm      1   0.12724 1.9227 0.00500 **
# CN_1cm    1   0.08536 1.2897 0.07000 . 
# IntO2     1   0.10806 1.6327 0.01167 * 
# NH4_Inv   1   0.18932 2.8607 0.00500 **
# Grain     1   0.09073 1.3709 0.05750 . 
# Residual 24   1.58836     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

test_each_CCAvegan3 <- t(apply (Env.vectors3, 2, FUN = function (x) as.matrix (anova (cca (Abundance_df ~ x))[1:4][1,])))
test_each_CCAvegan3 <- as.data.frame (test_each_CCAvegan3)
names (test_each_CCAvegan3) <- c("Df", "Variance", "F", "Pr(>F)")
test_each_CCAvegan3.adj <- test_each_CCAvegan3
test_each_CCAvegan3.adj$`Pr(>F)` <- p.adjust (test_each_CCAvegan3$`Pr(>F)`, method = 'BH')
test_each_CCAvegan3.adj[order (test_each_CCAvegan3.adj$Variance, decreasing = TRUE),]

# NEW RESULTS (JAN 20)

#         Df  Variance        F      Pr(>F)
# NH4_Inv  1 0.2883275 3.743387 0.001666667
# Farm     1 0.2528295 3.229358 0.001666667
# IntO2    1 0.2494075 3.180684 0.001666667
# Grain    1 0.1833972 2.270591 0.002500000
# CN_1cm   1 0.1427747 1.736465 0.006000000

# OLD RESULTS

#         Df  Variance        F  Pr(>F)
# NH4_Inv  1 0.2883275 3.743387 0.00125
# Farm     1 0.2528295 3.229358 0.00125
# IntO2    1 0.2494075 3.180684 0.00125
# Grain    1 0.1833972 2.270591 0.00125
# CN_1cm   1 0.1427747 1.736465 0.00300

######################################

# Second CCA ("monitoring" parameters)

# version of "monitoring" CCA with C, P, NOx, O2 and Temp
# NOTE: these not updated because no porosity correction required for BW values.

Env.CCA.mon1 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_cca_mon1.txt", header = TRUE)
Env.vectors.mon1 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_cca_mon1.txt", header = TRUE)

# old commands
# Env.CCA.mon1 <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris/micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_cca_mon1.txt", header = TRUE)
# Env.vectors.mon1 <- read.delim("/media/jesse/My Passport/Jobs/Helsinki Postdoc/jesse/jpharris/micca/Seili/16S/denovo_greedy_otus_L50M15_noprim/sampledata_vectors_mon1.txt", header = TRUE)

CCA.mon1 <- cca(Abundance_df ~ C_1cm + O2_BW + Temp_BW + NOX_BW + P_BW, data = Env.CCA.mon1)

# check variance inflation factors

vif.cca(CCA.mon1)
##     C_1cm     O2_BW   Temp_BW    NOX_BW      P_BW 
##  8.529165 20.266952 10.993087 16.020201 18.573223

set.seed(1)
drop1(CCA.mon1, test="perm")

#         Df    AIC      F Pr(>F)   
# <none>     25.881                 
# C_1cm    1 25.739 1.5331  0.015 * 
# O2_BW    1 25.401 1.2477  0.090 . 
# Temp_BW  1 27.096 2.7151  0.005 **
# NOX_BW   1 25.404 1.2495  0.125   
# P_BW     1 26.152 1.8873  0.005 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# get some high VIFs again (>20) - suspected collinearity between O2 and Temp to begin with.
# drop 1 analysis suggests that removing O2 would have less of an effect on the AIC that removing temp would.

# regardless... "sediment" cca also has O2 and there is a good reason to keep this parameter (since we know that organic
# loading has resulted in low O2 near the farm.) so let's take out temperature and keep O2.

# version of "monitoring" CCA with C, P, NOx and O2 (TEMP REMOVED)

Env.CCA.mon2 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_cca_mon2.txt", header = TRUE)
Env.vectors.mon2 <- read.delim("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/CCAfiles/16S_sampledata_vectors_mon2.txt", header = TRUE)

CCA.mon2 <- cca(Abundance_df ~ C_1cm + O2_BW + NOX_BW + P_BW, data = Env.CCA.mon2)

# check variance inflation factors
vif.cca(CCA.mon2) # higher than in the "sediment" CCA, but still seem acceptable (highest is 13.7 - people often use 20 as an upper limit)

##     C_1cm     O2_BW    NOX_BW      P_BW 
##  7.334638 13.709284 11.538491  8.494925

set.seed(1)
drop1(CCA.mon2, test="perm")

# <none>    27.096                 
# C_1cm   1 27.432 2.0243  0.005 **
# O2_BW   1 26.986 1.6253  0.010 **
# NOX_BW  1 27.158 1.7782  0.005 **
# P_BW    1 27.115 1.7399  0.005 **

head(summary(CCA.mon2)) # 27.7% constrained (in old data set, 32.95%)

with(Env.CCA.mon2, levels(Site))
## [1] "S1" "S2" "S3" "S4" "S5" "S6"

# Cairo(file = "Fig6_BacMon_Winter20.png", type = "png", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")

set.seed(1337)

par(mar=c(5,5,2,2))
plot(CCA.mon2 , display = "sites", type= "n", xlim = c(-1.25, 1.25), scaling = "sites", cex.lab = 1.5, cex.axis = 1.5)
with(Env.CCA.mon2, points(CCA.mon2 , display = "sites", col = colvec[Site],
                          scaling = "sites", pch = 21, bg = colvec[Site], cex = 2))
with(Env.CCA.mon2, legend("topleft", legend = levels(Site), bty = "n",
                          col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

fit.mon2 <- envfit(CCA.mon2, Env.vectors.mon2, permutations = 999, display = 'lc') # fit vectors
fit.mon2 # print

# ***VECTORS

# CCA1     CCA2     r2 Pr(>r)    
# O2_BW   0.99451  0.10463 0.9523  0.001 ***
# C_1cm  -0.86497  0.50182 0.9298  0.001 ***
# NOX_BW -0.81380 -0.58114 0.7437  0.001 ***
# P_BW   -0.93321 -0.35933 0.9687  0.001 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999

plot(fit.mon2, cex = 1.3)
dev.off()
dev.off()

# dev.off()

test_margin_CCAMon2 <- anova(CCA.mon2, by = 'margin', parallel = 4)
test_margin_CCAMon2$`Pr(>F)` <- p.adjust (test_margin_CCAMon2$`Pr(>F)`, method = 'BH')
test_margin_CCAMon2

# Permutation test for cca under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999

# NEW RESULTS (JAN 20)

#  Model: cca(formula = Abundance_df ~ C_1cm + O2_BW + NOX_BW + P_BW, data = Env.CCA.mon2)
#           Df ChiSquare      F   Pr(>F)   
#  C_1cm     1   0.14316 2.0243 0.004000 **
#  O2_BW     1   0.11495 1.6253 0.015000 * 
#  NOX_BW    1   0.12576 1.7782 0.004000 **
#  P_BW      1   0.12305 1.7399 0.005333 **
#  Residual 25   1.76805                   
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# OLD RESULTS

# Model: cca(formula = Abundance_df ~ C_1cm + O2_BW + NOX_BW + P_BW, data = Env.CCA.mon2)
# Df ChiSquare      F Pr(>F)   
# C_1cm     1   0.14316 2.0243  0.004 **
# O2_BW     1   0.11495 1.6253  0.010 **
# NOX_BW    1   0.12576 1.7782  0.004 **
# P_BW      1   0.12305 1.7399  0.004 **
# Residual 25   1.76805 

test_each_CCAmon2 <- t(apply (Env.vectors.mon2, 2, FUN = function (x) as.matrix (anova (cca (Abundance_df ~ x))[1:4][1,])))
test_each_CCAmon2 <- as.data.frame (test_each_CCAmon2)
names (test_each_CCAmon2) <- c("Df", "Variance", "F", "Pr(>F)")
test_each_CCAmon2.adj <- test_each_CCAmon2
test_each_CCAmon2.adj$`Pr(>F)` <- p.adjust (test_each_CCAmon2$`Pr(>F)`, method = 'BH')
test_each_CCAmon2.adj[order (test_each_CCAmon2.adj$Variance, decreasing = TRUE),]

#        Df  Variance        F Pr(>F)
# O2_BW   1 0.2864603 3.715928  0.001
# P_BW    1 0.2750509 3.549166  0.001
# C_1cm   1 0.2531561 3.234012  0.001
# NOX_BW  1 0.2172179 2.730144  0.001

#### BACTEROIDETES ##############################

# taxon bar plots suggest that members of the phylum Bacteroidetes might be of interest here.
# let's see what things look like after subsetting

# subset bacteroidetes

bacteroidetes <- subset_taxa(rawdata.2.hellinger, Rank2 == "Bacteroidetes")

bacteroidetes.psmelt <- psmelt(bacteroidetes) # melt to long format (also converts to data frame)
# View(bacteroidetes.psmelt)

# taxon bar plot for bacteroidetes

bacteroidetes.ra.r5 <- tax_glom(bacteroidetes, taxrank = rank_names(bacteroidetes)[5]) # genus level
bacteroidetes.ra.r5 <- psmelt(bacteroidetes.ra.r5) # melt to long format (also converts to data frame)

bacteroidetes.ra.r5.005 <- subset(bacteroidetes.ra.r5, Abundance > 0.005) # cut out OTUs less than 0.5%
sum(bacteroidetes.ra.r5.005$Abundance)/sum(bacteroidetes.ra.r5$Abundance) # retains 84.1% of original data

# Cairo(file = "Seili_16S_BacteroidetesBarR5NEW.tiff", type = "tiff", units = "cm", width = 23, height = 20, pointsize = 14, dpi = 300, bg= "white")
bacteroidetes.r5.plot <- ggplot(bacteroidetes.ra.r5.005, aes(x = Replicate, y = Abundance, fill = Rank5)) + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >0.5%)") +
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))
bacteroidetes.r5.plot
# dev.off()

# bacteroidetes is about 2x more abundant in S1 than other sites; also some limited enrichment in S2+3
# the genus Flavobacteriaceae is the most abundant Bacteroidetes genus in this data set

bacteroidetes.ra.r6 <- tax_glom(bacteroidetes, taxrank = rank_names(bacteroidetes)[6]) # spp level
bacteroidetes.ra.r6 <- psmelt(bacteroidetes.ra.r6) # melt to long format (also converts to data frame)

bacteroidetes.ra.r6.005 <- subset(bacteroidetes.ra.r6, Abundance > 0.005) # cut out OTUs less than 0.5%
sum(bacteroidetes.ra.r6.005$Abundance)/sum(bacteroidetes.ra.r6$Abundance) # retains 39.7% of original data

# Cairo(file = "Seili_16S_BacteroidetesBarR6NEW.tiff", type = "tiff", units = "cm", width = 23, height = 20, pointsize = 14, dpi = 300, bg= "white")
bacteroidetes.r6.plot <- ggplot(bacteroidetes.ra.r6.005, aes(x = Replicate, y = Abundance, fill = Rank6)) + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >0.5%)") +
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))
bacteroidetes.r6.plot
# dev.off()

# at least at level above >0.05% RA, see some unique spp. in S1 too

## BAR PLOT WITH 5 MOST ABUNDANT PHYLA MINUS PROTEOBACTERIA

# first remove Proteobacteria from the data set
minusproteo <- subset_taxa(rawdata.2.hellinger, Rank2 != "Proteobacteria")

# subset further to top 5 phyla
top5ph <- sort(tapply(taxa_sums(minusproteo), tax_table(minusproteo)[, "Rank2"], sum), 
               decreasing = TRUE)[1:5]
fivephyla <- subset_taxa(minusproteo, Rank2 %in% names(top5ph))

# agglomerate data to phylum level and order the phyla from highest to lowest abundance

fivephyla.glom <- tax_glom(fivephyla, taxrank = rank_names(fivephyla)[2], NArm = T,
                           bad_empty = c(NA,""," ","\t"))

fivephyla.glom.psmelt <- psmelt(fivephyla.glom)

fivephyla.glom.psmelt$Rank2 = factor(fivephyla.glom.psmelt$Rank2, 
                                           levels = c("Bacteroidetes","Planctomycetes","Acidobacteria","Actinobacteria", "Chloroflexi"))

# View(fivephyla.glom.psmelt)

# plot the data

# Cairo(file = "Fig5_Bac_Feb19.png", type = "png", units = "cm", width = 20, height = 18, pointsize = 14, dpi = 300, bg= "white")

scaleFUN <- function(x) sprintf("%.2f", x) # for rounding y axis to two decimals!

fivephyla.plot <- ggplot(fivephyla.glom.psmelt, aes(x = Replicate, y = Abundance, fill = Rank2)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(Rank2 ~ Site, scales = "free_y") +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (%)", labels = scaleFUN) + # scaleFUN is used here
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = -90, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.position="none")
fivephyla.plot

# dev.off()

## BAR PLOT FOCUSING ON PROTEOBACTERIA

# subset to proteobacteria
proteo <- subset_taxa(rawdata.2.hellinger, Rank2 == "Proteobacteria")

# agglomerate at rank3 (class)
proteo.glom <- tax_glom(proteo, taxrank = rank_names(proteo)[3], NArm = T,
                           bad_empty = c(NA,""," ","\t"))

proteo.glom.psmelt <- psmelt(proteo.glom)

proteo.plot <- ggplot(proteo.glom.psmelt, aes(x = Replicate, y = Abundance, fill = Rank3)) + 
  geom_bar(stat = "identity", position = "stack") +
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
proteo.plot

# this shows that there's no real difference going on... so better to focus on other groups

## SAVE THE DATA ########################################################################

# savedatafile = "Seili_16S_Tku"
# datestamp = gsub(":", "_", gsub("[[:space:]]+", "-", date()), fixed = TRUE)
# save.image(paste0(savedatafile, "-", datestamp, ".RData"))

## Statistical comparisons of alpha diversity data ## (Feb 20)

# NON-NORMALISED DATA

# assign the div values to an object
div_vals <- estimate_richness(rawdata, measures=c("Observed", "Chao1", "Shannon"))

# add site data
sitenames <- c("S1", "S2", "S3", "S4", "S5", "S6")
sitenames <- rep(sitenames, each = 5)
div_vals$site <- sitenames

# need to check out library sizes? 
# sample_sums(rawdata)

# compute the analysis of variance
div.aov <- aov(Shannon ~ site, data = div_vals)

# diagnostics
library(ggfortify)
autoplot(div.aov)

# summary of the analysis
summary(div.aov)

#              Df Sum Sq Mean Sq F value   Pr(>F)    
#  site         5 2.9245  0.5849   75.51 6.93e-14 ***
#  Residuals   24 0.1859  0.0077                     
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


div.aov.tukey <- TukeyHSD(div.aov)
multcompLetters(extract_p(div.aov.tukey$site))

#   S2   S3   S4   S5   S6   S1 
# "ab"  "a"  "a"  "b"  "c"  "d" 

# Cairo(file = "Bac_ShannonTukey_Winter20.png", type = "png", units = "cm", width = 15, height = 15, pointsize = 12, dpi = 300, bg= "white")
boxplot(Shannon ~ site, data = div_vals, ylim=c(6.1, 7.4), ylab = "Shannon index", xlab = "Site")
text(x= 1, y= 6.8, labels= "d")
text(x= 2, y= 7.3, labels= "ab")
text(x= 3, y= 7.3, labels= "a")
text(x= 4, y= 7.3, labels= "a")
text(x= 5, y= 7, labels= "b")
text(x= 6, y= 6.43, labels= "c")
# dev.off()
# dev.off()

# testing effect of library size

# create library size df with site names
librarysize <- data.frame(as.vector(sample_sums(rawdata)))
librarysize$site <- sitenames
library(dplyr)
librarysize <- librarysize %>% rename(size = as.vector.sample_sums.rawdata..)

# compute the analysis of variance
size.aov <- aov(Shannon ~ site, data = librarysize)




