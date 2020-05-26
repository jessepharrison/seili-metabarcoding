## Seili Data - Phyloseq Analysis (18S v4) ######################
## Jesse Harrison (June 2018) ####################################

## CLEAN-UP, LIBRARIES, PLOT THEME ###################

# clean up R
rm(list=ls())

# libraries
library(phyloseq)
library(ggplot2)
library(vegan)
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
library(indicspecies)
library(data.table)

# other attached packages:
# [1] indicspecies_1.7.6   multcompView_0.1-7   Cairo_1.5-10         RColorBrewer_1.1-2   knitr_1.25           microbiome_1.6.0     RVAideMemoire_0.9-73
# [8] cowplot_1.0.0        plyr_1.8.4           data.table_1.12.6    gridExtra_2.3        vegan_2.5-6          lattice_0.20-38      permute_0.9-5       
# [15] ggplot2_3.2.1        phyloseq_1.28.0  

# bw theme for plots
theme_set(theme_classic())

## WORKING DIR + DATA IMPORT ##########################################

setwd("/home/jharriso/Desktop/ISME/ISME winter 2019-2020/RData/")
load("Seili18s.RData")

# set working directory (based on location of micca output)
# setwd("/home/jesse/jpharris/micca/Seili/18S/denovo_greedy_otus_L250M85_trim_noprim") # primers removed in micca

# use import_biom() to import BIOM file, phylogenetic tree file and reference sequence file
# rawdata <- import_biom("tables_sampledata_export.biom", treefilename="tree_rooted.tree", refseqfilename="otus.fasta")
# the 'export' version of the BIOM contains minimal metadata - if want more, remove '_export' from file name

## ALPHA DIVERSITY ####################################################

# Cairo(file = "Seili_18S_alpha_rawdataNEW.tiff", type = "tiff", units = "cm", width = 22, height = 12, pointsize = 14, dpi = 300, bg= "white")
alphaplot <- plot_richness(rawdata, x="Site", measures=c("Observed", "Chao1", "Shannon"))
alphaplot <- alphaplot + geom_boxplot(data = alphaplot$data, aes(x=Site, color=NULL), alpha=0.1)
alphaplot
# dev.off()

## RAREFACTION WITHOUT EQUAL SAMPLING DEPTH ################################################################

# Cairo(file = "Seili_18S_rare_rawdataNEW.tiff", type = "tiff", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")
rarecurve(t(otu_table(rawdata)), step = 100,
          cex.lab = 1.5, cex.axis = 1.5, label = FALSE, ylab = "Phylotype")
# dev.off()

# trimmed data set looks better than original

## INSPECTING THE DATA AND FILTERING OUT ARTIFACTS / LOW-PREVALENCE PHYLA #################################

table(tax_table(rawdata)[, "Rank2"], exclude = NULL)

# Ambiguous_taxa      Amoebozoa Archaeplastida   Centrohelida  Cryptophyceae       Excavata     Haptophyta 
#              1              4             32              2              5              3              4 
# Incertae Sedis   Opisthokonta            SAR           <NA> 
#             22            154            455            225 

# one phylum that occur only once: "Ambiguous_taxa"
# there are also 225 NAs that are likely to be artifacts and should be removed.

rawdata.0 <- subset_taxa(rawdata, !is.na(Rank2) & !Rank2 %in% c("", "uncharacterized"))
table(tax_table(rawdata.0)[, "Rank2"], exclude = NULL) # exclude NAs

# next, we can inspect the data in more detail to see if certain phyla are comprised mostly of low-prevalence
# features (prevalence = no. of times an OTU is observed at least once).

# compute prevalence of each feature, store as data.frame
prevdf.raw = apply(X = otu_table(rawdata.0),
                   MARGIN = ifelse(taxa_are_rows(rawdata.0), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

# add taxonomy and total read counts to this data.frame
prevdf.raw = data.frame(Prevalence = prevdf.raw,
                        TotalAbundance = taxa_sums(rawdata.0),
                        tax_table(rawdata.0))

plyr::ddply(prevdf.raw, "Rank2", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# "Ambiguous_taxa" has low prevalence (2), so it should be safe to filte.

# define phyla to filter and remove them from the data set
filterPhyla.raw = c("Ambiguous_taxa")
rawdata.1 = subset_taxa(rawdata.0, !Rank2 %in% filterPhyla.raw)
rawdata.1

## ADDITIONAL WAYS TO INSPECT PREVALENCE / OTU COUNTS ##############################################

# source script with fast_melt() function
# ("melts" OTU table into format with three main columns: SampleID, TaxaID and count)
# taxa_summary.R can be downloaded from: http://evomics.org/phyloseq/taxa_summary-r/

source("/home/jesse/jpharris/micca/taxa_summary.R", local = TRUE)

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
prevdt.raw[(Prevalence <= 0), .N]
## [1] 0

# how many singletons?
prevdt.raw[(Prevalence <= 1), .N]
## [1] 133

# how many doubletons?
prevdt.raw[(Prevalence <= 2), .N]
## [1] 234

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

# Cairo(file = "Seili_18S_PrevFilterNEW.tiff", type = "tiff", units = "cm", width = 28, height = 20, pointsize = 14, dpi = 300, bg= "white")
prevdf.raw.phyla = subset(prevdf.raw, Rank2 %in% get_taxa_unique(rawdata.1, "Rank2"))
ggplot(prevdf.raw.phyla, aes(TotalAbundance, Prevalence / nsamples(rawdata.0),color=Rank2)) +
  # Include a guess for parameter (5 % of total samples)
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Rank2) + theme(legend.position="none")
# dev.off()

# some of the low-prevalence OTUs are also somewhat high in abundance (at least for Opisthokonta and SAR).
# try singleton/doubleton removal...

## EXAMINING DATA WITHOUT SINGLETONS / DOUBLETONS

rawdata.f <- filter_taxa(rawdata, function (x) {sum(x > 0) > 2}, prune = TRUE)
rarecurve(t(otu_table(rawdata.f)), step = 100)

table(tax_table(rawdata.f)[, "Rank2"], exclude = NULL)
# single-occurrence phyla: Amoebozoa and Haptophyta. NAs: 88

rawdata.f.0 <- subset_taxa(rawdata.f, !is.na(Rank2) & !Rank2 %in% c("", "uncharacterized"))
table(tax_table(rawdata.f.0)[, "Rank2"], exclude = NULL)

prevdf.raw.f = apply(X = otu_table(rawdata.f.0),
                   MARGIN = ifelse(taxa_are_rows(rawdata.f.0), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

prevdf.raw.f = data.frame(Prevalence = prevdf.raw.f,
                        TotalAbundance = taxa_sums(rawdata.f.0),
                        tax_table(rawdata.f.0))

plyr::ddply(prevdf.raw.f, "Rank2", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Amoebozoa: 3 and Haptophyta: 7, maybe leave as is. --> use rawdata.f.0

mdt.raw.f = fast_melt(rawdata.f.0)
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

prevdf.raw.phyla.f = subset(prevdf.raw.f, Rank2 %in% get_taxa_unique(rawdata.f.0, "Rank2"))
ggplot(prevdf.raw.phyla.f, aes(TotalAbundance, Prevalence / nsamples(rawdata.f.0),color=Rank2)) +
  # Include a guess for parameter (5 % of total samples)
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Rank2) + theme(legend.position="none")

# Cairo(file = "Seili_18S_rare_rawdata_nosingleNEW.tiff", type = "tiff", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")
rarecurve(t(otu_table(rawdata.f.0)), step = 100)
# dev.off()

# this is actually more stringent than the 5% prevalence filter.
# since used prevalence filter for 16S, process the 18S data in the same way.

prevalenceThreshold.raw = 0.05 * nsamples(rawdata.1)

# execute prevalence filter, using `prune_taxa()` function
keepTaxa.raw = rownames(prevdf.raw.phyla)[(prevdf.raw.phyla$Prevalence >= prevalenceThreshold.raw)]
rawdata.2 = prune_taxa(keepTaxa.raw, rawdata.1)

# Cairo(file = "Seili_18S_rare_afterprevfiltNEW.tiff", type = "tiff", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")
rarecurve(t(otu_table(rawdata.2)), step = 100)
# dev.off()

## HELLINGER TRANSFORMATION #######################################################################################

# Hellinger transformation for filtered "uneven" data set (rawdata.1)
# (= square root of relative abundance)

# sample_sums(rawdata.2)
# ntaxa(rawdata.2)

rawdata.2.hellinger <- rawdata.2
otu_table(rawdata.2.hellinger) <-otu_table(decostand(otu_table(rawdata.2.hellinger), 
                                                     method = "hellinger"), taxa_are_rows = TRUE)

rawdata.2.hellinger <- transform_sample_counts(rawdata.2.hellinger, function(x) x/sum(x)) # convert to proportions

# export the data as a CSV file
# rawdata.2.hellinger.psmelt <- psmelt(rawdata.2.hellinger) # melt to long format (also converts to data frame)
# fwrite(rawdata.2.hellinger.psmelt, "Seili_18S_Hellinger_RA_noprim.csv")

# if want to view in table format inside R: 
# View(rawdata.2.hellinger.psmelt)

# RA without Hellinger transformation? 
# rawdata.2.nohellinger <- transform_sample_counts(rawdata.2, function(x) x/sum(x)) # convert to proportions

## BAR PLOTS ######################################################################################################

rawdata.2.hellinger.psmelt <- psmelt(rawdata.2.hellinger) # melt to long format (also converts to data frame)
# View(rawdata.2.hellinger.psmelt)

colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
             "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "#DD776655")

# rank3 plot

raw2.ra.r3 <- tax_glom(rawdata.2.hellinger, taxrank = rank_names(rawdata.2.hellinger)[3], NArm = T,
                       bad_empty = c(NA,""," ","\t")) # filters out any junk (not necessarily required but keeping it there anyway)

raw2.ra.r3 <- psmelt(raw2.ra.r3) # melt to long format (also converts to data frame)
raw2.ra.r3.02 <- subset(raw2.ra.r3, Abundance > 0.02) # cut out OTUs less than 2%
sum(raw2.ra.r3.02$Abundance)/sum(raw2.ra.r3$Abundance) # retains 97.3% of original data

# Cairo(file = "Seili_18S_r3BarNEW.tiff", type = "tiff", units = "cm", width = 23, height = 20, pointsize = 14, dpi = 300, bg= "white")
r3.plot <- ggplot(raw2.ra.r3.02 , aes(x = Replicate, y = Abundance, fill = Rank3)) + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >2%)") +
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))
r3.plot
# dev.off()

# rank4 plot

raw2.ra.r4 <- tax_glom(rawdata.2.hellinger, taxrank = rank_names(rawdata.2.hellinger)[4], NArm = T,
                       bad_empty = c(NA,""," ","\t")) # filters out any junk (not necessarily required but keeping it there anyway)

raw2.ra.r4 <- psmelt(raw2.ra.r4) # melt to long format (also converts to data frame)
raw2.ra.r4.02 <- subset(raw2.ra.r4, Abundance > 0.02) # cut out OTUs less than 2%
sum(raw2.ra.r4.02$Abundance)/sum(raw2.ra.r4$Abundance) # retains 91.7% of original data

# Cairo(file = "Seili_18S_r4BarNEW.tiff", type = "tiff", units = "cm", width = 23, height = 20, pointsize = 14, dpi = 300, bg= "white")
r4.plot <- ggplot(raw2.ra.r4.02 , aes(x = Replicate, y = Abundance, fill = Rank4)) + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(name = "Relative abundance (taxa >2%)") +
  theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 10))
r4.plot
# dev.off()

## nMDS ORDINATION ########################################################################

set.seed(1)

rawdata.2.nmds <- ordinate(physeq = rawdata.2.hellinger, method = "NMDS", distance = "bray") 
# stress = 0.148

# Cairo(file = "Seili_18S_nMDSNEW.tiff", type = "tiff", units = "cm", width = 15, height = 15, pointsize = 14, dpi = 300, bg= "white")
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

# using Hellinger transformation gives a slight improvement in terms of separating sites

## PERMANOVA FOR FACTOR SITE ############################################################

## ADONIS() (PERMANOVA global test, Bray-Curtis distance matrix)

# make a data frame from the sample_data
rawdata.2.hellinger.df <- data.frame(sample_data(rawdata.2.hellinger))

# extract bray-curtis distances

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

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Site       5    1.9823 0.39646  2.4583 0.40577  0.001 ***
# Residuals 18    2.9030 0.16128         0.59423           
# Total     23    4.8853                 1.00000           
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# post-hoc pairwise comparison

set.seed(1)
pairwise.perm.manova(rawdata.2.hellinger.bray, rawdata.2.hellinger.df$Site, nperm = 999, p.method = "BH")

# Pairwise comparisons using permutation MANOVAs on a distance matrix 

# data:  rawdata.2.hellinger.bray by rawdata.2.hellinger.df$Site
# 999 permutations 

# S1    S2    S3    S4    S5   
# S2 0.036 -     -     -     -    
# S3 0.036 0.036 -     -     -    
# S4 0.036 0.036 0.036 -     -    
# S5 0.038 0.036 0.036 0.036 -    
# S6 0.036 0.036 0.036 0.036 0.036

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
# Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)   
# Groups     5 0.016340 0.0032679 5.763    999  0.004 **
# Residuals 18 0.010207 0.0005671                       
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# as with 16S data, there is a significant multivariate dispersion effect that needs to be taken into
# account when interpreting the PERMANOVA output.

# pairwise Tukey HSD comparisons

betaHSD <- TukeyHSD(beta)

# Cairo(file = "Seili_18S_Tukey1New.tiff", type = "tiff", units = "cm", width = 15, height = 15, pointsize = 12, dpi = 300, bg= "white")
psig = as.numeric(apply(betaHSD$`group`[,2:3],1,prod)>=0)+1
op = par(mar=c(4.2,5,3.8,2))
plot(betaHSD, col = psig, yaxt = "n")
for (j in 1:length(psig)){
  axis(2,at=j,labels=rownames(betaHSD$`group`)[length(psig)-j+1],
       las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
}
par(op)
# dev.off()

# significant differences in dispersion: S4 vs. S1, S4 vs. S5, S4 vs. S6

# or can extract letters like this:
multcompLetters(extract_p(betaHSD$group))
#  S2   S3   S4   S5   S6   S1 
# "ab" "ab"  "a"  "b"  "b"  "b" 

# box plot with letters

# Cairo(file = "Seili_18S_Tukey2NEW.tiff", type = "tiff", units = "cm", width = 15, height = 15, pointsize = 12, dpi = 300, bg= "white")
boxplot(beta, xlab = "Site", ylim=c(0.27, 0.42))
text(x= 1, y= 0.359, labels= "b")
text(x= 2, y= 0.3922, labels= "ab")
text(x= 3, y= 0.394, labels= "ab")
text(x= 4, y= 0.417, labels= "a")
text(x= 5, y= 0.35, labels= "b")
text(x= 6, y= 0.371, labels= "b")
# dev.off()

## CCA ######################################################################################################

# version with "sediment" parameters based on 16S analysis - grain size included

# import metadata - values are the same as for BIOM construction
Env.CCA1 <- read.delim("/home/jesse/jpharris/micca/Seili/18S/denovo_greedy_otus_L250M85_trim_noprim/sampledata_ccatrim3_porcorr.txt", header = TRUE)
Env.vectors1 <- read.delim("/home/jesse/jpharris/micca/Seili/18S/denovo_greedy_otus_L250M85_trim_noprim/sampledata_vectorstrim3_porcorr.txt", header = TRUE)
# separate "vectors" file without site + replicate data - want to plot only environmental variables as vectors!

# extract abundance matrix from the phyloseq object
Abundance = as(otu_table(rawdata.2.hellinger), "matrix")

# transpose the data and coerce to data frame
Abundance_tr <- t(Abundance)
Abundance_df <- as.data.frame(Abundance_tr)

# create ordination
CCA.vegan1 <- cca(Abundance_df ~ Farm + CN_1cm + IntO2 + NH4_Inv + Grain, data = Env.CCA1)

# check variance inflation factors
vif.cca(CCA.vegan1) # looks reasonable - max is 6.04
##     Farm   CN_1cm    IntO2  NH4_Inv    Grain 
## 5.350545 2.305570 6.036453 2.083147 5.476191 

set.seed(1)
drop1(CCA.vegan1, test="perm")

#         Df    AIC      F Pr(>F)   
# <none>     35.617                 
# Farm     1 34.995 1.0638  0.270   
# CN_1cm   1 35.296 1.3046  0.010 **
# IntO2    1 35.231 1.2524  0.025 * 
# NH4_Inv  1 35.182 1.2130  0.035 * 
# Grain    1 35.253 1.2699  0.015 * 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

head(summary(CCA.vegan1)) # 28.9% constrained

with(Env.CCA1, levels(Site))
## [1] "S1" "S2" "S3" "S4" "S5" "S6"

colvec <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD")

set.seed(1)
fit <- envfit(CCA.vegan1, Env.vectors1, permutations = 999) # fit vectors
fit # print

## 
## ***VECTORS
## 
##             CCA1     CCA2     r2 Pr(>r)    
## Farm     0.39170  0.92009 0.6862  0.001 ***
## CN_1cm  -0.93768  0.34749 0.4549  0.003 ** 
## NH4_Inv -0.93231 -0.36167 0.4599  0.011 *  
## IntO2   -0.86963 -0.49370 0.9166  0.001 ***
## Grain    0.91185 -0.41052 0.9525  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## Permutation: free
## Number of permutations: 999

# Cairo(file = "Seili_18S_CCANEW_porcorr.tiff", type = "tiff", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")
set.seed(1337)

par(mar=c(5,5,2,2))
plot(CCA.vegan1, display = "sites", type= "n", xlim = c(-1.5, 1.5), scaling = "sites", cex.lab = 1.5, cex.axis = 1.5)
with(Env.CCA1, points(CCA.vegan1, display = "sites", col = colvec[Site],
                      scaling = "sites", pch = 21, bg = colvec[Site], cex = 1.5))
with(Env.CCA1, legend("topleft", legend = levels(Site), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

plot(fit, p.max = 0.05) # limit plotting to most significant variables with argument p.max

# dev.off()

# version with "monitoring" parameters based on 16S analysis - temperature excluded

Env.CCA.mon2 <- read.delim("/home/jesse/jpharris/micca/Seili/18S/denovo_greedy_otus_L250M85_trim_noprim/sampledata_cca_mon2.txt", header = TRUE)
Env.vectors.mon2 <- read.delim("/home/jesse/jpharris/micca/Seili/18S/denovo_greedy_otus_L250M85_trim_noprim/sampledata_vectors_mon2.txt", header = TRUE)

CCA.mon2 <- cca(Abundance_df ~ C_1cm + O2_BW + NOX_BW + P_BW, data = Env.CCA.mon2)

# check variance inflation factors
vif.cca(CCA.mon2) # higher than in the "sediment" CCA, but still seem acceptable (highest is 13.7 - people often use 20 as an upper limit)
##     C_1cm     O2_BW    NOX_BW      P_BW 
##  7.334638 13.709284 11.538491  8.494925

set.seed(1)
drop1(CCA.mon2, test="perm") # all good
##        Df    AIC      F Pr(>F)  
# <none>    35.294                 
# C_1cm   1 35.009 1.4074  0.005 **
# O2_BW   1 34.821 1.2481  0.025 * 
# NOX_BW  1 34.962 1.3672  0.010 **
# P_BW    1 34.938 1.3468  0.010 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

head(summary(CCA.mon2)) # 23.7% constrained (in old data set, was 19.3%)

with(Env.CCA.mon2, levels(Site))
## [1] "S1" "S2" "S3" "S4" "S5" "S6"

# Cairo(file = "Seili_18S_CCAMonNEW.tiff", type = "tiff", units = "cm", width = 25, height = 25, pointsize = 14, dpi = 300, bg= "white")

set.seed(1337)

par(mar=c(5,5,2,2))
plot(CCA.mon2 , display = "sites", type= "n", xlim = c(-1.25, 1.1), 
     ylim = c(-1.1, 0.75), scaling = "sites", cex.lab = 1.5, cex.axis = 1.5)
with(Env.CCA.mon2, points(CCA.mon2 , display = "sites", col = colvec[Site],
                          scaling = "sites", pch = 21, bg = colvec[Site], cex = 1.5))
with(Env.CCA.mon2, legend("topleft", legend = levels(Site), bty = "n",
                          col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

set.seed(1)
fit.mon2 <- envfit(CCA.mon2, Env.vectors.mon2, permutations = 999) # fit vectors
fit.mon2 # print

## 
## ***VECTORS
## 
## CCA1     CCA2     r2 Pr(>r)    
## O2_BW   0.88474 -0.46609 0.5659  0.001 ***
## C_1cm  -0.99999 -0.00453 0.9593  0.001 ***
## NOX_BW -0.33514  0.94217 0.1072  0.308    
## P_BW   -0.72429  0.68949 0.4926  0.002 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## Permutation: free
## Number of permutations: 999

plot(fit.mon2, p.max = 0.05) # drops NOX_BW from the list

# dev.off()

## INDICATOR SPECIES ANALYSIS (METAZOA) ##############################################

# subset metazoa

metazoa <- subset_taxa(rawdata.2.hellinger, Rank4 == "Metazoa (Animalia)")

metazoa.psmelt <- psmelt(metazoa) # melt to long format (also converts to data frame)
# View(metazoa.psmelt)

# export the data as a CSV file
# fwrite(metazoa.psmelt, "Seili_18S_Metazoa_Hellinger_RA_noprim_July18.csv")

# taxon bar plot for Metazoa (rank 6)

metazoa.ra.r6 <- tax_glom(metazoa, taxrank = rank_names(metazoa)[6])
metazoa.ra.r6 <- psmelt(metazoa.ra.r6) # melt to long format (also converts to data frame)

metazoa.ra.r6.005 <- subset(metazoa.ra.r6, Abundance > 0.005) # cut out OTUs less than 0.5%
sum(metazoa.ra.r6.005$Abundance)/sum(metazoa.ra.r6$Abundance) # retains 92.6% of original data

# Cairo(file = "Seili_18S_MetazoaBarNEW.tiff", type = "tiff", units = "cm", width = 23, height = 20, pointsize = 14, dpi = 300, bg= "white")
meta.plot <- ggplot(metazoa.ra.r6.005, aes(x = Replicate, y = Abundance, fill = Rank6)) + geom_bar(stat = "identity", position = "stack") +
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
meta.plot
# dev.off()

## indicator species analysis

counts <- data.frame(t(as.matrix(otu_table(metazoa))))
groups <- as.vector(as.matrix(sample_data(metazoa)[,"Site"]))

set.seed(8046)

indicatorsp_18s <- multipatt(counts, groups, func = "r.g", control = how(nperm = 999))
summary(indicatorsp_18s, alpha = 1, indvalcomp = T) # not corrected for multiple testing

indic_df_18s <- indicatorsp_18s$sign
indic_df_18s # without multiple correction adjustment (output omitted)

indic_df_18s[9] <- lapply(indic_df_18s[9], p.adjust, method="BH")
indic_df_18s # with Benjamini-Hochberg adjustment

# the above has lots of OTUs with p > 0.05, let's create a list showing only those with p < 0.05

S1 <- as.matrix(indic_df_18s[which(indic_df_18s$s.S1 == 1 & indic_df_18s$p.value < 0.05),])
S2 <- as.matrix(indic_df_18s[which(indic_df_18s$s.S2 == 1 & indic_df_18s$p.value < 0.05),])
S3 <- as.matrix(indic_df_18s[which(indic_df_18s$s.S3 == 1 & indic_df_18s$p.value < 0.05),])
S4 <- as.matrix(indic_df_18s[which(indic_df_18s$s.S4 == 1 & indic_df_18s$p.value < 0.05),])
S5 <- as.matrix(indic_df_18s[which(indic_df_18s$s.S5 == 1 & indic_df_18s$p.value < 0.05),])
S6 <- as.matrix(indic_df_18s[which(indic_df_18s$s.S6 == 1 & indic_df_18s$p.value < 0.05),])

r_values_18s <- rbind(S1, S2, S3, S4, S5, S6)
colnames(r_values_18s)[1:6] <-c("S1","S2","S3","S4","S5","S6")
r_values_18s 

# there's a split in the data set - species that are either in S1-S3 or species that are in S4-S6.

#           S1 S2 S3 S4 S5 S6 index      stat    p.value
# DENOVO33   1  0  1  0  0  0     8 0.7994000 0.01288889
# DENOVO74   1  0  1  0  0  0     8 0.9019890 0.01288889
# DENOVO5    0  1  1  0  0  0    12 0.7935250 0.01288889
# DENOVO6    0  1  1  0  0  0    12 0.8383155 0.01288889
# DENOVO5    0  1  1  0  0  0    12 0.7935250 0.01288889
# DENOVO6    0  1  1  0  0  0    12 0.8383155 0.01288889
# DENOVO33   1  0  1  0  0  0     8 0.7994000 0.01288889
# DENOVO74   1  0  1  0  0  0     8 0.9019890 0.01288889
# DENOVO80   0  0  1  0  0  1    18 0.7323482 0.01933333
# DENOVO345  0  0  1  0  0  0     3 0.8377231 0.03866667
# DENOVO59   0  0  0  1  0  0     4 0.7356557 0.03314286
# DENOVO20   0  0  0  1  0  1    20 0.7761780 0.01933333
# DENOVO9    0  0  0  1  0  0     4 0.8970575 0.01288889
# DENOVO35   0  0  0  1  0  0     4 0.9270611 0.01288889
# DENOVO79   0  0  0  1  0  0     4 0.9339088 0.01288889
# DENOVO54   0  0  0  1  1  1    41 0.8795412 0.01288889
# DENOVO78   0  0  0  1  1  1    41 0.7004398 0.03123077
# DENOVO54   0  0  0  1  1  1    41 0.8795412 0.01288889
# DENOVO78   0  0  0  1  1  1    41 0.7004398 0.03123077
# DENOVO220  0  0  0  0  1  0     5 0.8820532 0.01740000
# DENOVO20   0  0  0  1  0  1    20 0.7761780 0.01933333
# DENOVO8    0  0  0  0  0  1     6 0.8375701 0.01288889
# DENOVO54   0  0  0  1  1  1    41 0.8795412 0.01288889
# DENOVO78   0  0  0  1  1  1    41 0.7004398 0.03123077
# DENOVO80   0  0  1  0  0  1    18 0.7323482 0.01933333

## Range of correlation coefficients
range(r_values_18s[,"stat"]) # [1] 0.7004398 0.9339088

## Total number of indicator OTUS
length(unique(rownames(r_values_18s))) # 15 total indicator taxa following BH adjustment of data

## OTU identities based on SILVA and BLAST

# Taxa only in S1-S3 (while also present in more than just one site): 

# DENOVO5  Copepoda;Calanoida;Acartia tonsa
# Acartia tonsa isolate MT00543 18S ribosomal RNA gene, partial sequence
# 100%

# DENOVO6  Copepoda;Calanoida;Ambiguous_taxa
# Eurytemora affinis isolate N4-2E_18S 18S ribosomal RNA gene, partial sequence
# 100%

# DENOVO33 Chromadorea;Monhysterida;Ambiguous_taxa
# Halomonhystera sp. n. JVC-2014 18S ribosomal RNA gene, partial sequence
# Geomonhystera disjuncta partial 18S rRNA gene, from Netherlands, Westerschelde
# 98% (both matches)

# DENOVO74 Chromadorea;Monhysterida;Ambiguous_taxa
# Daptonema hirsutum partial 18S rRNA gene, isolated from Tamar estuary
# 91%

# Taxa only in S4-S6 (while also present in more than just one site):

# DENOVO20 Chromadorea;Monhysterida;Ambiguous_taxa
# Zygonemella striata isolate Z3 18S ribosomal RNA gene, partial sequence
# Daptonema sp. PFN-2007 18S ribosomal RNA gene, partial sequence
# 97% (both matches)

# DENOVO54 Chromadorea;Chromadorida
# Cyatholaimus sp. BHMM-2005 18S small subunit ribosomal RNA gene, partial sequence
# 99%

# DENOVO78 Chromadorea;Araeolaimida;Ambiguous_taxa
# Leptolaimus sp. 1283 18S small subunit ribosomal RNA gene, partial sequence
# 99%

# Taxa unique to a single site:

# DENOVO345 (Site 3) Chromadorea;Desmodorida;Ambiguous_taxa
# DENOVO9 (Site 4) Podocopa;Podocopida 
# DENOVO35 (Site 4) Enoplia;Enoplida;Ambiguous_taxa
# DENOVO59 (Site 4) Chromadorea;Monhysterida
# DENOVO79 (Site 4) Enoplia;Enoplida;Ambiguous_taxa
# DENOVO220 (Site 5) Chromadorea;Monhysterida
# DENOVO8 (Site 6) Podocopa;Podocopida

# Most are specific to S4; also, there are no OTUs specific to S1 or S2

# “Misc” taxa – present in S3 and S6

# DENOVO80 Chromadorea;Araeolaimida
