library(tidyverse)
library(reshape2)
library(qiime2R)
library(phyloseq)
library(scales)
library(decontam)
library(formattable)
library(ggpubr)
library(ggrepel)

## see decontam vignette here: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
## see decontam paper too: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2

########################################################################
## data imports to create phyloseq object (useful for decontam functions)
########################################################################

## import metadata
metadata <- read_csv(file="~/Repos/nhguano/data/metadata/allbat_meta.csv")
metadata$is.neg <- ifelse(metadata$SampleType=="ncontrol", TRUE, FALSE)
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$SampleID
sam = sample_data(metadata)

## import taxonomy
taxonomy <- read_delim(file="~/Repos/nhguano/data/tax/tmp.raw_bigDB_NBtax.tsv", delim="\t")
taxonomy <- taxonomy %>% separate(., 
         col = Taxon, 
         sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% 
  select(-kingdom_name)
taxonomy <- as.data.frame(apply(taxonomy, 2, function(y) gsub(".__", "", y)))
taxonomy <- as.data.frame(apply(taxonomy, 2, function(y) gsub("^$|^ $", NA, y)))
colnames(taxonomy)[1] <- "ASVid"

## import ASV table; save as both physeq object (for Decontam) and long-format (for custom plotting later)
qzapath="/Users/do/Repos/nhguano/data/qiime_qza/ASVtable/tmp.raw_table.qza"
features <- read_qza(qzapath)
mat.tmp <- features$data
rm(features, qzapath)
OTU=otu_table(mat.tmp, taxa_are_rows = TRUE)
## create physeq object:
ps <- phyloseq(OTU, sam)
rm(OTU, sam)
## remove any ASVs that are only in one sample only; drop any samples that no longer have data
psf <- prune_taxa(taxa_sums(ps) > 1, ps)
  ntaxa(ps)   ## 10,797 ASVs
  ntaxa(psf)  ## 10,790 ASVs
psf <- prune_samples(sample_sums(psf) > 0, psf)
  nsamples(ps)  ## 3,622 samples
  nsamples(psf) ## 3,255 samples
  length(which(taxa_sums(psf) == 2))  ## 612 doubletons too...
  length(which(taxa_sums(psf) == 3))  ## 577 tripletons
  length(which(taxa_sums(psf) == 5))  ## 417 5-tons
  length(which(taxa_sums(psf) == 10))  ## 235 10-tons (just 235 ASVs in at least 10 samples)
  length(which(taxa_sums(psf) == 20))  ## 98 20-tons (just 98 ASVs in at least 20 samples, or about 1/10 of the data!)
  length(which(taxa_sums(psf) == 50))  ## 46 20-tons (just 46 ASVs in at least 50 samples, or about 1/4 of our data)
  
rm(ps)

## create long-formatted; add in metadata
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
long_df <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp)
colnames(long_df) <- c("ASVid", "SampleID", "Reads")
long_df <- merge(long_df, metadata)
long_df <- merge(long_df, taxonomy)
tmp1 <- long_df %>% group_by(ASVid) %>%  summarise(nReads=sum(Reads)) %>% arrange(-nReads) %>% mutate(ASValias=paste0("ASV-", row.names(.))) %>% select(-nReads)
long_df <- merge(long_df, tmp1)
rm(tmp1, metadata, taxonomy)

########################################################################
## plotting read abundances by SampleType
########################################################################
## This ggplot example shows how our contaminant samples generally hove lower abundances
## But the way it's ordered sort of obscures the image
  ## We'll make a custom plot of this below
df <- as.data.frame(sample_data(psf)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(psf)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot() + 
  geom_point(data=df %>% filter(SampleType == "sample"), aes(x=Index, y=LibrarySize), color="gray50") +
  geom_point(data=df %>% filter(SampleType == "mock"), aes(x=Index, y=LibrarySize), color="#caa102") +
  geom_point(data=df %>% filter(SampleType == "ncontrol"), aes(x=Index, y=LibrarySize), color="#512698") +
  theme_bw() +
  scale_y_log10()

## This custom plot will examine both the read abundances and the nASVs
## plot; save as 'contam_eval_allSampls_Counts-ASVs_scatterplot'; export at 800x400
df_sumry <- long_df %>% group_by(SampleType, SampleID) %>% summarise(ASVs=n_distinct(ASVid), Reads=sum(Reads))
ggplot() +
  geom_point(data = df_sumry %>% filter(SampleType == "sample"), aes(y=ASVs, x=Reads), color="gray70", shape=1) +
  geom_point(data = df_sumry %>% filter(SampleType == "mock"), aes(y=ASVs, x=Reads), color="#caa102", size=2.5, shape=15) +
  geom_point(data = df_sumry %>% filter(SampleType == "ncontrol"), aes(y=ASVs, x=Reads), color="#512698", size=1.75, shape=17) +
  theme_bw() +
  scale_x_log10(labels=comma) + annotation_logticks(sides = "b", colour = "grey20") +
  labs(x="sequence counts per Sample", y="ASVs per sample")
  ## here we see that the vast majority of ncontrols have less than 1000 reads; two are strikingly higher
  ## mock communities have only a few higher than expected ASVs (~24 ASVs expected, so cross-index contam likely only minor issue)

rm(df_sumry)
########################################################################
## decontam function to measure ASV prevalence
########################################################################

## Underlying hypothesis is that prevalence of contaminants is greater in ncontrols than in true samples ..
## ..due to the absence of competing DNA in the sequencing process (paraphrased from Paper!)
## Chi-squared 2x2 statistic for presence-absence performed on each feature
## We (below) also can add in batch classifications for things like the sequencing library or DNA plate a sample
## See this post about the interpretation of the pscore (it's not a p value): https://github.com/benjjneb/decontam/issues/28

## Going to test a range of prevalence thresholds first:
## Function will create a series of datasets:
basic_contam.function <- function(value, label){
  tmp <- isContaminant(psf, method="prevalence", neg="is.neg", threshold=value)
  tmp <- tmp %>% mutate(ASVid=row.names(tmp))
  tmp %>% mutate(threshold=label)
}

cp_basic01 <- basic_contam.function(0.1, "0.1")
cp_basic02 <- basic_contam.function(0.2, "0.2")
cp_basic03 <- basic_contam.function(0.3, "0.3")
cp_basic05 <- basic_contam.function(0.5, "0.5")
cp_basic07 <- basic_contam.function(0.7, "0.7")
basic_contam.prev <- rbind(cp_basic01, cp_basic02, cp_basic03, cp_basic05, cp_basic07)
rm(cp_basic01, cp_basic02, cp_basic03, cp_basic05, cp_basic07)

## We can see that the number of contaminant-suspected (FALSE) or non-suspected (TRUE) varies with threshold:
## Note that there are singleton ASVs that apparently aren't being filtered out by Phyloseq's function... labeleing those in this table
basic_contam.prev$contaminant <- ifelse(basic_contam.prev$prev==1, "singleton", basic_contam.prev$contaminant)
basic_contam_table <- basic_contam.prev %>% 
  group_by(threshold, contaminant) %>%
  tally() %>% 
  dcast(., threshold~contaminant, value.var = "n")
basic_contam_table <- basic_contam_table[c(1,2,4,3)]

## plot; save as 'decontam_prevTable_basic'; export at 480x240
formattable(basic_contam_table)
rm(basic_contam_table)
## plot; save as 'decontam_prevHist_basic'; export at 600x300
ggplot(data = basic_contam.prev %>% filter(threshold=="0.1"), aes(p)) + 
  geom_histogram(bins=100, color='black', fill="gray50") +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15)) +
  #geom_vline(xintercept = 0.2, color="red", linetype="dashed") +
  labs(x="decontam Score", y="Number ASVs")

## Can add in "batch" variable for each sequencing run or any other set of variables
## We have to drop out the mock samples for DNAplate-based analysis to work (they don't have a plate number)
psfm = subset_samples(psf, SampleType != "mock")                ## drops mock samples

batch_contam.function <- function(data, value, label, batchtype, batch){
  tmp <- isContaminant(seqtab=data, method="prevalence", neg="is.neg", threshold=value, batch=batch)
  tmp <- tmp %>% mutate(ASVid=row.names(tmp))
  tmp %>% 
    mutate(threshold=label) %>% 
    mutate(batch=batchtype)
}

## Calculate for each batch combination of threshold and batch type
cp_seq_01 <- batch_contam.function(psf, 0.1, "0.1", "SeqBatch", "SeqBatch")
cp_seq_02 <- batch_contam.function(psf, 0.2, "0.2", "SeqBatch", "SeqBatch")
cp_seq_03 <- batch_contam.function(psf, 0.3, "0.3", "SeqBatch", "SeqBatch")
cp_seq_05 <- batch_contam.function(psf, 0.5, "0.5", "SeqBatch", "SeqBatch")
cp_seq_07 <- batch_contam.function(psf, 0.7, "0.7", "SeqBatch", "SeqBatch")
cp_dna_01 <- batch_contam.function(psfm, 0.1, "0.1", "DNAplate", "DNAplate")
cp_dna_02 <- batch_contam.function(psfm, 0.2, "0.2", "DNAplate", "DNAplate")
cp_dna_03 <- batch_contam.function(psfm, 0.3, "0.3", "DNAplate", "DNAplate")
cp_dna_05 <- batch_contam.function(psfm, 0.5, "0.5", "DNAplate", "DNAplate")
cp_dna_07 <- batch_contam.function(psfm, 0.7, "0.7", "DNAplate", "DNAplate")

rm(psfm)

batch_contam.prev <- rbind(cp_seq_01, cp_seq_02, cp_seq_03, cp_seq_05, cp_seq_07, cp_dna_01, cp_dna_02, cp_dna_03, cp_dna_05, cp_dna_07)
rm(cp_seq_01, cp_seq_02, cp_seq_03, cp_seq_05, cp_seq_07, cp_dna_01, cp_dna_02, cp_dna_03, cp_dna_05, cp_dna_07)

## create table showing the TRUE/FALSE ASVs per prevalence and batch classes
batch_contam.prev$contaminant <- ifelse(batch_contam.prev$prev==1, "singleton", batch_contam.prev$contaminant)
batch_contam.table <- batch_contam.prev %>% 
  group_by(threshold, batch, contaminant) %>% 
  tally() %>% 
  spread(., contaminant, n) %>% arrange(batch, threshold)
batch_contam.table <- batch_contam.table[c(1,2,3,5,4)]

## plot; save as 'decontam_prevTable_batch'; export at 800x400
formattable(batch_contam.table)
rm(batch_contam.table)
## plot histogram; save as 'decontam_prevHist_batch'; export at 500x500
ggplot(data = batch_contam.prev %>% filter(threshold=="0.1"), aes(p, fill=batch)) + 
  geom_histogram(bins=100, color="black") +
  facet_grid(batch ~ threshold) +
  scale_fill_manual(values=c("grey25", "gray75")) +
  scale_y_continuous(breaks=c(0, 400, 800)) +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size=14), axis.title = element_text(size=15),
        strip.background.x = element_blank(), strip.text.x = element_blank(),
        strip.text.y = element_text(size=14)) +
  labs(x="decontam Score", y="Number ASVs")

rm(basic_contam.function, basic_contam.prev)
########################################################################
## decontam ASVs - figures plotting the number of samples and number of reads per SampleType and ContamType
## making 3 plots, then stitching together
########################################################################
contam_dna01_ASVs <- batch_contam.prev %>% 
  filter(threshold=="0.1" & batch=="DNAplate" & contaminant==TRUE) %>% 
  select(ASVid)
contam_dna02_ASVs <- batch_contam.prev %>% 
  filter(threshold=="0.2" & batch=="DNAplate" & contaminant==TRUE) %>% 
  select(ASVid)
contam_seq01_ASVs <- batch_contam.prev %>% 
  filter(threshold=="0.1" & batch=="SeqBatch" & contaminant==TRUE) %>% 
  select(ASVid)
contam_seq02_ASVs <- batch_contam.prev %>% 
  filter(threshold=="0.2" & batch=="SeqBatch" & contaminant==TRUE) %>% 
  select(ASVid)
contam_basic02_ASVs <- basic_contam.prev %>% 
  filter(threshold=="0.2" & contaminant == TRUE) %>% select(ASVid)

## not sure if it's better to combine these or not?
## notrun: all02_ASVs <- c(contam_basic02_ASVs$ASVid, contam_seq03_ASVs$ASVid, contam_dna03_ASVs$ASVid) %>% unique(.)

## This plot is going to just use the ASVs identified in the "SeqBatch" dataset.. 
## More controls per SeqBatch to better model the data; too few per DNAplate (and prevalence does better with more NTCs)
## Could combine from 'all02_ASVs' instead to be more conservative
contam_df_sumry <- long_df %>% 
  filter(ASVid %in% contam_seq01_ASVs$ASVid) %>%                  ## choose appropriate input here
  group_by(SampleType, ASVid) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% 
  mutate(Status="contam")
noncontam_df_sumry <- long_df %>% 
  filter(!ASVid %in% contam_seq01_ASVs$ASVid) %>%                  ## choose appropriate input here
group_by(SampleType, ASVid) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% 
  mutate(Status="non_contam")
contam_plot_sumry <- rbind(contam_df_sumry, noncontam_df_sumry)
#contam_plot_sumry$Status <- ifelse(is.na(contam_plot_sumry$Status) & contam_plot_sumry$SampleType=="ncontrol", "ncontrol_only", contam_plot_sumry$Status)
#contam_plot_sumry$Status <- ifelse(is.na(contam_plot_sumry$Status) & contam_plot_sumry$SampleType=="sample", "sample_only", contam_plot_sumry$Status)
#contam_plot_sumry$Status <- ifelse(is.na(contam_plot_sumry$Status) & contam_plot_sumry$SampleType=="mock", "mock_only", contam_plot_sumry$Status)
contam_plot_sumry <- contam_plot_sumry %>% filter(SampleType!="contaminant")
rm(contam_df_sumry, noncontam_df_sumry)

cpc <- ggplot(data = contam_plot_sumry %>% filter(SampleType=="ncontrol"), 
       aes(x=Reads, y=Samples, color=Status)) +
  geom_point(shape=17, size=3) +
  theme_bw() +
  scale_color_manual(values=c("firebrick", "dodgerblue2")) +
  scale_x_log10(labels=comma_format(accuracy = 1)) +
  scale_y_continuous(breaks = c(0, 25, 50)) +
  annotation_logticks(sides="b") +
  theme(legend.position = "right", legend.text = element_text(size=12), legend.title = element_blank())

cpm <- ggplot(data = contam_plot_sumry %>% filter(SampleType=="mock"), 
       aes(x=Reads, y=Samples, color=Status)) +
  geom_point(shape=15, size=2, alpha=0.8) +
  theme_bw() +
  scale_color_manual(values=c("firebrick", "dodgerblue2")) +
  scale_x_log10(labels=comma_format(accuracy = 1)) +
  scale_y_continuous(breaks = c(0,5,10)) +
  annotation_logticks(sides="b") +
  theme(legend.position = "right", legend.text = element_text(size=12), legend.title = element_blank())

cps <- ggplot(data = contam_plot_sumry %>% filter(SampleType=="sample"), 
       aes(x=Reads, y=Samples, color=Status)) +
  geom_point(shape=1) +
  theme_bw() +
  scale_color_manual(values=c("firebrick", "dodgerblue2")) +
  scale_x_log10(labels=comma_format(accuracy = 1)) +
  scale_y_continuous(breaks = c(0, 400, 800)) +
  annotation_logticks(sides="b") +
  theme(legend.position = "right", legend.text = element_text(size=12), legend.title = element_blank())

## plot; save as 'decontam_Reads-ASVs_ContamOrNot; export at 
ggarrange(cpc, cpm, cps, common.legend = FALSE,
          labels = c("A", "B", "C"), nrow = 3)
rm(cpc, cpm, cps)
  ## note that the biggest abundance Contamination ASVs aren't identified as likely contaminants
  ## I wonder if that contaminant that is in all the NTCs is also in most of the true samples (so you can't really see a big diff...)

########################################################################
## we can also plot the decontam-style figure showing the prevaelence of ncontrol vs. samples
########################################################################
psf.pa <- transform_sample_counts(psf, function(abund) 1*(abund>0))
psf.pa.neg <- prune_samples(sample_data(psf.pa)$is.neg == "TRUE", psf.pa)
psf.pa.pos <- prune_samples(sample_data(psf.pa)$is.neg == "FALSE", psf.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(psf.pa.pos), pa.neg=taxa_sums(psf.pa.neg),
                    contaminant=batch_contam.prev %>% filter(batch=="DNAplate" & threshold=="0.1") %>% select(contaminant))
df.pa$ASVid <- row.names(df.pa)
ASVkey <- long_df %>% group_by(ASVid, ASValias) %>% summarise(Samples=n_distinct(SampleID)) %>% select(-Samples)
df.pa <- merge(df.pa, ASVkey)
rm(ASVkey)

## plot scatterplot; save as 'decontam_Prevalence_ncontrolANDsample'; export at 700x700 (WARNING: changing dimensions messes with labels)
ggplot() + 
  geom_point(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) +
  scale_color_manual(values=c("dodgerblue2", "green", "firebrick")) +
  geom_label_repel(data=df.pa %>% filter(contaminant=="FALSE" & pa.neg >= 9), aes(x=pa.neg, y=pa.pos, label=ASValias),
                   force = 2, nudge_x = -5, nudge_y=95, segment.colour = "gray60", size=3) +
  geom_label_repel(data=df.pa %>% filter(contaminant=="TRUE" & pa.neg >= 9) %>% filter(pa.neg < 40), aes(x=pa.neg, y=pa.pos, label=ASValias),
                   force = 2, nudge_x = 4, nudge_y=-70, segment.colour = "gray60", size=3) +
  geom_label_repel(data=df.pa %>% filter(pa.neg > 40), aes(x=pa.neg, y=pa.pos, label=ASValias),
                   nudge_x = -5, segment.colour = "gray60", size=3) +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  theme_bw() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12),
        legend.position="top", legend.text = element_text(size=12), legend.title = element_blank())
rm(psf.pa, psf.pa.neg, psf.pa.pos, df.pa)
########################################################################
## future goals:
#1a. use mock data to highlight the bleed in (similar to tidybug plot); noramlize reads to per-sample
#1b. use mock data to highlight which mock samples are likely contaminants
#1c. drop those mock samples from true samples

#2. drop any NTCs unique to the negative controls
########################################################################

## append the `long_df` so that we set "contam" or not based on whether it was in the DNAplate-batched decontam with threshold 0.1...
long_df$ContamStatus <- ifelse(long_df$ASVid %in% contam_dna01_ASVs$ASVid, "contam", "noncontam")

## create a summary of all ASVs with total sequencing depth and number of detections; group by the sample type too
contamASVtaxa <- long_df %>% 
  filter(ContamStatus=="contam") %>% 
  group_by(ASValias, phylum_name, class_name, order_name, family_name, genus_name, species_name, SampleType) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads))

## Select the most prevalent ncontrol ASVs that are considered contaminants from the previous plot
## A few of these are highly abundant, but most occur in just a few samples
long_df %>% filter(SampleType=="ncontrol") %>% group_by(SampleID) %>% summarise(ASVs=n_distinct(ASVid)) %>% nrow()  ## there were 206 ncontrols to start with?
  ## pretty amazing there are rarely ever more than 5 ncontrol samples with the same ASV... really refutes the idea of pervasive contam

## plotting the per-sample abundances from the `decontam_Prevalence_ncontrolANDsample` plot:
prev_list <- c("ASV-1", "ASV-5", "ASV-12", "ASV-20",
               "ASV-2", "ASV-3", "ASV-7", "ASV-11")
prev_plotdat <- long_df %>% filter(ASValias %in% prev_list) %>% select(ASValias, Reads, ContamStatus, SampleType)
prev_plotdat_types <- c("ncontrol", "sample")
prev_plotdat$ASValias <- factor(prev_plotdat$ASValias, levels = c("ASV-1", "ASV-5", "ASV-12", "ASV-20", 
                                                                  "ASV-2", "ASV-3", "ASV-7", "ASV-11"))

## plot facets; save as 'decontam_selectASVabundance-perSample_contamComparison'; export at 
ggplot(prev_plotdat %>% filter(SampleType %in% prev_plotdat_types), 
       aes(x=SampleType, y=Reads, color=ContamStatus, shape=SampleType)) +
  geom_jitter() +
  facet_wrap(ContamStatus ~ ASValias, nrow = 2) +
  scale_color_manual(values=c("firebrick", "dodgerblue2")) +
  scale_shape_manual(values=c(15,1)) +
  scale_y_continuous(trans = 'log2', labels = comma_format(accuracy = 1)) +
  #scale_y_log10(labels=comma_format(accuracy = 1)) +
  theme_bw() +
  theme(legend.position = "none") + 
  labs(x="", y="sequence counts")

###########
## 1. are there any ASVs among all NTCs we sequenced that are unique to NTCs?
## 2. are there any ASVs among all samples that are unique to samples?
## 3. how many ASVs are in both samples and NTCs?
## 4. how many mock ASVs are in samples? NTCs?

## 1. 
mockASVs <- long_df %>% filter(SampleType=="mock") %>% select(ASVid)
ncontrolASVs <- long_df %>% filter(SampleType=="ncontrol") %>% select(ASVid)
sampleASVs <- long_df %>% filter(SampleType=="sample") %>% select(ASVid)

mns_all <- intersect(intersect(mockASVs$ASVid, ncontrolASVs$ASVid), sampleASVs$ASVid) ## just 20 ASVs in all 3?
mn_only <- intersect(mockASVs$ASVid, ncontrolASVs$ASVid) %>% setdiff(., sampleASVs$ASVid) ## just 1 ASV exclusive to mock and ncontrols
ms_only <- intersect(mockASVs$ASVid, sampleASVs$ASVid) %>% setdiff(., ncontrolASVs$ASVid) ## just 17 ASVs exclusive to mock and samples
ns_only <- intersect(ncontrolASVs$ASVid, sampleASVs$ASVid) %>% setdiff(., mockASVs$ASVid) ## 365 ASVs between negs and samples
  ## another big of evidence suggesting that cross talk is likely between plates, not between sequencing runs
  ## 
m_only
n_only
s_only