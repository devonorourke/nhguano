library(tidyverse)
library(reshape2)
library(qiime2R)
library(phyloseq)
library(scales)
library(decontam)
library(formattable)
library(ggpubr)
library(ggrepel)
library(ape)

## see decontam vignette here: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
## see decontam paper too: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2

## Data produced from this workflow is more fully described here: https://github.com/devonorourke/nhguano/blob/master/docs/decontam_workflow.md 

########################################################################
## data imports to create phyloseq object (useful for decontam functions)
########################################################################

## import metadata
metadata <- read_csv(file="https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/allbat_meta.csv")
metadata$is.neg <- ifelse(metadata$SampleType=="ncontrol", TRUE, FALSE)
metadata <- as.data.frame(metadata)

## create physeq metadata and taxonomy data for analyses:
row.names(metadata) <- metadata$SampleID
phy_meta <- sample_data(metadata)

## import ASV table; save as both physeq object (for Decontam) and long-format (for custom plotting later)
## import the original ASV table (not clustered):
## change path below to your own path!
# download.file(url = "https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/ASVtable/tmp.raw_table.qza",
#               destfile = "~/Desktop/tmp.raw_table.qza")

qzapath="~/github/nhguano/data/qiime_qza/ASVtable/tmp.raw_table.qza" ## setting file path to my local machine
features <- read_qza(qzapath)
mat.tmp <- features$data
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
OTU <- otu_table(df.tmp, taxa_are_rows = TRUE) ## import as physeq object 
## create phyloseq object
ps <- phyloseq(OTU, phy_meta)
rm(OTU, phy_meta)
## create long_df object
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
long_df <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp)
colnames(long_df) <- c("ASVid", "SampleID", "Reads")
long_df <- merge(long_df, metadata)
tmp1 <- long_df %>% group_by(ASVid) %>%  summarise(nReads=sum(Reads)) %>% arrange(-nReads) %>% mutate(ASValias=paste0("ASV-", row.names(.))) %>% select(-nReads)
long_df <- merge(long_df, tmp1)
rm(tmp1, features)

## How rare are these ASVs anyway?
nASVs <- long_df %>% group_by(ASVid) %>% summarise(Samples=n_distinct(SampleID)) %>% arrange(Samples)
nrow(nASVs %>% filter(Samples == 1))    ## singleton ASVs: 7,641
nrow(nASVs %>% filter(Samples >= 2))    ## singleton ASVs: 3,156
nrow(nASVs %>% filter(Samples >= 5))    ## singleton ASVs: 1,185
nrow(nASVs %>% filter(Samples >= 50))    ## singleton ASVs: 127
nrow(nASVs %>% filter(Samples >= 100))    ## singleton ASVs: 39
rm(nASVs)

## How many unique samples per Study type?
long_df %>% group_by(SampleType, SampleID) %>% summarise(Reads=sum(Reads)) %>% 
  group_by(SampleType) %>% summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% 
  mutate(pSamples=round(Samples/sum(Samples),4) * 100) %>% 
  mutate(pReads=round(Reads/sum(Reads),4) * 100) %>% 
  select(SampleType, Samples, pSamples, Reads, pReads)
## for all samples (could have as little as 1 read per sample!)..
## 183 ncontrols make up about 6% of all samples, yet contribute just 0.43% of all reads
## guano samples comprise majority of samples and sequence space

## If we filter by requiring at least 500 reads per sample (and thus drop out any sample that had very few overall reads)..
## ..which would be dropped by our rarefying sampling depth regardless, how many samples/reads per SampleType now?
long_df %>% group_by(SampleType, SampleID) %>% summarise(Reads=sum(Reads)) %>% filter(Reads >= 500) %>% 
  group_by(SampleType) %>% summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% 
  mutate(pSamples=round(Samples/sum(Samples),4) * 100) %>% 
  mutate(pReads=round(Reads/sum(Reads),4) * 100) %>% 
  select(SampleType, Samples, pSamples, Reads, pReads)
## Filtering demonstrates that among those samples with at least 500 reads...
## We've lost almost all of our control samples (and about 50% of guano samples!) so a read minmum disproportionately from negatives
## There are fewer overall reads remaining to control samples now than before (0.37% vs 0.43%), even after removing the "low abundant" control samples


########################################################################
## plotting read abundances by SampleType
########################################################################

## This plot will examine both the read abundances and the nASVs
## generate the data 
df_sumry <- long_df %>% group_by(SampleType, SampleID) %>% summarise(ASVs=n_distinct(ASVid), Reads=sum(Reads))
## set levels
df_sumry$SampleType <- factor(df_sumry$SampleType, levels = c("sample", "mock", "ncontrol"))

## plot; save as 'contam_eval_allSampls_Counts-ASVs_scatterplot'; export at 800x400
ggplot(data = df_sumry %>% filter(SampleType != "contaminant") %>% filter(Reads >= 500), 
       aes(y=ASVs, x=Reads, color=SampleType, shape=SampleType, label=SampleID)) +
  geom_point(data = df_sumry %>% filter(SampleType == "sample")) +
  geom_point(data = df_sumry %>% filter(SampleType == "mock"), size=2) +
  geom_point(data = df_sumry %>% filter(SampleType == "ncontrol")) +
  scale_shape_manual(values=c(15,17, 1)) + scale_color_manual(values = c("#caa102","#512698", "gray60")) +
  geom_vline(xintercept = 500, color="red", linetype="dotted") +
  theme_bw() + theme(legend.position = "top") + guides(label=FALSE) +
  scale_x_log10(limits=c(1, 1800000), labels=comma_format(accuracy = 1)) +  annotation_logticks(sides = "b", colour = "grey20") +
  labs(x="sequence counts per Sample", y="ASVs per sample", 
       subtitle = "Per-sample read abundances and ASV detections", caption = "Dotted line at 500 reads") +
  theme(legend.position = "top") +
  guides(label=FALSE)
## here we see that the vast majority of ncontrols have less than 1000 reads; two are strikingly higher
## mock communities have only a few higher than expected ASVs (~24 ASVs expected, so cross-index contam likely only minor issue)

## Same plot as above, but we remove any sample with fewer than 500 reads...
## plot; save as 'contam_eval_allSampls_Counts-ASVs_scatterplot_filt500ReadMin'; export at 600x600
ggplot(data = df_sumry %>% filter(SampleType != "contaminant") %>% filter(Reads >= 500), 
       aes(y=ASVs, x=Reads, color=SampleType, shape=SampleType, label=SampleID)) +
  geom_label_repel(data = df_sumry %>% filter(SampleID=="negoro14D04"), segment.alpha = 0.5, segment.colour = "purple", nudge_y = 70, size=3) +
  geom_label_repel(data = df_sumry %>% filter(SampleID=="negoro37F08"), segment.alpha = 0.5, segment.colour = "purple", nudge_x = .8, nudge_y = 50, size=3) +
  geom_label_repel(data = df_sumry %>% filter(SampleID=="negoro34A01"), segment.alpha = 0.5, segment.colour = "purple", nudge_y = 85, size=3, fill="white") +
  geom_label_repel(data = df_sumry %>% filter(SampleID=="negoro33G07"), segment.alpha = 0.5, segment.colour = "purple", nudge_x = -.15, nudge_y = 40, size=3, fill="white") +
  geom_label_repel(data = df_sumry %>% filter(SampleID=="negoro14G07"), segment.alpha = 0.5, segment.colour = "purple", nudge_x = .1, nudge_y = 70, size=3) +
  geom_point(data = df_sumry %>% filter(SampleType == "sample" & Reads >= 500)) +
  geom_point(data = df_sumry %>% filter(SampleType == "mock" & Reads >= 500), size=2) +
  geom_point(data = df_sumry %>% filter(SampleType == "ncontrol" & Reads >= 500)) +
  scale_shape_manual(values=c(15,17, 1)) + scale_color_manual(values = c("#caa102","#512698", "gray60")) +
  geom_vline(xintercept = 500, color="red", linetype="dotted") +
  theme_bw() + theme(legend.position = "top") + guides(label=FALSE) +
  scale_x_log10(labels=comma_format(accuracy = 1)) +  annotation_logticks(sides = "b", colour = "grey20") +
  labs(x="sequence counts per Sample", y="ASVs per sample",
       subtitle = "Per-sample read abundances and ASV detections\nSamples must have at least 500 sequence counts",
       caption = "Dotted line at 500 reads\nNTC with most total sequence counts labeled")

rm(df_sumry)
########################################################################
## decontam function to measure ASV prevalence
########################################################################

## Underlying hypothesis is that prevalence of contaminants is greater in ncontrols than in true samples ..
## ..due to the absence of competing DNA in the sequencing process (paraphrased from Paper!)
## Chi-squared 2x2 statistic for presence-absence performed on each feature
## We (below) also can add in batch classifications for things like the sequencing library or DNA plate a sample
## See this post about the interpretation of the pscore (it's not a p value): https://github.com/benjjneb/decontam/issues/28

## Filtering phyloseq object to remove any instances where ASV in just single sample
## function from here: https://github.com/joey711/phyloseq/issues/917
psf <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)

## Going to test a range of prevalence thresholds first:
## Function will create a series of datasets:
basic_contam.function <- function(threshold, label){
  #tmp <- isContaminant(psf, method="prevalence", neg="is.neg", threshold=threshold)
  tmp <- isContaminant(psf, method="prevalence", neg="is.neg", threshold=threshold)
  tmp <- tmp %>% mutate(ASVid=row.names(tmp))
  tmp %>% mutate(threshold=label) %>% mutate(batch="basic")
}

cp_basic01 <- basic_contam.function(0.1, "0.1")
cp_basic02 <- basic_contam.function(0.2, "0.2")
cp_basic03 <- basic_contam.function(0.3, "0.3")
cp_basic05 <- basic_contam.function(0.5, "0.5")
basic_contam.prev <- rbind(cp_basic01, cp_basic02, cp_basic03, cp_basic05)
rm(cp_basic01, cp_basic02, cp_basic03, cp_basic05)

## We can see that the number of contaminant-suspected (FALSE) or non-suspected (TRUE) varies with threshold:
basic_contam_table <- basic_contam.prev %>% 
  group_by(threshold, contaminant) %>%
  tally() %>% 
  dcast(., threshold~contaminant, value.var = "n")

## plot; save as 'decontam_prevTable_basic'; export at 300x225
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

## Of the ASVs flagged by the model, how many are private to the NTCs only? 
allASVs <- unique(long_df$ASVid)
cp_NTCbasic_ASVs_1 <- basic_contam.prev %>% filter(threshold=="0.1" & contaminant==TRUE) %>% pull(ASVid)
cp_NTCbasic_ASVs_2 <- basic_contam.prev %>% filter(threshold=="0.2" & contaminant==TRUE) %>% pull(ASVid)
cp_NTCbasic_ASVs_3 <- basic_contam.prev %>% filter(threshold=="0.3" & contaminant==TRUE) %>% pull(ASVid)
cp_NTCbasic_ASVs_5 <- basic_contam.prev %>% filter(threshold=="0.5" & contaminant==TRUE) %>% pull(ASVid)
uniqNTC_ASVs <- long_df %>% filter(SampleType=="ncontrol") %>% pull(ASVid) %>% unique()
uniqSamp_ASVs <- long_df %>% filter(SampleType=="sample") %>% pull(ASVid) %>% unique()
privateNTC_ASVS <- setdiff(uniqNTC_ASVs, uniqSamp_ASVs)
private_NTCbasic_ASVs_1 <- intersect(cp_NTCbasic_ASVs_1, privateNTC_ASVS)  ## 1 NTC flagged as "contaminant" isn't in true samples; present in 2 samples, generated just 6 reads
private_NTCbasic_ASVs_2 <- intersect(cp_NTCbasic_ASVs_2, privateNTC_ASVS)  ## same one NTC as above
private_NTCbasic_ASVs_3 <- intersect(cp_NTCbasic_ASVs_3, privateNTC_ASVS)   ## same one NTC as above
private_NTCbasic_ASVs_5 <- intersect(cp_NTCbasic_ASVs_5, privateNTC_ASVS)   ## same one NTC as above

## Of al the ASVs flagged by the model, how many are private to true samples?
cp_SAMPbasic_ASVs_1 <- basic_contam.prev %>% filter(threshold=="0.1" & contaminant==FALSE) %>% pull(ASVid)
cp_SAMPbasic_ASVs_2 <- basic_contam.prev %>% filter(threshold=="0.2" & contaminant==FALSE) %>% pull(ASVid)
cp_SAMPbasic_ASVs_3 <- basic_contam.prev %>% filter(threshold=="0.3" & contaminant==FALSE) %>% pull(ASVid)
cp_SAMPbasic_ASVs_5 <- basic_contam.prev %>% filter(threshold=="0.5" & contaminant==FALSE) %>% pull(ASVid)
privateSAMP_ASVS <- setdiff(uniqSamp_ASVs, uniqNTC_ASVs)
private_SAMPbasic_ASVs_1 <- intersect(cp_SAMPbasic_ASVs_1, privateSAMP_ASVS)
length(private_SAMPbasic_ASVs_1)  ## 2,762 ASVs are not found in any negative control
private_SAMPbasic_ASVs_2 <- intersect(cp_SAMPbasic_ASVs_2, privateSAMP_ASVS)
private_SAMPbasic_ASVs_3 <- intersect(cp_SAMPbasic_ASVs_3, privateSAMP_ASVS)
private_SAMPbasic_ASVs_5 <- intersect(cp_SAMPbasic_ASVs_5, privateSAMP_ASVS)
length(private_SAMPbasic_ASVs_2)  ## 2,762 ASVs are not found in any negative control
length(private_SAMPbasic_ASVs_3)  ## 2,762 ASVs are not found in any negative control
length(private_SAMPbasic_ASVs_5)  ## 2,762 ASVs are not found in any negative control

## Finally, how many ASVs are shared between the NTCs and true samples?
common_ASVS_NTCandSAMP <- intersect(uniqNTC_ASVs, uniqSamp_ASVs)
length(common_ASVS_NTCandSAMP)
## 382 ASVs are shared among these data types (or about 12% overall)

rm(list = ls(pattern = "cp_basic_ASVs_*"))
rm(list = ls(pattern = "cp_NTCbasic_ASVs_*"))
rm(list = ls(pattern = "cp_SAMPbasic_ASVs_*"))
rm(list = ls(pattern = "private_*"))
rm(common_ASVS_NTCandSAMP, uniqNTC_ASVs, uniqSamp_ASVs, allASVs)
rm(basic_contam.function, basic_contam_table)

########################################################################
## Batch effects with Decontam
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
cp_dna_01 <- batch_contam.function(psfm, 0.1, "0.1", "DNAplate", "DNAplate")
cp_dna_02 <- batch_contam.function(psfm, 0.2, "0.2", "DNAplate", "DNAplate")
cp_dna_03 <- batch_contam.function(psfm, 0.3, "0.3", "DNAplate", "DNAplate")
cp_dna_05 <- batch_contam.function(psfm, 0.5, "0.5", "DNAplate", "DNAplate")

rm(psfm)

batch_contam.prev <- rbind(cp_seq_01, cp_seq_02, cp_seq_03, cp_seq_05, cp_dna_01, cp_dna_02, cp_dna_03, cp_dna_05)
rm(cp_seq_01, cp_seq_02, cp_seq_03, cp_seq_05, cp_dna_01, cp_dna_02, cp_dna_03, cp_dna_05)

## create table showing the TRUE/FALSE ASVs per prevalence and batch classes
batch_contam.table <- batch_contam.prev %>% 
  group_by(threshold, batch, contaminant) %>% 
  tally() %>% 
  spread(., contaminant, n) %>% arrange(batch, threshold)

## plot; save as 'decontam_prevTable_batch'; export at 350x400
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

########################################################################
########################################################################
## How many ASVs are present in the following sets?
all.prev <- rbind(basic_contam.prev, batch_contam.prev)
rm(basic_contam.prev, batch_contam.prev)


## function counts the number of shared/distinct ASVs for the 7 sets among the 3 groups ("basic", "DNAplate",  or "SeqBatch")
## function works on a per-threshold basis
commonASVcounter <- function(ThisValue){
  basiccontamASVs <- all.prev %>% filter(threshold==ThisValue & batch=="basic" & contaminant==TRUE) %>% pull(ASVid)
  batchDNAcontamASVs <- all.prev %>% filter(threshold==ThisValue & batch=="DNAplate" & contaminant==TRUE) %>% pull(ASVid)
  batchSeqcontamASVs <- all.prev %>% filter(threshold==ThisValue & batch=="SeqBatch" & contaminant==TRUE) %>% pull(ASVid)
  all_shared <- data.frame(ASVs=length(intersect(intersect(basiccontamASVs, batchDNAcontamASVs), batchSeqcontamASVs))) %>% 
    mutate(Group="all shared", Threshold=ThisValue)
  only_basic_dna <- data.frame(ASVs=length(intersect(basiccontamASVs, batchDNAcontamASVs) %>% setdiff(., all_shared))) %>% 
    mutate(Group="only DNAplate + basic", Threshold=ThisValue)
  only_basic_seq <- data.frame(ASVs=length(intersect(basiccontamASVs, batchSeqcontamASVs) %>% setdiff(., all_shared))) %>% 
    mutate(Group="only SeqBatch + basic", Threshold=ThisValue)
  only_dna_seq <- data.frame(ASVs=length(intersect(batchDNAcontamASVs, batchSeqcontamASVs) %>% setdiff(., all_shared))) %>% 
    mutate(Group = "only DNAplate + SeqBatch", Threshold=ThisValue)
  only_basic <- data.frame(ASVs=length(setdiff(basiccontamASVs, batchDNAcontamASVs) %>% setdiff(., batchSeqcontamASVs))) %>% 
    mutate(Group = "only Basic", Threshold = ThisValue)
  only_dna <- data.frame(ASVs=length(setdiff(batchDNAcontamASVs, basiccontamASVs) %>% setdiff(., batchSeqcontamASVs))) %>% 
    mutate(Group = "only DNAplate", Threshold=ThisValue) 
  only_seq <- data.frame(ASVs=length(setdiff(batchSeqcontamASVs, basiccontamASVs) %>% setdiff(., batchDNAcontamASVs))) %>% 
    mutate(Group = "only SeqBatch", Threshold=ThisValue)
  
  rbind(all_shared, only_basic_dna, only_basic_seq, only_dna_seq, only_basic, only_dna, only_seq)
}

## use the function to get the number of shared ASVs per threshold
tmpASVcounts01 <- commonASVcounter("0.1")
tmpASVcounts02 <- commonASVcounter("0.2")
tmpASVcounts03 <- commonASVcounter("0.3")
tmpASVcounts05 <- commonASVcounter("0.5")

## combine the datasets into a single table
ASVcount_byThreshold_byGroup <- rbind(tmpASVcounts01, tmpASVcounts02, tmpASVcounts03, tmpASVcounts05) %>% 
  dcast(Group ~ Threshold, value.var = "ASVs")
rm(list=ls(pattern = 'tmpASVcounts*'))
ASVcount_byThreshold_byGroup$orderer <- c(1,5,6,2,3,7,4)
ASVcount_byThreshold_byGroup <- ASVcount_byThreshold_byGroup %>% arrange(orderer) %>% select(-orderer)
## plot; save as 'TableOf_ASVcounts_byThreshold_byGroup'; export at 350x300
formattable(ASVcount_byThreshold_byGroup)
## Most ASVs are shared either among all three filtering types or pairs of filtering types...
## Very few ASVs are flagged as contaminants only by one filtering type (with DNAplate batch type as most)

########################################################################
## decontam ASVs - figures plotting the number of samples and number of reads per SampleType and ContamType
## exploring how batch type selection identifies contaminants or not at common threshold (0.2)
########################################################################
contam_dna02_ASVs <- all.prev %>% 
  filter(threshold=="0.2" & batch=="DNAplate" & contaminant==TRUE) %>% 
  select(ASVid)
contam_seq02_ASVs <- all.prev %>% 
  filter(threshold=="0.2" & batch=="SeqBatch" & contaminant==TRUE) %>% 
  select(ASVid)

## not sure if it's better to combine these or not?
## notrun: all02_ASVs <- c(contam_basic02_ASVs$ASVid, contam_seq03_ASVs$ASVid, contam_dna03_ASVs$ASVid) %>% unique(.)

## This plot is going to just use the ASVs identified in the "SeqBatch" dataset.. 
## More controls per SeqBatch to better model the data; too few per DNAplate (and prevalence does better with more NTCs)
## Could combine from 'all02_ASVs' instead to be more conservative
dna_contam_df_sumry <- long_df %>% 
  filter(ASVid %in% contam_dna02_ASVs$ASVid) %>%  group_by(SampleType, ASVid, ) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% mutate(Status="contam") %>% mutate(Batch="DNAplate")
dna_noncontam_df_sumry <- long_df %>% 
  filter(!ASVid %in% contam_dna02_ASVs$ASVid) %>%  group_by(SampleType, ASVid) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% mutate(Status="non_contam") %>% mutate(Batch="DNAplate")
seq_contam_df_sumry <- long_df %>% 
  filter(ASVid %in% contam_seq02_ASVs$ASVid) %>%  group_by(SampleType, ASVid, ) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% mutate(Status="contam") %>% mutate(Batch="SeqBatch")
seq_noncontam_df_sumry <- long_df %>% 
  filter(!ASVid %in% contam_seq02_ASVs$ASVid) %>%  group_by(SampleType, ASVid) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)) %>% mutate(Status="non_contam") %>% mutate(Batch="SeqBatch")

contam_plot_sumry <- rbind(dna_contam_df_sumry, dna_noncontam_df_sumry, seq_contam_df_sumry, seq_noncontam_df_sumry)
contam_plot_sumry <- contam_plot_sumry %>% filter(SampleType!="contaminant")
rm(dna_contam_df_sumry, dna_noncontam_df_sumry, seq_contam_df_sumry, seq_noncontam_df_sumry)

## plot; save as 'decontam_Reads-ASVs_ContamOrNot; export at 1200x600
ggplot(data = contam_plot_sumry, aes(x=Reads, y=Samples, color=Status, shape=SampleType)) +
  geom_point() + 
  theme_bw() + 
  annotation_logticks(sides="b") +
  scale_color_manual(values=c("firebrick", "dodgerblue2")) + 
  scale_shape_manual(values=c(15,17,1)) +
  scale_x_log10(labels=comma_format(accuracy = 1)) +
  facet_wrap(Batch ~ SampleType, scales = "free_y", nrow=2) +
  theme(legend.position = "top", legend.text = element_text(size=12), legend.title = element_blank()) +
  labs(x="sequence counts per ASV", "Samples per ASV")
## note that the biggest abundance Contamination ASVs aren't identified as likely contaminants in the SeqBatch method, but are with DNAplate method

########################################################################
## we can also plot the decontam-style figure showing the prevaelence of ncontrol vs. samples for each batch type
## labeling a few of the most prevalent ASVs that are different or similar among the batch types
########################################################################

### one way to think about prevalence: which ASVs are:
### 1) shared among bat samples and negative control samples...
### 2) that are in at least 10 negative control samples...
### for each Sample Type, how many samples do we see those ASVs?
ncontrol_ASVids <- 
  long_df %>% 
  filter(SampleType == 'ncontrol') %>% 
  distinct(ASVid) %>% pull()

sample_ASVids <-
  long_df %>% 
  filter(SampleType == 'sample') %>% 
  distinct(ASVid) %>% pull()

common_ASVids <- intersect(ncontrol_ASVids, sample_ASVids)

prevalent_common_ASVids <- long_df %>% 
  filter(ASVid %in% common_ASVids) %>% 
  group_by(ASValias, ASVid, SampleType) %>% 
  summarise(nSamples = n_distinct(SampleID),
            sumReads = sum(Reads)) %>% 
  filter(SampleType == "ncontrol") %>% 
  filter(nSamples >= 10)

prev_common_ASVs_df <- long_df %>% 
  filter(ASVid %in% prevalent_common_ASVids$ASVid) %>% 
  filter(SampleType != "contaminant") %>% 
  filter(SampleType != "mock") %>% 
  group_by(ASValias, ASVid, SampleType) %>% 
  summarise(nSamples = n_distinct(SampleID),
            sumReads = sum(Reads)) %>% 
  ungroup() %>% 
  select(ASValias, SampleType, nSamples) %>% 
  pivot_wider(values_from = "nSamples", names_from = "SampleType", values_fill = 0) %>% 
  arrange(-ncontrol)

formattable(prev_common_ASVs_df)


###################


psf.pa <- transform_sample_counts(psf, function(abund) 1*(abund>0))
psf.pa.neg <- prune_samples(sample_data(psf.pa)$is.neg == "TRUE", psf.pa)
psf.pa.pos <- prune_samples(sample_data(psf.pa)$is.neg == "FALSE", psf.pa)

# Make data.frame of prevalence in positive and negative samples for each batch type
dna.df.pa <- data.frame(pa.pos=taxa_sums(psf.pa.pos), pa.neg=taxa_sums(psf.pa.neg),
                        contaminant=all.prev %>% filter(batch=="DNAplate" & threshold=="0.2") %>% select(contaminant)) %>% 
  mutate(ASVid=row.names(.)) %>% mutate(Batch="DNAplate")

seq.df.pa <- data.frame(pa.pos=taxa_sums(psf.pa.pos), pa.neg=taxa_sums(psf.pa.neg),
                        contaminant=all.prev %>% filter(batch=="SeqBatch" & threshold=="0.2") %>% select(contaminant)) %>% 
  mutate(ASVid=row.names(.)) %>% mutate(Batch="SeqBatch")
df.pa <- rbind(dna.df.pa, seq.df.pa)
ASVkey <- long_df %>% group_by(ASVid, ASValias) %>% summarise(Samples=n_distinct(SampleID)) %>% select(-Samples)
write_csv(ASVkey, file="~/github/nhguano/data/fasta/asvkey_allASVs.csv")
df.pa <- merge(df.pa, ASVkey)
rm(dna.df.pa, seq.df.pa)

## generate ASV list to use for selected points
## ASVs 2,5,7 are TRUE for DNAplate batch, but FALSE for SeqBatch method
## ASVs 8, 14, 16 are TRUE fpr both methods
## ASVs 15, 18, 20 are FALSE for both methods.
selectASVlist <- paste("ASV-", c(2,5,7,8,14,15,16,18,20), sep="")

## plot scatterplot; save as 'decontam_Prevalence_ncontrolANDsample'; export at 700x700 (WARNING: changing dimensions messes with labels)
ggplot() + 
  geom_point(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) +
  facet_wrap(~Batch) +
  scale_color_manual(values=c("dodgerblue2", "firebrick")) +
  geom_label_repel(data=df.pa %>% filter(ASValias %in% selectASVlist & pa.neg >= 40), aes(x=pa.neg, y=pa.pos, label=ASValias),
                   force = 2, nudge_x = -10, segment.colour = "gray60", size=3) +
  geom_label_repel(data=df.pa %>% filter(ASValias %in% selectASVlist & pa.neg < 40 & pa.pos > 200), aes(x=pa.neg, y=pa.pos, label=ASValias),
                   force = 2, nudge_x = 10, nudge_y = 25, segment.colour = "gray60", size=3) +
  geom_label_repel(data=df.pa %>% filter(ASValias %in% selectASVlist & pa.neg < 40 & pa.pos < 200), aes(x=pa.neg, y=pa.pos, label=ASValias),
                   force = 2, nudge_x = 20, segment.colour = "gray60", size=3) +
  labs(x="Prevalence (Negative Controls)", 
       y="Prevalence (True Samples)",
       caption = "FALSE indicates ASV not identified as contaminant\nTRUE indicates ASV identified as contaminant by the model") +
  theme_bw() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12),
        legend.position="top", legend.text = element_text(size=12), legend.title = element_blank(),
        strip.text = element_text(size=12))

rm(psf.pa, psf.pa.neg, psf.pa.pos, df.pa)
rm(contam_df_sumry, contam_dna02_ASVs, contam_seq02_ASVs, contam_plot_sumry)

########################################################################
## Plotting selected ASVs that are highly preavlent but..
## are either both TRUE, both FALSE, or differ, depending on the ..
## Decontam batch method selected
########################################################################

## append the `long_df` so that we set "contam" or not based on the status for each batch method
seq_contamASVs <- all.prev %>% filter(batch=="SeqBatch" & threshold=="0.2" & contaminant==TRUE) %>% pull(ASVid)
dna_contamASVs <- all.prev %>% filter(batch=="DNAplate" & threshold=="0.2" & contaminant==TRUE) %>% pull(ASVid)

## plotting the per-sample abundances from the `decontam_Prevalence_ncontrolANDsample` plot:
## Using the nine ASVs selected from previous plot to compare
keepsampletypelist <- c("ncontrol", "sample")
prev_plotdat <- long_df %>% filter(ASValias %in% selectASVlist & SampleType %in% keepsampletypelist) %>% 
  select(ASValias, Reads, SampleType, ASVid)
prev_plotdat$SeqBatch <- ifelse(prev_plotdat$ASVid %in% seq_contamASVs, TRUE, FALSE)
prev_plotdat$DNAplate <- ifelse(prev_plotdat$ASVid %in% dna_contamASVs, TRUE, FALSE)
prev_plotdat$ColorMark <- ""
prev_plotdat$ColorMark <- ifelse(prev_plotdat$DNAplate==TRUE & prev_plotdat$SeqBatch==TRUE, TRUE, prev_plotdat$ColorMark)
prev_plotdat$ColorMark <- ifelse(prev_plotdat$DNAplate==FALSE & prev_plotdat$SeqBatch==FALSE, FALSE, prev_plotdat$ColorMark)
prev_plotdat$ColorMark <- ifelse(prev_plotdat$DNAplate!=prev_plotdat$SeqBatch, 'DIFF', prev_plotdat$ColorMark)

## set levels plot; ordering by 'DIFF', then 'TRUE', then 'FALSE' bor ColorMark (batch comparison result)
prev_plotdat$ASValias <- factor(prev_plotdat$ASValias, levels = c(
  "ASV-2", "ASV-5", "ASV-7", 
  "ASV-8", "ASV-14", "ASV-16", 
  "ASV-15", "ASV-18", "ASV-20"))

ggplot(prev_plotdat,
       aes(x=SampleType, y=Reads, color=ColorMark, shape=SampleType, group=SampleType)) +
  geom_jitter(width = 0.3) +
  geom_boxplot(outlier.colour = NA, fill="white", alpha=0.5) +
  facet_grid(. ~ ASValias) +
  scale_color_manual(values=c("darkorchid4", "firebrick", "dodgerblue2"), name ="Decontam class") +
  scale_shape_manual(values=c(15, 1), name = "Sample type") +
  scale_y_continuous(trans = 'log2', labels = comma_format(accuracy = 1), breaks = c(16, 512, 16384), limits = c(0.99,155000)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x="", y="sequence counts")



########################################################################
########################################################################
## Diversity analyses to compare effect of NTC community comoposition by batch type (SeqBatch or DNAplate)
########################################################################
########################################################################

## See the `decontam_workflow.md` document for QIIME 2 code executed to generate the imported .qza files

## Import PCoA artifacts from QIIME
### download directly to local machine:
# download.file(url = 'https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/decontam_pcoa/ds_pcoa.qza', destfile = '~/Desktop/ds_pcoa.qza')
# download.file(url = 'https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/decontam_pcoa/bc_pcoa.qza', destfile = '~/Desktop/bc_pcoa.qza')
# download.file(url = 'https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/decontam_pcoa/uu_pcoa.qza', destfile = '~/Desktop/uu_pcoa.qza')
# download.file(url = 'https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/decontam_pcoa/wu_pcoa.qza', destfile = '~/Desktop/wu_pcoa.qza')

ds_pcoaqza_path="~/github/nhguano/data/qiime_qza/decontam_pcoa/ds_pcoa.qza"  # set path to your local machine !!!
bc_pcoaqza_path="~/github/nhguano/data/qiime_qza/decontam_pcoa/bc_pcoa.qza"
uu_pcoaqza_path="~/github/nhguano/data/qiime_qza/decontam_pcoa/uu_pcoa.qza"
wu_pcoaqza_path="~/github/nhguano/data/qiime_qza/decontam_pcoa/wu_pcoa.qza"

## Function to generate plot data from PCoA qza artifacts
pcoa2df <- function(Path, Metric){
  tmp_pcoaqza <- read_qza(Path)
  axisPC1label <- paste(round(tmp_pcoaqza$data$ProportionExplained[1, 1], 3)*100, "%)", sep = " ")
  axisPC2label <- paste(round(tmp_pcoaqza$data$ProportionExplained[1, 2], 3)*100, "%)", sep = " ")
  plotdat <- tmp_pcoaqza$data$Vectors[, 1:3] %>% 
    mutate(PC1label=paste("PC1 (", axisPC1label, sep = " "), 
           PC2label=paste("PC2 (", axisPC2label, sep = " "),
           Metric=Metric)
  as.data.frame(plotdat)
}

dspd <- pcoa2df(ds_pcoaqza_path, "Dice-Sorensen")
bcpd <- pcoa2df(bc_pcoaqza_path, "Bray-Curtis")
uupd <- pcoa2df(uu_pcoaqza_path, "Unweighted-UniFrac")
wupd <- pcoa2df(wu_pcoaqza_path, "Weighted-UniFrac")

allplotdat <- rbind(dspd, bcpd, uupd, wupd)
rm(dspd, bcpd, uupd, wupd)
allplotdat <- merge(allplotdat, metadata)
allplotdat <- allplotdat %>% filter(SampleType != "contaminant")
allplotdat$SeqBatch <- as.character(allplotdat$SeqBatch)
coordequalval <- max(abs(c(max(allplotdat$PC1), min(allplotdat$PC1), max(allplotdat$PC2), min(allplotdat$PC2))))


########################################################################
## first plot is to group data by SeqBatch
########################################################################
## notrun: custom 9-color palette for each SeqBatch 
## notrune: pal9 <- c("#cc53a2","#8dcf55","#784bc2","#cca552","#51314f","#8bcbae","#bf5243","#9b99c0","#525f3b")

## make 4 separate plots so we can include the % variance on each axis, then stitch into a single plot
pcoaplot_bySeq <- function(BetaTest){
  ggplot() +
    geom_text(data=allplotdat %>% filter(SampleType=='sample' & Metric==BetaTest), 
              aes(x=PC1, y=PC2, label=SeqBatch), color='gray50', alpha=0.6) +
    geom_text(data=allplotdat %>% filter(SampleType=="ncontrol" & Metric==BetaTest), 
              aes(x=PC1, y=PC2, label=SeqBatch), color='#512698', size=4) +
    scale_x_continuous(limits = c(-coordequalval, coordequalval)) + 
    scale_y_continuous(limits = c(-coordequalval, coordequalval)) + 
    labs(x=allplotdat %>% filter(Metric==BetaTest) %>% pull(PC1label) %>% unique(),
         y=allplotdat %>% filter(Metric==BetaTest) %>% pull(PC2label) %>% unique(),
         subtitle = BetaTest, 
         color = 'Sequencing Batch',
         shape = 'Contaminant') +
    theme_bw()
}

## collect data for each plot
dsplot_seq <- pcoaplot_bySeq("Dice-Sorensen")
bcplot_seq <- pcoaplot_bySeq("Bray-Curtis")
uuplot_seq <- pcoaplot_bySeq("Unweighted-UniFrac")
wuplot_seq <- pcoaplot_bySeq("Weighted-UniFrac")

## plot scatterplot; save as 'contam_pcoa_4metric_bySeq'; export at 750x750
ggarrange(dsplot_seq, bcplot_seq, uuplot_seq, wuplot_seq, 
          common.legend = TRUE, ncol = 2, nrow=2, labels = c("A", "B", "C", "D"))
## pretty clear from the weighted unifrac estimate that our NTCs are as likely to be in any SeqBatch as any other

########################################################################
## Next plot will group by DNAplate 
########################################################################

pcoaplot_byDNA <- function(BetaTest){
  ggplot() +
    geom_text(data=allplotdat %>% filter(SampleType=='sample' & Metric==BetaTest), 
              aes(x=PC1, y=PC2, label=DNAplate), color='gray50', alpha=0.6) +
    geom_text(data=allplotdat %>% filter(SampleType=="ncontrol" & Metric==BetaTest), 
              aes(x=PC1, y=PC2, label=DNAplate), color='#512698', size=4) +
    scale_x_continuous(limits = c(-coordequalval, coordequalval)) + 
    scale_y_continuous(limits = c(-coordequalval, coordequalval)) + 
    labs(x=allplotdat %>% filter(Metric==BetaTest) %>% pull(PC1label) %>% unique(),
         y=allplotdat %>% filter(Metric==BetaTest) %>% pull(PC2label) %>% unique(),
         subtitle = BetaTest, 
         color = 'Sequencing Batch',
         shape = 'Contaminant') +
    theme_bw()
}

## collect data for each plot
dsplot_dna <- pcoaplot_byDNA("Dice-Sorensen")
bcplot_dna <- pcoaplot_byDNA("Bray-Curtis")
uuplot_dna <- pcoaplot_byDNA("Unweighted-UniFrac")
wuplot_dna <- pcoaplot_byDNA("Weighted-UniFrac")

## plot scatterplot; save as 'contam_pcoa_4metric_byDNA'; export at 750x750
ggarrange(dsplot_dna, bcplot_dna, uuplot_dna, wuplot_dna, 
          common.legend = TRUE,
          ncol = 2, nrow=2, labels = c("A", "B", "C", "D"))
## here we find that DNAplates tend to be associated for unweighted abundance metrics (low abundnace reads are driving similarities)
##

######################################################################
## focus solely on DNAplate 14 - the one with two of the highest read abundances for ncontrol samples
######################################################################
pcoaplot_byDNAwell <- function(BetaTest){
  ggplot() +
    geom_text(data=allplotdat %>% filter(SampleType=='sample' & Metric==BetaTest & DNAplate=="33"), 
              aes(x=PC1, y=PC2, label=DNAwell), color='gray50', alpha=0.6) +
    geom_point(data=allplotdat %>% filter(SampleType=="ncontrol" & Metric==BetaTest & DNAplate=="33"), 
               aes(x=PC1, y=PC2), color='#512698', size=3) +
    geom_label_repel(data=allplotdat %>% filter(SampleType=="ncontrol" & Metric==BetaTest & DNAplate=="33"), 
                     aes(x=PC1, y=PC2, label=DNAwell), color='#512698', size=4, nudge_y = 0.5, force=5) +
    scale_x_continuous(limits = c(-0.3, 0.32)) + 
    scale_y_continuous(limits = c(-0.2, 0.5)) + 
    labs(x=allplotdat %>% filter(Metric==BetaTest) %>% pull(PC1label) %>% unique(),
         y=allplotdat %>% filter(Metric==BetaTest) %>% pull(PC2label) %>% unique(),
         subtitle = BetaTest, 
         color = 'Sequencing Batch',
         shape = 'Contaminant') +
    theme_bw()
}

dsplot_dnawell <- pcoaplot_byDNAwell("Dice-Sorensen")
bcplot_dnawell <- pcoaplot_byDNAwell("Bray-Curtis")
uuplot_dnawell <- pcoaplot_byDNAwell("Unweighted-Unifrac")
wuplot_dnawell <- pcoaplot_byDNAwell("Weighted-Unifrac")

## plot scatterplot; save as 'contam_pcoa_4metric_byDNAwell33only'; export at 750x750
ggarrange(dsplot_dnawell, bcplot_dnawell, uuplot_dnawell, wuplot_dnawell, 
          common.legend = TRUE,
          ncol = 2, nrow=2, labels = c("A", "B", "C", "D"))
## in this case, we see that the ncontrol wells don't correlate with surrounding wells 


## clenaup:
rm(list = ls(pattern = "dsplot_*"))
rm(list = ls(pattern = "bcplot_*"))
rm(list = ls(pattern = "uuplot_*"))
rm(list = ls(pattern = "wuplot_*"))
rm(coordequalval, allplotdat)
rm(list=ls(pattern="pcoa*"))

######################################################################
## Next section is to use mock data to visualize index bleed among the 9 libraries
######################################################################## 

## subset mock data 
mock <- long_df %>% filter(SampleType=="mock")
mock$SampleID <- as.character(mock$SampleID)

## rename mock sample names for plotting clarity 
mock$SampleID <- ifelse(mock$SampleID=="mockp1L2a", "mock_batch_1.2a", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockp1L2b", "mock_batch_1.2b", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockp3L1", "mock_batch_3.1", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockp3L2", "mock_batch_3.2", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockIM4p4L1", "mock_batch_4.1", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockIM4p4L2", "mock_batch_4.2", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockIM4p5L1", "mock_batch_5.1", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockIM4p5L2", "mock_batch_5.2", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockIM4p7L1", "mock_batch_7.1", mock$SampleID)
mock$SampleID <- ifelse(mock$SampleID=="mockIM4p7L2", "mock_batch_7.2", mock$SampleID)

## add in the "expected" or "unexpected" mock seq
### download to local machine:
#download.file(url='https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/mock/mock.expectSeqs.qza', destfile = '~/Desktop/mock.expectSeqs.qza')
#download.file(url='https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/mock/mock.notexpectSeqs.qza', destfile = '~/Desktop/mock.notexpectSeqs.qza')

expect_mockqza <- read_qza(file='~/github/nhguano/data/qiime_qza/mock/mock.expectSeqs.qza') ## set your path!!
expect_mockseqs <- data.frame(expect_mockqza$data) %>% mutate(ASVid=row.names(.)) %>% pull(ASVid) %>% unique()
length(expect_mockseqs)
## 374 seqs fit within 98% identity and 97% query coverage

notexpect_mockqza <- read_qza(file='~/github/nhguano/data/qiime_qza/mock/mock.notexpectSeqs.qza') ## set your path!!!
notexpect_mockseqs <- data.frame(notexpect_mockqza$data) %>% mutate(ASVid=row.names(.)) %>% pull(ASVid) %>% unique()
length(notexpect_mockseqs)  
## 10,423 don't fit this criteria !

## Use the ASVid names in the `expect_mockseqs` to classify a mock ASVid as mock expected or not
mock$MockASV <- ifelse(mock$ASVid %in% expect_mockseqs, TRUE, FALSE)

## plot the per-sample abundances of each ASV, faceted by mock sample, coloring the Mock ASV status
## plot scatterplot, save as 'contam_mockIndexBleed'; export at 
ggplot(mock, aes(x=SampleID, y=Reads, color=MockASV)) +
  geom_jitter(width = 0.3) +
  facet_wrap(~SampleID, nrow = 2, scales = "free_x") +
  scale_y_continuous(trans="log10") +
  labs(x="", y="sequence counts\n") +
  scale_color_manual(values=c("orangered3", "chartreuse4")) +
  theme_bw() +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())

## which mock ASVs are generally only detected in high proportions?
prevalent_ASVs_df <- mock %>% filter(MockASV == TRUE & Reads > 200)
prevalent_ASVs <- prevalent_ASVs_df %>% select(ASVid) %>% unique()
colnames(prevalent_ASVs) <- "featureid"
write_delim(prevalent_ASVs, file="~/github/nhguano/data/fasta/test_prevalentMockASVs.txt.gz")


rm(prevalent_ASVs, prevalent_ASVs_df, mmockexpectqza, expect_mockqza, bigmockASVs, notexpect_mockseqs, notexpect_mockqza,
   mockexpectqza, mock, p3L2asvs, qzapath, expect_mockseqs, expectseqs, oroList)