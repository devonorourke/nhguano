library(tidyverse)
library(reshape2)
library(qiime2R)

########################################################################
## data imports
########################################################################

## import metadata
meta <- read_delim(file="~/Repos/nhguano/data/metadata/nhbat_meta.tsv", delim = "\t")
meta$SampleID <- ifelse(meta$SampleID=="negoro35A01", "negoro35A01-contaminated", meta$SampleID)
meta$SampleType <- ifelse(meta$SampleType=="contaminant", "ncontrol", meta$SampleType)

## import ASV table
qzaimport.function <- function(qzapath){
  featuretable <- read_qza(qzapath)
  mat.tmp <- featuretable$data
  rm(featuretable)
  df.tmp <- as.data.frame(mat.tmp)
  rm(mat.tmp)
  df.tmp$OTUid <- rownames(df.tmp)
  rownames(df.tmp) <- NULL
  tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
  colnames(tmp) <- c("ASVid", "SampleID", "Reads")
  tmp
}

rawtable="/Users/do/Repos/nhguano/data/qiime_qza/ASVtable/all.raw_table.qza"
rawreads <- qzaimport.function(rawtable)
rm(rawtable, qzaimport.function)
## add ASValias; because all rarefied ASVs are present in nonrarefied set, we'll use the ..
## order of the ASV by total read abundance to set the labels (so ASV-1 has the most total reads, ASV-2 the second most, etc.):
tmp1 <- rawreads %>% group_by(ASVid) %>%  summarise(nReads=sum(Reads)) %>% arrange(-nReads) %>% mutate(ASValias=paste0("ASV-", row.names(.))) %>% select(-nReads)
rawreads <- merge(rawreads, tmp1)
rm(tmp1)

## merge with metadata:
nonfyd_reads <- merge(rawreads, meta, by='SampleID')
rm(meta, rawreads)


########################################################################
## quick summary of fraction of samples that are ncontrols vs. samples
########################################################################

## How many unique samples per Study type?
nonfyd_reads %>% 
  group_by(StudyID) %>% 
  summarise(nSamples=n_distinct(SampleID))

## How many samples are negative controls?
nonfyd_reads %>%
  group_by(SampleType) %>% 
  summarise(nSamples=n_distinct(SampleID)) %>% 
  mutate(pSamples=round(nSamples/sum(nSamples), 4) * 100)
    ## so controls make up ~6% of the total samples...

## What fraction of these reads per study type are ncontrols vs. samples?
nonfyd_reads %>% 
  group_by(SampleType) %>% 
  summarise(nSamples=n_distinct(SampleID), nRead=sum(Reads)) %>% 
  mutate(pRead = round(nRead/sum(nRead),4) * 100)
  ## negative controls make up ~0.4% of overall reads (despite making up for 6% of all samples)... good!

########################################################################
## how many of these ASVs are from host sequences?
## importing data from NaiveBayes classifier using db consisting only of host data
########################################################################

hosttax <- read_delim(file="~/Repos/nhguano/data/host/study.raw_hostDB_VStax.tsv",
                      delim="\t") 
hosttax <- hosttax %>% 
  separate(., 
           col = Taxon, 
           sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"))
hosttax <- as.data.frame(apply(hosttax, 2, function(y) gsub(".__", "", y)))
hosttax <- as.data.frame(apply(hosttax, 2, function(y) gsub("^$|^ $", NA, y)))
colnames(hosttax)[1] <- "ASVid"

hostASVs <- hosttax %>% filter(phylum_name == "Chordata") %>% select(ASVid, genus_name, species_name)
dirtyASVs <- merge(dirtyASVs, hostASVs, all.x = TRUE)
  ## no ncontrols have any suspected host ASVs
expectedASVs <- merge(expectedASVs, hostASVs, all.x = TRUE)
  ## the South American bats are persistent in low abundances only
  ## the North American bats are not frequent, but are  abundant

########################################################################
## Plot that illustrates that most of our samples have lots of ASVs and lots of reads..
## ..while ncontrols typically have very few
########################################################################

alldat <- nonfyd_reads %>% group_by(SampleID, SampleType) %>% 
  summarise(Reads=sum(Reads), ASVs=n_distinct(ASVid)) %>% 
  mutate(log2Reads=log2(Reads)) %>% 
  mutate(log10Reads=log10(Reads))

## Notice here how the majority of ncontrols have fewer reads and many fewer ASVs
## By subsampling at 6,000 sequences, we retain the most of our true samples and only a small number of NTCs to deal with
## plot; save as 'contam_eval_allSampls_Counts-ASVs_scatterplot'; export at 
ggplot(data=alldat, aes(x=Reads, y=ASVs)) +
  geom_point(data=alldat %>% filter(SampleType=="sample"), aes(x=Reads, y=ASVs), color="gray75") +
  geom_point(data=alldat %>% filter(SampleType=="contaminant"), aes(x=Reads, y=ASVs), color="purple", size=2) +
  geom_point(data=alldat %>% filter(SampleType=="ncontrol"), aes(x=Reads, y=ASVs), color="dodgerblue", size=2) +
  geom_point(data=alldat %>% filter(SampleType=="mock"), aes(x=Reads, y=ASVs), color="sienna2", size=3) +
  scale_x_continuous(trans='log2', labels=comma, breaks=c(4,32,256,2048,16384,131072,1048526)) +
  geom_vline(xintercept = 6000, linetype="dashed") +
  theme_bw() +
  labs(x="sequence counts", y="ASVs") +
  annotation_logticks(sides="b")

##################
## Remove samples with low abundances??
## How many true samples remain? How many NTCs remain?
## Focus only on NH guano samples
##################

nhdat <- nonfyd_reads %>% 
  group_by(SampleID, SampleType) %>% 
  summarise(nASVs=n_distinct(ASValias), nReads=sum(Reads))

## Drop all samples with less than N reads; require at least 2 ASVs
## 500 reads
nhdat %>% filter(nReads > 500 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
  ## 1791 samples, 62 NTCs

## 1000 reads
nhdat %>% filter(nReads > 1000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
  ## 1655 samples, 39 NTCs

## 2000 reads
nhdat %>% filter(nReads > 2000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
## 1535 samples, 27 NTCs

## 3000 reads
nhdat %>% filter(nReads > 3000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
## 1464 samples, 23 NTCs

## 4000 reads
nhdat %>% filter(nReads > 4000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
## 1410 samples, 19 NTCs

## 5000 reads
nhdat %>% filter(nReads > 5000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
## 1362 samples, 17 NTCs

## 6000 reads
nhdat %>% filter(nReads > 6000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
## 1316 samples, 10 NTCs

## 7000 reads
nhdat %>% filter(nReads > 7000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
## 1280 samples, 6 NTCs

## 10000 reads
nhdat %>% filter(nReads > 10000 & nASVs > 1) %>% group_by(SampleType) %>% summarise(nSamples=n_distinct(SampleID))
## 1190 samples, 2 NTCs

########################################################################
## Of the NTC's remaining...
## how many per SeqBatch?
########################################################################

## for at least 5,000 reads/sample... what are the SampleIDs
tmp <- nhdat %>% filter(nReads > 5000 & nASVs > 1 & SampleType=="ncontrol")
dirtysamps <- as.character(tmp$SampleID)
## and which ASVs do they contain
dirtyASVs <- nonfyd_reads %>% 
  filter(SampleID %in% dirtysamps) %>% 
  group_by(ASVid) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads))
## how many of these ASVs are seen in true samples?
## first, get list of expected samples with similar read abundance threshold
tmp2 <- nhdat %>% filter(nReads > 5000 & SampleType == "sample")
expectedSamps <- as.character(tmp2$SampleID)
## then get list of ASVs and compare how many are in common with ncontrols
expectedASVs <- nonfyd_reads %>% 
  filter(SampleID %in% expectedSamps) %>% 
  group_by(ASVid) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads))
commonDirtyASVs <- intersect(dirtyASVs$ASVid, expectedASVs$ASVid)
  ## 242 in common...
rm(tmp, tmp2)

ggplot(dirtyASVs, aes(x=)) +
  geom_point(data=expectedASVs %>% filter(!ASVid %in% commonDirtyASVs), aes(x=Samples, y=Reads), color="black") +
  geom_point(data=expectedASVs %>% filter(ASVid %in% commonDirtyASVs), aes(x=Samples, y=Reads), color="red") +
  geom_point(data=dirtyASVs %>% filter(!ASVid %in% commonDirtyASVs), aes(x=Samples, y=Reads), color="blue") +
  geom_point(data=dirtyASVs %>% filter(ASVid %in% commonDirtyASVs), aes(x=Samples, y=Reads), color="orange") +
  theme_bw() +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,500000))
