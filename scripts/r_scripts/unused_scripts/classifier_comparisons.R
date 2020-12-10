library(tidyverse)
library(reshape2)
library(qiime2R)
library(formattable)

################################################################################
## the same dataset was classified using either:
# A. VSEARCH (97% identity, 89% coverage)
# B. Naive Bayes classifier (default settings in QIIME 2)

## First two plots explore:
#1. How often do the two classifiers agree (same taxon name per level for each ASV)
#2. How often do the two classifiers assign some information or not (per ASV, per taxon level)
################################################################################

## import ASValias values used in inital plots (keeping consistent with earlier plots in decontam work)
asvkey <- read_csv(file="~/Repos/nhguano/data/tax/asvkey_allASVs.csv")

## import VSEARCH taxonomy
vstaxa <- read_delim(file="~/Repos/nhguano/data/tax/tmp.raw_bigDB_VStax_c94p97.tsv", delim="\t")
vstaxa <- vstaxa %>% separate(., col = Taxon, 
                              sep=';', into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub(".__", "", y)))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub("^$|^ $", NA, y)))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(vstaxa)[1] <- "ASVid"
vstaxa <- merge(vstaxa, asvkey)

## import Naive Bayes taxonomy
nbtaxa <- read_delim(file="~/Repos/nhguano/data/tax/tmp.raw_bigDB_nbtax.tsv", delim="\t")
nbtaxa <- nbtaxa %>% separate(., col = Taxon, 
                              sep=';', into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub(".__", "", y)))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub("^$|^ $", NA, y)))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(nbtaxa)[1] <- "ASVid"
nbtaxa <- merge(nbtaxa, asvkey)

rm(asvkey)
## common and different values per level:
matchfunction <- function(TaxonLevel){
  vstmp <- vstaxa %>% select(ASVid, TaxonLevel) %>% filter(complete.cases(.)) %>% unite(., "vsinput", c("ASVid", TaxonLevel)) %>% pull(vsinput)
  nbtmp <- nbtaxa %>% select(ASVid, TaxonLevel) %>% filter(complete.cases(.)) %>% unite(., "nbinput", c("ASVid", TaxonLevel)) %>% pull(nbinput)
  Match <- length(intersect(vstmp, nbtmp))
  Diff <- length(setdiff(vstmp, nbtmp))
  data.frame(TaxonLevel, Match, Diff)
}

## table calculating number of times they disagree/match, but note that these are ignoring comparisons where one might be NA
CompTable <- rbind(matchfunction("Kingdom"), matchfunction("Phylum"), matchfunction("Class"),
                   matchfunction("Order"), matchfunction("Family"), matchfunction("Genus"), matchfunction("Species"))

## plot table as image; save as 'taxcomp_matchStatus'; export at 275x300
formattable(CompTable)

## how often does a classifier assign a name or not?
emptyfunction <- function(TaxonLevel){
  vstmp <- vstaxa %>% select(ASVid, TaxonLevel) %>% rename(vtax=TaxonLevel) 
  nbtmp <- nbtaxa %>% select(ASVid, TaxonLevel) %>% rename(ntax=TaxonLevel)
  alltmp <- merge(vstmp, nbtmp)
  alltmp$Status <- ""
  alltmp$Status <- ifelse(is.na(alltmp$vtax) & is.na(alltmp$ntax), "both missing", alltmp$Status)
  alltmp$Status <- ifelse(is.na(alltmp$vtax) & !is.na(alltmp$ntax), "vsearch missing only", alltmp$Status)
  alltmp$Status <- ifelse(!is.na(alltmp$vtax) & is.na(alltmp$ntax), "nbayes missing only", alltmp$Status)
  alltmp$Status <- ifelse(!is.na(alltmp$vtax) & !is.na(alltmp$ntax), "neither missing", alltmp$Status)
  alltmp %>% group_by(Status) %>% tally() %>% spread(Status, n) %>% mutate(Level=TaxonLevel)
}

EmptyTable <- rbind(emptyfunction("Phylum"), emptyfunction("Class"),emptyfunction("Order"), 
                    emptyfunction("Family"), emptyfunction("Genus"), emptyfunction("Species"))
EmptyTable <- EmptyTable %>% select(Level, `both missing`, `neither missing`, `nbayes missing only`, `vsearch missing only`)
tmpk <- data.frame(emptyfunction("Kingdom"), "nbayes missing only"=0)
colnames(tmpk) <- c("both missing", "neither missing", "vsearch missing only", "Level", "nbayes missing only")
tmpk <- tmpk %>% select(Level, `both missing`, `neither missing`, `nbayes missing only`, `vsearch missing only`)
EmptyTable <- rbind(tmpk, EmptyTable)

## plot table as image; save as 'taxcomp_missingStatus'; export at 650x350
formattable(EmptyTable, 
            caption="Instances where one or both classifiers do not assign name to an ASV")

rm(CompTable, EmptyTable, tmpk)


################################################################################
## Classifier specifics: among our most prevalent ASVs, which are assigned names and which aren't?
################################################################################

## import metadata again
metadata <- read_csv(file="~/Repos/nhguano/data/metadata/allbat_meta.csv")
metadata$is.neg <- ifelse(metadata$SampleType=="ncontrol", TRUE, FALSE)
metadata <- as.data.frame(metadata)

## import ASV-filtered, no mock/ncontrol data this time
filtqzapath="/Users/do/Repos/nhguano/data/qiime_qza/ASVtable/sampleOnly_nobatASV_table.qza"
features <- read_qza(filtqzapath)
mat.tmp <- features$data
rm(features)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
flong_df <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp)
colnames(flong_df) <- c("ASVid", "SampleID", "Reads")
flong_df <- merge(flong_df, metadata, all.x = TRUE)
rm(metadata)

## merge separate dataasets for each taxonomy file:
vs_long <- merge(flong_df, vstaxa)
nb_long <- merge(flong_df, nbtaxa)

## How many ASVs contain at least Order/Family/Genus information?
vs_asvsumry <- as.data.frame(vs_long %>% group_by(ASVid, ASValias, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)))
nb_asvsumry <- as.data.frame(nb_long %>% group_by(ASVid, ASValias, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% summarise(Samples=n_distinct(SampleID), Reads=sum(Reads)))
## it looks like among the most prevalent ASVs (which are typically also those generating the most reads) both VSEARCH and Naive Bayes agree
## however there are a few instances in which Naive Bayes provides additional information where VSEARCH does not
## would be great to work out these differences in a hybrid classifier but need to test first before using on real dataset
## would need to include a least common ancestor approach for instances where VSEARCH and Naive Bayes disagree at some level

## Because there are a few highly prevalent ASVs that Naive Bayes appears to classify that VSEARCH missed, let's gather those

## Gathering ASVs that have less than Family-name information in VSEARCH method 
vs_miss_df <- as.data.frame(vs_asvsumry %>% filter(is.na(Family)))
  ## almost all of these contain ZERO information (Undefined through Class rank)
  ## those that contain Order information only tend to be lepidopterans, which we have found tend to cluster a lot with this amplicon
  ## these are often cases in which multipe Families are good hits within 97% identity, so the name gets collapsed to just the Order

## How many Naive Bayes-classified ASVs have some information for those VSEARCH ASVs that lack at least Family-name information?
nb_vsmiss_df <- as.data.frame(nb_asvsumry %>% filter(ASVid %in% vs_miss_df$ASVid))
  ## the most prevalent ASVs are frequently just one of a handful of taxa (ex. Phyllophaga hirsuta)..
  ## ..and previously unclassified taxa (by VSEARCH) now have information 

## How many of these ASVs aren't classified to at least Family-rank within Naive Bayes data too?
nb_vsmiss_df[3:8] <- lapply(nb_vsmiss_df[3:8], as.character)
nb_vsmiss_df <- nb_vsmiss_df %>% replace(is.na(.), "Undefined")
nb_vsmiss_withOrder_df <- nb_vsmiss_df %>% filter(Order != "Undefined")
  ## lots are. and we're retaining most of the reads; other ASVs we're missing tend to be in fewer samples or with lower total read coutns

rm(nb_vsmiss_df, nb_vsmiss_withOrder_df, vs_miss_df)
################################################################################
## Semi-hybrid classification
## Trusting the majority of VSEARCH information because it retains the LCA information when multiple best hits available
## (and is therefore more conservative)
## 1. Identify all taxa with at least Order-level information assigned by VSEARCH; keep these ASVs classified as is
## 2. For all other ASVs, replace with Naive Bayes information 
## 3. Identify if any remaining taxa that are highly prevalent are missing and search BOLD database manually
################################################################################

vsearch_asvs2keep <- vs_asvsumry %>% filter(!is.na(Family)) %>% pull(ASVid)
nbayes_asvs2use <- nb_asvsumry %>% filter(!ASVid %in% vsearch_asvs2keep) %>% filter(!is.na(Family)) %>% pull(ASVid)

## combine these datasets together:
vsearch_taxa_tmp <- vs_long %>% filter(ASVid %in% vsearch_asvs2keep)
nbayes_taxa_tmp <- nb_long %>% filter(!ASVid %in% vsearch_asvs2keep) %>% filter(!is.na(Family))
used_taxa_tmp <- rbind(vsearch_taxa_tmp, nbayes_taxa_tmp)
usedASVs <- used_taxa_tmp %>% pull(ASVid) %>% unique()
unused_taxa_tmp <- flong_df %>% filter(!ASVid %in% usedASVs)

## what kind of taxa are still unclassified? 
unused_taxa_sumry <- unused_taxa_tmp %>% group_by(ASVid) %>% 

unused_asvs <- vs_asvsumry %>% filter(!ASVid %in% vsearch_asvs2keep) %>% filter(!ASVid %in% nbayes_asvs2use) %>% pull(ASVid)
  ## 2766 ASVs still lack sufficient (Order name) information 

## how many of the remaining "unused" ASVs have many reads and/or high prevalence? 
unusedASV_sumry <- flong_df %>% filter(ASVid %in% unused_asvs) %>% group_by(ASVid) %>% summarise(Samples=n_distinct(SampleID), Reads=sum(Reads))
