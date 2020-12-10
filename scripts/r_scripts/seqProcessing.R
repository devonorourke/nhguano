library(tidyverse)
library(qiime2R)
library(magicfor)
library(SRS)
library(reshape2)
library(lubridate)

## goals in filtering with this script are to:
## 1. confirm that we can safely drop all mock-associated OTUs
## 2. drop all mock and negative control samples from dataset
## 3. identify all bat-host related OTUs from a series of classifier outputs and provide species IDs where possible
## 4. filter our raw dataset to retain only family-level Arthropod-classified OTUs, then create a new OTU table for diversity analyses
## 5&6. summarize the read depth of these remaining filtered taxa to determine sampling depth for subsampling for diversity estimates
## 7. filter the final samples at a particular read depth 

########################################
## part 1 - bring in the ASV table and .uc files to create a new OTU table (summarizing by ASVs grouped as centroids/OTUs from .uc fle)
########################################

# ## import abundance info
# ## import the original ASV table (not clustered):
  ## change path below to your own path!
    # download.file(url = "https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/ASVtable/tmp.raw_table.qza",
    #               destfile = "~/Desktop/tmp.raw_table.qza")
      
qzapath = "~/github/nhguano/data/qiime_qza/ASVtable/tmp.raw_table.qza"  ## amend this line to your own path!
featuretable <- read_qza(qzapath)
mat.asv <- featuretable$data
rm(featuretable, qzapath)

# ## convert the mat.tmp object into a long format
df.asv <- as.data.frame(mat.asv)
rm(mat.asv)
df.asv$ASVid <- rownames(df.asv)
rownames(df.asv) <- NULL
asv_table_long <- melt(df.asv, id = "ASVid") %>% filter(value != 0)
rm(df.asv)
colnames(asv_table_long) <- c("ASVid", "SampleID", "Reads")
#
# ## import the .uc file and reformat to identify the ASVs with their proper CentroidIDs
UC_df <- read_delim(file="https://github.com/devonorourke/nhguano/raw/master/data/ucfile/allSamps_clustered_p985.uc.gz",
                    delim = "\t", col_names = FALSE)
UC_df <- UC_df[,c(1,9,10)]
colnames(UC_df) <- c("RecordType", "ASVid", "CentroidID")
tmp.UC.hits <- UC_df %>%
  filter(RecordType == "H") %>%
  select(-RecordType)
tmp.UC.sentroids <- UC_df %>%
  filter(RecordType == "S") %>%
  mutate(CentroidID = ASVid) %>%
  select(-RecordType)
uc_data <- rbind(tmp.UC.sentroids, tmp.UC.hits) %>%
  separate(col = ASVid,
           into = c("ASVid", "ASVcounts"),
           sep = ";") %>%
  separate(col = CentroidID,
           into = c("OTUid", "OTUcounts"),
           sep = ";") %>%
  select(ASVid, OTUid)

rm(tmp.UC.hits, tmp.UC.sentroids, UC_df)

## merge with `asv_table_long`
asv_table_long_wCentroids <- merge(asv_table_long, uc_data, by="ASVid")

# ## group by the centroids, summing all ASVs that share a particular centroid
otu_table_long <- asv_table_long_wCentroids %>%
  group_by(OTUid, SampleID) %>%
  summarise(OTUreads = sum(Reads))
#
# ## write as long format:
write_csv(otu_table_long, quote=FALSE,
          path = "~/github/nhguano/data/text_tables/otu_tables/allSamples_OTUtable_long.csv")
# 
# ########################################
# ## part 2 - filtering mock and negative control samples
# ########################################
# 
# ## filtering out any OTUs associated with the 25 mock community ASVs
#   ## is it save to assume all OTUs with mock aren't absorbing a bunch of real diet components with similar % identity?
# 
## import mockASVs
mock_asvs <- read_delim(file = 'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/fasta/prevalentMockASVs.txt.gz',
           delim = '\t')
# ## filter the uc file to retain just those mock ASVs and their associated clusters
mock_uc_data <- uc_data %>%
  filter(ASVid %in% mock_asvs$featureid)
# ## are they all the centroids?
mock_uc_data$boolMatch <- mock_uc_data$ASVid == mock_uc_data$OTUid
  ## nope. 23 of 25 are though! good to see that our exact mock variants represent almost all of the OTU clusters.
  ## if centroids were something other than mock ASVs, we might suspect the same taxa were being consumed in our dataset (but they're not)

# ## For the boolean match == 'FALSE' term above... there are 2 such instances where the mock ASV is NOT the OTU centroid...
#   ##... how many times are those ASV observed:
  ASVs_notMockCentroids <- mock_uc_data %>% filter(boolMatch == FALSE) %>% select(OTUid)
#     ## 1) in mock samples?
    asv_table_long %>%
      filter(str_detect(SampleID, "^mock")) %>%
      filter(ASVid %in% ASVs_notMockCentroids$OTUid) %>%
      group_by(ASVid) %>%
      summarise(nSamples = n(),
                perASVreads = sum(Reads))
#         ## detected in 5 mock samples each, with about 1000 reads total per ASV...
#     ## 2) in true samples?
    asv_table_long %>%
      filter(!str_detect(SampleID, "^mock")) %>%
      filter(!str_detect(SampleID, "^neg")) %>%
      filter(ASVid %in% ASVs_notMockCentroids$OTUid) %>%
      group_by(ASVid) %>%
      summarise(nSamples = n(),
                perASVreads = sum(Reads))
#       ## detected in 1-2 samples each, with 6 or 12 total reads.
#     ## the fact that we see more of these 2 ASVs in the mocks than true samples is an indication ...
      ##... they are due to sequencing errors, or possibly if real are present only in the mock samples (and they get filtered out once we drop mock samples)
# 
# ## among these 25 mock ASVs, how many reads are there among the true samples in our dataset?
sumry_mockASVs_trueSamples <- asv_table_long %>%
  filter(!str_detect(SampleID, "^mock")) %>%
  filter(!str_detect(SampleID, "^neg")) %>%
  filter(ASVid %in% mock_uc_data$ASVid) %>%
  group_by(ASVid) %>%
  summarise(sumASVreads = sum(Reads),
            nSamples = n())
#   ## a bit of cross talk as expected... some ASVs in tens of samples, but overall read depth is a fraction of overall abundance (per mock ASV)
#   ## dropping these isn't biasing strongly to one kind of taxa, and is better to drop than retain and inflate richness values.
# 
# ## generate a list of mock-associated OTUs to drop from dataset
mock_OTU_list_toDrop <- mock_uc_data %>% select(OTUid) %>% pull()
# 
# ## drop thos mock OTUs from dataset, along with all negative control samples:
# ## import metadata to make sample filtering easy:
allMeta <- read_csv("https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/allbat_meta.csv")
# ## retain only samples relevant to this NH study:
NHsites <- c("HOP", "WLT", "GRN", "SWZ", "FOX", "WLD", "BRN", "MAP", "MAS", "GIL", "COR", "CNA", "ROL", "PEN", "EPS", "CNB", "HOL", "CHI", "MTV", "ALS")
NHsamples <- allMeta %>%
  filter(Site %in% NHsites) %>%
  filter(Date != "unknown") %>%
  select(SampleID) %>%
  pull()
# 
# ## drop the OTUs associated with the Mock ASVs/OTUs, and keep only those samples that were from NH study
otu_table_long_sampOnly <- otu_table_long %>%
  filter(!OTUid %in% mock_OTU_list_toDrop) %>%
  filter(SampleID %in% NHsamples) %>%
  filter(!grepl("contam", SampleID))

n_distinct(otu_table_long_sampOnly$SampleID)
#   ## 2,521 samples remaining with at least 1 read
# 
otu_table_long_sampOnly %>%
  group_by(SampleID) %>%
  summarise(sumReads = sum(OTUreads)) %>%
  filter(sumReads > 500) %>%
  nrow()
#   ## but only 1,317 of these have at least 500 reads! ... and of these reads, many can be bat related
# 
rm(asv_table_long, asv_table_long_wCentroids, otu_table_long,
   ASVs_notMockCentroids, mock_asvs, mock_uc_data, sumry_mockASVs_trueSamples, uc_data,
   mock_OTU_list_toDrop, NHsamples, NHsites)


# ########################################
# ## part 3 - identify which set of OTUs/ASVs were derived from what bat species
# ########################################
# ## make list of remaining OTUs present in true samples:
sampOTUs <- otu_table_long_sampOnly %>% distinct(OTUid) %>% pull()
# 
# ## import taxa files from BOLD database using NaiveBayes (NB) or VSEARCH (VS)
allSamps_NBtaxa_boldDB <- read_delim(file = "~/github/nhguano/data/taxonomy/NBtaxonomy.tsv.gz",
                              delim = "\t") %>%
  filter(`Feature ID` %in% sampOTUs) %>%
  rename(OTUid = `Feature ID`)
# 
allSamps_VStaxa_boldDB <- read_delim(file = "~/github/nhguano/data/taxonomy/VStaxonomy.tsv.gz",
                              delim = "\t") %>%
  filter(`Feature ID` %in% sampOTUs) %>%
  rename(OTUid = `Feature ID`)
# 
# ## find the bat reads assigned by either method:
NBbatOTUs <- allSamps_NBtaxa_boldDB %>%
  filter(grepl("Chiroptera", Taxon)) %>%
  select(OTUid)
nrow(NBbatOTUs)
#   ## 6 OTUs detected
# 
VSbatOTUs <- allSamps_VStaxa_boldDB %>%
  filter(grepl("Chiroptera", Taxon)) %>%
  select(OTUid)
nrow(VSbatOTUs)
#   ## 3 OTUs detected
# 
# ## are the 3 OTUs in VSEARCH a subset of those in the NaiveBayes classifier output?
intersect(VSbatOTUs$OTUid, NBbatOTUs$OTUid)
#   ## yep.
# 
# ## How many samples are these bat-classified taxa assigned to? How many total reads?
# ## generate a summary for both methods collectively.
# ##### first, using the NB classifier data:
NBbatOTUs_long <- otu_table_long_sampOnly %>%
  filter(OTUid %in% NBbatOTUs$OTUid)
NBbatOTUs_taxa <- allSamps_NBtaxa_boldDB %>%
  filter(OTUid %in% NBbatOTUs_long$OTUid) %>%
  select(-Confidence) %>%
  separate(col = Taxon, sep = ";",
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  select(OTUid, Species)
NBbatOTUs_long_wTaxa <- merge(NBbatOTUs_long, NBbatOTUs_taxa, by="OTUid")
NBbat_OTUcounts_wTaxa <- NBbatOTUs_long_wTaxa %>%
  group_by(SampleID, Species) %>%
  summarise(perbatSpeciesReadSums = sum(OTUreads)) %>%
  pivot_wider(names_from = "Species", values_from = "perbatSpeciesReadSums")
# ## how many of these are M. lucifugus?
NBbatOTUs_long_wTaxa %>% filter(Species == "s__Myotis lucifugus") %>% summarise(n_distinct(SampleID))
#   ## 579 samples have at least one OTU classified to MYLU
## which samples were these?
NBbatOTUs_sampleIDs <- NBbatOTUs_long_wTaxa %>% distinct(SampleID)
write_csv(NBbatOTUs_sampleIDs, "~/github/nhguano/data/taxonomy/NBclassifier_SampleIDs_wMYLU_assigned.csv")

# ##### second, using the VS classifier data:
VSbatOTUs_long <- otu_table_long_sampOnly %>%
  filter(OTUid %in% VSbatOTUs$OTUid)
VSbatOTUs_taxa <- allSamps_VStaxa_boldDB %>%
  filter(OTUid %in% VSbatOTUs_long$OTUid) %>%
  select(-Consensus) %>%
  separate(col = Taxon, sep = ";",
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  select(OTUid, Species)
VSbatOTUs_long_wTaxa <- merge(VSbatOTUs_long, VSbatOTUs_taxa, by="OTUid")
VSbat_OTUcounts_wTaxa <- VSbatOTUs_long_wTaxa %>%
  group_by(SampleID, Species) %>%
  summarise(perbatSpeciesReadSums = sum(OTUreads)) %>%
  pivot_wider(names_from = "Species", values_from = "perbatSpeciesReadSums")
# ## how many of these are M. lucifugus?
VSbatOTUs_long_wTaxa %>% filter(Species == "s__Myotis lucifugus") %>% summarise(n_distinct(SampleID))
#   ## 578 samples have at least one OTU classified to MYLU
## which samples were these?
VSbatOTUs_sampleIDs <- VSbatOTUs_long_wTaxa %>% distinct(SampleID)
write_csv(VSbatOTUs_sampleIDs, "~/github/nhguano/data/taxonomy/VSclassifier_SampleIDs_wMYLU_assigned.csv")


# #### we have a strong consensus from both classification methods that the vast majority...
#   ## ... of our bat COI fragments in our samples come from one species: Myotis lucifugus
# 
# ## generate summary table reporting these results:
vs_bat_sumry_tmp <- VSbatOTUs_long_wTaxa %>%
  group_by(Species) %>%
  summarise(distinctOTUs = n_distinct(OTUid),
            nReads_perSpecies = sum(OTUreads),
            nSamples_withOTU = n()) %>%
  pivot_longer(-Species, names_to = "Measure", values_to="Value") %>%
  mutate(ClassifierMethod = "vsearch")
# 
nb_bat_sumry_tmp <- NBbatOTUs_long_wTaxa %>%
  group_by(Species) %>%
  summarise(distinctOTUs = n_distinct(OTUid),
            nReads_perSpecies = sum(OTUreads),
            nSamples_withOTU = n()) %>%
  pivot_longer(-Species, names_to = "Measure", values_to="Value") %>%
  mutate(ClassifierMethod = "sklearn")
# 
batDataSumry <- rbind(vs_bat_sumry_tmp, nb_bat_sumry_tmp) %>%
  mutate(Species = case_when(Species == "s__Lasiurus borealis" ~ "Lasiurus borealis",
                             Species == "s__Myotis lucifugus" ~ "Myotis lucifugus",
                             Species == "s__Noctilio albiventris" ~ "Noctilio albiventris",
                             Species == "s__Rhogeessa tumida" ~ "Rhogeessa tumida")) %>%
  pivot_wider(names_from = "Species", values_from = "Value", values_fill = 0) %>%
  arrange(ClassifierMethod, Measure)
# 
# #write_csv(batDataSumry, path = "~/github/nhguano/supplementaryData/tableS4_batHost_summary.csv")
# 
rm(batDataSumry, nb_bat_sumry_tmp, vs_bat_sumry_tmp, NBbat_OTUcounts_wTaxa, VSbat_OTUcounts_wTaxa)
rm(list=ls(pattern="NBbatOTUs*"))
rm(list=ls(pattern="VSbatOTUs*"))
# 
# ########################################
# ## part 4 - filter our dataset to retain just those samples with at least family-level Arthropod classified OTUs
# ########################################
# 
# ## filter to retain only family-level Arthropod taxa for NB, among remaining Samples
arthOnly_NBtaxa <- allSamps_NBtaxa_boldDB %>%
  select(-Confidence) %>%
  filter(Taxon != "Unassigned") %>%
  separate(col = Taxon,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  mutate(Kingdom = gsub("tax=k__", "", Kingdom),
         Phylum = gsub("p__", "", Phylum),
         Class = gsub("c__", "", Class),
         Order = gsub("o__", "", Order),
         Family = gsub("f__", "", Family),
         Genus = gsub("g__", "", Genus)) %>%
  filter(Phylum == "Arthropoda") %>%
  mutate_all(list(~na_if(.,""))) %>%
  filter(!is.na(Family)) %>%
  filter(OTUid %in% sampOTUs)
# 
# ## filter to retain only family-level Arthropod taxa for VS
arthOnly_VStaxa <- allSamps_VStaxa_boldDB %>%
  select(-Consensus) %>%
  filter(Taxon != "Unassigned") %>%
  separate(col = Taxon,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  mutate(Kingdom = gsub("tax=k__", "", Kingdom),
         Phylum = gsub("p__", "", Phylum),
         Class = gsub("c__", "", Class),
         Order = gsub("o__", "", Order),
         Family = gsub("f__", "", Family),
         Genus = gsub("g__", "", Genus)) %>%
  filter(Phylum == "Arthropoda") %>%
  mutate_all(list(~na_if(.,""))) %>%
  filter(!is.na(Family)) %>%
  filter(OTUid %in% sampOTUs)
# 
# ## retain exact matches from VSEARCH first, and add in remaining taxa from naive Bayes classifier for remaining taxa:
vs_keepers_OTUs <- arthOnly_VStaxa %>% select(OTUid) %>% pull()
nb_keepers_toAdd <- arthOnly_NBtaxa %>%
  filter(!OTUid %in% vs_keepers_OTUs) %>%
  mutate(ClassifierMethod = "sklearn")
# 
# ## generate a final taxonomy dataset from these two OTUid lists:
arthOnly_mergedTaxa <- arthOnly_VStaxa %>%
  mutate(ClassifierMethod = "vsearch")
arthOnly_mergedTaxa <- rbind(nb_keepers_toAdd, arthOnly_mergedTaxa)
#   ## 3,186 OTUs remain

# 
# ### minor bit of reformatting/cleanup with these names:
# ## A. Two instances where Genus labels are actually ambiguous (that I could find) - both have numbers, so we'll drop those and leave as NA's
# ## B. Multiple instances where there are ambiguous things like 'sp.' or 'nr.' or 'gr.' in Species labels; drop those and replace with NA
# ## C. Multiple instnaces where other arbritrary labels are appended to things like "Janzen"; drop those too and replace with NA
# ## D. rename one odd instance where extra info added to species label (converting 'Coleophora duplicis group' --> 'Coleophora duplicis')
# ## E. rename a few other instances where species label is clear, but has additional number that isn't necessary
# 
arthOnly_mergedTaxa_cleaned <- arthOnly_mergedTaxa %>%
  mutate(Genus = ifelse(grepl('[0-9]', Genus), NA, Genus)) %>%
  mutate(Species = ifelse(grepl("s__Enicospilus purgatusDHJ02", Species), "Enicospilus purgatus", Species)) %>%
  mutate(Species = ifelse(grepl("s__Fissicrambus mutabilis PS3", Species), "Fissicrambus mutabilis", Species)) %>%
  mutate(Species = ifelse(grepl("s__Enicospilus flavostigmaDHJ658", Species), "Enicospilus flavostigma", Species)) %>%
  mutate(Species = ifelse(grepl("s__Phyllophaga elenansAS1", Species), "Phyllophaga elenans", Species)) %>%
  mutate(Species = ifelse(grepl("s__Chironomus nr. atroviridis 2i IP2013", Species), "Chironomus atroviridis", Species)) %>%
  mutate(Species = ifelse(grepl("s__Microthyris prolongalisDHJ02", Species), "Microthyris prolongalis", Species)) %>%
  mutate(Species = ifelse(grepl('[0-9]', Species), NA, Species)) %>%
  mutate(Species = ifelse(grepl('sp\\.', Species), NA, Species)) %>%
  mutate(Species = ifelse(grepl('nr\\.', Species), NA, Species)) %>%
  mutate(Species = ifelse(grepl('gr\\.', Species), NA, Species)) %>%
  mutate(Species = gsub("s__$", NA, Species)) %>%
  mutate(Species = gsub("s__", "", Species))
#          
# 
# ## merge this information with the long OTU table object and write to disk:
otu_table_long_sampOnly_wTaxa <- merge(otu_table_long_sampOnly, arthOnly_mergedTaxa_cleaned, by="OTUid")
# 
# 
# #################### sanity check == start ####################
# ##### sanity check: are any of these OTUs we _failed_ to classify present in many samples? with high read depth?
allSamps_NBtaxa_readCounts <- merge(allSamps_NBtaxa_boldDB, otu_table_long_sampOnly, by="OTUid")
allSamps_NBtaxa_OTUsumry <- allSamps_NBtaxa_readCounts %>%
  group_by(OTUid, Taxon) %>%
  summarise(nSamples = n(),
            totalReads = sum(OTUreads))
# ## what OTUs were left out from our final table that are present here?
NB_droppedOTUs <- setdiff(allSamps_NBtaxa_OTUsumry$OTUid, otu_table_long_sampOnly_wTaxa$OTUid)
#   ## 1,764 OTUs are dropped with that filtering!
# ## what were their persample/readdepth summaries?
dropped_NBtaxa_OTUsumry <- allSamps_NBtaxa_OTUsumry %>% filter(OTUid %in% NB_droppedOTUs)
#   ## just 29 of the OTUs we dropped are detected in more than 1% of our samples
#   ## manually BLASTing (NCBI webBlast) suggest these are degraded samples, as most % identity scores were well below 90% ID for anything
#   ## the largest OTU we're leaving out? the MYLU OTU!
# #################### sanity check == end ####################
# 
write_csv(otu_table_long_sampOnly_wTaxa, path = "~/github/nhguano/data/text_tables/otu_tables/allTrueSamps_OTUtable_long_wTaxa.csv")
# 
rm(allSamps_VStaxa_boldDB, arthOnly_mergedTaxa, arthOnly_NBtaxa, arthOnly_VStaxa,
   nb_keepers_toAdd, otu_table_long_sampOnly, sampOTUs, vs_keepers_OTUs, arthOnly_mergedTaxa_cleaned,
   NB_droppedOTUs, dropped_NBtaxa_OTUsumry, allSamps_NBtaxa_boldDB, allSamps_NBtaxa_readCounts, allSamps_NBtaxa_OTUsumry)

########################################
## part 5 - final filtering considerations: sampling depth among arth-specific OTUs
########################################
# ## can restart this section by loading two files:
# otu_table_long_sampOnly_wTaxa <- read_csv("https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/allTrueSamps_OTUtable_long_wTaxa.csv.gz")
# allMeta <- read_csv("https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/allbat_meta.csv")


## how many reads per sample?
perSampleOTUssumry <- otu_table_long_sampOnly_wTaxa %>% 
  group_by(SampleID) %>% 
  summarise(nOTUs = n_distinct(OTUid),
            sumReads = sum(OTUreads))

perSampleOTUssumry_wMeta <- merge(perSampleOTUssumry, allMeta, by='SampleID') %>% 
  select(-rBC, -fBC, -SeqBatch, -DNAplate, -DNAwell, -HostSpecies)

### run SRS on a range of sampling depths, then summarize number of OTUs per sample for each sampling depth:
magic_for(print, silent = TRUE)

OTUsumry_atSeqDepth <- for (NUM in seq(500,2500,500)){
  ## get the Samples with at least NUM reads
  tmp_sampleIDs_atSeqDepth <- perSampleOTUssumry %>% filter(sumReads >= NUM) %>% select(SampleID) %>% pull()
  ## reformat long to wide format for SRS sampling
  tmp_otu_table_wide <- otu_table_long_sampOnly_wTaxa %>% 
    filter(SampleID %in% tmp_sampleIDs_atSeqDepth) %>% 
    select(OTUid, OTUreads, SampleID) %>%
    dcast(data = ., formula = OTUid ~ SampleID, value.var='OTUreads', fill = 0)
  row.names(tmp_otu_table_wide) <- tmp_otu_table_wide$OTUid
  tmp_otu_table_wide$OTUid <- NULL
  tmp_otu_table_wide <- as.data.frame(tmp_otu_table_wide)
  ## get Cmin
  Cmin <- min(colSums(tmp_otu_table_wide))
  ## run SRS
  tmp_SRS_output <- SRS(tmp_otu_table_wide, Cmin)
  ## add featureIDs back into the SRSoutput table
  tmp_SRS_output$OTUid <- row.names(tmp_otu_table_wide)
  ## convert back from wide to long:
  tmp_SRS_long <- melt(tmp_SRS_output, id = "OTUid") %>% filter(value != 0)
  colnames(tmp_SRS_long) <- c("OTUid", "SampleID", "Reads")
  ## summarise number of OTUs per sample
  tmp_OTUs_per_sample_sumry <- tmp_SRS_long %>% 
    group_by(SampleID) %>% 
    summarise(nOTUs = n_distinct(OTUid)) %>% 
    mutate(minDepth = NUM)
  print(tmp_OTUs_per_sample_sumry)
}

magic_result_as_dataframe()
OTUsumry_atSeqDepth_df <- dplyr::bind_rows(OTUsumry_atSeqDepth)
rm(OTUsumry_atSeqDepth)

## plot the number of OTUs retained as a histogram, faceted by the 5 sampling depths:
ggplot(OTUsumry_atSeqDepth_df, aes(nOTUs)) + 
  geom_histogram(binwidth = 1, color="black", fill="gray85") +
  facet_wrap(~ minDepth, nrow=5) +
  labs(x = "OTUs detected", y = "number of samples") +
  geom_vline(xintercept = c(5,10,15,20,25), color="blue", linetype="dotted")

## a boxplot of the distribution of OTUs per sample, depending on sampling depths:
ggplot(OTUsumry_atSeqDepth_df, aes(x=minDepth, y=nOTUs, group=minDepth)) +
  #geom_boxplot(outlier.color = NA) + 
  geom_violin(color="black", fill="blue", alpha = 0.25) +
  geom_jitter(aes(x=minDepth, y=nOTUs), width = 25, color="gray15") +
  theme_minimal() +
  scale_x_discrete(limits = c(500, 1000, 1500, 2000, 2500)) +
  labs(x = "minimum seqs per sample", y = "OTUs per sample")
  
## how many distinct OTUs are retained in the global pool of data at each sampling depth?
## what fraction of reads does this global pool reflect?

magic_for(print, silent = TRUE)

global_OTUsumry_atSeqDepth <- for (NUM in seq(500,2500,500)){
  tmp_sampleIDs_atSeqDepth <- perSampleOTUssumry %>% filter(sumReads >= NUM) %>% select(SampleID) %>% pull()
  tmp_all_data_sampleSize <- otu_table_long_sampOnly_wTaxa %>% summarise(n_distinct(SampleID)) %>% pull()
  tmp_select_data_sampleSize <- otu_table_long_sampOnly_wTaxa %>% filter(SampleID %in% tmp_sampleIDs_atSeqDepth) %>% 
    summarise(n_distinct(SampleID)) %>% pull()
  tmp_all_data_globalOTUrichness <- otu_table_long_sampOnly_wTaxa %>%
    summarise(n_distinct(OTUid)) %>% pull()
  tmp_select_data_globalOTUrichness <- otu_table_long_sampOnly_wTaxa %>% filter(SampleID %in% tmp_sampleIDs_atSeqDepth) %>% 
    summarise(n_distinct(OTUid)) %>% pull()
  tmp_all_data_readDepth <- otu_table_long_sampOnly_wTaxa %>% summarise(sumReads = sum(OTUreads)) %>% pull()
  tmp_select_data_readDepth <- otu_table_long_sampOnly_wTaxa %>% 
    filter(SampleID %in% tmp_sampleIDs_atSeqDepth) %>% summarise(sumReads = sum(OTUreads)) %>% pull()
  tmp_global_OTU_richness <- otu_table_long_sampOnly_wTaxa %>%
    filter(SampleID %in% tmp_sampleIDs_atSeqDepth) %>%
    summarise(n_distinct(OTUid)) %>% pull()
  print(data.frame(Depth = NUM, 
                   samples_retained = tmp_select_data_sampleSize,
                   fraction_samples_retained = round((tmp_select_data_sampleSize/tmp_all_data_sampleSize),3),
                   distinct_OTUs_retained = tmp_select_data_globalOTUrichness,
                   fraction_distinct_OTUs_retained = round((tmp_select_data_globalOTUrichness/tmp_all_data_globalOTUrichness),3),
                   reads_retained = round((tmp_select_data_readDepth/tmp_all_data_readDepth),3)
                   ))
}

magic_result_as_dataframe()
global_OTUsumry_atSeqDepth_df <- dplyr::bind_rows(global_OTUsumry_atSeqDepth) %>% distinct()
rm(global_OTUsumry_atSeqDepth)

########################################
## part 6 - what questions can we answer assuming a read depth of 1000 sequences per sample?
########################################
  ## see this script if you want to convert WOY into even number of bins (of days):
  ## https://raw.githubusercontent.com/devonorourke/nhguano/4a07feca887bdfb6e531db4827771c18189f5f45/scripts/r_scripts/taxa_summaries.R

##########
## first question: what kinds of data do we have to look at diversity between years at same sites?
##########

## going to instead just define sampling collection by month of year
perSampleOTUssumry_wMeta$Month <- month(perSampleOTUssumry_wMeta$newDate)
perSampleOTUssumry_wMeta$Year <- year(perSampleOTUssumry_wMeta$newDate)

## regroup samples:
monthSumry <- perSampleOTUssumry_wMeta %>% 
  filter(sumReads >= 1000) %>% 
  group_by(Site, Month, Year) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>% 
  mutate(Year = as.factor(Year),
         Month = as.factor(Month))

monthSumry$Month <- factor(monthSumry$Month, levels = c(4,5,6,7,8,9,10))

ggplot(monthSumry, aes(x=Month, y=nSamples, group=Site, fill=Year)) +
  geom_col(position = position_dodge2()) +
  facet_grid(~Site, scales = "free_x", space = "free_x") +
  theme_minimal() +
  geom_hline(yintercept = c(2,4,6), linetype="dotted", color="navy")
    ## there are 5 sites where we have at least five samples collected for a given month. 
    ## these are all month 7 (July), at BRN, CRN, FOX, HOP, and MAP sites
  

##########
## second question: what kinds of data do we have to look at diversity changes within months for just 2016 samples?
##########
## drop instances where there are less than 5 samples per site month:
filtd_monthSumry <- monthSumry %>% 
  filter(nSamples >= 4)

ggplot(data = filtd_monthSumry %>% filter(Year == 2016), 
       aes(y=Month, x=Site, label=nSamples)) +
  geom_text() +
  theme_minimal()
  ## 7 sites have common sampling months in May/June/July: BRN, CNA, EPS, HOL, HOP, MAP, PEN

## how many samples are retained in each group at a sampling depth of 1000 sequences?
siteYearsumry <- perSampleOTUssumry_wMeta %>% 
  filter(sumReads >= 1000) %>% 
  group_by(StudyID, Site, WOY) %>% summarise(nSamples = n_distinct(SampleID))


########################################
## part 7 - final filtering
########################################
## get a list of just those samples with at least 1000 reads per sample of our arthropod-classified taxa
min1kreads_sampleIDs <- perSampleOTUssumry_wMeta %>% 
  filter(sumReads >= 1000) %>% 
  select(SampleID) %>% pull()

## filter data set with per-OTU read count info among these samples, and save for diversity assessments:
otu_table_long_sampOnly_wTaxa_min1kreads_perSample <- otu_table_long_sampOnly_wTaxa %>% 
  filter(SampleID %in% min1kreads_sampleIDs)

## provide OTU alias based on total abundance of reads across all samples for each OTU ("OTU-1", "OTU-2", etc.)
otu_Alias_df <- otu_table_long_sampOnly_wTaxa_min1kreads_perSample %>%
  group_by(OTUid) %>% 
  summarise(sumOTUreads = sum(OTUreads)) %>%
  arrange(-sumOTUreads) %>%
  mutate(OTUalias = paste0("OTU-", row.names(.))) %>% 
  ungroup %>% 
  select(-sumOTUreads)

otu_table_long_sampOnly_wTaxa_min1kreads_perSample <- merge(otu_table_long_sampOnly_wTaxa_min1kreads_perSample,otu_Alias_df)

## manually edit the single instance of OTU-1, which should be classified as 'Maladera castanea'
otu_table_long_sampOnly_wTaxa_min1kreads_perSample <- 
  otu_table_long_sampOnly_wTaxa_min1kreads_perSample %>% 
  mutate(Genus = case_when(OTUalias == "OTU-1" ~ 'Maladera',
                           TRUE ~ as.character(Genus)),
         Species = case_when(OTUalias == "OTU-1" ~ 'Maladera castanea',
                             TRUE ~ as.character(Species)))

## save the long-form OTU table of our filtered dataset
  ## will import this file in the diversityAnalyses.R code
write_csv(otu_table_long_sampOnly_wTaxa_min1kreads_perSample,
          path = "~/github/nhguano/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa.csv")


## get a list of just those OTUs that passed the filtering
  ## will import this list into QIIME to create the tree used in diversityAnalyses.R code
filtd_OTUid_names <- otu_table_long_sampOnly_wTaxa_min1kreads_perSample %>% 
  distinct(OTUid) %>%
  rename("Feature ID" = OTUid)
write_csv(filtd_OTUid_names,
          "~/github/nhguano/data/text_tables/otu_tables/filtd_sequence_featureIDs.csv")

###################################################################
############################# unused code
######################

## generate rarefaction curves, using SRS
# SRScurve(x = otu_table_wide[1:5], metric = "richness", step = 100, sample = 1000,
#          rarefy.comparison = T, rarefy.comparison.legend = T, ylab = "OTUs",
#          lty = c("solid","longdash"), xlim = c(0,10000))

# 
# ##########
# ## can also reformat the OTU table and check how many OTUs are dropped for a given sampling depth interactively
# ## Use SRS for subsampling
# otu_table_wide_forQIIME <- otu_table_long_sampOnly_wTaxa %>% 
#   filter(SampleID %in% sampls_min500reads) %>%
#   select(OTUid, OTUreads, SampleID) %>%
#   dcast(data = ., formula =  OTUid ~ SampleID, value.var='OTUreads', fill = 0)
# colnames(otu_table_wide_forQIIME)[1] <- "#OTU ID"
# write_delim(otu_table_wide_forQIIME, delim = "\t",
#           path="/Users/devonorourke/Documents/nau_projects/guano/NHguano_redux/redo_taxa_boldANML/otu_table_wide_forQIIME.tsv")
## export that file as QIIME format and view in view.qiime2.org to see how sampling depth changes influence # of OTUs retained in study
  ## 'biom convert -i otu_table_wide_forQIIME.tsv -o otu_table_wide_forQIIME.biom --to-hdf5'
  ## 'qiime tools import --input-path otu_table_wide_forQIIME.biom --output-path otu_table_NHsamplesMin500reads.qza --type 'FeatureTable[Frequency]''
  ## 'qiime diversity alpha-rarefaction --i-table otu_table_NHsamplesMin500reads.qza --p-metrics 'observed_otus' --o-visualization rarecurve_NHsamples_Min500Reads.qzv --p-min-depth 500 --p-max-depth 5000'
      ## the 'rarecurve_NHsamples_Min500Reads.qzv' is then viewed to produce per-sample rarefaction curves like with the SRS code above
  ## 'qiime feature-table summarize --i-table otu_table_NHsamplesMin500reads.qza --o-visualization otu_table_NHsamplesMin500reads_summary.qzv'
      ## the 'otu_table_NHsamplesMin500reads_summary.qzv' is viewed to determint tradeoffs in sampling depth and retained OTUs