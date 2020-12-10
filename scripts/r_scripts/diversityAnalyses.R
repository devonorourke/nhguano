library(tidyverse)
library(SRS)
library(reshape2)
library(lubridate)
library(svglite)
library(ape)
library(phyloseq)
library(btools)
library(Matrix)
library(multcompView)
library(FSA)
library(cowplot)
library(pairwiseAdonis)
library(ggrepel)
library(usedist)
library(ggpubr)
library(indicspecies)

########################################
## part0: importing filtered OTU table containing samples with min 1000 arthropod reads, and merging with metadata
## part1: evaluating overall richness of entire dataset (all sites, both years) for most frequently detected arthropod Orders
## part2: evaluating overall richness of entire dataset for most frequently detected OTUs, and grouping these by Genus-level
## part3: evaluating if alpha diversity and community composition change between years over similar time period
## part4: evaluating community composition

########################################
## part 0 - import filtered sequence data records and metadata, merge fields
########################################

## import file of counts of sequences per OTU per sample, filtered requiring a minimum of 1000 seqs per sample
otu_table_long_filtd <- read_csv("https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa.csv.gz")

## import and merge with metadata
allMeta <- read_csv("https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/allbat_meta.csv")

otu_table_long_filtd_wMeta <- merge(otu_table_long_filtd, allMeta, by="SampleID") %>%
  select(-DNAplate, -DNAwell, -SampleType, -Date, -HostSpecies, -rBC, -fBC, -SeqBatch, -ClassifierMethod)

## import list of sampleIDs where MYLU was identified as species:
## write a list of all SampleIDs classified as MYLU:
NBbat_sampleNames <- read_csv("https://raw.githubusercontent.com/devonorourke/nhguano/master/data/taxonomy/NBclassifier_SampleIDs_wMYLU_assigned.csv")
VSbat_sampleNames <- read_csv("https://raw.githubusercontent.com/devonorourke/nhguano/master/data/taxonomy/VSclassifier_SampleIDs_wMYLU_assigned.csv")
anyClassifier_MYLUassigned_sampleID <- unique(c(NBbat_sampleNames$SampleID, VSbat_sampleNames$SampleID))
otu_table_long_filtd_wMeta <- otu_table_long_filtd_wMeta %>%
  mutate(batTaxa = ifelse(SampleID %in% anyClassifier_MYLUassigned_sampleID, "MYLU", NA))

otu_table_long_filtd_wMeta %>% filter(batTaxa == 'MYLU') %>% summarise(n_distinct(SampleID))
  ## this shows that 249 of our 899 samples here were assigned 'MYLU' as bat taxonomy (so about 1/4 of samples)

## add separate Year and Month metadata labels
otu_table_long_filtd_wMeta <- otu_table_long_filtd_wMeta %>%
  mutate(Year = year(newDate),
         Month = month(newDate)) %>%
  select(-StudyID)

## instead of binning by Month, could also group into N-day intervals (e.x. 45-day windows)
  ### might be able to create N-day windows to regroup samples?
  ## https://raw.githubusercontent.com/devonorourke/nhguano/4a07feca887bdfb6e531db4827771c18189f5f45/scripts/r_scripts/taxa_summaries.R

## generate Day Of Year column
otu_table_long_filtd_wMeta$DOY <- yday(otu_table_long_filtd_wMeta$newDate)
## get range defined with the number of bins (devined in the "by" term of the 'breaks' argument)
window=37         ## 'window' represents number of days to span a given observation)
otu_table_long_filtd_wMeta$DayRange = cut(otu_table_long_filtd_wMeta$DOY, breaks = c(seq(from=1, to=365, by=window)), include.lowest = TRUE)
otu_table_long_filtd_wMeta$Window <- as.numeric(otu_table_long_filtd_wMeta$DayRange)
otu_table_long_filtd_wMeta$RangeStart <- as.character(otu_table_long_filtd_wMeta$DayRange)
otu_table_long_filtd_wMeta$RangeStart <- gsub("\\(", "", otu_table_long_filtd_wMeta$RangeStart)
otu_table_long_filtd_wMeta$RangeStart <- gsub("\\]", "", otu_table_long_filtd_wMeta$RangeStart)
otu_table_long_filtd_wMeta$RangeStart <- gsub("\\[", "", otu_table_long_filtd_wMeta$RangeStart)
otu_table_long_filtd_wMeta <- otu_table_long_filtd_wMeta %>% separate(., col = RangeStart, into=c("RangeStart", "RangeEnd"), sep = ",")
otu_table_long_filtd_wMeta$RangeStart <- as.numeric(otu_table_long_filtd_wMeta$RangeStart)
otu_table_long_filtd_wMeta$RangeEnd <- as.numeric(otu_table_long_filtd_wMeta$RangeEnd)
## convert RangeStart and RangeEnd back to dates
orival <- ymd("2016-01-01")
otu_table_long_filtd_wMeta$WindowStart <- as.Date(otu_table_long_filtd_wMeta$RangeStart, origin = orival)
otu_table_long_filtd_wMeta$WindowEnd <- as.Date(otu_table_long_filtd_wMeta$RangeEnd, origin = orival)
otu_table_long_filtd_wMeta <- otu_table_long_filtd_wMeta %>%
  select(-DayRange, -WindowStart, -WindowEnd)
# 
# ## write to disk with added date info:
write_csv(otu_table_long_filtd_wMeta,
        "~/github/nhguano/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz")
# 

## add a per-Window, per-Year, per-Site matrix of the number of samples collected for supplementary table in report:
nSamples_perSiteWindowYear_matrix <- otu_table_long_filtd_wMeta %>%
  group_by(Site, Year, Window) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>%
  mutate(WindowYear = paste0(Year, "-", Window)) %>%
  select(Site, WindowYear, nSamples) %>%
  ungroup() %>% 
  select(-Year) %>%
  arrange(WindowYear, Site) %>% 
  pivot_wider(values_from = nSamples, names_from = WindowYear, values_fill = 0) %>% 
  arrange(Site)

write_csv(nSamples_perSiteWindowYear_matrix, "~/github/nhguano/supplementaryData/tableS3_nSamples_perSiteWindowYear_matrix.csv")


# ## cleanup
rm(otu_table_long_filtd, NBbat_sampleNames, VSbat_sampleNames, anyClassifier_MYLUassigned_sampleID, orival, window, nSamples_perSiteWindowYear_matrix)
# 
# # ## windows to remember for manual plot labels, if necessary:
# # # Window WindowStart  WindowEnd
# # # 3  2016-03-16 2016-04-22
# # # 4  2016-04-22 2016-05-29
# # # 5  2016-05-29 2016-07-05
# # # 6  2016-07-05 2016-08-11
# # # 7  2016-08-11 2016-09-17
# # # 8  2016-09-17 2016-10-24
# # # 9  2016-10-24 2016-11-30

########################################
## part 1 - richness summaries 
### a. frequency of presence or absence of at least 1 OTU in a given arthropod Order/Family observed
### b. diversity of OTUs per arthropod order... how many OTUs per arth order?
########################################

## start with this file: 
#otu_table_long_filtd_wMeta <- read_csv('https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz')

# ### a.
total_nSamples = n_distinct(otu_table_long_filtd_wMeta$SampleID)

overall_arth_orders_detected_perSample <- otu_table_long_filtd_wMeta %>%
  select(SampleID, Order) %>%
  distinct() %>%
  group_by(Order) %>%
  summarise(nSamples_wOrder = n(),
            fracSamples = nSamples_wOrder/total_nSamples) %>%
  mutate(fracSamples = round(fracSamples, 3)) %>%
  arrange(-nSamples_wOrder)

### b.
overall_OTUs_per_arth_order <- otu_table_long_filtd_wMeta %>%
  group_by(Order) %>%
  summarise(distinct_OTUs = n_distinct(OTUid)) %>%
  arrange(-distinct_OTUs)

## merge into single data table:
overall_arth_order_sumry <- merge(overall_arth_orders_detected_perSample, overall_OTUs_per_arth_order, by="Order") %>%
  arrange(-nSamples_wOrder)

## collapse those orders with less than 10 samples detected into a single group:
rareOrderList <- overall_arth_order_sumry %>%
  filter(nSamples_wOrder < 10) %>%
  select(Order)

mainOrderdf <- overall_arth_order_sumry %>%
  filter(nSamples_wOrder >= 10)

rareOrderFrac <- otu_table_long_filtd_wMeta %>%
  filter(Order %in% rareOrderList$Order) %>% 
  summarise(nSamples_wOrder = n(),
            fracSamples = nSamples_wOrder/total_nSamples) %>%
  mutate(fracSamples = round(fracSamples, 3))

rareOrder_distinctOTUs <- otu_table_long_filtd_wMeta %>%
  filter(Order %in% rareOrderList$Order) %>% 
  summarise(distinct_OTUs = n_distinct(OTUid))

tmp_rare_df <- data.frame(Order = "other taxa", 
                          nSamples_wOrder = " < 10",
                          fracSamples = rareOrderFrac$fracSamples,
                          distinct_OTUs = rareOrder_distinctOTUs$distinct_OTUs)

overall_arth_order_sumry_abbreviated <- rbind(mainOrderdf, tmp_rare_df)

# ## save as Table 1 for manuscript
write_csv(overall_arth_order_sumry_abbreviated,
          path = "~/github/nhguano/tables/table1_overall_arthropodOrder_summaries.csv")

## save also the list of the "other taxa" Order names for legend for Table 1~
write_csv(rareOrderList, path = "~/github/nhguano/tables/table1_otherOrders_names_forLegend.csv")

# ## save as Table 1alt not for manuscript, but in case someone doesn't want these infrequent Orders collapsed together?
write_csv(overall_arth_order_sumry,
          path = "~/github/nhguano/tables/table1alt_overall_arthropodOrder_summaries_noInfrequentOrderCollapsing.csv")


# 
# ## cleanup
rm(overall_arth_orders_detected_perSample, overall_OTUs_per_arth_order, overall_arth_order_sumry, 
   rareOrderList, mainOrderdf, rareOrderFrac, rareOrder_distinctOTUs, tmp_rare_df)

########################################
## part 2 - core features
### a. which OTUs are observed in at least N% (or more) samples?
  ## requiring being detected in at least 5% min... that represents at least 45 of 899 samples...
### b. because multiple OTUs can be assigned the same species, analyzing these potential overlaps makes more sense to aggregate
  ## going to collapse OTUs with shared species (but again, will see that it's not perfect because it biases those not classified to species rank)
### c. instead of collapsing by shared species-labels, and after noting that 49 of 51 'common' OTUs were assigned at least to Genus...
### ... I finally grouped the analysis by shared Genus. The two outliers were counted individually and plotted among the grouped (by Genus) OTUs 

############333 !!!!!!!!!!!!!!!!!
##### in the final plot, note that we're grouping shared Genus ONLY AMONG THOSE CORE FEATURES (OTUs), not across all OTUs assigned to that Genus
############333 !!!!!!!!!!!!!!!!!
# 
# ### ai.
core_OTUs_df <- otu_table_long_filtd_wMeta %>%
  group_by(OTUid, OTUalias, Order, Family, Genus, Species) %>%
  summarise(nSamples_wOTU = n_distinct(SampleID),
            frac_Samples = nSamples_wOTU/total_nSamples) %>%
  filter(frac_Samples >= 0.05) %>% 
  arrange(Order, Genus, Species, -frac_Samples)

core_OTUs_df <- core_OTUs_df[c(3:8,2,1)]
write_csv(core_OTUs_df,
          path = "~/github/nhguano/supplementaryData/tableS5_coreOTU_summary.csv")

  ## note that there are multiple instances of the same species (ex. Phyllophaga hirsuta) having distinct cluster IDs
    ## perhaps not surprising if these OTUs were derived from different populations of these beetles?
  ## going to reformat these to cluster by species labels when available too...

### aii. Same as above, but further clustering if same species labels apply
not_NA_coreOTUs <- core_OTUs_df %>% filter(!is.na(Genus)) %>% filter(!is.na(Species)) ## list of OTUs we'll further cluster by species-level
NA_coreOTUs <- core_OTUs_df %>% filter(is.na(Genus) | is.na(Species)) %>% ungroup() %>% select(-OTUid) ## leave these alone, as NA's are ambiguous (could be diff species!)

tmp1 <- otu_table_long_filtd_wMeta %>%
  filter(OTUid %in% not_NA_coreOTUs$OTUid) %>%
  group_by(Order, Family, Genus, Species) %>%
  summarise(nSamples_wOTU = n_distinct(SampleID),
            frac_Samples = nSamples_wOTU/total_nSamples) %>%
  filter(frac_Samples >= 0.05)

core_OTUs_df_altClustBySpecies <- bind_rows(tmp1, NA_coreOTUs) %>%
  arrange(-nSamples_wOTU)

rm(tmp1, not_NA_coreOTUs, NA_coreOTUs)

## relabeling the empty Family/Genus labels
core_OTUs_df_altClustBySpecies <- core_OTUs_df_altClustBySpecies %>%
  mutate(altGenusLabel = ifelse(is.na(Genus), paste0("f.",Family, " ", OTUalias), Genus),
         altSpeciesLabel = ifelse(grepl("^f.", altGenusLabel), altGenusLabel, Species),
         altSpeciesLabel = ifelse(is.na(altSpeciesLabel), paste0("g.",altGenusLabel," ",OTUalias), Species),
         altSpeciesLabel = ifelse(is.na(altSpeciesLabel), altGenusLabel, altSpeciesLabel))
#   ## now there are only 37 core features that we know might be distinct at species level
#     ## 14 of 37 of these are NA at species-rank (and 1 of those 15 are NA for Genus level)
# 
# ### plotting these?
# ## how many different arthropod orders?
core_OTUs_df_altClustBySpecies %>% ungroup() %>% distinct(Order) %>% arrange(Order)
# 
# ## keep palette consistent with other figure of broad diet pattern (sections 3 and others later):
# ## 'gray25'  == Araneae
# ## 'darkgoldenrod' == Blattodea
# ## "orchid2" == Coleoptera
# ## "tan4"  == Diptera
# ## "#3E6C54" == Ephemeroptera
# ## "#FDAC9F" == Hemiptera
# ## "darkorange" == Hymenoptera
# ## "#808133" == Lepidoptera
# ## "cadetblue3" == Megaloptera
# ## "#e1b580" == Psocodoea
# ## "turquoise4" == Trichoptera
# ## 'gray50' == Trombidiformes
# 
coreOTUpal <- c('darkgoldenrod', "orchid2", "tan4", "#3E6C54", "#808133", "cadetblue3", "turquoise4")
# 
ggplot(core_OTUs_df_altClustBySpecies,
       aes(x = reorder(altSpeciesLabel, -frac_Samples),
           y = frac_Samples,
           fill = Order)) +
  geom_col() +
  scale_fill_manual(values = coreOTUpal) +
  facet_grid(~ Order, space = "free_x", scales = "free_x", shrink = TRUE) +
  theme_minimal() +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = "top") +
  labs(x = "", y="fraction of samples with OTU detected")
# 
# ### alternatively, could just resort to plotting the distinct Genus and combine OTUs within:
# ### using the original core_df data (which was done on a per-OTU level) to regroup according to shared Genus labels
#   ## not including instances where we lack the particular Genus label (2 of 51 OTUs were this case)
tmp_genera_notNA <- otu_table_long_filtd_wMeta %>%
  filter(OTUid %in% core_OTUs_df$OTUid) %>%
  filter(!is.na(Genus)) %>%
  group_by(Order, Family, Genus) %>%
  summarise(nSamples_wTaxa = n_distinct(SampleID),
            fracSamples_wTaxa = nSamples_wTaxa/total_nSamples) %>%
  filter(fracSamples_wTaxa >= 0.05)

## get a list of those 49 OTUs with at least Genus names (to find the 2 without)
tmp_list_genera_notNA <- otu_table_long_filtd_wMeta %>%
  filter(OTUid %in% core_OTUs_df$OTUid) %>%
  filter(!is.na(Genus)) %>%
  distinct(OTUid) %>% pull()

## pull out those two OTUs with unknown Genus labels
tmp_genera_isNA <- otu_table_long_filtd_wMeta %>%
  filter(OTUid %in% core_OTUs_df$OTUid) %>%
  filter(!OTUid %in% tmp_list_genera_notNA) %>%
  group_by(OTUid, OTUalias, Order, Family, Genus) %>%
  summarise(nSamples_wTaxa = n_distinct(SampleID),
            fracSamples_wTaxa = nSamples_wTaxa/total_nSamples) %>%
  filter(fracSamples_wTaxa >= 0.05) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family), Genus),
         Genus = paste0(Genus, " ", OTUalias))

core_genus_df <- bind_rows(tmp_genera_notNA, tmp_genera_isNA) %>%
  arrange(-nSamples_wTaxa)
# 
ggplot(core_genus_df,
       aes(x = reorder(Genus, -fracSamples_wTaxa),
           y = fracSamples_wTaxa,
           fill = Order)) +
  geom_col() +
  scale_fill_manual(values = coreOTUpal) +
  facet_grid(~ Order, space = "free_x", scales = "free_x", shrink = TRUE) +
  theme_classic() +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  #geom_hline(yintercept = 0.1) +
  labs(x = "", y="fraction of samples with taxa detected\n", fill = "Arthropod\norder") +
  guides(fill = guide_legend(nrow = 2))
# 
# ## save this one:
# ##### CAUTION: note that we're grouping shared Genus ONLY AMONG THOSE CORE FEATURES (OTUs), not across all OTUs assigned to that Genus
# 
ggsave("~/github/nhguano/figures/figure2_corefeatures_byGenus.png", height = 12, width = 17, units="cm")
ggsave("~/github/nhguano/figures/figure2_corefeatures_byGenus.pdf", height = 12, width = 17, units="cm")
#ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/corefeatures_byGenus.svg", height = 12, width = 17, units="cm")
# 

write_csv(core_genus_df, 
          file = "~/github/nhguano/data/text_tables/core_genus_detections.csv")

# ## cleanup
rm(core_OTUs_df_altClustBySpecies, core_genus_df,
   tmp_genera_isNA, tmp_genera_notNA, tmp_plot, tmp_list_genera_notNA)
# 
# ## how many of sites have at least 2 (or more) of these core OTUs?
coreOTU_ids <- unique(core_OTUs_df$OTUid)

tmp1_otu_bySite_sumry <- otu_table_long_filtd_wMeta %>%
  filter(OTUid %in% coreOTU_ids) %>%
  group_by(Site) %>%
  summarise(total_nSamples = n_distinct(SampleID))

tmp2_otu_bySite_sumry <- otu_table_long_filtd_wMeta %>%
  filter(OTUid %in% coreOTU_ids) %>%
  group_by(OTUalias, Site) %>%
  summarise(nSamples_wCoreOTU = n()) %>%
  filter(nSamples_wCoreOTU > 1)

core_OTUs_sumry <- merge(tmp2_otu_bySite_sumry, tmp1_otu_bySite_sumry, by="Site") %>%
  mutate(fractionSamples_wCoreOTU = nSamples_wCoreOTU/total_nSamples)
rm(tmp1_otu_bySite_sumry, tmp2_otu_bySite_sumry)

core_OTUs_bySite <- core_OTUs_sumry %>%
  group_by(Site) %>%
  summarise(nCoreOTUs= n_distinct(OTUalias),
            fract_CoreOTUs_perSite = nCoreOTUs / length(coreOTU_ids),
            fract_CoreOTUs_perSite = round(fract_CoreOTUs_perSite, 3)) %>% 
  arrange(fract_CoreOTUs_perSite)
    ## these core OTUs are very common... across the 19 sites surveyed...
      ##... all but one site have at least a quarter of the core OTUs (and it had just 2 samples total!)
      ##... 11 of the 19 sites surveyed have more than 50% of the core OTUs detected in multiple samples

## write the original sumry as matrix of OTU ~ Site, with elements representing frac Samples detected
core_OTUs_sumry_matrix_tmp1 <- core_OTUs_sumry %>%
  mutate(elementVal = (nSamples_wCoreOTU/total_nSamples),
         elementVal = round(elementVal, 2)) %>%
  select(Site, OTUalias, elementVal) %>%
  pivot_wider(names_from = "Site", values_from = "elementVal", values_fill=0)

core_OTUs_sumry_matrix_tmp2 <- core_OTUs_sumry %>%
  mutate(elementVal = paste0("n = ", total_nSamples)) %>%
  distinct(Site, elementVal) %>%
  mutate(OTUalias = "sampleSize") %>%
  pivot_wider(names_from = "Site", values_from = "elementVal") %>% 
  head(1)

core_OTUs_sumry_matrix <- rbind(core_OTUs_sumry_matrix_tmp2, core_OTUs_sumry_matrix_tmp1)
rm(core_OTUs_sumry_matrix_tmp2, core_OTUs_sumry_matrix_tmp1)
write_csv(core_OTUs_sumry_matrix,
          path = "~/github/nhguano/supplementaryData/tableS6_core_OTUs_sumry_matrix.csv")

## a few other summary notes:
## mean and median number of samples retained with at least 1 core OTU detected in a sample:
core_OTUs_sumry %>%
  distinct(Site, total_nSamples) %>%
  summarise(medianSamples = median(total_nSamples),
            meanSamples = mean(total_nSamples))
  ## 42 median, 46.2 mean
## mean and median number of samples across all samples (not just those with 1 or more core OTU)
otu_table_long_filtd_wMeta %>%
  group_by(Site) %>%
  summarise(nSamples = n_distinct(SampleID)) %>%
  summarise(meanSamplesTotal = mean(nSamples),
            medianSamplesTotal = median(nSamples))
  ## 43 median samples, 47.3 mean samples
    #### in other words, most samples have at least one core OTU!

rm(core_OTUs_df, core_OTUs_bySite, core_OTUs_sumry, core_OTUs_sumry_matrix, coreOTU_ids, coreOTUpal,
   tmp1_otu_bySite_sumry, tmp2_otu_bySite_sumry)

########################################
## part 3 - diversity estimates: single site
  ## simplest analysis: how does diversity change at a single site in a single year?
  ## FOX (Fox State Forest) best example with sampling depth across many windows
  ## using 37-day sampling window

## 0a. justify why we're using FOX
## 0b. normalize reads using SRS method
## 0c. create phyloseq object for diversity calculations

## 1a. calculate alpha diversity using 3 metrics (Richness, Shannon's, Faith's PD)
## 1b. compare sampling windows for each metric using Kruskal-Wallis and pairwise Wilcoxon

## 2a. calculate community composition distances using 4 metrics (Dice-Sorensen, Bray-Curtis, Unweighted/Weighted Unifrac)
## 2b. ordinate each distance for:
#  2bii. between-group medians (centroids)
# 2biii. within-group dispersions (via betadisper)
## 2c. ordinate each distance metric with PCoA
## 2d. pairwise adonis and heatmap for comparing differences between particular clusters

## 3. find the top 10 OTUs in each Window (in terms of SRS-normalized relative abundances)

########################################

## start with this file: 
# otu_table_long_filtd_wMeta <- read_csv("~/github/nhguano/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz")

# ##################
# ## part 0a: justification for why we're analyzing windows 3-8 at just FOX state:
otu_table_long_filtd_wMeta %>%
  filter(Year == 2016 & Site == "FOX") %>%
  group_by(Site, Window) %>%
  summarise(nSamples = n_distinct(SampleID)) %>%
  arrange(Window)
  ## drop out Window 9 (just 1 sample)
# 
# ## sample mean/median/sd for those windows at FOX in 2016
otu_table_long_filtd_wMeta %>%
  filter(Year == 2016 & Site == "FOX") %>%
  group_by(Site, Window) %>%
  summarise(nSamples = n_distinct(SampleID)) %>%
  filter(Window != 9) %>%
  arrange(Window) %>%
  summarise(mean_nSamples = mean(nSamples),
            median_nSamples = median(nSamples),
            sd_nSamples = sd(nSamples))
#   ## mean number of samples = 13.5, median 9.5, sd = 9.69
# 
# ##################
# ## part 0b: rarefy table using SRS for these FOX samples (windows 3-8, year 2016 only, FOX site only)
# 
# ## get metadata for these samples and then generate OTU table for alpha/adonis/betadisp/ordinations:
fox2016_meta <- otu_table_long_filtd_wMeta %>%
  mutate(Window = as.factor(Window)) %>%
  filter(Year == 2016 & Site == "FOX" & Window != 9) %>%
  select(SampleID, Window, newDate) %>%
  distinct()
#   
# ## select these FOX samples and generate SRS subsampled OTU table
# ## for SRS paper see: https://peerj.com/articles/9593/
fox_OTUtable_wide_raw <- as.data.frame(otu_table_long_filtd_wMeta %>%
  filter(Year == 2016 & Site == "FOX" & Window != 9) %>%
  select(SampleID, OTUreads, OTUid) %>%
  pivot_wider(values_from = "OTUreads", names_from = "SampleID", values_fill = 0))
row.names(fox_OTUtable_wide_raw) <- fox_OTUtable_wide_raw$OTUid
fox_OTUtable_wide_raw$OTUid <- NULL
# ## get Cmin for this subset
fox2016_Cmin <- min(colSums(fox_OTUtable_wide_raw))
# ## run SRS
fox_OTUtable_wide_norm <- SRS(fox_OTUtable_wide_raw, fox2016_Cmin, seed = 1)
# ## add featureIDs back into the SRSoutput table
row.names(fox_OTUtable_wide_norm) <- row.names(fox_OTUtable_wide_raw)
rm(fox_OTUtable_wide_raw, fox2016_Cmin)

## create binary version of same object
# fox_OTUtable_wide_norm_binary <- fox_OTUtable_wide_norm
# fox_OTUtable_wide_norm_binary[fox_OTUtable_wide_norm_binary > 0] <- 1

# 
# ##################
# ## part 0c: diversity estimates using SRS output calculated with Phyloseq to incorporate tree info for both alpha and distance calcs
# ## import tree object to for phylogenetic distance measures
# ## download qza file:
# # download.file("https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/trees/finalfiltd_rootedtree.nwk", "tree.nwk")
# # tree <- read.tree(file = "tree.nwk")
# 
# ## alternatively run from local:
tree <- read.tree(file = "~/github/nhguano/data/text_tables/trees/finalfiltd_rootedtree.nwk")
# 
# ## create phyloseq object
# ## import OTU data in phyloseq format
fox2016_phyOTU <- otu_table(fox_OTUtable_wide_norm, taxa_are_rows = TRUE)
# fox2016_phyOTU_binary <- otu_table(fox_OTUtable_wide_norm_binary, taxa_are_rows = TRUE)

## import taxonomy data in phyloseq format:
fox2016_taxdf <- as.data.frame(otu_table_long_filtd_wMeta %>% distinct(OTUid, Kingdom, Phylum, Class, Order, Family, Genus, Species))
row.names(fox2016_taxdf) <- fox2016_taxdf$OTUid
fox2016_taxdf$OTUid <- NULL
fox2016_taxmat <- as.matrix(fox2016_taxdf)
fox2016_phyTAX <- tax_table(fox2016_taxmat)
rm(fox2016_taxmat, fox2016_taxdf)

## get metadata:
fox2016_phyMETA <- as.data.frame(fox2016_meta)
row.names(fox2016_phyMETA) <- fox2016_phyMETA$SampleID
fox2016_phyMETA <- sample_data(fox2016_phyMETA)

# ## bundle OTU table, taxonomy, metadata, and tree info into single phyloseq objects
## for abundance-based OTU table
fox2016_phydat <- phyloseq(fox2016_phyOTU, fox2016_phyTAX)
fox2016_phydat <- merge_phyloseq(fox2016_phydat, tree, fox2016_phyMETA)
fox2016_phydat
## for binary OTU table
# fox2016_phydat_binary <- phyloseq(fox2016_phyOTU_binary, fox2016_phyTAX)
# fox2016_phydat_binary <- merge_phyloseq(fox2016_phydat_binary, tree, fox2016_phyMETA)
# fox2016_phydat_binary

rm(fox2016_phyOTU, fox2016_phyOTU_binary, fox2016_phyTAX)
# 
# ##################
# ## part 1a: alpha diversity estimates: observed OTUs, Shannon's, Faith's PD
# ## need to estimate Shannon's separate from Faith's/Observed
# ## get Observed and Faith's first
fox2016_faith_obs <- estimate_pd(fox2016_phydat)
# ## now get Shannon's entropy
fox2016_shannons <- estimate_richness(fox2016_phydat, measures = "Shannon") %>% rename("H" = Shannon)
# ## merge, then combine with metadata
fox2016_alpha_df <- data.frame(fox2016_faith_obs, fox2016_shannons) %>%
  mutate(SampleID = row.names(.))
rm(fox2016_shannons, fox2016_faith_obs)
fox2016_alpha_df <- merge(fox2016_alpha_df, fox2016_meta)
# ## long format for plotting ease and stat calcs:
fox2016_alpha_df_long <- fox2016_alpha_df %>%
  pivot_longer(cols = c("SR", "PD", "H"), values_to = "Alpha_value", names_to = "Metric")
# 
# 
# ##################
# ## part 1b: running KW tests for each alpha metric:
# ## why not just run ANOVA? Nothing is parametric about these data...
  shapiro.test(fox2016_alpha_df_long %>% filter(Metric == "SR") %>% select(Alpha_value) %>% pull())
#     ## W = 0.90191, p-value = 1.313e-05
  shapiro.test(fox2016_alpha_df_long %>% filter(Metric == "H") %>% select(Alpha_value) %>% pull())
#     ## W = 0.88998, p-value = 4.188e-06
  shapiro.test(fox2016_alpha_df_long %>% filter(Metric == "PD") %>% select(Alpha_value) %>% pull())
#     ## W = 0.93816, p-value = 0.0007057
# 
# ## KW tests:
kw_input_fox2016_SR <- fox2016_alpha_df_long %>% filter(Metric == "SR") %>% mutate(Window = as.character(Window)) %>% select(Window, Alpha_value)
kw_result_fox2016_SR <- kruskal.test(kw_input_fox2016_SR$Alpha_value ~ kw_input_fox2016_SR$Window)
kw_result_fox2016_SR
#   ## significant, but barely:
#   ## chi-squared = 13.351, df = 5, p-value = 0.0203
kw_input_fox2016_H <- fox2016_alpha_df_long %>% filter(Metric == "H") %>% mutate(Window = as.character(Window)) %>% select(Window, Alpha_value)
kw_result_fox2016_H <- kruskal.test(kw_input_fox2016_H$Alpha_value ~ kw_input_fox2016_H$Window)
#   ## not significant
#   ## chi-squared = 9.7861, df = 5, p-value = 0.08153
kw_input_fox2016_PD <- fox2016_alpha_df_long %>% filter(Metric == "PD") %>% mutate(Window = as.character(Window)) %>% select(Window, Alpha_value)
kw_result_fox2016_PD <- kruskal.test(kw_input_fox2016_PD$Alpha_value ~ kw_input_fox2016_PD$Window)
#   ## significant, but barely:
#   ## Kruskal-Wallis chi-squared = 11.989, df = 5, p-value = 0.03494

# ## cleanup
rm(kw_result_fox2016_H, kw_result_fox2016_PD, kw_result_fox2016_SR)
# 
# ### overview: might be group differences, but we follow up with pairwise test to identify which are significantly different:
# ### post hoc test using pairwise Dunn test:
# ### generate table of letters (for plots) and actual pairwise-pvals (For supplementary tables)
# ### letters first:
dunns_letters_fox2016 <- function(data, metric){
  tmp_dunn_list <- dunnTest(Alpha_value ~ as.factor(Window), data=data, method="bh")
  tmp_dunn_df <- data.frame(tmp_dunn_list[2]) %>% select(res.Comparison, res.P.unadj, res.P.adj) %>% 
    mutate(res.Comparison = str_replace_all(string = res.Comparison, pattern = " ", replacement = "")) %>% 
    arrange(res.Comparison)
  tmp_mat_un <- matrix(NA, nrow=6, ncol=6)
  tmp_mat_un[lower.tri(tmp_mat_un)] <- tmp_dunn_df$res.P.unadj
  row.names(tmp_mat_un) <- paste0("window", seq(3,8))
  colnames(tmp_mat_un) <- paste0("window", seq(3,8))
  tmp_mat_un <- Matrix::forceSymmetric(tmp_mat_un, uplo="L")
  tmp_lmat_un <- multcompLetters(tmp_mat_un, compare="<=", threshold=0.05, Letters=letters)
  tmp_labels_un <- data.frame(tmp_lmat_un$Letters) %>%
    rename(Letters = tmp_lmat_un.Letters) %>%
    mutate(Window = row.names(.),
           Window = str_remove(Window, "window"),
           Metric = metric,
           pType = "unadj") 
  tmp_mat_adj <- matrix(NA, nrow=6, ncol=6)
  tmp_mat_adj[lower.tri(tmp_mat_adj)] <- tmp_dunn_df$res.P.adj
  row.names(tmp_mat_adj) <- paste0("window", seq(3,8))
  colnames(tmp_mat_adj) <- paste0("window", seq(3,8))
  tmp_mat_adj <- Matrix::forceSymmetric(tmp_mat_adj, uplo="L")
  tmp_lmat_adj <- multcompLetters(tmp_mat_adj, compare="<=", threshold=0.05, Letters=letters)
  tmp_labels_adj <- data.frame(tmp_lmat_adj$Letters) %>%
    rename(Letters = tmp_lmat_adj.Letters) %>%
    mutate(Window = row.names(.),
           Window = str_remove(Window, "window"),
           Metric = metric,
           pType = "adj")
  rbind(tmp_labels_un, tmp_labels_adj)
}

fox2016_labels_SR <- dunns_letters_fox2016(kw_input_fox2016_SR, "SR")
fox2016_labels_H <- dunns_letters_fox2016(kw_input_fox2016_H, "H")
fox2016_labels_PD <- dunns_letters_fox2016(kw_input_fox2016_PD, "PD")
fox2016_LetterLabels_all <- rbind(fox2016_labels_SR, fox2016_labels_H, fox2016_labels_PD)
rm(fox2016_labels_SR, fox2016_labels_H, fox2016_labels_PD)

## we find that for adjusted p-values, no pairwise significant differences reported
## for p values NOT corrected for multiple tests, we see few differences still, generally between window8 and others
 
# ## use a modified function from above to retain the actual p-values for supplementary file:
dunns_pvals_fox2016 <- function(data, metric){
  tmp_dunn_list <- dunnTest(Alpha_value ~ as.factor(Window), data=data, method="bh")
  tmp_dunn_df <- data.frame(tmp_dunn_list[2]) %>% select(res.Comparison, res.P.unadj, res.P.adj) %>% 
    mutate(res.Comparison = str_replace_all(string = res.Comparison, pattern = " ", replacement = "")) %>% 
    arrange(res.Comparison) %>% 
    separate(col=res.Comparison, into = c("windowA", "windowB"), sep="-") %>% 
    mutate(Metric = metric,
           res.P.unadj = round(res.P.unadj, 3),
           res.P.adj = round(res.P.adj, 3))
  tmp_dunn_df
}
# 
fox2016_dunn_SR <- dunns_pvals_fox2016(kw_input_fox2016_SR, "SR")
fox2016_dunn_H <- dunns_pvals_fox2016(kw_input_fox2016_H, "H")
fox2016_dunn_PD <- dunns_pvals_fox2016(kw_input_fox2016_PD, "PD")
# 
fox2016_dunn_all <- rbind(fox2016_dunn_SR, fox2016_dunn_H, fox2016_dunn_PD)
rm(fox2016_dunn_SR, fox2016_dunn_H, fox2016_dunn_PD)

# write_csv(fox2016_dunn_all, path = "~/github/nhguano/data/text_tables/fox2016_dunn_allMetrics.csv")

#   ## no significant pairwise differences after correcting for multiple tests

# 
# ## plot to visualize these diversity metrics, grouped by sampling Window, faceted by Alpha metric
# ## set levels to order alpha metric facets
fox2016_alpha_df_long$Metric <- factor(fox2016_alpha_df_long$Metric, levels = c("SR", "H", "PD"))
# 
# ## and plot
fox2016_alpha_plotdat <- merge(fox2016_alpha_df_long, fox2016_LetterLabels_all, by=c("Window", "Metric"))

p_foxalpha_sr <- ggplot(fox2016_alpha_plotdat %>% filter(Metric == "SR" & pType == "adj"), aes(x=Window, y=Alpha_value, group=Window)) +
  geom_boxplot(outlier.color = NA, color = "black") +
  geom_jitter(alpha=0.55, width = 0.1, size=1.5) +
  labs(x="", y='Observed OTUs') +
  theme_classic() +
  theme(axis.text.x = element_blank())

p_foxalpha_h <- ggplot(fox2016_alpha_plotdat %>% filter(Metric == "H" & pType == "adj"), aes(x=Window, y=Alpha_value, group=Window)) +
  geom_boxplot(outlier.color = NA, color = "black") +
  geom_jitter(alpha=0.55, width = 0.1, size=1.5) +
  labs(x="", y='Shannon\'s H') +
  theme_classic() +
  theme(axis.text.x = element_blank())

p_foxalpha_pd <- ggplot(fox2016_alpha_plotdat %>% filter(Metric == "PD" & pType == "adj"), aes(x=Window, y=Alpha_value, group=Window)) +
  geom_boxplot(outlier.color = NA, color = "black") +
  geom_jitter(alpha=0.55, width = 0.1, size=1.5) +
  theme_classic() +
  labs(x="", y='Faith\'s PD')

ggarrange(p_foxalpha_sr, p_foxalpha_h, p_foxalpha_pd,
          labels = c("A", "B", "C"), ncol=1, nrow=3)

rm(p_foxalpha_sr, p_foxalpha_h, p_foxalpha_pd)
# 
# ##save
# ggsave("~/github/nhguano/supplementaryData/figureS1_fox2016_alpha_boxplots.png", height = 15, width = 10, units = 'cm', dpi=150)
# ggsave("~/github/nhguano/supplementaryData/figureS1_fox2016_alpha_boxplots.pdf", height = 15, width = 10, units = 'cm', dpi=300)
# 

## plot a heatmap of pvalues (adjusted and non-adjusted) for the pairwise comparisons of these Windows
fox2016_dunn_all$Metric <- factor(fox2016_dunn_all$Metric, levels = c("SR", "H", "PD"))
p_heatmap_fox2016_unadj <- 
  ggplot(fox2016_dunn_all, 
         aes(x=windowA, y=windowB, fill=res.P.unadj, label=res.P.unadj)) + 
    geom_tile(color="black") + 
    geom_text(size = 2.5) +
    facet_wrap(~Metric) + 
    coord_fixed() +
    scico::scale_fill_scico(palette = "bilbao", end = 0.65) +
    labs(x="\nsampling window", y="sampling window\n", fill="p.value") +
    theme_classic()

p_heatmap_fox2016_adj <- 
  ggplot(fox2016_dunn_all, 
         aes(x=windowA, y=windowB, fill=res.P.adj, label=res.P.adj)) + 
    geom_tile(color="black") + 
    geom_text(size = 2.5) +
    facet_wrap(~Metric) + 
    coord_fixed() +
    scico::scale_fill_scico(palette = "bilbao", end = 0.65) +
    labs(x="\nsampling window", y="sampling window\n", fill="p.value") +
    theme_classic()

ggarrange(p_heatmap_fox2016_adj, p_heatmap_fox2016_unadj,
          labels = c("A", "B"), nrow=2, ncol=1, common.legend = TRUE)

## save as supplementary figure 2
# ggsave("~/github/nhguano/supplementaryData/figureS2_fox2016_alpha_heatmap_pvals.png", height = 15, width = 23, units = 'cm', dpi=150)
# ggsave("~/github/nhguano/supplementaryData/figureS2_fox2016_alpha_heatmap_pvals.pdf", height = 15, width = 23, units = 'cm', dpi=300)

# ## alpha diversity cleanup:
rm(fox2016_alpha_df, fox2016_alpha_df_long,
   kw_input_fox2016_H, kw_input_fox2016_PD, kw_input_fox2016_SR, fox2016_LetterLabels_all, fox2016_dunn_all, 
   fox2016_alpha_plotdat, p_heatmap_fox2016_adj, p_heatmap_fox2016_unadj, dunns_letters_fox2016, dunns_pvals_fox2016)
# 
# 
# ##################
# ## part 2a: community composition
# ## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method
fox2016_dist_ds <- phyloseq::distance(fox2016_phydat, "bray", binary = TRUE)
fox2016_dist_uu <- phyloseq::distance(fox2016_phydat, "unifrac", weighted=FALSE)
# 
# ##################
# ## part 2b: community composition
# ## 2bi. run PERMANOVA (via Adonis) for each distance method to evaluate if samples within each Window (group) have different centroids
fox2016_adonis_function <- function(distanceData, distanceMetric){
  adonis_tmp <- adonis2(distanceData ~ Window, data = fox2016_meta, )
  data.frame(adonis_tmp) %>% mutate(Metric = distanceMetric,
                                            Class = row.names(.))
}
# 
fox2016_adonis_ds <- fox2016_adonis_function(fox2016_dist_ds, "Dice-Sorensen")
fox2016_adonis_uu <- fox2016_adonis_function(fox2016_dist_uu, "Unifrac Unweighted")
# 
fox2016_adonis_all <- rbind(fox2016_adonis_ds, fox2016_adonis_uu)
rm(fox2016_adonis_ds, fox2016_adonis_uu)
#   ## Window is significant main effect regardless of distance metric
# 
fox2016_adonis_all <- fox2016_adonis_all %>% 
 mutate(R2 = round(R2, 3),
        Sum.Sq = round(SumOfSqs, 3),
        F.value = round(F, 3)) %>% 
 select(-SumOfSqs, -F) %>% 
  relocate(Class, Df, Sum.Sq, F.value, R2, Pr..F., Metric)
write.csv(fox2016_adonis_all, quote=FALSE, row.names = FALSE,
          file = "~/github/nhguano/supplementaryData/tableS7_fox2016_adonis_allMetrics.csv")
# 
# ## part 2bii. Test for homogeneity of dispersion
fox2016_betadisper_function <- function(distanceData, distanceMetric){
 tmp_disper_list <- betadisper(d = distanceData, group =  fox2016_meta$Window, type = c("median"))
 tmp_disper_anova <- data.frame(anova(tmp_disper_list))
 tmp_disper_anova <- tmp_disper_anova %>% 
   mutate(Class = row.names(.)) %>% mutate(Metric = distanceMetric) %>% mutate(Factor = "Window")
   tmp_disper_anova[,c(6,1,2,3,4,5,8,7)]
 }
# 
fox2016_bdisp_ds <- fox2016_betadisper_function(fox2016_dist_ds, "Dice-Sorensen")
fox2016_bdisp_ds
#   ## non significant at p <= 0.05; p = 0.165 (i.e. we can't reject null that groups have same dispersions... good!)
fox2016_bdisp_uu <- fox2016_betadisper_function(fox2016_dist_uu, "Unifrac Unweighted")
fox2016_bdisp_uu
#   ## non significant but closer to p<0.05 threshold; p = 0.075
# 
fox2016_bdisp_all <- rbind(fox2016_bdisp_ds, fox2016_bdisp_uu)
rm(fox2016_bdisp_ds, fox2016_bdisp_uu)
fox2016_bdisp_all <- fox2016_bdisp_all %>% 
   mutate(Sum.Sq = round(Sum.Sq, 3),
          Mean.Sq = round(Mean.Sq, 3),
          F.value = round(F.value, 3),
          Pr..F. = round(Pr..F., 3)) %>% 
   select(-Factor)
write_csv(fox2016_bdisp_all,
          path = "~/github/nhguano/supplementaryData/tableS8_fox2016_betadisper_allMetrics.csv")
# 
# ##################
# ## part 2c: ordination of distance data sets
# ### ordinations:
fox2016_ordi_function <- function(distanceValues, distanceMetric){
  tmp_pcoa <- ordinate(fox2016_phydat, method="PCoA", distance = distanceValues)
  tmp_pcoa_list <- plot_ordination(fox2016_phydat, tmp_pcoa)
  tmp_pcoa_df <- tmp_pcoa_list$data %>%
    mutate(SampleID = row.names(.),
           Axis1lab = tmp_pcoa_list$labels$x,
           Axis2lab = tmp_pcoa_list$labels$y,
           Metric = distanceMetric)
  merge(tmp_pcoa_df, fox2016_meta) %>% mutate(Window = as.character(Window))
}
# 
fox2016_orddata_ds <- fox2016_ordi_function(fox2016_dist_ds, "Dice-Sorensen")
fox2016_orddata_uu <- fox2016_ordi_function(fox2016_dist_uu, "Unweighted-Unifrac")
# 
fox2016_orddata_all <- rbind(fox2016_orddata_ds, fox2016_orddata_uu)
rm(fox2016_orddata_ds, fox2016_orddata_uu)
# 
# ## plot with metadata
# #### generate palette that works with color blind data
# scico::scico(length(unique(fox2016_pcoa_df$Window)), palette = "batlow")
#   # original: "#001959" "#184E60" "#577646" "#B28C32" "#FCA68C" "#F9CCF9"
# # going to modify blue and green values to make easier to contrast Windows 6, 7, 8
fox2016_pal <- rev(c("#001447", "#267b97", "#557d3f", "#B28C32", "#FCA68C", "#F9CCF9"))
# 
# ## plot individually, then stitch together
fox2016_ordplot_function <- function(inputdata){
  ggplot(data = inputdata,
         aes(x=Axis.1, y=Axis.2)) +
    geom_point(aes(color = Window), size = 3) +
    stat_ellipse(aes(group = Window, color = Window), alpha = 0.5) +
    labs(x = unique(inputdata$Axis1lab),
         y = unique(inputdata$Axis2lab)) +
    theme_classic() +
    scale_color_manual(values = fox2016_pal) +
    coord_fixed() +
    guides(color=guide_legend(nrow=1,byrow=TRUE))
}
# 
p_fox2016_ord_ds <- fox2016_ordplot_function(fox2016_orddata_all %>% filter(Metric == "Dice-Sorensen"))
p_fox2016_ord_uu <- fox2016_ordplot_function(fox2016_orddata_all %>% filter(Metric == "Unweighted-Unifrac"))
 
ggarrange(p_fox2016_ord_ds, p_fox2016_ord_uu,
          labels = c("A", "B"),
          nrow = 1, ncol = 2,
          align = "hv", common.legend = TRUE)

ggsave("~/github/nhguano/figures/figure3a_fox2016_ordinations_all.png", width = 20, height = 10, units = 'cm', dpi=150)
ggsave("~/github/nhguano/figures/figure3a_fox2016_ordinations_all.pdf", width = 20, height = 10, units = 'cm', dpi=300)
  ## will end up including this image with other plots below to make a multifaceted plot 

rm(p_fox2016_ord_ds, p_fox2016_ord_uu)
# 
# ##################
# ## part 2d: pairwise Adonis to determine what, if any, community comopsition groups are different
# ## generate long form of all values for supplementary table
fox2016_pairwise_adonis_func <- function(distancevals, metric){
  pairwise.adonis(distancevals, fox2016_meta$Window, p.adjust.m = "BH") %>%
    mutate(Metric = metric) %>%
    separate(col = pairs, into = c("windowA", "windowB"), sep = "vs") %>%
    arrange(windowA, windowB) %>% 
    rename(Sum.Sq = SumsOfSqs, F.value = F.Model, BHadjusted.pvalue = "p.adjusted")
}
# 
fox2016_pwadonis_ds <- fox2016_pairwise_adonis_func(fox2016_dist_ds, "Dice-Sorensen")
fox2016_pwadonis_uu <- fox2016_pairwise_adonis_func(fox2016_dist_uu, "Unweighted-Unifrac")
fox2016_pwadonis_all <- rbind(fox2016_pwadonis_ds, fox2016_pwadonis_uu)
rm(fox2016_pwadonis_ds, fox2016_pwadonis_uu)

write_csv(fox2016_pwadonis_all, 
          "~/github/nhguano/data/text_tables/fox2016_pwAdonis_allmetrics.csv")


# 
# ## cleanup
rm(fox2016_adonis_all, fox2016_bdisp_all, fox2016_orddata_all, fox2016_pwadonis_all)
rm(list = ls(pattern = "fox2016_.*._func*"))
rm(list = ls(pattern = "p_fox2016.*"))

# 
# ##################
# ## part 3: we find that early/late windows are different in terms of community composition... so what are those particular taxa?
# 
# ## pivot wide to long for SRS output, merge with taxonomy and metadata info
fox_OTUtable_wide_norm$OTUid <- row.names(fox_OTUtable_wide_norm)
fox_OTUtable_long_norm <- fox_OTUtable_wide_norm %>%
  pivot_longer(-OTUid, names_to = "SampleID", values_to = "SRSreads") %>%
  filter(SRSreads > 0)
fox_taxa <- otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% fox_OTUtable_long_norm$SampleID) %>%
  distinct(OTUid, OTUalias, Class, Order, Family, Genus, Species)
fox_OTUtable_long_norm_wTaxa <- merge(fox_OTUtable_long_norm, fox_taxa)
fox_OTUtable_long_norm_wTaxa_wMeta <- merge(fox_OTUtable_long_norm_wTaxa, fox2016_meta)
rm(fox_OTUtable_wide_norm, fox_OTUtable_long_norm, fox_OTUtable_long_norm_wTaxa)
# 
# ## calculate fraction of reads at Order and Genus-levels (separately)... but will focus here strictly on # detections

###### first plot is at Order level...
## for all data, across all sampling windows, how many detections total?
fox2016_nGlobalDetections <- fox_OTUtable_long_norm_wTaxa_wMeta %>% 
  summarise(nDetections_perOrder = n()) %>% pull()

## for all data, across all sampling windows, what is the fraction of detections per Order?
## retain only those Orders detected in at least 1% of all samples
fox2016_topOrders <- fox_OTUtable_long_norm_wTaxa_wMeta %>% 
  group_by(Class, Order) %>% 
  summarise(nDetections_perOrder = n()) %>% 
  mutate(pDetections_perOrder = nDetections_perOrder / fox2016_nGlobalDetections) %>% 
  filter(pDetections_perOrder > 0.01)

## for all data, how many detections occur per Order across all samples in a given sampling window?
## note we're selecting just those top Orders from above
fox2016_nDetections_perOrder_perWindow <- fox_OTUtable_long_norm_wTaxa_wMeta %>% 
  filter(Order %in% fox2016_topOrders$Order) %>% 
  group_by(Window, Class, Order) %>% 
  summarise(nDetections_perOrder = n())

## summarise the total number of detections occurring in a sampling window across all samples and all taxa
fox2016_nDetections_allOrders_perWindow <- fox2016_nDetections_perOrder_perWindow %>% 
  filter(Order %in% fox2016_topOrders$Order) %>% 
  group_by(Window) %>% 
  summarise(nDetections_allOrders = sum(nDetections_perOrder))

## get the fraction of detections each order is represented as per window
fox2016_orderDetectionsPerWindow <- merge(fox2016_nDetections_perOrder_perWindow, fox2016_nDetections_allOrders_perWindow, by = "Window") %>%
  mutate(pDetections_perOrder = nDetections_perOrder / nDetections_allOrders,
         infrequentCase = ifelse(pDetections_perOrder < 0.01, TRUE, FALSE)) ## sanity check about filtering above

## plot
# ## keep palette consistent with other figure of broad diet pattern:
# ## 'gray25'  == Araneae; 'darkgoldenrod' == Blattodea; "orchid2" == Coleoptera; "tan4"  == Diptera
# ## "#3E6C54" == Ephemeroptera; "#FDAC9F" == Hemiptera; "darkorange" == Hymenoptera; "#808133" == Lepidoptera
# ## "cadetblue3" == Megaloptera; "#e1b580" == Psocodoea; "turquoise4" == Trichoptera; 'gray50' == Trombidiformes

#10 Orders included are: "Araneae" "Coleoptera" "Diptera" "Ephemeroptera" "Hemiptera" "Hymenoptera" "Lepidoptera" "Psocodoea" "Trichoptera" "Trombidiformes"
fox2016_orderPal <- c("gray25", "orchid2", "tan4", "#3E6C54", "#FDAC9F", 'darkorange', "#808133", "#e1b580", "turquoise4", 'gray50')

## retain same color scheme as earlier
p_3b <- ggplot(data = fox2016_orderDetectionsPerWindow,
       aes(x=Window, y=pDetections_perOrder, fill=Order)) +
  geom_col() +
  scale_fill_manual(values=fox2016_orderPal) +
  labs(x="\nWindow", y="fraction of detections per window") +
  theme_classic()
p_3b

ggsave("~/github/nhguano/figures/figure3b_fox2016_detectionsPerWindowPerOrder.png", height=5, width=5, dpi=150)
ggsave("~/github/nhguano/figures/figure3b_fox2016_detectionsPerWindowPerOrder.pdf", height=5, width=5, dpi=300)
  ## will include this with ordination plot earlier, as well as genus-level below below

###### second plot is at Genus level...
### a bit more refined: summarise the abundant taxa in each sampling window by arthropod Genus
### still analyzing in terms of detections (not read counts), but we do filter using reads and detections
fox2016_sumry_PerWindow <- fox_OTUtable_long_norm_wTaxa_wMeta %>%
  group_by(Window) %>% 
  summarise(nSamples_perWindow = n_distinct(SampleID),
            nReads_perWindow_global = sum(SRSreads))

fox2016_topGenus_perWindow <- fox_OTUtable_long_norm_wTaxa_wMeta %>%
  mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>%
  group_by(Window, Order, Genus) %>%
  summarise(nReads_perWindow_perGenus = sum(SRSreads),
            nSamples_perWindow_perGenus = n_distinct(SampleID)) %>%
  merge(., fox2016_sumry_PerWindow) %>%
  mutate(pReads_perWindow_perGenus = nReads_perWindow_perGenus / nReads_perWindow_global,
         pSamples_perWindow_perGenus = nSamples_perWindow_perGenus / nSamples_perWindow) %>%
  mutate(pReads_perWindow_perGenus = round(pReads_perWindow_perGenus, 3),
         pSamples_perWindow_perGenus = round(pSamples_perWindow_perGenus, 3)) %>%
  mutate(Window = as.numeric(as.character(Window))) %>%
  filter(pSamples_perWindow_perGenus >= 0.2 & pReads_perWindow_perGenus >= 0.005) %>%
  select(Window, Order, Genus, pSamples_perWindow_perGenus)

## which Orders are include in the plot?
fox2016_topGenus_perWindow %>% arrange(Order) %>% distinct(Order) %>% pull()
## 7 Orders represented... keep same color palette per Order
#"Coleoptera" "Diptera" "Ephemeroptera" "Hymenoptera" "Lepidoptera" "Psocodoea" "Trichoptera"
fox2016_genusPal <- c("orchid2", "tan4", "#3E6C54", 'darkorange', "#808133", "#e1b580", "turquoise4")

## plot
p_3c <- ggplot() +
  geom_point(data=fox2016_topGenus_perWindow,
             aes(x=Window, y=pSamples_perWindow_perGenus),
             show.legend = TRUE) +
  geom_label_repel(data=fox2016_topGenus_perWindow,
                  aes(x=Window, y=pSamples_perWindow_perGenus, label=Genus, 
                      fill=Order, fontface = "bold",
                      color = ifelse(Order == "Psocodea", "gray90", "gray20")),
                  direction = "y", nudge_x = -0.05,
                  hjust = 1,
                  segment.size = 0.2,
                  segment.colour = "black",
                  min.segment.length = 0,
                  size = 4) +
  scale_fill_manual(values = fox2016_genusPal) +
  scale_color_manual(values = c("white", "black")) +
  scale_x_continuous(limits = c(2.6,8), breaks = seq(3,8), labels = seq(3,8 )) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(0,1)) +
  labs(x = "\n Window", y="fraction samples with taxa", color = "Order") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  guides(color = FALSE, 
         fill = guide_legend(override.aes = aes(color = NA)))

p_3c
ggsave("~/github/nhguano/figures/figure3c_fox2016_detectionsPerWindowPerGenus.png", height=9, width=15, dpi=150)
ggsave("~/github/nhguano/figures/figure3c_fox2016_detectionsPerWindowPerGenus.pdf", height=9, width=15, dpi=300)


ggarrange(p_3b, p_3c,
          nrow=1, ncol=2)

## end up stitching together the ordination (3a), order-level (3b), and genus-level (3c) observations into single plot
  ## in Illustrator...

# ## cleanup:
rm(list = ls(pattern = "^fox*"))
rm(p_3b, p_3c)

########################################
## part 4 - diversity estimates: multipe sites, same year
## how does diversity change between sites and sampling windows in a single year?
## using same 37-day sampling window

## 0a. justify why we're using sites selected
## 0b. normalize reads using SRS method for just these samples in these sites, then create phyloseq object for diversity calculations

## 1a. calculate alpha diversity using 3 metrics (Richness, Shannon's, Faith's PD)
## 1b. compare sampling windows for each metric using Kruskal-Wallis and pairwise Dunn's

## 2a. calculate community composition distances using 4 metrics (Dice-Sorensen, Unweighted Unifrac)
## 2b. ordinate each distance for:
#  2bii. between-group medians (centroids)
# 2biii. within-group dispersions (via betadisper)
## 2c. ordinate each distance metric with PCoA
## 2d. pairwise adonis and heatmap for comparing differences between particular clusters
## 2e. get the top number of detections per site+window, per taxa grouped by order

## 3. indicator species work

########################################

## start with this files: 
# otu_table_long_filtd_wMeta <- read_csv("~/github/nhguano/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz")
# import tree from local:
# tree <- read.tree(file = "~/github/nhguano/data/text_tables/trees/finalfiltd_rootedtree.nwk")

##################
## part 0a: justification for why we're analyzing particular site-windows:
otu_table_long_filtd_wMeta %>% 
  filter(Year == 2016) %>% 
  group_by(Site, Window) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>%
  arrange(Window) %>% 
  pivot_wider(names_from = "Window", values_from = "nSamples")
## 2 different window spans possible to investigate; all site-windows have minimum 6 samples per site-window: 
  ## windows 4/5/6: FOX, HOP, CNA, EPS, HOL, MAP, PEN
    ## calling these "win456"
  ## windows 5/6/7: FOX, CNB, CHI, MTV
    ## calling these "win567"

## basic stats for these two collections of samples in these two site-windows
win456_sitenames <- c("FOX", "HOP", "CNA", "EPS", "HOL", "MAP", "PEN")
win567_sitenames <- c("FOX", "CNB", "CHI", "MTV")
win456_windownames <- c(4,5,6)
win567_windownames <- c(5,6,7)

## ensure we use only Samples that have at least 2 OTUs per sample (no samples with just one OTU)
win456_samplenames <- otu_table_long_filtd_wMeta %>% 
  filter(Year == 2016 & Site %in% win456_sitenames & Window %in% win456_windownames) %>% 
  group_by(SampleID) %>% 
  summarise(nOTUs = n_distinct(OTUid)) %>%
  arrange(nOTUs) %>%
  filter(nOTUs > 1) %>% 
  select(SampleID) %>% pull()

win567_samplenames <- otu_table_long_filtd_wMeta %>% 
  filter(Year == 2016 & Site %in% win567_sitenames & Window %in% win567_windownames) %>% 
  group_by(SampleID) %>% 
  summarise(nOTUs = n_distinct(OTUid)) %>%
  arrange(nOTUs) %>%
  filter(nOTUs > 1) %>% 
  filter(SampleID != 'oro163923') %>% ## dropping this one sample because after SRS subsampling there is just 1 OTU left
  select(SampleID) %>% pull()

## for win456
otu_table_long_filtd_wMeta %>% 
  filter(SampleID %in% win456_samplenames) %>% 
  group_by(Site, Window) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>%
  arrange(Window) %>% 
  filter(nSamples >= 6 & Window %in% win456_windownames) %>%
  ungroup() %>% 
  summarise(mean_nSamples = mean(nSamples), median_nSamples = median(nSamples), sd_nSamples = sd(nSamples))
    ## mean 15.8 samples per site window; median 15, sd 9.14

## for win567
otu_table_long_filtd_wMeta %>% 
  filter(SampleID %in% win567_samplenames) %>% 
  group_by(Site, Window) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>%
  arrange(Window) %>% 
  filter(nSamples >= 6 & Window %in% win567_windownames) %>%
  ungroup() %>% 
  summarise(mean_nSamples = mean(nSamples), median_nSamples = median(nSamples), sd_nSamples = sd(nSamples))
## mean 13.8 samples per site window; median 14, sd 4.71

### going to focus on 'win456' group, as it has more sites (7 vs. 4), and more samples per site
### retained all script sfor 'win567' at "unused section" code at end of script


##################
## part 0b: rarefy table using SRS for samples in each group
## then generate phyloseq object by importing with metadata and tree file to generate phyloseq object

## get metadata for these samples and then generate OTU table for alpha/adonis/betadisp/ordinations:
win456_meta <- otu_table_long_filtd_wMeta %>% 
  filter(SampleID %in% win456_samplenames) %>% 
  select(SampleID, Window, Site, newDate) %>% 
  mutate(Window = as.factor(Window)) %>% 
  distinct()
  ## 331 samples total

## create phyloseq object:
windows2016_getPhyloseq_function <- function(metadata_file){
  tmp_sampleIDs <- metadata_file %>% select(SampleID) %>% pull()
  tmp_OTUtable_wide_raw <- as.data.frame(
    otu_table_long_filtd_wMeta %>% 
      filter(SampleID %in% tmp_sampleIDs) %>% 
      select(SampleID, OTUreads, OTUid) %>%
      pivot_wider(values_from = "OTUreads", names_from = "SampleID", values_fill = 0))
  
  row.names(tmp_OTUtable_wide_raw) <- tmp_OTUtable_wide_raw$OTUid
  tmp_OTUtable_wide_raw$OTUid <- NULL
  tmp_Cmin <- min(colSums(tmp_OTUtable_wide_raw))
  tmp_OTUtable_wide_norm <- SRS(tmp_OTUtable_wide_raw, tmp_Cmin, set_seed = 1)
  row.names(tmp_OTUtable_wide_norm) <- row.names(tmp_OTUtable_wide_raw)

  tmp_phyOTU <- otu_table(tmp_OTUtable_wide_norm, taxa_are_rows = TRUE)
  tmp_phydat <- phyloseq(tmp_phyOTU, metadata_file)
  merge_phyloseq(tmp_phydat, tree)
}

win456_phydat <- windows2016_getPhyloseq_function(win456_meta)

##################
## part 1a: alpha diversity: species richness, Shannon's entropy, Faith's PD...
windows2016_getAlphaVals_function <- function(phydat, metadata){
  tmp_faith_obs <- estimate_pd(phydat)
  tmp_shannons <- estimate_richness(phydat, measures = "Shannon") %>% rename("H" = Shannon)
  tmp_alpha_df <- data.frame(tmp_faith_obs, tmp_shannons) %>% 
    mutate(SampleID = row.names(.)) %>% 
    pivot_longer(cols = c("SR", "PD", "H"), values_to = "Alpha_value", names_to = "Metric")
  merge(tmp_alpha_df, metadata, by="SampleID")
}

win456_alpha_df <- windows2016_getAlphaVals_function(win456_phydat, win456_meta)

##################
## part 1b: applying KW test for global differences in distinct diversity metrics for a given Site+Window group
## ...then post hoc Dunn's test for pairwise diffs in alpha vals by Site+Group

windows2016_getAlphaStats_function <- function(alphadat, metric){
  tmp_kw_input <- alphadat %>% 
    mutate(Grouper = as.character(paste0(Window, Site))) %>% 
    filter(Metric == metric) %>% 
    select(Grouper, Alpha_value)
  tmp_kw_output <- kruskal.test(tmp_kw_input$Alpha_value ~ tmp_kw_input$Grouper)
  data.frame(tmp_kw_output[1], tmp_kw_output[2], tmp_kw_output[3]) %>% 
    mutate(Metric = metric,
           statistic = round(statistic, 4),
           p.value = round(p.value, 5)) %>% 
    rename("Chi_sq" = statistic,
           "df" = parameter)
}

win456_KW_result_SR <- windows2016_getAlphaStats_function(win456_alpha_df, "SR")
win456_KW_result_SR
  ##significant; chi-squared = 49.586, df = 20, p-value = 0.00025
win456_KW_result_H <- windows2016_getAlphaStats_function(win456_alpha_df, "H")
win456_KW_result_H  
##significant; chi-squared = 39.058, df = 20, p-value = 0.00656
win456_KW_result_PD <- windows2016_getAlphaStats_function(win456_alpha_df, "PD")
win456_KW_result_PD
  ##significant; chi-squared = 91.1182, df = 20, p-value = 0

## combine into supplementary table and report?
windows2016_kW_result_SR <- rbind(win456_KW_result_SR, win456_KW_result_H, win456_KW_result_PD)
rm(win456_KW_result_SR, win456_KW_result_H, win456_KW_result_PD)

write_csv(windows2016_kW_result_SR,
          file = "~/github/nhguano/data/text_tables/windows2016_KWvalues_all.csv")

rm(windows2016_kW_result_SR)

## get Dunn's pairwise vals for each group (Site+Window):
dunns_pvals_windows2016 <- function(data, metric, windowgroup){
  tmp_dunn_list <- dunnTest(Alpha_value ~ as.factor(Grouper), 
                            data=data %>% 
                              filter(Metric == metric) %>% 
                              mutate(Grouper = paste0(Window, Site)), 
                            method="bh")
  tmp_dunn_df <- data.frame(tmp_dunn_list[2]) %>% select(res.Comparison, res.P.unadj, res.P.adj) %>% 
    mutate(res.Comparison = str_replace_all(string = res.Comparison, pattern = " ", replacement = "")) %>% 
    arrange(res.Comparison) %>% 
    separate(col=res.Comparison, into = c("windowA", "windowB"), sep="-") %>% 
    mutate(Metric = metric,
           Windowgroup = windowgroup,
           res.P.unadj = round(res.P.unadj, 3),
           res.P.adj = round(res.P.adj, 3))
  tmp_dunn_df
}

win456_dunn_SR <- dunns_pvals_windows2016(win456_alpha_df, "SR", "456")
  ## among pairwise comps with pval <= 0.05 after adjusted only 10 shown, and 9 are either 4|5 EPS...
win456_dunn_H <- dunns_pvals_windows2016(win456_alpha_df, "H", "456")
  ## only 3 pairwise comps with pval <= 0.05 after adjusted... all 3 are either 5|6 HOL...
win456_dunn_PD <- dunns_pvals_windows2016(win456_alpha_df, "PD", "456")
  ## 4MAP (9) and 6HOL (8) account for half of the 29 pairwise values with p<0.05 after BH correction
    ## note the increase in number of sig pairwise diffs with PD than SR|H... why is phylogenetic diversity more different?
    ## also note that none of these sig diffs are within same sites... always between sites
    ## it is the case, however, that within Windows are represented in these sig sites... so ...
    ## seems like major difference is SITE; not Window...

windows2016_dunn_all <- rbind(win456_dunn_SR, win456_dunn_H, win456_dunn_PD)
rm(win456_dunn_SR, win456_dunn_H, win456_dunn_PD)

## save as long .csv format, but then make plot for easier viz in pub
windows2016_dunn_all <- windows2016_dunn_all %>% 
  rename("pvalue" = res.P.unadj, "BHadjusted_pvalue" = res.P.adj) %>% 
  mutate(Windowgroup = paste0("window", Windowgroup))

write_csv(windows2016_dunn_all,
          file = "~/github/nhguano/data/text_tables/windows2016_dunn_all.csv")

## easier to plot these pairwise comps than to list all of them?
windows2016_dunn_all$Metric <- factor(windows2016_dunn_all$Metric, levels = c(
  "SR", "H", "PD"))
## add asterisk to signifiy "significant" vals with p.adj <= 0.05 instead of plotting full value (too many boxes!)
windows2016_dunn_all <- windows2016_dunn_all %>% 
  mutate(sigVal = ifelse(BHadjusted_pvalue <= 0.05, "*", ""))

## plot separate windowgroups, then stitch together
ggplot(windows2016_dunn_all %>% filter(Windowgroup == "window456"), 
                     aes(x=windowB, y=windowA, fill=BHadjusted_pvalue, label=sigVal)) +
  geom_tile(color="black", size=0.25) +
  geom_text() +
  facet_grid(Windowgroup ~ Metric, space = "free") +
  scico::scale_fill_scico(palette = "bilbao") +
  coord_fixed() +
  labs(x="", y="", fill = "BH-corrected\np-value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        strip.text.y = element_blank(),
        legend.position = "top")

## save plot
ggsave("~/github/nhguano/supplementaryData/figureS5_windows2016_dunn_heatmap_pvalues_all.png", width = 21, height = 12, units = 'cm')
ggsave("~/github/nhguano/supplementaryData/figureS5_windows2016_dunn_heatmap_pvalues_all.pdf", width = 20, height = 12, units = 'cm')

## plot these groups as boxplot too, and generate alphabetical labels to show pairwise differences...
### first, generate the labels
windows2016_DunnLetters_function <- function(data, metric, windowgroup){
  tmp_dunn_list <- dunnTest(Alpha_value ~ as.factor(Grouper), method="bh",
                            data=data %>% filter(Metric == metric) %>% mutate(Grouper = paste0(Window, Site)))
  tmp_dunn_df <- data.frame(tmp_dunn_list[2]) %>% select(res.Comparison, res.P.adj) %>% 
    mutate(res.Comparison = str_replace_all(string = res.Comparison, pattern = " ", replacement = "")) %>% 
    arrange(res.Comparison)
  ngroups <- tmp_dunn_df %>% 
    mutate(Rower = row.names(.)) %>% 
    select(res.Comparison, Rower) %>% 
    separate(col=res.Comparison, into=c('windowA', 'windowB'), sep="-") %>%
    pivot_longer(-Rower) %>%
    summarise(n_distinct(value)) %>% pull()
  tmp_mat <- matrix(NA, nrow=ngroups, ncol=ngroups)
  tmp_mat[lower.tri(tmp_mat)] <- tmp_dunn_df$res.P.adj
  tmp_row1name <- tmp_dunn_df %>% 
    head(ngroups-1) %>%
    separate(col=res.Comparison, into=c('windowA', 'windowB'), sep="-") %>% 
    head(1) %>% select(windowA) %>% pull()
  tmp_rowOthersNames <- tmp_dunn_df %>% 
    head(ngroups-1) %>%
    separate(col=res.Comparison, into=c('windowA', 'windowB'), sep="-") %>% select(windowB) %>% pull()
  row.names(tmp_mat) <- paste0("window", c(tmp_row1name, tmp_rowOthersNames))
  colnames(tmp_mat) <- row.names(tmp_mat)
  tmp_mat <- Matrix::forceSymmetric(tmp_mat, uplo="L")
  tmp_lmat <- multcompLetters(tmp_mat, compare="<=", threshold=0.05, Letters=letters)
  data.frame(tmp_lmat$Letters) %>%
    rename(Letters = tmp_lmat.Letters) %>%
    mutate(Grouper = row.names(.),
           Grouper = str_remove(Grouper, "window"),
           Metric = metric,
           Window = substr(Grouper, 1, 1),
           Site = substr(Grouper, 2, 4)) %>% 
    select(-Grouper)
}
  
  
win456_dunnLetters_SR <- windows2016_DunnLetters_function(win456_alpha_df, "SR", "456")
win456_dunnLetters_H <- windows2016_DunnLetters_function(win456_alpha_df, "H", "456")
win456_dunnLetters_PD <- windows2016_DunnLetters_function(win456_alpha_df, "PD", "456")

windows2016_dunnLetters_all456 <- rbind(win456_dunnLetters_SR, win456_dunnLetters_H, win456_dunnLetters_PD)
rm(win456_dunnLetters_SR, win456_dunnLetters_H, win456_dunnLetters_PD)

#### next, create a boxplot using those Lettered labels to indicate different Dunn PW comp groups
## merge letter data (per group) with individual data (per sample) by site and window (the group!):
win456_alpha_plotdat <- 
  merge(win456_alpha_df, windows2016_dunnLetters_all456, by=c("Site", "Window", "Metric")) %>% 
  mutate(Grouper = paste0(Window, Site),
         FacePlotLabel = case_when(
           Metric == "SR" ~ "Observed OTUs", 
           Metric == "H" ~ "Shannon\'s H", 
           Metric == "PD" ~ "Faith\'s PD"))

win456_alpha_plotdat <- win456_alpha_plotdat %>% 
  group_by(Metric) %>% 
  summarise(LabelAlpha_value = max(Alpha_value) * 1.1) %>% 
  merge(., win456_alpha_plotdat, by="Metric")


## set levels
win456_alpha_plotdat$FacePlotLabel <- factor(win456_alpha_plotdat$FacePlotLabel, levels=c("Observed OTUs", "Shannon\'s H", "Faith\'s PD"))


## plot
ggplot() +
  geom_boxplot(data = win456_alpha_plotdat,
               aes(x=Window, y=Alpha_value),
               outlier.color = NA) +
  geom_jitter(data = win456_alpha_plotdat,
              aes(x=Window, y=Alpha_value),
              alpha=0.65, width = 0.1) +
  geom_text(data=win456_alpha_plotdat, 
            aes(x=Window, y=LabelAlpha_value, label=Letters),
            #angle=45, hjust=1,
            size = 2) +
  facet_grid(FacePlotLabel ~ Site, scales = "free", space = "free_x", switch = "y") +
  labs(x="", y="") +
  theme_classic() +
  theme(panel.spacing.x = unit(0.5, "lines"), panel.spacing.y = unit(3, "lines"),
        strip.placement.y = "outside",
        #strip.text.x = element_blank(), ## remove this line if you want to show per site facets at top of plot
        strip.background.y = element_blank())
        
  
ggsave("~/github/nhguano/supplementaryData/figureS4_windows2016_boxplot_alphavals_all.png", width = 27, height = 14, units = 'cm', dpi=300)
ggsave("~/github/nhguano/supplementaryData/figureS4_windows2016_boxplot_alphavals_all.pdf", width = 27, height = 14, units = 'cm', dpi=450)


##################
## part 2a:

## 2a. calculate community composition distances using 2 metrics (Dice-Sorensen, Bray-Curtis, Unweighted/Weighted Unifrac)
## 2b. ordinate each distance for:
#  2bii. between-group medians (centroids)
# 2biii. within-group dispersions (via betadisper)
## 2c. ordinate each distance metric with PCoA
## 2d. pairwise adonis and heatmap for comparing differences between particular clusters

## part 2a: community composition
## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method, per window (456 or 567)
win456_dist_ds <- phyloseq::distance(win456_phydat, "bray", binary = TRUE)
win456_dist_uu <- phyloseq::distance(win456_phydat, "unifrac", weighted=FALSE)


##################
## part 2b: community composition
## 2bi. run PERMANOVA (via Adonis) for each distance method to evaluate if samples within each Site-Window group have different centroids
windows2016_adonis_function <- function(distanceData, distanceMetric, metadata, windowgroup){
  tmp_metadata <- metadata
  tmp_metadata$Window <- as.factor(tmp_metadata$Window)
  adonis_tmp <- adonis2(distanceData ~ Site * Window, data = tmp_metadata)
  data.frame(adonis_tmp) %>% 
    mutate(Metric = distanceMetric, Class = row.names(.),
           WindowGroup = windowgroup,
           R2 = round(R2, 3),
           SumsOfSqs = round(SumOfSqs, 3),
           R2 = round(R2, 4),
           Fmodel = round(`F`, 3)) %>% 
    select(Class, Df, SumsOfSqs, R2, Fmodel, `Pr..F.`, Metric, WindowGroup) %>% 
    rename(Sum.Sq = "SumsOfSqs", F.value = "Fmodel")
}

win456_adonis_ds <- windows2016_adonis_function(win456_dist_ds, "Dice-Sorensen", win456_meta, "windows456")
win456_adonis_uu <- windows2016_adonis_function(win456_dist_uu, "UniFrac-Unweighted", win456_meta, "windows456")
win456_adonis_all <- rbind(win456_adonis_ds, win456_adonis_uu)
rm(win456_adonis_ds, win456_adonis_uu)

write_csv(win456_adonis_all,
          path = "~/github/nhguano/supplementaryData/tableS9_windows2016_adonis_allMetrics_allWindows.csv")

rm(win456_adonis_all)

##################
## part 2bii. Test for homogeneity of dispersion
  ## have to test for Site and Window main effects separately for each WindowGroup&DistanceMetric combo...
window2016_betadisper_function <- function(distanceData, distanceMetric, metadata, windowgroup){
  tmp_metadata <- metadata
  tmp_disper_list_window <- betadisper(d = distanceData, group =  tmp_metadata$Window, type = c("median"))
  tmp_disper_anova_window <- data.frame(anova(tmp_disper_list_window)) %>% 
    mutate(Class = row.names(.),
           Metric = distanceMetric, 
           Factor = "Window", 
           WindowGroup = windowgroup)
  tmp_disper_list_site <- betadisper(d = distanceData, group =  tmp_metadata$Site, type = c("median"))
  tmp_disper_anova_site <- data.frame(anova(tmp_disper_list_site)) %>% 
    mutate(Class = row.names(.),
           Metric = distanceMetric, 
           Factor = "Site", 
           WindowGroup = windowgroup)
  rbind(tmp_disper_anova_window, tmp_disper_anova_site) %>% 
    rename(SumsOfSqs = `Sum.Sq`, Fmodel = `F.value`) %>% 
    select(-`Mean.Sq`) %>% 
    relocate(Class, Factor, Df, SumsOfSqs, Fmodel, `Pr..F.`, Metric, WindowGroup) %>% 
    mutate(SumsOfSqs = round(SumsOfSqs, 3),
           Fmodel = round(Fmodel, 4),
           `Pr..F.` = round(`Pr..F.`, 5)) %>% 
    rename(Sum.Sq = "SumsOfSqs", F.value = "Fmodel")
}

win456_bdisp_ds <- window2016_betadisper_function(win456_dist_ds, "Dice-Sorensen", win456_meta, "windows456")
win456_bdisp_uu <- window2016_betadisper_function(win456_dist_uu, "Unifrac-Unweighted", win456_meta, "windows456")

windows2016_bdisp_all <- rbind(win456_bdisp_ds, win456_bdisp_uu)
rm(win456_bdisp_ds, win456_bdisp_uu)
write_csv(windows2016_bdisp_all,
          path = "~/github/nhguano/supplementaryData/tableS10_windows2016_betaDispersion_allMetrics_allWindows.csv")


##################
## part 2c: ordinations
## gather data for ordination per distance metric and window group
windows2016_ordi_function <- function(distanceValues, distanceMetric, phydat, windowgroup, metadata){
  tmp_metadata <- metadata
  tmp_pcoa <- ordinate(phydat, method="PCoA", distance = distanceValues)
  tmp_pcoa_list <- plot_ordination(phydat, tmp_pcoa)
  tmp_pcoa_df <- tmp_pcoa_list$data %>%
    mutate(SampleID = row.names(.),
           Axis1lab = tmp_pcoa_list$labels$x,
           Axis2lab = tmp_pcoa_list$labels$y,
           Metric = distanceMetric)
  merge(tmp_pcoa_df, tmp_metadata) %>% 
    mutate(Window = as.character(Window),
           WindowGroup = windowgroup)
}

win456_orddata_ds <- windows2016_ordi_function(win456_dist_ds, "Dice-Sorensen", win456_phydat, "win456", win456_meta)
win456_orddata_uu <- windows2016_ordi_function(win456_dist_uu, "Unweighted-UniFrac", win456_phydat, "win456", win456_meta)

## for all plot types, set up these color/shape parameters
## palette to match windows from other plots!
win456_pal <- rev(c("#557d3f", "#B28C32", "#FCA68C"))

## need to define shape point types; ensure overlapping Sites have same point value
win456_shapes <- c(9, 3, 17, 15, 8, 16, 11)

### 1) single plots per WindowGroup & Distance estimate
window2016_plotfunction <- function(ordinationdata, shapepal, colorpal){
  
  ordinationdata$Metric <- factor(ordinationdata$Metric, levels = c(
    "Dice-Sorensen", "Unweighted-UniFrac"))
  
  ggplot(data = ordinationdata,
         aes(x=Axis.1, y=Axis.2, color=Window, shape = Site)) +
    geom_point(aes(color = Window, shape = Site), size = 1.5) +
    stat_ellipse(aes(group = Window, color=Window), alpha = 0.5) +
    labs(x = unique(ordinationdata$Axis1lab), y = unique(ordinationdata$Axis2lab)) +
    theme_classic() +
    facet_wrap(~Metric) +
    scale_color_manual(values = colorpal) +
    scale_shape_manual(values = shapepal) +
    coord_fixed() +
    guides(shape = guide_legend(order = 0), colour = guide_legend(order = 1))
}

p_win456_ord_ds <- window2016_plotfunction(win456_orddata_ds, win456_shapes, win456_pal)
p_win456_ord_ds
ggsave(filename = "~/github/nhguano/figures/figure4ai_pcoa_win456_ds.png", width = 17, height = 15, unit = "cm", dpi=150)
ggsave(filename = "~/github/nhguano/figures/figure4ai_pcoa_win456_ds.pdf", width = 17, height = 15, unit = "cm", dpi=300)

p_win456_ord_uu <- window2016_plotfunction(win456_orddata_uu, win456_shapes, win456_pal)
p_win456_ord_uu
ggsave(filename = "~/github/nhguano/figures/figure4aii_pcoa_win456_uu.png", width = 17, height = 15, unit = "cm", dpi=150)
ggsave(filename = "~/github/nhguano/figures/figure4aii_pcoa_win456_uu.pdf", width = 17, height = 15, unit = "cm", dpi=300)


### 2) combined plots with different distance estimates for data with shared WindowGroups
tmp_win456_pcoa_legend <- cowplot::get_legend(p_win456_ord_ds)

p_windows2016_pcoa_left <- plot_grid(
  p_win456_ord_ds + theme(legend.position = "none"), 
  p_win456_ord_uu + theme(legend.position = "none"),
  nrow = 1, ncol = 2, labels = c("A", "B"),
  align = "hv")

plot_grid(p_windows2016_pcoa_left, tmp_win456_pcoa_legend, ncol=2, nrow=1, 
          rel_widths = c(0.9, 0.1))

ggsave(filename = "~/github/nhguano/figures/figure4a_pcoa_windows2016_allMetrics.png", width = 27, height = 12, unit = "cm", dpi=150)
ggsave(filename = "~/github/nhguano/figures/figure4a_pcoa_windows2016_allMetrics.pdf", width = 27, height = 12, unit = "cm", dpi = 300)
ggsave(filename = "~/github/nhguano/figures/figure4a_pcoa_windows2016_allMetrics.svg", width = 27, height = 12, unit = "cm", dpi = 300)

rm(tmp_win456_pcoa_legend, p_windows2016_pcoa_left)

##################
## part 2d: pairwise adonis

## testing for groups as distinct site+window combinations like with Dunn's for alpha diversity above...
window2016_pairwise_adonis_func <- function(distancevals, metric, windowgroup, metadata){
  tmp_pairwiseadonis_window <- 
    metadata_vars <- paste0(metadata$Window, metadata$Site)
    pairwise.adonis(distancevals, 
                    metadata_vars,
                    p.adjust.m = "BH") %>%
    mutate(Metric = metric, WindowGroup = windowgroup, Factor = "WindowSite") %>%
    separate(col = pairs, into = c("groupA", "groupB"), sep = "vs") %>%
    arrange(groupA, groupB)
}

win456_pwadonis_ds <- window2016_pairwise_adonis_func(win456_dist_ds, "Dice-Sorensen", "windows456", win456_meta)
win456_pwadonis_uu <- window2016_pairwise_adonis_func(win456_dist_uu, "Unweighted-UniFrac", "windows456", win456_meta)

windows2016_pwadonis_all <- rbind(win456_pwadonis_ds, win456_pwadonis_uu)
rm(win456_pwadonis_ds, win456_pwadonis_uu)
write_csv(windows2016_pwadonis_all,
          path = "~/github/nhguano/data/text_tables/windows2016_pwAdonis_allmetrics.csv")

## cleanup
rm(list=ls(pattern = "^win.*._dist.*"))
rm(list=ls(pattern = ".*._sitenames"))
rm(list=ls(pattern = ".*._samplenames"))
rm(list=ls(pattern = ".*._windownames"))
rm(list=ls(pattern = "^window"))
rm(dunns_pvals_windows2016)

##################
## part 2e: what taxa are different between these Site+Window groups?
## identify the proportion of detections BY arthropod ORDER

## gather data from same srs sampled dataset - makes directly comparable to PCoA input (distance matrix calcualted from this data)
# ## pivot wide to long for SRS output, merge with taxonomy and metadata info
get_windows2016_srsdata_wtaxa_function <- function(phydat_data, metadata){
  tmp_otutable_wide <- as.data.frame(otu_table(phydat_data))
  tmp_otutable_wide$OTUid <- row.names(tmp_otutable_wide)
  tmp_otutable_long <- tmp_otutable_wide %>% 
    pivot_longer(-OTUid, values_to = 'SRSreads', names_to = 'SampleID') %>% 
    filter(SRSreads > 0)
  tmp_taxa <- otu_table_long_filtd_wMeta %>%
    filter(SampleID %in% tmp_otutable_long$SampleID) %>%
    distinct(OTUid, OTUalias, Class, Order, Family, Genus, Species)
  tmp_otutable_long_wTaxa <- merge(tmp_otutable_long, tmp_taxa)
  merge(tmp_otutable_long_wTaxa, metadata)
}

win456_otutable_long_wtaxa_wmeta <- get_windows2016_srsdata_wtaxa_function(win456_phydat, win456_meta)

####### first plot is at Order level...
windows2016_getTopOrderData_function <- function(otutable_long){
  ## for all data, across all sampling windows, how many detections total PER SITE?
  tmp_nGlobalDetections <- otutable_long %>% group_by(Site) %>% summarise(nDetections_perSite = n())
  ## retain only those arth Orders where fraction of detections >= 2% in at least 1 site (arth Order can exist only in single site!)
  tmp_nDetections_perSite_perOrder_all <- otutable_long %>% 
    group_by(Site, Class, Order) %>% 
    summarise(nDetections_perOrder = n()) %>% 
    merge(., tmp_nGlobalDetections) %>% 
    mutate(pDetections_perSite = nDetections_perOrder / nDetections_perSite)
  tmp_topOrders_detected <- tmp_nDetections_perSite_perOrder_all %>% 
    filter(pDetections_perSite >= 0.02) %>% 
    distinct(Order) %>% 
    pull()
  
  ## for each site+window combination, how many detections total?
  tmp_nSiteWindow_detections_global <- otutable_long %>% group_by(Site, Window) %>% summarise(nDetections_perSiteWindow = n())
  ## figure out fraction of orders per site+window group, aggregating those "other" orders less than 2% of detections globally
  otutable_long %>% 
    mutate(Order = case_when(!Order %in% tmp_topOrders_detected ~ "other", ## this is the aggregation step for other Orders
                             TRUE ~ as.character(Order))) %>% 
    group_by(Site, Window, Order) %>% 
    summarise(nDetections_perOrder_perSiteWindow = n()) %>% 
    merge(., tmp_nSiteWindow_detections_global, by = c("Site", "Window")) %>% 
    mutate(pDetections_perSiteWindow = nDetections_perOrder_perSiteWindow / nDetections_perSiteWindow) %>% 
    select(Site, Window, Order, pDetections_perSiteWindow) %>% 
    mutate(Grouper = paste0(Window, Site))
}

win456_OrderDetect_plotdat <- windows2016_getTopOrderData_function(win456_otutable_long_wtaxa_wmeta)

## keep order palete colors consistent between plots!
win456_topOrders <- win456_OrderDetect_plotdat %>% distinct(Order) %>% arrange(Order) %>% pull
#10+1 Orders included are: "Araneae" "Coleoptera" "Diptera" "Ephemeroptera" "Hemiptera" "Lepidoptera" "Megaloptera" "Psocodoea" "Trichoptera" "Trombidiformes" 'other'
win456_orderPal <- c("gray25", "orchid2", "tan4", "#3E6C54", "#FDAC9F", "#808133", 'cadetblue3', "#e1b580", "turquoise4", 'gray50', '#e9ee7a')

## order these Orders to place the aggregated 'other' at end of legend:
win456_OrderDetect_plotdat$Order <- factor(win456_OrderDetect_plotdat$Order, levels = c(
  "Araneae", "Coleoptera", "Diptera", "Ephemeroptera", "Hemiptera", "Lepidoptera", 
  "Megaloptera", "Psocodea", "Trichoptera", "Trombidiformes", "other"))

# plot, faceting site groups together and ordered by sampling window within facet
ggplot(win456_OrderDetect_plotdat,
       aes(x=Window, y=pDetections_perSiteWindow, fill=Order)) + 
  geom_col() + 
  scale_fill_manual(values = win456_orderPal) +
  facet_grid(~Site, space="free_x", scales = "free_x") +
  labs(x="\n sampling Window", y = "fraction of detections\n") +
  theme_classic() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow=2))

ggsave("~/github/nhguano/figures/figure4b_windows2016_topOrder.png", height = 10, width = 20, dpi = 150, units="cm")
ggsave("~/github/nhguano/figures/figure4b_windows2016_topOrder.pdf", height = 10, width = 20, dpi = 300, units="cm")


##################
## part 2e:
## indicator species work: 

##setting up a palette for our association statistic values upfront:
### keep scal consistent across all different plot types
windows2016_assocPal <- scico::scico(3, palette = 'lajolla')
windows2016_assocPal

## what are the orders we're going to keep for the indispec analyses?
## focusing on the 10 orders listed in fig2d - the ones with at least 2% detections across any site+window
top456Orders <- win456_OrderDetect_plotdat %>% distinct(Order) %>% filter(Order != "other") %>% pull()

get_indispecDataTable_function <- function(metadata, orderlist){
  ## collapse read counts to shared genus labels, then reformat for SRS sampling
  tmp_forSRS_df <- as.data.frame(otu_table_long_filtd_wMeta %>% 
                                   filter(SampleID %in% metadata$SampleID) %>% ## retain same samples as in Physeq object
                                   filter(Order %in% orderlist) %>% 
                                   mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>% ## rename genus labels to avoid any NA aggregating
                                   group_by(SampleID, Order, Genus) %>% 
                                   summarise(taxaReads = sum(OTUreads)) %>%
                                   mutate(taxaName = paste0(Order, "-", Genus)) %>% 
                                   ungroup() %>% 
                                   select(SampleID, taxaName, taxaReads) %>% 
                                   pivot_wider(names_from = "SampleID", values_from = "taxaReads", values_fill = 0))
  
  ## perform srs sampling    
  row.names(tmp_forSRS_df) <- tmp_forSRS_df$taxaName
  tmp_forSRS_df$taxaName <- NULL
  tmp_Cmin <- min(colSums(tmp_forSRS_df))
  tmp_SRSout_norm <- SRS(tmp_forSRS_df, tmp_Cmin, set_seed = 1)
  row.names(tmp_SRSout_norm) <- row.names(tmp_forSRS_df)
  
  ## flip axes of matrix for indispecies program
  tmp_indi_mat <- t(tmp_SRSout_norm)
}  

win456_indiDataTable <- get_indispecDataTable_function(win456_meta, top456Orders)

### perform indicator species analysis three different ways:
### (A) by Windows, (B) by Site, or (C) by Site+Window combination
### can modify this function below, but our main plot will focus just on those singleton groups...
  ### ... that is, those indicator species exclusive to a particular Site, or Window, or Site+Window
  ### if interested to explore further group combinations, just modify the `maxorderval` in the argument to generate the data

#### by Windows first:
doIndiSpec_function <- function(datatable, metadata, metagroup, maxorderval){
  tmp_indval_list <- multipatt(datatable, 
                               metadata %>% 
                                 select(metagroup) %>% pull(), 
                               control = how(nperm=999),
                               max.order = maxorderval)
}

win456_indispec_list_windows <- doIndiSpec_function(win456_indiDataTable, win456_meta, "Window", 1)
win456_indispec_plotdat_windows <- win456_indispec_list_windows$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(win456_indispec_list_windows$comb)) %>%
          mutate(index = row.names(.)))

## set levels for plot, ordering group that spans windows between distinct windows; 
## plot
p_win456_indispec_windows <- ggplot(win456_indispec_plotdat_windows,
                                    aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(. ~ Order, scales = "free_x", switch = "y", space="free") +
  labs(y="window", x="", fill = "association\nstatistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(angle = 45, vjust=0.5),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"))
p_win456_indispec_windows
ggsave("~/github/nhguano/figures/figure4c_specindi_byWindow_win456.png", width = 250, height = 100, units="mm", dpi=150)
ggsave("~/github/nhguano/figures/figure4c_specindi_byWindow_win456.pdf", width = 250, height = 100, units="mm", dpi=300)

## repeat for Sites (instead of window groups)
win456_indispec_list_sites <- doIndiSpec_function(win456_indiDataTable, win456_meta, "Site", 1)
win456_indispec_plotdat_sites <- win456_indispec_list_sites$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(win456_indispec_list_sites$comb)) %>% ## gets the index names to match group levels
          mutate(index = row.names(.)))

# ## set levels for x axis:
win456_indispec_plotdat_sites$indexName <-
  factor(win456_indispec_plotdat_sites$indexName,
         levels = win456_indispec_plotdat_sites %>% distinct(indexName, index) %>% arrange(index) %>% pull(indexName))


## plot
p_win456_indispec_sites <- ggplot(win456_indispec_plotdat_sites,
                                  aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(. ~ Order, scales = "free_x", switch = "y", space="free") +
  labs(y="window", x="", fill = "association\nstatistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(angle = 45, vjust=0.5),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"))
p_win456_indispec_sites
ggsave("~/github/nhguano/figures/figure4c_specindi_bySite_win456.png", width = 250, height = 125, units="mm", dpi=150)
ggsave("~/github/nhguano/figures/figure4c_specindi_bySite_win456.pdf", width = 250, height = 125, units="mm", dpi=300)

## let's try it for the various Site+Window groups too!
win456_indispec_list_sitewindows <- doIndiSpec_function(win456_indiDataTable, 
                                                        win456_meta %>% mutate(SiteWindow = paste0(Site,Window)), 
                                                        "SiteWindow", 1)
win456_indispec_plotdat_sitewindows <- win456_indispec_list_sitewindows$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(win456_indispec_list_sitewindows$comb)) %>% ## gets the index names to match group levels
          mutate(index = row.names(.))) %>% 
  mutate(Site = substr(indexName, 1, 3),
         Window = substr(indexName, 4, 4))

# ## set levels for x axis:
win456_indispec_plotdat_sitewindows$indexName <-
  factor(win456_indispec_plotdat_sitewindows$indexName,
         levels = win456_indispec_plotdat_sitewindows %>% distinct(indexName, index) %>% arrange(index) %>% pull(indexName))


## plot
p_win456_indispec_sitewindows <- ggplot(win456_indispec_plotdat_sitewindows,
                                        aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(Site ~ Order, scales = "free", switch = "y", space="free") +
  labs(y="window", x="", fill = "association\nstatistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(angle = 45, vjust=0.5),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"))
p_win456_indispec_sitewindows
ggsave("~/github/nhguano/figures/figure4c_specindi_bySiteWindow_win456.png", width = 250, height = 200, units="mm", dpi=150)
ggsave("~/github/nhguano/figures/figure4c_specindi_bySiteWindow_win456.pdf", width = 250, height = 200, units="mm", dpi=300)

###############
## grouping all these plots together so they share a single x-axis (common taxa/genus labels)

## first, gather the metadata needed from the output of the individual indispec analyses and merge into single data frame
win456_winfilt <- win456_meta %>% distinct(Window) %>% pull()
win456_sitefilt <- win456_meta %>% distinct(Site) %>% pull()
win456_sitewindow_filt <- win456_meta %>% mutate(SiteWindow = paste0(Site,Window)) %>% distinct(SiteWindow) %>% pull()

win456_indispec_plotdat_allGroups <- 
  rbind(win456_indispec_plotdat_windows %>% filter(indexName %in% win456_winfilt) %>% mutate(groupLabel = "window"),
        win456_indispec_plotdat_sites %>% filter(indexName %in% win456_sitefilt) %>% mutate(groupLabel = "site"),
        win456_indispec_plotdat_sitewindows %>% filter(indexName %in% win456_sitewindow_filt) %>% select(-Site, -Window) %>% mutate(groupLabel = "site+window"))

## order these "grouping" levels of 'window', 'site', and 'window+site' groups...
win456_indispec_plotdat_allGroups$groupLabel <- factor(
  win456_indispec_plotdat_allGroups$groupLabel,
  levels = c("window", "site", "site+window"))

# plot
p_win456_indispec_alldat_wide <- ggplot(win456_indispec_plotdat_allGroups,
                                        aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_x_discrete(position = "bottom") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(groupLabel ~ Order, scales = "free", space="free", switch = "y") +
  labs(x="", y="", fill = "association\nstatistic") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.placement = "outside",
        aspect.ratio = 1,
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"),
        axis.text.x = element_text(angle=90, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.spacing.y = unit(0.5, "lines"))
p_win456_indispec_alldat_wide
## goint to need to stagger the x.axis facet labels so they don't get cut off between the panels...
ggsave("~/github/nhguano/figures/figure4c_specindi_allGroups_win456_WIDE.png", width = 250, height = 150, units="mm", dpi=150)
ggsave("~/github/nhguano/figures/figure4c_specindi_allGroups_win456_WIDE.pdf", width = 250, height = 150, units="mm", dpi=300)


### as a supplementary file, going to also do this species indicator analysis to include a max of 2 layers
  ## going to investigate how many significant hits for these site+window combos - do most sitll aggregate to HOL?
win456_indispec_list_window_multigroup <- doIndiSpec_function(win456_indiDataTable, win456_meta, "Window", 2)
win456_indispec_list_site_multigroup <- doIndiSpec_function(win456_indiDataTable, win456_meta, "Site", 2)
win456_indispec_list_sitewindows_multigroupOkay <- doIndiSpec_function(win456_indiDataTable, 
                                                        win456_meta %>% mutate(SiteWindow = paste0(Site,Window)), 
                                                        "SiteWindow", 2)

getplotdat_indispec_function <- function(indiplotdat_list, grouplabeler){
  indiplotdat_list$sign %>% 
    mutate(taxaName = row.names(.)) %>% 
    separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
    filter(p.value <= 0.05) %>% 
    select(index, stat, Order, Genus) %>% 
    merge(., 
          data.frame(indexName = colnames(indiplotdat_list$comb)) %>%
            mutate(index = row.names(.))) %>% 
    mutate(GroupLabel = grouplabeler)
}

win456_indispec_plotdat_window_multigroupOkay <- getplotdat_indispec_function(win456_indispec_list_window_multigroup, "Window")
win456_indispec_plotdat_site_multigroupOkay <- getplotdat_indispec_function(win456_indispec_list_site_multigroup, "Site")
win456_indispec_plotdat_sitewindows_multigroupOkay <- getplotdat_indispec_function(win456_indispec_list_sitewindows_multigroupOkay, "Site+Window")

## combine and retain only those hits that are a combination:
win456_indispec_plotdat_allGroups_multigroupOnly <- 
  bind_rows(win456_indispec_plotdat_window_multigroupOkay,
            win456_indispec_plotdat_site_multigroupOkay,
            win456_indispec_plotdat_sitewindows_multigroupOkay) %>% 
  filter(grepl("\\+", indexName))

## now plot for supplementary file
win456_indispec_plotdat_allGroups_multigroupOnly$GroupLabel <-
  factor(win456_indispec_plotdat_allGroups_multigroupOnly$GroupLabel,
         levels = c("Window", "Site", "Site+Window"))

ggplot(win456_indispec_plotdat_allGroups_multigroupOnly,
       aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_x_discrete(position = "bottom") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(GroupLabel ~ Order, scales = "free", space="free", switch = "y") +
  labs(x="", y="", fill = "association\nstatistic") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10, angle=90, hjust=1),
        strip.text.y.left = element_text(size = 10, angle=0, hjust=1),
        strip.placement = "outside",
        aspect.ratio = 1,
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"),
        axis.text.x = element_text(angle=90, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.spacing.y = unit(0.5, "lines"))
ggsave("~/github/nhguano/supplementaryData/figureS6_specindi_multigroupComps_win456_WIDE.png", width = 200, height = 150, units="mm", dpi=150)
ggsave("~/github/nhguano/supplementaryData/figureS6_specindi_multigroupComps_win456_WIDE.pdf", width = 200, height = 150, units="mm", dpi=300)


#### one last little summary - group these significant indicator species results by the Group (window | site | site+window)
### tally up the number of significant indicators distinct to a given factor within the group
win456_indispec_plotdat_allGroups_multigroupOkay <- 
  bind_rows(win456_indispec_plotdat_window_multigroupOkay,
            win456_indispec_plotdat_site_multigroupOkay,
            win456_indispec_plotdat_sitewindows_multigroupOkay)

win456_indiDataTable_summaryTopHits_perGroupLabel <- win456_indispec_plotdat_allGroups_multigroupOkay %>% 
  group_by(GroupLabel, indexName) %>% 
  summarise(nTotalSigTaxa = n()) %>% 
  arrange(GroupLabel, -nTotalSigTaxa) %>% 
  rename(GroupType = "GroupLabel", Group = 'indexName')

write_csv(win456_indiDataTable_summaryTopHits_perGroupLabel,
          file = "~/github/nhguano/supplementaryData/tableS11_indiDataTable_summaryTopHitsPerGroupLabel.csv")

## we can see here that:
## for Window group... window6 has the bulk of distinct genus-level taxa (25 of 37 identified)
## for Site group .... HOL has bulk of distinct genus-level taxa (19 of 42)
## Site+Window groups mostly unique to HOL ...

## what fraction of detections do these indicator species make up?
## just because these taxa are associated with a group, are they still pretty rare?

### function to grab these data:
get_longSRSoutput_forindicTest <- function(metadata, orderlist){
  tmp_forSRS_df <- as.data.frame(otu_table_long_filtd_wMeta %>% 
                                   filter(SampleID %in% metadata$SampleID) %>% ## retain same samples as in Physeq object
                                   filter(Order %in% orderlist) %>% 
                                   mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>% ## rename genus labels to avoid any NA aggregating
                                   group_by(SampleID, Order, Genus) %>% 
                                   summarise(taxaReads = sum(OTUreads)) %>%
                                   mutate(taxaName = paste0(Order, "-", Genus)) %>% 
                                   ungroup() %>% 
                                   select(SampleID, taxaName, taxaReads) %>% 
                                   pivot_wider(names_from = "SampleID", values_from = "taxaReads", values_fill = 0))
  
  ## perform srs sampling, then pivot back to long format 
  row.names(tmp_forSRS_df) <- tmp_forSRS_df$taxaName
  tmp_forSRS_df$taxaName <- NULL
  tmp_Cmin <- min(colSums(tmp_forSRS_df))
  tmp_SRSout_norm <- SRS(tmp_forSRS_df, tmp_Cmin, set_seed = 1)
  row.names(tmp_SRSout_norm) <- row.names(tmp_forSRS_df)
  tmp_SRSout_norm$TaxaLabel <- row.names(tmp_SRSout_norm)
  tmp_SRSout_long <- tmp_SRSout_norm %>% 
    pivot_longer(-TaxaLabel, names_to = "SampleID", values_to = "TaxaReads") %>% 
    filter(TaxaReads > 0) %>% 
    merge(., metadata %>% select(-newDate))
}

win456_srsLong_genusTaxaCounts <- get_longSRSoutput_forindicTest(win456_meta, win456_topOrders)

## figure out, for each GroupType, what fraction of samples are in that state
### get the taxa names per GroupType (window, site, site+windnow), then count what fraction of detections those taxa constitute in a their respective group
win456_indispec_windowHitTaxaNames <- win456_indispec_plotdat_allGroups_multigroupOkay %>% 
  filter(GroupLabel == "Window") %>% 
  distinct(Order, Genus) %>% 
  mutate(TaxaLabel = paste0(Order,"-",Genus)) %>% 
  select(TaxaLabel) %>% pull()

win456_indispec_siteHitTaxaNames <- win456_indispec_plotdat_allGroups_multigroupOkay %>% 
  filter(GroupLabel == "Site") %>% 
  distinct(Order, Genus) %>% 
  mutate(TaxaLabel = paste0(Order,"-",Genus)) %>% 
  select(TaxaLabel) %>% pull()

win456_indispec_sitewindowHitTaxaNames <- win456_indispec_plotdat_allGroups_multigroupOkay %>% 
  filter(GroupLabel == "Site+Window") %>% 
  distinct(Order, Genus) %>% 
  mutate(TaxaLabel = paste0(Order,"-",Genus)) %>% 
  select(TaxaLabel) %>% pull()

### then apply this to create a new DataFrame
win456_srsLong_genusTaxaCounts <- win456_srsLong_genusTaxaCounts %>%
  mutate(sigfor_window = case_when(TaxaLabel %in% win456_indispec_windowHitTaxaNames ~ "TRUE",
                                   TRUE ~ as.character("FALSE")),
         sigfor_site = case_when(TaxaLabel %in% win456_indispec_siteHitTaxaNames ~ "TRUE",
                                   TRUE ~ as.character("FALSE")),
         sigfor_sitewindow = case_when(TaxaLabel %in% win456_indispec_sitewindowHitTaxaNames ~ "TRUE",
                                   TRUE ~ as.character("FALSE")))

### now use the detections (not TaxaRead sums) to aggregate numbers of detections for each GroupType
### 

#sigfor_label <- sigfor_window
#grouptype = "Window" (or "site", etc.)

tmp_detectionsPerGroup_window <-
  win456_srsLong_genusTaxaCounts %>%
  group_by(sigfor_window, Window) %>% ## adjust for status
  summarise(nDetections_perGroup = n(),
            pDetections_perGroup = nDetections_perGroup/nrow(win456_srsLong_genusTaxaCounts)) %>% 
  mutate(pDetections_perGroup = round(pDetections_perGroup, 2)) %>% 
  select(-nDetections_perGroup) %>% 
  rename(indicTaxa = 'sigfor_window', Group = "Window") %>% ## adjust
  arrange(Group, rev(indicTaxa))

tmp_detectionsPerGroup_site <-
  win456_srsLong_genusTaxaCounts %>%
  group_by(sigfor_site, Site) %>% ## adjust for status
  summarise(nDetections_perGroup = n(),
            pDetections_perGroup = nDetections_perGroup/nrow(win456_srsLong_genusTaxaCounts)) %>% 
  mutate(pDetections_perGroup = round(pDetections_perGroup, 2)) %>% 
  select(-nDetections_perGroup) %>% 
  rename(indicTaxa = 'sigfor_site', Group = "Site") %>% ## adjust
  arrange(Group, rev(indicTaxa))

tmp_detectionsPerGroup_sitewindow <-
  win456_srsLong_genusTaxaCounts %>%
  mutate(SiteWindow = paste0(Site, Window)) %>% 
  group_by(sigfor_sitewindow, SiteWindow) %>% ## adjust for status
  summarise(nDetections_perGroup = n(),
            pDetections_perGroup = nDetections_perGroup/nrow(win456_srsLong_genusTaxaCounts)) %>% 
  mutate(pDetections_perGroup = round(pDetections_perGroup, 2)) %>% 
  select(-nDetections_perGroup) %>% 
  rename(indicTaxa = 'sigfor_sitewindow', Group = "SiteWindow") %>% ## adjust
  arrange(Group, rev(indicTaxa))

#### overall appears pretty balanced between TRUE/FALSE indicator species taxa vs. non
  ### HOL appears to have more detections overall associated with that site though.

#### could also look at this just in terms of TRUE/FALSE, without further differentiating by the particular window/site/site+window
win456_srsLong_genusTaxaCounts %>%
  group_by(sigfor_window) %>% ## adjust for status
  summarise(nDetections_perGroup = n(),
            pDetections_perGroup = nDetections_perGroup/nrow(win456_srsLong_genusTaxaCounts))

win456_srsLong_genusTaxaCounts %>%
  group_by(sigfor_site) %>% ## adjust for status
  summarise(nDetections_perGroup = n(),
            pDetections_perGroup = nDetections_perGroup/nrow(win456_srsLong_genusTaxaCounts))

win456_srsLong_genusTaxaCounts %>%
  group_by(sigfor_sitewindow) %>% ## adjust for status
  summarise(nDetections_perGroup = n(),
            pDetections_perGroup = nDetections_perGroup/nrow(win456_srsLong_genusTaxaCounts))


## cleanup:
rm(list = ls(pattern = "^win456"))
rm(list = ls(pattern = "^p_win456"))
rm(list = ls(pattern = "^sigIndi"))
rm(indispec_distinctHits_perGroupFactor, top456Orders, windows2016_assocPal)
rm(doIndiSpec_function, get_indispecDataTable_function, 
   get_longSRSoutput_forindicTest, get_windows2016_srsdata_wtaxa_function,
   getplotdat_indispec_function, windows2016_getTopOrderData_function)


########################################
## part 5 - evaluate between between sites in single window of sampling season between 2015 and 2016 years
########################################

## start with this files: 
# otu_table_long_filtd_wMeta <- read_csv("~/github/nhguano/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz")
# import tree from local:
# tree <- read.tree(file = "~/github/nhguano/data/text_tables/trees/finalfiltd_rootedtree.nwk")


# ##################
# ## part 5a: In comparing between 2015 qnd 2016... Which sites make sense? Which months?
otu_table_long_filtd_wMeta %>%
  group_by(Site, Year, Window) %>%
  summarise(nSamples = n_distinct(SampleID)) %>%
  pivot_wider(names_from = "Year", values_from = "nSamples") %>%
  filter(!is.na(`2015`)) %>%
  filter(!is.na(`2016`)) %>%
  arrange(Window, Site)

## focusing on just 3 sites at single sampling window: COR and HOP and MAP
  ## Site  Window  `2016` `2015`
  ## COR        6     20     12
  ## HOP        6     18     13
  ## MAP        6     12      9

# ##################
# ## part 5b: select those sites we'll analyze and run SRS to noramlize counts
multiyear_site_names <- c("COR", "HOP", "MAP")
multiyear_otu_table_long <- otu_table_long_filtd_wMeta %>% 
  filter(Site %in% multiyear_site_names & Window == 6) %>% 
  filter(SampleID != 'oro160680') %>% filter(SampleID != 'oro160679') %>% ## dropping these two samples because they are filled almost entirely with mite DNA
  filter(SampleID != 'oro151711') ## dropping this one sample because after SRS sampling we have just a single OTU detected 

# ## get metadata just for those sites:
multiyear_meta <- multiyear_otu_table_long %>% select(SampleID, Site, Year, Window, SiteLat, SiteLong) %>% distinct()
#

### how many samples per group? (group == site+year)
multiyear_meta %>% 
  group_by(Site, Year) %>% 
  tally() %>% 
  ungroup() %>% 
  summarise(meanSamps = mean(n), medianSamps = median(n), sdSamps = sd(n))
## mean=13.5; median=12; sd=4.46 samples per group

# ## calculating alpha and beta diversity measures; what sampling depth to use? keep at 1000 seqs? increase to min of these samps?
multiyear_otu_table_long %>%
  group_by(SampleID, Site, Year) %>%
  summarise(readsPerSample = sum(OTUreads)) %>% 
  arrange(readsPerSample) %>% head()
  # ## yep; 1001 reads is lowest, with multiple samples in that range of abundances... comparable to all other tests

## use SRS to subsample data prior to diversity tests
multiyear_otu_table_wide_raw <- multiyear_otu_table_long %>%
  select(OTUid, OTUreads, SampleID) %>%
  dcast(data = ., formula = OTUid ~ SampleID, value.var='OTUreads', fill = 0)
# 
row.names(multiyear_otu_table_wide_raw) <- multiyear_otu_table_wide_raw$OTUid
multiyear_otu_table_wide_raw$OTUid <- NULL
multiyear_otu_table_wide_raw <- as.data.frame(multiyear_otu_table_wide_raw)
## get Cmin
multiyear_Cmin <- min(colSums(multiyear_otu_table_wide_raw))
# ## run SRS
multiyear_otu_table_wide_norm <- SRS(multiyear_otu_table_wide_raw, multiyear_Cmin, set_seed = 1)
## add featureIDs back into the SRSoutput table
row.names(multiyear_otu_table_wide_norm) <- row.names(multiyear_otu_table_wide_raw)
rm(multiyear_otu_table_wide_raw, multiyear_Cmin, multiyear_otu_table_long)
## create long-formatted info for later:
multiyear_otu_table_long_norm <- 
  multiyear_otu_table_wide_norm %>% 
  mutate(OTUid = row.names(.)) %>% 
  pivot_longer(-OTUid, values_to = "OTUreads", names_to = "SampleID") %>% 
  filter(OTUreads > 0)
## add in metadata for those subset of samples:
multiyear_otu_table_long_norm_wMeta <- merge(multiyear_otu_table_long_norm, multiyear_meta)
## add in taxa for those subset of samples:
multiyear_taxa_tmp <- otu_table_long_filtd_wMeta %>% 
  filter(SampleID %in% multiyear_otu_table_long_norm$SampleID) %>% 
  select(OTUid, Class, Order, Family, Genus, Species, OTUalias) %>% 
  distinct()
multiyear_otu_table_long_norm_wMeta_andTaxa <- merge(multiyear_otu_table_long_norm_wMeta, multiyear_taxa_tmp)
rm(multiyear_taxa_tmp, multiyear_otu_table_long_norm_wMeta, multiyear_otu_table_long_norm)

# 
# ##################
# ## part 5c: diversity estimates using SRS output calculated with Phyloseq to incorporate tree info for both alpha and distance calcs
# ## create phyloseq object
# ## import OTU data in phyloseq format
multiyear_phyOTU <- otu_table(multiyear_otu_table_wide_norm, taxa_are_rows = TRUE)
## bundle OTU table, metadata, and tree info into single phyloseq object
myltiyear_phymeta <- as.data.frame(multiyear_meta)
myltiyear_phymeta$SampleID -> row.names(myltiyear_phymeta)
myltiyear_phymeta <- sample_data(myltiyear_phymeta)
multiyear_phydat <- phyloseq(multiyear_phyOTU, tree)
multiyear_phydat <- merge_phyloseq(multiyear_phydat, myltiyear_phymeta)
rm(multiyear_phyOTU, multiyear_otu_table_wide_norm, myltiyear_phymeta)
# 
# ##################
# ## part 1a: alpha diversity estimates: observed OTUs, Shannon's, Faith's PD
# ## need to estimate Shannon's separate from Faith's/Observed
# ## get Observed and Faith's first
multiyear_faith_obs <- estimate_pd(multiyear_phydat)
## now get Shannon's entropy
multiyear_shannons <- estimate_richness(multiyear_phydat, measures = "Shannon") %>% rename("H" = Shannon)
## merge, then combine with metadata
multiyear_alpha_df <- data.frame(multiyear_faith_obs, multiyear_shannons) %>%
  mutate(SampleID = row.names(.))
rm(multiyear_shannons, multiyear_faith_obs)
multiyear_alpha_df <- merge(multiyear_alpha_df, multiyear_meta)
## long format for plotting ease and stat calcs:
multiyear_alpha_df_long <- multiyear_alpha_df %>%
  pivot_longer(cols = c("SR", "PD", "H"), values_to = "Alpha_value", names_to = "Metric")
# 
# 
# ##################
## part 1b: running KW tests for each alpha metric... we're testing for differences between YEAR (given it's within the same sampling Window)
## KW tests:
kw_input_multiyear_SR <- multiyear_alpha_df_long %>% filter(Metric == "SR") %>% mutate(KWgroup = paste(Site, Year, sep="-")) %>% select(KWgroup, Alpha_value)
kw_result_multiyear_SR <- kruskal.test(kw_input_multiyear_SR$Alpha_value ~ kw_input_multiyear_SR$KWgroup)
kw_result_multiyear_SR
  ## not significant for Site-Year group
  ## chi-squared = 6.5427, df = 5, p-value = 0.2569

kw_input_multiyear_H <- multiyear_alpha_df_long %>% filter(Metric == "H") %>% mutate(KWgroup = paste(Site, Year, sep="-")) %>% select(KWgroup, Alpha_value)
kw_result_multiyear_H <- kruskal.test(kw_input_multiyear_H$Alpha_value ~ kw_input_multiyear_H$KWgroup)
kw_result_multiyear_H
  ## not significant
  ## chi-squared = 2.514, df = 5, p-value = 0.7744

kw_input_multiyear_PD <- multiyear_alpha_df_long %>% filter(Metric == "PD") %>% mutate(KWgroup = paste(Site, Year, sep="-")) %>% select(KWgroup, Alpha_value)
kw_result_multiyear_PD <- kruskal.test(kw_input_multiyear_PD$Alpha_value ~ kw_input_multiyear_PD$KWgroup)
kw_result_multiyear_PD
  ## not significant
  ## chi-squared = 3.5567, df = 5, p-value = 0.6148

# ## cleanup
rm(list=ls(pattern="^kw_"))

### quick sanity check - do most of our samples contain at least a handful of OTUs? 
  multiyear_alpha_df_long %>% filter(Metric == "SR" & Alpha_value > 4) %>% nrow()
  ## 68 of 81 (~ 84 %) samples have at least 5 OTUs
  multiyear_alpha_df_long %>% filter(Metric == "SR" & Alpha_value > 9) %>% nrow()
  ## 36 of 81 (~ 47 %) samples have at least 10 OTUs

## add grouper term to simplify plotting
multiyear_alpha_df_long$SiteYear <- paste(multiyear_alpha_df_long$Site, multiyear_alpha_df_long$Year, sep = "\n")
## add a term for facet labels to simplify plotting
multiyear_alpha_df_long <- 
  multiyear_alpha_df_long %>% 
  mutate(MetricStripLabel = case_when(Metric == "SR" ~ "Observed OTUs",
                                      Metric == "H" ~ "Shannon\'s H",
                                      Metric == "PD" ~ "Faith\'s PD"))

## visualize these data in boxplots, grouping data by sampling site, with it's pair of years next to each other
### plotting by individual Metrics, then stitching together into multifaceted plot
multiyear_alpha_df_long$MetricStripLabel <- 
  factor(multiyear_alpha_df_long$MetricStripLabel, levels = c("Observed OTUs", "Shannon\'s H", "Faith\'s PD"))

multiyear_alphaviz_function <- function(metric){
  ggplot(multiyear_alpha_df_long %>% filter(Metric == metric), 
         aes(x=SiteYear, y=Alpha_value)) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.15, alpha=0.5, size=2) +
    facet_wrap(~MetricStripLabel, scales = "free_y", strip.position = "left") +
    labs(x="", y="") +
    theme_classic() +
    theme(strip.placement = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size=14),
          axis.text = element_text(size=12))
}

p_mySR <- multiyear_alphaviz_function("SR")
p_myH <- multiyear_alphaviz_function("H")
p_myPD <- multiyear_alphaviz_function("PD")

ggarrange(p_mySR, p_myH, p_myPD, ncol=1, nrow=3)

ggsave("~/github/nhguano/supplementaryData/figureS6_multiyearAlphaBoxplot_AllMetrics.png", height=21, width=14, units = "cm", dpi=150)
ggsave("~/github/nhguano/supplementaryData/figureS6_multiyearAlphaBoxplot_AllMetrics.pdf", height=21, width=14, units = "cm", dpi=300)

## cleanup:
rm(list=ls(pattern = "^p_my"))
rm(multiyear_alpha_df, multiyear_alpha_df_long)

# ### doesn't appear that there's any story in terms of differences in #OTUs detected per sample among these year/site groupings
# ## but are there any compositional changes?

# ##################
# ## part 2b: community composition
# ## 2bi. run PERMANOVA (via Adonis) for each distance method to evaluate if samples have different centroids among Sites or Years

# ## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method
multiyear_dist_ds <- phyloseq::distance(multiyear_phydat, "bray", binary = TRUE)
multiyear_dist_uu <- phyloseq::distance(multiyear_phydat, "unifrac", weighted=FALSE)

multiyear_adonis_function <- function(distanceData, distanceMetric){
  adonis_tmp <- adonis2(distanceData ~ Site * Year, data = multiyear_meta, )
  data.frame(adonis_tmp) %>% mutate(Metric = distanceMetric,
                                    Class = row.names(.))
}

multiyear_adonis_ds <- multiyear_adonis_function(multiyear_dist_ds, "Dice-Sorensen")
multiyear_adonis_uu <- multiyear_adonis_function(multiyear_dist_uu, "Unifrac Unweighted")

multiyear_adonis_all <- rbind(multiyear_adonis_ds, multiyear_adonis_uu)
rm(multiyear_adonis_ds, multiyear_adonis_uu)
## SITE and YEAR main effects are significant
  ## interaction term for SITE*YEAR is significant too...

multiyear_adonis_all <- 
  multiyear_adonis_all[,c(7,1,2,3,4,5,6)] %>% 
  rename(Fmodel = `F`) %>% 
  mutate(SumOfSqs = round(SumOfSqs, 3),
         R2 = round(R2, 3),
         `Pr..F.` = round(`Pr..F.`, 4),
         Fmodel = round(Fmodel, 3)) %>% 
  rename(Sum.Sq = SumOfSqs, F.value = Fmodel)
write_csv(multiyear_adonis_all,
          file = "~/github/nhguano/supplementaryData/tableS12_multiyear_Adonis_AllMetrics.csv")
 
# 
# ## part 2bii. Test for homogeneity of dispersion
  ## need to test for Site and Year separately; then bind into single data frame (per distance Metric)
multiyear_betadisper_function <- function(distanceData, distanceMetric){
  tmp_disper_list_site <- betadisper(d = distanceData, group =  multiyear_meta$Site, type = c("median"))
  tmp_disper_anova_site <- data.frame(anova(tmp_disper_list_site)) %>% 
    mutate(Class = row.names(.), Metric = distanceMetric, GroupFactor = "Site")
  tmp_disper_list_year <- betadisper(d = distanceData, group =  multiyear_meta$Year, type = c("median"))
  tmp_disper_anova_year <- data.frame(anova(tmp_disper_list_year)) %>% 
    mutate(Class = row.names(.), Metric = distanceMetric, GroupFactor = "Year")
  rbind(tmp_disper_anova_site, tmp_disper_anova_year) %>% 
    select(-`Mean.Sq`) %>% 
    rename(SumOfSqs = "Sum.Sq", Fmodel = `F.value`, pValue = "Pr..F.") %>% 
    relocate(Class, GroupFactor, Df, SumOfSqs, Fmodel, pValue, Metric) %>% 
    mutate(SumOfSqs = round(SumOfSqs, 3), 
           Fmodel = round(Fmodel, 3),
           pValue = round(pValue, 3))
}
 
multiyear_bdisp_ds <- multiyear_betadisper_function(multiyear_dist_ds, "Dice-Sorensen")
multiyear_bdisp_ds
## non significant for either GroupFactor: Site p = 0.641; Year p = 0.382

multiyear_bdisp_uu <- multiyear_betadisper_function(multiyear_dist_uu, "Unifrac Unweighted")
multiyear_bdisp_uu
## non significant for either GroupFactor: Site p = 0.338; Year p = 0.360

multiyear_bdisp_all <- rbind(multiyear_bdisp_ds, multiyear_bdisp_uu)
rm(multiyear_bdisp_ds, multiyear_bdisp_uu)
write_csv(multiyear_bdisp_all,
          file = "~/github/nhguano/supplementaryData/tableS13_multiyear_Betadisper_AllMetrics.csv")

### with group dipsersions being similar, our observed differences for Site and/or Year (and/or interaction)...
### ...are driven by group median diff's between groups!

# ##################
# ## part 2c: ordination of distance data sets
# ### ordinations:
multiyear_ordi_function <- function(distanceValues, distanceMetric){
  tmp_pcoa <- ordinate(multiyear_phydat, method="PCoA", distance = distanceValues)
  tmp_pcoa_list <- plot_ordination(multiyear_phydat, tmp_pcoa)
  tmp_pcoa_df <- tmp_pcoa_list$data %>%
    mutate(SampleID = row.names(.),
           Axis1lab = tmp_pcoa_list$labels$x,
           Axis2lab = tmp_pcoa_list$labels$y,
           Metric = distanceMetric)
  merge(tmp_pcoa_df, multiyear_meta) %>% mutate(Window = as.character(Window))
}

multiyear_orddata_ds <- multiyear_ordi_function(multiyear_dist_ds, "Dice-Sorensen")
multiyear_orddata_uu <- multiyear_ordi_function(multiyear_dist_uu, "Unweighted-UniFrac")

## plot with metadata
#### use shapes to distinguish SITE, and color for YEAR
#### Point shapes to match other common sites in other figures:
  ### HOP == asterisk == 8; MAP == circle == 16; COR == (something not yet used...diamond!) == 18

## plot individually, then stitch together
multiyear_ordplot_function <- function(inputdata){
  inputdata$Year = as.character(inputdata$Year)
  ggplot(data = inputdata,
         aes(x=Axis.1, y=Axis.2, color = Year, shape = Site)) +
    geom_point(size = 3, alpha=0.8) +
    stat_ellipse(aes(group = Year), alpha = 0.5) +
    coord_fixed() +
    labs(x = unique(inputdata$Axis1lab),
         y = unique(inputdata$Axis2lab)) +
    facet_wrap(~Metric) +
    scale_color_manual(values = c('dodgerblue', '380282')) +
    scale_shape_manual(values = c(18, 8, 16)) +
    theme_classic()
}

p_multiyear_ord_ds <- multiyear_ordplot_function(multiyear_orddata_ds)
p_multiyear_ord_ds

p_multiyear_ord_uu <- multiyear_ordplot_function(multiyear_orddata_uu)
p_multiyear_ord_uu

p_multiyear_ord_allMetrics <- ggarrange(p_multiyear_ord_ds, p_multiyear_ord_uu, 
          common.legend = TRUE, nrow=1, ncol=2, align = "h")

### edit save path!~~
ggsave("~/github/nhguano/figures/figure5a_multiyear_ordinations_all.png", width = 20, height = 11, units = 'cm', dpi = 150)
ggsave("~/github/nhguano/figures/figure5a_multiyear_ordinations_all.pdf", width = 20, height = 11, units = 'cm', dpi = 300)

## cleanup
#rm(p_multiyear_ord_ds, p_multiyear_ord_uu)

# 
# ##################
# ## part 2d: pairwise Adonis to determine what, if any, community comopsition groups are different
# ## generate long form of all values for supplementary table
multiyear_pairwise_adonis_func <- function(distancevals, metric){
  pairwise.adonis(distancevals, paste0(multiyear_meta$Site, "-", multiyear_meta$Year), p.adjust.m = "BH") %>%
    mutate(Metric = metric) %>%
    separate(col = pairs, into = c("Group_A", "Group_B"), sep = "vs") %>%
    arrange(p.adjusted)
}

multiyear_pwadonis_ds <- multiyear_pairwise_adonis_func(multiyear_dist_ds, "Dice-Sorensen")
multiyear_pwadonis_uu <- multiyear_pairwise_adonis_func(multiyear_dist_uu, "Unweighted-UniFrac")
multiyear_pwadonis_all <- rbind(multiyear_pwadonis_ds, multiyear_pwadonis_uu)
write_csv(multiyear_pwadonis_all,
          file = "~/github/nhguano/data/text_tables/multiyear_pairwiseAdonis_all.csv")


## plotting these pairwise pvalues:
multiyear_pairwise_adonis_heatmap_function <- function(metric){
  tmp <- multiyear_pwadonis_all %>% 
    filter(Metric == metric) %>% 
    select(Group_A, Group_B, p.adjusted) %>% 
    dcast(Group_A ~ Group_B, value.var = 'p.adjusted')
  tmp_pwadonis_df <- tmp[c(3,1,4,5,2),c(1,2,5,6,3,4)]
  row.names(tmp_pwadonis_df) <- tmp_pwadonis_df$Group_A
  tmp_pwadonis_df$Group_A <- NULL
  tmp_pwadonis_mat <- as.matrix(tmp_pwadonis_df)
  tmp_pwadonis_mat_newRow <- rep(NA, ncol(tmp_pwadonis_mat))
  tmp_pwadonis_mat <- rbind(tmp_pwadonis_mat, tmp_pwadonis_mat_newRow)
  row.names(tmp_pwadonis_mat) <- c(head(str_trim(rownames(tmp_pwadonis_mat), "both"),-1), "HOP-2016")
  tmp_pwadonis_mat_newCol <- rep(NA, nrow(tmp_pwadonis_mat))
  tmp_pwadonis_mat <- cbind(tmp_pwadonis_mat_newCol,tmp_pwadonis_mat)
  colnames(tmp_pwadonis_mat)[1] <- "HOP-2015"
  tmp_pwadonis_mat <- Matrix::forceSymmetric(tmp_pwadonis_mat, uplo = "U")
  row.names(tmp_pwadonis_mat) <- str_trim(rownames(tmp_pwadonis_mat))
  colnames(tmp_pwadonis_mat) <- str_trim(colnames(tmp_pwadonis_mat), "both")
  tmp_pwadonis_df_resolved <- as.data.frame(as.matrix(tmp_pwadonis_mat))
  tmp_pwadonis_df_resolved$Group_A <- row.names(tmp_pwadonis_df_resolved)
  tmp_pwadonis_df_plotdat <- 
    tmp_pwadonis_df_resolved %>% 
    pivot_longer(-Group_A, names_to = "Group_B", values_to = "pvalue", values_drop_na = TRUE) %>% 
    mutate(Metric = metric,
           pvalue = round(pvalue, 3)) %>% 
    mutate(Group_A = gsub("-", "\n", Group_A),
           Group_B = gsub("-", "\n", Group_B))
  
  tmp_pwadonis_df_plotdat <- 
    tmp_pwadonis_df_plotdat %>% 
    mutate(DropStatus = case_when(
      Group_A == "COR\n2016" & Group_B == "COR\n2015" ~ "drop",
      Group_A == "HOP\n2015" & Group_B == "COR\n2015" ~ "drop",
      Group_A == "HOP\n2016" & Group_B == "COR\n2015" ~ "drop",
      Group_A == "MAP\n2015" & Group_B == "COR\n2015" ~ "drop",
      Group_A == "MAP\n2016" & Group_B == "COR\n2015" ~ "drop",
      Group_A == "HOP\n2015" & Group_B == "COR\n2016" ~ "drop",
      Group_A == "HOP\n2016" & Group_B == "COR\n2016" ~ "drop",
      Group_A == "MAP\n2015" & Group_B == "COR\n2016" ~ "drop",
      Group_A == "MAP\n2016" & Group_B == "COR\n2016" ~ "drop",
      Group_A == "HOP\n2016" & Group_B == "HOP\n2015" ~ "drop",
      Group_A == "MAP\n2015" & Group_B == "HOP\n2015" ~ "drop",
      Group_A == "MAP\n2016" & Group_B == "HOP\n2015" ~ "drop",
      Group_A == "MAP\n2015" & Group_B == "HOP\n2016" ~ "drop",
      Group_A == "MAP\n2016" & Group_B == "HOP\n2016" ~ "drop",
      Group_A == "MAP\n2016" & Group_B == "MAP\n2015" ~ "drop",
      TRUE ~ "keep"))
  
  
  tmp_pwadonis_df_plotdat$Group_A <- 
    factor(tmp_pwadonis_df_plotdat$Group_A,
           levels=c("COR\n2015", "COR\n2016", "HOP\n2015", "HOP\n2016", "MAP\n2015", "MAP\n2016"))
  
  tmp_pwadonis_df_plotdat$Group_B <-
    factor(tmp_pwadonis_df_plotdat$Group_B,
           levels=c("COR\n2015", "COR\n2016", "HOP\n2015", "HOP\n2016", "MAP\n2015", "MAP\n2016"))
  
  ggplot(tmp_pwadonis_df_plotdat %>% filter(DropStatus == "keep"),
         aes(Group_B, Group_A, fill=pvalue, label=pvalue)) +
    geom_tile(color="black") +
    geom_text(size = 4.5) +
    coord_fixed() +
    scico::scale_fill_scico(begin = 0, breaks = c(0.05, 0.125, 0.2), limits = c(0,0.25)) +
    facet_wrap(~Metric, nrow=2) +
    labs(x = "", y = "", fill = "BH-adjusted\np-value") +
    theme_classic() +
    theme(legend.position = "right",
          axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          strip.text = element_text(size=14))
}

p_multiyear_pwadonis_dice <- multiyear_pairwise_adonis_heatmap_function("Dice-Sorensen")
p_multiyear_pwadonis_dice

p_multiyear_pwadonis_uuni <- multiyear_pairwise_adonis_heatmap_function("Unweighted-UniFrac")
p_multiyear_pwadonis_uuni

### stitch together into one faceted plot:
ggarrange(p_multiyear_pwadonis_dice, p_multiyear_pwadonis_uuni,
          align="hv", nrow=1, common.legend=TRUE, legend = "right")

#### adjust path
ggsave("~/github/nhguano/figures/unused_multiyear_pairwiseAdonis_heatmap_pvalues_all.png", width = 22, height = 12, units = 'cm', dpi=150)
ggsave("~/github/nhguano/figures/unused_multiyear_pairwiseAdonis_heatmap_pvalues_all.pdf", width = 22, height = 12, units = 'cm', dpi=300)

## cleanup
rm(multiyear_OTUtable_wide_raw, legend_multiyear, tmp_multiyear_plot)
rm(list=ls(pattern="^p_multiyear"))

###### 
## identify indicator taxa (@ genus level not species) to each site+year group, as well as for site (only), and year (only)
### for indicator species analysis, retain just those orders with at least 2% detection in at least one site+year group
multiyear_getTopOrders_function <- function(inputdata){
  ## for all data, across all sampling windows, how many detections total PER YEAR?
  tmp_nGlobalDetections <- inputdata %>% group_by(Site, Year) %>% summarise(nDetections_perSite = n())
  ## retain only those arth Orders where fraction of detections >= 2% in at least 1 site (arth Order can exist only in single site!)
  tmp_nDetections_perSite_perOrder_all <- inputdata %>% 
    group_by(Site, Year, Class, Order) %>% 
    summarise(nDetections_perOrder = n()) %>% 
    merge(., tmp_nGlobalDetections) %>% 
    mutate(pDetections_perSite = nDetections_perOrder / nDetections_perSite)
  tmp_nDetections_perSite_perOrder_all %>% 
    filter(pDetections_perSite >= 0.02) %>% 
    distinct(Order) %>% 
    pull()
}

multiyear_topOrderNames <- multiyear_getTopOrders_function(multiyear_otu_table_long_norm_wMeta_andTaxa)
## 9 orders to consider when filtering from raw data

## get indicator species data by aggregating across shared genus labels
## indicator "taxa" are genus-level not species level!
multiyear_get_indispecData <- function(metadata){
  ## collapse read counts to shared genus labels, then reformat for SRS sampling
  tmp_forSRS_df <- as.data.frame(otu_table_long_filtd_wMeta %>% 
                                   filter(SampleID %in% metadata$SampleID) %>% ## retain same samples as in Physeq object
                                   filter(Order %in% multiyear_topOrderNames) %>% 
                                   mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>% ## rename genus labels to avoid any NA aggregating
                                   group_by(SampleID, Order, Genus) %>% 
                                   summarise(taxaReads = sum(OTUreads)) %>%
                                   mutate(taxaName = paste0(Order, "-", Genus)) %>% 
                                   ungroup() %>% 
                                   select(SampleID, taxaName, taxaReads) %>% 
                                   pivot_wider(names_from = "SampleID", values_from = "taxaReads", values_fill = 0))
  
  ## perform srs sampling    
  row.names(tmp_forSRS_df) <- tmp_forSRS_df$taxaName
  tmp_forSRS_df$taxaName <- NULL
  tmp_Cmin <- min(colSums(tmp_forSRS_df))
  tmp_SRSout_norm <- SRS(tmp_forSRS_df, tmp_Cmin, set_seed = 1)
  row.names(tmp_SRSout_norm) <- row.names(tmp_forSRS_df)
  
  ## flip axes of matrix for indispecies program
  tmp_indi_mat <- t(tmp_SRSout_norm)
}  

multiyear_indiDataTable <- multiyear_get_indispecData(multiyear_meta)

########### indicator species work
#### perform indicator species for three groups: Site (both years), Year (all sites), Site+Year:
doIndiSpec_function <- function(datatable, metadata, metagroup, maxorderval){
  tmp_indval_list <- multipatt(datatable, 
                               metadata %>% 
                                 select(metagroup) %>% pull(), 
                               control = how(nperm=999),
                               max.order = maxorderval)
}


########## first for Sites (both years)
multiyear_indispec_list_site <- doIndiSpec_function(multiyear_indiDataTable, multiyear_meta, "Site", 1)
multiyear_indispec_plotdat_site <- multiyear_indispec_list_site$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(multiyear_indispec_list_site$comb)) %>%
          mutate(index = row.names(.)))

## plot
p_my_indispec_site <- ggplot(multiyear_indispec_plotdat_site,
                             aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(. ~ Order, scales = "free_x", switch = "y", space="free") +
  labs(y="Site", x="", fill = "association\nstatistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(angle = 45, vjust=0.5),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"))

p_my_indispec_site

########## next for Year (all sites)
multiyear_indispec_list_year <- doIndiSpec_function(multiyear_indiDataTable, multiyear_meta, "Year", 1)
multiyear_indispec_plotdat_year <- multiyear_indispec_list_year$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(multiyear_indispec_list_year$comb)) %>%
          mutate(index = row.names(.)))

## plot
p_my_indispec_year <- ggplot(multiyear_indispec_plotdat_year,
                             aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(. ~ Order, scales = "free_x", switch = "y", space="free") +
  labs(y="Year", x="", fill = "association\nstatistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(angle = 45, vjust=0.5),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"))
p_my_indispec_year

########## lastly for particular Site+Year groups
multiyear_indispec_list_siteyear <- doIndiSpec_function(multiyear_indiDataTable, 
                                                        multiyear_meta %>% mutate(SiteYear = paste0(Site, Year)), 
                                                        "SiteYear", 1)
multiyear_indispec_plotdat_siteyear <- multiyear_indispec_list_siteyear$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(multiyear_indispec_list_siteyear$comb)) %>%
          mutate(index = row.names(.)))

## set levels for plot, ordering group that spans windows between distinct windows; 
## plot
p_my_indispec_siteyear <- ggplot(multiyear_indispec_plotdat_siteyear,
                                 aes(y=indexName, x=Genus, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(. ~ Order, scales = "free_x", switch = "y", space="free") +
  labs(y="Site+Year", x="", fill = "association\nstatistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(angle = 45, vjust=0.5),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"))
p_my_indispec_siteyear


### aggregate all these data together to show common genus labels across Site, Year, and Site+Year groups in one plot
multiyear_indispec_plotdat_alldat <- 
  bind_rows(multiyear_indispec_plotdat_siteyear %>% mutate(GroupLabel = "Site+Year"), 
            multiyear_indispec_plotdat_site %>% mutate(GroupLabel = "Site"), 
            multiyear_indispec_plotdat_year %>% mutate(GroupLabel = "Year"))

## order levels of y axis facet label
multiyear_indispec_plotdat_alldat$GroupLabel <-
  factor(multiyear_indispec_plotdat_alldat$GroupLabel, 
         levels = c("Site", "Year", "Site+Year"))

ggplot(multiyear_indispec_plotdat_alldat,
       aes(x=Genus, y=indexName, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid(GroupLabel ~ Order, scales = "free", switch = "y", space="free") +
  labs(y="Site+Year", x="", fill = "association\nstatistic") +
  labs(y="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(angle = 45, vjust=0.5, size=12),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"),
        panel.spacing.y = unit(2, "lines"))


########### allow for combinations of groups in analyses
multiyear_indispec_list_site_comboOkay <- doIndiSpec_function(multiyear_indiDataTable, multiyear_meta, "Site", 3)
multiyear_indispec_plotdat_site_comboOkay <- multiyear_indispec_list_site_comboOkay$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(multiyear_indispec_list_site_comboOkay$comb)) %>%
          mutate(index = row.names(.)))

multiyear_indispec_list_year_comboOkay <- doIndiSpec_function(multiyear_indiDataTable, multiyear_meta, "Year", 2)
multiyear_indispec_plotdat_year_comboOkay <- multiyear_indispec_list_year_comboOkay$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(multiyear_indispec_list_year_comboOkay$comb)) %>%
          mutate(index = row.names(.)))

multiyear_indispec_list_siteyear_comboOkay <- doIndiSpec_function(multiyear_indiDataTable, 
                                                                  multiyear_meta %>% mutate(SiteYear = paste0(Site, Year)), 
                                                                  "SiteYear", 6)
multiyear_indispec_plotdat_siteyear_comboOkay <- multiyear_indispec_list_siteyear_comboOkay$sign %>% 
  mutate(taxaName = row.names(.)) %>% 
  separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
  filter(p.value <= 0.05) %>% 
  select(index, stat, Order, Genus) %>% 
  merge(., 
        data.frame(indexName = colnames(multiyear_indispec_list_siteyear_comboOkay$comb)) %>%
          mutate(index = row.names(.)))

multiyear_indispec_plotdat_siteyear_comboOkay$indexName <- 
  factor(multiyear_indispec_plotdat_siteyear_comboOkay$indexName, 
         levels = c("COR2016", "COR2015+COR2016", "MAP2015", "MAP2016", "COR2016+MAP2015", "COR2015+HOP2015"))

p_multiyear_indispec_plotdat_siteyear_comboOkay <- ggplot(multiyear_indispec_plotdat_siteyear_comboOkay,
       aes(x=Genus, y=indexName, fill=stat)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
                       breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
                       limits = c(0,1)) +
  facet_grid( ~ Order, scales = "free", switch = "y", space="free") +
  labs(y="Site+Year", x="", fill = "association\nstatistic") +
  labs(y="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(angle = 45, vjust=0.5, size=8),
        strip.placement = "outside",
        aspect.ratio = 1,
        strip.background = element_rect(color = NA, fill="gray90"),
        panel.grid.major.y = element_line(color="gray85"),
        panel.grid.major.x = element_line(color="gray85"),
        panel.spacing.y = unit(2, "lines"))

p_multiyear_indispec_plotdat_siteyear_comboOkay
ggsave("~/github/nhguano/figures/figure5b_multiyear_specIndia_allGroups_heatmap.png", height=10, width=15, units="cm", dpi=150)
ggsave("~/github/nhguano/figures/figure5b_multiyear_specIndia_allGroups_heatmap.pdf", height=10, width=15, units="cm", dpi=300)
  ## will add silhouettes and resolve Order names for strip.x panels (in Illustrator) for publication

################## 
## Last section just to compare why our PERMANOVA is significant, but PCoA isn't grouping into obvious cluster...
### and why our indicator species analysis didn't turn up many significant shared genus labels
#### assumption is the difference is driven by the levels of analysis: PCoA and PERMANOVA at OTU-level, indispec at Genus level

### First, how many detections total, at OTU level?
multiyear_totalOTUs_topOrders <- otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  filter(Order %in% multiyear_topOrderNames) %>%
  nrow()
multiyear_totalOTUs_topOrders
## 1032 total detections among all these multiyear samples

#### Second, identify the distinct OTUs in dataset, and tally how many distinct Site+Year groups they occur in
multiyear_getOTUsSumryTable <- 
  otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  filter(Order %in% multiyear_topOrderNames) %>%
  group_by(OTUid) %>% 
  summarise(detections_perGroup = n_distinct(Site, Year)) %>%
  group_by(detections_perGroup) %>% 
  tally()
## 298 OTUs seen in just a single Site+Year
## 7 OTUs seen in all 6 (11 seen in 5 of 6)

#### Now go back and make that original list of how many Site+Year groups a given OTU exists:
multiyear_getOTUsFreqTable <- 
  otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  filter(Order %in% multiyear_topOrderNames) %>%
  group_by(OTUid) %>% 
  summarise(detections_perGroup = n_distinct(Site, Year))

### get a list of OTUs seen in a single Site+Year, and count up how many detections they represent across the entire study
multiyear_OTUlist_1site <-
  multiyear_getOTUsFreqTable %>% 
  filter(detections_perGroup == 1) %>% 
  select(OTUid) %>% pull()

### same list idea, but for OTUs in all 6 sites
multiyear_OTUlist_AllSites <- 
  multiyear_getOTUsFreqTable %>% 
  filter(detections_perGroup == 6) %>% 
  select(OTUid) %>% pull()

#### now count the number of detections for each of those lists:
##### in just 1 site:
otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  filter(OTUid %in% multiyear_OTUlist_1site) %>% 
  nrow()
### represent 337 total detections (or 1032 total, or 32.7%)

##### in all 6 sites:
otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  filter(OTUid %in% multiyear_OTUlist_AllSites) %>% 
  nrow()
### 173 detections (16.8%)

##### so OTUs detected in just a single site account for 32.7% of all detections, across 298 OTUs
##### yet OTUs detected in all six sites account for 16.8% of all detections, among just 7 OTUs!

########################
### could also look at number of detections among OTUs with those shared genus labels?
#### first recognize the genus-labels that are detected in 1 site+year group only
#### then recognize the genus-labels detected in all 6 site+year groups
#### then filter at OTU level by those genus labels (but group Order-Genus "TaxaLabel"), and determine fraction

#### identify the distinct Genus labels in dataset, and tally how many distinct Site+Year groups they occur in
multiyear_getGenusFreqTable <- 
  otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  filter(Order %in% multiyear_topOrderNames) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>% ## rename genus labels to avoid any NA aggregating
  mutate(TaxaLabel = paste(Order, Genus, sep="-")) %>%
  group_by(TaxaLabel) %>% 
  summarise(detections_perGroup = n_distinct(Site, Year))

#### generate 2 lists of genus labels: one for those that exist in all 6 site+year groups, and one for genus labels in only 1 site+year group
### create a list of these Genus names that exist in all 6 site+year groups
multiyear_genusList_allSites <- 
  multiyear_getGenusFreqTable %>% 
  filter(detections_perGroup == 6) %>% 
  select(TaxaLabel) %>% pull()
### create a list of these Genus names that exist in only 1 site+year group
multiyear_genusList_singleSites <- 
  multiyear_getGenusFreqTable %>% 
  filter(detections_perGroup == 1) %>% 
  select(TaxaLabel) %>% pull()

#### calculate fraction of OTU detections for just those particular categores (all 6 groups, or singleton group taxa)
### for all 6 site+year group genus-labels, how many detections exist?
### detections are still at OTU level, just organized by shared genus labels
otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>% ## rename genus labels to avoid any NA aggregating
  mutate(TaxaLabel = paste(Order, Genus, sep="-")) %>%
  filter(TaxaLabel %in% multiyear_genusList_allSites) %>% 
  nrow()
## 646 total detections among genera with shared genus labels (so, 65.3% detections fall into these shared genus labels)

otu_table_long_filtd_wMeta %>%
  filter(SampleID %in% multiyear_meta$SampleID) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>% ## rename genus labels to avoid any NA aggregating
  mutate(TaxaLabel = paste(Order, Genus, sep="-")) %>%
  filter(TaxaLabel %in% multiyear_genusList_singleSites) %>% 
  nrow()
## 121 detections among the genus labels unique to single sites (just 11.7% of shared genus-labels)

#### how many distinct genus labels are there in 1 site only, 2 sites, ... all 6 sites?
multiyear_getGenusFreqTable %>% 
  group_by(detections_perGroup) %>% 
  tally()
### 90 distinct genus labels in singleton site+year groups
### 10 distinct genus labels in ALL site+year groups

############## very last plot looking at number of taxa detected in a particular order for each site+year group... simple barplot...
### for each Site+Year, what fraction of detections occurred for those top 9 orders used in the indispec analysis?
### show a barplot of # detections (of OTUs), per each Order...

multiyear_perSiteYear_detections <- multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  filter(Order %in% multiyear_topOrderNames) %>%
  group_by(Site, Year) %>% 
  summarise(nDetections_perSiteYear=n())

multiyear_perSiteYearOrder_detections <- multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  filter(Order %in% multiyear_topOrderNames) %>%
  group_by(Site, Year, Order) %>% 
  summarise(nDetections_perSiteYearOrder=n()) %>% 
  merge(., multiyear_perSiteYear_detections, by=c("Site", "Year")) %>%
  mutate(pDetections = nDetections_perSiteYearOrder/nDetections_perSiteYear,
         GroupLabeler = paste0(Site, "\n", Year))

### keeping consistent color palette:
multiyear_orderPal <- c('gray25', 'orchid2', 'tan4', '#3E6C54', '#FDAC9F', '#808133', 'cadetblue3', 'gray75', 'turquoise4')

ggplot(multiyear_perSiteYearOrder_detections,
       aes(x=GroupLabeler, y=pDetections, fill=Order)) +
  geom_col() +
  labs(x="", y="fraction of detections") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = multiyear_orderPal) +
  facet_grid(~Site, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave("~/github/nhguano/figures/figure5c_multiyear_detectionsPerOrderBarplot.png", height=10, width=15, units="cm", dpi=150)
ggsave("~/github/nhguano/figures/figure5c_multiyear_detectionsPerOrderBarplot.pdf", height=10, width=15, units="cm", dpi=300)

## cleanup
rm(list=ls(pattern = "^multiyear"))
rm(list=ls(pattern = "^p_m"))



############# UNUSED ####################
# ## for all data, how many detections occur per Genus across all samples in a given sampling window?
#   ## note we're selecting just those top Orders from above
# fox2016_nDetections_perGenus_perWindow <- fox_OTUtable_long_norm_wTaxa_wMeta %>% 
#   filter(Order %in% fox2016_topOrders$Order) %>% 
#   mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family), Genus)) %>% ## this groups all ambiguous Genus labels for a given Order+Family label
#   group_by(Window, Order, Family, Genus) %>% 
#   summarise(nDetections_perGenus = n())
# ## summarise the total number of detections occurring in a sampling window across all samples and all taxa
# fox2016_nDetections_allGenera_perWindow <- fox2016_nDetections_perGenus_perWindow %>% 
#   group_by(Window) %>% 
#   summarise(nDetections_allGenera = sum(nDetections_perGenus))
# ## get the fraction of detections each order is represented as per window
# fox2016_genusDetectionsPerWindow <- merge(fox2016_nDetections_perGenus_perWindow, fox2016_nDetections_allGenera_perWindow, by = "Window") %>%
#   mutate(pDetections_perGenus = nDetections_perGenus / nDetections_allGenera,
#          randomLabeler = paste0("label", row_number())) ## for filtering
# ## get the top taxa per Window, requiring detections in multiple sapmles per Window
# fox2016_TOP_genusDetectionsPerWindow <- fox2016_genusDetectionsPerWindow %>% 
#   group_by(Window) %>% 
#   slice_max(order_by = pDetections_perGenus, n=10) %>% 
#   filter(nDetections_perGenus > 1)
# ## aggregate all the other taxa not included in these TOP taxa to get the bar plot to add to 100%
# fox2016_infrequent_genusDetectionsPerWindow <- fox2016_genusDetectionsPerWindow %>% 
#   filter(!randomLabeler %in% fox2016_TOP_genusDetectionsPerWindow$randomLabeler) %>% 
#   group_by(Window) %>% 
#   summarise(pDetections_perGenus = sum(pDetections_perGenus)) %>% 
#   mutate(Genus = "other taxa",
#          Order = NA)
# ## merge top taxa with aggregated taxa for plot:
# fox2016_genusDetectionsPerWindow_aggregated <- fox2016_TOP_genusDetectionsPerWindow %>% 
#   select(Window, Genus, Order, pDetections_perGenus) %>% 
#   bind_rows(., fox2016_infrequent_genusDetectionsPerWindow)
# ## plot:
# ## 8 Orders represented... keep same color palette per Order
# #"Araneae" "Coleoptera" "Diptera" "Ephemeroptera" "Hymenoptera" "Lepidoptera" "Psocodoea" "Trichoptera"
# fox2016_genusPal <- c("gray25", "orchid2", "tan4", "#3E6C54", 'darkorange', "#808133", "#D49347", "turquoise4")
# #ggplot(data=fox2016_genusDetectionsPerWindow_aggregated,
# ggplot(data=fox2016_TOP_genusDetectionsPerWindow,
#        aes(x=Window, y=pDetections_perGenus)) +
#   geom_point() +
#   geom_text_repel(data=fox2016_TOP_genusDetectionsPerWindow %>% filter(Order == "Diptera" | Order == "Coleoptera"),
#                   aes(x=Window, y=pDetections_perGenus, label=Genus, color=Order),
#                   direction = "y", nudge_x = 0.05,
#                   hjust = 0,
#                   segment.size = 0.2) +
#   geom_text_repel(data=fox2016_TOP_genusDetectionsPerWindow %>% filter(Order != "Diptera" & Order != "Coleoptera"),
#                   aes(x=Window, y=pDetections_perGenus, label=Genus, color=Order),
#                   direction = "y", nudge_x = -0.05,
#                   hjust = 1,
#                   segment.size = 0.2) +
#   scale_color_manual(values = fox2016_genusPal) +
#   theme_classic()

# ## this generates a data table that is easy to share (wide format)
# fox2016_topGenera_byWindows_grid <- fox_OTUtable_long_norm_wTaxa_wMeta %>%
#   mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>%
#   group_by(Window, Order, Genus) %>%
#   summarise(nReads_perWindow_perGenus = sum(SRSreads),
#             nSamples_perWindow_perGenus = n_distinct(SampleID)) %>%
#   merge(., fox2016_sumry_PerWindow) %>%
#   mutate(pReads_perWindow_perGenus = nReads_perWindow_perGenus / nReads_perWindow_global,
#          pSamples_perWindow_perGenus = nSamples_perWindow_perGenus / nSamples_perWindow) %>%
#   mutate(pReads_perWindow_perGenus = round(pReads_perWindow_perGenus, 3),
#          pSamples_perWindow_perGenus = round(pSamples_perWindow_perGenus, 3)) %>%
#   mutate(Window = as.numeric(as.character(Window))) %>%
#   filter(pSamples_perWindow_perGenus >= 0.2 & pReads_perWindow_perGenus >= 0.005) %>%
#   select(Window, Order, Genus, pReads_perWindow_perGenus, nSamples_perWindow_perGenus) %>%
#   mutate(pReads_perWindow_perGenus = pReads_perWindow_perGenus * 100) %>%
#   arrange(Order) %>%
#   pivot_wider(values_from = "pReads_perWindow_perGenus", names_from = "Window", values_fill = 0) %>% 
#   rename()
# # 
# # ## save to disk
# write_csv(fox2016_topGenera_byWindows_grid,
#           path = "~/github/nhguano/fox2016_topGenera_byWindows.csv")

# 
# ### same data as the wide format above, but retained in long form for plotting
# fox2016_topGenera_byWindows <- fox_OTUtable_long_norm_wTaxa_wMeta %>% 
#   mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>% 
#   group_by(Window, Order, Genus) %>% 
#   summarise(nReads_perWindow_perGenus = sum(SRSreads),
#             nSamples_perWindow_perGenus = n_distinct(SampleID)) %>% 
#   merge(., fox2016_sumry_PerWindow) %>% 
#   mutate(pReads_perWindow_perGenus = nReads_perWindow_perGenus / nReads_perWindow_global,
#          pSamples_perWindow_perGenus = nSamples_perWindow_perGenus / nSamples_perWindow) %>% 
#   mutate(pReads_perWindow_perGenus = round(pReads_perWindow_perGenus, 3),
#          pSamples_perWindow_perGenus = round(pSamples_perWindow_perGenus, 3)) %>% 
#   mutate(Window = as.numeric(as.character(Window)),
# ## used to avoid jitter/label mapping issue in plot
#          WindowPlotval = Window + (runif(nrow(.), min=-0.2, max=0.2))) %>%
# ## taxa must be present in at least 20% of samples AND have read abundance > half a percent (per sampling window) 
#   filter(pSamples_perWindow_perGenus >= 0.2 & pReads_perWindow_perGenus >= 0.005) %>% 
#   select(Window, WindowPlotval, Order, Genus, pReads_perWindow_perGenus, pSamples_perWindow_perGenus) %>%
#   arrange(Order)
# 
# ## plotting in terms of samples first
#   ## keep palette consistent with other figure of broad diet pattern:
#   ## 'gray25'  == Araneae; 'darkgoldenrod' == Blattodea; "orchid2" == Coleoptera; "tan4"  == Diptera
#   ## "#3E6C54" == Ephemeroptera; "#FDAC9F" == Hemiptera; "darkorange" == Hymenoptera; "#808133" == Lepidoptera
#   ## "cadetblue3" == Megaloptera; "#D49347" == Psocodoea; "turquoise4" == Trichoptera; 'gray50' == Trombidiformes
# 
# fox2016_topGenera_bySamples_pal <- c("orchid2", "tan4", "#3E6C54", "darkorange", "#808133", "#D49347", "turquoise4")
# 
# ggplot(fox2016_topGenera_byWindows, 
#        aes(x=WindowPlotval, y=pSamples_perWindow_perGenus, color=Order, label=Genus)) +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 2.55, xmax=3.45, ymin=0, ymax=1), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 3.65, xmax=4.45, ymin=0, ymax=1), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 4.65, xmax=5.45, ymin=0, ymax=1), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 5.65, xmax=6.45, ymin=0, ymax=1), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 6.65, xmax=7.45, ymin=0, ymax=1), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 7.65, xmax=8.45, ymin=0, ymax=1), fill="gray85") +
#   geom_point() +
#   scale_color_manual(values = fox2016_topGenera_bySamples_pal) +
#   geom_label_repel(direction = "y", hjust = 0, nudge_x = 0.05, fill = alpha("white", 0.75), show.legend = FALSE) +
#   scale_x_continuous(breaks = seq(3,8)) +
#   labs(x = "\nWindow", y = "Samples detected\n", color = "Arthropod\norder") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15),
#         legend.title = element_text(size = 14), legend.text = element_text(size = 12))
# 
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_dotplot_arthGenus_perSAMPLE_byWindows.png", width = 35, height = 20, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_dotplot_arthGenus_perSAMPLE_byWindows.svg", width = 35, height = 20, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_dotplot_arthGenus_perSAMPLE_byWindows.pdf", width = 35, height = 20, units = 'cm')
# 
# ## and plot in terms of read abundance:
# ggplot(fox2016_topGenera_byWindows, 
#        aes(x=WindowPlotval, y=pReads_perWindow_perGenus, color=Order, label=Genus)) +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 2.55, xmax=3.45, ymin=0, ymax=max(pReads_perWindow_perGenus)), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 3.65, xmax=4.45, ymin=0, ymax=max(pReads_perWindow_perGenus)), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 4.65, xmax=5.45, ymin=0, ymax=max(pReads_perWindow_perGenus)), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 5.65, xmax=6.45, ymin=0, ymax=max(pReads_perWindow_perGenus)), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 6.65, xmax=7.45, ymin=0, ymax=max(pReads_perWindow_perGenus)), fill="gray85") +
#   geom_rect(inherit.aes = FALSE, aes(xmin = 7.65, xmax=8.45, ymin=0, ymax=max(pReads_perWindow_perGenus)), fill="gray85") +
#   geom_point() +
#   scale_color_manual(values = fox2016_topGenera_bySamples_pal) +
#   geom_label_repel(force = 1, size = 3,
#     hjust = 0, nudge_x = 0.05, fill = alpha("white", 0.75), show.legend = FALSE) +
#   scale_x_continuous(breaks = seq(3,8)) +
#   labs(x = "\nWindow", y = "Samples detected\n", color = "Arthropod\norder") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15),
#         legend.title = element_text(size = 14), legend.text = element_text(size = 12))
# 
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_dotplot_arthGenus_perREADS_byWindows.png", width = 35, height = 20, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_dotplot_arthGenus_perREADS_byWindows.svg", width = 35, height = 20, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_dotplot_arthGenus_perREADS_byWindows.pdf", width = 35, height = 20, units = 'cm')
#   ## MIGHT be worth moving around some labels in Illustrator to make the plot cleaner
# 

## create heatmap matrix of these values to report both adjusted and unadjusted pvalues
# fox2016_pwadonis_all$Metric <- factor(fox2016_pwadonis_all$Metric, levels = c(
#   'Dice-Sorensen', 'Unweighted-Unifrac'))
# 
# p_fox2016pwadonis_adj <- ggplot(fox2016_pwadonis_all %>% mutate(p.adjusted = round(p.adjusted, 2)), 
#                                 aes(x=windowA,y=windowB,fill=p.adjusted,label = p.adjusted)) +
#   geom_tile(color = "black") +
#   geom_text(size=2.5) +
#   coord_fixed() +
#   #scale_fill_viridis_c(option = "viridis", direction = -1, alpha = 0.5) +
#   scico::scale_fill_scico(palette = "bilbao", end = 0.65) +
#   facet_wrap(~Metric) +
#   labs(x="\nsampling window", y="sampling window\n", fill="p.value") +
#   theme_classic()
# 
# p_fox2016pwadonis_unadj <- ggplot(fox2016_pwadonis_all %>% mutate(p.value = round(p.value, 2)), 
#                                   aes(x=windowA,y=windowB,fill=p.value, label = p.value)) +
#   geom_tile(color = "black") +
#   geom_text(size=2.5) +
#   coord_fixed() +
#   #scale_fill_viridis_c(option = "viridis", direction = -1, alpha = 0.5) +
#   scico::scale_fill_scico(palette = "bilbao", end = 0.65) +
#   facet_wrap(~Metric) +
#   labs(x="\nsampling window", y="sampling window\n", fill="p.value") +
#   theme_classic()
# 
# ggarrange(p_fox2016pwadonis_adj, p_fox2016pwadonis_unadj, 
#           labels = c("A", "B"), nrow=2, ncol=1, common.legend = TRUE)
# 
# ggsave("~/github/nhguano/supplementaryData/figureS3_fox2016_pairwiseAdonis_heatmap_pvalues_all.png", width = 16, height = 16, units = 'cm', dpi=150)
# ggsave("~/github/nhguano/supplementaryData/figureS3_fox2016_pairwiseAdonis_heatmap_pvalues_all.pdf", width = 16, height = 16, units = 'cm', dpi=300)

##################
## part 2e: pairwise distance boxplots
## ended up running a separate series of calculations using QIIME in a shell script...
## available here: https://raw.githubusercontent.com/devonorourke/nhguano/master/scripts/shell_scripts/get_beta_group_significance_windows2016data.sh
## the result is a pair of files (win456 and win567) containing all pairwise comparisons of distances for each of the four distance metrics

## note that the original tree .nwk file (imported as 'tree' in R) was already available for use in the shell script above
## First, this is what was exported for use in that shell script:
#### OTUtable for win456 data
# win456_otutable <- as.data.frame(as.matrix(otu_table(win456_phydat)))
# win456_otutable$`OTU id` <- row.names(win456_otutable)
# win456_otutable <- win456_otutable %>% 
#   select(`OTU id`, everything())
# write_tsv(win456_otutable, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win456_OTUtable.tsv")
# rm(win456_otutable)
# #### OTUtable for win567 data
# win567_otutable <- as.data.frame(as.matrix(otu_table(win567_phydat)))
# win567_otutable$`OTU id` <- row.names(win567_otutable)
# win567_otutable <- win567_otutable %>% 
#   select(`OTU id`, everything())
# write_tsv(win567_otutable, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win567_OTUtable.tsv")
# rm(win567_otutable)
# ## Metadata for win456 and win567 datasets
# win456_meta_alt <- win456_meta %>% 
#   mutate(Window = paste0("window_", Window))
# write_tsv(win456_meta_alt, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win456meta_forQIIME.txt")
# win567_meta_alt <- win567_meta %>% 
#   mutate(Window = paste0("window_", Window))
# write_tsv(win567_meta_alt, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win567meta_forQIIME.txt")
# 
# ## cleanup:
# rm(list = ls(pattern = "^win"))
# 
# ##################
# ## imports data from output of `.sh` script...
# ## imports metadata and merges with distance data
# ## calculates the bottom 5% distance value for a given Metric among within site points 
# 
# betagroupsig_plotdat_function <- function(dist_url, meta_url, metric){
#   tmp_betagroupsig_allmetrics <- 
#     read_delim(delim = "\t", dist_url) %>% 
#     mutate(MetricLabel = case_when(Metric == "dice" ~ 'Dice-Sorensen',
#                                    Metric == "bray" ~ 'Bray-Curtis',
#                                    Metric == "uuni" ~ 'Unweighted-UniFrac',
#                                    Metric == "wuni" ~ 'Weighted-UniFrac'),
#            Group1 = str_replace(Group1, "window_", "Window "),
#            Group2 = str_replace(Group2, "window_", "Window ")) %>% 
#     filter(FactorGroup == "site" & Metric == metric)
#   
#   win567_meta_tmp <- read_tsv(meta_url) %>% select(SampleID, Window) %>% mutate(Window = str_replace(Window, "window_", "Window "))
#   
#   tmp_betagroupsig_allmetrics <- 
#     merge(tmp_betagroupsig_allmetrics, win567_meta_tmp, by.x = "SubjectID1", by.y = "SampleID", all.x = TRUE) %>% 
#     rename(WindowSubject1 = "Window")
#   
#   tmp_betagroupsig_allmetrics <- 
#     merge(tmp_betagroupsig_allmetrics, win567_meta_tmp, by.x = "SubjectID2", by.y = "SampleID", all.x = TRUE) %>% 
#     rename(WindowSubject2 = "Window")
#   
#   tmp_betagroupsig_allmetrics <- 
#     tmp_betagroupsig_allmetrics %>% 
#     mutate(WindowCompareTest = (WindowSubject1 == WindowSubject2),
#            WindowLabel = case_when(WindowCompareTest == TRUE ~ paste0('within ',WindowSubject1),
#                                    WindowCompareTest == FALSE ~ paste0(WindowSubject1, ' and ', WindowSubject2)))
#   
#   tmp_betagroupsig_allmetrics <- tmp_betagroupsig_allmetrics %>% 
#     #group_by(grp = paste(pmax(SubjectID1, SubjectID2), pmin(SubjectID1, SubjectID2), sep = "_")) %>% 
#     #slice(1) %>% 
#     #ungroup() %>% 
#     #select(-grp) %>% 
#     mutate(WindowLabel = case_when(WindowLabel == "Window 6 and Window 4" ~ "Window 4 and Window 6",
#                                    WindowLabel == "Window 6 and Window 5" ~ "Window 5 and Window 6",
#                                    WindowLabel == "Window 5 and Window 4" ~ "Window 4 and Window 5",
#                                    WindowLabel == "Window 7 and Window 5" ~ "Window 5 and Window 7",
#                                    WindowLabel == "Window 7 and Window 6" ~ "Window 6 and Window 7",
#                                    TRUE ~ as.character(.$WindowLabel)),
#            SiteLabel = paste0(Group1, "_", Group2))
#   
#   dist_sumry <- tmp_betagroupsig_allmetrics %>% 
#     group_by(Metric) %>% 
#     summarise(DistanceBar = quantile(Distance, 0.05, q=0.05))
#   
#   merge(tmp_betagroupsig_allmetrics, dist_sumry)
# }
# 
# win456_betagroupsig_plotdat_dice <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
#   'dice')
# 
# win456_betagroupsig_plotdat_bray <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
#   'bray')
# 
# win456_betagroupsig_plotdat_uuni <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
#   'uuni')
# 
# win456_betagroupsig_plotdat_wuni <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
#   'wuni')
# 
# 
# win567_betagroupsig_plotdat_dice <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
#   'dice')
# 
# win567_betagroupsig_plotdat_bray <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
#   'bray')
# 
# win567_betagroupsig_plotdat_uuni <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
#   'uuni')
# 
# win567_betagroupsig_plotdat_wuni <- betagroupsig_plotdat_function(
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
#   'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
#   'wuni')
# 
# 
# ## generate data for plotting function that will plot each Metric (for each windowgroup) separately
# ## use value in 'DistanceBar' to highlight the bottom 5% of distance values per Metric...
# ## helps illustrate which Sites are more similar to each other within each Window
# 
# ## droplist for the grid to focus on just one of the two reversible pairs
# dropsitelabels_win567 <- c("CNB_CHI", "FOX_CHI", "MTV_CHI", "FOX_CNB", "MTV_CNB", "MTV_FOX")
# 
# dropsitelabels_win456 <- c("EPS_CNA", "FOX_CNA", "HOL_CNA", "HOP_CNA", "MAP_CNA", "PEN_CNA",
#                            "FOX_EPS", "HOL_EPS", "HOP_EPS", "MAP_EPS", "PEN_EPS",
#                            "HOL_FOX", "HOP_FOX", "MAP_FOX", "PEN_FOX",
#                            "HOP_HOL", "MAP_HOL", "PEN_HOL",
#                            "MAP_HOP", "PEN_HOP",
#                            "PEN_MAP")
# 
# betagroupsig_plotfunction_win567 <- function(plotdat, dropsitelabels){
#   
#   plotdat$WindowLabel <-
#     factor(plotdat$WindowLabel, 
#            levels = c("within Window 5", "within Window 6", "within Window 7",
#                       "Window 5 and Window 6", "Window 6 and Window 7", "Window 5 and Window 7"))
#   ggplot(data = plotdat %>% filter(!SiteLabel %in% dropsitelabels),
#          aes(x=Group1, 
#              y=Distance,
#              color=WindowLabel)) +
#     geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1),
#                alpha = 0.5) +
#     geom_boxplot(position = position_dodge(1), outlier.shape=NA,
#                  fill="white", alpha=0.5) +
#     #geom_hline(yintercept = unique(plotdat$DistanceBar), color="gray20", alpha = 0.4) +
#     scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.75)) +
#     scale_color_viridis_d(option = 'magma', begin = 0, end = 0.9) +
#     labs(x="", y="Distance\n", fill="Sample pair") +
#     facet_grid(Group2 ~ Group1, scales = "free_x", space = "free_x") +
#     theme_classic() +
#     theme(panel.spacing.y = unit(2, "lines"))
# }
# 
# # p_betasig_win567_dice <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_dice, dropsitelabels_win567)
# # p_betasig_win567_dice
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_dice.png", height = 15, width = 22, units = 'cm', dpi=150)
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_dice.pdf", height = 15, width = 22, units = 'cm', dpi=300)
# # 
# # p_betasig_win567_bray <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_bray, dropsitelabels_win567)
# # p_betasig_win567_bray
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_bray.png", height = 15, width = 22, units = 'cm', dpi=150)
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_bray.pdf", height = 15, width = 22, units = 'cm', dpi=300)
# # 
# # p_betasig_win567_uuni <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_uuni, dropsitelabels_win567)
# # p_betasig_win567_uuni
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_uuni.png", height = 15, width = 20, units = 'cm', dpi=150)
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_uuni.pdf", height = 15, width = 20, units = 'cm', dpi=300)
# # 
# # p_betasig_win567_wuni <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_wuni, dropsitelabels_win567)
# # p_betasig_win567_wuni
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_wuni.png", height = 15, width = 20, units = 'cm', dpi=150)
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_wuni.pdf", height = 15, width = 20, units = 'cm', dpi=300)
# # 
# # ggpubr::ggarrange(p_betasig_win567_dice, p_betasig_win567_bray, p_betasig_win567_uuni, p_betasig_win567_wuni,
# #                   common.legend = TRUE, align = "hv", labels = c("A", "B", "C", "D"))
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_allMetrics.png", height = 20, width = 20, units = 'cm', dpi=150)
# # ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_allMetrics.pdf", height = 20, width = 20, units = 'cm', dpi=300)
# 
# rm(list = ls(pattern = "^p_betasig_win567"))
# rm(list = ls(pattern = "^win567"))
# 
# 
# betagroupsig_plotfunction_win456 <- function(plotdat, dropsitelabels){
#   
#   plotdat$WindowLabel <-
#     factor(plotdat$WindowLabel, 
#            levels = c("within Window 4", "within Window 5", "within Window 6",
#                       "Window 4 and Window 5", "Window 5 and Window 6", "Window 4 and Window 6"))
#   ggplot(data = plotdat %>% filter(!SiteLabel %in% dropsitelabels),
#          aes(x=Group1, 
#              y=Distance,
#              color=WindowLabel)) +
#     geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1),
#                alpha = 0.5) +
#     geom_boxplot(position = position_dodge(1), outlier.shape=NA,
#                  fill="white", alpha=0.5) +
#     #geom_hline(yintercept = unique(plotdat$DistanceBar), color="gray20", alpha = 0.4) +
#     scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.75)) +
#     scale_color_viridis_d(option = 'magma', begin = 0, end = 0.9) +
#     labs(x="", y="Distance\n", fill="Sample pair") +
#     facet_grid(Group2 ~ Group1, scales = "free_x", space = "free_x") +
#     theme_classic() +
#     theme(panel.spacing.y = unit(2, "lines"))
# }

# p_betasig_win456_dice <- betagroupsig_plotfunction_win456(win456_betagroupsig_plotdat_dice, dropsitelabels_win456)
# p_betasig_win456_dice
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_dice.png", height = 20, width = 25, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_dice.pdf", height = 20, width = 25, units = 'cm', dpi=300)
# 
# p_betasig_win456_bray <- betagroupsig_plotfunction_win456(win456_betagroupsig_plotdat_bray, dropsitelabels_win456)
# p_betasig_win456_bray
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_bray.png", height = 20, width = 25, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_bray.pdf", height = 20, width = 25, units = 'cm', dpi=300)
# 
# p_betasig_win456_uuni <- betagroupsig_plotfunction_win456(win456_betagroupsig_plotdat_uuni, dropsitelabels_win456)
# p_betasig_win456_uuni
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_uuni.png", height = 20, width = 25, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_uuni.pdf", height = 20, width = 25, units = 'cm', dpi=300)
# 
# p_betasig_win456_wuni <- betagroupsig_plotfunction_win456(win456_betagroupsig_plotdat_wuni, dropsitelabels_win456)
# p_betasig_win456_wuni
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_wuni.png", height = 20, width = 25, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_wuni.pdf", height = 20, width = 25, units = 'cm', dpi=300)
# 
# ggpubr::ggarrange(p_betasig_win456_dice, p_betasig_win456_bray, p_betasig_win456_uuni, p_betasig_win456_wuni,
#                   common.legend = TRUE, align = "hv", labels = c("A", "B", "C", "D"))
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_allMetrics.png", height = 30, width = 20, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win456_betagroupsig_allMetrics.pdf", height = 30, width = 20, units = 'cm', dpi=300)
# 
# rm(list = ls(pattern = "^p_betasig_win456"))
# rm(list = ls(pattern = "^win456"))
# 
# whitenerlist = c("Araneae", "Diptera", "Ephemeroptera", "Trombidiformes")
# 
# p_genera456dots <- ggplot() +
#   geom_point(data=win456_GenusDetect_plotdat,
#              aes(x=Window, y=pSamples_perSiteWindow_perGenus),
#              show.legend = TRUE) +
#   geom_label_repel(data=win456_GenusDetect_plotdat, 
#                    aes(x=Window, y=pSamples_perSiteWindow_perGenus, label=Genus, 
#                        fill=Order, fontface = "bold",
#                        color = ifelse(Order %in% whitenerlist, "gray90", "gray20")),
#                    direction = "y", nudge_x = -0.05,
#                    hjust = 1,
#                    segment.size = 0.2,
#                    segment.colour = "black",
#                    size = 2,
#                    min.segment.length = 0) +
#   scale_fill_manual(values = win456_generaPal) +
#   scale_color_manual(values = c("black", "white")) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(0,1)) +
#   scale_x_continuous(breaks = c(4, 5, 6), limits = c(3.3, 6.2)) +
#   facet_wrap(~Site, scales = "free_x", ncol=2) +
#   labs(x = "\n Window", y="fraction samples with taxa", color = "Order") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.position = "top") +
#   guides(color = FALSE, 
#          fill = guide_legend(override.aes = aes(color = NA)))

## for win567, plotting all as single row (4 wide) to match spacing with win456
# p_genera567dots <- ggplot() +
#   geom_point(data=win567_GenusDetect_plotdat,
#              aes(x=Window, y=pSamples_perSiteWindow_perGenus),
#              show.legend = TRUE) +
#   geom_label_repel(data=win567_GenusDetect_plotdat, 
#                    aes(x=Window, y=pSamples_perSiteWindow_perGenus, label=Genus, 
#                        color = ifelse(Order %in% whitenerlist, "gray90", "gray20"),
#                        fill=Order, fontface = "bold"),
#                    direction = "y", nudge_x = -0.05,
#                    hjust = 1,
#                    segment.size = 0.2,
#                    segment.colour = "black",
#                    size = 2,
#                    force
#                    min.segment.length = 0) +
#   scale_fill_manual(values = win567_generaPal) +
#   scale_color_manual(values = c("black", "white")) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(0,1)) +
#   scale_x_continuous(breaks = c(5, 6, 7), limits = c(4.3, 7.2)) +
#   facet_wrap(~Site, scales = "free_x", nrow=2, ncol=2) +
#   labs(x = "\n Window", y="fraction samples with taxa", color = "Order") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.position = "top") +
#   guides(color = FALSE, 
#          fill = guide_legend(override.aes = aes(color = NA)))
# 
# ## plot to stitch together, to get the heights to be equivalent per facet
# windows2016_generaLayout <- "
# AAB
# "
# p_5c <- p_genera456dots + p_genera567dots + plot_layout(design = windows2016_generaLayout)
# p_5c


# ## part 0b: rarefy table using SRS for samples in each group
# ## then generate phyloseq object by importing with metadata and tree file to generate phyloseq object
# win567_meta <- otu_table_long_filtd_wMeta %>% 
#   filter(SampleID %in% win567_samplenames) %>% 
#   select(SampleID, Window, Site, newDate) %>% 
#   mutate(Window = as.factor(Window)) %>% 
#   distinct()
# ## 166 samples total
# 
# ## create phyloseq object:
# win567_phydat <- windows2016_getPhyloseq_function(win567_meta)
# 
# ##################
# ## part 1a: alpha diversity: species richness, Shannon's entropy, Faith's PD...
# win567_alpha_df <- windows2016_getAlphaVals_function(win567_phydat, win567_meta)
# 
# ##################
# ## part 1b: applying KW test for global differences in distinct diversity metrics for a given Site+Window group
# ## ...then post hoc Dunn's test for pairwise diffs in alpha vals by Site+Group
# win567_KW_result_SR <- windows2016_getAlphaStats_function(win567_alpha_df, "SR")
# win567_KW_result_SR
# ##significant; chi-squared = 46.508, df = 11, p-value = 0
# win567_KW_result_H <- windows2016_getAlphaStats_function(win567_alpha_df, "H")
# win567_KW_result_H
# ## NOT! significant; chi-squared = 19.035, df = 11, p-value = 0.0605
# win567_KW_result_PD <- windows2016_getAlphaStats_function(win567_alpha_df, "PD")
# win567_KW_result_PD
# ##significant; chi-squared = 45.168, df = 11, p-value = 0
# 
# ## get Dunn's pairwise vals for each group (Site+Window):
# win567_dunn_SR <- dunns_pvals_windows2016(win567_alpha_df, "SR", "567")
# ## chichester seems very different here - 14 of 16 sig sites have to do with CHI 5|6|7...
# win567_dunn_H <- dunns_pvals_windows2016(win567_alpha_df, "H", "567")
# ## 5CHI again. only 2 with pval <= 0.05. CHI5:CHI6, and CHI5:CHI7! ... so why does CHI have such a huge drop at from 5 to windows 6|7?
# win567_dunn_PD <- dunns_pvals_windows2016(win567_alpha_df, "PD", "567")
# ## CHI site involved in 10 of 11 sig pairwise diffs... 
# #### the story of win567 is CHICHESTER only...
# 
# ## plot separate windowgroups, then stitch together
# p_dunn_567 <- ggplot(windows2016_dunn_all %>% filter(Windowgroup == "window567"), 
#                      aes(x=windowB, y=windowA, fill=BHadjusted_pvalue, label=sigVal)) +
#   geom_tile(color="black", size=0.25) +
#   geom_text() +
#   facet_grid(Windowgroup ~ Metric, space = "free") +
#   scico::scale_fill_scico(palette = "bilbao") +
#   coord_fixed() +
#   labs(x="", y="", fill = "BH-corrected\np-value") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle=90, hjust=1),
#         strip.text.y = element_blank())
# 
# ## plot these groups as boxplot too, and generate alphabetical labels to show pairwise differences...
# ### first, generate the labels
# win567_dunnLetters_SR <- windows2016_DunnLetters_function(win567_alpha_df, "SR", "567")
# win567_dunnLetters_H <- windows2016_DunnLetters_function(win567_alpha_df, "H", "567")
# win567_dunnLetters_PD <- windows2016_DunnLetters_function(win567_alpha_df, "PD", "567")
# windows2016_dunnLetters_all567 <- rbind(win567_dunnLetters_SR, win567_dunnLetters_H, win567_dunnLetters_PD)
# 
# ## merge letter data (per group) with individual data (per sample) by site and window (the group!):
# win567_alpha_plotdat <- 
#   merge(win567_alpha_df, windows2016_dunnLetters_all567, by=c("Site", "Window", "Metric")) %>% 
#   mutate(Grouper = paste0(Window, Site),
#          FacePlotLabel = case_when(
#            Metric == "SR" ~ "Observed OTUs", 
#            Metric == "H" ~ "Shannon\'s H", 
#            Metric == "PD" ~ "Faith\'s PD"))
# 
# win567_alpha_plotdat <- win567_alpha_plotdat %>% 
#   group_by(Metric) %>% 
#   summarise(LabelAlpha_value = max(Alpha_value) * 1.1) %>% 
#   merge(., win567_alpha_plotdat, by="Metric")
# 
# ## set levels
# win567_alpha_plotdat$FacePlotLabel <- factor(win567_alpha_plotdat$FacePlotLabel, levels=c("Observed OTUs", "Shannon\'s H", "Faith\'s PD"))
# 
# ## plot
# ggplot() +
#   geom_boxplot(data = win567_alpha_plotdat,
#                aes(x=Window, y=Alpha_value),
#                outlier.color = NA) +
#   geom_jitter(data = win567_alpha_plotdat,
#               aes(x=Window, y=Alpha_value),
#               alpha=0.65, width = 0.1) +
#   geom_text(data=win567_alpha_plotdat, 
#             aes(x=Window, y=LabelAlpha_value, label=Letters),
#             #angle=45, hjust=1,
#             size = 2) +
#   facet_grid(FacePlotLabel ~ Site, scales = "free", space = "free_x", switch = "y") +
#   labs(x="", y="") +
#   theme_classic() +
#   theme(panel.spacing.x = unit(0.5, "lines"), panel.spacing.y = unit(3, "lines"),
#         strip.placement.y = "outside",
#         #strip.text.x = element_blank(), ## remove this line if you want to show per site facets at top of plot
#         strip.background.y = element_blank())
# 
# 
# 
# # plot_grid(p_boxplot_win456,p_boxplot_win567,
# #           labels = c("A", "B"),
# #           nrow=1, ncol=2, 
# #           rel_widths = c(7/11, 4/11))
# 
# 
# ##################
# ## part 2a: community composition
# ## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method, per window (456 or 567)
# win567_dist_ds <- phyloseq::distance(win567_phydat, "bray", binary = TRUE)
# win567_dist_uu <- phyloseq::distance(win567_phydat, "unifrac", weighted=FALSE)
# 
# ##################
# ## part 2b: community composition
# ## 2bi. run PERMANOVA (via Adonis) for each distance method to evaluate if samples within each Site-Window group have different centroids
# win567_adonis_ds <- windows2016_adonis_function(win567_dist_ds, "Dice-Sorensen", win567_meta, "windows567")
# win567_adonis_uu <- windows2016_adonis_function(win567_dist_uu, "Unifrac-Unweighted", win567_meta, "windows567")
# win567_adonis_all <- rbind(win567_adonis_ds, win567_adonis_uu)
# rm(win567_adonis_ds, win567_adonis_uu)
# windows2016_adonis_allMetrics_allWindows <- rbind(win456_adonis_all, win567_adonis_all)
# 
# ##################
# ## part 2bii. Test for homogeneity of dispersion
# ## have to test for Site and Window main effects separately for each WindowGroup&DistanceMetric combo...
# win567_bdisp_ds <- window2016_betadisper_function(win567_dist_ds, "Dice-Sorensen", win567_meta, "windows567")
# win567_bdisp_uu <- window2016_betadisper_function(win567_dist_uu, "Unifrac-Unweighted", win567_meta, "windows567")
# 
# 
# ##################
# ## part 2c: ordinations
# win567_orddata_ds <- windows2016_ordi_function(win567_dist_ds, "Dice-Sorensen", win567_phydat, "win567", win567_meta)
# win567_orddata_uu <- windows2016_ordi_function(win567_dist_uu, "Unweighted-UniFrac", win567_phydat, "win567", win567_meta)
# 
# ## for all plot types, set up these color/shape parameters
# ## palette to match windows from other plots!
# win567_pal <- rev(c("#267b97", "#557d3f", "#B28C32"))
# 
# ## need to define shape point types; ensure overlapping Sites have same point value
# win567_shapes <- c(0, 1, 17, 4)
# ## note the orders here: the third values need to be the same so that "FOX" site is common in both 'win456' and 'win567' plots
# p_win567_ord_ds <- window2016_plotfunction(win567_orddata_ds, win567_shapes, win567_pal)
# 
# ### 1) single plots per WindowGroup & Distance estimate
# p_win567_ord_ds
# ggsave(filename = "~/github/nhguano/figures/figure4aiii_pcoa_win567_ds.png", width = 17, height = 15, unit = "cm", dpi=150)
# ggsave(filename = "~/github/nhguano/figures/figure4aiii_pcoa_win567_ds.pdf", width = 17, height = 15, unit = "cm", dpi=300)
# 
# p_win567_ord_uu <- window2016_plotfunction(win567_orddata_uu, win567_shapes, win567_pal)
# p_win567_ord_uu
# ggsave(filename = "~/github/nhguano/figures/figure4aiv_pcoa_win567_uu.png", width = 17, height = 15, unit = "cm", dpi=150)
# ggsave(filename = "~/github/nhguano/figures/figure4aiv_pcoa_win567_uu.pdf", width = 17, height = 15, unit = "cm", dpi=300)
# 
# 
# ### 2) combined plots with different distance estimates for data with shared WindowGroups
# tmp_win567_pcoa_legend <- cowplot::get_legend(p_win567_ord_ds)
# 
# ##################
# ## part 2d: pairwise adonis
# win567_pwadonis_ds <- window2016_pairwise_adonis_func(win567_dist_ds, "Dice-Sorensen", "windows567", win567_meta)
# win567_pwadonis_uu <- window2016_pairwise_adonis_func(win567_dist_uu, "Unweighted-UniFrac", "windows567", win567_meta)
# 
# ##################
# ## part 2e: what taxa are different between these Site+Window groups?
# ## identify the proportion of detections BY arthropod ORDER
# win567_otutable_long_wtaxa_wmeta <- get_windows2016_srsdata_wtaxa_function(win567_phydat, win567_meta)
# 
# ####### first plot is at Order level...
# win567_OrderDetect_plotdat <- windows2016_getTopOrderData_function(win567_otutable_long_wtaxa_wmeta)
# 
# ## keep order palete colors consistent between plots!
# win567_OrderDetect_plotdat %>% distinct(Order) %>% arrange(Order) %>% pull()
# ##8+1 orders included are: Coleoptera" "Diptera" "Ephemeroptera" "Hemiptera" "Hymenoptera" "Lepidoptera" "Megaloptera" "other" "Trichoptera"
# win567_orderPal <- c("orchid2", "tan4", "#3E6C54", "#FDAC9F", "darkorange", "#808133", 'cadetblue3', '#e9ee7a', "turquoise4")
# 
# # plot, faceting site groups together and ordered by sampling window within facet
# ggplot(win567_OrderDetect_plotdat,
#        aes(x=Window, y=pDetections_perSiteWindow, fill=Order)) + 
#   geom_col() + 
#   scale_fill_manual(values = win567_orderPal) +
#   facet_grid(~Site, space="free_x", scales = "free_x") +
#   labs(x="\n sampling Window", y = "fraction of detections\n") +
#   theme_classic() +
#   theme(legend.position = "top") +
#   guides(fill = guide_legend(nrow=2))
# 
# 
# # ## plot to stitch together, same width in export as 5c!
# # windows2016_orderLayout <- "
# # AAAAAAABBBB"
# # p_5b <- p_order456bar + p_order567bar + plot_layout(design = windows2016_orderLayout)
# # p_5b
# 
# 
# ####### second plot is at Genus level...
# ## find taxa frequently detected with shared genus labels
# windows2016_getTopGeneraData_function <- function(otutable_long){
#   ## gather per-SiteWindow sample size and total read counts for filtering in next step
#   tmp_sumry_PerSiteWindow <- otutable_long %>%
#     group_by(Site, Window) %>% 
#     summarise(nSamples_perSiteWindow = n_distinct(SampleID),
#               nReads_perSiteWindow_global = sum(SRSreads))
#   ## group by shared genus label, requiring only taxa in at least 20% of samples with a minimum of 0.5% of reads
#   tmp_topGenus_perSiteWindow <- otutable_long %>%
#     mutate(Genus = ifelse(is.na(Genus), paste0("f. ", Family, " ", OTUalias), Genus)) %>%
#     group_by(Site, Window, Order, Genus) %>%
#     summarise(nReads_perSiteWindow_perGenus = sum(SRSreads),
#               nSamples_perSiteWindow_perGenus = n_distinct(SampleID)) %>%
#     merge(., tmp_sumry_PerSiteWindow, by = c("Site", "Window")) %>%
#     mutate(pReads_perSiteWindow_perGenus = nReads_perSiteWindow_perGenus / nReads_perSiteWindow_global,
#            pSamples_perSiteWindow_perGenus = nSamples_perSiteWindow_perGenus / nSamples_perSiteWindow) %>%
#     mutate(Window = as.numeric(as.character(Window))) %>%
#     filter(pSamples_perSiteWindow_perGenus >= 0.2 & pReads_perSiteWindow_perGenus >= 0.005) %>%
#     select(Site, Window, Order, Genus, pSamples_perSiteWindow_perGenus)
# }
# 
# win456_GenusDetect_plotdat <- windows2016_getTopGeneraData_function(win456_otutable_long_wtaxa_wmeta)
# win567_GenusDetect_plotdat <- windows2016_getTopGeneraData_function(win567_otutable_long_wtaxa_wmeta)
# 
# ## what orders remain in each of these plots?
# win456_GenusDetect_plotdat %>% distinct(Order) %>% arrange(Order) %>% pull()
# #"Araneae" "Coleoptera" "Diptera" "Ephemeroptera" "Hemiptera" "Hymenoptera" "Lepidoptera" "Megaloptera" "Neuroptera" "Sarcoptiformes" "Trichoptera" "Trombidiformes"
# win456_generaPal <- c('gray25', "orchid2", "tan4", "#3E6C54", "#FDAC9F", "darkorange", "#808133", "cadetblue3", '#32CD32', 'gray75', "turquoise4", 'gray50')
# 
# win567_GenusDetect_plotdat %>% distinct(Order) %>% arrange(Order) %>% pull()
# #"Araneae" "Blattodea" "Coleoptera" "Diptera" "Ephemeroptera" "Hemiptera" "Hymenoptera" "Lepidoptera" "Megaloptera" "Mesostigmata" "Psocodea" "Trichoptera" "Trombidiformes"
# win567_generaPal <- c('gray25', 'darkgoldenrod', "orchid2", "tan4", "#3E6C54", "#FDAC9F", "darkorange", "#808133", "cadetblue3", 'gray90', '#e1b580', "turquoise4", 'gray50')
# 
# 
# ## plotting win456 in a 2-row setup because 7 sites is just too bunched together to plot all the genus labels
# p_genera456text <- ggplot() +
#   geom_point(data=win456_GenusDetect_plotdat,
#              aes(x=Window, y=pSamples_perSiteWindow_perGenus),
#              show.legend = TRUE) +
#   geom_text_repel(data=win456_GenusDetect_plotdat, 
#                   aes(x=Window, y=pSamples_perSiteWindow_perGenus, label=Genus, 
#                       color=Order, fontface = "bold"),
#                   direction = "y", nudge_x = -0.05,
#                   hjust = 1,
#                   segment.size = 0.2,
#                   segment.colour = "black",
#                   size = 2,
#                   min.segment.length = 0) +
#   scale_color_manual(values = win456_generaPal) +
#   scale_y_continuous(breaks = c(0, 0.5, 1.0), limits = c(0,1)) +
#   scale_x_continuous(breaks = c(4, 5, 6), limits = c(3.3, 6.2)) +
#   facet_wrap(~Site, scales = "free_x", ncol=4) +
#   labs(x = "\n Window", y="fraction samples with taxa", color = "Order") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.position = "none")
# 
# ## for win567, plotting all as single row (4 wide) to match spacing with win456
# p_genera567text <- ggplot() +
#   geom_point(data=win567_GenusDetect_plotdat,
#              aes(x=Window, y=pSamples_perSiteWindow_perGenus),
#              show.legend = TRUE) +
#   geom_text_repel(data=win567_GenusDetect_plotdat, 
#                   aes(x=Window, y=pSamples_perSiteWindow_perGenus, label=Genus, 
#                       color=Order, fontface = "bold"),
#                   direction = "y", nudge_x = -0.05,
#                   hjust = 1,
#                   segment.size = 0.2,
#                   segment.colour = "black",
#                   size = 2,
#                   min.segment.length = 0) +
#   scale_color_manual(values = win567_generaPal) +
#   scale_y_continuous(breaks = c(0, 0.5, 1.0), limits = c(0,1)) +
#   scale_x_continuous(breaks = c(5, 6, 7), limits = c(4.3, 7.2)) +
#   facet_wrap(~Site, scales = "free_x", nrow=2, ncol=2) +
#   labs(x = "\n Window", y="fraction samples with taxa", color = "Order") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.position = "none")
# 
# ## plot to stitch together, to get the heights to be equivalent per facet
# windows2016_generaLayout <- "
# AAB
# "
# p_5c <- p_genera456text + p_genera567text + plot_layout(design = windows2016_generaLayout)
# p_5c
# 
# ggsave("~/github/nhguano/figures/figure4c_windows2016_topGenera.png", height = 30, width = 50, dpi = 150, units="cm")
# ggsave("~/github/nhguano/figures/figure4c_windows2016_topGenera.pdf", height = 30, width = 50, dpi = 300, units="cm")
# 
# 
# ## get a legend for these plots...
# windows2016_generaPal <- c('gray25',
#                            'darkgoldenrod',
#                            "orchid2", "tan4", "#3E6C54", "#FDAC9F", "darkorange", "#808133", "cadetblue3", 
#                            'gray90',
#                            '#32CD32', 
#                            'orchid2',
#                            'gray75', "turquoise4", 'gray50')
# 
# p_windows2016_generaLegend <- rbind(win456_GenusDetect_plotdat, win567_GenusDetect_plotdat) %>% 
#   ggplot() +
#   geom_point(aes(x=Site, y=Window, color=Order)) +
#   scale_color_manual(values = windows2016_generaPal)
# 
# p_windows2016_generaLegend
# ggsave("~/github/nhguano/figures/figure4legend_windows2016_topGenera_legendPalleteOnly.png", height = 15, width = 5, dpi = 150, units="cm")
# ggsave("~/github/nhguano/figures/figure4legend_windows2016_topGenera_legendPalleteOnly.pdf", height = 15, width = 5, dpi = 300, units="cm")
# 
# 
# ## part 2e:
# ## indicator species work: 
# ## focusing on the 10 orders listed in fig2d - the ones with at least 2% detections across any site+window
# top567Orders <- win567_OrderDetect_plotdat %>% distinct(Order) %>% filter(Order != "other") %>% pull()
# 
# win567_indiDataTable <- get_indispecDataTable_function(win567_meta, top567Orders)
# 
# ### perform indicator species analysis three different ways:
# ### (A) by Windows, (B) by Site, or (C) by Site+Window combination
# #### by Windows first:
# ## generate data, filter for pvalue sig threshold p<=0.05
# win567_indispec_list_windows <- doIndiSpec_function(win567_indiDataTable, win567_meta, "Window", 3)
# win567_indispec_plotdat_windows <- win567_indispec_list_windows$sign %>% 
#   mutate(taxaName = row.names(.)) %>% 
#   separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
#   filter(p.value <= 0.05) %>% 
#   select(index, stat, Order, Genus) %>% 
#   merge(., 
#         data.frame(indexName = colnames(win567_indispec_list_windows$comb)) %>% ## gets the index names to match group levels
#           mutate(index = row.names(.)))
# 
# ## set levels for plot, ordering group that spans windows between distinct windows; 
# ## put shared begin/end windows at end of plot
# win567_indispec_plotdat_windows$indexName <- factor(win567_indispec_plotdat_windows$indexName, levels = c(
#   "5", "5+6", "6", "6+7", "7", "5+7"))
# ## plot
# p_win567_indispec_windows <- ggplot(win567_indispec_plotdat_windows,
#                                     aes(x=indexName, y=Genus, fill=stat)) +
#   geom_tile(color="black") +
#   scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
#                        breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
#                        limits = c(0,1)) +
#   facet_grid(Order ~ ., scales = "free_y", switch = "y", space="free") +
#   scale_x_discrete(labels = c("5", "5\n+\n6", "6", "6\n+\n7", "7", "5\n+\n7")) +
#   labs(x="window", y="", fill = "association\nstatistic") +
#   theme_classic() +
#   theme(strip.text.y.left = element_text(angle = 0),
#         strip.placement = "outside",
#         aspect.ratio = 0.9,
#         strip.background = element_rect(color = NA, fill="gray90"),
#         panel.grid.major.y = element_line(color="gray85"),
#         panel.grid.major.x = element_line(color="gray85"))
# p_win567_indispec_windows
# ggsave("~/github/nhguano/figures/figure4c_specindi_byWindow_win567.png", width = 150, height = 250, units="mm", dpi=150)
# ggsave("~/github/nhguano/figures/figure4c_specindi_byWindow_win567.pdf", width = 150, height = 250, units="mm", dpi=300)
# 
# ## repeat for Sites (instead of window groups)
# win567_indispec_list_sites <- doIndiSpec_function(win567_indiDataTable, win567_meta, "Site", 3)
# win567_indispec_plotdat_sites <- win567_indispec_list_sites$sign %>% 
#   mutate(taxaName = row.names(.)) %>% 
#   separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
#   filter(p.value <= 0.05) %>% 
#   select(index, stat, Order, Genus) %>% 
#   merge(., 
#         data.frame(indexName = colnames(win567_indispec_list_sites$comb)) %>% ## gets the index names to match group levels
#           mutate(index = row.names(.)))
# 
# ## what are our levels to order this time?
# win567_indispec_plotdat_sites %>% distinct(indexName, index) %>% arrange(index) 
# 
# # ## set levels for x axis:
# win567_indispec_plotdat_sites$indexName <-
#   factor(win567_indispec_plotdat_sites$indexName,
#          levels = win567_indispec_plotdat_sites %>% distinct(indexName, index) %>% arrange(index) %>% pull(indexName))
# 
# 
# ## plot
# p_win567_indispec_sites <- ggplot(win567_indispec_plotdat_sites,
#                                   aes(x=indexName, y=Genus, fill=stat)) +
#   geom_tile(color="black") +
#   scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
#                        breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
#                        limits = c(0,1)) +
#   facet_grid(Order ~ ., scales = "free_y", switch = "y", space="free") +
#   labs(x="site", y="", fill = "association\nstatistic") +
#   theme_classic() +
#   theme(strip.text.y.left = element_text(angle = 0),
#         strip.placement = "outside",
#         aspect.ratio = 1,
#         strip.background = element_rect(color = NA, fill="gray90"),
#         axis.text.x = element_text(angle = 90, hjust=1),
#         panel.grid.major.y = element_line(color="gray85"),
#         panel.grid.major.x = element_line(color="gray85"))
# p_win567_indispec_sites
# ggsave("~/github/nhguano/figures/figure4c_specindi_bySite_win567.png", width = 150, height = 250, units="mm", dpi=150)
# ggsave("~/github/nhguano/figures/figure4c_specindi_bySite_win567.pdf", width = 150, height = 250, units="mm", dpi=300)
# 
# ## let's try it for the various Site+Window groups too!
# ### note that we're limiting our observations to ONLY singleton groups... that is...
# ## the site+window combination where an indicator species is ONLY associated with that one site+window
# ## combinations beyond this get pretty unweildy, but between 2-3 levels are still manageable if reviewers are concerned
# win567_indispec_list_sitewindows <- doIndiSpec_function(win567_indiDataTable, 
#                                                         win567_meta %>% mutate(SiteWindow = paste0(Site,Window)), 
#                                                         "SiteWindow", 1)
# win567_indispec_plotdat_sitewindows <- win567_indispec_list_sitewindows$sign %>% 
#   mutate(taxaName = row.names(.)) %>% 
#   separate(col=taxaName, into=c('Order', 'Genus'), sep="-", extra = "merge") %>% 
#   filter(p.value <= 0.05) %>% 
#   select(index, stat, Order, Genus) %>% 
#   merge(., 
#         data.frame(indexName = colnames(win567_indispec_list_sitewindows$comb)) %>% ## gets the index names to match group levels
#           mutate(index = row.names(.))) %>% 
#   mutate(Site = substr(indexName, 1, 3),
#          Window = substr(indexName, 4, 4))
# 
# # ## set levels for x axis:
# win567_indispec_plotdat_sitewindows$indexName <-
#   factor(win567_indispec_plotdat_sitewindows$indexName,
#          levels = win567_indispec_plotdat_sitewindows %>% distinct(indexName, index) %>% arrange(index) %>% pull(indexName))
# 
# 
# ## plot
# p_win567_indispec_sitewindows <- ggplot(win567_indispec_plotdat_sitewindows,
#                                         aes(x=Window, y=Genus, fill=stat)) +
#   geom_tile(color="black") +
#   scale_fill_gradientn(colours = scico::scico(30, palette = "lajolla"), na.value = "transparent",
#                        breaks = c(0.2, 0.5, 0.8), labels = c(0.2, 0.5, 0.8),
#                        limits = c(0,1)) +
#   facet_grid(Order ~ Site, scales = "free", switch = "y", space="free") +
#   labs(x="window", y="", fill = "association\nstatistic") +
#   theme_classic() +
#   theme(strip.text.y.left = element_text(angle = 0),
#         strip.placement = "outside",
#         aspect.ratio = 1,
#         strip.background = element_rect(color = NA, fill="gray90"),
#         panel.grid.major.y = element_line(color="gray85"),
#         panel.grid.major.x = element_line(color="gray85"))
# p_win567_indispec_sitewindows
# ggsave("~/github/nhguano/figures/figure4c_specindi_bySiteWindow_win567.png", width = 150, height = 250, units="mm", dpi=150)
# ggsave("~/github/nhguano/figures/figure4c_specindi_bySiteWindow_win567.pdf", width = 150, height = 250, units="mm", dpi=300)


