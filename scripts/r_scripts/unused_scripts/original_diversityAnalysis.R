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
# ## "#D49347" == Psocodoea
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
# ############333 !!!!!!!!!!!!!!!!!
# ##### CAUTION: note that we're grouping shared Genus ONLY AMONG THOSE CORE FEATURES (OTUs), not across all OTUs assigned to that Genus
# 
ggsave("~/github/nhguano/figures/figure2_corefeatures_byGenus.png", height = 12, width = 17, units="cm")
ggsave("~/github/nhguano/figures/figure2_corefeatures_byGenus.pdf", height = 12, width = 17, units="cm")
#ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/corefeatures_byGenus.svg", height = 12, width = 17, units="cm")
# 
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
# ## bundle OTU table, metadata, and tree info into single phyloseq object
fox2016_phydat <- phyloseq(fox2016_phyOTU, fox2016_meta)
fox2016_phydat <- merge_phyloseq(fox2016_phydat, tree)
rm(fox2016_phyOTU)
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

write_csv(fox2016_dunn_all,
          path = "~/github/nhguano/data/text_tables/fox2016_dunn_allMetrics.csv")

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
ggsave("~/github/nhguano/supplementaryData/figureS1_fox2016_alpha_boxplots.png", height = 15, width = 10, units = 'cm', dpi=150)
ggsave("~/github/nhguano/supplementaryData/figureS1_fox2016_alpha_boxplots.pdf", height = 15, width = 10, units = 'cm', dpi=300)
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
ggsave("~/github/nhguano/supplementaryData/figureS2_fox2016_alpha_heatmap_pvals.png", height = 15, width = 23, units = 'cm', dpi=150)
ggsave("~/github/nhguano/supplementaryData/figureS2_fox2016_alpha_heatmap_pvals.pdf", height = 15, width = 10, units = 'cm', dpi=300)

# ## alpha diversity cleanup:
rm(fox2016_alpha_df, fox2016_alpha_df_long, fox2016_labels_H, fox2016_labels_PD, fox2016_labels_SR,
   kw_input_fox2016_H, kw_input_fox2016_PD, kw_input_fox2016_SR, fox2016_LetterLabels_all, fox2016_dunn_all, 
   fox2016_alpha_plotdat, p_heatmap_fox2016_adj, p_heatmap_fox2016_unadj, dunns_letters_fox2016, dunns_pvals_fox2016)
# 
# 
# ##################
# ## part 2a: community composition
# ## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method
# fox2016_dist_ds <- phyloseq::distance(fox2016_phydat, "bray", binary = TRUE)
# fox2016_dist_bc <- phyloseq::distance(fox2016_phydat, "bray", binary = FALSE)
# fox2016_dist_uu <- phyloseq::distance(fox2016_phydat, "unifrac", weighted=FALSE)
# fox2016_dist_wu <- phyloseq::distance(fox2016_phydat, "wunifrac")
# 
# ##################
# ## part 2b: community composition
# ## 2bi. run PERMANOVA (via Adonis) for each distance method to evaluate if samples within each Window (group) have different centroids
# fox2016_adonis_function <- function(distanceData, distanceMetric){
#   adonis_tmp <- adonis2(distanceData ~ Window, data = fox2016_meta, )
#   data.frame(adonis_tmp) %>% mutate(Metric = distanceMetric,
#                                             Class = row.names(.))
# }
# 
# fox2016_adonis_ds <- fox2016_adonis_function(fox2016_dist_ds, "Dice-Sorensen")
# fox2016_adonis_bc <- fox2016_adonis_function(fox2016_dist_bc, "Bray-Curtis")
# fox2016_adonis_uu <- fox2016_adonis_function(fox2016_dist_uu, "Unifrac Unweighted")
# fox2016_adonis_wu <- fox2016_adonis_function(fox2016_dist_wu, "Unifrac Weighted")
# 
# fox2016_adonis_all <- rbind(fox2016_adonis_ds, fox2016_adonis_bc, fox2016_adonis_uu, fox2016_adonis_wu)
# rm(fox2016_adonis_ds, fox2016_adonis_bc, fox2016_adonis_uu, fox2016_adonis_wu)
#   ## Window is significant main effect regardless of distance metric
# 
# fox2016_adonis_all <- fox2016_adonis_all[,c(7,1,2,3,4,5,6)] %>% 
#   mutate(R2 = round(R2, 3),
#          SumOfSqs = round(SumOfSqs, 3),
#          F = round(F, 3)) %>% 
#   rename(SS = "SumOfSqs")
# # write.csv(fox2016_adonis_all, quote=FALSE, row.names = FALSE,
# #           file = "~/Documents/nau_projects/guano/NHguano_redux/new_data/fox2016_adonis_allMetrics.csv")
# 
# ## part 2bii. Test for homogeneity of dispersion
# fox2016_betadisper_function <- function(distanceData, distanceMetric){
#   tmp_disper_list <- betadisper(d = distanceData, group =  fox2016_meta$Window, type = c("median"))
#   tmp_disper_anova <- data.frame(anova(tmp_disper_list))
#   tmp_disper_anova <- tmp_disper_anova %>% 
#     mutate(Class = row.names(.)) %>% mutate(Metric = distanceMetric) %>% mutate(Factor = "Window")
#   tmp_disper_anova[,c(6,1,2,3,4,5,8,7)]
# }
# 
# fox2016_bdisp_ds <- fox2016_betadisper_function(fox2016_dist_ds, "Dice-Sorensen")
#   ## non significant; p = 0.163 (i.e. we can't reject null that groups have same dispersions... good!)
# fox2016_bdisp_bc <- fox2016_betadisper_function(fox2016_dist_bc, "Bray-Curtis")
#   ## non significant; p = 0.089
# fox2016_bdisp_uu <- fox2016_betadisper_function(fox2016_dist_uu, "Unifrac Unweighted")
#   ## non significant; p = 0.079
# fox2016_bdisp_wu <- fox2016_betadisper_function(fox2016_dist_wu, "Unifrac Weighted")
#   ## SIGNIFICANT!; p = 0.0082 ... need to treat this interpretation more carefully
# 
# fox2016_bdisp_all <- rbind(fox2016_bdisp_ds, fox2016_bdisp_bc, fox2016_bdisp_uu, fox2016_bdisp_wu)
# rm(fox2016_bdisp_ds, fox2016_bdisp_bc, fox2016_bdisp_uu, fox2016_bdisp_wu)
# fox2016_bdisp_all <- fox2016_bdisp_all %>% 
#   mutate(Sum.Sq = round(Sum.Sq, 3),
#          Mean.Sq = round(Mean.Sq, 3),
#          F.value = round(F.value, 3),
#          Pr..F. = round(Pr..F., 3)) %>% 
#   select(-Factor)
# # write_csv(fox2016_bdisp_all, 
# #           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/fox2016_betadisper_allMetrics.csv")
# 
# ##################
# ## part 2c: ordination of distance data sets
# ### ordinations:
# fox2016_ordi_function <- function(distanceValues, distanceMetric){
#   tmp_pcoa <- ordinate(fox2016_phydat, method="PCoA", distance = distanceValues)
#   tmp_pcoa_list <- plot_ordination(fox2016_phydat, tmp_pcoa)
#   tmp_pcoa_df <- tmp_pcoa_list$data %>% 
#     mutate(SampleID = row.names(.),
#            Axis1lab = tmp_pcoa_list$labels$x,
#            Axis2lab = tmp_pcoa_list$labels$y,
#            Metric = distanceMetric)
#   merge(tmp_pcoa_df, fox2016_meta) %>% mutate(Window = as.character(Window))
# }
# 
# fox2016_orddata_ds <- fox2016_ordi_function(fox2016_dist_ds, "Dice-Sorensen")
# fox2016_orddata_bc <- fox2016_ordi_function(fox2016_dist_bc, "Bray-Curtis")
# fox2016_orddata_uu <- fox2016_ordi_function(fox2016_dist_uu, "Unweighted-Unifrac")
# fox2016_orddata_wu <- fox2016_ordi_function(fox2016_dist_wu, "Weighted-Unifrac")
# 
# fox2016_orddata_all <- rbind(fox2016_orddata_ds, fox2016_orddata_bc, fox2016_orddata_uu, fox2016_orddata_wu)
# rm(fox2016_orddata_ds, fox2016_orddata_bc, fox2016_orddata_uu, fox2016_orddata_wu)
# 
# ## plot with metadata
# #### generate palette that works with color blind data
# scico::scico(length(unique(fox2016_pcoa_df$Window)), palette = "batlow")
#   # original: "#001959" "#184E60" "#577646" "#B28C32" "#FCA68C" "#F9CCF9"
# # going to modify blue and green values to make easier to contrast Windows 6, 7, 8
# fox2016_pal <- rev(c("#001447", "#267b97", "#557d3f", "#B28C32", "#FCA68C", "#F9CCF9"))
# 
# ## plot individually, then stitch together
# fox2016_ordplot_function <- function(inputdata){
#   ggplot(data = inputdata,
#          aes(x=Axis.1, y=Axis.2)) +
#     geom_point(aes(color = Window), size = 3) +
#     stat_ellipse(aes(group = Window, color = Window), alpha = 0.5) +
#     labs(x = unique(inputdata$Axis1lab),
#          y = unique(inputdata$Axis2lab)) +
#     theme_classic() +
#     scale_color_manual(values = fox2016_pal) +
#     coord_fixed()
# }
# 
# p_fox2016_ord_ds <- fox2016_ordplot_function(fox2016_orddata_all %>% filter(Metric == "Dice-Sorensen"))
# p_fox2016_ord_bc <- fox2016_ordplot_function(fox2016_orddata_all %>% filter(Metric == "Bray-Curtis"))
# p_fox2016_ord_uu <- fox2016_ordplot_function(fox2016_orddata_all %>% filter(Metric == "Unweighted-Unifrac"))
# p_fox2016_ord_wu <- fox2016_ordplot_function(fox2016_orddata_all %>% filter(Metric == "Weighted-Unifrac"))
# 
# legend_fox2016 <- get_legend(p_fox2016_ord_ds + guides(color = guide_legend(nrow = 1)))
# 
# tmp_fox2016_plot <- plot_grid(
#           p_fox2016_ord_ds + theme(legend.position = "none"), 
#           p_fox2016_ord_bc + theme(legend.position = "none"),
#           p_fox2016_ord_uu + theme(legend.position = "none"), 
#           p_fox2016_ord_wu + theme(legend.position = "none"),
#           labels = c("A", "B", "C", "D"),
#           nrow = 2, ncol = 2,
#           align = "hv")
# 
# plot_grid(legend_fox2016, tmp_fox2016_plot,
#           nrow = 2, ncol = 1, 
#           rel_heights = c(.2, 1))
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_ordinations_all.png", width = 20, height = 21, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_ordinations_all.svg", width = 20, height = 21, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_ordinations_all.pdf", width = 20, height = 21, units = 'cm')
# 
# rm(p_fox2016_ord_ds, p_fox2016_ord_bc, p_fox2016_ord_uu, p_fox2016_ord_wu)
# 
# ##################
# ## part 2d: pairwise Adonis to determine what, if any, community comopsition groups are different
# ## generate long form of all values for supplementary table
# fox2016_pairwise_adonis_func <- function(distancevals, metric){
#   tmp_pairwiseadonis <- pairwise.adonis(distancevals, fox2016_meta$Window, p.adjust.m = "bonferroni") %>% 
#     mutate(Metric = metric) %>% 
#     separate(col = pairs, into = c("Window_A", "Window_B"), sep = "vs") %>% 
#     arrange(p.adjusted)
# }
# 
# fox2016_pwadonis_ds <- fox2016_pairwise_adonis_func(fox2016_dist_ds, "Dice-Sorensen")
# fox2016_pwadonis_bc <- fox2016_pairwise_adonis_func(fox2016_dist_bc, "Bray-Curtis")
# fox2016_pwadonis_uu <- fox2016_pairwise_adonis_func(fox2016_dist_uu, "Unweighted-Unifrac")
# fox2016_pwadonis_wu <- fox2016_pairwise_adonis_func(fox2016_dist_wu, "Weighted-Unifrac")
# fox2016_pwadonis_all <- rbind(fox2016_pwadonis_ds, fox2016_pwadonis_bc, fox2016_pwadonis_uu, fox2016_pwadonis_wu)
# rm(fox2016_pwadonis_ds, fox2016_pwadonis_bc, fox2016_pwadonis_uu, fox2016_pwadonis_wu)
# # write_csv(fox2016_pwadonis_all, 
# #           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/fox2016_pairwiseAdonis_all.csv")
# 
# fox2016_pwadonis_all$Metric <- factor(fox2016_pwadonis_all$Metric, levels = c(
#   'Dice-Sorensen', 'Bray-Curtis', 'Unweighted-Unifrac', 'Weighted-Unifrac'))
# ggplot(fox2016_pwadonis_all %>% mutate(p.adjusted = round(p.adjusted, 2)), aes(x=Window_B,
#                                  y=Window_A,
#                                  fill=p.adjusted,
#                                  label = p.adjusted)) +
#   geom_tile(color = "black") +
#   geom_text() +
#   coord_fixed() +
#   scale_fill_viridis_c(option = "viridis", direction = -1, alpha = 0.5) +
#   facet_wrap(~Metric, nrow=2) +
#   labs(x = "Window", y = "Window", fill = "Bonferonni-adjusted\np-value") +
#   theme_classic() +
#   theme(legend.position = "top")
#   
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_pairwiseAdonis_heatmap_pvalues_all.png", width = 16, height = 16, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_pairwiseAdonis_heatmap_pvalues_all.svg", width = 16, height = 16, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_pairwiseAdonis_heatmap_pvalues_all.pdf", width = 16, height = 16, units = 'cm')
# 
# ## cleanup
# # rm(list = ls(pattern = "^fox2016"))
# # rm(fox_OTUtable_wide_raw, legend_fox2016, tmp_fox2016_plot)
# 
# ##################
# ## part 3: we find that early/late windows are different in terms of community composition... so what are those particular taxa?
# 
# ## create fox2016_meta again:
# fox2016_meta <- otu_table_long_filtd_wMeta %>% 
#   mutate(Window = as.factor(Window)) %>% 
#   filter(Year == 2016 & Site == "FOX" & Window != 9) %>% 
#   select(SampleID, Window, newDate) %>% 
#   distinct()
# 
# ## pivot wide to long for SRS output, merge with taxonomy and metadata info
# fox_OTUtable_wide_norm$OTUid <- row.names(fox_OTUtable_wide_norm)
# fox_OTUtable_long_norm <- fox_OTUtable_wide_norm %>% 
#   pivot_longer(-OTUid, names_to = "SampleID", values_to = "SRSreads") %>% 
#   filter(SRSreads > 0)
# fox_taxa <- otu_table_long_filtd_wMeta %>% 
#   filter(SampleID %in% fox_OTUtable_long_norm$SampleID) %>% 
#   distinct(OTUid, OTUalias, Class, Order, Family, Genus, Species)
# fox_OTUtable_long_norm_wTaxa <- merge(fox_OTUtable_long_norm, fox_taxa)
# fox_OTUtable_long_norm_wTaxa_wMeta <- merge(fox_OTUtable_long_norm_wTaxa, fox2016_meta)
# 
# rm(fox_OTUtable_wide_norm, fox_OTUtable_long_norm, fox_OTUtable_long_norm_wTaxa)
# 
# ## calculate fraction of reads at Order and Genus-levels (separately)...
# ## could also do this strictly for detections...
# ### first, gather the total number of reads per sampling window (accounting for different sample sizes in each window)
# fox2016_sumry_PerWindow <- 
#   fox_OTUtable_long_norm_wTaxa_wMeta %>%
#   group_by(Window) %>%
#   summarise(nReads_perWindow_global = sum(SRSreads),
#             nSamples_perWindow = n_distinct(SampleID))
# 
# ### next, summarise the total number of reads per sampling window for each arthropod Order
# fox2016_sumry_ReadsPerWindow_perOrder <- 
#   fox_OTUtable_long_norm_wTaxa_wMeta %>% 
#   group_by(Window, Class, Order) %>% 
#   summarise(nReads_perWindow_perOrder = sum(SRSreads)) %>% 
#   merge(., fox2016_sumry_PerWindow) %>% 
#   mutate(pReads_perWindow_perOrder = nReads_perWindow_perOrder / nReads_perWindow_global) %>% 
#   mutate(pReads_perWindow_perOrder = round(pReads_perWindow_perOrder, 3))
# 
# ## keep palette consistent with other figure of broad diet pattern:
# ## 'gray25'  == Araneae; 'darkgoldenrod' == Blattodea; "orchid2" == Coleoptera; "tan4"  == Diptera
# ## "#3E6C54" == Ephemeroptera; "#FDAC9F" == Hemiptera; "darkorange" == Hymenoptera; "#808133" == Lepidoptera
# ## "cadetblue3" == Megaloptera; "#D49347" == Psocodoea; "turquoise4" == Trichoptera; 'gray50' == Trombidiformes
# 
# # fox2016_sumyr_ReadsPerWindow_perOrder %>% filter(pReads_perWindow_perOrder >= 0.02) %>% distinct(Order) %>% arrange(Order) %>% pull()
# # > "Coleoptera"     "Diptera"        "Ephemeroptera"  "Hemiptera"      "Hymenoptera"    "Lepidoptera"    "Trichoptera"   "Trombidiformes"
# 
# fox2016_orderPal <- c("orchid2", "tan4", "#3E6C54", "#FDAC9F", 'darkorange', "#808133", "turquoise4", 'gray50')
# 
# ## and plot
# ggplot(fox2016_sumry_ReadsPerWindow_perOrder %>% filter(pReads_perWindow_perOrder >= 0.02),
#        aes(x = Window, 
#            y = pReads_perWindow_perOrder, 
#            fill = Order)) +
#   geom_col(color="gray25") +
#   theme_classic() +
#   labs(x = '\nWindow', y = 'fraction of sequences\n', fill = 'Arthropod\norder') +
#   scale_fill_manual(values = fox2016_orderPal)
# 
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_barplot_arthOrders_byWindows.png", width = 10, height = 8, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_barplot_arthOrders_byWindows.svg", width = 10, height = 8, units = 'cm')
# # ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/fox2016_barplot_arthOrders_byWindows.pdf", width = 10, height = 8, units = 'cm')
# 
# ### a bit more refinfed: summarise the abundant taxa in each sampling window by arthropod Genus
#   ### can look at it in terms of fraction of samples with taxa detected per window, or...
#   ###... by the fraction of reads per taxa per window
# 
# ## going to summarise these data by finding the top10 values for each window, then creating a data table combining all those values
#   ### note we're including distinct OTUs with "NA" in this output, reformatting them to their particular OTU and family name
#   ### (this only occurs in two of the 38 elements)
# 
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
#   select(Window, Order, Genus, pReads_perWindow_perGenus) %>% 
#   mutate(pReads_perWindow_perGenus = pReads_perWindow_perGenus * 100) %>% 
#   arrange(Order) %>% 
#   pivot_wider(values_from = "pReads_perWindow_perGenus", names_from = "Window", values_fill = 0)
# 
# ## save to disk
# # write_csv(fox2016_topGenera_byWindows_grid,
# #           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/fox2016_topGenera_byWindows.csv")
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
# ## cleanup:
# rm(list = ls(pattern = "^fox*"))

########################################
## part 4 - diversity estimates: multipe sites, same year
## how does diversity change between sites and sampling windows in a single year?
## using same 37-day sampling window

## 0a. justify why we're using sites selected
## 0b. normalize reads using SRS method for just these samples in these sites, then create phyloseq object for diversity calculations

## 1a. calculate alpha diversity using 3 metrics (Richness, Shannon's, Faith's PD)
## 1b. compare sampling windows for each metric using Kruskal-Wallis and pairwise Wilcoxon
## also going to generate a heatmap of these pvalues given there are so many pairwise comps with Wilcoxon test
## also going to provide Lettered labels for pairwise comp differences between groups

## 2a. calculate community composition distances using 4 metrics (Dice-Sorensen, Bray-Curtis, Unweighted/Weighted Unifrac)
## 2b. ordinate each distance for:
#  2bii. between-group medians (centroids)
# 2biii. within-group dispersions (via betadisper)
## 2c. ordinate each distance metric with PCoA
## 2d. pairwise adonis and heatmap for comparing differences between particular clusters

########################################

##################
## part 0a: justification for why we're analyzing particular site-windows:
otu_table_long_filtd_wMeta %>% 
  filter(Year == 2016) %>% 
  group_by(Site, Window) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>%
  arrange(Window) %>% 
  pivot_wider(names_from = "Window", values_from = "nSamples")
## selecting 2 different window spans; all site-windows have minimum 6 samples per site-window: 
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


##################
## part 0b: rarefy table using SRS for samples in each group (win456 and win567), ...
## then generate phyloseq object by importing with metadata and tree file to generate phyloseq object

## get metadata for these samples and then generate OTU table for alpha/adonis/betadisp/ordinations:
win456_meta <- otu_table_long_filtd_wMeta %>% 
  filter(SampleID %in% win456_samplenames) %>% 
  select(SampleID, Window, Site, newDate) %>% 
  mutate(Window = as.factor(Window)) %>% 
  distinct()
## 331 samples total

win567_meta <- otu_table_long_filtd_wMeta %>% 
  filter(SampleID %in% win567_samplenames) %>% 
  select(SampleID, Window, Site, newDate) %>% 
  mutate(Window = as.factor(Window)) %>% 
  distinct()
## 166 samples total

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
  tmp_OTUtable_wide_norm <- SRS(tmp_OTUtable_wide_raw, tmp_Cmin)
  row.names(tmp_OTUtable_wide_norm) <- row.names(tmp_OTUtable_wide_raw)
  
  tmp_phyOTU <- otu_table(tmp_OTUtable_wide_norm, taxa_are_rows = TRUE)
  tmp_phydat <- phyloseq(tmp_phyOTU, metadata_file)
  merge_phyloseq(tmp_phydat, tree)
}

win456_phydat <- windows2016_getPhyloseq_function(win456_meta)
win567_phydat <- windows2016_getPhyloseq_function(win567_meta)

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
win567_alpha_df <- windows2016_getAlphaVals_function(win567_phydat, win567_meta)

##################
## part 1b: applying KW test for global differences in distinct diversity metrics for a given Site+Window group
## ...then Wilcoxon rank sum for pairwise diffs in alpha vals by Site+Group

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
##significant; chi-squared = 49.586, df = 20, p-value = 0.00025
win456_KW_result_H <- windows2016_getAlphaStats_function(win456_alpha_df, "H")
##significant; chi-squared = 39.058, df = 20, p-value = 0.00656
win456_KW_result_PD <- windows2016_getAlphaStats_function(win456_alpha_df, "PD")
##significant; chi-squared = 91.158, df = 20, p-value = 0

win567_KW_result_SR <- windows2016_getAlphaStats_function(win567_alpha_df, "SR")
##significant; chi-squared = 46.508, df = 11, p-value = 0
win567_KW_result_H <- windows2016_getAlphaStats_function(win567_alpha_df, "H")
## NOT! significant; chi-squared = 19.035, df = 11, p-value = 0.0605
win567_KW_result_PD <- windows2016_getAlphaStats_function(win567_alpha_df, "PD")
##significant; chi-squared = 45.035, df = 11, p-value = 0

## combine into supplementary table and report?
windows2016_kW_result_SR <- rbind(win456_KW_result_SR, win456_KW_result_H, win456_KW_result_PD, 
                                  win567_KW_result_SR, win567_KW_result_H, win567_KW_result_PD)
rm(win456_KW_result_SR, win456_KW_result_H, win456_KW_result_PD, 
   win567_KW_result_SR, win567_KW_result_H, win567_KW_result_PD)

# write_csv(windows2016_kW_result_SR,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/windows2016_KWvalues_all.csv")

rm(windows2016_kW_result_SR)

## get wilcoxon pairwise vals for each group (Site+Window):
windows2016_wilcoxon_function <- function(alphadat, metric, windowgroup){
  tmp_wilcoxon_input <- alphadat %>% 
    mutate(Grouper = (paste0(Window, Site))) %>% 
    filter(Metric == metric) %>% 
    select(Grouper, Alpha_value)
  attach(tmp_wilcoxon_input)
  tmp_tri <- pairwise.wilcox.test(Alpha_value, Grouper, p.adjust.method = "BH")
  detach()
  h0_mat <- data.frame(tmp_tri$p.value)
  h0_mat$Pair <- row.names(h0_mat)
  h0_mat %>% 
    pivot_longer(-Pair, names_to="Names", values_to = "BHpval") %>% 
    filter(!is.na(`BHpval`)) %>%
    mutate(Names = sub("X", "", Names),
           Metric = metric,
           Windows = windowgroup,
           `BHpval` = round(`BHpval`, 3)) %>% 
    arrange(Pair, Names)
}

win456_wilcox_SR <- windows2016_wilcoxon_function(win456_alpha_df, "SR", "456")
## only 5 pairwise comps with pval <= 0.05! (note 5EPS site as culprit?)
win456_wilcox_H <- windows2016_wilcoxon_function(win456_alpha_df, "H", "456")
## only 4 pairwise comps with pval <= 0.05! (same 5EPS group...)
win456_wilcox_PD <- windows2016_wilcoxon_function(win456_alpha_df, "PD", "456")
## 41! pairwise comps with pval <= 0.05... why the big diff?

win567_wilcox_SR <- windows2016_wilcoxon_function(win567_alpha_df, "SR", "567")
## chichester seems very different here
win567_wilcox_H <- windows2016_wilcoxon_function(win567_alpha_df, "H", "567")
## 5CHI again. only 3 with pval <= 0.05
win567_wilcox_PD <- windows2016_wilcoxon_function(win567_alpha_df, "PD", "567")
## similar trend of many more sig differences in groups with Faith's PD calc...

windows2016_wilcox_all <- rbind(win456_wilcox_SR, win456_wilcox_H, win456_wilcox_PD,
                                win567_wilcox_SR, win567_wilcox_H, win567_wilcox_PD)
#rm(win456_wilcox_SR, win456_wilcox_H, win456_wilcox_PD, win567_wilcox_SR, win567_wilcox_H, win567_wilcox_PD)

windows2016_wilcox_all <- windows2016_wilcox_all %>% 
  rename("Window_A" = Pair, "Window_B" = Names, "pval BHcorrected" = BHpval)
# write_csv(windows2016_wilcox_all,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/windows2016_wilcox_all.csv")

## easier to plot these pairwise comps than to list all of them?
windows2016_wilcox_all$Metric <- factor(windows2016_wilcox_all$Metric, levels = c(
  "SR", "H", "PD"))

p_wilcox_456 <- ggplot(windows2016_wilcox_all %>% 
                         filter(Windows == "456"), 
                       aes(x=Window_B, y=Window_A, fill=`pval BHcorrected`)) +
  geom_tile(color="black") +
  facet_grid(Windows ~ Metric, space = "free") +
  scale_fill_viridis_c(option = "viridis", direction = -1, alpha = 0.7) +
  coord_fixed() +
  labs(x="", y="", fill = "BH-corrected\np-value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        strip.text.y = element_blank())

p_wilcox_567 <- ggplot(windows2016_wilcox_all %>% 
                         filter(Windows == "567") %>% 
                         mutate(Windows = "5 | 6 | 7"), 
                       aes(x=Window_B, y=Window_A, fill=`pval BHcorrected`)) +
  geom_tile(color="black") +
  facet_grid(Windows ~ Metric, space = "free") +
  scale_fill_viridis_c(option = "viridis", direction = -1, alpha = 0.7) +
  coord_fixed() +
  labs(x="", y="", fill = "BH-corrected\np-value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        strip.text.y = element_blank())

windows2016_wilcox_legend <- get_legend(p_wilcox_456 + theme(legend.position = "top"))

plot_grid(windows2016_wilcox_legend,
          p_wilcox_456 + theme(legend.position = "none"),
          p_wilcox_567 + theme(legend.position = "none"),
          nrow = 3, ncol = 1,
          labels = c(NA, "A", "B"),
          rel_heights = c(0.2, 1, 1))

## save plot
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/windows2016_wilcoxon_heatmap_pvalues_all.png", width = 19, height = 19, units = 'cm')
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/windows2016_wilcoxon_heatmap_pvalues_all.svg", width = 19, height = 19, units = 'cm')
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/windows2016_wilcoxon_heatmap_pvalues_all.pdf", width = 19, height = 19, units = 'cm')

rm(p_wilcox_456, p_wilcox_567, windows2016_wilcox_legend)

## plot these groups as boxplot too, and generate alphabetical labels to show significant pairwise differences...
win456_alpha_df$Metric <- factor(win456_alpha_df$Metric, levels=c("SR", "H", "PD"))
win567_alpha_df$Metric <- factor(win567_alpha_df$Metric, levels=c("SR", "H", "PD"))

p_boxplot_win456 <- ggplot(win456_alpha_df %>% mutate(Grouper = paste0(Window, Site)), 
                           aes(x=Grouper, y=Alpha_value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(alpha=0.5, width = 0.1) +
  facet_grid(Metric ~ Site, scales = "free", space="free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x="", y="Alpha value")

p_boxplot_win567 <- ggplot(win567_alpha_df %>% mutate(Grouper = paste0(Window, Site)), 
                           aes(x=Grouper, y=Alpha_value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(alpha=0.5, width = 0.1) +
  facet_grid(Metric ~ Site, scales = "free", space="free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x="", y="Alpha value")

plot_grid(p_boxplot_win456,p_boxplot_win567,
          labels = c("A", "B"),
          nrow=2, ncol=1)

# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/windows2016_boxplot_alphavals_all.png", width = 25, height = 25, units = 'cm')
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/windows2016_boxplot_alphavals_all.pdf", width = 25, height = 25, units = 'cm')
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/windows2016_boxplot_alphavals_all.svg", width = 25, height = 25, units = 'cm')

rm(p_boxplot_win456, p_boxplot_win567)

## generate pariwise Wilcoxon matrix of differences according to lettered labels (instead of raw pvalues)
windows2016_wilcoxonLetters_function <- function(alphadat, metric, windowgroup){
  tmp_wilcoxon_input <- alphadat %>% 
    mutate(Grouper = (paste0(Window, Site))) %>% 
    filter(Metric == metric) %>% 
    select(Grouper, Alpha_value)
  attach(tmp_wilcoxon_input)
  tmp_tri <- pairwise.wilcox.test(Alpha_value, Grouper, p.adjust.method = "BH")
  detach()
  h0_mat <- tmp_tri$p.value
  newcol_vals <- as.matrix(rep(NA, length(row.names(h0_mat))))
  h1_mat <- cbind(h0_mat, newcol_vals)
  colnames(h1_mat)[ncol(h1_mat)] <- row.names(h0_mat) %>% tail(1)
  newrow_vals <- t(as.matrix(rep(NA, length(colnames(h1_mat)))))
  h2_mat <- rbind(newrow_vals, h1_mat)
  row.names(h2_mat)[1] <- colnames(h2_mat)[1]
  h2_mat <- Matrix::forceSymmetric(h2_mat,uplo="L")
  h2_mat <- as.matrix(h2_mat)
  h2_mat[is.na(h2_mat)] <- 1
  tmp_lmat <- multcompLetters(h2_mat, compare="<=", threshold=0.05, Letters=letters)
  data.frame(tmp_lmat$Letters) %>% 
    rename(Letters = tmp_lmat.Letters) %>% 
    mutate(Window = row.names(.),
           Metric = metric,
           Windows = windowgroup) %>%
    arrange(Letters)
}

win456_wilcoxLetters_SR <- windows2016_wilcoxonLetters_function(win456_alpha_df, "SR", "456")
win456_wilcoxLetters_H <- windows2016_wilcoxonLetters_function(win456_alpha_df, "H", "456")
win456_wilcoxLetters_PD <- windows2016_wilcoxonLetters_function(win456_alpha_df, "PD", "456")

win567_wilcoxLetters_SR <- windows2016_wilcoxonLetters_function(win567_alpha_df, "SR", "567")
win567_wilcoxLetters_H <- windows2016_wilcoxonLetters_function(win567_alpha_df, "H", "567")
win567_wilcoxLetters_PD <- windows2016_wilcoxonLetters_function(win567_alpha_df, "PD", "567")

windows2016_wilcoxLetters_all <- rbind(win456_wilcoxLetters_SR, win456_wilcoxLetters_H, win456_wilcoxLetters_PD,
                                       win567_wilcoxLetters_SR, win567_wilcoxLetters_H, win567_wilcoxLetters_PD)

# write_csv(windows2016_wilcoxLetters_all, 
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/windows2016_wilcoxLetters_all.csv")

rm(win456_wilcoxLetters_SR, win456_wilcoxLetters_H, win456_wilcoxLetters_PD,
   win567_wilcoxLetters_SR, win567_wilcoxLetters_H, win567_wilcoxLetters_PD,
   windows2016_wilcoxLetters_all)

rm(win456_alpha_df, win567_alpha_df, windows2016_wilcox_all, windows2016_wilcox_legend,
   windows2016_getAlphaStats_function, windows2016_getAlphaVals_function, 
   windows2016_wilcoxon_function, windows2016_wilcoxonLetters_function)

##################
## part 2a:

## 2a. calculate community composition distances using 4 metrics (Dice-Sorensen, Bray-Curtis, Unweighted/Weighted Unifrac)
## 2b. ordinate each distance for:
#  2bii. between-group medians (centroids)
# 2biii. within-group dispersions (via betadisper)
## 2c. ordinate each distance metric with PCoA
## 2d. pairwise adonis and heatmap for comparing differences between particular clusters

## part 2a: community composition
## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method, per window (456 or 567)
win456_dist_ds <- phyloseq::distance(win456_phydat, "bray", binary = TRUE)
win456_dist_bc <- phyloseq::distance(win456_phydat, "bray", binary = FALSE)
win456_dist_uu <- phyloseq::distance(win456_phydat, "unifrac", weighted=FALSE)
win456_dist_wu <- phyloseq::distance(win456_phydat, "wunifrac")

win567_dist_ds <- phyloseq::distance(win567_phydat, "bray", binary = TRUE)
win567_dist_bc <- phyloseq::distance(win567_phydat, "bray", binary = FALSE)
win567_dist_uu <- phyloseq::distance(win567_phydat, "unifrac", weighted=FALSE)
win567_dist_wu <- phyloseq::distance(win567_phydat, "wunifrac")


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
    select(Class, Df, SumsOfSqs, R2, Fmodel, `Pr..F.`, Metric, WindowGroup)
}

win456_adonis_ds <- windows2016_adonis_function(win456_dist_ds, "Dice-Sorensen", win456_meta, "456")
win456_adonis_bc <- windows2016_adonis_function(win456_dist_bc, "Bray-Curtis", win456_meta, "456")
win456_adonis_uu <- windows2016_adonis_function(win456_dist_uu, "Unifrac-Unweighted", win456_meta, "456")
win456_adonis_wu <- windows2016_adonis_function(win456_dist_wu, "Unifrac-Weighted", win456_meta, "456")
win456_adonis_all <- rbind(win456_adonis_ds, win456_adonis_bc, win456_adonis_uu, win456_adonis_wu)
rm(win456_adonis_ds, win456_adonis_bc, win456_adonis_uu, win456_adonis_wu)
win567_adonis_ds <- windows2016_adonis_function(win567_dist_ds, "Dice-Sorensen", win567_meta, "567")
win567_adonis_bc <- windows2016_adonis_function(win567_dist_bc, "Bray-Curtis", win567_meta, "567")
win567_adonis_uu <- windows2016_adonis_function(win567_dist_uu, "Unifrac-Unweighted", win567_meta, "567")
win567_adonis_wu <- windows2016_adonis_function(win567_dist_wu, "Unifrac-Weighted", win567_meta, "567")
win567_adonis_all <- rbind(win567_adonis_ds, win567_adonis_bc, win567_adonis_uu, win567_adonis_wu)
rm(win567_adonis_ds, win567_adonis_bc, win567_adonis_uu, win567_adonis_wu)
windows2016_adonis_allMetrics_allWindows <- rbind(win456_adonis_all, win567_adonis_all)

# write_csv(windows2016_adonis_allMetrics_allWindows, 
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/windows2016_adonis_allMetrics_allWindows.csv")

rm(win456_adonis_all, win567_adonis_all, windows2016_adonis_allMetrics_allWindows)
## Window is significant main effect regardless of distance metric or window... 
## greatest proportion of variation at SITE captured in win456 group, using Weighted-Unifrac measure (almost 2x as other metrics)
## proportion variation for Windows not as sig as SITE for either windowGroup, but generally a bit higher for group "win567" in comparable distance metric
## Appears that SITE main effect is driving biggest fraction of community composition variation for "win456"; 
## Appears WINDOW is marginally higher than SITE for main effect variation in "win567" group
## caveats: different window overlaps!
## caveats: different sites - potentially different environmental conditions and arthropod community availability!!

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
           `Pr..F.` = round(`Pr..F.`, 5))
}

win456_bdisp_ds <- window2016_betadisper_function(win456_dist_ds, "Dice-Sorensen", win456_meta, "win456")
win456_bdisp_bc <- window2016_betadisper_function(win456_dist_bc, "Bray-Curtis", win456_meta, "win456")
win456_bdisp_uu <- window2016_betadisper_function(win456_dist_uu, "Unifrac-Unweighted", win456_meta, "win456")
win456_bdisp_wu <- window2016_betadisper_function(win456_dist_wu, "Unifrac-Weighted", win456_meta, "win456")

win567_bdisp_ds <- window2016_betadisper_function(win567_dist_ds, "Dice-Sorensen", win567_meta, "win567")
win567_bdisp_bc <- window2016_betadisper_function(win567_dist_bc, "Bray-Curtis", win567_meta, "win567")
win567_bdisp_uu <- window2016_betadisper_function(win567_dist_uu, "Unifrac-Unweighted", win567_meta, "win567")
win567_bdisp_wu <- window2016_betadisper_function(win567_dist_wu, "Unifrac-Weighted", win567_meta, "win567")

windows2016_bdisp_all <- rbind(win456_bdisp_ds, win456_bdisp_bc, win456_bdisp_uu, win456_bdisp_wu,
                               win567_bdisp_ds, win567_bdisp_bc, win567_bdisp_uu, win567_bdisp_wu)
## all beta dispersions for each SITE or WINDOW effect (for each distance metric and window group) is SIGNIFICANT! ... 
## need to treat this interpretation more carefully - if variance within groups (individ. samples relative to group centroid)...
## ..differs, then the group difference can be driven due to differences in dispersion, not because of differences in community represenation itself
## plotting these data with ellipses helps tease apart the trends... see next section!
# write_csv(windows2016_bdisp_all,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/windows2016_betaDispersion_allMetrics_allWindows.csv")

rm(win456_bdisp_ds, win456_bdisp_bc, win456_bdisp_uu, win456_bdisp_wu,
   win567_bdisp_ds, win567_bdisp_bc, win567_bdisp_uu, win567_bdisp_wu)



##################
## part 2c: ordinations
## going to plot a few different ways: 
## 1) individual metrics for each windowGroup
## 2) one large plot per WindowGroup for all 4 distance estimates

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
win456_orddata_bc <- windows2016_ordi_function(win456_dist_bc, "Bray-Curtis", win456_phydat, "win456", win456_meta)
win456_orddata_uu <- windows2016_ordi_function(win456_dist_uu, "Unweighted-UniFrac", win456_phydat, "win456", win456_meta)
win456_orddata_wu <- windows2016_ordi_function(win456_dist_wu, "Weighted-UniFrac", win456_phydat, "win456", win456_meta)
win567_orddata_ds <- windows2016_ordi_function(win567_dist_ds, "Dice-Sorensen", win567_phydat, "win567", win567_meta)
win567_orddata_bc <- windows2016_ordi_function(win567_dist_bc, "Bray-Curtis", win567_phydat, "win567", win567_meta)
win567_orddata_uu <- windows2016_ordi_function(win567_dist_uu, "Unweighted-UniFrac", win567_phydat, "win567", win567_meta)
win567_orddata_wu <- windows2016_ordi_function(win567_dist_wu, "Weighted-UniFrac", win567_phydat, "win567", win567_meta)

## for all plot types, set up these color/shape parameters
## palette to match windows from other plots!
win456_pal <- rev(c("#557d3f", "#B28C32", "#FCA68C"))
win567_pal <- rev(c("#267b97", "#557d3f", "#B28C32"))

## need to define shape point types; ensure overlapping Sites have same point value
win456_shapes <- c(15, 16, 17, 9, 10, 11, 12)
win567_shapes <- c(0, 1, 17, 4)
## note the orders here: the third values need to be the same so that "FOX" site...
## ... common in both 'win456' and 'win567' are using the same shape types

### 1) single plots per WindowGroup & Distance estimate
window2016_plotfunction <- function(ordinationdata, shapepal, colorpal){
  
  ordinationdata$Metric <- factor(ordinationdata$Metric, levels = c(
    "Dice-Sorensen", "Bray-Curtis", "Unweighted-UniFrac", "Weighted-UniFrac"))
  
  ggplot(data = ordinationdata,
         aes(x=Axis.1, y=Axis.2, color=Window, shape = Site)) +
    geom_point(aes(color = Window, shape = Site), size = 3) +
    stat_ellipse(aes(group = Window, color = Window), alpha = 0.5) +
    labs(x = unique(ordinationdata$Axis1lab), y = unique(ordinationdata$Axis2lab)) +
    theme_classic() +
    facet_wrap(~Metric) +
    scale_color_manual(values = colorpal) +
    scale_shape_manual(values = shapepal) +
    coord_fixed()
}

p_win456_ord_ds <- window2016_plotfunction(win456_orddata_ds, win456_shapes, win456_pal)
p_win456_ord_ds
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_ds.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_ds.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_ds.svg", width = 17, height = 15, unit = "cm")

p_win456_ord_bc <- window2016_plotfunction(win456_orddata_bc, win456_shapes, win456_pal)
p_win456_ord_bc
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_bc.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_bc.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_bc.svg", width = 17, height = 15, unit = "cm")

p_win456_ord_uu <- window2016_plotfunction(win456_orddata_uu, win456_shapes, win456_pal)
p_win456_ord_uu
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_uu.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_uu.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_uu.svg", width = 17, height = 15, unit = "cm")

p_win456_ord_wu <- window2016_plotfunction(win456_orddata_wu, win456_shapes, win456_pal)
p_win456_ord_wu
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_wu.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_wu.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_wu.svg", width = 17, height = 15, unit = "cm")

p_win567_ord_ds <- window2016_plotfunction(win567_orddata_ds, win567_shapes, win567_pal)
p_win567_ord_ds
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_ds.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_ds.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_ds.svg", width = 17, height = 15, unit = "cm")

p_win567_ord_bc <- window2016_plotfunction(win567_orddata_bc, win567_shapes, win567_pal)
p_win567_ord_bc
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_bc.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_bc.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_bc.svg", width = 17, height = 15, unit = "cm")

p_win567_ord_uu <- window2016_plotfunction(win567_orddata_uu, win567_shapes, win567_pal)
p_win567_ord_uu
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_uu.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_uu.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_uu.svg", width = 17, height = 15, unit = "cm")

p_win567_ord_wu <- window2016_plotfunction(win567_orddata_wu, win567_shapes, win567_pal)
p_win567_ord_wu
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_wu.png", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_wu.pdf", width = 17, height = 15, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_wu.svg", width = 17, height = 15, unit = "cm")


### 2) combined plots with different distance estimates for data with shared WindowGroups
ggarrange(p_win456_ord_ds, p_win456_ord_bc, p_win456_ord_uu, p_win456_ord_wu,
          common.legend = TRUE, nrow = 2, ncol = 2, align = "hv")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_allMetrics.png", width = 25, height = 25, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_allMetrics.pdf", width = 25, height = 25, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win456_allMetrics.svg", width = 25, height = 25, unit = "cm")


ggarrange(p_win567_ord_ds, p_win567_ord_bc, p_win567_ord_uu, p_win567_ord_wu,
          common.legend = TRUE, nrow = 2, ncol = 2, align = "hv")

# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_allMetrics.png", width = 25, height = 25, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_allMetrics.pdf", width = 25, height = 25, unit = "cm")
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/pcoa_win567_allMetrics.svg", width = 25, height = 25, unit = "cm")

rm(list = ls(pattern = "^win456_orddata"))
rm(list = ls(pattern = "^win567_orddata"))
rm(list = ls(pattern = "^p_win"))
rm(win567_pal, win456_pal, win456_shapes, win567_shapes, windows2016_orddate_all)

##################
## part 2d: pairwise adonis
## didn't ultimately use these values - shown here anyway in case someone wanted to run this test

### here we're testing the distances individually for each metadata variable (Window and Site)
window2016_pairwise_adonis_func <- function(distancevals, metric, windowgroup, metadata){
  tmp_pairwiseadonis_window <- 
    pairwise.adonis(distancevals, metadata$Window, p.adjust.m = "bonferroni") %>%
    mutate(Metric = metric, WindowGroup = windowgroup, Factor = "Window") %>%
    separate(col = pairs, into = c("Factor_A", "Factor_B"), sep = "vs") %>%
    arrange(p.adjusted)
  tmp_pairwiseadonis_site <- 
    pairwise.adonis(distancevals, metadata$Site, p.adjust.m = "bonferroni") %>%
    mutate(Metric = metric, WindowGroup = windowgroup, Factor = "Site") %>%
    separate(col = pairs, into = c("Factor_A", "Factor_B"), sep = "vs") %>%
    arrange(p.adjusted)
  rbind(tmp_pairwiseadonis_window, tmp_pairwiseadonis_site)
}

win456_pwadonis_ds <- window2016_pairwise_adonis_func(win456_dist_ds, "Dice-Sorensen", "win456", win456_meta)
win456_pwadonis_bc <- window2016_pairwise_adonis_func(win456_dist_bc, "Bray-Curtis", "win456", win456_meta)
win456_pwadonis_uu <- window2016_pairwise_adonis_func(win456_dist_uu, "Unweighted-Unifrac", "win456", win456_meta)
win456_pwadonis_wu <- window2016_pairwise_adonis_func(win456_dist_wu, "Weighted-Unifrac", "win456", win456_meta)
win567_pwadonis_ds <- window2016_pairwise_adonis_func(win567_dist_ds, "Dice-Sorensen", "win567", win567_meta)
win567_pwadonis_bc <- window2016_pairwise_adonis_func(win567_dist_bc, "Bray-Curtis", "win567", win567_meta)
win567_pwadonis_uu <- window2016_pairwise_adonis_func(win567_dist_uu, "Unweighted-Unifrac", "win567", win567_meta)
win567_pwadonis_wu <- window2016_pairwise_adonis_func(win567_dist_wu, "Weighted-Unifrac", "win567", win567_meta)

windows2016_pwadonis_all <- rbind(win456_pwadonis_ds, win456_pwadonis_bc, win456_pwadonis_uu, win456_pwadonis_wu,
                                  win567_pwadonis_ds, win567_pwadonis_bc, win567_pwadonis_uu, win567_pwadonis_wu)

rm(win456_pwadonis_ds, win456_pwadonis_bc, win456_pwadonis_uu, win456_pwadonis_wu,
   win567_pwadonis_ds, win567_pwadonis_bc, win567_pwadonis_uu, win567_pwadonis_wu)
# write_csv(windows2016_pwadonis_all,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/windows2016_pairwiseAdonis_all.csv")

## if we want to build a heatmap of these values, we're going to need to rebuild a symmetric matrix for both factors (Site and Window)
## how about using the distances and creating boxplots like in QIIME's beta group significance?
## couldn't for the life of me figure out how to get the pairwise distance calculations to work right
## this got me close, but still didn't generate all the pairwise comparisons needed (there were a few empty box plots?)
# tmp_meta <- win567_meta
# tmp_samplist <- tmp_meta %>% select(SampleID) %>% pull()
# tmp_subset_dist <- dist_subset(win567_dist_ds, tmp_samplist)
# tmp_dist_names <- labels(tmp_subset_dist)
# tmp_dist_site_perSample <- tmp_meta %>% filter(SampleID %in% tmp_dist_names) %>% select(Site) %>% pull()
# tmp_pairwise <- dist_groups(tmp_subset_dist, tmp_dist_site_perSample)
# ggplot(tmp_pairwise, aes(x=Label, y=Distance)) + 
#   geom_boxplot(outlier.colour = NA) + 
#   geom_point(position = position_jitter(width = 0.25), alpha = 0.5) + 
#   facet_grid(Group1 ~ .) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


##################
## part 2e: pairwise distance boxplots
## ended up running a separate series of calculations using QIIME in a shell script...
## available here: https://raw.githubusercontent.com/devonorourke/nhguano/master/scripts/shell_scripts/get_beta_group_significance_windows2016data.sh
## the result is a pair of files (win456 and win567) containing all pairwise comparisons of distances for each of the four distance metrics

## note that the original tree .nwk file (imported as 'tree' in R) was already available for use in the shell script above
## First, this is what was exported for use in that shell script:
#### OTUtable for win456 data
win456_otutable <- as.data.frame(as.matrix(otu_table(win456_phydat)))
win456_otutable$`OTU id` <- row.names(win456_otutable)
win456_otutable <- win456_otutable %>% 
  select(`OTU id`, everything())
write_tsv(win456_otutable, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win456_OTUtable.tsv")
rm(win456_otutable)
#### OTUtable for win567 data
win567_otutable <- as.data.frame(as.matrix(otu_table(win567_phydat)))
win567_otutable$`OTU id` <- row.names(win567_otutable)
win567_otutable <- win567_otutable %>% 
  select(`OTU id`, everything())
write_tsv(win567_otutable, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win567_OTUtable.tsv")
rm(win567_otutable)
## Metadata for win456 and win567 datasets
win456_meta_alt <- win456_meta %>% 
  mutate(Window = paste0("window_", Window))
write_tsv(win456_meta_alt, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win456meta_forQIIME.txt")
win567_meta_alt <- win567_meta %>% 
  mutate(Window = paste0("window_", Window))
write_tsv(win567_meta_alt, path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/win567meta_forQIIME.txt")

## cleanup:
rm(list = ls(pattern = "^win"))

##################
## imports data from output of `.sh` script...
## imports metadata and merges with distance data
## calculates the bottom 5% distance value for a given Metric among within site points 

betagroupsig_plotdat_function <- function(dist_url, meta_url, metric){
  tmp_betagroupsig_allmetrics <- 
    read_delim(delim = "\t", dist_url) %>% 
    mutate(MetricLabel = case_when(Metric == "dice" ~ 'Dice-Sorensen',
                                   Metric == "bray" ~ 'Bray-Curtis',
                                   Metric == "uuni" ~ 'Unweighted-UniFrac',
                                   Metric == "wuni" ~ 'Weighted-UniFrac'),
           Group1 = str_replace(Group1, "window_", "Window "),
           Group2 = str_replace(Group2, "window_", "Window ")) %>% 
    filter(FactorGroup == "site" & Metric == metric)
  
  win567_meta_tmp <- read_tsv(meta_url) %>% select(SampleID, Window) %>% mutate(Window = str_replace(Window, "window_", "Window "))
  
  tmp_betagroupsig_allmetrics <- 
    merge(tmp_betagroupsig_allmetrics, win567_meta_tmp, by.x = "SubjectID1", by.y = "SampleID", all.x = TRUE) %>% 
    rename(WindowSubject1 = "Window")
  
  tmp_betagroupsig_allmetrics <- 
    merge(tmp_betagroupsig_allmetrics, win567_meta_tmp, by.x = "SubjectID2", by.y = "SampleID", all.x = TRUE) %>% 
    rename(WindowSubject2 = "Window")
  
  tmp_betagroupsig_allmetrics <- 
    tmp_betagroupsig_allmetrics %>% 
    mutate(WindowCompareTest = (WindowSubject1 == WindowSubject2),
           WindowLabel = case_when(WindowCompareTest == TRUE ~ paste0('within ',WindowSubject1),
                                   WindowCompareTest == FALSE ~ paste0(WindowSubject1, ' and ', WindowSubject2)))
  
  tmp_betagroupsig_allmetrics <- tmp_betagroupsig_allmetrics %>% 
    #group_by(grp = paste(pmax(SubjectID1, SubjectID2), pmin(SubjectID1, SubjectID2), sep = "_")) %>% 
    #slice(1) %>% 
    #ungroup() %>% 
    #select(-grp) %>% 
    mutate(WindowLabel = case_when(WindowLabel == "Window 6 and Window 4" ~ "Window 4 and Window 6",
                                   WindowLabel == "Window 6 and Window 5" ~ "Window 5 and Window 6",
                                   WindowLabel == "Window 5 and Window 4" ~ "Window 4 and Window 5",
                                   WindowLabel == "Window 7 and Window 5" ~ "Window 5 and Window 7",
                                   WindowLabel == "Window 7 and Window 6" ~ "Window 6 and Window 7",
                                   TRUE ~ as.character(.$WindowLabel)),
           SiteLabel = paste0(Group1, "_", Group2))
  
  dist_sumry <- tmp_betagroupsig_allmetrics %>% 
    group_by(Metric) %>% 
    summarise(DistanceBar = quantile(Distance, 0.05, q=0.05))
  
  merge(tmp_betagroupsig_allmetrics, dist_sumry)
}

win456_betagroupsig_plotdat_dice <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
  'dice')

win456_betagroupsig_plotdat_bray <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
  'bray')

win456_betagroupsig_plotdat_uuni <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
  'uuni')

win456_betagroupsig_plotdat_wuni <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win456_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win456meta_forQIIME.txt',
  'wuni')


win567_betagroupsig_plotdat_dice <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
  'dice')

win567_betagroupsig_plotdat_bray <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
  'bray')

win567_betagroupsig_plotdat_uuni <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
  'uuni')

win567_betagroupsig_plotdat_wuni <- betagroupsig_plotdat_function(
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/beta_group_sig/win567_betasig_allMetrics_data.tsv.gz',
  'https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/win567meta_forQIIME.txt',
  'wuni')


## generate data for plotting function that will plot each Metric (for each windowgroup) separately
## use value in 'DistanceBar' to highlight the bottom 5% of distance values per Metric...
## helps illustrate which Sites are more similar to each other within each Window

## droplist for the grid to focus on just one of the two reversible pairs
dropsitelabels_win567 <- c("CNB_CHI", "FOX_CHI", "MTV_CHI", "FOX_CNB", "MTV_CNB", "MTV_FOX")

dropsitelabels_win456 <- c("EPS_CNA", "FOX_CNA", "HOL_CNA", "HOP_CNA", "MAP_CNA", "PEN_CNA",
                           "FOX_EPS", "HOL_EPS", "HOP_EPS", "MAP_EPS", "PEN_EPS",
                           "HOL_FOX", "HOP_FOX", "MAP_FOX", "PEN_FOX",
                           "HOP_HOL", "MAP_HOL", "PEN_HOL",
                           "MAP_HOP", "PEN_HOP",
                           "PEN_MAP")

betagroupsig_plotfunction_win567 <- function(plotdat, dropsitelabels){
  
  plotdat$WindowLabel <-
    factor(plotdat$WindowLabel, 
           levels = c("within Window 5", "within Window 6", "within Window 7",
                      "Window 5 and Window 6", "Window 6 and Window 7", "Window 5 and Window 7"))
  ggplot(data = plotdat %>% filter(!SiteLabel %in% dropsitelabels),
         aes(x=Group1, 
             y=Distance,
             color=WindowLabel)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1),
               alpha = 0.5) +
    geom_boxplot(position = position_dodge(1), outlier.shape=NA,
                 fill="white", alpha=0.5) +
    #geom_hline(yintercept = unique(plotdat$DistanceBar), color="gray20", alpha = 0.4) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.75)) +
    scale_color_viridis_d(option = 'magma', begin = 0, end = 0.9) +
    labs(x="", y="Distance\n", fill="Sample pair") +
    facet_grid(Group2 ~ Group1, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(panel.spacing.y = unit(2, "lines"))
}

# p_betasig_win567_dice <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_dice, dropsitelabels_win567)
# p_betasig_win567_dice
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_dice.png", height = 15, width = 22, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_dice.pdf", height = 15, width = 22, units = 'cm', dpi=300)
# 
# p_betasig_win567_bray <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_bray, dropsitelabels_win567)
# p_betasig_win567_bray
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_bray.png", height = 15, width = 22, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_bray.pdf", height = 15, width = 22, units = 'cm', dpi=300)
# 
# p_betasig_win567_uuni <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_uuni, dropsitelabels_win567)
# p_betasig_win567_uuni
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_uuni.png", height = 15, width = 20, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_uuni.pdf", height = 15, width = 20, units = 'cm', dpi=300)
# 
# p_betasig_win567_wuni <- betagroupsig_plotfunction_win567(win567_betagroupsig_plotdat_wuni, dropsitelabels_win567)
# p_betasig_win567_wuni
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_wuni.png", height = 15, width = 20, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_wuni.pdf", height = 15, width = 20, units = 'cm', dpi=300)
# 
# ggpubr::ggarrange(p_betasig_win567_dice, p_betasig_win567_bray, p_betasig_win567_uuni, p_betasig_win567_wuni,
#                   common.legend = TRUE, align = "hv", labels = c("A", "B", "C", "D"))
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_allMetrics.png", height = 20, width = 20, units = 'cm', dpi=150)
# ggsave(filename = "~/Documents/nau_projects/guano/NHguano_redux/new_figures/win567_betagroupsig_allMetrics.pdf", height = 20, width = 20, units = 'cm', dpi=300)

rm(list = ls(pattern = "^p_betasig_win567"))
rm(list = ls(pattern = "^win567"))


betagroupsig_plotfunction_win456 <- function(plotdat, dropsitelabels){
  
  plotdat$WindowLabel <-
    factor(plotdat$WindowLabel, 
           levels = c("within Window 4", "within Window 5", "within Window 6",
                      "Window 4 and Window 5", "Window 5 and Window 6", "Window 4 and Window 6"))
  ggplot(data = plotdat %>% filter(!SiteLabel %in% dropsitelabels),
         aes(x=Group1, 
             y=Distance,
             color=WindowLabel)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1),
               alpha = 0.5) +
    geom_boxplot(position = position_dodge(1), outlier.shape=NA,
                 fill="white", alpha=0.5) +
    #geom_hline(yintercept = unique(plotdat$DistanceBar), color="gray20", alpha = 0.4) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.75)) +
    scale_color_viridis_d(option = 'magma', begin = 0, end = 0.9) +
    labs(x="", y="Distance\n", fill="Sample pair") +
    facet_grid(Group2 ~ Group1, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(panel.spacing.y = unit(2, "lines"))
}

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


########################################
## part 5 - evaluate between between sites in single window of sampling season between 2015 and 2016 years
########################################

# ##################
# ## part 5a: In comparing between 2015 qnd 2016... Which sites make sense? Which months?
otu_table_long_filtd_wMeta %>%
  group_by(Site, Year, Window) %>%
  summarise(nSamples = n_distinct(SampleID)) %>%
  pivot_wider(names_from = "Year", values_from = "nSamples") %>%
  filter(!is.na(`2015`)) %>%
  filter(!is.na(`2016`)) %>%
  arrange(Window, Site)

## focusing on just 3 sites at single samplign window: COR and HOP and MAP
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
multiyear_otu_table_wide_norm <- SRS(multiyear_otu_table_wide_raw, multiyear_Cmin)
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
write_csv(multiyear_otu_table_long_norm_wMeta_andTaxa, 
          path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/multiyear_otu_table_long_norm_wMeta_andTaxa.csv")

# 
# ##################
# ## part 5c: diversity estimates using SRS output calculated with Phyloseq to incorporate tree info for both alpha and distance calcs
# ## create phyloseq object
# ## import OTU data in phyloseq format
multiyear_phyOTU <- otu_table(multiyear_otu_table_wide_norm, taxa_are_rows = TRUE)
## bundle OTU table, metadata, and tree info into single phyloseq object
multiyear_phydat <- phyloseq(multiyear_phyOTU, multiyear_meta)
multiyear_phydat <- merge_phyloseq(multiyear_phydat, tree)
rm(multiyear_phyOTU, multiyear_otu_table_wide_norm)
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
## chi-squared = 3.5474, df = 5, p-value = 0.6162

# ## cleanup
rm(list=ls(pattern="^kw_"))

### quick sanity check - do most of our samples contain at least a handful of OTUs? 
multiyear_alpha_df_long %>% filter(Metric == "SR" & Alpha_value > 4) %>% nrow()
## 68 of 81 (~ 84 %) samples have at least 5 OTUs
multiyear_alpha_df_long %>% filter(Metric == "SR" & Alpha_value > 9) %>% nrow()
## 36 of 81 (~ 47 %) samples have at least 10 OTUs

## add grouper term to simplify plotting
multiyear_alpha_df_long$SiteYear <- paste(multiyear_alpha_df_long$Site, multiyear_alpha_df_long$Year, sep = "\n")

## visualize these data in boxplots, grouping data by sampling site, with it's pair of years next to each other
### plotting by individual Metrics, then stitching together into multifaceted plot
multiyear_alpha_df_long$Metric <- 
  factor(multiyear_alpha_df_long$Metric, levels = c("SR", "H", "PD"))

multiyear_alphaviz_function <- function(metric, yaxislab){
  ggplot(data = multiyear_alpha_df_long %>% filter(Metric == metric), 
         aes(x=SiteYear, y=Alpha_value, fill=SiteYear)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
    geom_boxplot(outlier.color = NA, color="gray50", fill="white", width = 0.8, alpha=0.25) +
    facet_grid(Metric~Site, scales = "free") +
    labs(x="", y=yaxislab) +
    theme_classic() +
    theme(strip.text = element_blank(),
          legend.position = "none")
}
p_multiyear_alpha_sr <- multiyear_alphaviz_function("SR", "Observed OTUs")
p_multiyear_alpha_sr
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_sr.png", height=10, width=15, units="cm", dpi=150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_sr.pdf", height=10, width=15, units="cm", dpi=300)

p_multiyear_alpha_h <- multiyear_alphaviz_function("H", "H")
p_multiyear_alpha_h
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_h.png", height=10, width=15, units="cm", dpi=150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_h.pdf", height=10, width=15, units="cm", dpi=300)

p_multiyear_alpha_pd <- multiyear_alphaviz_function("PD", "Faith's PD diversity")
p_multiyear_alpha_pd
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_pd.png", height=10, width=15, units="cm", dpi=150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_pd.pdf", height=10, width=15, units="cm", dpi=300)

ggarrange(p_multiyear_alpha_sr, p_multiyear_alpha_h, p_multiyear_alpha_pd,
          labels = c("A", "B", "C"), align = "hv", ncol=1, nrow = 3)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_allMetrics.png", height=20, width=15, units = "cm", dpi=150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_alpha_boxplot_allMetrics.pdf", height=20, width=15, units = "cm", dpi=300)

## cleanup:
rm(list=ls(pattern = "^p_multiyear"))
rm(multiyear_alpha_df, multiyear_alpha_df_long)

# ### doesn't appear that there's any story in terms of differences in #OTUs detected per sample among these year/site groupings
# ## but are there any compositional changes?

# ##################
# ## part 2b: community composition
# ## 2bi. run PERMANOVA (via Adonis) for each distance method to evaluate if samples have different centroids among Sites or Years

# ## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method
multiyear_dist_ds <- phyloseq::distance(multiyear_phydat, "bray", binary = TRUE)
multiyear_dist_bc <- phyloseq::distance(multiyear_phydat, "bray", binary = FALSE)
multiyear_dist_uu <- phyloseq::distance(multiyear_phydat, "unifrac", weighted=FALSE)
multiyear_dist_wu <- phyloseq::distance(multiyear_phydat, "wunifrac")

multiyear_adonis_function <- function(distanceData, distanceMetric){
  adonis_tmp <- adonis2(distanceData ~ Site * Year, data = multiyear_meta, )
  data.frame(adonis_tmp) %>% mutate(Metric = distanceMetric,
                                    Class = row.names(.))
}

multiyear_adonis_ds <- multiyear_adonis_function(multiyear_dist_ds, "Dice-Sorensen")
multiyear_adonis_bc <- multiyear_adonis_function(multiyear_dist_bc, "Bray-Curtis")
multiyear_adonis_uu <- multiyear_adonis_function(multiyear_dist_uu, "Unifrac Unweighted")
multiyear_adonis_wu <- multiyear_adonis_function(multiyear_dist_wu, "Unifrac Weighted")

multiyear_adonis_all <- rbind(multiyear_adonis_ds, multiyear_adonis_bc, multiyear_adonis_uu, multiyear_adonis_wu)
rm(multiyear_adonis_ds, multiyear_adonis_bc, multiyear_adonis_uu, multiyear_adonis_wu)
## what we find is that when we weight data with abundance information, SITE, but not YEAR is significant
## when data abundance information is not weighted, both SITE and YEAR main effects are significant
## interaction term for SITE*YEAR is significant for all 4 metrics, so only some Site differences are significant... which of the 3 are different from each other??

multiyear_adonis_all <- 
  multiyear_adonis_all[,c(7,1,2,3,4,5,6)] %>% 
  rename(Fmodel = `F`) %>% 
  mutate(SumOfSqs = round(SumOfSqs, 3),
         R2 = round(R2, 3),
         `Pr..F.` = round(`Pr..F.`, 4),
         Fmodel = round(Fmodel, 3))
# write.csv(multiyear_adonis_all, quote=FALSE, row.names = FALSE,
#           file = "~/Documents/nau_projects/guano/NHguano_redux/new_data/multiyear_adonis_allMetrics.csv")



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

multiyear_bdisp_bc <- multiyear_betadisper_function(multiyear_dist_bc, "Bray-Curtis")
multiyear_bdisp_bc
## non significant for either GroupFactor: Site p = 0.859; Year p = 0.620

multiyear_bdisp_uu <- multiyear_betadisper_function(multiyear_dist_uu, "Unifrac Unweighted")
multiyear_bdisp_uu
## non significant for either GroupFactor: Site p = 0.338; Year p = 0.360

multiyear_bdisp_wu <- multiyear_betadisper_function(multiyear_dist_wu, "Unifrac Weighted")
multiyear_bdisp_wu
## non significant for either GroupFactor: Site p = 0.448; Year p = 0.674

multiyear_bdisp_all <- rbind(multiyear_bdisp_ds, multiyear_bdisp_bc, multiyear_bdisp_uu, multiyear_bdisp_wu)
rm(multiyear_bdisp_ds, multiyear_bdisp_bc, multiyear_bdisp_uu, multiyear_bdisp_wu)
# write_csv(multiyear_bdisp_all,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/multiyear_betadisper_allMetrics.csv")


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
multiyear_orddata_bc <- multiyear_ordi_function(multiyear_dist_bc, "Bray-Curtis")
multiyear_orddata_uu <- multiyear_ordi_function(multiyear_dist_uu, "Unweighted-Unifrac")
multiyear_orddata_wu <- multiyear_ordi_function(multiyear_dist_wu, "Weighted-Unifrac")

multiyear_orddata_all <- rbind(multiyear_orddata_ds, multiyear_orddata_bc, multiyear_orddata_uu, multiyear_orddata_wu)
rm(multiyear_orddata_ds, multiyear_orddata_bc, multiyear_orddata_uu, multiyear_orddata_wu)

## plot with metadata
#### use shapes to distinguish years, and color to distinguish sites
## plot individually, then stitch together
multiyear_ordplot_function <- function(inputdata){
  inputdata$Year = as.character(inputdata$Year)
  ggplot(data = inputdata,
         aes(x=Axis.1, y=Axis.2, color = Site, shape = Year)) +
    geom_point(size = 3, alpha=0.8) +
    stat_ellipse(aes(group = Site), alpha = 0.5) +
    labs(x = unique(inputdata$Axis1lab),
         y = unique(inputdata$Axis2lab)) +
    scale_color_manual(values = c('dodgerblue', '#ffb07c', '380282')) +
    scale_shape_manual(values = c(16, 17)) +
    coord_fixed() +
    theme_classic() +
    theme(legend.position = "top", legend.box = "horizontal", legend.margin = margin())
}

p_multiyear_ord_ds <- multiyear_ordplot_function(multiyear_orddata_all %>% filter(Metric == "Dice-Sorensen"))
p_multiyear_ord_ds
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_dice.png", height = 12, width = 12, units = "cm", dpi = 150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_dice.pdf", height = 12, width = 12, units = "cm", dpi = 300)

p_multiyear_ord_bc <- multiyear_ordplot_function(multiyear_orddata_all %>% filter(Metric == "Bray-Curtis"))
p_multiyear_ord_bc
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_bray.png", height = 12, width = 12, units = "cm", dpi = 150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_bray.pdf", height = 12, width = 12, units = "cm", dpi = 300)

p_multiyear_ord_uu <- multiyear_ordplot_function(multiyear_orddata_all %>% filter(Metric == "Unweighted-Unifrac"))
p_multiyear_ord_uu
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_uuni.png", height = 12, width = 12, units = "cm", dpi = 150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_uuni.pdf", height = 12, width = 12, units = "cm", dpi = 300)

p_multiyear_ord_wu <- multiyear_ordplot_function(multiyear_orddata_all %>% filter(Metric == "Weighted-Unifrac"))
p_multiyear_ord_wu
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_wuni.png", height = 12, width = 12, units = "cm", dpi = 150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pcoa_wuni.pdf", height = 12, width = 12, units = "cm", dpi = 300)

legend_multiyear <- get_legend(p_multiyear_ord_ds)

tmp_multiyear_plot <- plot_grid(
  p_multiyear_ord_ds + theme(legend.position = "none"),
  p_multiyear_ord_bc + theme(legend.position = "none"),
  p_multiyear_ord_uu + theme(legend.position = "none"),
  p_multiyear_ord_wu + theme(legend.position = "none"),
  labels = c("A", "B", "C", "D"),
  nrow = 2, ncol = 2,
  align = "hv")

plot_grid(legend_multiyear, tmp_multiyear_plot,
          nrow = 2, ncol = 1,
          rel_heights = c(.1, 1))
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_ordinations_all.png", width = 20, height = 21, units = 'cm', dpi = 150)
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_ordinations_all.pdf", width = 20, height = 21, units = 'cm', dpi = 300)

## cleanup
rm(p_multiyear_ord_ds, p_multiyear_ord_bc, p_multiyear_ord_uu, p_multiyear_ord_wu)

# 
# ##################
# ## part 2d: pairwise Adonis to determine what, if any, community comopsition groups are different
# ## generate long form of all values for supplementary table
multiyear_pairwise_adonis_func <- function(distancevals, metric){
  pairwise.adonis(distancevals, paste0(multiyear_meta$Site, "-", multiyear_meta$Year), p.adjust.m = "bonferroni") %>%
    #pairwise.adonis(distancevals, paste0(multiyear_meta$Site, "-", multiyear_meta$Year), p.adjust.m = "BH") %>%
    mutate(Metric = metric) %>%
    separate(col = pairs, into = c("Window_A", "Window_B"), sep = "vs") %>%
    arrange(p.adjusted)
}

multiyear_pwadonis_ds <- multiyear_pairwise_adonis_func(multiyear_dist_ds, "Dice-Sorensen")
multiyear_pwadonis_bc <- multiyear_pairwise_adonis_func(multiyear_dist_bc, "Bray-Curtis")
multiyear_pwadonis_uu <- multiyear_pairwise_adonis_func(multiyear_dist_uu, "Unweighted-Unifrac")
multiyear_pwadonis_wu <- multiyear_pairwise_adonis_func(multiyear_dist_wu, "Weighted-Unifrac")
multiyear_pwadonis_all <- rbind(multiyear_pwadonis_ds, multiyear_pwadonis_bc, multiyear_pwadonis_uu, multiyear_pwadonis_wu)
rm(multiyear_pwadonis_ds, multiyear_pwadonis_bc, multiyear_pwadonis_uu, multiyear_pwadonis_wu)
# write_csv(multiyear_pwadonis_all,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/multiyear_pairwiseAdonis_all.csv")

## abundance unweighted metrics agree with each other: between YEAR differences in composition for COR only, not for HOP or MAP
## abundances weighted metrics agree with each other: between YEAR differences for MAP only; not for HOP nor for COR
## so the weird thing here is that HOP has no year over year composition differences by any metric, ...
## but COR only does if you use abundance weighted metrics, while MAP only does if you use unweighted metrics!
## note that the same trends are observed regardless of whether or not phylogenetic info is used - it just depends on whether abundances are included
## ALSO note that if we switch up the multiple significance test measure (say from bonferroni to benjamini-hochberg), the trends are similar, but not identical:
## same: HOP has no p-value < 0.05 for any metric
## same: COR is sig for unweighted, not sig for weighted
## diff: MAP is significant for dice, as well as bray/wuni... just not for uuni. 


## plotting these pairwise pvalues:
multiyear_pairwise_adonis_heatmap_function <- function(metric){
  tmp <- multiyear_pwadonis_all %>% 
    filter(Metric == metric) %>% 
    select(Window_A, Window_B, p.adjusted) %>% 
    dcast(Window_A ~ Window_B, value.var = 'p.adjusted')
  tmp_pwadonis_df <- tmp[c(3,1,4,5,2),c(1,2,5,6,3,4)]
  row.names(tmp_pwadonis_df) <- tmp_pwadonis_df$Window_A
  tmp_pwadonis_df$Window_A <- NULL
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
  tmp_pwadonis_df_resolved$Window_A <- row.names(tmp_pwadonis_df_resolved)
  tmp_pwadonis_df_plotdat <- 
    tmp_pwadonis_df_resolved %>% 
    pivot_longer(-Window_A, names_to = "Window_B", values_to = "pvalue", values_drop_na = TRUE) %>% 
    mutate(Metric = metric) %>% 
    mutate(Window_A = gsub("-", "\n", Window_A),
           Window_B = gsub("-", "\n", Window_B))
  
  tmp_pwadonis_df_plotdat <- 
    tmp_pwadonis_df_plotdat %>% 
    mutate(DropStatus = case_when(
      Window_A == "COR\n2016" & Window_B == "COR\n2015" ~ "drop",
      Window_A == "HOP\n2015" & Window_B == "COR\n2015" ~ "drop",
      Window_A == "HOP\n2016" & Window_B == "COR\n2015" ~ "drop",
      Window_A == "MAP\n2015" & Window_B == "COR\n2015" ~ "drop",
      Window_A == "MAP\n2016" & Window_B == "COR\n2015" ~ "drop",
      Window_A == "HOP\n2015" & Window_B == "COR\n2016" ~ "drop",
      Window_A == "HOP\n2016" & Window_B == "COR\n2016" ~ "drop",
      Window_A == "MAP\n2015" & Window_B == "COR\n2016" ~ "drop",
      Window_A == "MAP\n2016" & Window_B == "COR\n2016" ~ "drop",
      Window_A == "HOP\n2016" & Window_B == "HOP\n2015" ~ "drop",
      Window_A == "MAP\n2015" & Window_B == "HOP\n2015" ~ "drop",
      Window_A == "MAP\n2016" & Window_B == "HOP\n2015" ~ "drop",
      Window_A == "MAP\n2015" & Window_B == "HOP\n2016" ~ "drop",
      Window_A == "MAP\n2016" & Window_B == "HOP\n2016" ~ "drop",
      Window_A == "MAP\n2016" & Window_B == "MAP\n2015" ~ "drop",
      TRUE ~ "keep"))
  
  
  tmp_pwadonis_df_plotdat$Window_A <- 
    factor(tmp_pwadonis_df_plotdat$Window_A,
           levels=c("COR\n2015", "COR\n2016", "HOP\n2015", "HOP\n2016", "MAP\n2015", "MAP\n2016"))
  
  tmp_pwadonis_df_plotdat$Window_B <-
    factor(tmp_pwadonis_df_plotdat$Window_B,
           levels=c("COR\n2015", "COR\n2016", "HOP\n2015", "HOP\n2016", "MAP\n2015", "MAP\n2016"))
  
  ggplot(tmp_pwadonis_df_plotdat %>% filter(DropStatus == "keep"),
         aes(Window_B, Window_A, fill=pvalue, label=pvalue)) +
    geom_tile() +
    geom_text() +
    coord_fixed() +
    scale_fill_viridis_c(option = "viridis", direction = -1, alpha = 0.5, begin = 0, end = 1, limits = c(0,1), breaks = c(0.25, 0.5, 0.75)) +
    facet_wrap(~Metric, nrow=2) +
    labs(x = "", y = "", fill = "Bonferonni-adjusted\np-value") +
    theme_classic() +
    theme(legend.position = "top")
  
}

p_multiyear_pwadonis_dice <- multiyear_pairwise_adonis_heatmap_function("Dice-Sorensen")
p_multiyear_pwadonis_dice
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_dice.png", width = 10, height=12, units = "cm", dpi=150)
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_dice.pdf", width = 10, height=12, units = "cm", dpi=300)

p_multiyear_pwadonis_bray <- multiyear_pairwise_adonis_heatmap_function("Bray-Curtis")
p_multiyear_pwadonis_bray
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_bray.png", width = 10, height=12, units = "cm", dpi=150)
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_bray.pdf", width = 10, height=12, units = "cm", dpi=300)


p_multiyear_pwadonis_uuni <- multiyear_pairwise_adonis_heatmap_function("Unweighted-Unifrac")
p_multiyear_pwadonis_uuni
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_uuni.png", width = 10, height=12, units = "cm", dpi=150)
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_uuni.pdf", width = 10, height=12, units = "cm", dpi=300)

p_multiyear_pwadonis_wuni <- multiyear_pairwise_adonis_heatmap_function("Weighted-Unifrac")
p_multiyear_pwadonis_wuni
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_wuni.png", width = 10, height=12, units = "cm", dpi=150)
ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pwadonis_heatmap_wuni.pdf", width = 10, height=12, units = "cm", dpi=300)

### stitch together into one faceted plot:
ggarrange(p_multiyear_pwadonis_dice, p_multiyear_pwadonis_bray, p_multiyear_pwadonis_uuni, p_multiyear_pwadonis_wuni,
          labels = c("A", "B", "C", "D", align="hv", nrow=2, ncol=1), common.legend=TRUE)

# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pairwiseAdonis_heatmap_pvalues_all.png", width = 22, height = 24, units = 'cm', dpi=150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_pairwiseAdonis_heatmap_pvalues_all.pdf", width = 22, height = 16, units = 'cm', dpi=300)

## cleanup
rm(multiyear_OTUtable_wide_raw, legend_multiyear, tmp_multiyear_plot)
rm(list=ls(pattern="^p_multiyear"))

### last part: 
## identifying what particular taxa are driving these differences in years for MAP (for bray/wuni) or COR (for dice/uuni)...
## consider working backwards: find the particular OTUs that are different, but then figure out if these are redundant at all among other taxa present (say at Family level?)

#multiyear_otu_table_long_norm_wMeta_andTaxa <- read_csv("~/Documents/nau_projects/guano/NHguano_redux/new_data/multiyear_otu_table_long_norm_wMeta_andTaxa.csv")

## for each site/year group, what are the topN OTUs we detect... are these the same in most site/year groups?
## first, how many samples are in each Site/Year group?
multiyear_sampsPerGroup <- multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  mutate(GrouperTerm = paste0(Site, "-", Year)) %>% 
  group_by(GrouperTerm) %>% 
  summarise(nSamplesInGroup = n_distinct(SampleID))

## next, which OTUs in each Site/Group are detected in at least N% of samples (set here at 33%, so 1/3 of samples...) PER SITE/YEAR group
multiyear_frequentOTUs_bySiteYear_list <- multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  mutate(GrouperTerm = paste0(Site, "-", Year)) %>% 
  group_by(GrouperTerm, OTUalias) %>% 
  summarise(nSamplesWithOTU = n()) %>% 
  merge(., multiyear_sampsPerGroup, by="GrouperTerm") %>% 
  mutate(pSamplesWithOTU = nSamplesWithOTU/nSamplesInGroup) %>% 
  ungroup() %>% 
  group_by(GrouperTerm) %>% 
  slice_max(order_by = pSamplesWithOTU, n=10) %>%
  filter(pSamplesWithOTU > 0.32) %>% 
  ungroup() %>% 
  distinct(OTUalias)

## using just these OTUs, determine how many samples had each of these OTUs for a given Site/Year group
multiyear_frequentOTUs_bySiteYear_matrix <- multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  mutate(GrouperTerm = paste0(Site, "-", Year)) %>% 
  filter(OTUalias %in% multiyear_frequentOTUs_bySiteYear_list$OTUalias) %>% 
  group_by(GrouperTerm, OTUalias) %>% 
  summarise(nSamplesWithOTU = n_distinct(SampleID)) %>% 
  merge(., multiyear_sampsPerGroup, by="GrouperTerm") %>% 
  mutate(pSamplesWithOTU = nSamplesWithOTU/nSamplesInGroup) %>% 
  ungroup() %>% 
  group_by(GrouperTerm) %>% 
  #slice_max(order_by = pSamplesWithOTU, n=10) %>%
  #mutate(printVal = paste0(nSamplesWithOTU, " / ", nSamplesInGroup)) %>% 
  #select(printVal, GrouperTerm, OTUalias) %>% 
  #pivot_wider(names_from="GrouperTerm", values_from="printVal")
  select(pSamplesWithOTU, GrouperTerm, OTUalias) %>% 
  mutate(pSamplesWithOTU = round(pSamplesWithOTU, 3)) %>% 
  pivot_wider(names_from="GrouperTerm", values_from="pSamplesWithOTU")

## add in taxa information to these OTUs
multiyear_frequentOTUs_bySiteYear_matrix_wTaxa <- 
  multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  filter(OTUalias %in% multiyear_frequentOTUs_bySiteYear_matrix$OTUalias) %>% 
  distinct(OTUalias, Class, Order, Family, Genus, Species) %>% 
  merge(., multiyear_frequentOTUs_bySiteYear_matrix)

# write_csv(multiyear_frequentOTUs_bySiteYear_matrix_wTaxa,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/multiyear_frequentOTUs_bySiteYear_matrix_wTaxa.csv")

### plot this as heatmap showing the fraction of these OTUs detected 
multiyear_frequentOTUs_bySiteYear_long_wTaxa <- multiyear_frequentOTUs_bySiteYear_matrix_wTaxa %>% 
  pivot_longer(cols = c("COR-2015", "COR-2016", "HOP-2015", "HOP-2016", "MAP-2015", "MAP-2016"), 
               values_to = "pSamplesWithOTU", names_to = "GrouperTerm") %>%
  filter(pSamplesWithOTU > 0) %>% 
  mutate(GenusLabel = ifelse(is.na(Genus), paste0(OTUalias, " f. ", Family), Genus),
         SpeciesLabel = ifelse(is.na(Species), paste0(OTUalias, " g. ", Genus), Species),
         SpeciesLabel = ifelse(grepl("g. NA$", SpeciesLabel), GenusLabel, SpeciesLabel),
         Site = substr(GrouperTerm, 1, 3))

ggplot(multiyear_frequentOTUs_bySiteYear_long_wTaxa, 
       aes(x=SpeciesLabel, y=GrouperTerm, fill=pSamplesWithOTU)) +
  geom_tile(color="black") +
  scico::scale_fill_scico(begin = 0.2, breaks = c(0.25, 0.5, 0.75), limits = c(0,1)) +
  facet_grid(Site ~ Order, scales = "free", space="free") +
  theme_classic() +
  labs(x = '', y='', fill='fraction\nsamples\nwith taxa\ndetected') +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        strip.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
        strip.text.y = element_blank())

# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_frequentOTUs_heatmap.png", height=10, width=15, units="cm", dpi=150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_frequentOTUs_heatmap.svg", height=10, width=15, units="cm", dpi=300)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_frequentOTUs_heatmap.pdf", height=10, width=15, units="cm", dpi=300)
## will need to manually edit the .pdf to make the Order labels fit better. stagger, then add silhouettes!

## might also try focusing on these OTUs at genus level...
## could collapse OTUs into shared genus-level (works for all but 2 OTUs)
## would then look across ENTIRE dataset and aggregate detections for those genus labels
## manually add back in the two other OTUs (1 and 377) also, though their numbers won't change
## can then maybe focus on whether there is any difference in our overall trend when examining at OTU vs. genus-level

## gathers the 14 defined Genera shared among those top OTUs; it's missing the two OTUs without any known genus
multiyear_frequentOTUs_Genus_list <- multiyear_frequentOTUs_bySiteYear_long_wTaxa %>% 
  distinct(Genus, OTUalias) %>% 
  filter(!is.na(Genus)) %>% 
  distinct(Genus)

## add in those two Genus-level data later:
multiyear_frequentOTUs_OTUnoGenus_list <- multiyear_frequentOTUs_bySiteYear_long_wTaxa %>% 
  distinct(Genus, OTUalias) %>% 
  filter(is.na(Genus)) %>% 
  select(OTUalias)

## let's collect the number of samples with those 14 genera identified
## this 'might' combine multiple OTUs with shared genus labels

## gather those data per genus label:
multiyear_frequentTaxa_bySiteYear_matrix_wTaxa_tmp = multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  mutate(GrouperTerm = paste0(Site, "-", Year)) %>% 
  filter(Genus %in% multiyear_frequentOTUs_Genus_list$Genus) %>%
  group_by(GrouperTerm, Order, Genus) %>% 
  summarise(nSamplesWithTaxa = n_distinct(SampleID)) %>% 
  merge(., multiyear_sampsPerGroup, by="GrouperTerm") %>% 
  mutate(pSamplesWithTaxa = nSamplesWithTaxa/nSamplesInGroup) %>% 
  ungroup() %>% 
  group_by(GrouperTerm) %>% 
  select(pSamplesWithTaxa, GrouperTerm, Order, Genus) %>% 
  mutate(pSamplesWithTaxa = round(pSamplesWithTaxa, 3)) %>% 
  pivot_wider(names_from="GrouperTerm", values_from="pSamplesWithTaxa")

multiyear_frequentTaxa_bySiteYear_matrix_wTaxa_tmp2 <- multiyear_otu_table_long_norm_wMeta_andTaxa %>% 
  mutate(GrouperTerm = paste0(Site, "-", Year)) %>% 
  filter(OTUalias %in% multiyear_frequentOTUs_OTUnoGenus_list$OTUalias) %>%
  group_by(GrouperTerm, Order, Family, OTUalias) %>% 
  summarise(nSamplesWithTaxa = n_distinct(SampleID)) %>% 
  merge(., multiyear_sampsPerGroup, by="GrouperTerm") %>% 
  mutate(pSamplesWithTaxa = nSamplesWithTaxa/nSamplesInGroup) %>% 
  ungroup() %>% 
  group_by(GrouperTerm) %>% 
  select(pSamplesWithTaxa, GrouperTerm, Order, Family, OTUalias) %>% 
  mutate(pSamplesWithTaxa = round(pSamplesWithTaxa, 3)) %>% 
  pivot_wider(names_from="GrouperTerm", values_from="pSamplesWithTaxa") %>% 
  mutate(Genus = paste0("f. ", Family, " ", OTUalias)) %>% 
  select(-Family, -OTUalias) %>% 
  relocate(Order, Genus, everything())

multiyear_frequentTaxa_bySiteYear_matrix_wTaxa <- rbind(multiyear_frequentTaxa_bySiteYear_matrix_wTaxa_tmp, multiyear_frequentTaxa_bySiteYear_matrix_wTaxa_tmp2)
# write_csv(multiyear_frequentTaxa_bySiteYear_matrix_wTaxa,
#           path = "~/Documents/nau_projects/guano/NHguano_redux/new_data/multiyear_frequentTaxa_bySiteYear_matrix_wTaxa.csv")

## same plot, just grouped by genus labels instead of OTUalias:
### plot this as heatmap showing the fraction of these OTUs detected 
multiyear_frequentTaxa_bySiteYear_long_wTaxa <- multiyear_frequentTaxa_bySiteYear_matrix_wTaxa %>% 
  pivot_longer(cols = c("COR-2015", "COR-2016", "HOP-2015", "HOP-2016", "MAP-2015", "MAP-2016"), 
               values_to = "pSamplesWithTaxa", names_to = "GrouperTerm") %>%
  filter(pSamplesWithTaxa > 0) %>% 
  mutate(Site = substr(GrouperTerm, 1, 3))

ggplot(multiyear_frequentTaxa_bySiteYear_long_wTaxa, 
       aes(x=Genus, y=GrouperTerm, fill=pSamplesWithTaxa)) +
  geom_tile(color="black") +
  scico::scale_fill_scico(begin = 0.2, breaks = c(0.25, 0.5, 0.75), limits = c(0,1)) +
  facet_grid(Site ~ Order, scales = "free", space="free") +
  labs(x = '', y='', fill='fraction\nsamples\nwith taxa\ndetected') +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        strip.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
        strip.text.y = element_blank())

# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_frequentTaxa_heatmap.png", height=10, width=15, units="cm", dpi=150)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_frequentTaxa_heatmap.svg", height=10, width=15, units="cm", dpi=300)
# ggsave("~/Documents/nau_projects/guano/NHguano_redux/new_figures/multiyear_frequentTaxa_heatmap.pdf", height=10, width=15, units="cm", dpi=300)
## will add silhouettes and resolve Order names for strip.x panels (in Illustrator) for publication

## cleanup
rm(list=ls(pattern = "^multiyear"))
