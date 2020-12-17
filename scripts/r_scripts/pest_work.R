library(tidyverse)

## import file of counts of sequences per OTU per sample, filtered requiring a minimum of 1000 seqs per sample
read_dat <- read_csv("https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz")

################################################################################
## USDA comparisons
################################################################################

## import pest data frames 
## file initially dowonloaded from: https://www.aphis.usda.gov/aphis/ourfocus/planthealth/import-information/rppl/rppl-table
## converted into .csv file
usda_pest <- read_csv(file = "https://github.com/devonorourke/nhguano/raw/master/data/pests/usda_aphis_pestlist.csv") %>% 
  filter(Order %in% unique(read_dat$Order)) %>% 
  select(Order, Family, Genus, Species) %>% 
  distinct()

## how many are exact species matches between USDA pest list and our bat samples?
usda_match_species <- intersect(unique(read_dat$Species), usda_pest$Species) %>% na.omit()
usda_match_species[1]
## just one match: Lymantria dispar (gypsy moth)

## how many are matching common genera?
usda_match_genus <- intersect(unique(read_dat$Genus), usda_pest$Genus)
length(usda_match_genus)
  ## 94 match at genus level... but they're not lining up as exact species-level matches....

## how many samples and sites for each species here?
read_dat %>% 
  filter(Species %in% usda_match_species) %>% 
  group_by(Species) %>% 
  summarise(Samples=n_distinct(SampleID), Site=n_distinct(Site))
  ## just 4 samples at 4 sites with this gypsy moth...
  ## likely aren't going to find perfect matches; lots of USDA matches at Genus work, but are listed in their database as "Genus sp."
  ## ex. Maladera is "Maladera sp.", but we have Maladera castanea .... arg.

## how many of our OTUs match these 94 genera from above?
usda_mg <- read_dat %>% 
  filter(Genus %in% usda_match_genus) %>% 
  group_by(Order, Family, Genus) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(OTUreads))
  ## Melanotus is really preavlent (in 204 samples) - click beetles are ag. pests...
  ## Phyllophaga (469 samples), Maladera (330 samples), and several other Genera listed with 10's of samples

################################################################################
## Perdue comparisons
################################################################################

## import Perdue genera; unique Genera adapted from their "Priority Pest Lists" Excel file
## see: http://caps.ceris.purdue.edu/pest-surveillance-guidelines/priority-pest-list-excel/2019
purdue_genus <- read_csv(file="https://github.com/devonorourke/nhguano/raw/master/data/pests/perdue.pestlist_genera.csv", col_names = FALSE) %>% 
  rename(Genus = X1)

## what Genera in our dataset are in this list?
perdue_genMatch <- intersect(unique(read_dat$Genus), purdue_genus$Genus)
perdue_genMatch
  ## four genera listed:  
    ## - "Monochamus" (beetle)
    ## - "Dendroctonus" (beetle) 
    ## - "Crocidosema" (moth)
    ## - "Lymantria" (moth)

## import species list to see what Species these could be:
perdue_species <- read_delim(file = "https://github.com/devonorourke/nhguano/raw/master/data/pests/perdue.pestlist_species.csv", 
                             delim = " ", col_names = FALSE) %>% 
  rename(Genus = X1, Species = X2)
perdue_species$Species <- gsub("_", " ", perdue_species$Species)
perdue_match_species <- perdue_species %>% filter(Genus %in% perdue_genMatch)
perdue_match_species
  ## 11 possible species in these 4 genera

## are any of these species in our dataset?
read_dat %>% 
  filter(Species %in% perdue_species$Species) %>% 
  group_by(Species) %>% 
  summarise(Samples=n_distinct(SampleID), Site=n_distinct(Site))
  ## just one exactly: Lymantria dispar (gypsy moth)
    ## however, we have a bunch of 'Dendroctonus' with no Species name - not sure what it is... 'Dendroctonus sp.' isn't helpful...
    ## likewise, we have a few 'Monochamus' species named, but not the species names in the Perdue list!


################################################################################
## USFS comparisons
## these insects aren't invasive, but are known forest pests
################################################################################

## file initially created from USFS site: 
## https://www.fs.fed.us/foresthealth/publications.shtml 
## --> select "Forest Insect & Disease Leaflets"
## --> select "FIDLs by scientific name" for list

## converted into .tsv file; import
usfs_pest <- read_delim(file = "https://github.com/devonorourke/nhguano/raw/master/data/pests/usfs.pests.tsv", delim="\t",
                        col_names = FALSE) %>% 
  rename(FIDL=X1, Species=X2, CommonName=X3, Year=X4) %>% 
  mutate(Genus=Species) %>% 
  separate(., col=Genus, into=c("Genus", "delete"), sep = " ") %>% 
  select(-delete)

## any Species in our dataset?
usfs_matches <- read_dat %>% 
  filter(Species %in% usfs_pest$Species) %>% 
  group_by(Species) %>% 
  summarise(Samples=n_distinct(SampleID), Site=n_distinct(Site)) %>% 
  arrange(-Samples)
  ## everything super rare... just six distinct species named

## try searching just for Genus; noticed that some species names had typos (dissitri instead of dissitria)
usfs_genus_match <- read_dat %>% 
  filter(Genus %in% usfs_pest$Genus) %>% 
  group_by(Order, Genus, Species) %>% 
  summarise(Samples=n_distinct(SampleID), Sites=n_distinct(Site), nReads = sum(OTUreads)) %>% 
  arrange(-Samples)
usfs_genus_match
  ## Phyllophaga by far most prevalent; 
  ## forest tent caterpillar (Malacosoma disstria) in 25 samples at 13 sites... so certainly real, but not frequent
  ## northeastern pine sawyer (Monochamus notatus) same story: 36 samples; 12 sites...

################################################################################
### MERGE IT ALL
## group all the three datasets together: USDA, USFWS, Perdue, and see what we can gather:
################################################################################

usfs_pest_tmp <- usfs_pest %>% select(Genus, Species)
all_pest_tmp <- bind_rows(usfs_pest_tmp, perdue_species) %>%
  distinct()
rm(usfs_pest_tmp)

## remove duplicates
notDupPestNames <- all_pest_tmp %>% filter(!Species %in% usda_pest$Species) %>% select(Species) %>% pull()
all_pest_tmp2 <- all_pest_tmp %>% 
  filter(Species %in% notDupPestNames)

all_pest_df <- bind_rows(all_pest_tmp2, usda_pest) %>% 
  relocate(Order, Family, Genus, Species) %>% 
  filter(!is.na(Genus)) %>%  ## drop any instances of a pest listed with Family taxonomic label but lacking Genus
  mutate(Species = ifelse(Species == 'Malacosoma disstri', 'Malacosoma disstria', Species)) ## fix label mistake to match our taxa label

## now perform search at genus level, finding which of our taxa might fit:
GenusSelectorLabelPestDF <- all_pest_df %>% 
  mutate(GenusSelectorLabel = paste0(Order,"-",Genus)) %>% 
  select(GenusSelectorLabel) %>% pull()

tmp <- read_dat %>% 
  mutate(GenusSelectorLabel = paste0(Order,"-",Genus)) %>% 
  filter(GenusSelectorLabel %in% GenusSelectorLabelPestDF) %>% 
  group_by(Order, Genus) %>% 
  summarise(Samples=n_distinct(SampleID), Sites=n_distinct(Site), nReads = sum(OTUreads))

# ## which of these were detected in at least 5% of samples and at least 2 or more sites?
# GenusSelectorLabel_freqHits <- tmp %>%
#   filter(Samples > 45 & Sites > 2) %>%
#   mutate(OrderGenusLabel = paste0(Order,"-",Genus)) %>%
#   select(OrderGenusLabel) %>% pull()

# what particular taxa do we have for those Genera?
## summarise at the per-species level (note this will increase counts for any ambiguous species by aggregating them!)
tmp2 <- read_dat %>%
  mutate(GenusSelectorLabel = paste0(Order,"-",Genus)) %>%
  filter(GenusSelectorLabel %in% paste0(tmp$Order,"-",tmp$Genus)) %>%
  filter(OTUalias != 'OTU-1480') %>%  ## dropping this one weird instance of a mislabeled Genus
  group_by(Order, Family, Genus, Species) %>%
  summarise(nSamples = n_distinct(SampleID),
            nSites = n_distinct(Site),
            nDistinctYears = n_distinct(Year),
            nDistinctOTUs = n_distinct(OTUalias)) %>% 
  filter(nSamples > 9) %>% ## retain only those present in at least 10 samples
  ungroup() %>% 
  arrange(-nSamples) %>% 
  mutate(Species = gsub(".+? ", "", Species),
         Species = ifelse(is.na(Species), "sp.", Species))

## summarise only at the shared genus label here...         
tmp3 <- read_dat %>%
  mutate(GenusSelectorLabel = paste0(Order,"-",Genus)) %>%
  filter(GenusSelectorLabel %in% paste0(tmp$Order,"-",tmp$Genus)) %>%
  filter(OTUalias != 'OTU-1480') %>%  ## dropping this one weird instance of a mislabeled Genus
  group_by(Order, Family, Genus) %>%
  summarise(nSamples = n_distinct(SampleID),
            nSites = n_distinct(Site)) %>% 
  filter(Genus %in% tmp2$Genus) ## retain only those Genera already filtered from table above in tmp2


## create a data frame that summarizes the species that are contained within a particular Order/Genus label
  ## add the number of samples ID'd for each of these species in parentheses
pest_detect_sumry_tmp <- tmp2 %>% 
  mutate(PrintVal = paste0(Species, " (", nSamples, ")")) %>%
  select(Order, Genus, Species, PrintVal) %>% 
  group_by(Order, Genus) %>%
  summarise(DetectedSpecies = toString(PrintVal)) %>%
  ungroup()

pest_detect_sumry <- merge(pest_detect_sumry_tmp, tmp3, by=c("Order", "Genus")) %>%
  relocate(Order, Family, Genus, nSamples, nSites, DetectedSpecies) %>% 
  arrange(Order, -nSamples)

write_delim(pest_detect_sumry, 
          path = "~/github/nhguano/tables/table2_pestDetectedSumry.txt",
          delim = ";")
  ## note we'll updated this table2 in Excel to make a touch fancier and to highlight the exact ones that are pests,
  ## ... and distinguish those outside of NH/northeast

