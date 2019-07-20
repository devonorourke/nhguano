library(tidyverse)

## importing read counts from all samples sequenced in study 
## this was created in the decontam_efforts.R script and modified in the taxa_summaries script to include 10-day windows
## see script here: https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/taxa_summaries.R

## import read data (importing just the read-filtered, min-sample data)
read_dat = read_csv("https://github.com/devonorourke/nhguano/raw/master/data/filtered_dataset_wWindows_min1000Reads_min2SamplesperSiteWindow.csv")

## import pest data frames 
## file initially dowonloaded from:https://www.aphis.usda.gov/aphis/ourfocus/planthealth/import-information/rppl/rppl-table
## converted into .csv file
usda_pest <- read_csv(file = "~/Repos/nhguano/data/pests/usda_aphis_pestlist.csv") %>% 
  filter(Order %in% unique(read_dat$Order)) %>% 
  select(Order, Family, Genus, Species) %>% 
  distinct()
#usda_pest <- read_csv(file = "https://github.com/devonorourke/nhguano/raw/master/data/pests/usda_aphis_pestlist.csv")


################################################################################
## USDA comparisons
################################################################################
## how many are exact species matches between USDA pest list and our bat samples?
usda_match_species <- intersect(unique(read_dat$Species), usda_pest$Species) %>% na.omit()
  ## just three matches 

## how many are matching common genera?
usda_match_genus <- intersect(unique(read_dat$Genus), usda_pest$Genus)
## 101 match at genus level!?

## how many samples and sites for each species here?
read_dat %>% 
  filter(Species %in% usda_match_species) %>% 
  group_by(Species) %>% 
  summarise(Samples=n_distinct(SampleID), Site=n_distinct(Site))


## likely aren't going to find perfect matches; lots of USDA matches at Genus work, but are listed in their database as "Genus sp."
## ex. Maladera is "Maladera sp.", but we have Maladera castanea
## how many of our ASVs match these 101 genera?
usda_mg <- read_dat %>% 
  filter(Genus %in% usda_match_genus) %>% 
  group_by(Order, Family, Genus) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads))

## Melanotus hyslopi is really preavlent, as is unnamed species.
## USDA has "Melanotus sp." named ... 


## import Perdue genera; unique Genera adapted from their "Priority Pest Lists" Excel file
## see: http://caps.ceris.purdue.edu/pest-surveillance-guidelines/priority-pest-list-excel/2019
purdue_genus <- read_csv(file="~/Repos/nhguano/data/pests/perdue.pestlist_genera.csv", col_names = FALSE) %>% 
  rename(Genus = X1)

## what Genera in our dataset are in this list?
perdue_genMatch <- intersect(unique(read_dat$Genus), purdue_genus$Genus)
  ## four genera listed 

## import species list to see what Species these could be:
perdue_species <- read_delim(file = "~/Repos/nhguano/data/pests/perdue.pestlist_species.csv", 
                             delim = " ", col_names = FALSE) %>% 
  rename(Genus = X1, Species = X2)
perdue_species$Species <- gsub("_", " ", perdue_species$Species)
perdue_match_species <- perdue_species %>% filter(Genus %in% perdue_genMatch)

## are any of these species in our dataset?
read_dat %>% 
  filter(Species %in% perdue_species$Species) %>% 
  group_by(Species) %>% 
  summarise(Samples=n_distinct(SampleID), Site=n_distinct(Site))
  ## just one: Lymantria dispar


################################################################################
## USFS comparisons
## these insects aren't invasive, but are known forest pests
################################################################################

## file initially created from USFS site: https://www.fs.fed.us/foresthealth/
## converted into .tsv file; import
usfs_pest <- read_delim(file = "~/Repos/nhguano/data/pests/usfs.pests.tsv", delim="\t",
                        col_names = FALSE) %>% 
  rename(FIDL=X1, Species=X2, CommonName=X3, Year=X4) %>% 
  mutate(Genus=Species) %>% 
  separate(., col=Genus, into=c("Genus", "delete"), sep = " ") %>% 
  select(-delete)

## any Species in our dataset?
usfs_matches <- read_dat %>% 
  filter(Species %in% usfs_pest$Species) %>% 
  group_by(Species) %>% 
  summarise(Samples=n_distinct(SampleID), Site=n_distinct(Site))
## Monochamus (white-spotted sawyer) is the most prevalent 
## Enaphalodes (red oak borer)

## try searching just for Genus; noticed that some species names had typos (dissitri instead of dissitria)
usfs_genus_match <- read_dat %>% 
  filter(Genus %in% usfs_pest$Genus) %>% 
  group_by(Order, Genus, Species) %>% 
  summarise(Samples=n_distinct(SampleID), Sites=n_distinct(Site))
