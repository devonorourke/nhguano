library(ggmap)
library(tidyverse)
library(urbnmapr)
library(cowplot)
library(lubridate)

## import site late/long vals only for samples used in diet analysis (not all samples collected)
site_data <- as.data.frame(read_csv(file="https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz")) %>% 
  distinct(Site, SiteLat, SiteLong) %>% 
  rename(lat = 'SiteLat', lon = 'SiteLong')

## define boundaries for base map
## could use just those areas we sampled within NH, but get's kinda too zoomed in...
# maxlat = max(metadata$SiteLat) + 0.15
# minlat = min(metadata$SiteLat) - 0.15
# maxlon = max(metadata$SiteLong) + 0.15
# minlon = min(metadata$SiteLong) - 0.15
# nhlatlonbox <- c(left = minlon, bottom = minlat, right = maxlon, top = maxlat)

## use a broader boundary to show most of NH instead of zooming into just southern half (makes inset map easier to understand)
alt_NHwide_boundary <- c(left = -72.6, bottom = 42.5, right = -70.5, top = 44.5)

## generate base map for plot
p_baseNHmap <- get_stamenmap(bbox = alt_NHwide_boundary,
              zoom = 10, 
              #maptype = "toner-background", ## shows roads and lakes ok, but would be better without roads
              maptype = "terrain-background",
              color = "color") %>% 
  ggmap(darken = c(0.2, "white")) ## increase number for more transparent color


## add in metadata to base map:
closeSites = c("CNA", "CNB", "WLT", "MAP", "BRN")
pNHwithlabs <- p_baseNHmap + 
  geom_label(data = site_data %>% filter(!Site %in% closeSites), aes(x=lon, y=lat, label=Site), size=3.5) +
  geom_label_repel(data = site_data %>% filter(Site %in% closeSites), aes(x=lon, y=lat, label=Site), size=3.5, direction = "x", segment.colour = "gray30")

pNHwithlabs
ggsave(filename = "~/github/nhguano/figures/figure1_nhmapNoSampleSizes.png", width=8, height=11, dpi=150)
ggsave(filename = "~/github/nhguano/figures/figure1_nhmapNoSampleSizes.svg", width=8, height=11, dpi=300)

#### get map of northeast USA to pin NH into this map at top right:
urbnmapr_data = urbnmapr::states
stateNamesToKeep <- c('Maine', 'New Hampshire', 'Vermont', 'New York', 'Massachusetts', 'Connecticut')
new_data <- urbnmapr_data %>% filter(state_name %in% stateNamesToKeep)

p_northeast <- ggplot() + 
  geom_polygon(data = new_data, mapping = aes(x = long, y = lat, group = group),
               fill = 'gray50', color = 'white') +
  geom_polygon(data = new_data %>% filter(state_name == "New Hampshire"), 
               mapping = aes(x = long, y = lat, group = group), fill = 'black', color = 'white') +
  #coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  theme_nothing() +
  theme(plot.background = element_rect(fill = "gray80", color = "black"))

#### now stitch those two together into a single plot:
ggdraw() +
  draw_plot(pNHwithlabs) +
  draw_plot(plot = p_northeast,
            x = 0.5,
            y = 0.65, 
            width = 0.42,
            height = 0.42,
            scale = 0.5)

ggsave(filename = "~/github/nhguano/figures/figure1_nhmapNoSampleSizes_withInset.png", height = 8, width = 11, dpi=150)
ggsave(filename = "~/github/nhguano/figures/figure1_nhmapNoSampleSizes_withInset.pdf", height = 8, width = 11, dpi=300)

## this image will be used in Figure 1, as panel "A", in addition to the landscape cover barplot generated in `spatial_work_NHguano.R`

#############
## let's also include a summary of the data collected vs. processed:
## using raw collection records to determine number of guano samples collected (not necessarily sequenced!)
## and certainly not necessarily sequenced AND generated enough data to be included in analaysis

## import 2015 and 2016 data
nhsites2016 <- c("ALS", "ROL", "COR", "MAP", "BRN", "MAS", "GIL", "MTV", "FOX", "CNB", "HOL", "PEN", 'CHI', 'EPS', 'CNA', 'HOP')

rawmeta2015 <- read_csv(file="https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/NHraw2015collectiondata.csv") %>% 
  rename(Site = "Location name", CollectionDate = "Collection Date", SampleID = `Sample ID`) %>% 
  dplyr::select(-Box, -`WEEK OF YEAR`, -status, -PlateNumber, -X8) %>% 
  mutate(Site = case_when(Site == "maple hill" ~ "MAP",
                          Site == "Hopkinton" ~ "HOP",
                          Site == "massabesic" ~ "MAS",
                          Site == "willard" ~ "WLD", 
                          Site == "brown lane" ~ "BRN",
                          Site == "fox state" ~ "FOX",
                          Site == "wilton" ~ "WLT",
                          Site == "Cornish" ~ "COR",
                          Site == "Greenfield" ~ "GRN",
                          Site == "Swanzey" ~ "SWZ",
                          Site == "gilsum" ~ "GIL",
                          TRUE ~ as.character(Site))) %>% 
  filter(Site != "unknown") %>% filter(!is.na(Site)) %>% ## drop any samples with unknown/missing Site info
  mutate(Date = mdy(CollectionDate),
         Year = year(Date),
         SampleID = paste0("oro15", SampleID)) %>% 
  dplyr::select(-CollectionDate)

rawmeta2016 <- read_csv(file="https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/NHraw2016collectiondata.csv") %>% 
  dplyr::select(StudyID, SampleID, LocationName, CollectionDate) %>% 
  mutate(Date = mdy(CollectionDate),
         Year = year(Date)) %>% 
  dplyr::select(-CollectionDate) %>% 
  filter(LocationName %in% nhsites2016) %>% 
  filter(StudyID == "oro16") %>% 
  filter(!grepl("contaminated", SampleID)) %>%  ## discard the few contaminated samples
  mutate(SampleID = paste0("oro16", SampleID),
         Year = ifelse(is.na(Year), 2016, Year)) %>% 
  rename(Site = "LocationName") %>% 
  dplyr::select(-StudyID)

rawmetaAll <- rbind(rawmeta2015, rawmeta2016)
rm(rawmeta2015, rawmeta2016)

## count how many samples were collected in each dataset:
allmeta_sampleSumry <- rawmetaAll %>% 
  group_by(Year, Site) %>% 
  summarise(SamplesCollected = n_distinct(SampleID)) %>% 
  pivot_wider(names_from = Year, values_from = SamplesCollected) %>% 
  rename(collected2015 = `2015`, collected2016 = `2016`)

## bring in seq data relevant to these samples... how many of these samples ended up generating sufficient seq data used in our diet analyses?
allrawSamplesWithSeqData <- read_csv("https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa_wDateBins.csv.gz")
allmeta_seqSumry <- allrawSamplesWithSeqData %>% 
  group_by(Year, Site) %>% 
  summarise(SamplesAnalyzed = n_distinct(SampleID)) %>% 
  pivot_wider(names_from = Year, values_from = SamplesAnalyzed) %>% 
  rename(evaluated2015 = `2015`, evaluated2016 = `2016`)

allmeta_fullSumry <- merge(allmeta_sampleSumry, allmeta_seqSumry, all = TRUE) %>% 
  relocate(Site, collected2015, evaluated2015, collected2016, evaluated2016) %>% 
  mutate_all(funs(replace_na(.,"..")))

write_csv(allmeta_fullSumry,
          path="~/github/nhguano/supplementaryData/tableS2_samplesPerSiteSumry.csv")


rm(allrawSamplesWithSeqData)
rm(allmeta_sampleSumry, allmeta_seqSumry)


#### unused code
# ## summarise number of samples collected at each site in each year in wide format for exportable .csv
# collection_sumry_wide <- metadata %>% 
#   group_by(Site, Year) %>% 
#   summarise(Samples=n_distinct(SampleID)) %>% 
#   arrange(Site, Year) %>% 
#   ungroup() %>%
#   pivot_wider(values_from = "Samples", names_from="Year", values_fill=0) %>% 
#   relocate(Site, `2015`, `2016`)

# get_stamenmap(bbox = alt_NHwide_boundary,
#               zoom = 10, 
#               #maptype = "toner-background", ## shows roads and lakes ok, but would be better without roads
#               maptype = "toner-background",
#               color = "color") %>% ggmap()


## summarise number of samples collected at each site in each year in long format for plotting
# collection_sumry_long <- metadata %>% 
#   group_by(Site, Year, SiteLat, SiteLong) %>% 
#   summarise(Samples=n_distinct(SampleID)) %>% 
#   arrange(Site, Year) %>% 
#   rename(lon = "SiteLong", lat = "SiteLat") %>% 
#   ungroup()

# tmp2 <- metadata %>% 
#   group_by(Site, Year, SiteLat, SiteLong) %>% 
#   summarise(Samples=n_distinct(SampleID)) %>% 
#   arrange(Site, Year) %>% 
#   rename(lon = "SiteLong", lat = "SiteLat") %>% 
#   ungroup() %>% 
#   mutate(Year = as.character(Year)) %>% 
#   mutate(lon = ifelse(Year == "2015", lon+0.05, lon-0.05))

## relabel metadata for better plot labels:
# tmp <- metadata %>% 
#   rename(lat = "SiteLat", lon = "SiteLong") %>% 
#   group_by(Site, Year, lat, lon) %>% 
#   summarise(Samples=n_distinct(SampleID)) %>% 
#   arrange(Site, Year) %>% 
#   pivot_wider(values_from = "Samples", names_from="Year") %>% 
#   #mutate(Site = as.factor(Site)) %>% 
#   mutate(LabelVal = case_when(is.na(`2016`) ~ paste0(Site, "\n", `2015`,"*"),
#                               is.na(`2015`) ~ paste0(Site, "\n", `2016`, "**"),
#                               TRUE ~ paste0(Site, "\n", `2015`, "*|", `2016`,"**"))) %>% 
#   select(-`2016`, -`2015`)
