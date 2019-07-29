## ggmap for plotting diet info:
## see here for ggmap program: https://github.com/dkahle/ggmap
## see here for useful static and animated ggmap applications:
  ## http://www.bernhardlearns.com/2018/11/animating-r-plot-with-gganimate.html
## needed to register to access Google API to run ggmap; start here: https://console.cloud.google.com
  ## once API key created, in R runonce: register_google(key = "{enter_your_long_key}", account_type = "standard", write = TRUE)

library(ggmap)
library(tidyverse)
library(formattable)
library(htmltools)
library(webshot)
library(ggrepel)

## import metadata; import the full read dataset as well as the min-read/min-samp-per-Site/Window filtered dataset
alldat <- as.data.frame(read_csv(file="~/Repos/nhguano/data/filtered_dataset_wWindows.csv"))
filtdat <- as.data.frame(read_csv(file="~/Repos/nhguano/data/filtered_dataset_wWindows_min1000Reads_min2SamplesperSiteWindow.csv"))

################################################################################
## part 1 - summary table of collections at each site
################################################################################

## summarise number of samples collected at each site
collection_sumry_alldat <- as.data.frame(alldat %>% 
  group_by(Site, StudyID, SiteLat, SiteLong) %>% 
  summarise(SamplesCollected=n_distinct(SampleID), SamplingWindowsSurveyed=n_distinct(Window)) %>% 
  arrange(Site, StudyID)) %>% 
  mutate(Year=ifelse(StudyID=="oro15", "2015", "2016"))

collection_sumry_filtdat <- as.data.frame(filtdat %>% 
                                           group_by(Site, StudyID, SiteLat, SiteLong) %>% 
                                           summarise(SamplesAnalyzed=n_distinct(SampleID), SamplingWindowsAnalyzed=n_distinct(Window)) %>% 
                                           arrange(Site, StudyID)) %>% 
  mutate(Year=ifelse(StudyID=="oro15", "2015", "2016")) %>% 
  select(Site, StudyID, SamplesAnalyzed, SamplingWindowsAnalyzed)

## merge data together
collection_sumry <- merge(collection_sumry_alldat, collection_sumry_filtdat, by=c('Site', 'StudyID'), all.x = TRUE)
## fill missing weeks that were not longer in filtered sample (but were in the original unfiltered dataset)
collection_sumry[is.na(collection_sumry)] <- 0

## formattable object to save
csumtable <- formattable(collection_sumry %>% 
              select(Site, Year, SamplesCollected, SamplesAnalyzed, SamplingWindowsSurveyed, SamplingWindowsAnalyzed),
            list(SamplesCollected = color_tile("white", "darkorange"),
                 SamplesAnalyzed = color_tile("white", "lightgoldenrod2"),
                 SamplingWindowsSurveyed = color_tile("white", "dodgerblue2"),
                 SamplingWindowsAnalyzed = color_tile("white", "lightblue2")))

## export formattable object because table is too long for basic R export
## save table as 'alldata_collectionSummary'; 
  ## modify the "width = 30%" value if you want to stretch/shrink the spacing between columns
export_formattable <- function(f, file, width = "75%", height = NULL,
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

## export here; save as 'alldata_collectionSummary'
setwd("~/Repos/nhguano/figures/")
export_formattable(csumtable,"collectionSummary.png")

class(csumtable)

################################################################################
## part 2 - map of collection locations
################################################################################

## import metadata
metadata <- as.data.frame(read_csv(file="https://github.com/devonorourke/nhguano/raw/master/data/metadata/nhbat_meta.csv"))

## define boundaries?
latbot <- min(metadata$SiteLat) - .25
lattop <- max(metadata$SiteLat) + .25
latmid <- lattop - ((lattop - latbot) / 2)

lonleft <- min(metadata$SiteLong) - .25
lonright <- max(metadata$SiteLong) + .25
lonmid <- lonleft - ((lonleft - lonright) / 2)

nh <- c(left = lonleft, bottom = latbot, right = lonright, top = lattop)

## get directrly from Google
nhmap_bw <- get_googlemap(center = c(lonmid, latmid), zoom = 9, color = "bw")
nhmap_color <- get_googlemap(center = c(lonmid, latmid), zoom = 9, color = "color")

## see also: markers (add metadata markers with Site codes); path (layer data)

keepstudy1_sites <- c("FOX", "HOL", "MAP", "MTV", "CHI")
othermap_sites <- c("WLD", "SWZ", "WLT", "GRN")

## saving both colors; 
## save as 'collection_map_color'; export at 800x800
ggmap(nhmap_bw) +  ## switch to "nhmap_bw" to generate black/white plot
  geom_point(data = collection_sumry %>% filter(Year=="2016"), 
             aes(x=SiteLong, y=SiteLat, size=Samples, shape=Year, color=Year)) +
  geom_point(data = collection_sumry %>% filter(Year=="2015"), 
             aes(x=SiteLong, y=SiteLat, size=Samples, shape=Year, color=Year)) +
  scale_color_manual(values=c("#f1a340", "#998ec3")) +
  geom_label_repel(data = collection_sumry %>% filter(Site %in% keepstudy1_sites & Year == "2016"), 
#                   aes(x=SiteLong, y=SiteLat, label=Site), color="gray25", nudge_x = .15, nudge_y = -.05) +
                   aes(x=SiteLong, y=SiteLat, label=Site), color="gray25", nudge_x = .15, nudge_y = -.05, fill="#fbeeac") +
  geom_label_repel(data = collection_sumry %>% filter(!Site %in% keepstudy1_sites & Year == "2016" & SiteLong > -72), 
                   aes(x=SiteLong, y=SiteLat, label=Site), color="gray25", nudge_x = -.15, nudge_y = .05, label.size = 0.15) +
  geom_label_repel(data = collection_sumry %>% filter(!Site %in% keepstudy1_sites & Year == "2016" & SiteLong < -72),
                 aes(x=SiteLong, y=SiteLat, label=Site), color="gray25", nudge_x = .1, nudge_y = .0, label.size = 0.15) +
  geom_label_repel(data = collection_sumry %>% filter(Site %in% othermap_sites),
                   aes(x=SiteLong, y=SiteLat, label=Site), color="gray25", nudge_x = .0, nudge_y = .1, label.size = 0.15) +
  theme_bw(base_size = 16) +
  labs(x="", y="") +
  guides(color = guide_legend(override.aes = list(size = 6)))
   