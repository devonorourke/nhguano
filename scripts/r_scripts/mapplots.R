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

## import metadata
metadata <- as.data.frame(read_csv(file="~/Repos/nhguano/data/metadata/nhbat_meta.csv"))

################################################################################
## part 1 - summary table of collections at each site
################################################################################

## summarise number of samples collected at each site
collection_sumry <- as.data.frame(metadata %>% 
  group_by(Site, StudyID, SiteLat, SiteLong) %>% 
  summarise(Samples=n_distinct(SampleID), Weeks=n_distinct(WOY)) %>% 
  arrange(Site, StudyID)) %>% 
  mutate(Year=ifelse(StudyID=="oro15", "2015", "2016"))

## formattable object to save
csumtable <- formattable(collection_sumry %>% 
              select(Site, Year, Weeks, Samples),
            list(Samples = color_tile("white", "orange"),
                 Weeks = color_tile("white", "dodgerblue2")))

## export formattable object because table is too long for basic R export
## save table as 'alldata_collectionSummary'; 
  ## modify the "width = 30%" value if you want to stretch/shrink the spacing between columns
export_formattable <- function(f, file, width = "30%", height = NULL,
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
export_formattable(csumtable,"alldata_collectionSummary.png")

################################################################################
## part 2 - map of collection locations
################################################################################


## define boundaries?
nh <- c(left = -73, bottom = 42.5, right = -70.5, top = 45.5)

## get directrly from Google
get_googlemap(center = c(-71.768521,43.395290), 
              zoom = 9, color = "color") %>% ggmap()

## see also: markers (add metadata markers with Site codes); path (layer data)
