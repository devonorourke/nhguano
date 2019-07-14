library(tidyverse)
library(qiime2R)
library(gganimate)
library(lubridate)
library(ggrepel)
library(ggpubr)
library(scales)
library(formattable)
library(htmltools)
library(webshot)

## importing read counts from all samples sequenced in study 
## this was the last file output from the 'decontam' R script
## see script here:
## see documentation of that process here:

## import read data; this was derived from the decontam_efforts.R script
## see: https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R
read_dat = read_csv("~/Repos/nhguano/data/filtered_dataset.csv")

## going to lump collections into 10-day bins instead of defining collection date by WOY
## because we want to keep window start/end dates consistent, going to use what was prevoiusly..
  ##..determined in earlier work to define these boundaries

## for each collection year...
mindate_2015 = read_dat %>% filter(StudyID=="oro15") %>% summarise(date=min(newDate)) %>% pull(date)
mindate_2016 = read_dat %>% filter(StudyID=="oro16") %>% summarise(date=min(newDate)) %>% pull(date)
  ## started earlier in 2016 than in 2015

## what's the latest in the season?
maxdate_2015 = read_dat %>% filter(StudyID=="oro15") %>% summarise(date=max(newDate)) %>% pull(date)
maxdate_2016 = read_dat %>% filter(StudyID=="oro16") %>% summarise(date=max(newDate)) %>% pull(date)
  ## collected later into 2016 too

## going to define sampling windows as 10-day bins, starting on January 1 of each year
## generate Day Of Year column
read_dat$DOY <- yday(read_dat$newDate)
## get range defined with the number of bins (devined in the "by" term of the 'breaks' argument)
window=14         ## 'window' represents number of days to span a given observation)
read_dat$DayRange = cut(read_dat$DOY, breaks = c(seq(from=1, to=365, by=window)), include.lowest = TRUE)
read_dat$Window <- as.numeric(read_dat$DayRange)
read_dat$RangeStart <- as.character(read_dat$DayRange)
read_dat$RangeStart <- gsub("\\(", "", read_dat$RangeStart)
read_dat$RangeStart <- gsub("\\]", "", read_dat$RangeStart)
read_dat$RangeStart <- gsub("\\[", "", read_dat$RangeStart)
read_dat <- read_dat %>% separate(., col = RangeStart, into=c("RangeStart", "RangeEnd"), sep = ",")
read_dat$RangeStart <- as.numeric(read_dat$RangeStart)
read_dat$RangeEnd <- as.numeric(read_dat$RangeEnd)
## convert RangeStart and RangeEnd back to dates 
orival <- ymd("2016-01-01")
read_dat$WindowStart <- as.Date(read_dat$RangeStart, origin = orival)
read_dat$WindowEnd <- as.Date(read_dat$RangeEnd, origin = orival)

## write this file to disk
write.csv(read_dat, file="~/Repos/nhguano/data/filtered_dataset_wWindows.csv", row.names = FALSE)

## filter data so that we retain only those samples with > 1000 reads (these should be the same samples as in the rarefied dataset)
## what samples have > 1000 reads?
min1000_samps <- read_dat %>% 
  group_by(SampleID) %>% 
  summarise(Reads=sum(Reads)) %>% 
  filter(Reads >= 1000) %>% 
  pull(SampleID)

tmp <- read_dat %>% 
  filter(SampleID %in% min1000_samps) %>% 
  group_by(StudyID, Window, Site) %>% 
  summarise(Samples=n_distinct(SampleID)) %>% 
  mutate(Labeler=paste0(Site,StudyID))

## filter any instances where there is just 1 sample per site/Window
tmp_filt <- tmp %>% filter(Samples > 1)

## filter read data so that we retain just those samples with > 1000 reads in Site/Windows with > 1 sample
tmp_dat_minreads <- read_dat %>% 
  mutate(Labeler=paste0(Site,StudyID)) %>% 
  filter(Labeler %in% tmp$Labeler & SampleID %in% min1000_samps)

tmp_dat_minreads_minSamples <- read_dat %>% 
  mutate(Labeler=paste0(Site,StudyID)) %>% 
  filter(Labeler %in% tmp_filt$Labeler & SampleID %in% min1000_samps)

## write each of these files to disk:
write.csv(tmp_dat_minreads, file="~/Repos/nhguano/data/filtered_dataset_wWindows_min1000Reads.csv", row.names = FALSE)
write.csv(tmp_dat_minreads_minSamples, file="~/Repos/nhguano/data/filtered_dataset_wWindows_min1000Reads_min2SamplesperSiteWindow.csv", row.names = FALSE)

################################################################################
## estimate the prevalent taxa across all sites
################################################################################
## generate dataset summarizing number of distinct (named) Family/Genus/Species..
  ##..as well as number of unique ASVs in that Order across entire dataset
## very different comparisons if we include ALL data (read_dat) vs. the sample filtered data
Order_sumry_alldat <- read_dat %>% 
  group_by(Order) %>% 
  summarise(Samples=n_distinct(SampleID),Family=n_distinct(Family),Genus=n_distinct(Genus),Species=n_distinct(Species),ASV=n_distinct(ASVid)) %>% 
  arrange(-ASV)

Order_sumry_minReads <- tmp_dat_minreads %>% 
  group_by(Order) %>% 
  summarise(Samples=n_distinct(SampleID), Family=n_distinct(Family), Genus=n_distinct(Genus), Species=n_distinct(Species), ASV=n_distinct(ASVid)) %>% 
  arrange(-ASV)

Order_sumry_minReads_minSamples <- tmp_dat_minreads_minSamples %>% 
  group_by(Order) %>% 
  summarise(Samples=n_distinct(SampleID), Family=n_distinct(Family), Genus=n_distinct(Genus), Species=n_distinct(Species), ASV=n_distinct(ASVid)) %>% 
  arrange(-ASV)

## plot tables; 
## save as 'Order_summary_allSites_allDat'; export at 500x470
formattable(Order_sumry_alldat, list(
  ASV = color_tile("white", "#b1d1fc"),
  Species = color_tile("white", "#76cd26"),
  Genus = color_tile("white", "#fcb001"),
  Family = color_tile("white", "#caa0ff")
))

## save as 'Order_summary_min1000reads_anySiteWindow'; export at 500x470
formattable(Order_sumry_minReads, list(
  ASV = color_tile("white", "#b1d1fc"),
  Species = color_tile("white", "#76cd26"),
  Genus = color_tile("white", "#fcb001"),
  Family = color_tile("white", "#caa0ff")
))


## save as 'Order_summary_min1000reads_atleast2sampesPerSiteWindow'; export at 500x470
formattable(Order_sumry_minReads_minSamples, list(
  ASV = color_tile("white", "#b1d1fc"),
  Species = color_tile("white", "#76cd26"),
  Genus = color_tile("white", "#fcb001"),
  Family = color_tile("white", "#caa0ff")
))

## even after dropping over 1000 samples, we don't see a big shift in the number of Species or ASVs detected
## the good news is that we're likely capturing the majority of the diversity in the dataset
## the bad news is that by requiring at least 1000 reads and 2 samples per Site/Window, we're losing a lot of samples

################################################################################
## highlight the particularly prevalent (most detected) taxa by Genus
################################################################################

## as before, show the difference between using the entire dataset versus a filtered one...
## all data:
genus_sumry_alldat <- as.data.frame(read_dat %>% 
                                      mutate(YearLabel=paste0(StudyID, Window)) %>% 
                                      group_by(Order, Family, Genus) %>% 
                                      summarise(ASVs=n_distinct(ASVid), Sites=n_distinct(Site), CollectionWindows=n_distinct(YearLabel),Samples=n_distinct(SampleID)) %>% 
                                      arrange(-ASVs))

genus_sumry_min1000reads <- as.data.frame(tmp_dat_minreads %>% 
                                      mutate(YearLabel=paste0(StudyID, Window)) %>% 
                                      group_by(Order, Family, Genus) %>% 
                                      summarise(ASVs=n_distinct(ASVid), Sites=n_distinct(Site), CollectionWindows=n_distinct(YearLabel),Samples=n_distinct(SampleID)) %>% 
                                      arrange(-ASVs))

genus_sumry_min1000reads_minSamples <- as.data.frame(tmp_dat_minreads_minSamples %>% 
                                            mutate(YearLabel=paste0(StudyID, Window)) %>% 
                                            group_by(Order, Family, Genus) %>% 
                                            summarise(ASVs=n_distinct(ASVid), Sites=n_distinct(Site), CollectionWindows=n_distinct(YearLabel),Samples=n_distinct(SampleID)) %>% 
                                            arrange(-ASVs))

## write full table to disk:
write.csv(genus_sumry_alldat, file="~/Repos/nhguano/data/diet_summary_tables/genus_summary_alldata.csv", row.names = FALSE, quote = FALSE)
write.csv(genus_sumry_min1000reads, file="~/Repos/nhguano/data/diet_summary_tables/genus_summary_min1000reads.csv", row.names = FALSE, quote = FALSE)
write.csv(genus_sumry_min1000reads_minSamples, file="~/Repos/nhguano/data/diet_summary_tables/genus_summary_min1000reads_min2SamplesPerSiteWindow.csv", row.names = FALSE, quote = FALSE)


## filter just those instances where a particular Genera are observed in at least 15 sites (of 20 total)
genus_fmtbl_alldat <- formattable(genus_sumry_alldat %>% filter(Sites >= 15) %>% arrange(-Sites, -Samples), 
                           list(ASVs = color_tile("white", "dodgerblue"),Samples = color_tile("white", "red")))

genus_fmtbl_minReads <- formattable(genus_sumry_min1000reads %>% filter(Sites >= 15) %>% arrange(-Sites, -Samples), 
                                  list(ASVs = color_tile("white", "dodgerblue"),Samples = color_tile("white", "red")))

genus_fmtbl_minReadsminSamps <- formattable(genus_sumry_min1000reads_minSamples %>% filter(Sites >= 15) %>% arrange(-Sites, -Samples), 
                                  list(ASVs = color_tile("white", "dodgerblue"),Samples = color_tile("white", "red")))

## export formattable object because table is too long for basic R export
## modify the "width = 30%" value if you want to stretch/shrink the spacing between columns
export_formattable <- function(f, file, width = "100%", height = NULL,
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

## export here:
setwd("~/Repos/nhguano/figures/")
export_formattable(genus_fmtbl_alldat,"Genus_summary_alldata.png")
export_formattable(genus_fmtbl_minReads,"Genus_summary_min1000Reads.png")
export_formattable(genus_fmtbl_minReadsminSamps,"Genus_summary_min1000readsmin2SamplesperSiteWindow.png")

################################################################################
## there appears to be a lot of shared taxa being detected at most sites...
## generate a few summary values
## focusing on just the filtered data here
################################################################################

## how many Genera are detected in at least 10 Sites (50%?)
genus_sumry_min1000reads_minSamples %>% 
  filter(!is.na(Genus) & Sites >= 10) %>% 
  group_by(Order) %>% tally() %>% arrange(-n)
  ## lots of different Lepidoptera (19), Coleoptera (17), Diptera (17), 
  ## but also several Ephemeroptera (5), Hemiptera (3), all others 2 or less

## how about 16 sites (80%?)
genus_sumry_min1000reads_minSamples %>% 
  filter(!is.na(Genus) & Sites >= 16) %>% 
  group_by(Order) %>% tally() %>% arrange(-n)
  ## fewer now... but still many Coleoptera (7) and some Diptera(4)...
  ## but Lepidoptera(2), Megaloptera (1) and Trichoptera (1) are only others

## what ASVs are the most prevalent? generate data:
ASVsumry <- tmp_dat_minreads_minSamples %>% 
  group_by(ASValias, Order, Genus, Species) %>% 
  summarise(Reads=sum(Reads), Samples=n_distinct(SampleID))

## 12-order palette used in other scripts
pal12 <- c('#fffe71', '#CE9834', '#9A6600', '#9BCC94', '#405E00', '#993303', 
           '#FF6501', '#336799', '#9ACEFF', 'gray75', '#D14A89', 'gray25')

## list of ASVs to highlight for discussion?
# notrun: prevASVs <- c("ASV-1", "ASV-2", "ASV-5", "ASV-7", "ASV-11", "ASV-12", "ASV-8", "ASV-4", "ASV-61", "ASV-20", "ASV-18", "ASV-27")

## plot scatterplot; save as "ASVsumry_alldata_Reads-and-Samples'; export at 900x450
ggplot(ASVsumry, aes(x=Samples, y=Reads, color=Order)) + 
  #geom_label_repel(data = ASVsumry %>% filter(ASValias %in% prevASVs), 
  #                 aes(x=Samples, y=Reads, color=Order, label=ASValias),
  #                 size=2.5, nudge_x = 10) +
  geom_point() +
  scale_y_continuous(trans = "log2", labels=comma_format(accuracy = 1)) +
  #scale_x_continuous(trans="log2", breaks = c(1,4,16,64,256)) +
  #annotation_logticks(sides = "b") +
  scale_color_manual(values=pal12) +
  theme_bw(base_size = 16) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(x="Samples", y="Sequence counts")


