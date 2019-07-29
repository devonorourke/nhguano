library(tidyverse)
library(reshape2)
library(qiime2R)
library(gganimate)
library(lubridate)
library(ggpubr)
library(scales)
library(geofacet)
library(vegan)

## importing read counts from all samples sequenced in study 
## this was created in the decontam_efforts.R script and modified in the taxa_summaries script to include 10-day windows
## see script here: https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/taxa_summaries.R

## using alpha diversity estimates calculated from QIIME 
## see https://github.com/devonorourke/nhguano/blob/master/docs/diversity_analyses.md

## import read data (importing just the read-filtered, min-sample data)
read_dat = read_csv("https://github.com/devonorourke/nhguano/raw/master/data/filtered_dataset_wWindows_min1000Reads_min2SamplesperSiteWindow.csv")

################################################################################
################################################################################
## Alpha diversity assessments
################################################################################
################################################################################

## import alpha diversity values from QIIME .qza artifacts
alphaimport <- function(urlpath, alphametric){
  tmp.qza <- read_qza(urlpath)
  tmp.dat <- tmp.qza$data %>% mutate(SampleID=row.names(.)) %>% mutate(Metric=alphametric)
  colnames(tmp.dat)[1] <- "Value"
  tmp.dat
}

## download files from url provided; running in this script from local path
## example not run: download.file(url = 'https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/alpha/alpha.vals_ob.qza', destfile = ob.qza)

ob_url <- "~/Repos/nhguano/data/qiime_qza/alpha/alpha.vals_ob.qza"
sh_url <- "~/Repos/nhguano/data/qiime_qza/alpha/alpha.vals_sh.qza"
fp_url <- "~/Repos/nhguano/data/qiime_qza/alpha/alpha.vals_fp.qza"

ob_dat <- alphaimport(ob_url, "observed")
sh_dat <- alphaimport(sh_url, "shannon")
fp_dat <- alphaimport(fp_url, "faithpd")

alpha_dat <- rbind(ob_dat, sh_dat, fp_dat)
rm(ob_url, fp_url, sh_url, ob_dat, sh_dat, fp_dat)

################################################################################
## going to split up dataset into 2016 and 2015 components and focus on selected sites
## Data is too discontinuous to leave all Sites for year/over/year comparisons
## Some sites too underrepresented to include for meaningful statistical analyses
## Subsetting 2016 and specific sites with majority of sampling first
################################################################################

## import metadata again
metadata <- as.data.frame(read_csv(file="https://github.com/devonorourke/nhguano/raw/master/data/metadata/nhbat_meta.csv"))
tmp_alphadat <- merge(alpha_dat, metadata) %>% select(-SampleType, -Date) ## selects just those samples present in alpha calculations (i.e. the samples that were rarefied)
tmp_alphadat <- tmp_alphadat %>% mutate(Year=ifelse(StudyID=="oro15", "2015", "2016"))
metadata <- metadata %>% filter(SampleID %in% tmp_alphadat$SampleID)

## Collective summary of all samples (to determine which sites to keep for 2016, etc.)
alpha_sumry <- tmp_alphadat %>% 
  filter(Metric=="observed") %>% 
  group_by(Site, Year) %>% 
  summarise(Samples=n_distinct(SampleID))

## Collecitve summary across _all_ sites (not just subset)
alpha_stats <- tmp_alphadat %>% 
  group_by(Site, Year, Metric) %>% 
  summarise(meanValue=mean(Value), medianValue=median(Value), sdValue=sd(Value)) %>% 
  mutate(meanValue=round(meanValue, 2), sdValue=round(sdValue, 2)) %>% 
  arrange(Metric, Site, 
          factor(Year, levels = c("2015", "2016")))

write.csv(alpha_stats, file="~/Repos/nhguano//data/stats/alpha/alpha_stats.csv",
          row.names = FALSE, quote = FALSE)
rm(alpha_stats)

## focusing on 2016 samples at these sites:
study1list <- c("MAP", "FOX", "BRN", "HOP", "HOL", "CNB", "CNA", "PEN", "MTV")
## update metadata to include just these Sites, 
study1meta <- metadata %>% filter(Site %in% study1list & StudyID == "oro16")
## now and add in Window information
tmp_metamrg <- read_dat %>% 
  filter(SampleID %in% study1meta$SampleID) %>% 
  select(SampleID, Window, WindowStart) %>% 
  distinct()
study1meta <- merge(study1meta, tmp_metamrg)
rm(tmp_metamrg, tmp_alphadat, alpha_sumry)

################################################################################
## Summary comparing WOY vs. Window strategy
################################################################################
study1_WOY_sumry <- read_dat %>% 
  filter(SampleID %in% study1meta$SampleID) %>% 
  group_by(Site, WOY) %>% 
  summarise(Samples = n_distinct(SampleID))

study1_Window_sumry <- read_dat %>% 
  filter(SampleID %in% study1meta$SampleID) %>% 
  group_by(Site, Window) %>% 
  summarise(Samples = n_distinct(SampleID))

## here's what WOY looks like
p1_woy <- ggplot(data = study1_WOY_sumry, aes(x=WOY, y=Site, label=Samples)) + geom_text() + theme_bw()
## here's what Window-approach (10 day bins) looks like
p1_win <- ggplot(data = study1_Window_sumry, aes(x=Window, y=Site, label=Samples)) + geom_text() + theme_bw()

## plot together; save as 'select2016Samples_windowAndWOY_comparison'; export at 1600x600
ggarrange(p1_woy, p1_win, ncol=2, labels=c("A", "B"))
  ## the Window approach has fewer gaps in data, and generally has a few more samples per observation
  ## given the low number of samples per site we'll want this window approach (lose the fine-scale temporal resolution, but at least we have the samples)

rm(metadata, p1_win, p1_win_f, p1_woy, p1_woy_f, 
   study1_win_sumry_filt, study1_Window_sumry, study1_WOY_sumry, study1_WOY_sumry_filt,
   filt10_samples, filt2_samples, filt5_samples, SiteWindowKeep, alphaimport)

################################################################################
## Are there distinct differences in alpha diversity per site per Window of time?
## Run ANOVA to test for main effects of Site and Time (Window)
################################################################################

## merge updated metadata (including Window) with alpha data
alpha_study1dat <- merge(study1meta, alpha_dat)

## set our categorical data for ANOVA as factors
alpha_study1dat$aovWindow <- as.factor(alpha_study1dat$Window)
alpha_study1dat$aovSite <- as.factor(alpha_study1dat$Site)

## collect just these names for future use in data filtering (in other scripts):
alpha_study1dat_names <- alpha_study1dat %>% select(SampleID) %>% unique()
colnames(alpha_study1dat_names) <- "#OTU ID"
write.table(alpha_study1dat_names, file="~/Repos/nhguano/data/metadata/alpha_study1names.txt", row.names = FALSE, quote = FALSE)

alpha_aov_ob <- aov(Value ~ aovSite * aovWindow, data=alpha_study1dat %>% filter(Metric == "observed"))
capture.output(summary(alpha_aov_ob),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_ob.txt")  ## Site sig, Window sig, interaction barely
alpha_aov_sh <- aov(Value ~ aovSite * aovWindow, data=alpha_study1dat %>% filter(Metric == "shannon"))
capture.output(summary(alpha_aov_sh),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_sh.txt")  ## Site sig, Window sig, interaction barely
alpha_aov_fp <- aov(Value ~ aovSite * aovWindow, data=alpha_study1dat %>% filter(Metric == "faithpd"))
capture.output(summary(alpha_aov_fp),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_fp.txt")  ## Site sig, Window barely, interaction barely
  ## seems like location and time of sampling both relevant, but interaction component suggests that it depends on which time/site we're looking at...

## Tukey posthoc summaries
tukeyfunction <- function(anovafile, whichfactor){
  posthoc_tmp_site <- TukeyHSD(x = anovafile, whichfactor)
  data.frame(posthoc_tmp_site[[1]]) %>% mutate(Pairs=row.names(.)) %>% arrange(p.adj) %>% mutate(p.adj=round(p.adj, 3)) %>% 
    mutate(diff=round(diff,2)) %>% mutate(lwr=round(lwr, 3)) %>% mutate(upr=round(upr,3))
}

tuk_ob_site <- tukeyfunction(alpha_aov_ob, 'aovSite')
tuk_ob_window <- tukeyfunction(alpha_aov_ob, 'aovWindow')
write.table(tuk_ob_site, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_ob_site.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 
write.table(tuk_ob_window, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_ob_window.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

tuk_sh_site <- tukeyfunction(alpha_aov_sh, 'aovSite')
tuk_sh_window <- tukeyfunction(alpha_aov_sh, 'aovWindow')
write.table(tuk_sh_site, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_sh_site.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 
write.table(tuk_sh_window, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_sh_window.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

tuk_fp_site <- tukeyfunction(alpha_aov_fp, 'aovSite')
tuk_fp_window <- tukeyfunction(alpha_aov_fp, 'aovWindow')
write.table(tuk_fp_site, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_fp_site.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 
write.table(tuk_fp_window, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_fp_window.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

rm(list=ls(pattern="tuk_*"))
rm(list=ls(pattern="alpha_aov*"))

################################################################################
## Let's plot how these are changing by Site and Time (for each Window)
## Generating a static plot is going to be tough to see all at once..
## ..so we'll also generate a matching animation for the static one
################################################################################

## static boxplot for each Site+Time window for selected 2016 sites
## order by geographic proximity (relatively anyway)
alpha_study1dat$Site <- factor(alpha_study1dat$Site, levels = c(
  "HOL", "FOX", "HOP", "PEN", "CNA", "CNB", "MTV", "MAP", "BRN"))

## plot boxplot; save as 'alpha_ob_boxplotPerWindowPerSite_select2016'; export at 500x1000
ggplot(data=alpha_study1dat %>% filter(Metric=="observed"), 
       aes(x=Window, y=Value, group=Window), fill="gray50") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  facet_grid(Site ~ .) +
  scale_x_continuous(breaks=c(7,10,13,16,19,22),
                     labels=c('2016-03-26', '2016-05-07', '2016-06-18', '2016-07-30', '2016-09-10', '2016-10-22')) +
  #scale_y_continuous(trans="log2", labels = comma_format(accuracy = 1)) + 
  labs(x="", y="Observed richness") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))

## plot boxplot; save as 'alpha_sh_boxplotPerWindowPerSite_select2016'; export at 500x1000
psh <- ggplot(data=alpha_study1dat %>% filter(Metric=="shannon"), 
       aes(x=Window, y=Value, group=Window), fill="gray50") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  facet_grid(Site ~ .) +
  scale_x_continuous(breaks=c(7,10,13,16,19,22),
                     labels=c('2016-03-26', '2016-05-07', '2016-06-18', '2016-07-30', '2016-09-10', '2016-10-22')) +
  labs(x="", y="Shannon's entropy") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))


## plot boxplot; save as 'alpha_fp_boxplotPerWindowPerSite_select2016'; export at 800x1200
pfp <- ggplot(data=alpha_study1dat %>% filter(Metric=="faithpd"), 
       aes(x=Window, y=Value, group=Window), fill="gray50") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  facet_grid(Site ~ .) +
  scale_x_continuous(breaks=c(7,10,13,16,19,22),
                     labels=c('2016-03-26', '2016-05-07', '2016-06-18', '2016-07-30', '2016-09-10', '2016-10-22')) +
  labs(x="", y="Faith's PD") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))

## save as 'alpha_shanfaith_combo'; export at 
ggarrange(psh, pfp, ncol=2, labels=c("A", "B"))

## animate the observed ASVs; note we've switched the x axis to Sites, and are transitioning through the Windows of date
## first one is for observed ASVs
ob_ani1 <- ggplot(data=alpha_study1dat %>% filter(Metric=="observed"), 
       aes(x=Site, y=Value, group=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_y_continuous(trans="log2", labels = comma_format(accuracy = 1)) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(x="", y="Observed richness", title = 'Date: {closest_state}') +
  ease_aes('quadratic-in-out') +
  transition_states(WindowStart, 1, 3)
ob_rendered <- animate(ob_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/alpha_ob.gif", animation = ob_rendered)
rm(ob_ani1, ob_rendered)

## next one for Shannon's entropy alpha vals
sh_ani1 <- ggplot(data=alpha_study1dat %>% filter(Metric=="shannon"), 
                  aes(x=Site, y=Value, group=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(x="", y="Shannon's entropy", title = 'Date: {closest_state}') +
  ease_aes('quadratic-in-out') +
  transition_states(WindowStart, 1, 3)
sh_rendered <- animate(sh_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/alpha_sh.gif", animation = sh_rendered)
rm(sh_ani1, sh_rendered)

## last one for Faiths's PD alpha vals
fp_ani1 <- ggplot(data=alpha_study1dat %>% filter(Metric=="faithpd"), 
                  aes(x=Site, y=Value, group=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(x="", y="Faith's PD", title = 'Date: {closest_state}') +
  ease_aes('quadratic-in-out') +
  transition_states(WindowStart, 1, 3)
fp_rendered <- animate(fp_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/alpha_fp.gif", animation = fp_rendered)
rm(fp_ani1, fp_rendered)


################################################################################
################################################################################
## Relative abundance stacked barplots
## first plot shows all data; second plot is just for selected 2016 sites
## Producing both static and animation equivalent for 2016 only
################################################################################
################################################################################
## use same 12-color palette:
pal12 <- c('#fffe71', '#CE9834', '#9A6600', '#9BCC94', '#405E00', '#993303', 
           '#FF6501', '#336799', '#9ACEFF', 'gray75', '#D14A89', 'gray25')

## all data:
barpoints <- read_dat %>% 
  group_by(StudyID, Site, Window) %>% 
  summarise(Samples=n_distinct(SampleID)) %>% 
  mutate(Year=ifelse(StudyID=="oro15", "2015", "2016")) %>% 
  mutate(Labeler=paste(Site,Year,sep="-")) %>% 
  filter(Samples > 1)

barplot_plotdat <- read_dat %>% 
  mutate(Year=ifelse(StudyID=="oro15", "2015", "2016"), Labeler=paste(Site,Year,sep="-")) %>% 
  filter(Labeler %in% barpoints$Labeler) %>% 
  group_by(Site, Year, Labeler, Window, WindowStart, Order) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(pReads = Reads/sum(Reads))

## rearrange the levels for the facet wrap:
barplot_plotdat$Labeler <- factor(barplot_plotdat$Labeler, levels = c(
  "BRN-2015","COR-2015","FOX-2015","GIL-2015","HOP-2015",
  "BRN-2016","COR-2016","FOX-2016","GIL-2016","HOP-2016",
  "MAP-2015","MAS-2015","SWZ-2015","WLD-2015","WLT-2015",
  "MAP-2016","MAS-2016","ALS-2016","CHI-2016","CNA-2016",
  "EPS-2016","HOL-2016","MTV-2016","PEN-2016","CNB-2016"))

## plot barplot; save as 'Order_relAbund_stackedBars_allDat'; export at 950x950
ggplot(data = barplot_plotdat,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  scale_x_continuous(breaks=c(7,10,13,16,19,22),
                     labels=c('Mar-26', 'May-07', 'Jun-18', 'Jul-30', 'Sep-10', 'Oct-22')) +
  scale_y_continuous(breaks = c(0,.5, 1), labels=c(0, 50, 100)) +
  facet_wrap(~ Labeler, nrow=5) +
  labs(x="", y="Relative Abundance (%)", fill="Order") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=12), legend.position = "top") +
  guides(fill = guide_legend(nrow = 2))

########## can also plot the frequency of occurrence, not abundance of reads; just one minor switch here:
barplot_plotdat2 <- read_dat %>% 
  mutate(Year=ifelse(StudyID=="oro15", "2015", "2016"), Labeler=paste(Site,Year,sep="-")) %>% 
  filter(Labeler %in% barpoints$Labeler) %>% 
  group_by(Site, Year, Labeler, Window, WindowStart, Order) %>% 
  mutate(Reads=ifelse(Reads > 0, 1, 0)) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(pReads = Reads/sum(Reads))

## rearrange the levels for the facet wrap:
barplot_plotdat2$Labeler <- factor(barplot_plotdat2$Labeler, levels = c(
  "BRN-2015","COR-2015","FOX-2015","GIL-2015","HOP-2015",
  "BRN-2016","COR-2016","FOX-2016","GIL-2016","HOP-2016",
  "MAP-2015","MAS-2015","SWZ-2015","WLD-2015","WLT-2015",
  "MAP-2016","MAS-2016","ALS-2016","CHI-2016","CNA-2016",
  "EPS-2016","HOL-2016","MTV-2016","PEN-2016","CNB-2016"))

## plot barplot; save as 'Order_relOccur_stackedBars_allDat'; export at 950x800
ggplot(data = barplot_plotdat2,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  scale_x_continuous(breaks=c(7,10,13,16,19,22),
                     labels=c('Mar-26', 'May-07', 'Jun-18', 'Jul-30', 'Sep-10', 'Oct-22')) +
  scale_y_continuous(breaks = c(0,.5, 1), labels=c(0, 50, 100)) +
  facet_wrap(~ Labeler, nrow=5) +
  labs(x="", y="Proportion of detections (%)", fill="Order") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=12), legend.position = "top") +
  guides(fill = guide_legend(nrow = 2))


## for 2016 data now... calculate per-Order read abundance; grouping all samples with shared Site+Window
barplot_study1_plotdat <- read_dat %>% 
  filter(SampleID %in% alpha_study1dat_names$`#OTU ID`) %>% 
  group_by(Site, Window, WindowStart, Order) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(pReads = Reads/sum(Reads))

## 12 color palette for the unique Orders here
pal12 <- c('#fffe71', '#CE9834', '#9A6600', '#9BCC94', '#405E00', '#993303', 
           '#FF6501', '#336799', '#9ACEFF', 'gray75', '#D14A89', 'gray25')

## values of colors are (yellow, lightbrown, darkbrown, lightgreen, darkgreen, darkred, 
##                       orange, darkblue, lightblue, lightgray, purple, darkgray)

## rearrange the levels for the facet wrap here to match previous ordering by Site
barplot_study1_plotdat$Site <- factor(barplot_study1_plotdat$Site, levels = c(
  "HOL", "FOX", "HOP", "PEN", "CNA", "CNB", "MTV", "MAP", "BRN"))

## plot is faceting Sites, showing each Window as separate stacked bars of Order relative abundances
## plot barplot; save as 'Order_relAbund_stackedBars_select2016'; export at 800x800
ggplot(data = barplot_study1_plotdat,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  scale_x_continuous(breaks=c(7,10,13,16,19,22),
                     labels=c('2016-03-26', '2016-05-07', '2016-06-18', '2016-07-30', '2016-09-10', '2016-10-22')) +
  scale_y_continuous(breaks = c(0,.5, 1)) +
  facet_grid(Site ~ .) +
  labs(x="", y="Relative Abundance (%)", fill="Order") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


################################################################################
## Animation of above plot; transitioning frames with Windows
## changing x-axis to Site instead of time (like with other alpha div plots)
################################################################################
bar_ani1 <- ggplot(data = barplot_study1_plotdat, aes(x=Site, y=pReads, group=Site, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  labs(x="", y="Relative Abundance (%)", fill="Taxonomic\nOrder",
       title = 'Sampling Date: {closest_state}') +
  theme_bw(base_size = 16) +
  ease_aes('quadratic-in-out') +
  transition_states(WindowStart, 1, 3)
bar_rendered <- animate(bar_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/taxaOrder_stackedBars_select2016.gif", animation = bar_rendered)
rm(bar_ani1, bar_rendered)

################################################################################
## It's possible to set up these bar charts in relative positions as a grid
## We can create a custom plot using geofacet library.
## This wasn't run because there weren't as many sites that needed to be compared
## Code here is preserved to just use as an example in future
## See: https://hafen.github.io/grid-designer/
################################################################################

## this script was used to figure out how to use existing lat/lon for each site and convert to a grid
## it maintains a 1/10th of a degree as a cell, so spacing in resulting geofacet grid is still basically proportional
# notrun: read_dat %>% select(Site, SiteLat, SiteLong) %>% distinct() %>% mutate(tmprow=SiteLat*10, tmpcol=-(SiteLong)*10) %>% 
  #       mutate(row=abs(ceiling(tmprow-(max(tmprow)))) + 1) %>% mutate(col=(ceiling(max(tmpcol)-tmpcol))+1) %>% select(Site, row, col)

## I looked at the output of the above code, pasted in a spreadsheet to figure out which row/col's were duplicate, and adjusted
## the following row/cols were derived from that spreadsheet, but I modified two pairs of sites (BRN/MAP and CNA/CNB) that had overlapping grids

# mygrid <- data.frame(
#  name = c("Holderness", "Cornish", "CanterburyW", "CanterburyE", "Penacook", "Chichester", "Alsted", "Hillsborough","Hopkinton", 
#           "Epsom", "Rollins", "Gilsum", "Antrim", "Massabesic", "Swanzey", "Greenfield", "MontVernon","Wilton", "HollisN", "HollisS"),
#  code = c("HOL", "COR", "CNA", "CNB", "PEN", "CHI", "ALS", "FOX", "HOP", 
#           "EPS", "ROL", "GIL", "WLD", "MAS", "SWZ", "GRN", "MTV", "WLT", "BRN", "MAP"),
#  row = c(1, 3, 4, 4, 5, 5, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9, 9, 10, 10, 11),
#  col = c(9, 1, 10, 11, 9, 12, 2, 6, 8, 12, 17, 2, 5, 12, 3, 7, 9, 8, 9, 9), stringsAsFactors = FALSE)

## For the 2016 dataset, the plots were too spread out to see the barplots easily, so I modified the distances between sites
## these are NOT well attributed to their geographic coordinates, but their relative positions are maintained
my2016grid <- data.frame(
  name = c("Holderness", "CanterburyW", "CanterburyE", "Penacook", "Hillsborough","Hopkinton", "MontVernon", "HollisN", "HollisS"),
  code = c("HOL", "CNA", "CNB", "PEN", "FOX", "HOP", "MTV", "BRN", "MAP"),
  row = c(1,2,2,3,4,4,5,6,6), col = c(3,3,4,3,1,2,2,2,3))
  
## add the 'code' name to each 
barplot_study1_plotdat$code <- barplot_study1_plotdat$Site

## save static plot as 'geofacet_barplots_2016sites'; export at 800x800
ggplot(data = barplot_study1_plotdat,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  scale_y_continuous(breaks = c(0,.5, 1)) +
  facet_geo(~ code, grid = my2016grid) +
  labs(x="", y="Relative Abundance (%)", fill="Order") +
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "top") +
  guides(fill = guide_legend(nrow = 2))


rm(barplot_study1_plotdat, alpha_dat, alpha_study1dat, read_dat, alpha_study1dat_names)

################################################################################
################################################################################
## Beta diversity assessments for 2016 selected sites
################################################################################
################################################################################

## import distance values from QIIME for Adonis:
## import diversity files
adonisfunction <- function(urlpath, betametric){
  tmp.qza <- read_qza(urlpath)
  tmp.dist <- tmp.qza$data
  tmp.adonis <- adonis(formula = tmp.dist ~ Site * Window, data = study1meta)
  tmp.adonis_data <- tmp.adonis$aov.tab
  tmp.adonis_data %>% 
    mutate(Terms=row.names(.), Metric=betametric) %>% 
    select(Terms, Df, SumsOfSqs, MeanSqs, F.Model, R2, `Pr(>F)`, Metric)
}

## download files from url provided; running in this script from local path
## example not run: download.file(url = 'https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/distmat/ds_dist.qza', destfile = ds_dist.qza)

ds.dist_url <- "~/Repos/nhguano/data/qiime_qza/distmat/select2016/s16_dist_ds.qza"
bc.dist_url <- "~/Repos/nhguano/data/qiime_qza/distmat/select2016/s16_dist_bc.qza"
uu.dist_url <- "~/Repos/nhguano/data/qiime_qza/distmat/select2016/s16_dist_uu.qza"
wu.dist_url <- "~/Repos/nhguano/data/qiime_qza/distmat/select2016/s16_dist_wu.qza"
dc.dist_url <- "~/Repos/nhguano/data/qiime_qza/distmat/select2016/deicode_distmat.qza"
dc_min10.dist_url <- "~/Repos/nhguano/data/qiime_qza/distmat/select2016/deicode_distmat_min10feat.qza"

ds_adonis_df <- adonisfunction(ds.dist_url, "ds")
bc_adonis_df <- adonisfunction(bc.dist_url, "bc")
uu_adonis_df <- adonisfunction(uu.dist_url, "uu")
wu_adonis_df <- adonisfunction(wu.dist_url, "wu")
dc_adonis_df <- adonisfunction(dc.dist_url, "dc")
dcmin10_adonis_df <- adonisfunction(dc_min10.dist_url, "dc_min10")

## write files to disk:
write_csv(ds_adonis_df, path = "~/Repos/nhguano/data/stats/beta/ds_adonis_select2016.csv", quote=FALSE)
write_csv(bc_adonis_df, path = "~/Repos/nhguano/data/stats/beta/bc_adonis_select2016.csv", quote=FALSE)
write_csv(uu_adonis_df, path = "~/Repos/nhguano/data/stats/beta/uu_adonis_select2016.csv", quote=FALSE)
write_csv(wu_adonis_df, path = "~/Repos/nhguano/data/stats/beta/wu_adonis_select2016.csv", quote=FALSE)
write_csv(dc_adonis_df, path = "~/Repos/nhguano/data/stats/beta/dc_adonis_select2016.csv", quote=FALSE)
write_csv(dcmin10_adonis_df, path = "~/Repos/nhguano/data/stats/beta/dc_min10_adonis_select2016.csv", quote=FALSE)

rm(ds.dist_url, bc.dist_url, uu.dist_url, wu.dist_url, dc.dist_url, dc_min10.dist_url,
   ds_adonis_df, bc_adonis_df, uu_adonis_df, wu_adonis_df, dc_adonis_df, dcmin10_adonis_df)

## import PCoA files for plotting
ds.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_ds.qza"
bc.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_bc.qza"
uu.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_uu.qza"
wu.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_wu.qza"

pcoaplotfunction <- function(urlpath, betametric){
  tmp.qza <- read_qza(urlpath)
  tmp.data <- tmp.qza$data$Vectors[,1:3]
  tmp.lab.pc1 <- tmp.qza$data$ProportionExplained[1] %>% mutate(PC1=PC1*100) %>% pull() %>% round(., 2) %>% paste("PC1", ., sep = " ") %>% paste0(., "%")
  tmp.lab.pc2 <- tmp.qza$data$ProportionExplained[2] %>% mutate(PC2=PC2*100) %>% pull() %>% round(., 2) %>% paste("PC2", ., sep = " ") %>% paste0(., "%")
  data.frame(tmp.data, tmp.lab.pc1, tmp.lab.pc2, betametric)
}

ds_pcoa_df <- pcoaplotfunction(ds.pcoa_url, "ds")
bc_pcoa_df <- pcoaplotfunction(bc.pcoa_url, "bc")
uu_pcoa_df <- pcoaplotfunction(uu.pcoa_url, "uu")
wu_pcoa_df <- pcoaplotfunction(wu.pcoa_url, "wu")

## plot separately to include each unique x/y axis PC % variance label, then pull together in single plot
pcoaplotvizfunction <- function(data, CaptionName){
  tmp_df <- merge(data, study1meta)
  pc1lab <- tmp_df %>% distinct(tmp.lab.pc1)
  pc1lab <- as.character(pc1lab$tmp.lab.pc1)
  pc2lab <- tmp_df %>% distinct(tmp.lab.pc2)
  pc2lab <- as.character(pc2lab$tmp.lab.pc2)
  ggplot(data = tmp_df, aes(x=PC1, y=PC2, label=Site, color=Window)) +
    geom_text(size = 2.5) + 
    scale_x_continuous(limits = c(-.63,.63)) +
    scale_y_continuous(limits = c(-.63,.63)) +
    scale_color_distiller(palette = "Spectral",
                          name = "Sampling Window",
                          breaks = c(7, 10, 13, 16, 19, 22), 
                          labels = c("Apr-7", "May-20", "June-18", "July-30", "Sept-10", "Oct-22")) +
    labs(x=pc1lab, y=pc2lab, subtitle = CaptionName) +
    theme(legend.position = "top") +
    theme(panel.background = element_rect(fill = "gray40"),
          plot.background = element_rect(fill = "gray95"),
          legend.key.width = unit(3.5, "cm"))
}

## gather individual plots
ds_pcoa_plot <- pcoaplotvizfunction(ds_pcoa_df, "Dice-Sorensen")
bc_pcoa_plot <- pcoaplotvizfunction(bc_pcoa_df, "Bray-Curtis")
uu_pcoa_plot <- pcoaplotvizfunction(uu_pcoa_df, "unweighted Unifrac")
wu_pcoa_plot <- pcoaplotvizfunction(wu_pcoa_df, "weighted Unifrac")

## plot as one big set; save as 'pcoa_select2016_all4'; export at 925x925
ggarrange(ds_pcoa_plot, bc_pcoa_plot, uu_pcoa_plot, wu_pcoa_plot,
          ncol=2, nrow=2, common.legend = TRUE)

## plot each separately too; export each at 900x900
## save as 'pcoa_ds_select2016'
ds_pcoa_plot

## save as 'pcoa_bc_select2016'
bc_pcoa_plot

## save as 'pcoa_uu_select2016'
uu_pcoa_plot

## save as 'pcoa_wu_select2016'
wu_pcoa_plot

rm(list=ls(pattern="_pcoa_*"))
rm(list=ls(pattern="*_url"))
rm(list=ls(pattern="*lab*"))

################################################################################
## Biplots to add informative species?
## performing just for weighted Unifrac as an example
## adapted from: https://forum.qiime2.org/t/how-to-make-pcoa-biplot-in-r-using-q2-deicode-ordination/8377/5
################################################################################

## read in data
wu_biplot_url <- "/Users/do/Repos/nhguano/data/qiime_qza/biplots/s16_pcoabiplot_wu.qza"
wu_biplot.qza <- read_qza(wu_biplot_url)

## merge metadata with ordi data for plot
wu_fullplotdata <- wu_biplot.qza$data$Vectors %>% left_join(study1meta)

## get taxonomy information
tax <- read_csv(file = "https://github.com/devonorourke/nhguano/raw/master/data/filtered_dataset_wWindows_min1000Reads_min2SamplesperSiteWindow.csv") %>% 
  distinct(ASVid, Class, Order, Family, Genus, Species, ASValias) %>% 
  rename(FeatureID=ASVid)

## save as 'biplot_wu_select2016'; export at 925x925
bigplot <- ggplot() +
  geom_text(data=wu_fullplotdata, aes(x=PC1, y=PC2, label=Site, color=Window), size=2.5) +
  scale_color_distiller(palette = "Spectral", 
                        name = "Sampling Window",
                        breaks = c(7, 10, 13, 16, 19, 22), 
                        labels = c("Apr-7", "May-20", "June-18", "July-30", "Sept-10", "Oct-22")) +
  theme(panel.background = element_rect(fill = "gray40"),
        plot.background = element_rect(fill = "gray95"),
        legend.position = "top",
        legend.key.width = unit(3.5, "cm")) +
  xlab(paste(paste("PC1", round(100*wu_biplot.qza$data$ProportionExplained[1],2),"%"))) +
  ylab(paste(paste("PC2", round(100*wu_biplot.qza$data$ProportionExplained[2],2),"%"))) +
  geom_segment(data=wu_biplot.qza$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                 top_n(8, a) %>% #keep 8 furthest away points
                 mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% # scale arrows linearly... is this ok? 
                 left_join(tax),
               aes(x=0, xend=PC1, y=0, yend=PC2),
               arrow = arrow(length = unit(0.3,"cm")))

## save as 'biplot_inset_wu_select2016'; export at 660x660
insetplot <- ggplot() + 
  theme(panel.background = element_rect(fill = "gray90"), 
        plot.background = element_rect(fill = "gray90"),
        legend.position = "top") +
  geom_segment(data=wu_biplot.qza$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% 
                 top_n(8, a) %>% 
                 mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% 
                 left_join(tax),
               aes(x=0, xend=PC1, y=0, yend=PC2, color=Order),
               size=2,
               arrow = arrow(length = unit(0.3,"cm"))) +
  geom_label_repel(data=wu_biplot.qza$data$Species %>% 
               mutate(a=sqrt(PC1^2+PC2^2)) %>% 
               top_n(8, a) %>% 
               mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% 
               select(FeatureID, PC1, PC2) %>% 
               left_join(tax),
            aes(x=PC1, y=PC2, fill=Order, label=ASValias), force=5, 
            segment.alpha = 0.1, color="black", size=5) +
  scale_fill_manual(values = c('#9A6600', '#9BCC94', '#405E00')) + ## to match other barplots!
  scale_color_manual(values = c('#9A6600', '#9BCC94', '#405E00')) +
  labs(x="", y="")

insetplot

## can't get both to set into each other for some reason... tried ggpubr, cowplot, ... maybe because of x/y negative dimensions?
## save individual plots and stitch together in Illustator :(

################################################################################
## Biplot from DEICODE
## see program info here: https://github.com/biocore/DEICODE
## see tutorial here: https://forum.qiime2.org/t/robust-aitchison-pca-beta-diversity-with-deicode/8333
## relies on non-rarefied data; 
################################################################################

## read in data
dc_biplot_url <- "/Users/do/Repos/nhguano/data/qiime_qza/biplots/deicode_biplot.qza"
dc_biplot.qza <- read_qza(dc_biplot_url)

## create data frame for species plot
dc_fullplotdata <- dc_biplot.qza$data$Vectors %>% left_join(study1meta)

## save as 'biplot_dc_select2016'; export at 925x925
dc_bigplot <- ggplot() +
  geom_text(data=dc_fullplotdata, aes(x=PC1, y=PC2, label=Site, color=Window), size=2.5) +
  scale_color_distiller(palette = "Spectral", 
                        name = "Sampling Window",
                        breaks = c(7, 10, 13, 16, 19, 22), 
                        labels = c("Apr-7", "May-20", "June-18", "July-30", "Sept-10", "Oct-22")) +
  theme(panel.background = element_rect(fill = "gray40"),
        plot.background = element_rect(fill = "gray95"),
        legend.position = "top",
        legend.key.width = unit(3.5, "cm")) +
  xlab(paste(paste("PC1", round(100*dc_biplot.qza$data$ProportionExplained[1],2),"%"))) +
  ylab(paste(paste("PC2", round(100*dc_biplot.qza$data$ProportionExplained[2],2),"%"))) +
  geom_segment(data=dc_biplot.qza$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                 top_n(8, a) %>% #keep 8 furthest away points
                 mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% # scale arrows linearly... is this ok? 
                 left_join(tax),
               aes(x=0, xend=PC1, y=0, yend=PC2),
               arrow = arrow(length = unit(0.3,"cm")))

## save as 'biplot_inset_dc_select2016'; export at 660x660
dc_insetplot <- ggplot() + 
  theme(panel.background = element_rect(fill = "gray90"), 
        plot.background = element_rect(fill = "gray90"),
        legend.position = "top") +
  geom_segment(data=dc_biplot.qza$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% 
                 top_n(8, a) %>% 
                 mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% 
                 left_join(tax),
               aes(x=0, xend=PC1, y=0, yend=PC2, color=Order),
               size=2,
               arrow = arrow(length = unit(0.3,"cm"))) +
  geom_label_repel(data=dc_biplot.qza$data$Species %>% 
                     mutate(a=sqrt(PC1^2+PC2^2)) %>% 
                     top_n(8, a) %>% 
                     mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% 
                     select(FeatureID, PC1, PC2) %>% 
                     left_join(tax),
                   aes(x=PC1, y=PC2, fill=Order, label=ASValias), force=5, 
                   segment.alpha = 0.1, color="black", size=5) +
  scale_fill_manual(values = c('#9A6600', '#9BCC94', '#405E00')) + ## to match other barplots!
  scale_color_manual(values = c('#9A6600', '#9BCC94', '#405E00')) +
  labs(x="", y="")

################################################################################
## repeat for plots that are filtered to require at least 10 reads per sample
################################################################################

## read in data
dcmin10_biplot_url <- "/Users/do/Repos/nhguano/data/qiime_qza/biplots/deicode_biplot_min10feat.qza"
dcmin10_biplot.qza <- read_qza(dcmin10_biplot_url)

## create data frame for species plot
dcmin10_fullplotdata <- dcmin10_biplot.qza$data$Vectors %>% left_join(study1meta)

## save as 'biplot_dcmin10_select2016'; export at 925x925
dcmin10_bigplot <- ggplot() +
  geom_text(data=dcmin10_fullplotdata, aes(x=PC1, y=PC2, label=Site, color=Window), size=2.5) +
  scale_color_distiller(palette = "Spectral", 
                        name = "Sampling Window",
                        breaks = c(7, 10, 13, 16, 19, 22), 
                        labels = c("Apr-7", "May-20", "June-18", "July-30", "Sept-10", "Oct-22")) +
  theme(panel.background = element_rect(fill = "gray40"),
        plot.background = element_rect(fill = "gray95"),
        legend.position = "top",
        legend.key.width = unit(3.5, "cm")) +
  xlab(paste(paste("PC1", round(100*dcmin10_biplot.qza$data$ProportionExplained[1],2),"%"))) +
  ylab(paste(paste("PC2", round(100*dcmin10_biplot.qza$data$ProportionExplained[2],2),"%"))) +
  geom_segment(data=dcmin10_biplot.qza$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                 top_n(8, a) %>% #keep 8 furthest away points
                 mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% # scale arrows linearly... is this ok? 
                 left_join(tax),
               aes(x=0, xend=PC1, y=0, yend=PC2),
               arrow = arrow(length = unit(0.3,"cm")))

## save as 'biplot_inset_dcmin10_select2016'; export at 660x660
dcmin10_insetplot <- ggplot() + 
  theme(panel.background = element_rect(fill = "gray90"), 
        plot.background = element_rect(fill = "gray90"),
        legend.position = "top") +
  geom_segment(data=dcmin10_biplot.qza$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% 
                 top_n(8, a) %>% 
                 mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% 
                 left_join(tax),
               aes(x=0, xend=PC1, y=0, yend=PC2, color=Order),
               size=2,
               arrow = arrow(length = unit(0.3,"cm"))) +
  geom_label_repel(data=dcmin10_biplot.qza$data$Species %>% 
                     mutate(a=sqrt(PC1^2+PC2^2)) %>% 
                     top_n(8, a) %>% 
                     mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>% 
                     select(FeatureID, PC1, PC2) %>% 
                     left_join(tax),
                   aes(x=PC1, y=PC2, fill=Order, label=ASValias), force=5, 
                   segment.alpha = 0.1, color="black", size=5) +
  scale_fill_manual(values = c('#9A6600', '#9BCC94', '#405E00')) + ## to match other barplots!
  scale_color_manual(values = c('#9A6600', '#9BCC94', '#405E00')) +
  labs(x="", y="")


################################################################################
## Assessing overlaps in diet trends in 2015 and 2016
## Fewer sites sampled in both years, with data in fewer overlapping sampling windows
################################################################################

## import read data (importing just the read-filtered, min-sample data)
read_dat = read_csv("https://github.com/devonorourke/nhguano/raw/master/data/filtered_dataset_wWindows_min1000Reads_min2SamplesperSiteWindow.csv")

## What sites have best overlaps?
study2_Window_sumry <- read_dat %>% 
  mutate(Year=ifelse(StudyID == "oro16", "2016", "2015")) %>% 
  mutate(Labeler=paste(Site,Year,sep="-")) %>% 
  group_by(Labeler, Window) %>% 
  summarise(Samples = n_distinct(SampleID))

ggplot(data = study2_Window_sumry, aes(x=Window, y=Labeler, label=Samples)) + geom_text() + theme_bw()
  ## just three sites have datasets where there are >1 samples per Site/Window continually

## make list of site weeks to extract:
foxwindow <- c(14:17)
corwindow <- c(14:17)
hopwindow <- c(12:15)

foxlist <- read_dat %>% filter(Site=="FOX" & Window %in% foxwindow) %>% pull(SampleID) %>% unique()
  ## just 29 samples over those 4 weeks...
corlist <- read_dat %>% filter(Site == "COR" & Window %in% corwindow) %>% pull(SampleID) %>% unique()
  ## 45 samples 
hoplist <- read_dat %>% filter(Site == "HOP" & Window %in% hopwindow) %>% pull(SampleID) %>% unique()
  ## 51 samples 
study2list <- c(foxlist, corlist, hoplist)

## filter dataset to include just those samples
study2_dat <- read_dat %>% filter(SampleID %in% study2list)

#rm(read_dat, study2_Window_sumry)
rm(list=ls(pattern="*list*"))
rm(list=ls(pattern="*window*"))

################################################################################
## bar plot for overlapping sites and windows
################################################################################

barplot_study2_plotdat <- study2_dat %>% 
  mutate(Year=ifelse(StudyID == "oro16", "2016", "2015")) %>% 
  mutate(Labeler=paste(Site,Year,sep="-")) %>% 
  group_by(Labeler, Window, WindowStart, Order) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(pReads = Reads/sum(Reads))

## 12 color palette for the unique Orders here (there are 12 Orders again just like in 2016)
pal12 <- c('#fffe71', '#CE9834', '#9A6600', '#9BCC94', '#405E00', '#993303', 
           '#FF6501', '#336799', '#9ACEFF', 'gray75', '#D14A89', 'gray25')

## values of colors are (yellow, lightbrown, darkbrown, lightgreen, darkgreen, darkred, 
##                       orange, darkblue, lightblue, lightgray, purple, darkgray)

## rearrange the levels for the facet wrap here to match previous ordering by Site
barplot_study2_plotdat$Site <- factor(barplot_study1_plotdat$Site, levels = c(
  "HOL", "FOX", "HOP", "PEN", "CNA", "CNB", "MTV", "MAP", "BRN"))

## plot is faceting Sites, showing each Window as separate stacked bars of Order relative abundances
## plot barplot; save as 'Order_relAbund_stackedBars_select2016'; export at 800x800
ggplot(data = barplot_study2_plotdat,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  #scale_x_continuous(breaks=c(7,10,13,16,19,22),
  #                   labels=c('2016-03-26', '2016-05-07', '2016-06-18', '2016-07-30', '2016-09-10', '2016-10-22')) +
  scale_y_continuous(breaks = c(0,.5, 1)) +
  facet_grid(Labeler ~ .) +
  labs(x="", y="Relative Abundance (%)", fill="Order") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


################################################################################
## there are so few samples per Site/Week overlap that a lot of these differences..
## ..are likely due to undersampling.
## What about evaluating how many taxa are detected year over year, regardless of sampling Window?
################################################################################

## First, examining how many overlapping taxa per Site, year over year
oro15sites <- read_dat %>% filter(StudyID=="oro15") %>% distinct(Site)
oro16sites <- read_dat %>% filter(StudyID=="oro16") %>% distinct(Site)
sitesbothyears <- intersect(oro15sites, oro16sites) %>% pull()
rm(oro15sites, oro16sites)

## get list of ASVs for each distinct Site and Year, then find intersetions of ASVs for each Site between years
siteyearasv <- function(location){
  tmp15.asvs <- read_dat %>% filter(Site == location, StudyID == "oro15")
  tmp16.asvs <- read_dat %>% filter(Site == location, StudyID == "oro16")
  commonASVs = intersect(tmp15.asvs$ASValias, tmp16.asvs$ASValias) %>% data.frame(commonASVs = .)
  diffASVs = setdiff(tmp15.asvs$ASValias, tmp16.asvs$ASValias) %>% data.frame(diffASVs = .)
  nSamples = read_dat %>% filter(Site == location) %>% summarise(Samples=n_distinct(SampleID)) %>% pull(Samples)
  common_Order = read_dat %>% filter(ASValias %in% commonASVs$commonASVs & Site == location) %>% 
    group_by(Order, ASValias) %>% summarise(Reads=sum(Reads), 
                                            Samples=n_distinct(SampleID)) %>% 
    mutate(pReads = round(100*(Reads/sum(Reads)),2), 
           pSamples=round(100*(Samples/nSamples),2),
           Status = 'shared')
  diff_Order = read_dat %>% filter(ASValias %in% diffASVs$diffASVs & Site == location) %>% 
    group_by(Order, ASValias) %>% summarise(Reads=sum(Reads), 
                                            Samples=n_distinct(SampleID)) %>% 
    mutate(pReads = round(100*(Reads/sum(Reads)),2), 
           pSamples=-round(100*(Samples/nSamples),2),
           Status = 'distinct')
  common_Genus = read_dat %>% filter(ASValias %in% commonASVs$commonASVs & Site == location) %>% 
    group_by(Genus, ASValias) %>% summarise(Reads=sum(Reads), Samples=n_distinct(SampleID)) %>% 
    mutate(pReads = round(100*(Reads/sum(Reads)),2),
           pSamples=round(100*(Samples/nSamples),2),
           Status = 'shared')
  diff_Genus = read_dat %>% filter(ASValias %in% diffASVs$diffASVs & Site == location) %>% 
    group_by(Genus, ASValias) %>% summarise(Reads=sum(Reads), Samples=n_distinct(SampleID)) %>% 
    mutate(pReads = round(100*(Reads/sum(Reads)),2),
           pSamples=-round(100*(Samples/nSamples),2),
           Status = 'distinct')
  tmp_Order = rbind(common_Order, diff_Order) %>% mutate(Site = location)
  tmp_Genus = rbind(common_Genus, diff_Genus) %>% mutate(Site = location)
  list(tmp_Order, tmp_Genus)
}

hopcomps <- siteyearasv("HOP")
foxcomps <- siteyearasv("FOX")
mapcomps <- siteyearasv("MAP")
brncomps <- siteyearasv("BRN")
mascomps <- siteyearasv("MAS")
gilcomps <- siteyearasv("GIL")
corcomps <- siteyearasv("COR")

Order_plotdat <- rbind(hopcomps[[1]], foxcomps[[1]], mapcomps[[1]], brncomps[[1]], mascomps[[1]], gilcomps[[1]], corcomps[[1]])
Genus_plotdat <- rbind(hopcomps[[2]], foxcomps[[2]], mapcomps[[2]], brncomps[[2]], mascomps[[2]], gilcomps[[2]], corcomps[[2]])

## just keep the select Orders that we repeatedly discuss in the paper (also the most prevalent)
keepOrders <- c("Coleoptera", "Diptera", "Ephemeroptera", "Lepidoptera", "Trichoptera")

## drop MAS site because it's undersampled:
Order_plotdat <- Order_plotdat %>% filter(Site != "MAS")

## set all Samples for distinct to negative values for barplot
Order_plotdat$Samples = ifelse(Order_plotdat$Status=="distinct", -Order_plotdat$Samples, Order_plotdat$Samples)

## ordering remaining sites by geographic position:
Order_plotdat$Site <- factor(Order_plotdat$Site, levels = c("COR", "GIL", "FOX", "HOP", "MAP", "BRN"))

## plot barchart; save as 'ASVs_OrdersObserved_yearOverYear'; export at 750x750
ggplot(Order_plotdat %>% filter(Order %in% keepOrders), 
       aes(x=Order, y=Samples, color=Status, label = ASValias)) +
  geom_jitter(width = 0.2) +
  facet_wrap(~ Site, nrow=3) +
  scale_y_continuous(limits = c(-40, 54), 
                     breaks = c(-40, 0, 40), 
                     labels = c("40", '0', "40")) +
  labs(x="", y="Samples with ASV", color="ASVs observed in 2015 and 2016") +
  scale_color_manual(values = c('#dfc27d', '#018571')) +
  theme_bw(base_size = 16) +
  theme(legend.position = "top", axis.text.x = element_text(size = 11, angle = 22.5, hjust=1))

################################################################################
## among those highly prevalent ASVs; what's their per-sample read distribution look like?
################################################################################
## summarise the number of samples and reads per Order
ASVsumry <- read_dat %>% 
  group_by(ASValias, Order, Family, Genus, Species) %>% 
  summarise(Samples=n_distinct(SampleID), Reads=sum(Reads))

## select the top 5 most prevalent ASVs for each Order
top5ASVs_perOrder <- ASVsumry %>% 
  group_by(Order) %>% 
  top_n(5, Samples) %>% 
  filter(Samples > 10)

## calculate the faction of reads (rather than plotting absolute abundances)
read_dat <- read_dat %>% 
  group_by(SampleID) %>% 
  mutate(fracReads=100*(Reads/sum(Reads)))

## use same 12-color palette:
pal12 <- c('#fffe71', '#CE9834', '#9A6600', '#9BCC94', '#405E00', '#993303', 
           '#FF6501', '#336799', '#9ACEFF', 'gray75', '#D14A89', 'gray25')

## plot the read abundances of those 72 ASVs, faceting by Order 
## save as 'ASVabundances_perSample_perOrder_top5detected_notrarefied'; export at 1050x800
ggplot(read_dat %>% filter(ASValias %in% top5ASVs_perOrder$ASValias),
       aes(x=ASValias, y=fracReads, color=Order)) +
  geom_jitter(width = 0.2) +
  facet_wrap(~Order, nrow=3, scales = 'free_x') +
  scale_color_manual(values = pal12) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(size =7.5), legend.position = "top") +
  labs(x="", y="Sequence counts (%)", color="", caption = "sequence counts calculated using raw data (not rarefied)")

## do we see the same patterns for rarefied data?
rare_qza <- read_qza("~/Repos/nhguano/data/qiime_qza/ASVtable/sampleOnly_rfyd_table.qza")
mat.tmp <- rare_qza$data
rm(rare_qza)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
rare_df <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
#rm(df.tmp)
colnames(rare_df) <- c("ASVid", "SampleID", "Reads")
rare_df <- rare_df %>% 
  filter(ASVid %in% read_dat$ASVid) %>% 
  filter(SampleID %in% read_dat$SampleID)
## apply taxonomy info 
taxmrg <- read_dat %>% select(ASVid, ASValias, Order) %>% distinct(ASVid, ASValias, Order)
rare_df <- merge(rare_df, taxmrg, all.x = TRUE)
rare_df$fracReads = 100*(rare_df$Reads/1000)

## save as 'ASVabundances_perSample_perOrder_top5detected_rarefiedData'; export at 1050x800
ggplot(rare_df %>% filter(ASValias %in% top5ASVs_perOrder$ASValias),
       aes(x=ASValias, y=fracReads, color=Order)) +
  geom_jitter(width = 0.2) +
  facet_wrap(~Order, nrow=3, scales = 'free_x') +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_color_manual(values = pal12) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(size =7.5), legend.position = "top") +
  labs(x="", y="Sequence counts (%)", color="", caption = "sequence counts calculated using rarefied data")

## now we'll select the top 50 most deeply sequenced samples
sampleSumry <- read_dat %>% 
  group_by(SampleID, Site, Window) %>% 
  summarise(Reads=sum(Reads))

topSamples <- sampleSumry %>% ## selecting the three most abundant samples per Site from raw sequence counts
  arrange(-Reads) %>% 
  group_by(Site, Window) %>% 
  top_n(1, Reads) %>%   ## sample with most reads per site per week
  group_by(Site) %>% 
  top_n(3, Reads)   ## top 3 most abundant samples per site
  ## min has > 20k reads. okay to use as is.

## show the proportions of reads rather than absolute abundances instead
read_dat <- read_dat %>% 
  group_by(SampleID) %>% 
  mutate(fracReads=Reads/sum(Reads))

## save as 'ASVabundances_perSample__top3perSite'; export at 800x800
ggplot(read_dat %>% filter(SampleID %in% topSamples$SampleID),
                               aes(x=SampleID, y=Reads, color=Order)) +
  geom_jitter(data=read_dat %>% filter(SampleID %in% topSamples$SampleID & Order == "Coleoptera"), ## plot beetles first so we don't hide other bugs
              aes(x=SampleID, y=fracReads, color=Order), width=0.2) +
  geom_jitter(data=read_dat %>% filter(SampleID %in% topSamples$SampleID & Order != "Coleoptera"),
              aes(x=SampleID, y=fracReads, color=Order), width=0.2) +
  scale_color_manual(values = pal12) +
  scale_y_continuous(breaks = c(0, .5, 1), labels = c(0, 50, 100)) +
  theme_bw(base_size = 16) +
  facet_wrap(~Site, ncol=6, scales = "free_x") +
  theme(legend.position = "top", axis.text.x = element_blank()) +
  labs(x="Samples", y="Sequence counts (%)", color="")
  

