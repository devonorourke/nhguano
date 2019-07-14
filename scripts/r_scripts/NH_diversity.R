library(tidyverse)
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

## import read data
read_dat = read_csv("https://github.com/devonorourke/nhguano/raw/master/data/filtered_dataset_wWindows.csv")

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
metadata <- metadata %>% filter(SampleID %in% tmp_alphadat$SampleID)

## Collective summary of all samples (to determine which sites to keep for 2016, etc.)
alpha_sumry <- tmp_alphadat %>% 
  filter(Metric=="observed") %>% 
  group_by(Site, StudyID) %>% 
  summarise(Samples=n_distinct(SampleID))
  ## 5 sites have > 50 samples/site for 2016; selecting these for initial alpha diversity analyses 
  ## See subset below

## Collecitve summary across _all_ sites (not just subset)
alpha_stats <- tmp_alphadat %>% 
  group_by(Site, StudyID, Metric) %>% 
  summarise(meanValue=mean(Value), medianValue=median(Value), sdValue=sd(Value)) %>% 
  mutate(meanValue=round(meanValue, 2), sdValue=round(sdValue, 2)) %>% 
  arrange(Metric, Site, 
          factor(StudyID, levels = c("oro15", "oro16")))

write.csv(alpha_stats, file="~/Repos/nhguano//data/stats/alpha/alpha_stats.csv",
          row.names = FALSE, quote = FALSE)
rm(alpha_stats)

## focusing on 2016 samples at these sites:
study1list <- c("FOX", "HOL", "MAP", "HOP", "BRN", "MTV")
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
  ## given the low overall diversity, having more samples per bin is likely helpful at capturing more diversity per timepoint than WOY

## What if we filter out samples with very few ASVs?
alpha_dat %>% filter(Metric=="observed" & Value > 0) %>% nrow() ## there are 911 to start
alpha_dat %>% filter(Metric=="observed" & Value > 1) %>% nrow() ## filtering out singleton ASVs drops it down to 905 samples
alpha_dat %>% filter(Metric=="observed" & Value > 2) %>% nrow() ## no doubletons means 880
alpha_dat %>% filter(Metric=="observed" & Value > 3) %>% nrow() ## 851
alpha_dat %>% filter(Metric=="observed" & Value > 4) %>% nrow() ## 806
alpha_dat %>% filter(Metric=="observed" & Value > 5) %>% nrow() ## 760
alpha_dat %>% filter(Metric=="observed" & Value > 10) %>% nrow() ## 491 when we require at least 10 ASVs per sample

## replot the same data where we're removing samples with just 1 ASV
filt2_samples <- alpha_dat %>% filter(Metric=="observed" & Value > 1) %>% pull(SampleID)

study1_WOY_sumry_filt <- read_dat %>% 
  filter(SampleID %in% study1meta$SampleID & SampleID %in% filt2_samples) %>%   ## vary this line for one of above 3 filt groups
  group_by(Site, WOY) %>% 
  summarise(Samples = n_distinct(SampleID))
  
study1_win_sumry_filt <- read_dat %>% 
  filter(SampleID %in% study1meta$SampleID & SampleID %in% filt2_samples) %>%   ## vary this line for one of above 3 filt groups
  group_by(Site, Window) %>% 
  summarise(Samples = n_distinct(SampleID))

p1_woy_f <- ggplot(data = study1_WOY_sumry_filt, aes(x=WOY, y=Site, label=Samples)) + geom_text() + theme_bw()
p1_win_f <- ggplot(data = study1_win_sumry_filt, aes(x=Window, y=Site, label=Samples)) + geom_text() + theme_bw()
ggarrange(p1_woy_f, p1_win_f, ncol=2, labels=c("A", "B"))

## try plotting again:
p1_woy_f <- ggplot(data = study1_WOY_sumry_filt %>% filter(!Site %in% dropstudy1_sites), 
                   aes(x=WOY, y=Site, label=Samples)) + geom_text() + theme_bw()
p1_win_f <- ggplot(data = study1_Win_sumry_filt %>% filter(!Site %in% dropstudy1_sites),  
                   aes(x=Window, y=Site, label=Samples)) + geom_text() + theme_bw()
ggarrange(p1_woy_f, p1_win_f, ncol=2, labels=c("A", "B"))

## need to drop instances where there is just one sample in a given Site+Window:
SiteWindowKeep <- study1_win_sumry_filt %>% 
  mutate(SiteWindow = paste0(Site,Window)) %>% 
  filter(Samples > 1) %>% 
  pull(SiteWindow)

study1meta <- study1meta %>% 
  mutate(SiteWindow = paste0(Site,Window)) %>% 
  filter(SiteWindow %in% SiteWindowKeep)

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
capture.output(summary(alpha_aov_ob),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_ob.txt")  ## Site sig, Window not, interaction sig
alpha_aov_sh <- aov(Value ~ aovSite * aovWindow, data=alpha_study1dat %>% filter(Metric == "shannon"))
capture.output(summary(alpha_aov_sh),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_sh.txt")  ## ## Site sig, Window not, interaction not
alpha_aov_fp <- aov(Value ~ aovSite * aovWindow, data=alpha_study1dat %>% filter(Metric == "faithpd"))
capture.output(summary(alpha_aov_fp),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_fp.txt")  ## ## Site sig, Window not, interaction sig
  ## Window significant for observed, but not for Shannon's or Faith's... 

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
  "HOL", "FOX", "HOP", "MTV", "MAP", "BRN"))

## plot boxplot; save as 'alpha_ob_boxplotPerWindowPerSite_select2016'; export at 500x1000
ggplot(data=alpha_study1dat %>% filter(Metric=="observed"), 
       aes(x=Window, y=Value, group=Window), fill="gray50") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  facet_grid(Site ~ .) +
  scale_x_continuous(breaks=c(10,13,16,19,22,25,28,31),
                     labels=c('2016-04-01', '2016-05-01', '2016-05-31', '2016-06-30', '2016-07-30', '2016-08-29', '2016-09-28', '2016-10-28')) +
  scale_y_continuous(trans="log2", labels = comma_format(accuracy = 1)) + 
  labs(x="", y="Observed richness") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))

## plot boxplot; save as 'alpha_sh_boxplotPerWindowPerSite_select2016'; export at 500x1000
ggplot(data=alpha_study1dat %>% filter(Metric=="shannon"), 
       aes(x=Window, y=Value, group=Window), fill="gray50") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  facet_grid(Site ~ .) +
  scale_x_continuous(breaks=c(10,13,16,19,22,25,28,31),
                     labels=c('2016-04-01', '2016-05-01', '2016-05-31', '2016-06-30', '2016-07-30', '2016-08-29', '2016-09-28', '2016-10-28')) +
  labs(x="", y="Shannon's entropy") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))


## plot boxplot; save as 'alpha_fp_boxplotPerWindowPerSite_select2016'; export at 800x1200
ggplot(data=alpha_study1dat %>% filter(Metric=="faithpd"), 
       aes(x=Window, y=Value, group=Window), fill="gray50") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  facet_grid(Site ~ .) +
  scale_x_continuous(breaks=c(10,13,16,19,22,25,28,31),
                     labels=c('2016-04-01', '2016-05-01', '2016-05-31', '2016-06-30', '2016-07-30', '2016-08-29', '2016-09-28', '2016-10-28')) +
  labs(x="", y="Faith's PD") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))

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
## Generating stacked bar plots for these 2016 sites
## Producing both static and animation equivalent
################################################################################
################################################################################

## calculate per-Order read abundance; grouping all samples with shared Site+Window

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
  "HOL", "FOX", "HOP", "MTV", "MAP", "BRN"))

## plot is faceting Sites, showing each Window as separate stacked bars of Order relative abundances
## plot barplot; save as 'Order_relAbund_stackedBars_select2016'; export at 800x800
ggplot(data = barplot_study1_plotdat,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  scale_x_continuous(breaks=c(10,13,16,19,22,25,28,31),
                     labels=c('2016-04-01', '2016-05-01', '2016-05-31', '2016-06-30', '2016-07-30', '2016-08-29', '2016-09-28', '2016-10-28')) +
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
# my2016grid <- data.frame(
#  name = c("Holderness", "CanterburyE", "Chichester", "Hillsborough","Hopkinton", "Epsom", "MontVernon", "HollisN", "HollisS"),
#  code = c("HOL", "CNB", "CHI", "FOX", "HOP", "EPS", "MTV", "BRN", "MAP"),
#  row = c(1,3,3,4,4,4,5,6,6), col = c(3,3,4,1,2,4,3,3,4))
  
## add the 'code' name to each 
# barplot_study1_plotdat$code <- barplot_study1_plotdat$Site

## save static plot as 'geofacet_barplots_2016sites'; export at 
#ggplot(data = barplot_study1_plotdat,
#       aes(x=Window, y=pReads, group=Window, fill=Order)) +
#  geom_bar(stat="identity") +
#  scale_fill_manual(values=pal12) +
#  scale_y_continuous(breaks = c(0,.5, 1)) +
#  facet_geo(~ code, grid = my2016grid) +
#  labs(x="", y="Relative Abundance (%)", fill="Order") +
#  theme_bw() +
#  theme(axis.text.x = element_blank(), legend.position = "top") +
#  guides(fill = guide_legend(nrow = 2))


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

ds_adonis_df <- adonisfunction(ds.dist_url, "ds")
bc_adonis_df <- adonisfunction(bc.dist_url, "bc")
uu_adonis_df <- adonisfunction(uu.dist_url, "uu")
wu_adonis_df <- adonisfunction(wu.dist_url, "wu")

## write files to disk:
write_csv(ds_adonis_df, path = "~/Repos/nhguano/data/stats/beta/ds_adonis_select2016.csv", quote=FALSE)
write_csv(bc_adonis_df, path = "~/Repos/nhguano/data/stats/beta/bc_adonis_select2016.csv", quote=FALSE)
write_csv(uu_adonis_df, path = "~/Repos/nhguano/data/stats/beta/uu_adonis_select2016.csv", quote=FALSE)
write_csv(wu_adonis_df, path = "~/Repos/nhguano/data/stats/beta/wu_adonis_select2016.csv", quote=FALSE)

rm(ds.dist_url, bc.dist_url, uu.dist_url, wu.dist_url, 
   ds_adonis_df, bc_adonis_df, uu_adonis_df, wu_adonis_df)

## import PCoA files for plotting
ds.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_ds.qza"
bc.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_bc.qza"
uu.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_uu.qza"
wu.pcoa_url <- "~/Repos/nhguano/data/qiime_qza/pcoa/select2016/s16_pcoa_wu.qza"

pcoaplotfunction <- function(urlpath, betametric){
  tmp.qza <- read_qza(urlpath)
  tmp.data <- tmp.qza$data$Vectors[,1:3]
  tmp.lab.pc1 <- tmp.qza$data$ProportionExplained[1] %>% pull() %>% round(., 2) %>% paste("PC1", ., sep = " ") %>% paste0(., "%")
  tmp.lab.pc2 <- tmp.qza$data$ProportionExplained[2] %>% pull() %>% round(., 2) %>% paste("PC2", ., sep = " ") %>% paste0(., "%")
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
  ggplot(data = tmp_df, aes(x=PC1, y=PC2, color=Site, label=Window)) +
    geom_text() + 
    scale_x_continuous(limits = c(-.63,.63)) +
    scale_y_continuous(limits = c(-.63,.63)) +
    labs(x=pc1lab, y=pc2lab, subtitle = CaptionName) +
    theme_bw(base_size = 16) + theme(legend.position = "top")
}

## gather individual plots
ds_pcoa_plot <- pcoaplotvizfunction(ds_pcoa_df, "Dice-Sorensen")
bc_pcoa_plot <- pcoaplotvizfunction(bc_pcoa_df, "Bray-Curtis")
uu_pcoa_plot <- pcoaplotvizfunction(uu_pcoa_df, "unweighted Unifrac")
wu_pcoa_plot <- pcoaplotvizfunction(wu_pcoa_df, "weighted Unifrac")

## plot as one big set; save as 'pcoa_select2016'; export at 900x900
ggarrange(ds_pcoa_plot, bc_pcoa_plot, uu_pcoa_plot, wu_pcoa_plot,
          ncol=2, nrow=2, common.legend = TRUE)

rm(list=ls(pattern="*_pcoa_*"))
rm(list=ls(pattern="*_url"))
