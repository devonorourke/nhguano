library(tidyverse)
library(qiime2R)
library(gganimate)
library(lubridate)
library(ggpubr)
library(scales)
library(geofacet)

## importing read counts from all samples sequenced in study 
## this was the last file output from the 'decontam' R script
## see script here:
## see documentation of that process here:

## import read data
read_dat = read_csv("~/Repos/nhguano/data/filtered_dataset_wWindows.csv")

################################################################################
################################################################################
## Alpha diversity assessments
################################################################################
################################################################################


## using alpha diversity estimates calculated from QIIME 
## see https://github.com/devonorourke/nhguano/blob/master/docs/diversity_analyses.md

## import diversity files
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

rm(ob_url, fp_url, sh_url)
################################################################################
## going to split up dataset into 2016 and 2015 components and focus on selected sites
## Data is too discontinuous to leave all Sites for year/over/year comparisons
## Some sites too underrepresented to include for meaningful statistical analyses
## Subsetting 2016 and specific sites with majority of sampling first
################################################################################

## import metadata again
metadata <- as.data.frame(read_csv(file="~/Repos/nhguano/data/metadata/nhbat_meta.csv"))
tmp_alphadat <- merge(alpha_dat, metadata) %>% select(-SampleType, -Date)
meta <- metadata %>% filter(SampleID %in% tmp_alphadat$SampleID)
rm(tmp_alphadat)

## Collective summary of all samples (to determine which sites to keep for 2016, etc.)
alpha_sumry <- tmp_alphadat %>% 
  filter(Metric=="observed") %>% 
  group_by(Site, StudyID) %>% 
  summarise(Samples=n_distinct(SampleID))
  ## 9 sites have > 150 samples/site for 2016; selecting these for initial alpha diversity analyses 
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

## focusing on 2016 samples at these sites:
study1list <- c("FOX", "HOL", "MAP", "HOP", "CHI", "CNB", "MTV", "EPS", "BRN")
## update metadata to include just these Sites, 
study1meta <- meta %>% filter(Site %in% study1list & StudyID == "oro16")
## now and add in Window information
tmp_metamrg <- read_dat %>% 
  filter(SampleID %in% study1meta$SampleID) %>% 
  select(SampleID, Window, WindowStart) %>% 
  distinct()
study1meta <- merge(study1meta, tmp_metamrg)
rm(tmp_metamrg)
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
alpha_dat %>% filter(Metric=="observed" & Value > 0) %>% nrow() ## there are 2,489 to start
alpha_dat %>% filter(Metric=="observed" & Value > 1) %>% nrow() ## filtering out singleton ASVs drops it down to 2,215 samples
alpha_dat %>% filter(Metric=="observed" & Value > 2) %>% nrow() ## no doubletons means 2,056
alpha_dat %>% filter(Metric=="observed" & Value > 3) %>% nrow() ## 1,904
alpha_dat %>% filter(Metric=="observed" & Value > 4) %>% nrow() ## 1,743
alpha_dat %>% filter(Metric=="observed" & Value > 5) %>% nrow() ## 1,573
alpha_dat %>% filter(Metric=="observed" & Value > 10) %>% nrow() ## 866 when we require at least 10 ASVs per sample

filt2_samples <- alpha_dat %>% filter(Metric=="observed" & Value > 1) %>% pull(SampleID)
filt5_samples <- alpha_dat %>% filter(Metric=="observed" & Value >= 5) %>% pull(SampleID)
filt10_samples <- alpha_dat %>% filter(Metric=="observed" & Value >= 10) %>% pull(SampleID)

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

## there are several samples with just a single ASV; dropping these samples but retaining any other sample with > 1 ASV
## Because we subsampled, many of the rare variants are likely lost; these wouldn't contribute to the overall read abundances in bar plots...
## ...but likely would have big impacts in diversity measures that don't incorporate abundance information

##notrun: keepstudy1_sites <- c("FOX", "HOL", "MAP", "MTV", "CHI")
rm(p1_woy_f, p1_woy, p1_win, p1_win_f, study1_Window_sumry, study1_WOY_sumry, rangematcher_df,
   study1_win_sumry_filt, study1_WOY_sumry_filt)

################################################################################
## Are there distinct differences in alpha diversity per site per Window of time?
## Run ANOVA to test for main effects of Site and Time (Window)
################################################################################

## merge updated metadata (including Window) with alpha data
alpha_study1dat <- merge(study1meta, alpha_dat)
## filtering out singletons; keeping all 9 2016 study sites of interest:
alpha_study1dat_filt_noSingles <- alpha_study1dat %>% filter(SampleID %in% filt2_samples & Site %in% study1list)
## set our categorical data for ANOVA as factors
alpha_study1dat_filt_noSingles$Window <- as.factor(alpha_study1dat_filt_noSingles$Window)
alpha_study1dat_filt_noSingles$Site <- as.factor(alpha_study1dat_filt_noSingles$Site)
class(alpha_study1dat_filt_noSingles$Window)

## collect just these names for future use in data filtering (in other scripts):
alpha_study1dat_names <- alpha_study1dat_filt_noSingles %>% select(SampleID) %>% unique()
colnames(alpha_study1dat_names) <- "#OTU ID"
write.table(alpha_study1dat_names, file="~/Repos/nhguano/data/metadata/alpha_study1names.txt", row.names = FALSE, quote = FALSE)

alpha_aov_ob <- aov(Value ~ Site * Window, data=alpha_study1dat_filt_noSingles %>% filter(Metric == "observed"))
capture.output(summary(alpha_aov_ob),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_ob.txt")  ## Site sig, Window not, interaction sig
alpha_aov_sh <- aov(Value ~ Site * Window, data=alpha_study1dat_filt_noSingles %>% filter(Metric == "shannon"))
capture.output(summary(alpha_aov_sh),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_sh.txt")  ## ## Site sig, Window not, interaction not
alpha_aov_fp <- aov(Value ~ Site * Window, data=alpha_study1dat_filt_noSingles %>% filter(Metric == "faithpd"))
capture.output(summary(alpha_aov_fp),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_fp.txt")  ## ## Site sig, Window not, interaction sig
  ## Window significant for observed, but not for Shannon's or Faith's... 

## Tukey posthoc summaries
tukeyfunction <- function(anovafile, whichfactor){
  posthoc_tmp_site <- TukeyHSD(x = anovafile, whichfactor)
  data.frame(posthoc_tmp_site[[1]]) %>% mutate(Pairs=row.names(.)) %>% arrange(p.adj) %>% mutate(p.adj=round(p.adj, 3)) %>% 
    mutate(diff=round(diff,2)) %>% mutate(lwr=round(lwr, 3)) %>% mutate(upr=round(upr,3))
}

tuk_ob_site <- tukeyfunction(alpha_aov_ob, 'Site')
tuk_ob_window <- tukeyfunction(alpha_aov_ob, 'Window')
write.table(tuk_ob_site, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_ob_site.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 
write.table(tuk_ob_window, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_ob_window.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

tuk_sh_site <- tukeyfunction(alpha_aov_sh, 'Site')
tuk_sh_window <- tukeyfunction(alpha_aov_sh, 'Window')
write.table(tuk_sh_site, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_sh_site.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 
write.table(tuk_sh_window, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_sh_window.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

tuk_fp_site <- tukeyfunction(alpha_aov_fp, 'Site')
tuk_fp_window <- tukeyfunction(alpha_aov_fp, 'Window')
write.table(tuk_fp_site, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_fp_site.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 
write.table(tuk_fp_window, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_fp_window.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

rm(list=ls(pattern="tuk_*"))
rm(list=ls(pattern="alpha_aov*"))
rm(ob_dat, sh_dat, fp_dat)
rm(filt2_samples, filt5_samples, filt10_samples)
rm(alphaaov, alphaimport, alpha_sumry)
################################################################################
## Let's plot how these are changing by Site and Time (for each Window)
## Generating a static plot is going to be tough to see all at once..
## ..so we'll also generate a matching animation for the static one
################################################################################

## static boxplot for each Site+Time window for selected 2016 sites
## color palette:
#pal5 <- c("#ba495b","#56ae6c","#9350a1","#b0923b","#697cd4")
# pal9 <- ??

## plot boxplot; save as 'alpha_ob_boxplotPerWindowPerSite_select2016'; export at 800x1200
ggplot(data=alpha_study1dat_filt_noSingles %>% filter(Metric=="observed"), 
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

## plot boxplot; save as 'alpha_sh_boxplotPerWindowPerSite_select2016'; export at 800x1200
ggplot(data=alpha_study1dat_filt_noSingles %>% filter(Metric=="shannon"), 
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
ggplot(data=alpha_study1dat_filt_noSingles %>% filter(Metric=="faithpd"), 
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
ob_ani1 <- ggplot(data=alpha_study1dat_filt_noSingles %>% filter(Metric=="observed"), 
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
sh_ani1 <- ggplot(data=alpha_study1dat_filt_noSingles %>% filter(Metric=="shannon"), 
                  aes(x=Site, y=Value, group=Site, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(x="", y="Shannon's entropy", title = 'Date: {closest_state}') +
  ease_aes('quadratic-in-out') +
  transition_states(WindowStart, 1, 3)
sh_rendered <- animate(sh_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/alpha_sh.gif", animation = sh_rendered)
rm(sh_ani1, sh_rendered)

## last one for Faiths's PD alpha vals
fp_ani1 <- ggplot(data=alpha_study1dat_filt_noSingles %>% filter(Metric=="faithpd"), 
                  aes(x=Site, y=Value, group=Site, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
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

## import data with read counts; retain only relevant samples to 2016
barplot_study1_dat <- read_csv(file="~/Repos/nhguano/data/filtered_dataset_wWindows.csv") %>% 
  filter(SampleID %in% alpha_study1dat$SampleID)

## calculate per-Order read abundance; grouping all samples with shared Site+Window
barplot_study1_plotdat <- barplot_study1_dat %>% 
  group_by(Site, Window, WindowStart, Order) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(pReads = Reads/sum(Reads))

## 12 color palette for the unique Orders here
pal12 <- c('#fffe71', '#CE9834', '#9A6600', '#9BCC94', '#405E00', '#993303', 
           '#FF6501', '#336799', '#9ACEFF', 'gray75', '#D14A89', 'gray25')

## values of colors are (yellow, lightbrown, darkbrown, lightgreen, darkgreen, darkred, 
##                       orange, darkblue, lightblue, lightgray, purple, darkgray)

## can't exactly order with coordinates, but we can rearrange the levels for the facet wrap here..
## ..to better match their relative positions in geographic space
barplot_study1_plotdat$Site <- factor(barplot_study1_plotdat$Site, levels = c(
  "HOL", "FOX", "HOP", "CNB", "CHI", "EPS", "MTV", "BRN", "MAP"))


## plot is faceting Sites, showing each Window as separate stacked bars of Order relative abundances
## plot barplot; save as 'Order_relAbund_stackedBars_select2016'; export at 800x1000
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
## this works really well as it shows how often the same sites have similar abundances across sampling windows...
## going to make another animation that shows the random assortment of these sites..
  ## ..to use in presentation (this emphasizes that ordering by location makes interpreting these data a lot easier!)

barplot_study1_plotdat_falseorder <- barplot_study1_plotdat
barplot_study1_plotdat_falseorder$Site <- factor(barplot_study1_plotdat_falseorder$Site, levels = c(
  "BRN","CNB","CHI","EPS", "FOX", "HOL", "HOP", "MAP", "MTV"))
bar_ani2 <- ggplot(data = barplot_study1_plotdat_falseorder, aes(x=Site, y=pReads, group=Site, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  labs(x="", y="Relative Abundance (%)", fill="Taxonomic\nOrder",
       title = 'Sampling Date: {closest_state}') +
  theme_bw(base_size = 16) +
  ease_aes('quadratic-in-out') +
  transition_states(WindowStart, 1, 3)

bar_rendered2 <- animate(bar_ani2, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/taxaOrder_stackedBars_select2016_randomOrder.gif", animation = bar_rendered2)
rm(bar_ani2, bar_rendered2)

################################################################################
## Setting up similar animation, but overlaying these bar charts on map
################################################################################

## creating a custom plot using geofacet library.
## generating the custom grid using this site as reference: https://hafen.github.io/grid-designer/

## this script was used to figure out how to use existing lat/lon for each site and convert to a grid
## it maintains a 1/10th of a degree as a cell, so spacing in resulting geofacet grid is still basically proportional
read_dat %>% 
  select(Site, SiteLat, SiteLong) %>% 
  distinct() %>% 
  mutate(tmprow=SiteLat*10, tmpcol=-(SiteLong)*10) %>% 
  mutate(row=abs(ceiling(tmprow-(max(tmprow)))) + 1) %>% 
  mutate(col=(ceiling(max(tmpcol)-tmpcol))+1) %>% 
  select(Site, row, col)

## I looked at the output of the above code, pasted in a spreadsheet to figure out which row/col's were duplicate, and adjusted
## the following row/cols were derived from that spreadsheet, but I modified two pairs of sites (BRN/MAP and CNA/CNB) that had overlapping grids

mygrid <- data.frame(
  name = c("Holderness", "Cornish", "CanterburyW", "CanterburyE", "Penacook", "Chichester", "Alsted", "Hillsborough","Hopkinton", 
           "Epsom", "Rollins", "Gilsum", "Antrim", "Massabesic", "Swanzey", "Greenfield", "MontVernon","Wilton", "HollisN", "HollisS"),
  code = c("HOL", "COR", "CNA", "CNB", "PEN", "CHI", "ALS", "FOX", "HOP", 
           "EPS", "ROL", "GIL", "WLD", "MAS", "SWZ", "GRN", "MTV", "WLT", "BRN", "MAP"),
  row = c(1, 3, 4, 4, 5, 5, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9, 9, 10, 10, 11),
  col = c(9, 1, 10, 11, 9, 12, 2, 6, 8, 12, 17, 2, 5, 12, 3, 7, 9, 8, 9, 9),
  stringsAsFactors = FALSE
)

## For the 2016 dataset, the plots were too spread out to see the barplots easily, so I modified the distances between sites
## these are NOT well attributed to their geographic coordinates, but their relative positions are maintained
my2016grid <- data.frame(
  name = c("Holderness", "CanterburyE", "Chichester", "Hillsborough","Hopkinton", "Epsom", "MontVernon", "HollisN", "HollisS"),
  code = c("HOL", "CNB", "CHI", "FOX", "HOP", "EPS", "MTV", "BRN", "MAP"),
  row = c(1,3,3,4,4,4,5,6,6),
  col = c(3,3,4,1,2,4,3,3,4))
  
## add the 'code' name to each 
barplot_study1_plotdat$code <- barplot_study1_plotdat$Site

## save static plot as 'geofacet_barplots_2016sites'; export at 
ggplot(data = barplot_study1_plotdat,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  #scale_x_continuous(breaks=c(10,13,16,19,22,25,28,31),
  #                   labels=c('2016-04-01', '2016-05-01', '2016-05-31', '2016-06-30', '2016-07-30', '2016-08-29', '2016-09-28', '2016-10-28')) +
  scale_y_continuous(breaks = c(0,.5, 1)) +
  facet_geo(~ code, grid = my2016grid) +
  labs(x="", y="Relative Abundance (%)", fill="Order") +
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "top") +
  guides(fill = guide_legend(nrow = 2))


ggplot() +
    geom_bar(data = barplot_study1_plotdat, 
           aes(x=Site, y=pReads, group=Site, fill=Order),
           stat="identity") +
  scale_fill_manual(values=pal12) +
  labs(x="", y="", fill="Order",
       title = 'Sampling Date: {closest_state}') +
  theme_bw(base_size = 16) +
  ease_aes('quadratic-in-out') +
  transition_states(WindowStart, 1, 3)

rm(barplot_study1_dat, barplot_study1_plotdat, barplot_study1_plotdat_falseorder, mygrid, my_grid, my2016grid, select2016grid)
rm(alpha_dat, alpha_study1dat, alpha_study1dat_filt_noSingles, alpha_study1dat_names, bar_ani2)

################################################################################
################################################################################
## Beta diversity assessments for 2016 selected sites
################################################################################
################################################################################


  
