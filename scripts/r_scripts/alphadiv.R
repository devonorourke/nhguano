library(tidyverse)
library(qiime2R)
library(gganimate)
library(lubridate)
library(ggpubr)
library(scales)

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
## subset metadata to include just these values
study1meta <- meta %>% filter(Site %in% study1list & StudyID == "oro16")

## will summarize by WOY, but also wanted to explore if we could better define things..
## ..by a Window span of something other than defining WOY as "starting on Sunday"..
## because of discontinous collection by volunteers (some started Sunday, others Tuesday, others Friday, etc)
## Defining a Window of 10 days to capture the most number of overlapping samples seemed to work best

## summary comparing these two approaches follows after this section

################################################################################
## date manipulations
## do these before any substantial plots to define what the Window span will be for dates (if not using WOY)
## 10-day spans seem to work pretty well...
## for whatever metadata file you start with, you can manipulate the date range as follows
## using just the select 9 sites in 2016 data as example
################################################################################

## what's the earliest time we recorded info for 2016 data?
mindate = study1meta %>% summarise(date=min(newDate)) %>% pull(date)
## what's the latest in the season?
maxdate = study1meta %>% summarise(date=max(newDate)) %>% pull(date)
## what's that span of dates?
dateinterval <- interval(mindate, maxdate)
as.duration(dateinterval) ## just over 29 weeks studied in 2016 (more than half the year!)
## how many days (or seconds/minutes, etc...) between min/max
time_length(dateinterval, unit = "day") ## 204 days

## if you already subsetted your data (like in the "selectmeta" object)
## breaking into N-day length bins
## use the min/maxdate interval information above to determine bins
## setting at 15 bins (so about every 2 weeks)
##notrun: study1meta <- meta %>% filter(Site %in% study1list & StudyID == "oro16")
study1meta$DOY <- yday(study1meta$newDate)
minDOY = min(study1meta$DOY)  ## 98
maxDOY = max(study1meta$DOY)  ## 302
allDOY=maxDOY-minDOY+1 
allDOY ## 205 days total (a span of 204 days)

## get range defined with the number of bins (devined in the "by" term of the 'breaks' argument)
window=10         ## 'window' represents number of days to span a given observation)
buckets = ceiling(allDOY/window)    ## 'buckets' represents number of bins needed to span duration of study
study1meta$DayRange = cut(study1meta$DOY, breaks = c(seq(from=minDOY, to=maxDOY+window, by=window)), include.lowest = TRUE)
study1meta$Window <- as.numeric(study1meta$DayRange)

study1meta$RangeStart <- as.character(study1meta$DayRange)
study1meta$RangeStart <- gsub("\\(", "", study1meta$RangeStart)
study1meta$RangeStart <- gsub("\\]", "", study1meta$RangeStart)
study1meta$RangeStart <- gsub("\\[", "", study1meta$RangeStart)
study1meta <- study1meta %>% separate(., col = RangeStart, into=c("RangeStart", "RangeEnd"), sep = ",")

## create a dataframe of the range days and convert back into a true date to merge with metadata
possibleRangeVals <- c(study1meta$RangeStart, study1meta$RangeEnd) %>% unique()
rangematcher_df <- data.frame(calendardate=as.Date(mindate) + 0:(allDOY+window-1), WindowRange=seq(minDOY, (maxDOY+window)))
## merge into metadata to define dates of window span
study1meta <- merge(study1meta, rangematcher_df, by.x='RangeStart', by.y='WindowRange')
study1meta <- study1meta %>% select(-RangeStart) %>% rename(DateWindowStart=calendardate)
study1meta <- merge(study1meta, rangematcher_df, by.x='RangeEnd', by.y='WindowRange')
study1meta <- study1meta %>% select(-RangeEnd) %>% rename(DateWindowEnd=calendardate)
## write this file to disk
write.csv(study1meta, file="~/Repos/nhguano/data/metadata/study1meta.csv", row.names = FALSE)

################################################################################
## Summary comparing WOY vs. Window strategy
################################################################################
study1_WOY_sumry <- study1meta %>% select(Site, WOY, SampleID) %>% 
  group_by(Site, WOY) %>% summarise(Samples=n_distinct(SampleID))

study1_Win_sumry <- study1meta %>% select(Site, Window, SampleID) %>%
  group_by(Site, Window) %>% summarise(Samples=n_distinct(SampleID))

## here's what WOY looks like
p1_woy <- ggplot(data = study1_WOY_sumry, aes(x=WOY, y=Site, label=Samples)) + geom_text() + theme_bw()
## here's what Window-approach (10 day bins) looks like
p1_win <- ggplot(data = study1_Win_sumry, aes(x=Window, y=Site, label=Samples)) + geom_text() + theme_bw()

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

study1_WOY_sumry_filt <- study1meta %>% select(Site, WOY, SampleID) %>% filter(SampleID %in% filt5_samples) %>% 
  group_by(Site, WOY) %>% summarise(Samples=n_distinct(SampleID))
study1_Win_sumry_filt <- study1meta %>% select(Site, Window, SampleID) %>% filter(SampleID %in% filt5_samples) %>% 
  group_by(Site, Window) %>% summarise(Samples=n_distinct(SampleID))

p1_woy_f <- ggplot(data = study1_WOY_sumry_filt, aes(x=WOY, y=Site, label=Samples)) + geom_text() + theme_bw()
p1_win_f <- ggplot(data = study1_Win_sumry_filt, aes(x=Window, y=Site, label=Samples)) + geom_text() + theme_bw()
ggarrange(p1_woy_f, p1_win_f, ncol=2, labels=c("A", "B"))

## it looks like when we have a 5-ASV requirement, two of the original nine sites loose a lot of data; dropping these
dropstudy1_sites <- c("BRN", "HOP")

## try plotting again:
p1_woy_f <- ggplot(data = study1_WOY_sumry_filt %>% filter(!Site %in% dropstudy1_sites), 
                   aes(x=WOY, y=Site, label=Samples)) + geom_text() + theme_bw()
p1_win_f <- ggplot(data = study1_Win_sumry_filt %>% filter(!Site %in% dropstudy1_sites),  
                   aes(x=Window, y=Site, label=Samples)) + geom_text() + theme_bw()
ggarrange(p1_woy_f, p1_win_f, ncol=2, labels=c("A", "B"))

## The most continuous overlaps in sampling dates exist for 5 sites: HOL, FOX, MAP, CHI, and MTV
## Most sites > 20 miles apart (just one pair are > 10 miles apart (MTV and MAP)); can rule out overlapping populations
## These are the best candidates for our timeseries studies so we'll test these first

keepstudy1_sites <- c("FOX", "HOL", "MAP", "MTV", "CHI")

rm(p1_woy_f, p1_woy, p1_win, p1_win_f, study1_Win_sumry, study1_WOY_sumry, rangematcher_df)
################################################################################
## Are there distinct differences in alpha diversity per site per Window of time?
## Run ANOVA to test for main effects of Site and Time (Window)
################################################################################

## merge updated metadata (including Window) with alpha data
alpha_study1dat <- merge(study1meta, alpha_dat)
## filtering out singletons and dropping two sites with largely discontinuous datasets (BRN and HOP)
alpha_study1dat_filt <- alpha_study1dat %>% filter(SampleID %in% filt5_samples & Site %in% keepstudy1_sites)
## collect just these names for future use in data filtering (in other scripts):
alpha_study1dat_names <- alpha_study1dat_filt %>% select(SampleID) %>% unique()
colnames(alpha_study1dat) <- "#OTU ID"
write.table(alpha_study1dat, file="~/Repos/nhguano/data/metadata/alpha_study1names.txt", row.names = FALSE, quote = FALSE)

alphaaov <- function(alphametric) {
  aov.out <- aov(Value ~ Site * Window, data=alpha_study1dat_filt %>% filter(Metric == alphametric))
}

alpha_aov_ob <- aov(Value ~ Site * Window, data=alpha_study1dat_filt %>% filter(Metric == "observed"))
capture.output(summary(alpha_aov_ob),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_ob.txt")  ## Site sig, Window not, interaction sig
alpha_aov_sh <- aov(Value ~ Site * Window, data=alpha_study1dat_filt %>% filter(Metric == "shannon"))
capture.output(summary(alpha_aov_sh),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_sh.txt")  ## ## Site sig, Window not, interaction not
alpha_aov_fp <- aov(Value ~ Site * Window, data=alpha_study1dat_filt %>% filter(Metric == "faithpd"))
capture.output(summary(alpha_aov_fp),file="~/Repos/nhguano/data/stats/alpha/aov_alpha_aov_fp.txt")  ## ## Site sig, Window not, interaction sig
  ## Window significant for observed, but not for Shannon's or Faith's... 

## Tukey posthoc summaries
posthoc_ob <- TukeyHSD(x = alpha_aov_ob, 'Site')
tmp_tuk_ob <- data.frame(posthoc_ob[[1]]) %>% mutate(Pairs=row.names(.)) %>% arrange(p.adj) %>% mutate(p.adj=round(p.adj, 3)) %>% 
  mutate(diff=round(diff,2)) %>% mutate(lwr=round(lwr, 3)) %>% mutate(upr=round(upr,3))
write.table(tmp_tuk_ob, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_ob.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

posthoc_sh <- TukeyHSD(x = alpha_aov_sh, 'Site')
tmp_tuk_sh <- data.frame(posthoc_sh[[1]]) %>% mutate(Pairs=row.names(.)) %>% arrange(p.adj) %>% mutate(p.adj=round(p.adj, 3)) %>% 
  mutate(diff=round(diff,2)) %>% mutate(lwr=round(lwr, 3)) %>% mutate(upr=round(upr,3))
write.table(tmp_tuk_sh, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_sh.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

posthoc_fp <- TukeyHSD(x = alpha_aov_fp, 'Site')
tmp_tuk_fp <- data.frame(posthoc_fp[[1]]) %>% mutate(Pairs=row.names(.)) %>% arrange(p.adj) %>% mutate(p.adj=round(p.adj, 3)) %>% 
  mutate(diff=round(diff,2)) %>% mutate(lwr=round(lwr, 3)) %>% mutate(upr=round(upr,3))
write.table(tmp_tuk_fp, file="~/Repos/nhguano/data/stats/alpha/alpha_tuk_fp.tsv", row.names = FALSE, quote=FALSE, sep = "\t") 

## ANOVA suggests no significant effect for temporal variation, but Site is significant, though so is interaction..
  ##..so the impact of Site still depends on the Window (time)
## Pairwise Tukey tests suggest that MAP is different from other sites for "observed" and "faith" but not "shannon"..
  ##..and this likely is because of a few samples in early May that had a much higher diversity than the two other..
  ##..locations sampled at that time (while the other 2 weren't sampled at all)
  ## Therefore this effect is likely localized
## Need to vissualize to make this make sense (see next section)

################################################################################
## Let's plot how these are changing by Site and Time (for each Window)
## Generating a static plot is going to be tough to see all at once..
## ..so we'll also generate a matching animation for the static one
################################################################################

## static boxplot for each Site+Time window for selected 2016 sites
## making separate plots to use distinct axes

## color palette:
pal5 <- c("#ba495b","#56ae6c","#9350a1","#b0923b","#697cd4")

## plot boxplot; save as 'alpha_ob_boxplotPerWindowPerSite_select2016'; export at 800x800
ggplot(data=alpha_study1dat_filt %>% filter(Metric=="observed"), 
             aes(x=Window, y=Value, group=Window, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
  facet_wrap(~ Site, ncol=1) +
  scale_x_continuous(breaks=c(1,6,11,16,21),
                     labels=c("2016-04-07", "2016-05-27", "2016-07-16", "2016-09-04", "2016-10-24	")) +
  scale_y_continuous(trans="log2", labels = comma_format(accuracy = 1)) + 
  labs(x="", y="Observed richness") +
  theme_bw() + theme(legend.position = "none")

## plot boxplot; save as 'alpha_sh_boxplotPerWindowPerSite_select2016'; export at 800x800
ggplot(data=alpha_study1dat_filt %>% filter(Metric=="shannon"), 
       aes(x=Window, y=Value, group=Window, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
  facet_wrap(~ Site, ncol=1) +
  scale_x_continuous(breaks=c(1,6,11,16,21),
                     labels=c("2016-04-07", "2016-05-27", "2016-07-16", "2016-09-04", "2016-10-24	")) +
  labs(x="", y="Shannon's entropy") +
  theme_bw() + theme(legend.position = "none")

## plot boxplot; save as 'alpha_fp_boxplotPerWindowPerSite_select2016'; export at 800x800
ggplot(data=alpha_study1dat_filt %>% filter(Metric=="faithpd"), 
       aes(x=Window, y=Value, group=Window, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
  facet_wrap(~ Site, ncol=1) +
  scale_x_continuous(breaks=c(1,6,11,16,21),
                     labels=c("2016-04-07", "2016-05-27", "2016-07-16", "2016-09-04", "2016-10-24	")) +
  labs(x="", y="Faith's PD") +
  theme_bw() + theme(legend.position = "none")

## animate the observed ASVs; note we've switched the x axis to Sites, and are transitioning through the Windows of date
## first one is for observed ASVs
ob_ani1 <- ggplot(data=alpha_study1dat_filt %>% filter(Metric=="observed"), 
       aes(x=Site, y=Value, group=Site, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
  scale_y_continuous(trans="log2", labels = comma_format(accuracy = 1)) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(x="", y="Observed richness", title = 'Date: {closest_state}') +
  ease_aes('quadratic-in-out') +
  transition_states(DateWindowStart, 1, 3)
ob_rendered <- animate(ob_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/alpha_ob.gif", animation = ob_rendered)
rm(ob_ani1, ob_rendered)

## next one for Shannon's entropy alpha vals
sh_ani1 <- ggplot(data=alpha_study1dat_filt %>% filter(Metric=="shannon"), 
                  aes(x=Site, y=Value, group=Site, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(x="", y="Shannon's entropy", title = 'Date: {closest_state}') +
  ease_aes('quadratic-in-out') +
  transition_states(DateWindowStart, 1, 3)
sh_rendered <- animate(sh_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/alpha_sh.gif", animation = sh_rendered)
rm(sh_ani1, sh_rendered)

## last one for Faiths's PD alpha vals
fp_ani1 <- ggplot(data=alpha_study1dat_filt %>% filter(Metric=="faithpd"), 
                  aes(x=Site, y=Value, group=Site, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.15) +
  scale_color_manual(values=pal5) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(x="", y="Faith's PD", title = 'Date: {closest_state}') +
  ease_aes('quadratic-in-out') +
  transition_states(DateWindowStart, 1, 3)
fp_rendered <- animate(fp_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/alpha_fp.gif", animation = fp_rendered)
rm(fp_ani1, fp_rendered)

rm(alpha_dat)
################################################################################
## Generating stacked bar plots for these 5 sites
## Examining change in abundance among taxa; grouped by taxonomic Order
## Producing both static and animation equivalent
################################################################################

## generate dataset:
## don't care about minimum ASV filtering here; including any sample in Site+Window
## get names of 5 sites
barplot_study1_samps <- alpha_study1dat %>% filter(Site %in% keepstudy1_sites) %>% pull(SampleID) %>% unique()
## import data with read counts
read_dat <- read_csv(file="~/Repos/nhguano/data/filtered_dataset.csv")
barplot_study1_dat <- read_dat %>% filter(SampleID %in% barplot_study1_samps)
rm(read_dat)
## add the date Window information here...
study1meta_reduced <- study1meta %>% select(SampleID, DOY, Window, DateWindowStart, DateWindowEnd)
barplot_study1_dat <- merge(barplot_study1_dat, study1meta_reduced, by="SampleID")
rm(study1meta_reduced)

## calculate per-Order read abundance; grouping all samples with shared Site+Window
barplot_study1_plotdat <- barplot_study1_dat %>% 
  group_by(Site, Window, DateWindowStart, Order) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(pReads = Reads/sum(Reads))

## 12 color palette for the unique Orders here
pal12 <- c('#fffe71', '#CE9834', '#9A6600', '#9BCC94', '#405E00', '#993303', 
           '#FF6501', '#336799', '#9ACEFF', 'gray75', '#D14A89', 'gray25')

## values of colors are (yellow, lightbrown, darkbrown, lightgreen, darkgreen, darkred, 
##                       orange, darkblue, lightblue, lightgray, purple, darkgray)

## plot is faceting Sites, showing each Window as separate stacked bars of Order relative abundances
## plot barplot; save as 'Order_relAbund_stackedBars_select2016'; export at 
ggplot(data = barplot_study1_plotdat,
       aes(x=Window, y=pReads, group=Window, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal12) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15,17,19,21),
                     labels=c("2016-04-07", "2016-04-27", "2016-05-17", "2016-06-06", "2016-06-26", "2016-07-16", "2016-08-05", "2016-08-25", "2016-09-14", "2016-10-04", "2016-10-24")) +
  facet_wrap(~Site, ncol=1) +
  labs(x="", y="Relative Abundance (%)", fill="Taxonomic\nOrder") +
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
  transition_states(DateWindowStart, 1, 3)

bar_rendered <- animate(bar_ani1, duration = 21, fps = 9)
anim_save(filename = "~/Repos/nhguano/figures/anis/taxaOrder_stackedBars_select2016.gif", animation = bar_rendered)
rm(bar_ani1, bar_rendered)
