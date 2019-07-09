library(tidyverse)
library(formattable)

################################################################################
## the same dataset was classified using either:
# A. VSEARCH (97% identity, 89% coverage)
# B. Naive Bayes classifier (default settings in QIIME 2)

## First two plots explore:
#1. How often do the two classifiers agree (same taxon name per level for each ASV)
#2. How often do the two classifiers assign some information or not (per ASV, per taxon level)
################################################################################

## how often do the two bigDB classifiers match for the same taxon?
## import VSEARCH taxonomy
vstaxa <- read_delim(file="~/Repos/nhguano/data/tax/tmp.raw_bigDB_VStax.tsv", delim="\t")
vstaxa <- vstaxa %>% separate(., col = Taxon, 
                              sep=';', into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub(".__", "", y)))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub("^$|^ $", NA, y)))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
vstaxa <- as.data.frame(apply(vstaxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(vstaxa)[1] <- "ASVid"

## import Naive Bayes taxonomy
nbtaxa <- read_delim(file="~/Repos/nhguano/data/tax/tmp.raw_bigDB_nbtax.tsv", delim="\t")
nbtaxa <- nbtaxa %>% separate(., col = Taxon, 
                              sep=';', into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub(".__", "", y)))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub("^$|^ $", NA, y)))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
nbtaxa <- as.data.frame(apply(nbtaxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(nbtaxa)[1] <- "ASVid"

## common and different values per level:
matchfunction <- function(TaxonLevel){
  vstmp <- vstaxa %>% select(ASVid, TaxonLevel) %>% filter(complete.cases(.)) %>% unite(., "vsinput", c("ASVid", TaxonLevel)) %>% pull(vsinput)
  nbtmp <- nbtaxa %>% select(ASVid, TaxonLevel) %>% filter(complete.cases(.)) %>% unite(., "nbinput", c("ASVid", TaxonLevel)) %>% pull(nbinput)
  Match <- length(intersect(vstmp, nbtmp))
  Diff <- length(setdiff(vstmp, nbtmp))
  data.frame(TaxonLevel, Match, Diff)
}

## table calculating number of times they disagree/match, but note that these are ignoring comparisons where one might be NA
CompTable <- rbind(matchfunction("Kingdom"), matchfunction("Phylum"), matchfunction("Class"),
                   matchfunction("Order"), matchfunction("Family"), matchfunction("Genus"), matchfunction("Species"))

## plot table as image; save as 'taxcomp_matchStatus'; export at 275x300
formattable(CompTable)

## how often does a classifier assign a name or not?
emptyfunction <- function(TaxonLevel){
  vstmp <- vstaxa %>% select(ASVid, TaxonLevel) %>% rename(vtax=TaxonLevel) 
  nbtmp <- nbtaxa %>% select(ASVid, TaxonLevel) %>% rename(ntax=TaxonLevel)
  alltmp <- merge(vstmp, nbtmp)
  alltmp$Status <- ""
  alltmp$Status <- ifelse(is.na(alltmp$vtax) & is.na(alltmp$ntax), "both missing", alltmp$Status)
  alltmp$Status <- ifelse(is.na(alltmp$vtax) & !is.na(alltmp$ntax), "vsearch missing only", alltmp$Status)
  alltmp$Status <- ifelse(!is.na(alltmp$vtax) & is.na(alltmp$ntax), "nbayes missing only", alltmp$Status)
  alltmp$Status <- ifelse(!is.na(alltmp$vtax) & !is.na(alltmp$ntax), "neither missing", alltmp$Status)
  alltmp %>% group_by(Status) %>% tally() %>% spread(Status, n) %>% mutate(Level=TaxonLevel)
}

EmptyTable <- rbind(emptyfunction("Phylum"), emptyfunction("Class"),emptyfunction("Order"), 
                    emptyfunction("Family"), emptyfunction("Genus"), emptyfunction("Species"))
EmptyTable <- EmptyTable %>% select(Level, `both missing`, `neither missing`, `nbayes missing only`, `vsearch missing only`)
tmpk <- data.frame(emptyfunction("Kingdom"), "nbayes missing only"=0)
colnames(tmpk) <- c("both missing", "neither missing", "vsearch missing only", "Level", "nbayes missing only")
tmpk <- tmpk %>% select(Level, `both missing`, `neither missing`, `nbayes missing only`, `vsearch missing only`)
EmptyTable <- rbind(tmpk, EmptyTable)

## plot table as image; save as 'taxcomp_missingStatus'; export at 650x350
formattable(EmptyTable, 
            caption="Instances where one or both classifiers do not assign name to an ASV")

rm(CompTable, EmptyTable, tmpk)

################################################################################
## For all the instances in which Naive Bayes assigns a name where VSEARCH doesn't, what was the % identity for that VSEARCH alignment?
## In other words, how frequently was VSEARCH really close to making a call but was discarded because it was below the % identity threshold?
################################################################################

## subset just those instnaces where Nbayes makes a determination and VSEARCH does not
## focusing on only those instances where Naive Bayes had at least Family-name assigned (don't care if it only assigned to Phylum, or Class, or Order)
nbtmp <- nbtaxa %>% filter(!is.na(Family)) %>% select(ASVid) %>% pull()
vstmp <- vstaxa %>% filter(is.na(Family)) %>% select(ASVid) %>% pull()
asvqueries <- data.frame(intersect(vstmp, nbtmp))
colnames(asvqueries) <- "featureid"
write.table(asvqueries, file="~/Repos/nhguano/data/tax/vsearch_missingFaminfo_asvs.txt", 
            col.names = TRUE, quote = FALSE, row.names = FALSE)
## using that `asvqueries` object in standalone vsearch to generate the % alignments for each of these ASVs
## ran this code:
    ## vsearch --usearch_global $READS --db $REFS --id 0.8 --query_cov 0.89 --strand both \
    ## --maxaccepts 100 --threads 24 --blast6out vsearch_missingFam_vsearchOut.tsv
## importing that data output here to filter

famtaxtest <- read_delim(file='https://github.com/devonorourke/nhguano/raw/master/data/tax/vsearch_missingFam_vsearchOut.tsv.gz', 
                         delim='\t', col_names = FALSE)
colnames(famtaxtest) <- c("ASVid", "sid", "pid", "alnlen", "qlo", "qhi")
famtaxtest$qlen <- abs(famtaxtest$qlo - famtaxtest$qhi) +1  ## checking to see if all qcov are > 0.89
famtaxtest$qcov <- famtaxtest$alnlen / famtaxtest$qlen      ## yep.
famtaxtest$sid <- ifelse(is.na(famtaxtest$sid), "missingsomething", famtaxtest$sid) ## some subject ID's not getting reported here...

tmp <- famtaxtest %>% group_by(ASVid) %>% top_n(1, pid) %>% top_n(1, alnlen)


group %>% group_by(Subject) %>% top_n(1, pt)