library(qiime2R)
library(tidyverse)
library(reshape2)

# ## import abundance info
# ## import the original ASV table (not clustered):
## change path below to your own path!
# download.file(url = "https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/ASVtable/tmp.raw_table.qza",
#               destfile = "~/Desktop/tmp.raw_table.qza")

qzapath = "~/github/nhguano/data/qiime_qza/ASVtable/tmp.raw_table.qza"  ## amend this line to your own path!
featuretable <- read_qza(qzapath)
mat.asv <- featuretable$data
rm(featuretable, qzapath)

## convert the mat.tmp object into a long format
df.asv <- as.data.frame(mat.asv)
rm(mat.asv)
df.asv$ASVid <- rownames(df.asv)
rownames(df.asv) <- NULL
asv_table_long <- melt(df.asv, id = "ASVid") %>% filter(value != 0)
rm(df.asv)
colnames(asv_table_long) <- c("ASVid", "SampleID", "Reads")

## create new header per ASV:
asv_sumry_long <- asv_table_long %>% 
  group_by(ASVid) %>% 
  summarise(sumReads = sum(Reads)) %>% 
  mutate(newLabel = paste0(ASVid, ";size=", sumReads))

## import the .tsv version of the fasta 
## modify this path to update github path!!!
fasta_tsv <- read_delim(file = "~/github/nhguano/data/text_tables/asv_data/allSamps_ASVseqs.tsv.gz",
                        delim = "\t", col_names = FALSE)
colnames(fasta_tsv) <- c("ASVid", "sequence")
## merge with new ASV header 
asv_fasta_renamed <- merge(asv_sumry_long, fasta_tsv, by="ASVid") %>% 
  select(-ASVid, -sumReads) %>% 
  mutate(newLabel = paste0(">", newLabel))
## write as .csv file and convert into fasta file for reuse with vsearch 
write_csv(asv_fasta_renamed,
          file = "~/github/nhguano/data/text_tables/asv_data/allSamps_ASVseqs_wSizes.csv.gz")
