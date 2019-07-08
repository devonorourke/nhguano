
# Overview
We initially processed all samples following steps described in the [sequence_processing](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md) document, and performed a detailed analysis of potential contamination among negative control and positive control samples as outlined in the [decontam_workflow](https://github.com/devonorourke/nhguano/blob/master/docs/decontam_workflow.md) file. These results suggest that a very minor amount of contamination exists, though it isn't clear whether this is due to either the DNA extraction or sequencing processes (or a combination of both). Fortunately our analyses suggest that these are extremely rare events, occur randomly with respect to which ASV is detected in a particular well, and occur with such low abundances that any diversity estimate that utilizes abundance information will not likely be biased due to a contamination event. However, it is also true that we would expect unweighted abundance metrics to be more sensitive to these false positive events. Within-sample observed richness is likely to be marginally higher, and between-sample dissimilarities are likely to be lower than a dataset that would otherwise have not generated any contaminants.  We therefore acknowledge that some amount of contamination is possible among samples, but it is likely minor and proceeded without removing any particular ASVs other than the ASVs known to be assigned to the positive control samples (i.e. the biological mock community).

The control-removed ASV table was then filtered so that bat (host) DNA was removed. We used a pair of databases and classified data using an alignment approach in VSEARCH.


1c. Removing selected ASVs
1ci. Are these all M. lucifigus? Create table of putative North American bat when species name is provided. Hoping to make statment that we suspect our analyses  are not influenced by other species.
1cii. Find any supporting evidence of previous info about emergence counts in those homes?

2ai. Filtering data for taxonomic specificity. Do in QIIME with some sort of feature filter thing. Create new filtered table and repseq file with arthonly data.
2aii. Create the summarize table at the same time to quickly figure out what a range of sampling depths are going to work best to preserve the most samples. Then run rarefying alpha viz to figure out what level to choose.

2b. Import that table into R and run ??

# Determining sampling depth for normalization
4. Figure out sampling depth with alpha rarefaction viz. Using default 10 iterations, but setting 1000 reads as minimum.

```
qiime diversity alpha-rarefaction \
  --p-metrics observed_otus --p-min-depth 1000 --p-max-depth 10000 --p-iterations 10 \
  --i-table study.raw_table.qza --o-visualization study.raw.alphaRareViz.qzv
```
