Filtering out non NH-bat samples, NTCs, mock samples.
```
## filtering table
qiime feature-table filter-samples \
  --i-table all.raw_table.qza --o-filtered-table study.raw_table.qza \
  --m-metadata-file "$META" --p-where "SampleID='subject-1'"

## filtering repSeqs
qiime feature-table filter-seqs \
  --i-data all.raw_repSeqs.qza --i-table study.raw_table.qza --o-filtered-data study.raw_repSeqs.qza
```

# Overview
1. Remove bat host dna, keep only arthropod data with family level info

1a. Setting up host db, classifying with VSEARCH
1b. Setting up big db, classifying with VSEARCH and again with Naive Bayes
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
