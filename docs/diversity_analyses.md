# Table of Contents

- [Overview](#overview)
- [Clustering representative sequences (ASVs) into representative clusters (OTUs)](#clustering-representative-sequences--asvs--into-representative-clusters--otus-)
- [Classifying OTUs and filtering samples](#classifying-otus-and-filtering-samples)
  * [Database development](#database-development)
  * [Training a naive Bayes classifier](#training-a-naive-bayes-classifier)
  * [Classifying representative sequences](#classifying-representative-sequences)
  * [Filtering data for bat and arthropod-specific representative sequences](#filtering-data-for-bat-and-arthropod-specific-representative-sequences)
- [Diversity analysis summaries](#diversity-analysis-summaries)
  * [Overview](#overview-1)
  * [New Hampshire-wide analyses (all samples)](#new-hampshire-wide-analyses--all-samples-)
    + [Frequently detected arthropod orders](#frequently-detected-arthropod-orders)
    + [Frequently detected OTUs](#frequently-detected-otus)
    + [Pest analyses](#pest-analyses)
  * [Spatiotemporal analyses (select samples)](#spatiotemporal-analyses--select-samples-)
    + [Single site, single year, multiple sampling windows](#single-site--single-year--multiple-sampling-windows)
    + [Multiple sites, single year, multiple sampling windows](#multiple-sites--single-year--multiple-sampling-windows)
    + [Multiple sites, multiple years, single sampling window](#multiple-sites--multiple-years--single-sampling-window)

> <small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


# Overview
We initially processed all samples following steps described in the [Sequence Processing](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md) document, and evaluated potential sources of contamination among negative control and positive control samples as outlined in the [contamination workflow](https://github.com/devonorourke/nhguano/blob/master/docs/contamination_evaluation.md) file. We determined that ASVs specific to the mock community were to be removed from this dataset, while both positive and negative control are also to be discarded. This document outlines the remaining sequence processing leading into the diversity investigations, as well as brief summaries of the various components of the diversity work itself. Notably, within those brief summaries are links to the outputs of these steps (though all data produced is available within various directories in this repository).

The work presented herein:
1. Clustering representative sequences (ASVs) into representative clusters (OTUs)
2. Classifying OTUs, including:
  - Creating a custom database from BOLD references
  - Training a naive Bayes classifier
  - Assigning taxonomic information to representative sequences with VSEARCH and naive Bayes classifiers
  - Identifying bat host species among samples
  - Filtering OTUs to retain only arthropod-specific diet components
3. Diversity analysis summaries

The R script [seqProcessing.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/seqProcessing.R) was used to tie together some of the outputs from steps 1-4 above. Likewise, the [diversityAnalyses.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/diversityAnalyses.R) R script was used to complete almost all sections of step 5 for the diversity analysis summaries. Additional R scripts are noted when appropriate for each of these steps in the following sections.

# Clustering representative sequences (ASVs) into representative clusters (OTUs)
We clustered the exact sequence variants ([tmp.raw_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza), produced at the end of the [Sequence Processing](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md) document) into representative sequence clusters. Our previous work with this and other bat diet datasets had shown that many ASVs often were assigned the exact same taxonomic information, which tends to increase the species richness estimates, and can lead to differences in community composition between groups when there likely is none (in our opinion, a bat doesn't know the differences between June bugs with a 1 base pair difference in their COI sequence). Clustering certainly reduces these diversity estimates, and may hide potential group differences, but this more conservative approach was a tradeoff we felt appropriate for our questions motivated by bat diet differences in space and time: if OTUs were different (rather than ASVs) it is very likely those differences are meaningful in terms of the kinds of taxa being detected.  

Because QIIME 2 did not have the functionality to cluster ASVs into OTUs using abundance information _and_ produce a temporary `.uc` file as output, we had to perform a little bit of reformatting to run VSEARCH manually. This involved exporting the ASV table, generating a plain text fasta file with abundance information (per ASV), then clustering with VSEARCH and generating a `.uc` file that identified the ASV to the group Centroid (OTU). That `.uc` file was then used in the [seqProcessing.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/seqProcessing.R) script to aggregate the raw reads from each ASV per sample into OTUs per sample. This process worked as follows:

First, we applied an R script [addingAbundanceInfoToDerepASVfasta.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/addingAbundanceInfoToDerepASVfasta.R) to read in the asv table QIIME object [tmp.raw_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza), then aggregate the reads per ASV. This resulted in a tab-separated file, [allSamps_ASVseqs_wSizes.csv.gz](https://github.com/devonorourke/nhguano/blob/master/data/text_tables/asv_data/allSamps_ASVseqs_wSizes.csv.gz), which served as the input to modify to create a fasta file for VSEARCH clustering. Because it was a column-based format instead of the traditional fasta file, we modified that file as follows:

```
zcat allSamps_ASVseqs_wSizes.csv.gz | \
awk 'NR > 1 {print $0}' | tr ',' '\n' | gzip \
> allSamps_ASVseqs_wSizes.fasta.gz
```

The [allSamps_ASVseqs_wSizes.fasta.gz](https://github.com/devonorourke/nhguano/blob/master/data/text_tables/asv_data/allSamps_ASVseqs_wSizes.fasta.gz) file then served as input into VSEARCH to cluster at a 98.5% identity. To collapse ASVs into OTUs and generate the `.uc` file:
```
vsearch --cluster_size allSamps_ASVseqs_wSizes.fasta.gz \
--id 0.985 --qmask none --xsize --threads 10 --minseqlength 1 --fasta_width 0 \
--centroids allSamps_clustered_p985.fasta \
--uc allSamps_clustered_p985.uc

gzip *.uc
```

This produced the , which served as input to the [seqProcessing.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/seqProcessing.R) script to collapse each ASV into it's proper OTU (and aggregate read counts).  

Next, we applied the [seqProcessing.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/seqProcessing.R) R script to first aggregate the sequence count information within the [tmp.raw_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) ASV table into an OTU table using the [allSamps_clustered_p985.uc.gz](https://github.com/devonorourke/nhguano/tree/master/data/ucfile/allSamps_clustered_p985.uc.gz) file.  In that same R script, we next removed mock community samples, exact mock community ASVs, and negative control samples. The next step was to identify which OTUs were associated with bat diet (we restricted those OTUs classified to arthropods with at least family-level information), and which OTUs were likely derived from a bat host. This required classifying each OTU, explained in the next section.

# Classifying OTUs and filtering samples
## Database development
Representative sequences were classified using a custom database curated with reference sequences and taxonomic information obtained from the Barcode of Life Database ([BOLD](https://v4.boldsystems.org/index.php)), and used in an earlier study [Robeson et al., in press 2020](https://www.biorxiv.org/content/10.1101/2020.10.05.326504v1.full). We relied in particular on a series of tools available in the [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) QIIME 2 plugin to filter the raw reference sequences, as well as [Seqkit](https://github.com/shenwei356/seqkit). Briefly, we first obtained the BOLD COI sequences using a custom R script, [bold_datapull_byGroup.R](https://github.com/devonorourke/COIdatabases/blob/master/scripts/bold_datapull_byGroup.R), which queried the BOLD API using the ['bold'](https://github.com/ropensci/bold) R package. Reference sequences were filtered for nucleotide ambiguity and homopolymer runs (no more than 5 degenerate bases per sequence or 12 homopolymers) with 'qiime rescript cull-seqs' and length (min 250 bp max 1600 bp) with 'qiime rescript filter-seqs-length', and dereplicated by applying a Least Common Ancestor method that gave preference to sequences represented most frequently with 'qiime rescript dereplicate --p-mode 'super' --p-derep-prefix'. Remaining sequences were trimmed to boundaries defined by our COI primer sequences by performing multiple sequence alignment of reference and primer sequences [MAFFT](https://mafft.cbrc.jp/alignment/software/). The remaining primer-trimmed, nucleotide quality and length-filtered references were once more filtered and gaps removed with 'qiime rescript degap-seqs' (retaining only reference sequences with a minimum length of 170 bp), then dereplicated a second time with the same LCA method used earlier with 'qiime rescript dereplicate'. The final database consists of 739,345 unique COI reference sequences and taxonomic labels.

The workflow is fully detailed in a QIIME 2 forum post [available here](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references), with all supporting scripts made available in a separate GitHub repository, [COIdatabases](https://github.com/devonorourke/COIdatabases). Database files are hosted through Open Science Framework in this [project folder](https://osf.io/d4jra/).  

## Training a naive Bayes classifier
We used a hybrid approach in classifying representative sequences, prioritizing exact matches from VSEARCH first, then retaining Naive Bayes classifications with sufficient information. We first trained the Naive Bayes classifier required to assign taxonomic information to representative sequences with the RESCRIPt command 'qiime rescript evaluate-fit-classifier', producing the QIIME object required for Naive Bayes classification. The sequence (bold_anml_seqs.qza) and taxonomy (bold_anml_taxa.qza) artifacts used as inputs in classifier training are available via Open Science Framework in this [project folder](https://osf.io/d4jra/), as is the 'bold_anml_classifier.qza' classifier object used as input for classifying the representative sequences in the next step of the workflow:
```
qiime rescript evaluate-fit-classifier \
  --i-sequences bold_anml_seqs.qza \
  --i-taxonomy bold_anml_taxa.qza \
  --p-reads-per-batch 6000 \
  --p-n-jobs 6 \
  --output-dir fitClassifier_boldANML
```
> The `bold_anml_classifier.qza` file is contained within the `fitClassifier_boldANML` output directory and is used in NBayes classification next.

## Classifying representative sequences

ASV sequences were classified using VSEARCH and naive Bayes separately. For VSEARCH, we required 100% identity across 94% query coverage as follows:
```
qiime feature-classifier classify-consensus-vsearch \
   --i-query tmp.raw_repSeqs.qza \
   --o-classification allASVs_VSp100c94_taxa.qza \
   --i-reference-reads bold_anml_seqs.qza \
   --i-reference-taxonomy bold_anml_taxa.qza \
   --p-perc-identity 1.0 --p-query-cov 0.94 --p-strand both --p-threads 20 --verbose
```

Default parameters were used for naive Bayes classification:
```
qiime feature-classifier classify-sklearn \
  --i-reads tmp.raw_repSeqs.qza \
  --i-classifier bold_anml_classifier.qza \
  --p-n-jobs 1 --p-reads-per-batch 2000 \
  --o-classification allASVs_NB_taxa.qza
```

Both the VSEARCH and Naive Bayes representative sequence taxonomic information outputs, [allASVs_VSp100c94_taxa.qza](https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/taxonomy/allASVs_VSp100c94_taxa.qza) and [allASVs_NB_taxa.qza](https://github.com/devonorourke/nhguano/raw/master/data/qiime_qza/taxonomy/allASVs_NB_taxa.qza), respectively, were exported as .tsv files (available as [allASVs_VS_taxa.tsv.g](https://github.com/devonorourke/nhguano/raw/master/data/taxonomy/allASVs_VS_taxa.tsv.gz) and [allASVs_NB_taxa.tsv.gz](https://github.com/devonorourke/nhguano/raw/master/data/taxonomy/allASVs_NB_taxa.tsv.gz)). These files were used as inputs in the [seqProcessing.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/seqProcessing.R) script, to filter remaining ASVs (focusing on only the centroids after clustering at 98.5%, called "OTUs" hereafter) for a minimum amount of taxonomic information.

## Filtering data for bat and arthropod-specific representative sequences
All filtering took place within the [seqProcessing.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/seqProcessing.R) script. A brief summary of these actions and their resulting outputs are highlighted below:

- Sequence counts per ASV per sample were aggregated for representative centroids (OTUs) - see [allSamples_OTUtable_long.csv.gz](https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/allSamples_OTUtable_long.csv.gz).
- Mock-specific ASVs were removed from the dataset; these particular sequences are listed in the file ['prevalentMockASVs.txt.gz'](https://raw.githubusercontent.com/devonorourke/nhguano/master/data/fasta/prevalentMockASVs.txt.gz). In addition, only New Hampshire guano samples were retained (thus dropping negative and positive control samples).
- All bat-associated sequences classified were identified separately with naive Bayes and VSEARCH, with both read abundances and sample occurrence summarized in [Table S4](https://github.com/devonorourke/nhguano/blob/master/supplementaryData/tableS4_batHost_summary.csv) of the manuscript. These data indicate that nearly all of our bat-classified COI data are derived from _Myotis lucifugus_ (Little brown bat), with at least 578 samples having at least one OTU classified to _M. lucifugs_ (578 for VSEARCH-classified OTUs, 579 for naive Bayes).
- For both VSEARCH and naive Bayes-classified OTUs, we separately filtered the dataset to retain only those OTUs with taxonomic family information assigned to the Phylum "Arthropoda". Thus, an OTU may be included that lacked genus or species labels, provided it retained an unambiguous family label. To select a final list of OTUs, we first retained VSEARCH-classified OTUs that fit these filtering criteria (representing near exact matches). Among VSEARCH-classified OTUs that did not pass this filtering threshold, we then selected from the naive Bayes-classified OTUs that met the same standard, if possible. In all, 559 OTUs were selected using the VSEARCH method, and 2,627 from naive Bayes.
- Ambiguous species labels were removed for taxonomic labels like "sp.", "nr.", or "gr.". Additionally, any alphanumeric suffixes in species labels were removed, thus a label like "Enicospilus purgatusDHJ02" was truncated to "Enicospilus purgatus". These data are available in the ['allTrueSamps_OTUtable_long_wTaxa.csv.gz'](https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/allTrueSamps_OTUtable_long_wTaxa.csv.gz) file.
- Among all OTUs remaining across all samples, we next noramlized samples using a method of scaling with ranked subsampling using the [SRS](https://peerj.com/articles/9593/) function ['SRS'](https://cran.r-project.org/web/packages/SRS/index.html). We retained only those samples with a per-sample minimum of 1,000 arthropod-classified reads. The resulting file ['min1kseqs_Samps_OTUtable_long_wTaxa.csv.gz'](https://github.com/devonorourke/nhguano/raw/master/data/text_tables/otu_tables/min1kseqs_Samps_OTUtable_long_wTaxa.csv.gz) was used for all subsequent diversity analyses processed in the [diversityAnalyses.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/diversityAnalyses.R) script.


# Diversity analysis summaries

## Overview
Diversity analyses were conducted using the entire collection of all samples and all sites, as well as selected samples at particular locations and dates. The majority of these investigations were completed using the [diversityAnalyses.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/diversityAnalyses.R) script. In addition, geographic location and land cover information at each site are shown in [Figure 1](https://raw.githubusercontent.com/devonorourke/nhguano/master/figures/figure1_nhmap_wBarplotCoverclass.png) of the manuscript, created using the R script [mapplots.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/mapplots.R). The pest analysis was completed using the R script [pest_work.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/pest_work.R).

Among all sampling locations and dates (899 samples total), we evaluated:  
- The most frequently detected arthropod orders
- The most frequently detected arthropod sequence clusters (OTUs), grouped by shared genus labels
- Potential taxa classified as pests by the US Forest Service (USFS) or US Department of Agriculture (USDA)  

Three separate analyses investigated how bat diets varied with temporal and spatial factors as follows: comparing effect of bat diet on samping date only, by comparing samples collected at a single site across multiple sampling windows in a single year; comparing effects between multiple sites and sampling windows in a single year; comparing effects of multiple sites at a single sampling date between multiple years.

- A single site, single year, multiple sampling windows (81 samples total):
  - Site: Fox State Forest, Hillsboro, New Hampshire USA (FOX)
  - Year: 2016
  - Sampling windows: 3 through 8
- Multiple sites, single year, multiple sampling windows (331 samples total):
  - 7 sites:
    - Cornish, New Hampshire, USA (COR)
    - Epsom, New Hampshire, USA (EPS)
    - Epsom, New Hampshire, USA (EPS)
    - Fox State Forest, Hillsboro, New Hampshire USA (FOX)
    - Holderness, New Hampshire, USA (HOL)
    - Maple Hill barn, Antrim, New Hampshire, USA (MAP)
    - Penacook, New Hampshire, USA (PEN)
  - Year: 2016
  - Sampling windows: 4 through 6
- Multiple sites, multiple years, single sampling windows (81 samples total):
  - 3 sites:
    - Cornish, New Hampshire, USA (COR)
    - Hopkington, New Hampshire, USA (HOP)
    - Maple Hill barn, Antrim, New Hampshire, USA (MAP)
  - Years: 2015 and 2016
  - Sampling window: 6 only

> Prefixes for site names match those presented in [Figure 1](https://raw.githubusercontent.com/devonorourke/nhguano/master/figures/figure1_nhmap_wBarplotCoverclass.png) and are applied to all other figures and tables throughout the manuscript:


  For estimating changes in richness and diet composition over time, samples were grouped into 37-day sampling windows. This retained the greatest number of samples in the smallest range of time shared across the most locations that generated sufficient sequence data. The particular sampling window periods are as follows:

| sampling window | start date | end date |
| --- | --- | --- |
| 3 | March 16 | April 22
| 4 | April 23 | May 29
| 5 | May 30 | July 05
| 6 | July 06 | August 11
| 7 | August 12 | September 17
| 8 | September 18 | October 24
| 9 | October 24 | November 30


## New Hampshire-wide analyses (all samples)
Analyses focused on those samples filtered to retain at least 1,000 arthropod-classified sequences, as previously described in the [Classifying OTUs and filtering samples](#classifying-otus-and-filtering-samples) section. 

### Frequently detected arthropod orders
### Frequently detected OTUs
### Pest analyses

## Spatiotemporal analyses (select samples)
### Single site, single year, multiple sampling windows
### Multiple sites, single year, multiple sampling windows
### Multiple sites, multiple years, single sampling window
