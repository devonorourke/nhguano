# Contamination overview
Raw sequence data was processed as described in the [sequence processing](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md) document. Demultiplexed paired-end read data from multiple libraries were trimmed with Cutadapt, denoised with DADA2, and merged into a single QIIME-formatted ASV table and fasta file. These data&mdash;the QIIME-formatted [ASV table](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) and the [ASV sequence](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/repSeqs/tmp.raw_repSeqs.qza) files merged from all libraries&mdash;served as input for these analyses.

We generated sequencing data from 9 Illumina libraries. These data consisted of true samples, positive biological mock samples, and negative controls (DNA extraction blanks that were amplified in conjunction with true samples). Because negative control samples (herein termed NTCs, for Negative Template Controls) generated some sequencing data we explored whether there was evidence for pervasive reagent contamination, specific sequencing batch contamination, or particular amplicon sequence variants (ASVs) that were likely contaminants. We focused on contamination originating in two ways:
1. **Platform-based contamination** Also known as "index-bleed", or "cross-talk", these data are a result of the sequencing of the libraries whereby the indexes assigned to one sample are mistakenly assigned to a different sample. As a result, the ASVs present in one sample are derived from a different sample. Because we included a biological mock community in every library, we were able to empirically evaluate the amount of cross talk for each library.
2. **Wet-bench contamination** can arise from the various molecular workflows conducted by the researcher in preparing the libraries: DNA extraction and PCR amplification are the principal sources of potential contamination in our experiments. DNA was extracted from samples using a 96-well plate format that has multiple opportunities for contamination: the initial bead-beating step requires the user to place individual guano samples in a 96-well plate with liquid reagents and shake the plate up to 30 Hz for 20 minutes. It's possible that the silicone mat used to secure and separate individual samples within the wells leaks liquid among the wells. It's also possible that some of this liquid and debris can be mixed between wells during the step following the shaking where the user removes the silicone mat. In addition, multiple rounds of pipetting various supernatants into and out of the wells poses opportunities for multichannel tips to drag across non-target wells. We later amplified our DNA via PCR, which required additional pipetting and added another layer of reagent (principally primers and PCR master mix) contamination. Finally, reagent contamination from the kits directly (also known as "kitomes") are possible, though less likely for our system given that these are arthropod amplicons and most kitomes are suspected to have microbial contaminants.

Our contamination evaluation explores the distributions of read abundances per sample and per ASV (i.e. how many sequence counts were present per sample? per ASV among all samples?), the distribution of unique ASVs per sample, and the prevalence of individual ASVs among all samples. We applied the R library [decontam](https://github.com/benjjneb/decontam) to evaluate the likelihood of ASV contamination, as describe below:

# Negative control samples generally produce few ASVs or sequence counts per sample

While a complete [script for this workflow](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R) defines all programs and code required for the following analysis, we wanted to provide a more extended visualization and commentary here.

Data was processed in RStudio v-1.2.1335 with R v-3.5.3. We imported the [ASV table](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) with the [qiime2R](https://github.com/jbisanz/qiime2R) package, and created a [phyloseq](https://joey711.github.io/phyloseq/) object that contained sequence counts per ASV per sample, as well as [metadata](https://github.com/devonorourke/nhguano/blob/master/data/metadata/allbat_meta.csv) information such as the sequencing library (`SeqBatch`) a sample was derived, as well as the DNA plate the guano sample was extracted (`DNAplate`). The `SampleType` was defined as either **sample** (guano sample), **ncontrol** (NTCs), or **mock** (biological mock sample).

While there were 10,797 unique sequence variants identified across all nine libraries, the majority of these were observed in just a single sample (**7,641** ASVs were singletons - a sequence variant detected in just one sample). Relatively few ASVs were observed in many samples (one indication of the lack of particular ASVs contamination): we observed just **39** ASVs in at least 100 samples (note that many sequencing runs contained ~400 samples). This is important in our understanding of contamination in general: we included 205 negative control samples among the 9 sequenced libraries, and most (183 samples) generated at least one sequence following DADA2 denoising. But this does not in and of itself indicate that these samples were contaminated with guano! We would expect some degree of cross-talk from the sequencing platform (for extended discussion of the various contamination sources and more links on the subject, check out [this blog post](http://fiererlab.org/2018/08/15/garbage-in-garbage-out-wrestling-with-contamination-in-microbial-sequencing-projects/)). What we want to understand isn't whether or not an ASV is observed in a **ncontrol** sample. What we're interested is whether the prevalence of a particular ASV is more likely to be detected in a control than a true sample. Furthermore, we'd like to better understand the distribution of read abundances of those ASVs in **ncontrol** samples: very low abundances are potentially a result of cross-talk. They may also be a result of very low-level contamination, but given that our diversity estimates include measures that incorporate read abundances, these rare and low abundant sequences are likely not have a major impact in our downstream analyses.

Therefore, we want to determine if:
1. The abundant ASVs are more likely to be observed in the **ncontrol** than guano **sample** `SampleType` data
2. The **ncontrol** samples generated relatively more or less sequence counts than true samples
3. The relative proportion of sequence counts in **mock** samples are of the _expected_ sequence types versus any unexpected ASVs

We first partitioned the samples into their respective `SampleType` groups (**sample**, **ncontrol**, and **mock**) to determine how the per-sample ASV prevalence (number of ASVs detected per sample) and per-sample read abundances were distributed. As seen in the figure below, the observed richness of ASVs (y axis) and the number of raw sequence counts (x axis) is typically far greater among true guano **samples** than the majority of **ncontrol** samples. In addition, the biological **mock** samples produce relatively high read abundances (reflecting the fact that purified samples as opposed to guano samples often amplify better) and have most sequence richness compared to many other guano samples. This is expected because mock samples were supposed to only have ~ 20 unique sequence variants.

![imagehere: contam_eval_allSampls_Counts-ASVs_scatterplot](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/contam_eval_allSampls_Counts-ASVs_scatterplot.png)

A few negative control samples are concerning in that they produce large numbers of total reads per sample, however we see that the vast majority of negative control samples are producing relatively few total reads, and are generally containing fewer than 10 unique ASVs (contrasted relative to the true samples, which contain many instances of samples with greater sequence counts). In other words, the negative controls represent the left tail of the read depth and ASV prevalence distributions. In fact, if we filter _out_ all samples that have fewer than 500 total sequences per sample, we find that the majority of negative control samples are discarded: just **24** of the original **~200** control samples generated at least 500 reads. The following plot below labels a few of the remaining NTC samples:
> all points shown are from same dataset; we simply dropped the samples with less than 500 total reads

![imagehere: contam_eval_allSampls_Counts-ASVs_scatterplot_filt500ReadMin](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/contam_eval_allSampls_Counts-ASVs_scatterplot_filt500ReadMin.png)

The labels for each negative control reflect the DNA plate number and DNA well position. For example, `negoro14D04` indicates that this particular NTC was extracted from plate **14** in well **D04**. It's possible that these few outlier NTCs were the result of trace amounts of bat guano being loaded into their respective wells: a piece of guano is the size of a large rice grain, and it's quite plausible among the >2000 guano samples that a few pieces broke into a well and were not detected by the researcher. Interestingly, _plate 14_ has two of the top 5 most abundant negative control samples (and there were 38 unique DNA plates with control samples), and perhaps this is an example in which the bead beating portion of the extraction process contributed to some contamination. One strong indication of this would be whether the community composition of the observed sequence variants are more similar to each other on this plate than other NTC samples. We'll explicitly test this in this workflow later.

Overall, this broad perspective illustrates that while we have many negative control samples producing some amount of sequence data, the vast majority of these **ncontrol** samples produce very few reads or ASVs. In fact, among the 9 libraries of data, the 183 NTC samples with sequence data analyzed constitute just 0.43% of the overall sequence data, despit accounting for about 6% of the total number of samples processed. Thus, if contamination is persistent, it is generally at very low abundances.

# Relatively few ASVs are identified by Decontam as likely contaminants

One approach to classifying a particular ASV as a suspected contaminant would be to highlight those sequence variants that are more prevalent in control samples than in true samples. The [decontam paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2) suggests that the reason is because the lack of competing DNA in control samples would allow for particular ASVs to be repeatedly amplified. In other words, the lack of variation in control samples would tend to produce the same ASVs among all control samples, while the highly variable sequence variants among true samples would tend to compete for sequence space, and thus a particular ASV would not be identified in proportionally as many true samples as in a control.   

The decontam function creates a score by generating a chi-square statistic from a 2x2 presence-absence table for each feature. The score statistic reflects the one-sided tail of the chi-square distribution, where smaller decontam scores indicate the contaminant model (an ASV is likely a contaminant) is more likely. We followed the authors suggestion of testing multiple thresholds with which an ASV would be flagged as a contaminant; a threshold of 0.5 would indicate ASVs that are more prevalent in a higher proportion of negative controls than the guano samples. As depicted in the histogram below, we found that the distribution of these decontam scores reflected that few sequence variants were likely contaminants:

![imagehere:decontam_prevHist_basic](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/decontam_prevHist_basic.png)

There are indeed _some_ ASVs in the left tail of that distribution, indicating that the program has targeted a few sequence variants as likely contaminants. We can see how many are flagged as contaminants at a given threshold in the following table:

![imagehere:decontam_prevTable_basic](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/decontam_prevTable_basic.png)

A few notes about this table:
1. Those marked as `FALSE` may be a situation where an ASV is present in both negative and true samples, but is better reflected by the non-contaminat model. It's also possible, however, that the ASV isn't detected in a negative control sample at all. These certainly wouldn't be contaminant ASVs then! In fact, there are **2,762** ASVs that are private to the true samples (that is, they are not detected in any negative control sample).
2. Those marked as `TRUE` may be a situation where an ASVs is present in both negative and true samples, but is better reflected by the contaminant model. It's also possible, however, that an ASV isn't detected in a true sample at all. These certainly wouldn't be contaminants _either_ because they aren't in our samples of concern! There was, in fact, just 1 ASV marked as `TRUE` by the decontam program that fit this description (that the ASV was a suspected contaminant by the program, but was not present in any of our guano samples), with just 6 reads among 2 samples.

Notably, the Decontam package also contains a parameter to control for batch effects. In our case we wanted to test for both DNA plate ("DNAplate") and sequencing ("SeqBatch") batch effects. We found that there were similar distributions of the left tail of Decontam scores for each batch type, but unique profiles to the right of the distribution:

![imagehere:decontam_prevHist_batch](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/decontam_prevHist_batch.png)

The overall number of ASVs identified as suspected contaminants were fairly similar among the `DNAplate` and `SeqBatch` batch types:

![imagehere:decontam_prevTable_batch](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/decontam_prevTable_batch.png)

Nevertheless, the distributional differences for each batch type had me wondering whether or not the _same_ ASVs were flagged as contaminants or not. We wouldn't expect that to necessarily be true: as noted by the Decontam authors, this program is not intended to model cross-talk, which is most likely the driving force behind any `SeqBatch` effect. Instead, the `DNAplate` model is better suited to model the ASV prevalence within a set of samples that shared a common DNA extraction condition. However, the statistical power to detect such differences is diminished when grouping by this batch type: there were between 1-10 (mean 5.42, median 6) NTC samples per DNA plate analyzed.

The following table demonstrates how many ASVs are shared or are distinct to a given set. With three different _'batch'_ types to consider in the Decontam model ("basic" == no batch, "DNAplate" for the DNA extraction plate type, and "SeqBatch" for the Illumina sequencing run) there seven distinct sets (listed in the table as a "Group"). We find that most ASV identified as a contaminant are generally shared among the three different batch types at a given threshold:

![imagehere:TableOf_ASVcounts_byThreshold_byGroup.png](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/TableOf_ASVcounts_byThreshold_byGroup.png)

Because there were distinct ASVs being flagged depending on the batch type, I wanted to explore _which_ ASVs those where, and how frequently they were detected among all the samples. I chose a threshold of 0.2 for the ASVs identified as contaminants in either the `SeqBatch` or `DNAplate` batch method. The following plot illustrates that the particular batch type we're employing makes a big difference if we were going to blindly remove a particular ASV from the dataset that was flagged by Decontam:

![decontam_Reads-ASVs_ContamOrNot](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/decontam_Reads-ASVs_ContamOrNot.png)

A few observations from this figure:
1. The most prevalent ASVs in our dataset are routinely identified as contaminants by the DNAplate method, but are typically not by the sequencing batch method. This makes me suspicious that there is not a sufficient statistical power to support the DNAplate-based batch filter. For example, it could be that a partiuclar ASV is marked as a contaminant in 3 of 40 DNA plates, but not so in the other 37. Nevertheless, we're seeing some ASVs in hundreds of samples, and these are likely _not_ contaminants: they're the kinds of sequence variants we expect in the bat diet (each ASV was classified and assigned taxonomic identity, explained in the last section of this document).
2. The mock ASVs sequenced were often flagged as likely contaminants themselves. Given that these samples were not part of the DNA extraction process, the more abundant samples are only a result of sequencing cross-talk (and we know Decontam isn't built for that). Given that the mock samples often generated among the greatest per-sample read abundances per Illumina run, it's not surprising that we might be seeing the mock ASVs showing up in negative control samples at low abundances. This is explored in more detail in the next section.
3. There are negative control samples that are _not_ identified as contaminants. These are clearly not of concern. It's unclear to me how an ASV would be unique to a control sample and not present in any mock sample or guano sample (perhaps it's a chimeric read that DADA2 missed, or a sequencing error that wasn't filtered out?), but these are not of a concern with our data and will be filtered out.

Because the Decontam package works by a presence-absence detection method, I wondered whether these ASVs being marked as contaminants were potentially overly sensitive, and whether accounting for differences in read abundances may be worth investigating. A negative control sample may have low read abundances whereas a true guano sample may have higher read abundances, for instance, yet these low reads are still "detected" in many control samples. As mentioned in the beginning of this document, we expected contaminants to occur more often in negative control samples than true samples, but we're not accounting for any sort of read abundances. If it turns out that most negative controls have very low sequence abundances of these ASVs that are highly prevalent, then perhaps these aren't very concerning. Likewise, it may be that when we visualize the per-sample ASV prevalence (like in the plot above) and find that a few ASVs have large read abundances, it could be that just one or two negative control samples generated that majority of those sequences. Thus, I wanted to examine the read counts on a per-sample basis, and selected a few ASVs that were highly prevalent, and marked either as contaminants or not among the two batch methods used in Decontam. As you can see in the plot below, the ASVs that are highly prevalent in our dataset are often considered contaminants (i.e. `TRUE` in the plot below) when using the `DNAplate` batch method in Decontam, but are `FALSE` when using the `SeqBatch` method:

![imagehere:decontam_Prevalence_ncontrolANDsample](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/decontam_Prevalence_ncontrolANDsample.png)

I've highlighted 9 different ASVs above:
- The contamination status is different for 3 ASVs depending on the batch type (ASVs 2, 5, and 7)
- Three ASVs are both `FALSE` (not identified as contaminant): ASVs 15, 18, and 20
- Three ASVs are both `TRUE` (identified as contaminants in both batch types): ASVs 8, 14, and 16

It's interesting to see that there is a generally linear trend among the presence of an ASV in both control and true samples. This indicates that we are more likely to detect an ASV in a control sample when it is found in many true samples. On one hand, that would be expected if it was a pervasive contaminant, but interestingly, these ASVs are rarely detected in many negative control samples. Among the 382 ASVs detected in these control samples also present in true samples, just 13 occur in at least 10 negative control samples (out of a possible 183 control samples with sequence data), suggesting that a study-wide contaminant is highly unlikely. Rather, it's likely those commonly observed ASVs are already in the 96-well plate and have a greater likelihood of being picked up in a negative control well during DNA extraction - indeed, that's exactly what we discovered, given that these 13 ASVs in at least 10 negative control samples are sequence variants. Here's a brief table summarizing the number of samples each of these prevalent ASVs are detected, according to sample type:

> the table below shows the number of samples in which a particular ASV is detected:
![imagehere: prevalent_ASVs_inGreatesNum_negControlSamps](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/prevalent_ASVs_inGreatesNum_negControlSamps.png)

Using those same 9 ASVs mentioned above, let's look at the distributions of the per-sample read abundances each of those ASVs in the among negative control samples relative to true samples:

> note the differences in the Y-axis for sequence counts between negative control and true samples!

![imagehere:decontam_selectASVabundance-perSample_contamComparison](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/decontam_selectASVabundance-perSample_contamComparison.png)

1. There doesn't appear to be a large difference in read abundances among the ASVs that are suspected as contaminants using the `DNAplate` batch method in Decontam, but not identified as contaminants using the `SeqBatch` class as the batch factor (top row of "A" in panel, purple square, ASVs 2,5,7). The sequence variants that were `TRUE` for both like ASVs 8 or 14 (middle row, blue square) appear to have a few more detections than the ASVs that were considered not to be indicators of contamination for either batch group (bottom row, red square), but compared with the sheer number of true samples ("B" panel in figure) these differences are questionable. 
2. Among control samples, we rarely detect many samples with high read abundances, but these same ASVs frequently are detected in higher amounts in true samples. Overall, we tend to see ASVs in our negative controls with the highest prevalence when they are also among the most prevalent ASVs in true samples. What's interesting about the above plot is that it suggests that for those ASVs which were represented in many controls samples like ASV2 (57 samples), ASV5 (33 samples, and ASV7 (30 samples), we rarely see many samples with high read abundances, but there are a few specific samples that generate the majority of these reads. If there were consistent and persistent contamination, I would have suspected that ASVs would have generated more similar patterns of sequence counts. However, elevated ASVs randomly increased in single negative control samples is exactly what you'd expect in situations where an empty well has no alternative DNA tempaltes to compete with during DNA amplification.

Overall, the analyses thus far have suggested that:
1. ASVs are most prevalent in negative control samples when they are also highly prevalent among true samples.
2. Negative controls rarely produce substantial read abundances collectively (per sample) or individually per sequence variant.
3. There are rare cases in which a particular ASV generates a substantial number of reads or is highly prevalent, but even among those instances the number of samples that have substantial read counts is very rare.
4. There is little evidence that there are specific sequencing run (SeqBatch) or DNA extraction (DNAplate) batch effects. Most ASVs identified by Decontam as contaminants are shared among both of these batch methods, though the DNAplate batch type generated the most number of suspected ASVs.


We wanted to next explore whether the composition of the sequence variants was associated with these batch effects. This required rarefying the data to conduct a series of distance estimates using nonphylogenetic and phylogenetic metrics.

# Contamination events are likley localized to individual plates, random across a given plate, and not pervasive within a plate
## Generating the required data: rarefying data, calculating distances, and PCoA
We'll first rarefy our data prior to any diversity evaluations. QIIME 2 has an [alpha rarefaction function](https://docs.qiime2.org/2019.4/plugins/available/diversity/alpha-rarefaction/) whereby we can visualize how the observed richness of a sample (or other alpha diversity metrics) changes with sampling depth. Because this analysis is focused on looking at how negative control samples may or may not associate with a particular group (i.e. DNA extraction plate, or Sequencing batch) we care more about preserving as many NTC samples as possible. As a result, we're first going to filter out all other non-NTC samples from the `tmp.raw_table.qza` initially imported in the [decontam R script](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R), then run the function on those NTC samples for a range of sampling depths. We'll then choose a single depth, rarefy _all_ samples at that depth, then use that rarefied table to conduct the distance estimates.  

We know that most NTCs have low abundances, so we'll investigate the observed richness of ASVs in a narrow range: 100 to 1000 sequence counts:

> `$TABLE` refers to the [tmp.raw_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) table
> `$META` refers to the QIIME-formatted metadata file [qiime_allbat_meta.tsv](https://github.com/devonorourke/nhguano/blob/master/data/metadata/qiime_allbat_meta.tsv)

```
## Filtering data to retain only NTC samples
qiime feature-table filter-samples \
  --i-table $TABLE \
  --m-metadata-file $META \
  --p-where "SampleType='ncontrol'" \
  --o-filtered-table tmp.neg_table.qza

## Summary of remaining samples:
qiime feature-table summarize --i-table tmp.neg_table.qza --o-visualization tmp.neg_sumrytable.qzv

## Alpha diversity visualization
qiime diversity alpha-rarefaction \
--i-table tmp.neg_table.qza --o-visualization tmp.neg_alphaviz.qza\
--p-metrics observed_otus --p-max-depth 1000 --p-min-depth 100
```

The [tmp.neg_sumrytable.qzv](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/tmp.neg_sumrytable.qzv) file can be loaded into the [QIIME viewer online](view.qiime2.org) and there is an interactive feature that illustrates how the number of features and samples are lost as a result of the sampling depth. The [tmp.neg_alphaviz.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/tmp.neg_alphaviz.qza.qzv) file contains a different summary of the number of observed ASVs at each of the specified sampling depths.  

The `tmp.neg_sumrytable.qzv` file (visualized in the table below) suggests that a sampling depth between 400-600 sequences will retain a balance between preserving the most number of sequence variants and the greatest number of samples. However this is only relative to the negative control samples; clearly, we'd want to have a higher number of sequences for the true samples.

| Sampling depth (bp) | % ASVs remaining | # Samples remaining |
| -------------- | ------ | --------- |
| 100 | 7% | 78 |
| 200 | 8.3% | 47 |
| 300 | 9.5% | 36 |
| 400 | 11.7% | 33 |
| 500 | 10.6% | 24 |
| 600 | 11.3% | 21 |
| 700 | 11.3% | 18 |
| 800 | 9.1% | 13 |
| 900 | 8.8% | 11 |
| 1000 | 8% | 9 |

We can see that the overall observed richness doesn't change much across that range of 400-600 bases in the `tmp.neg_alphaviz.qza` file. Here's a screenshot of what that file looks like in the viewer:

![imagehere:contam_alpharareViz](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/contam_alpharareViz.png)

Given these data, we selected a sampling depth of 500 reads in attempts to preserve diversity without losing too many NTC samples. We next rarefied the _entire_ dataset at that depth in QIIME. We also created another summary visualization of the remaining table to clarify which samples were retained in the analysis.

```
qiime feature-table rarefy \
  --i-table $TABLE \
  --p-sampling-depth 500 \
  --o-rarefied-table contam_rfyd_table.qza

qiime feature-table summarize --i-table contam_rfyd_table.qza --o-visualization contam_rfyd_sumrytable.qzv
```

As we can see in the [contam_rfyd_sumrytable.qzv](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/contam_rfyd_sumrytable.qzv) file, there are 1,515 samples remaining at this sampling depth, 24 of which are NTCs. The resulting [contam_rfyd_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/contam_rfyd_table.qza) ASV table contains the rarefied data, but we're going to filter out any samples with just a single ASV before we can use the data to calculate distances and ordinate the data:

```
qiime feature-table filter-samples --p-min-features 2 \
  --i-table contam_rfyd_table.qza --o-filtered-table contam_rfyd_filtd_table.qza

qiime feature-table summarize --i-table contam_rfyd_filtd_table.qza --o-visualization contam_rfyd__filtd_sumrytable.qzv
```

As we can see in the [contam_rfyd_filtd_sumrytable.qzv](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/contam_rfyd_filtd_sumrytable.qzv) visualization file, the filtered [contam_rfyd_filtd_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/contam_rfyd_filtd_table.qza) ASVtable now has just **1,498** samples instead of the original **1,515** samples. Notably, all 24 NTCs remain.  

Next, we're going to be using Unifrac to estimate community compositional differences between samples, so we need to create a tree using the ASV sequences in the dataset. We'll again rely on a QIIME function to [align ASVs to build a tree](https://docs.qiime2.org/2019.4/plugins/available/phylogeny/align-to-tree-mafft-fasttree/). This is done with the original [tmp.raw_repSeqs.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/repSeqs/tmp.raw_repSeqs.qza) QIIME-formatted fasta file equivalent.

> `$READS` refers to the `tmp.raw_repSeqs.qza` file linked above

```
## create tree with FastTree
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences $READS \
--p-n-threads 24 \
--i-sequences tmp.raw_repSeqs.qza \
--p-n-threads 24 \
--o-alignment raw.ASV_alignment.qza --o-masked-alignment raw.ASV_alignment_masked.qza \
--o-tree raw.ASVtree_unrooted.qza --o-rooted-tree raw.ASVtree_rooted.qza
```

The rooted [raw.ASVtree_rooted.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/trees) file is used for the Unifrac distance estimates next.

QIIME has a pair of beta diversity functions for [non-phylogenetic](https://docs.qiime2.org/2019.4/plugins/available/diversity/beta/) and [phylogenetic](https://docs.qiime2.org/2019.4/plugins/available/diversity/beta-phylogenetic/) metrics. We're using four metrics here in this comparison:
- Dice-Sorensen (a.k.a. _Observed otus_): unweighted abundance, unweighted phylogenetic
- Bray-Curtis: weighted abundance, unweighted phylogenetic
- unweighted Unifrac: unweighted abundance, weighted phylogenetic
- weighted Unifrac: weighted abundance, weighted phylogenetic  

Each distance estimate is calculated as follows:
> `$TABLE` refers to the `contam_rfyd_filtd_table.qza` file
> `$TREE` refers to the `raw.ASVtree_rooted.qza` file

```
## non-phylogenetic
qiime diversity beta --i-table $TABLE --p-metric dice --o-distance-matrix contam_ds_dist.qza
qiime diversity beta --i-table $TABLE --p-metric braycurtis --o-distance-matrix contam_bc_dist.qza

## phylogenetic
qiime diversity beta-phylogenetic --i-table $TABLE --i-phylogeny "$TREE" --p-metric unweighted_unifrac --o-distance-matrix contam_uu_dist.qza
qiime diversity beta-phylogenetic --i-table $TABLE --i-phylogeny "$TREE" --p-metric weighted_unifrac --o-distance-matrix contam_wu_dist.qza
```

We then use each of these distance estimates in a Principal Components Analysis:
> `$PCOADIR` refers to the [output director where each PCoA artifact is found](https://github.com/devonorourke/nhguano/data/qiime_qza/pcoa/contam_evals) - the `*pcoa.qza` files

```
qiime diversity pcoa --i-distance-matrix contam_ds_dist.qza --o-pcoa ds_pcoa.qza
qiime diversity pcoa --i-distance-matrix contam_bc_dist.qza --o-pcoa bc_pcoa.qza
qiime diversity pcoa --i-distance-matrix contam_uu_dist.qza --o-pcoa uu_pcoa.qza
qiime diversity pcoa --i-distance-matrix contam_wu_dist.qza --o-pcoa wu_pcoa.qza
```

These PCoA `.qza` files were imported into the [R script for this workflow](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R) to produce the subsequent ordinations presented in the following section. The `*dist.qza` objects were used to conduct a permutational analysis of variance using the [QIIME adonis plugin](https://docs.qiime2.org/2019.4/plugins/available/diversity/adonis/), which is a wrapper for the [Vegan Adonis R script](http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html) that allows us to assess the proportion of variation associated with a series of main effect variables. For our analysis we examined the model of (y ~ SampleType * SeqBatch * DNAplate).  

> `$META` refers to the QIIME-formatted metadata file [qiime_allbat_meta.tsv](https://github.com/devonorourke/nhguano/blob/master/data/metadata/qiime_allbat_meta.tsv)
> `$ADONISDIR` refers to the directory with the [output of the Adonis directory]((https://github.com/devonorourke/nhguano/data/qiime_qzv/contam_analysis/adonis)) containing the `.qzv` visualization files

```
qiime diversity adonis \
  --i-distance-matrix contam_ds_dist.qza --o-visualization "$ADONISDIR"/contam_ds_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
qiime diversity adonis \
  --i-distance-matrix contam_bc_dist.qza --o-visualization "$ADONISDIR"/contam_bc_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
qiime diversity adonis \
  --i-distance-matrix contam_uu_dist.qza --o-visualization "$ADONISDIR"/contam_uu_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
qiime diversity adonis \
  --i-distance-matrix contam_wu_dist.qza --o-visualization "$ADONISDIR"/contam_wu_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
```

## Diversity analyses
PERMANOVA reports for each of the four distance metrics are available as `.qzv` files in the [Adonis directory of the Repo]((https://github.com/devonorourke/nhguano/data/qiime_qzv/contam_analysis/adonis)) and can be viewed in the [QIIME viewer online](view.qiime2.org). We find that while there are significant main effects for all three groups (**SampleType**, **SeqBatch**, and **DNAplate**), the strength of these relationships are extremely low with a single exception: the main effect of DNAplate. This finding suggests that when contamination does persist, it is a limited to the samples within the single plate itself, thus project-wide pervasive contamination is not of a great concern.  

The per-sample dissimilarities were subsequently ordinated and visualized in the following plots. Both present the first two principle component axes for each of the four distance measures. Each plot shows negative control samples (purple color) and positive control samples (grey), with the proportion of variation captured by each fo the first two PCs shown on axis labels:

![imagehere:contam_pcoa_4metric_bySeq](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/contam_pcoa_4metric_bySeq.png)

Little variation is explained by Sequencing Run (by library), though it is strongest for weighted-UniFrac (~31.4%). Because sequencing batches generally had samples that were extracted in similar locations and dates we would expect there to be some minor relationship to sequencing batches even among negative controls. However when sequence abundances are taken into account, or phylogenetic information is added, or both, very little clustering by sequencing run group is evident. With so little overall variation explained by any of the metrics, we find no reason to be concerned about contamination due to the sequencing batch.  

The next plot highlights how negative control samples relate to their respective positive controls with respect to the DNA plate a sample was extracted from (again, purple indicates a negative control sample, and a gray sample is a guano sample):  

![imagehere:contam_pcoa_4metric_byDNA](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/contam_pcoa_4metric_byDNA.png)

We see that negative control and true guano samples are more frequently clustered by DNAplate numbers than with sequencing batch numbers. This is precicely what we would expect if there was a minor amount of variation occurring during DNA extraction - negative control samples that looked more like _other_ plates would be an indication of reagent contamination occurring across multiple extraction experiments and would be even more concerning. Thus we'd expect both guano and negative control samples to cluster together given that the guano samples were generally extracted in batches related to the site and week they were obtained. Nevertheless, while some samples are more similar to each other within the same DNAplate, the similarity of multiple negative control samples within a single DNAplate often vary in the same ordination space. If the entirety of an extraction was contaminated, we would expect all the negative control samples to cluster together more tightly to each other than they do to the other true samples, yet we don't have any strong evidence for that.

To show one example of this, we selected a single DNAplate (**DNAplate 33**) and illustrate that the specific well location of a negative control sample often does not match up in composition to the expected neighboring wells:

![imagehere:contam_pcoa_4metric_byDNAwel33only](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/contam_pcoa_4metric_byDNAwell33only.png)

These data for just a single DNAplate indicate two further pieces of evidence suggesting that we do not need to be concerned about particular ASVs contaminating entire DNAplates:
1. We do not see that the NTC samples are clustering together, thus there are no specific ASVs that appear to be dominating the negative control composition  
2. The NTCs do not necessarily contain similar compositions to all the neighboring wells where they were extracted. If there was local contamination, we would expect that a given NTC sample would be similar to the surrounding wells. For example, the NTC sample in well location **B02** would be surrounded by plate well locations: A01, B01, C01, A02, C02, A03, B03, and C03. We find some evidence for that (ex. C02) but other instances where it's not apparent at all (ex. C01, A03, B03).

Collectively these analyses point to random and infrequent instances where some minor amount of contamination may have occurred during DNA extraction, we expect that these instances are rare and will likely not contribute to any batch effects when analyzing the various sets of samples that were extracted across different plates. Furthermore, we find little evidence for any sequencing group bias. Overall we see no reason to remove any ASVs from these samples, and will simply remove the negative control samples from subsequent analysis.


# Mock community samples suggest there is some degree of cross-talk, but it is extremely low
We included a positive control in each sequencing run. These biological mock samples consisted of about 20 unique sequence variants spanning 10 distinct arthropod Orders as described in Michelle Jusino's [paper](https://doi.org/10.1111/1755-0998.12951). We were interested in using these mock samples to evaluate two different contamination features:
1. Reagent contamination. We do not expect the same contaminants identified in the Decontam workflow (see above) to be present in the mock samples because these mock samples were not extracted in conjunction with negative controls or true samples (the mock samples consist of equimolar pools of plasmids containing the individual sequence variants).  In addition, mock samples were amplified using similar primer constructs as the guano samples, except these reactions took place using batches of reagents that were separate from those used in the true samples.  
2. Cross talk. We _would_ expect some unexpected sequence variants to be present in the mock samples because of the sequencing process itself. These would typically be low abundance reads in the mock sample, though the particular unexpected ASVs in these mocks would most likely be derived from the most highly abundant ASVs in other samples.  

We applied the following code to perform the necessary tasks to complete this evaluation:
1. Subset the original ASV table (consisting of DADA2-filtered reads from all samples) to create a  mock-only ASV table  
2. Use QIIME's [quality control]() plugin to assign ASVs in the mock ASV table as either "expected" or "unexpected". We aligned all ASV sequences observed in each mock samle to a fasta file containing the known mock sequences, [partialCOImock.fasta](https://raw.githubusercontent.com/devonorourke/nhguano/master/data/qiime_qza/mock/partialCOImock.fasta), and considered any sequence within 99% identity as "expected", while any other sequence as "unexpected". We then created a list of the expected and unexpected ASVid's, and added applied that information within the [decontam_efforts.R](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R) script to generate the data necessary to create the resulting visualizations (last section, starting on line 660 of the R script). In addition, the samples that were verified as expected mock samples were ultimately removed from analysis (as mock sample ASVs can also result in cross-talk _into_ true samples).  
> `$ALLTABLE` refers to the initial [ASV table](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) imported in this workflow containing all DADA2-processed samples
> `$ALLSEQS` refers to the initial [ASV fasta](https://github.com/devonorourke/nhguano/data/qiime_qza/repSeqs/tmp.raw_repSeqs.qza) file associated with the ASV table
> `$META` refers to the QIIME-formatted metadata file [qiime_allbat_meta.tsv](https://github.com/devonorourke/nhguano/blob/master/data/metadata/qiime_allbat_meta.tsv)
> `$MOCKSEQ` refers to the QIIME-formatted mock sequences file [partialCOImock.seqs.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/mock/partialCOImock.seqs.qza)

```
## filtering data table to retain only mock samples and sequences
qiime feature-table filter-samples \
  --i-table "$TABLE" \
  --m-metadata-file "$META" --p-where "SampleType='mock'" \
  --o-filtered-table mock.table.qza

## filter the original fasta file to retain only the ASVs remaining in that mock table
qiime feature-table filter-seqs \
  --i-data tmp.raw_repSeqs.qza \
  --i-table mock.table.qza \
  --o-filtered-data mock.seqs.qza

## run VSEARCH to align all ASVs in mock samples to known mock reference sequences
MOCKSEQ=/mnt/lustre/macmaneslab/devon/guano/mockFastas/partialCOImock.seqs.qza
ALLSEQS=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/tmp.raw_repSeqs.qza
qiime quality-control exclude-seqs --p-method vsearch \
  --i-query-sequences $ALLSEQS \
  --i-reference-sequences $MOCKSEQ \
  --p-perc-identity 0.98 \
  --o-sequence-hits mock.expectSeqs.qza \
  --o-sequence-misses mock.notexpectSeqs.qza
```

The ASVs identified in the [mock.expectSeqs.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/mock/mock.expectSeqs.qza) file served as input for the final section of the [R script for this workflow](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R). We assigned any ASVid identified by the above alignment as a MockASV in the following plot.

![imagehere:contam_mockIndexBleed](https://github.com/devonorourke/nhguano/blob/master/figures/figs_for_docs/contam_mockIndexBleed.png)

The plot above indicates that ASVs not matching at least 99% similar to the expected mock sequences (the _unexpected_ sequences, listed in red in the plot above) were extremely low in terms of sequence abundance (per ASV, per library) relative to the expected mock sequences: the highest _unexpected_ ASV contained just 14 reads, with a mean of just 8.3 reads per unexpected ASV. This is exactly what we would expect: a very small proportion of cross talk in terms of read abundance. On a per-sample basis, there are 321 NH guano samples that contain at least one of these mock-associated ASVs, but there are always very few sequence counts attributed to any one of those expected mock sequences. While there was a maximum of 249 sequences detected for a single mock-associated ASV in a single sample, the average number with which any one mock-associated ASV is detected in a sample is just 2 sequences. Among all 25 of these ASVs, we detect a mean of 27 sequences and a median of 15 sequences. These data suggest that cross-talk may be routinely occurring at extremely low abundances from our mock samples into the broader population of sequence samples, a phenomenon that was expected given that these mock samples were typically more deeply sequenced per run than other guano samples.

In addition, we can clearly see that just a few of these ASVs are repeatedly detected within the mock samples themselves, and these are very likely our expected mock sequences. Just 24 ASVs are present in at least 5 of 10 mock samples, and there are just 25 distinct ASVs with at least 200 reads per sample. In fact, these sequence variants represent those directly from the mock sample and should be removed from the final dataset. We also queried these 25 ASV sequences with the nr Database in NCBI using their online BLAST application; as expected, nearly all sequences were matches for the expected taxonomies, and the two instances where they were not observed to be matches were when no suitable match was found (thus they were not the _incorrect_ assignment). The list of these [25 ASV identifiers is found here](https://github.com/devonorourke/nhguano/data/fasta/prevalentMockASVs.txt).

Collectively we find that  there is clear evidence for per-library cross-talk. though the abundances with which those sequence-based contamination events occur is very low. This suggests that we are best suited to utilize diversity metrics that incorporate abundance information in our analyses, otherwise we will likely be inflating our alpha diversity estimates, and likely inflating the similarity between samples that are otherwise less related. We therefore proceed by removing all mock and negative control samples, as well as removing the 25 ASVs identified as likely mock sequences in out dataset.

# Final thoughts

These analyses provided us with the confidence to proceed with removing all negative control samples without discarding any particular ASVs asociated with those negative controls. However, the expected mock sequences were identified as being clearly associated with the mock samples, and to avoid having those ASVs be incorporated to inflate subsequent species richness metrics, we did recommend dropping those 25 ASVs. While this may reduce some potential taxa that are true biological diet components (as the mock community was itself made of invertebrate representative sequences that _might_ exist in a bat diet), this felt like the more conservative and appropriate tradeoff.

## Next steps in analysis
We next proceed to the main analysis portion of the experiment, beginnning by clustering the exact sequence variants (ASVs) into representative variants sharing (at least) 98.5% similarity (OTUs). All negative and positive (mock) control samples are discarded, as well as the specific mock-associated ASVs identified in this analysis. These OTUs are then classified using both alignment and kmer-based approaches, allowing the identification of bat-host species identification, and arthropod prey items. We then proceed to evaluating dietary richness and composition among samples collected in particular spatial and temporal groups. All of these analyses are further described in the [diversity workflow](https://github.com/devonorourke/nhguano/blob/master/docs/diversity_analyses.md) document. 
